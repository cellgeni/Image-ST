#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

import numpy as np
import pandas as pd
import fire
from avg_spot_profile import main as average_spot_profiles
from decoding_functions import decoding_function, decoding_output_to_dataframe
import logging

# from prepare_ISS import prepare_iss
import os
from codebook_qc import load_codebook, qc_codebook, to_starfish_codebook
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

VERSION = "0.1.0"


# def ReadPrepCodebook_ISS(codebook_path):
#     '''
#     ISS stands for In Situ Sequencing. It is a method for encoding ISS data.
#     '''
#     #I consider single channel!
#     codebook_in = load_codebook(codebook_path)
#     codes = codebook_in['code']
#     n_genes = len(codes); n_rounds = len(str(codes[0]))
#     n_channel = -1
#     for channel in codebook_in.columns:
#         if 'channel' in channel:
#             channel_index = int(channel.split('_')[1].replace('channel', ''))
#             if channel_index > n_channel:
#                 n_channel = channel_index
#     n_channel -= 1 #remove DAPI/Hoechst channel
#     codebook_3d = np.zeros((n_genes, n_channel, n_rounds), dtype =  'uint8')
#     for ng in range(n_genes):
#         for nr in range(n_rounds):
#             for nch in range(n_channel):
#                 codebook_3d[ng, nch, nr] = int(codebook_in['channel_' + str(nch+1)][ng][nr])
#                 # codebook_3d[ng, 0, nr] = int(str(codes[ng])[nr])
#     gene_list_obj = np.array(codebook_in['gene'], dtype = object)
#     return gene_list_obj, codebook_3d, n_genes


# def ReadPrepCodebook_MER(codebook_path, N_readouts):
#     '''
#     MER stands for Multiplex Error Robust. It is a method for encoding MERFISH data.
#     '''
#     codebook_in = load_codebook(codebook_path)
#     print(codebook_in)
#     n_genes = codebook_in.shape[0]; n_rounds = N_readouts
#     codebook_3d = np.zeros((n_genes, 1, n_rounds), dtype =  'uint8')
#     for ng in range(n_genes):
#         for nr in range(n_rounds):
#             column_name = 'Readout_' + str(nr+1)
#             codebook_3d[ng, 0, nr] = int(codebook_in[column_name][ng])
#     gene_list_obj = np.array(codebook_in['gene'], dtype = object)
#     return gene_list_obj, codebook_3d, n_genes


def decode(
    spot_locations_p: str,
    spot_profile_p: str,
    codebook_p: str,
    out_name: str,
    readouts_csv: str = None,
    keep_noises=True,
    min_prob=0.95,
    R: int = None,
    codebook_target_col: str = "Gene",
    codebook_code_col: str = "code",
    coding_col_prefix: str = "cycle\d_channel\d_+",
) -> pd.DataFrame:
    """
    Decodes spots using the Postcode algorithm.

    Args:
        spot_locations_p (str): A file path to pandas DataFrame containing the spot locations.
        spot_profile_p (str): A file path to numpy array containing the spot profiles (N x C x R).
        codebook (str): Cortana-like codebook with only one channel and number of rounds (readouts)
        out_name (str): name of the output file
        readouts_csv (str, optional): csv file with table which describes link between cycle-channels and readouts
        keep_noises (bool, optional): Whether to keep spots that were classified as 'background' or 'infeasible'.
        min_prob: [0,1] - value of minimum allowed probability of decoded spot
            Defaults to True.
        R (int): Number of rounds. Defaults to None.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the decoded spots and their locations.
    """
    is_merfish = readouts_csv and os.path.getsize(readouts_csv) != 0
    codebook = load_codebook(codebook_p, codebook_code_col)
    qc_codebook(codebook, codebook_code_col, coding_col_prefix)
    starfish_book = to_starfish_codebook(
        codebook,
        target_col=codebook_target_col,
        code_col=codebook_code_col,
        is_merfish=is_merfish,
    )
    starfish_book.to_json(f"{Path(codebook_p).stem}_starfish_codebook.json")
    codebook_arr = np.array(starfish_book).transpose(0, 2, 1)
    gene_list = np.array(starfish_book.target)
    K = len(starfish_book.target)

    spot_profile = np.load(spot_profile_p)

    if len(spot_profile.shape) == 2 and R:
        # if the spot_profile is two dimensional, it is assumed that the spot_profile is in
        # the shape of (n_channel*n_cycle, n_spot). Then reshape it.
        n_ch, n_spots = spot_profile.shape
        # fine DAPI/Hoechst channel indexes and remove them from the profile
        n_chs_per_cycle = n_ch // R
        coding_mask_pre_cycle = np.ones(n_chs_per_cycle)
        coding_mask_pre_cycle[0] = 0  # remove DAPI/Hoechst channel
        coding_ch_mask = np.array(list(coding_mask_pre_cycle) * R, dtype=bool)
        spot_profile = spot_profile[coding_ch_mask].reshape(
            R, n_chs_per_cycle - 1, n_spots
        )
        spot_profile = spot_profile.transpose(2, 1, 0)
        # np.save(f"{stem}_reshaped_spot_profile.npy", spot_profile)
        print(
            spot_profile.shape,
            "\n",
            spot_profile[0],
            type(spot_profile[0]),
            "\n",
            spot_profile[0].dtype,
        )

    if is_merfish:
        spot_profile, _ = average_spot_profiles(spot_profile, readouts_csv)
    # if os.path.getsize(readouts_csv) != 0: # MERFISH-like data, this file should be provided
    #     spot_profile, N_readouts = average_spot_profiles(spot_profile, readouts_csv) # Average is chosen over max for MERFISH-like profiles
    #     gene_list, codebook_arr, K = ReadPrepCodebook_MER(codebook_p, N_readouts)
    # else:
    #     codebook_arr, spot_profile, gene_list, K = prepare_iss(codebook_p, spot_profile_p, **prepare_iss_kwargs)

    spot_locations = pd.read_csv(spot_locations_p)
    if len(spot_locations.columns) == 2:
        spot_locations["spot_id"] = spot_locations.index
    elif len(spot_locations.columns) == 3:
        spot_locations.columns = ["spot_id", "y", "x"]
    else:
        raise ValueError("spot_locations_p should have 2 or 3 columns")

    assert spot_locations.shape[0] == spot_profile.shape[0], (
        "Number of spots in spot_locations and spot_profile do not match"
    )
    # Decode using postcode
    out = decoding_function(spot_profile, codebook_arr, print_training_progress=False)

    # Reformat output into pandas dataframe
    df_class_names = np.concatenate((gene_list, ["infeasible", "background", "nan"]))
    if is_merfish:
        barcodes_0123_str = ["".join(k) for k in codebook_arr[:, 0, :].astype(str)]
    else:
        barcodes_0123_str = [
            "".join(np.argmax(k, axis=0).astype(str)) for k in codebook_arr.astype(str)
        ]
        # barcodes_0123_str = ["".join(k) for k in codebook_arr.astype(str)]
    # barcodes_0123 = codebook_arr[:,0,:]
    # barcodes_AGCT = np.empty(K, dtype='object')
    # for k in range(K):
    #     barcodes_AGCT[k] = ''.join(list(barcodes_0123[k, :].astype(str)))
    df_class_codes = np.concatenate((barcodes_0123_str, ["NA", "0000", "NA"]))
    decoded_spots_df = decoding_output_to_dataframe(out, df_class_names, df_class_codes)

    decoded_df_s = pd.concat([decoded_spots_df, spot_locations], axis=1)
    decoded_df_s = decoded_df_s[decoded_df_s["Probability"] > min_prob]
    decoded_df_s["Z"] = 0

    if keep_noises:
        decoded_df_s.to_csv(out_name, index=False)
    else:
        # Remove infeasible and background codes
        decoded_df_s[
            ~np.isin(decoded_df_s["Name"], ["background", "infeasible"])
        ].reset_index(drop=True).to_csv(out_name, index=False)


if __name__ == "__main__":
    options = {"run": decode, "version": VERSION}
    fire.Fire(options)
