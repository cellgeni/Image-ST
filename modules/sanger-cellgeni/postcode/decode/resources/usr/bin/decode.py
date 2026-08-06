#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

import numpy as np
import pandas as pd
import fire
from decoding_functions import train_and_save_model
from starfish.core.codebook.codebook import Codebook
import logging
import torch


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

VERSION = "0.2.0"


def decode(
    spot_locations_p: str,
    spot_profile_p: str,
    starfish_codebook_p: str,
    output_path: str,
    up_prc_to_remove: float = 99.95,
    num_iter: int = 60,
    zero_threshold=0.5,
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
        up_prc_to_remove: spots with any feature above this percentile (per
            feature, across all spots) are dropped before training. Pass 100 to
            disable this filtering entirely (e.g. for a small/pre-cleaned
            profile where the percentile itself may be unstable).
        num_iter: number of SVI training iterations.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the decoded spots and their locations.
    """
    starfish_book = Codebook.open_json(starfish_codebook_p)
    codebook_arr = np.array(starfish_book).transpose(0, 2, 1)
    # gene_list = np.array(starfish_book.target)
    # K = len(starfish_book.target)

    # float16's max finite value (~65504) is below raw 16-bit intensity range
    # (up to 65535) - any saturated pixel would silently overflow to inf here
    # and poison the log-normalization downstream into NaN. torch_format()
    # already upcasts to float32 before training, so keep full precision until
    # then rather than truncating range for no benefit.
    spot_profile = np.load(spot_profile_p).astype(np.float32)
    spot_locations = pd.read_csv(spot_locations_p)
    assert spot_locations.shape[0] == spot_profile.shape[0], (
        "Number of spots in spot_locations and spot_profile do not match"
    )

    if len(spot_locations.columns) == 2:
        spot_locations["spot_id"] = spot_locations.index
    elif len(spot_locations.columns) == 3:
        spot_locations.columns = ["spot_id", "y", "x"]
    else:
        raise ValueError("spot_locations_p should have 2 or 3 columns")

    # Check if there are too many zeros in the profile
    zero_proportion = np.mean(spot_profile == 0, axis=(1, 2))

    if np.all(zero_proportion > zero_threshold):
        logger.info("Too many zeros in the profile, skipping decoding.")
        torch.save({}, output_path)
        # spot_locations["Name"] = "high_zero_proportion"
        # spot_locations.to_csv(out_name, index=False)
        # return spot_locations

    logger.info(spot_profile)
    to_serialize_params = train_and_save_model(
        spot_profile,
        codebook_arr,
        print_training_progress=False,
        up_prc_to_remove=up_prc_to_remove,
        num_iter=num_iter,
    )

    # Serialize parameters, losses, and data_norm to disk
    torch.save(to_serialize_params, output_path)

    # Decode using postcode
    # try:
    #     out = decoding_function(
    #         spot_profile, codebook_arr, print_training_progress=False, **decoding_params
    #     )
    # except ValueError as e:
    #     logger.error(f"Decoding failed: {e}")
    #     spot_locations["Name"] = "decoding_failed"
    #     spot_locations.to_csv(out_name, index=False)
    #     raise e

    # # Reformat output into pandas dataframe
    # df_class_names = np.concatenate((gene_list, ["infeasible", "background", "nan"]))
    # df_class_codes = np.concatenate(
    #     (barcodes_0123_str, ["infeasible", "background", "NA"])
    # )
    # decoded_spots_df = decoding_output_to_dataframe(out, df_class_names, df_class_codes)

    # decoded_df_s = pd.concat([decoded_spots_df, spot_locations], axis=1)
    # decoded_df_s.to_csv(out_name, index=False)


if __name__ == "__main__":
    options = {"run": decode, "version": VERSION}
    fire.Fire(options)
