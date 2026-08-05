#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute
# Author: Stanislaw Markarchuk

import numpy as np
import pandas as pd
from tifffile import TiffFile
from ome_types import from_xml


def extract_channel_names(img_path):
    with TiffFile(img_path) as tif:
        ome_meta = tif.ome_metadata
    ome = from_xml(ome_meta)
    ch_names_list = []
    for i in range(len(ome.images[0].pixels.channels)):
        ch_names_list.append(
            ome.images[0].pixels.channels[i].name.replace(" ", "-")
        )  # here I replace with "-" to not mess with "_" as separator for cycle/channel
    return ch_names_list


def prepare_codebook_ISS(
    codebook_path, gene_name="gene", cyc_name="cycle", ch_name="channel", separator="_"
):
    # firstly find number of cycles and channels from column names
    codebook = pd.read_csv(codebook_path)
    list_of_col_names = [s for s in codebook.columns if cyc_name in s]
    n_cycles = [
        int(col_name.split(cyc_name, 1)[1][0]) for col_name in list_of_col_names
    ]
    Ncycles = np.max(n_cycles)
    n_channels = [
        int(col_name.split(ch_name, 1)[1][0]) for col_name in list_of_col_names
    ]
    Nchannels = np.max(n_channels)
    Ngenes = codebook.shape[0]
    codebook_3d = np.zeros((Ngenes, Nchannels, Ncycles), dtype="uint8")
    print("Found " + str(Ncycles) + " cycles and " + str(Nchannels) + " channels")
    channel_order = []
    for ng in range(Ngenes):
        for nch in range(Nchannels):
            for ncyc in range(Ncycles):
                column_subname = (
                    cyc_name + str(ncyc + 1) + separator + ch_name + str(nch + 1)
                )
                column_name = [s for s in list_of_col_names if column_subname in s][0]
                codebook_3d[ng, nch, ncyc] = int(codebook[column_name][ng])
                if ng == 0 and ncyc == 0:
                    channel_order.append(column_name.split("_")[2])
    gene_list_obj = np.array(codebook[gene_name], dtype=object)
    return gene_list_obj, codebook_3d, Ngenes, Ncycles, Nchannels, channel_order


def prepare_spot_profile_ISS(
    spot_profile_path,
    Ncycles,
    Nchannels,
    order_channels_spot_profile,
    order_channels_codebook,
    start_cycle=2,
):
    spot_profile = np.load(spot_profile_path, allow_pickle=True)
    Nspots = spot_profile.shape[1]
    assert (Ncycles + start_cycle - 1) * Nchannels == spot_profile.shape[
        0
    ], "Dimensions of spot profile array is not matching codebook. Please check number of cycle, channels and also number of cycles to skip in spot profile array"
    spot_profile_3d = np.zeros((Nspots, Nchannels, Ncycles), dtype=spot_profile.dtype)
    for ns in range(Nspots):
        for nch in range(Nchannels):
            nch_spot_profile = order_channels_spot_profile.index(
                order_channels_codebook[nch]
            )
            for ncyc in range(Ncycles):
                spot_profile_3d[ns, nch, ncyc] = spot_profile[
                    (ncyc + start_cycle - 1) * Nchannels + nch_spot_profile, ns
                ]
    return spot_profile_3d


def delete_channel(codebook_3d, spot_profile_3d, channel_order, channels_to_delete):
    nums = [channel_order.index(element) for element in channels_to_delete]
    codebook_3d = np.delete(codebook_3d, nums, axis=1)
    spot_profile_3d = np.delete(spot_profile_3d, nums, axis=1)
    return codebook_3d, spot_profile_3d


def prepare_iss(
    codebook_path,
    spot_profile_path,
    image_path=None,
    channel_names=None,
    start_cycle=2,
    channels_to_delete=["DAPI"],
):
    assert image_path != "None", "Please specify path to the original image!"
    gene_list_obj, codebook_3d, Ngenes, Ncycles, Nchannels, channel_order = (
        prepare_codebook_ISS(codebook_path)
    )
    if image_path != None:
        channel_names = extract_channel_names(image_path)

    if channel_names == None:
        raise ValueError(
            "Please specify path to the original image or provide channel names!"
        )

    spot_profile_3d = prepare_spot_profile_ISS(
        spot_profile_path, Ncycles, Nchannels, channel_names, channel_order, start_cycle
    )
    if len(channels_to_delete) > 0:
        codebook_3d, spot_profile_3d = delete_channel(
            codebook_3d, spot_profile_3d, channel_order, channels_to_delete
        )
    return codebook_3d, spot_profile_3d, gene_list_obj, Ngenes
