#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2025 Wellcome Sanger Institute

import numpy as np
import pandas as pd
import fire
from avg_spot_profile import main as average_spot_profiles
import logging

import os
from codebook_qc import (
    load_tabular_codebook,
    qc_codebook,
    to_starfish_codebook,
    filter_columns_by_regex,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

VERSION = "0.2.0"


def reorder_profile(profile, channel_orders=[0, 1, 2, 3, 4], n_cycles=None):
    """
    Reorder hyperstack based on channel orders for each cycle

    Parameters:
    -----------
    profile : np.ndarray
        Input hyperstack array with shape (n_channels * n_cycle, ...) or (n_cycles * n_channels, ...)
    channel_orders : list, optional
        List containing channel orders for each cycle.
    n_cycles : int, optional
        Number of cycles in the hyperstack. If not provided, it will be inferred from the
        shape of the profile.

    Returns:
    --------
    np.ndarray
        Reordered hyperstack with correct channel ordering
    """
    assert (
        len(profile.shape) == 2
    ), "Profile should be a 2D array with shape (n_channels * n_cycles, n_features)"
    n_total_channels = profile.shape[0]
    channel_orders = np.array(channel_orders, dtype=int)
    n_channel = len(channel_orders)
    cycle_channel_mask = channel_orders != 0
    logger.info(cycle_channel_mask)
    if n_channel != n_total_channels:
        assert (
            n_total_channels % n_channel == 0
        ), f"Number of channels ({n_total_channels}) must be a multiple of the number of channel ({n_channel})"
        n_cycles = n_total_channels // n_channel if n_cycles is None else n_cycles
        channel_orders = np.vstack(
            [
                (i * n_channel + channel_orders) * cycle_channel_mask
                for i in range(n_cycles)
            ]
        )
        channel_orders = channel_orders[channel_orders != 0] - 1
    else:
        assert (
            n_cycles is not None
        ), "n_cycles should be set if complete channel_orders is provided"
        n_channel = n_total_channels // n_cycles
        # User gave the full channel orders, should be 0-based
        channel_orders = channel_orders
    logger.info(channel_orders)
    # Reorder the profile
    return (
        profile[channel_orders]
        .reshape(n_channel, n_cycles, profile.shape[1])
        .transpose(2, 0, 1)
    )


def decode(
    profile: str,
    tabular_codebook: str,
    out_starfish_codebook: str,
    out_reformatted_profile: str,
    readouts_csv: str = None,
    R: int = None,
    codebook_targer_col: str = "Gene",
    codebook_code_col: str = "code",
    coding_col_prefix: str = "cycle\d_channel\d_+",
    channel_orders: list = [0, 1, 2, 3, 4],
) -> pd.DataFrame:
    """
    Decodes spots using the Postcode algorithm.

    Args:
        profile (str): A file path to numpy array containing the spot profiles (N x C x R).
        tabular_codebook (str): Cortana-like codebook with only one channel and number of rounds (readouts)
        out_starfish_codebook (str): name of the output starfish codebook file
        out_reformatted_profile (str): name of the output reformatted profile file
        readouts_csv (str, optional): csv file with table which describes link between cycle-channels and readouts
        R (int): Number of rounds. Defaults to None.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the decoded spots and their locations.
    """
    is_merfish = readouts_csv and os.path.getsize(readouts_csv) != 0
    codebook = load_tabular_codebook(tabular_codebook, codebook_code_col)
    qc_codebook(codebook, codebook_code_col, coding_col_prefix)

    # Channels per round can't be inferred from the code digits alone (e.g. a
    # 2-channel selector code and a single-channel binary-presence code both
    # only ever use digits {0, 1}) - derive it from the codebook's own coding
    # columns instead, the same ones qc_codebook already validated against.
    n_channel = None
    if not is_merfish:
        n_round = codebook[codebook_code_col].str.len().max()
        n_coding_cols = filter_columns_by_regex(codebook, coding_col_prefix).shape[1]
        assert n_coding_cols % n_round == 0, (
            f"Number of coding columns ({n_coding_cols}) is not a multiple of "
            f"the code length ({n_round}) - cannot infer channels per round."
        )
        n_channel = n_coding_cols // n_round

    starfish_book = to_starfish_codebook(
        codebook,
        target_col=codebook_targer_col,
        code_col=codebook_code_col,
        is_merfish=is_merfish,
        n_channel=n_channel,
    )
    starfish_book.to_json(out_starfish_codebook)

    spot_profile = np.load(profile)

    if len(spot_profile.shape) == 2 and R:
        reordered_profile = reorder_profile(
            spot_profile, channel_orders=channel_orders, n_cycles=R
        )
    else:
        reordered_profile = spot_profile

    if is_merfish:
        processed_spot_profile, _ = average_spot_profiles(
            reordered_profile, readouts_csv
        )
    else:
        processed_spot_profile = reordered_profile
    # Convert codebook_arr to the form of barcodes_0123_str
    np.save(out_reformatted_profile, processed_spot_profile)


if __name__ == "__main__":
    options = {"run": decode, "version": VERSION}
    fire.Fire(options)
