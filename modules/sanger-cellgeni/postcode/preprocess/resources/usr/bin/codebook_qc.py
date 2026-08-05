#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2025 Wellcome Sanger Institute

import pandas as pd
from starfish.core.codebook.codebook import Codebook
from starfish.types import Axes, Features
import numpy as np
import re


def to_starfish_codebook(
    codebook: pd.DataFrame,
    target_col: str,
    code_col: str,
    is_merfish: bool = False,
    n_channel: int = None,
) -> Codebook:
    """Build a starfish Codebook from a tabular (gene, code, ...) codebook.

    n_channel: number of real channels per round. For a non-MERFISH, multi-
        channel codebook (code digit selects which of n_channel channels is
        active per round), this MUST be passed explicitly - it cannot be
        reliably inferred from the code string's digit range alone, since e.g.
        a 2-channel selector code and a single-channel binary-presence code
        both only ever use digits {0, 1}. When not passed, falls back to the
        old max-digit-based heuristic (only correct for single-channel codes).
    """
    all_codes = codebook[code_col]
    all_codes_lengths = [len(str(code)) for code in all_codes]
    n_round = max(all_codes_lengths)
    assert len(set(all_codes_lengths)) == 1, (
        "All codes in the codebook must have the same length. "
        f"Found lengths: {set(all_codes_lengths)}"
    )
    all_codes_str = "".join(all_codes)
    max_digit = max(int(digit) for digit in all_codes_str)
    min_digit = min(int(digit) for digit in all_codes_str)
    is_zero_based = min_digit == 0
    if is_merfish:
        n_channel = 1
    elif n_channel is None:
        n_channel = max_digit if is_zero_based else max_digit + 1
    mappings = []
    for _, row in codebook.iterrows():
        mapping = {}
        mapping[Features.TARGET] = row[target_col]
        codeward = []
        for r, c in enumerate(str(row[code_col])):
            if is_merfish:
                # single channel per round; the digit itself is the presence value
                code_value = c
            else:
                # digit selects which channel is active this round - the digit's
                # numeric value is a channel index, not an intensity, so the
                # active channel is always marked present (1), regardless of
                # which digit (e.g. "0") happened to select it.
                code_value = 1
            codeward.append(
                {
                    Axes.ROUND.value: r,
                    Axes.CH.value: 0 if is_merfish else int(c),
                    Features.CODE_VALUE: code_value,
                }
            )
        mapping[Features.CODEWORD] = codeward
        mappings.append(mapping)
    return Codebook.from_code_array(mappings, n_round=n_round, n_channel=n_channel)


def filter_columns_by_regex(df: pd.DataFrame, pattern: str) -> pd.DataFrame:
    regex = re.compile(pattern)
    return df.loc[:, df.columns.to_series().apply(lambda x: bool(regex.match(x)))]


def load_tabular_codebook(codebook_path: str, codebook_code_col: str) -> pd.DataFrame:
    """
    Read codebook from a file.

    Args:
        codebook_path (str): A file path to the codebook.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the codebook.
    """
    if codebook_path.endswith(".xlsx"):
        codebook = pd.read_excel(codebook_path, dtype={codebook_code_col: str})
    elif codebook_path.endswith(".csv"):
        codebook = pd.read_csv(codebook_path, dtype={codebook_code_col: str})
    elif codebook_path.endswith(".tsv"):
        codebook = pd.read_csv(codebook_path, sep="\t", dtype={codebook_code_col: str})
    else:
        raise ValueError("Codebook file must be in .xlsx, .csv or .tsv format")
    return codebook.dropna()
    # codebook.to_json(f"{Path(codebook_p).stem}.json")


def filter_columns_by_regex(df: pd.DataFrame, pattern: str) -> pd.DataFrame:
    regex = re.compile(pattern)
    return df.loc[:, df.columns.to_series().apply(lambda x: bool(regex.match(x)))]


def qc_codebook(
    codebook, code_col: str, coding_col_prefix: str = "cycle\d_channel\d_*"
):
    df = filter_columns_by_regex(codebook, coding_col_prefix).astype(int)
    _, n_code = df.shape
    if "channel" not in coding_col_prefix:
        for i, code in enumerate(codebook[code_col]):
            assert code == "".join(
                df.iloc[i, :].values.astype(str)
            ), "ISS-like Codebook is not consistent with the coding columns"
    else:
        for i, code in enumerate(codebook[code_col]):
            current_code = np.zeros(n_code, dtype=int)
            for j, char in enumerate(code):
                current_code[j * n_code // len(code) + int(char)] = 1
            assert "".join(df.iloc[i, :].values.astype(str)) == "".join(
                current_code.astype(str)
            ), "Merifish-like Codebook is not consistent with the coding columns"
