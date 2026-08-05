#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

import numpy as np
import pandas as pd
import fire


def get_n_readout(string):
    ss = string.split("Readout ")
    if len(ss) > 1:
        n = ss[1]
    else:
        n = -1  # when the name is not starting with "Readout", for example "Anchor"
    return int(n)


def find_max_n_readouts(table):
    Nmax = 0
    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            N = get_n_readout(table.iloc[i, j])
            if N > Nmax:
                Nmax = N
    return Nmax


def prep_averaged_spot_profiles(spot_profiles, readout_table, Nmax):
    avg_spot_profile = np.zeros(
        (spot_profiles.shape[0], 1, Nmax), dtype=spot_profiles.dtype
    )
    for n_r in range(Nmax):
        list_of_cyc_ch = []
        for i in range(readout_table.shape[0]):
            for j in range(readout_table.shape[1]):
                nn = get_n_readout(readout_table.iloc[i, j])
                if nn > -1:
                    if nn == n_r + 1:
                        list_of_cyc_ch.append(np.array([i, j]))
        for nl in range(len(list_of_cyc_ch)):
            avg_spot_profile[:, 0, n_r] += spot_profiles[
                :, list_of_cyc_ch[nl][1], list_of_cyc_ch[nl][0]
            ]
        avg_spot_profile[:, 0, n_r] = (
            avg_spot_profile[:, 0, n_r] / len(list_of_cyc_ch)
        ).astype(spot_profiles.dtype)
    return avg_spot_profile


def main(spot_profiles, path_readouts_csv):
    ## this function takes computed for all cycle-channels spot profiles and create an average version for each readout
    ## readouts_csv should contain N rows and M columns, where N - is number of cycles and M is number of channels
    ## cells should be filled like "Readout 1", "Readout 2" etc starting from 1 not 0
    readouts_table = pd.read_csv(path_readouts_csv)
    assert (
        spot_profiles.shape[1] == readouts_table.shape[1]
    ), "number of channels is different in npy file and in readouts csv!"
    assert (
        spot_profiles.shape[2] == readouts_table.shape[0]
    ), "number of cycles is different in npy file and in readouts csv!"
    Nmax = find_max_n_readouts(readouts_table)
    print("There is " + str(Nmax) + " readouts in total")
    spot_profiles_avg = prep_averaged_spot_profiles(spot_profiles, readouts_table, Nmax)
    return spot_profiles_avg, Nmax


if __name__ == "__main__":
    fire.Fire(main)
