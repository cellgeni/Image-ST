#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

""" """

import fire
from spotiflow.model import Spotiflow
from imagetileprocessor import slice_and_crop_image
import numpy as np
import csv


def main(
    image_path: str,
    output_name: str,
    C: int,
    x_min: int,
    x_max: int,
    y_min: int,
    y_max: int,
    model_name: str = "general",
    Z: int = 0,
):
    crop = slice_and_crop_image(image_path, x_min, x_max, y_min, y_max, Z, C, 0)
    crop = np.squeeze(crop)
    assert len(crop.shape) == 2, "Crop should be a 2D array"
    model = Spotiflow.from_pretrained(model_name)
    peaks, details = model.predict(crop)

    with open(output_name, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["y", "x"])  # write column names
        if len(peaks) > 0:
            # Serialize peaks to disk as CSV
            peaks[:, 0] = peaks[:, 0] + y_min
            peaks[:, 1] = peaks[:, 1] + x_min
            writer.writerows(peaks)


if __name__ == "__main__":
    options = {"run": main, "version": "0.1.0"}
    fire.Fire(options)
