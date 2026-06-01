#!/usr/bin/env python3

import argparse
import sys

import numpy as np
import tifffile


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an int16, uint16, or int32 TIFF label image to uint32 NumPy array."
    )
    parser.add_argument(
        "--input", required=True, help="Input int16, uint16, or int32 TIFF label image."
    )
    parser.add_argument(
        "--output", required=True, help="Output uint32 .npy label array."
    )
    parser.add_argument(
        "--allow-negative",
        action="store_true",
        help="Allow negative signed integer values and cast them to uint32.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    image = tifffile.imread(args.input)
    if image.dtype not in (np.int16, np.uint16, np.int32):
        sys.exit(f"Expected int16, uint16, or int32 input, found {image.dtype}.")

    if not args.allow_negative and np.any(image < 0):
        sys.exit("Input contains negative labels; refusing to wrap them to uint32.")

    np.save(args.output, image.astype(np.uint32, copy=False))


if __name__ == "__main__":
    main()
