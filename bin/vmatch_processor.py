#!/usr/bin/env python3

import argparse
import re
import sys

import pandas as pd


def process_vmatch_output(vmatch_file, result_file):
    """Process vmatch output to identify direct repeats"""

    repeats_counter = 0
    with open(vmatch_file, "r") as input_file, open(result_file, "w") as output_file:
        for line in input_file:
            if not line.startswith("#"):
                repeats_counter += 1
                lines = line.strip().split()
                ice_element = lines[1]
                start_1 = str(int(lines[2]) + 1)
                end_1 = str(int(lines[2]) + int(lines[0]))
                start_2 = str(int(lines[6]) + 1)
                end_2 = str(int(lines[6]) + int(lines[4]))
                DR = "\t".join(
                    [
                        ice_element,
                        start_1,
                        end_1,
                        start_2,
                        end_2,
                    ]
                )
                output_file.write(DR + "\n")

    return repeats_counter


def main():
    parser = argparse.ArgumentParser(
        description="Process vmatch output to identify direct repeats"
    )
    parser.add_argument("vmatch_file", help="Input vmatch output file")
    parser.add_argument(
        "output_file", help="Output TSV file with processed direct repeats"
    )
    args = parser.parse_args()

    try:
        num_repeats = process_vmatch_output(
            args.vmatch_file,
            args.output_file,
        )
        print(f"Successfully found {num_repeats} potential direct repeats")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
