#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import migfinder as mf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate the input.yaml for the CWL pipeline")
    parser.add_argument("-f", dest="fasta",
                        required=True,
                        help="Fasta file")
    args = parser.parse_args()

    # The output file will be stored in the same path as the provided fasta file
    mf.main(args.fasta)
