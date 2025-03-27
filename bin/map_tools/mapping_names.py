#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
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


def names_map(map_file):
    names_equiv = {}
    with open(map_file, "r") as input_map:
        for line in input_map:
            new_name, old_name = line.rstrip().split("\t")
            names_equiv[new_name.replace(">", "")] = old_name
    inv_names_equiv = {v: k for k, v in names_equiv.items()}

    return (names_equiv, inv_names_equiv)
