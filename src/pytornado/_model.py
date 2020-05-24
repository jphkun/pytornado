#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------
# Copyright 2019-2020 Airinnova AB and the PyTornado authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ----------------------------------------------------------------------

# Author: Aaron Dettmann

"""
Frame model generator
"""

from numbers import Number

from mframework import FeatureSpec, ModelSpec

from ._run import run_model


class S:
    any_int = {'type': int}
    any_num = {'type': Number}
    pos_int = {'type': int, '>': 0}
    pos_number = {'type': Number, '>': 0}
    string = {'type': str, '>': 0}
    vector3x1 = {'type': list, 'min_len': 3, 'max_len': 3, 'item_types': Number}
    vector6x1 = {'type': list, 'min_len': 6, 'max_len': 6, 'item_types': Number}


# =================
# ===== MODEL =====
# =================
mspec = ModelSpec()

fspec = FeatureSpec()
fspec.add_prop_spec('segment_vertices', {'type': dict}, singleton=False, required=True, doc="Add a wing segment")
mspec.add_feature_spec('wing', fspec, singleton=False, required=True, doc="Add a wing")

# ===================
# ===== RESULTS =====
# ===================
rspec = ModelSpec()
mspec.results = rspec


# ===== MODEL =====
class Model(mspec.user_class):
    def run(self):
        super().run()
        run_model(self)
        return self.results
