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
Run the model
"""

from . import MODULE_NAME
from . import _solve as sol
from . import _log as log
from . import _meshing as mesh
from . import _plot as plot
from .__version__ import __version__

logger = log.logger


def run_model(m):
    """
    Run the complete model analysis

    Args:
        :m: instance of model
    """

    logger.info(f"===== {MODULE_NAME} {__version__} =====")
    ...
