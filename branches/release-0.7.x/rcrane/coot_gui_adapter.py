#!/usr/bin/env python
"""Allows functions from coot_gui to be imported."""

# Copyright 2012 Kevin Keating
# 
# Licensed under the Educational Community License, Version 2.0 (the
# "License"); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
# 
# http://www.osedu.org/licenses/ECL-2.0
# 
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS"
# BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
# or implied. See the License for the specific language governing
# permissions and limitations under the License.

#"import coot_gui" results in an error, so this module is required to retrieve
#functions that are defined in coot_gui

import os, sys
from os.path import exists, join
from coot import *

#search the Python path for coot_utils
for curpath in sys.path:
    abspath = join(curpath, "coot_gui.py")
    if exists(abspath):
        #when we find it, exec it
        execfile(abspath)
        break