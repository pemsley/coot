#!/usr/bin/env python
"""Allows functions from coot_utils to be imported"""

# Copyright 2011 Kevin Keating
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

#"import coot_utils" results in an error, so this module is required to retrieve
#functions that are defined in coot_utils

import os, sys
from os.path import exists, join
from coot import *

use_gui_qm = False #coot_utils requires this variable to be defined

#search the Python path for coot_utils
for curpath in sys.path:
    abspath = join(curpath, "coot_utils.py")
    if exists(abspath):
        #when we find it, exec it
        execfile(abspath)
        break