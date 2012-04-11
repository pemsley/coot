#!/usr/bin/env python
"""This module is a simple wrapper for RCrane modules to make RCrane easy to import and to avoid polluting Coot's Python namespace"""

# Copyright 2010, 2011, 2012 Kevin Keating
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

from traceGui import TraceGui as newTrace
from about import createAboutDialog
from rotamerize import pickAtoms as newRotamerize
from rotamerize import storeMenuItem as storeRotamerizeMenuItem
from rotamerize import setRotamerizeMaxNucleotides, getRotamerizeMaxNucleotides
from citationPopup import dontShowPopup as dontShowCitationPopup
from calcCoords import disablePhenixRestraints, enablePhenixRestraints, usingPhenixRestraints
from menu import createRCraneMenu

try:
    from debug import enableDebugging
except ImportError:
    #if we can't find the debug module, then just replace the enableDebugging function with a stub
    #(the debug module isn't included in the distributable)
    def enableDebugging(): return False

newTrace5to3 = lambda: newTrace(direction = 5)
    #a convenience function for users who are setting keyboard shortcuts
newTrace3to5 = newTrace
    #a convenience function to allow for consistant naming between different trace directions
newRotamerizeWithoutDensity = lambda: newRotamerize(ignoreDensity = True)
    #a convenience function to rotamerize without density