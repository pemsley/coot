#!/usr/bin/env python
"""Ensure that Coot is able to properly run RCrane."""

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

from coot import svn_revision
from coot import monomer_restraints_py as monomer_restraints

def checkVersion():
    """Make sure that Coot is a new enough version
    
    ARGUMENTS:
        None
    RETURNS:
        True if Coot is SVN revision 3562 or newer (corresponds to the 0.6.2 release build)
        False otherwise
    """
    
    #this function could be (and previously was) implement by parsing coot_version(), which gives us more info
    #but testing svn_revision is far safer as it doesn't depend on the precise string formatting of coot_version()'s output
    
    if svn_revision() >= 3562: #3562 is the 0.6.2 release
        return True
    else:
        return False


def checkMonomerLibrary():
    """Make sure that the monomer library is up-to-date
    
    ARGUMENTS:
        None
    RETURNS:
        True if Coot is using an appropriate monomer library (the one included with Coot 0.6.2 or CCP4 6.2)
        False otherwise
    """
    
    if monomer_restraints("A")['_chem_comp'][3] == "RNA":
        return True
    else:
        return False