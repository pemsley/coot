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

import gtk
from coot import svn_revision, coot_version
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


def checkCootAndReportErrors():
    """Run both pre-launch checks (Coot version and monomer library) and report any errors that occur
    
    ARGUMENTS:
        None
    RETURNS:
        True if the Coot version is recent enough to run RCrane and Coot is using an appropriate monomer library
        False otherwise
    EFFECTS:
        If Coot fails a pre-launch check, an error will be printed to standard output and a pop-up
        will appear explaining the error
    """
    
    if not checkVersion():
        #if this version of Coot isn't new enough for RCrane, print an error message to the terminal and pop up an error dialog
        print "RCrane requires Coot version 0.6.2 or newer."
        print "Current Coot version is " + str(coot_version())
        print "RCrane not launched."
        
        errorDialog = gtk.MessageDialog(type = gtk.MESSAGE_ERROR, buttons = gtk.BUTTONS_OK)
        errorDialog.set_title("RCrane Error")
        errorDialog.set_markup("RCrane requires Coot version 0.6.2 or newer.  RCrane not launched.")    
        errorDialog.run()
        errorDialog.destroy()
        return False
        
    elif not checkMonomerLibrary():
        #if Coot's monomer library isn't correct, print an error message to the terminal and pop up an error dialog
        print "Monomer library is out of date.  RCrane not launched."
        
        errorDialog = gtk.MessageDialog(type = gtk.MESSAGE_ERROR, buttons = gtk.BUTTONS_OK)
        errorDialog.set_title("RCrane Error")
        errorDialog.set_markup("""RCrane was not launched because Coot appears to be using an out-of-date monomer library.  Coot cannot refine RNA unless the monomer library is updated.

This problem is typically caused by running Coot on a computer with an old version of CCP4 (&lt;6.2) installed.  If this is the case, there are two solutions:

1) Set $COOT_REFMAC_LIB_DIR to
/path/to/Coot/share/coot/lib before running Coot.
Under Bash, this can be done with the command:
  export COOT_REFMAC_LIB_DIR=/path/to/coot/share/coot/lib
with /path/to/coot replaced by the appropriate Coot directory.

2) Upgrade your CCP4 installation to 6.2 or newer.""")
        
        errorDialog.run()
        errorDialog.destroy()
        return False
        
    else:
        return True