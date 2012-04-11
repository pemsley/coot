# Copyright 2010, 2011, 2102 Kevin Keating
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

import inspect, sys, os.path, gtk
import imp

#first, find the location of this script (i.e. the launch.py file) so we know where the RCrane package is located
#note that we can't use __file__ because of the way Coot calls the script, so we use the inspect module instead

fileloc = inspect.currentframe().f_code.co_filename
fileloc = os.path.abspath(fileloc)
rcranePath = os.path.dirname(fileloc)

#import the preLaunchCheck module so we can make sure Coot is able to run RCrane
#we use imp to do the importing in case the RCrane directory name has a period in it
preLaunchCheckPath = os.path.join(rcranePath, "preLaunchCheck.py")
preLaunchCheckFileObject = open(preLaunchCheckPath, "r")
try:
    preLaunchCheck = imp.load_module("rcranePreLaunchCheck", preLaunchCheckFileObject, preLaunchCheckPath, ('.py', 'r', imp.PY_SOURCE))
finally:
    #if the import fails, we still want to close the file object
    if preLaunchCheckFileObject: preLaunchCheckFileObject.close()

#make sure Coot is able to run RCrane
if not(preLaunchCheck.checkVersion()):
    #if this version of Coot isn't new enough for RCrane, print an error message to the terminal and pop up an error dialog
    print "RCrane requires Coot version 0.6.2 or newer."
    print "Current Coot version is " + str(coot_version())
    print "RCrane not launched."
    
    errorDialog = gtk.MessageDialog(type = gtk.MESSAGE_ERROR, buttons = gtk.BUTTONS_OK)
    errorDialog.set_title("RCrane Error")
    errorDialog.set_markup("RCrane requires Coot version 0.6.2 or newer.  RCrane not launched.")    
    errorDialog.run()
    errorDialog.destroy()
    
elif not(preLaunchCheck.checkMonomerLibrary()):
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
    
else:
    print "Launching RCrane"
    rcrane = imp.load_module("rcrane", None, rcranePath, ('', '', imp.PKG_DIRECTORY))
        #we use imp to import the rcrane module in case the directory name has a period in it
        #this also allows us to re-import the module later without knowing the directory name
    rcrane.createRCraneMenu()

#destroy any variables that we've created, since the variables will be visible from the Coot Python window
del fileloc
del rcranePath
del preLaunchCheckPath
del preLaunchCheckFileObject
del preLaunchCheck

#errorDialog only exists if RCrane wasn't launched
try:
    del errorDialog
except NameError:
    pass