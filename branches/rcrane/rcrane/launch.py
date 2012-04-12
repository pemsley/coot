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

#if Coot is able to run RCrane, then load RCrane and create the menu
if preLaunchCheck.checkCootAndReportErrors():
    print "Launching RCrane"
    rcrane = imp.load_module("rcrane", None, rcranePath, ('', '', imp.PKG_DIRECTORY))
        #we use imp to import the rcrane module in case the directory name has a period in it
        #this also allows us to re-import the module later without knowing the directory name
    rcrane.createRCraneMenu()

#destroy any variables that we've created, since the variables will be visible from the Coot Python scripting window
del fileloc
del rcranePath
del preLaunchCheckPath
del preLaunchCheckFileObject
del preLaunchCheck
