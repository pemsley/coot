#!/usr/bin/env python
"""Store and check RCrane settings"""

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

import os, os.path, re

__filename = None

rcranePath = os.path.dirname(os.path.abspath(__file__))


def settingsFilename():
    """Return the filename of the RCrane launch/settings script
    
    ARGUMENTS:
        None
    RETURNS:
        None
    EXCEPTIONS:
        RCraneSettingsError - if the appropriate filename cannot be determined
                              or if the .coot-preferences directory does not exist and cannot be created
    """
    
    #cache the filename results, since it can take a fair number of hard drive accesses to determine it
    global __filename
    if __filename is not None:
        return __filename
    
    homedir = None
    
    #check for a HOME environment variable first
    #this should only exist under a posix environment (i.e. Linux, OS X, msys)
    #This may break if we're running under Windows and the user has manually set a %HOME% environment variable
    #but that seems unlikely.  We could avoid that problem by checking os.name first (i.e. don't look for HOME
    #if os.name == 'nt'), but that would fail under MSys/MinGW
    try:
        homedir = os.environ["HOME"]
    except KeyError:
        pass
    
    #if we can find HOME, look for COOT_HOME (which is used in WinCoot)
    if homedir == None:
        try:
            homedir = os.environ["COOT_HOME"]
        except KeyError:
            pass
    
    if homedir == None:
        #if we can't find HOME or COOT_HOME, then the user is probably running this from the command line in Windows
        #in a last ditch effort, look for the WinCoot binary
        
        if os.path.getsize("C:\WinCoot\bin\coot") and os.path.getsize("C:\WinCoot\bin\coot-real.exe"):
            #first check C:\WinCoot, the default install location
            homedir == "C:\WinCoot"
        elif os.path.getsize("C:\Program Files\WinCoot\bin\coot") and os.path.getsize("C:\Program Files\WinCoot\bin\coot-real.exe"):
            #then check C:\Program Files\WinCoot, which could be a common alternate path
            homedir == "C:\Program Files\WinCoot"
        elif os.path.getsize("C:\Program Files (x86)\WinCoot\bin\coot") and os.path.getsize("C:\Program Files (x86)\WinCoot\bin\coot-real.exe"):
            #then check C:\Program Files (x86)\WinCoot (we're really grasping at straws here)
            homedir == "C:\Program Files (x86)\WinCoot"
        #if none of those paths work, then give up
        #short of searching the entire hard drive, we're not going to find the Coot installation directory
        
    if homedir == None:
        #if we still haven't found the appropriate directory, give up
        #this probably means we're running in Windows from the command line and WinCoot is installed in a non-standard directory
        raise RCraneSettingsError("Can't determine home directory ($HOME or %COOT_HOME%).")
        return
    
    
    #make sure the .coot-preferences directory exists, and try to create it if it doesn't
    preferencesDir = os.path.join(homedir, ".coot-preferences")
    
    if os.path.exists(preferencesDir):
        if not os.path.isdir(preferencesDir):
            raise RCraneSettingsError("Could not create Coot preferences directory: " + preferencesDir + "\nA file of the same name already exists")
    else:
        try:
            os.makedirs(preferencesDir)
        except OSError, error:
            raise RCraneSettingsError("Could not create Coot preferences directory: " + preferencesDir + "\nError: " + str(error))
    
    #now that we've figured out the filename, cache it
    filename = os.path.join(homedir, ".coot-preferences", "rcrane.py")
    __filename = filename
    
    #return the name of the file to create
    return filename


def addSetting(setting):
    """Adds the specified setting to the end of the RCrane configuration file
    
    ARGUMENTS:
        setting - the setting to add to the configuration file
    RETURNS:
        True if the setting was added succesfully
        False if the setting already existed
    """
    
    #read in the current settings file
    settingsFile = open(settingsFilename(), 'a+')
    
    curline = None
    for curline in settingsFile:
        #make sure the line doesn't already exist
        if re.match(setting, curline):
            settingsFile.close()
            return False
    
    #find out if the last line of the file ends with a newline
    if curline is None or curline[-1] == "\n":
        newline = ""
    else:
        newline = "\n"

    settingsFile.write(newline + setting + "\n")
    settingsFile.close()
    


class RCraneSettingsError(Exception):
    """Errors associated with setting or reading RCrane settings"""
    
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value
