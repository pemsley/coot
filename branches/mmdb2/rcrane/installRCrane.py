#!/usr/bin/env python
"""Install RCrane so that it runs every time Coot launches."""

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

import sys
import os
import os.path
from os.path import basename, abspath, dirname
import inspect
from shutil import copy

def main():
    """Install RCrane so that it runs every time Coot launches"""
    
    try:
        filename = determineFilename()
        checkFile(filename)
        launchLoc = writeFile(filename)
    except RCraneInstallError, error:
        print error
        
        if withinCoot():
            #if we're running in Coot, use GTK to display any error messages in a message box
            import gtk
            errorDialog = gtk.MessageDialog(type=gtk.MESSAGE_ERROR, buttons=gtk.BUTTONS_OK, message_format=error)
            errorDialog.set_title("RCrane Installation Error")
            errorDialog.run()
            errorDialog.destroy()
    else:
        successMessage = "RCrane installed successfully.  " + filename + " created."
        print successMessage
        
        if withinCoot():
            #if we're running in Coot
            
            #if RCrane hasn't yet been launched, launch it
            if not sys.modules.has_key("rcrane"):
                run_python_script(launchLoc)
                successMessage += "\n\nRcrane launched."
        
            #report successful install (and possibly launch)
            import gtk
            msg = gtk.MessageDialog(type=gtk.MESSAGE_INFO, buttons=gtk.BUTTONS_OK, message_format=successMessage)
            msg.set_title("RCrane Installation")
            msg.run()
            msg.destroy()


def determineFilename():
    """Determine the appropriate filename for an RCrane launch script
    
    ARGUMENTS:
        None
    RETURNS:
        None
    EXCEPTIONS:
        RCraneInstallError - if the appropriate filename cannot be determined
                             or if the .coot-preferences directory does not exist and cannot be created
    """
    
    #NOTE: this function is nearly identical to settingsFilename() in settings.py, but in order to import that
    #function in place of this one, we'd need to know the RCrane directory
    #and most of this function is dedicated to finding the RCrane directory, so we'd still need this
    #function in order to import that function (i.e. we're not actually going to save any code by doing it that way)
    
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
    
    #if we can find HOME, look for COOT_HOME (which is used in WinCoot <= 0.6.2)
    if homedir == None:
        try:
            homedir = os.environ["COOT_HOME"]
        except KeyError:
            pass
    
    if homedir == None:
        #if we can't find HOME or COOT_HOME, then the user is probably running this from the command line in Windows
        #(or possibly we're running WinCoot > 0.6.2, but that doesn't exist yet)
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
        raise RCraneInstallError("Can't determine home directory ($HOME or %COOT_HOME%).  Cannot install RCrane.")
        return
    
    
    #make sure the .coot-preferences directory exists, and try to create it if it doesn't
    preferencesDir = os.path.join(homedir, ".coot-preferences")
    
    if os.path.exists(preferencesDir):
        if not os.path.isdir(preferencesDir):
            raise RCraneInstallError("Could not create Coot preferences directory: " + preferencesDir + "\nA file of the same name already exists")
    else:
        try:
            os.mkdir(preferencesDir)
        except OSError, error:
            raise RCraneInstallError("Could not create Coot preferences directory: " + preferencesDir + "\nError: " + str(error))
    
    #return the name of the file to create
    return os.path.join(homedir, ".coot-preferences", "rcrane.py")
    

def checkFile(filename):
    """Make sure that it's safe to overwrite any file at filename
    
    ARGUMENTS:
        filename - the filename to check
    RETURNS:
        None
    EFFECTS:
        will create a backup of filename if neccesary (if filename contains something other than the standard launch script)
    EXCEPTIONS:
        RCraneInstallError - if a backup is necessary and cannot be completed
    """
    
    if not os.path.exists(filename):
        return
    
    try:
        input = open(filename)
    except IOError, error:
        raise RCraneInstallerError("Could not open " + filename + " for reading\nError: " + str(error))
    
    for curline in input:
        import re
        #figure out if this file contains anything other than the launch command (or blank lines)
        if not re.match(r"""^run_python_script\(["'].*launch.py["']\)\s*$""", curline) or re.match("^\s*$", curline):
            break
    else:
        return
    
    #if we get here, then the file contains something we didn't recognize, so create a backup of it
    backupFilename = filename + ".bak"
    
    #try to find an rcrane.py.bak filename that we can use
    for curSuffix in [""] + map(str, range(1,99)):
        #print "trying suffix:", curSuffix
        if not os.path.exists(backupFilename + curSuffix):
            break
    else:
        raiseRCraneInstallError("Could not find a suitable backup filename for the existing rcrane.py.\nPlease check " + dirname(filename))
    backupFilename += curSuffix
    
    try:
        copy(filename, backupFilename)
    except IOError, error:
        raise RCraneInstallError("Could not backup existing " + filename + "\nError: " + str(error))


def writeFile(filename):
    """Create an RCrane launch script
    
    ARGUMENTS:
        filename - the filename for the launch script
    RETURNS:
        launchLoc - the location of RCrane's launch.py
    EXCEPTIONS:
        RCraneInstallError - if the launch script cannot be created
    """
    
    try:
        out = open(filename, 'w')
    except IOError, error:
        raise RCraneInstallerError("Could not open " + filename + "\nError: " + str(error))
    
    
    #figure out where the RCrane directory is (i.e. the directory where this file is located)
    #since this script may be run from within Coot, we must use the inspect module instead of something simpler
    fileloc = inspect.currentframe().f_code.co_filename
    fileloc = abspath(fileloc)
    rcranedir = dirname(fileloc)
    
    launchLoc = os.path.join(rcranedir, "launch.py")
    
    out.write('run_python_script(r"' + launchLoc + '")\n')
        #the path must be a raw string if we are on Windows and any of the directories in the launchLoc path
        #start with an "n" or a "t" (or any other escapable character)
    out.close()
    
    return launchLoc


def withinCoot():
    """Figure out if we're running from within Coot
    
    ARGUMENTS:
        None
    RETURNS:
        None
    """
    
    #do this by looking for the coot_version function, which won't exist if we're running this script from the command line
    #we could also look at sys.argv[0] and make sure we're running coot-real
    
    try:
        coot_version()
    except NameError:
        return False
    else:
        return True


class RCraneInstallError(Exception):
    """Errors associated with the installation of RCrane"""
    
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value


main()