
#load the rcrane module
import rcrane

def import_rcrane_wrapper():
   
   #create the menu entry in the Extensions menu
   extMenu = coot_menubar_menu("Calculate")
   rcraneLaunch = gtk.MenuItem("RCrane launch")
   rcraneLaunch.connect("activate", lambda x: rcrane.createRCraneMenuWithCheck())
   rcraneLaunch.show()
   extMenu.insert(rcraneLaunch, len(extMenu.get_children())-1)

   #destroy any variables that we've created, since the variables will be visible from the
   # Coot Python scripting window
   del extMenu
   del rcraneLaunch

   #check to see if the user has an old version of RCrane installed
   #if they do, stop it from running
   def checkForOldRCraneVersions():
       from rcrane.settings import settingsFilename
       rcraneSettingsFilename = settingsFilename()

       #try to open the RCrane settings file
       try:
           rcraneSettings = open(rcraneSettingsFilename, 'r')
       except IOError:
           #if the file doesn't exist, then the user doesn't have an old version of RCrane installed
           return False

       for curline in rcraneSettings:
           if curline[0:18] == "run_python_script(":
               rcraneSettings.close()
               break
       else:
           #the rcrane.py file exists, but it doesn't have a run_python_script( line,
           #so the user doesn't have an old version of RCrane installed
           rcraneSettings.close()
           return False

       #if we reach here, then the user has an old version of RCrane installed
       from shutil import move
       from os.path import dirname, join

       #move the user's rcrane.py to rcrane.py.bak, since it's no longer necessary
       try:
           move(rcraneSettingsFilename, join(dirname(rcraneSettingsFilename), "rcrane.py.bak"))
       except IOError:
           print("Cannot move ~/.coot-preferences/rcrane.py to ~/.coot-preferences/rcrane.py.bak")
           return False

       #let the user know what we just did
       noticeDialog = gtk.MessageDialog(type = gtk.MESSAGE_WARNING, buttons = gtk.BUTTONS_OK)
       noticeDialog.set_title("Old version of RCrane detected")
       noticeDialog.set_markup("An old version of RCrane was detected.  Coot now includes "
           + "RCrane 1.1, which may be accessed via Extensions -> RCrane launch.  To prevent the "
           + "old version of RCrane from running, your "
           + "~/.coot-preferences/rcrane.py file was moved to ~/.coot-preferences/rcrane.py.bak.")
       noticeDialog.run()
       noticeDialog.destroy()

       return True


   checkForOldRCraneVersions()
   del checkForOldRCraneVersions #destroy the function as well, since we don't need it after it's been run

# do this now if we load stuff from guile-gtk

# print '===== in rcrane loader: use_gui_qm is ', use_gui_qm

if (use_gui_qm != 2):
    import_rcrane_wrapper()
