#load the rcrane module
import rcrane

#create the menu entry in the Extensions menu
extMenu = coot_menubar_menu("Extensions")
rcraneLaunch = gtk.MenuItem("RCrane launch")
rcraneLaunch.connect("activate", lambda x: rcrane.createRCraneMenuWithCheck())
rcraneLaunch.show()
extMenu.insert(rcraneLaunch, len(extMenu)-1)

#destroy any variables that we've created, since the variables will be visible from the Coot Python scripting window
del extMenu
del rcraneLaunch