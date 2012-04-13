#!/usr/bin/env python
"""Functions for the RCrane menu."""

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

import gtk, os.path, imp
from coot_gui_adapter import coot_menubar_menu, add_simple_coot_menu_menuitem
import rcrane


def createRCraneMenu():
    """Create the RCrane menu
    
    ARGUMENTS:
        None
    RETURNS:
        True if the RCrane menu was created
        False if the RCrane menu already had entries and no new entries were created
    """
    
    #create a separate menu for RCrane (default)
    menu = coot_menubar_menu("_RCrane")
    
    #if the RCrane menu already has entries, then give up
    if len(menu.get_children()) > 0:
        return False

    add_simple_coot_menu_menuitem(menu, "New trace 5'->3'...", lambda x: rcrane.newTrace(direction=3))
    add_simple_coot_menu_menuitem(menu, "New trace 3'->5'...", lambda x: rcrane.newTrace(direction=5))

    
    #add_simple_coot_menu_menuitem(menu, "Rotamerize existing structure...", lambda x: rcrane.rotamerize())
    rotamerize = gtk.MenuItem("Rotamerize existing structure...")
    rotamerize.connect("activate", lambda x: rcrane.newRotamerize())
    rotamerize.show()
    menu.append(rotamerize)
    rcrane.storeRotamerizeMenuItem(rotamerize)
    
    #if Coot is newer than 3728 and has the regularize_zone_with_score_py function, then add a menu option for Rotamerize without density
    try:
        from coot import regularize_zone_with_score_py

        rotamerizeNoDensity = gtk.MenuItem("Rotamerize without density...")
        rotamerizeNoDensity.connect("activate", lambda x: rcrane.newRotamerize(ignoreDensity = True))
        rotamerizeNoDensity.show()
        menu.append(rotamerizeNoDensity)
        rcrane.storeRotamerizeMenuItem(rotamerizeNoDensity, ignoreDensity = True)
    except ImportError:
        print "Coot r3728 (0.7-pre) required for Rotamerize without density.\n\tRotamerize without density not available."
    
    add_simple_coot_menu_menuitem(menu, "About RCrane...", lambda x: rcrane.createAboutDialog())
    
    return True
