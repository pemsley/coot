# extensions.py
# Copyright 2007, 2008 by Bernhard Lohkamp
# Copyright 2006, 2007, 2008 by The University of York
# Copyright 2015 by Medical Research Council
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

# import pygtk, gtk, pango

import coot
import gi
gi.require_version("Gtk", "4.0") 
from gi.repository import Gtk, GObject
from gi.repository import Gio
from gi.repository import GLib
import coot_gui_api # built into coot binary, not an extension
import time
import numbers
import coot_utils
import coot_gui
import gui_add_linked_cho
import refmac
import gui_prosmart
import fitting
import shelx
import shelx_extensions
import find_baddies
import jligand_gui
import parse_pisa_xml
import ncs

def test_func_for_menus(menumodel: Gio.MenuModel, menu_label: str, submenu_label: str):
    print("Finding", menu_label, submenu_label)
    print("menumodel dir", dir(menumodel))
    if menumodel == False:
        print("********** in my_menumodel_get_menu_item_index()", submenu_label, "null menumodel")
    n = menumodel.get_n_items()
    quoted_menu_label = "'" + menu_label + "'"
    quoted_submenu_label = "'" + submenu_label + "'"
    for i in range(n):
        mmml = menumodel.get_item_attribute_value(i, "label", GLib.VariantType.new("s"))

    idx = 0
    for i in range(n):
        idx = i
        thing = menumodel.get_item_link(idx, "label")
        print("thing", thing, "dir thing", dir(thing))
        #  name = thing.get_name()
        # print("name:", name)
        # l = thing.get_property("label")
        # v = thing.get_value()
        # print("v:", v)

    
menumodel = coot_gui_api.main_menumodel()
if menumodel:
    test_func_for_menus(menumodel, "Edit", "Setting...")
