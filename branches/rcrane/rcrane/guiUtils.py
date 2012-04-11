#!/usr/bin/env python
"""Common constants and functions used in the various GUI modules"""

# Copyright 2010, 2011 Kevin Keating
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

import os.path
import gtk

#spacing parameters for GTK boxes
VBOX_SPACING = 5
HBOX_SPACING = 10
BUTTON_SPACING = 2

rcranePath = os.path.dirname(os.path.abspath(__file__))
    
def buttonWithIcon(label, image):
    """Create a button object with the specified label and GTK stock image.
    
    ARGUMENTS:
        label  - the text label for the button
        image  - the name of the GTK stock image
    RETURNS:
        button - the constructed button object
    NOTE:
        The buttons created here are different than the stock buttons created by gtk.Button().
        Buttons created by gtk.Button() ignore the label if it does not match the stock image and
        resize differently than buttons created here.
    """
    
    button = gtk.Button()
    hbox = gtk.HBox(False, 0)
    #hbox.set_border_width(2)
    buttonImage = gtk.image_new_from_stock(image, gtk.ICON_SIZE_BUTTON)
    buttonLabel = gtk.Label(label)
    hbox.pack_start(buttonImage, False, False, 3)
    hbox.pack_start(buttonLabel, False, False, 3)
    button.add(hbox)
    
    return button


def createRCraneWindowObject(title = "RCrane"):
    """Create a new window object with the appropriate name and icon
    
    OPTIONAL ARGUMENTS:
        title - the title for the window.  Defaults to "RCrane"
    RETURNS:
        window - a window object
    """
    
    window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    window.set_title(title)
    window.set_icon_from_file(os.path.join(rcranePath, "icon.svg"))
    return window


def createRCraneErrorDialog(errorText):
    """Create a RCrane error dialog popup
    
    ARGUMENTS:
        errorText - the text of the error message to report to the user
    RETURNS:
        None
    """
    
    errorDialog = gtk.MessageDialog(type = gtk.MESSAGE_ERROR, buttons = gtk.BUTTONS_OK)
    errorDialog.set_title("RCrane Error")
    errorDialog.set_markup(errorText)
    errorDialog.run()
    errorDialog.destroy()