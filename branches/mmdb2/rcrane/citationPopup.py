#!/usr/bin/env python
"""Display a pop-up that informs the user about RCrane's citations"""

# Copyright 2011, 2012 Kevin Keating
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

from settings import addSetting

__dontShowAgain = False
__dontShowAgainButton = None

def createCitationPopup():
    """Creates and displays a popup informing the user about RCrane's citations.
    
    ARGUMENTS:
        None
    RETURNS:
        None
    EFFECTS:
        creates and displays the citation dialog
    NOTE:
        This function returns either:
            a) immediately, if the user has already checked the "Don't show this again" option
            or b) after the user has closed the dialog
    """
    
    global __dontShowAgainButton
    if __dontShowAgain:
        return
    
    citeDialog = gtk.Dialog("RCrane", None, gtk.DIALOG_MODAL | gtk.DIALOG_NO_SEPARATOR)
    citeDialog.set_resizable(False)
    citeDialog.connect("delete_event", __citationPopupCloseWin)
    
    citeLabel = gtk.Label()
    citeLabel.set_markup("  RCrane partially automates the RNA model building process.  All  \n" +
                         "  publications resulting from the use of RCrane must acknowledge:  \n\n" + 
                         "  Keating KS and Pyle AM. Semiautomated model building for RNA  \n" +
                         "  \tcrystallography using a directed rotameric approach. <span style='italic'>Proc  \n" +
                         "  \tNatl Acad USA</span>, <span weight='bold'>107</span> (2010) 8177-8182.  ")
    citeLabel.show()
    citeDialog.vbox.pack_start(citeLabel, True, True, 7)
    
    #change the layout so the check button appears on the left
    citeDialog.action_area.set_layout(gtk.BUTTONBOX_EDGE)
    
    dontShowAgainButton = gtk.CheckButton("Don't show this again")
    dontShowAgainButton.show()
    citeDialog.action_area.pack_start(dontShowAgainButton, False, False, 0)
    
    #store the Don't Show Again button object so we can later check to see if it's been checked
    __dontShowAgainButton = dontShowAgainButton
    
    okButton = gtk.Button("_OK", gtk.STOCK_OK)
    okButton.connect("clicked", __citationPopupOK, citeDialog)
    okButton.show_all()
    citeDialog.action_area.pack_end(okButton, False, False, 0)
    
    citeDialog.run()
    citeDialog.destroy()


def __citationPopupCloseWin(window, widget):
    """Respond to the user closing the citation pop-up window
    
    ARUGMENTS:
        window - the pop-up window object
        widget - the widget used to close the window
    NOTES:
        This function simply calls __citationPopupOK with a reversed argument list
    """
    
    __citationPopupOK(widget, window)


def __citationPopupOK(widget, window):
    """Respond to the user clicking OK
    
    ARGUMENTS:
        widget - the OK button
        window - the pop-up window object
    RETURNS:
        None
    """
    
    global __dontShowAgainButton
    
    #if the user has clicked Don't show again
    if __dontShowAgainButton.get_active():
        #then store that setting in a module variable
        global __dontShowAgain
        __dontShowAgain = True
        
        #and store the setting in the RCrane settings file
        addSetting("rcrane.dontShowCitationPopup()")
            #if RCrane isn't installed, this command will produce an error the next time Coot is launched
            #but the error shouldn't cause any problems (other than possibly confusing the user)
    
    #clear the Don't show again button object
    __dontShowAgainButton = None
    
    #close the window
    window.destroy()


def dontShowPopup():
    """Set the flag so that the citation pop-up is not shown
    
    ARGUMENTS:
        None
    RETURNS:
        None
    NOTES:
        This function is intended to be called in the RCrane settings file
    """
    
    global __dontShowAgain
    __dontShowAgain = True
