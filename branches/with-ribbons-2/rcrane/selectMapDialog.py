#!/usr/bin/env python
"""Create a select map dialog box, similar to Coot's show_select_map_dialog(), but allow for a callback."""

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

from coot import molecule_name, graphics_n_molecules, is_valid_map_molecule, imol_refinement_map, set_imol_refinement_map
import gtk


def selectMapDialog(callback = None):
    """Display a GUI for selecting the refinement map.
    
    OPTIONAL ARGUMENTS:
        callback - a function to call when the user clicks OK
    RETURNS:
        None
    """
    
    #create a list of map molecules and their names
    mapList = []
    for imol in xrange(graphics_n_molecules()):
        if is_valid_map_molecule(imol):
            mapList.append(imol)
    
    #if we have no maps, then display an error and don't invoke the callback
    if len(mapList) == 0:
        warningMessage = gtk.MessageDialog(type=gtk.MESSAGE_ERROR, buttons=gtk.BUTTONS_OK, message_format="Please load a map before using RCrane.")
        warningMessage.set_title("RCrane Error")
        warningMessage.run()
        warningMessage.destroy()
    
    else:
        #if we have one or more maps
        dialog = gtk.Dialog(title="Select Map for Fitting", buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK))
        
        #create a combo box that lists maps
        mapListStore = gtk.ListStore(int, str)
        for imol in mapList:
            mapListStore.append([imol, str(imol) + " " + molecule_name(imol)])
        mapCombo = gtk.ComboBox(mapListStore)
        mapRender = gtk.CellRendererText()
        mapCombo.pack_start(mapRender, True)
        mapCombo.add_attribute(mapRender, 'text', 1)
        
        #select the current refinemnt map, or the first map if no map has been selected
        currentMap = imol_refinement_map()
        if currentMap == -1:
            set_imol_refinement_map(mapListStore[0][0])
            mapCombo.set_active(0)
        else:
            #figure out which combobox row corresponds to the currently selected map
            comboIndex = [i[0] for i in mapListStore].index(currentMap)
            mapCombo.set_active(comboIndex)
        
        #connect the changed signal for the combobox
        mapCombo.connect("changed", __changeMap)
        
        #add the label and frame
        chooseLabel = gtk.Label("Choose a Map:")
        frame = gtk.Frame()
        frameBox = gtk.VBox(False, 0)
        frame.add(frameBox)
        frameBox.pack_start(chooseLabel, False, False, 2)
        frameBox.pack_start(mapCombo, False, False, 2)
        
        #pack things into the top area of the dialog box
        dialog.vbox.pack_start(frame, False, False, 2)
        dialog.vbox.show_all()
        
        #connect the callback function
        if callback is not None:
            #if there is a callback function and the user clicks OK, then call the callback
            #if the user clicks the close button, then don't call the callback
            def checkResponse(dialog, responseID, callback):
                if responseID == gtk.RESPONSE_OK: callback()
            dialog.connect("response", checkResponse, callback)
        
        #create the dialog box
        dialog.run()
        dialog.destroy()
        

def __changeMap(mapCombo):
    """Respond to the user selecting a map
        
    ARGUMENTS:
        mapCombo - the ComboBox object for map selection
    RETURNS:
        None
    EFFECTS:
        sets the refinement map to the newly selected map
    """
    
    #figure out the newly selected map
    mapListStore  = mapCombo.get_model()
    selectedMapIter = mapCombo.get_active_iter()
    selectedImol = mapListStore.get_value(selectedMapIter, 0)
    
    #set the refinement map
    set_imol_refinement_map(selectedImol)
