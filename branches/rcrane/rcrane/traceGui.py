#!/usr/bin/env python
"""A class for the graphical interface used to trace the RNA backbone."""

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

from coot import imol_refinement_map, new_generic_object_number, set_display_generic_object, close_generic_object, graphics_draw, set_rotation_centre, to_generic_object_add_line, rotation_centre_position, is_closed_generic_object_p, show_select_map_dialog, zoom_factor, screen_vectors_py, add_status_bar_text
screen_vectors = screen_vectors_py
def rotation_centre():   return [rotation_centre_position(0), rotation_centre_position(1), rotation_centre_position(2)]
    #rotation_centre is defined in coot_utils.py, but you have to exec coot_utils.py (rather than import it)
    #since we only need this one function, it's easier to redefine it here
import gtk
import os
from copy import deepcopy
#from pprint import pprint   #for debugging
#from cPickle import dump    #for debugging
#from time import sleep      #for debugging

from nextPhos import NextPhos
from strucCalc import plus, minus, rotateAtoms
from rotamerSeq import determineRotamerSeq, determinePucker
from calcCoords import calcCoords
from reviewSuitesGui import ReviewSuitesGui
from guiUtils import buttonWithIcon, HBOX_SPACING, VBOX_SPACING, BUTTON_SPACING, createRCraneWindowObject
from pseudoMolecule import PseudoMolecule
from selectMapDialog import selectMapDialog
from citationPopup import createCitationPopup

srcPath = os.path.dirname(os.path.abspath(__file__))

CROSS_LENGTH_PHOS = 0.3
CROSS_COLOR_PHOS  = "green"
CROSS_LENGTH_C1   = 0.2
CROSS_COLOR_C1    = "cyan"
CANDIDATE_PHOS_TITLE = "Candidate phosphate positions"
CANDIDATE_C1_TITLE = "Candidate C1' positions"
BOX_LENGTH = 0.3
BATON_COLOR = "orange"
BATON_TITLE = "Currently Selected Peak"
ADJUST_PHOS_TITLE = "Manually adjusted phosphate position"

class TraceGui:
    """A class for the graphical interface used to trace the RNA backbone."""
    
    def __init__(self, direction = 3, curCoords = None):
        """Initialize a TraceGui object and create a GUI for tracing a molecule        
        
        OPTIONAL ARGUMENTS:
            direction      - which direction to trace the chain: 3 implies 5'->3'
                                                                 5 implies 3'->5'
                             defaults to 3 (5'->3')
            curCoords      - the current (x,y,z) coordinates for where to start building
                             if not given, the screen center is used
        RETURNS:
            an initialized gui object
        EFFECTS:
            Creates and displays an RCrane dialog box
        """
        
        #create the citation pop-up if necessary
        createCitationPopup()
        
        if curCoords is None:
            curCoords = rotation_centre()

        #initialize the object
        self.__pseudoMolecule = PseudoMolecule()
        self.__direction = direction
        self.__nextPhos = NextPhos(phosDistFilename  = os.path.join(srcPath, "data/phosDistSmoothed.csv"),
                                   phosAngleFilename = os.path.join(srcPath, "data/phosAngleSmoothed.csv"),
                         sugarPhosSugarAngleFilename = os.path.join(srcPath, "data/sugarPhosSugarAngleSmoothed.csv"),
                                       basesFilename = os.path.join(srcPath, "data/bases.pdb"),
                                   pseudoChiFilename = os.path.join(srcPath, "data/pseudoChiSmoothed.csv"))
        
        self.__candidatePhosphates   = None
        self.__currentPhosphateIndex = 0
        
        self.__candidateC1s     = None
        self.__currentC1Index   = 0
        
        #a list of base locations in the format self.__candidateBases[phosphateIndex][C1Index] = [base type, base coords]
        #base type is stored for each coordinate because bases are mutated lazily
        #(i.e we only mutate the base that is currently displayed, so when we display a base, we may have to mutate it
        #if the base type in the self.__candidateBases array is different than self.__baseType)
        self.__candidateBases = None
        self.__baseType = "C"
        
        #Coot drawing objects
        self.__candidatePhosObject  = None
        self.__candidateC1Object    = None
        self.__batonObject          = None
        
        #coordinates that have been set by the user (via "Manually Adjust" or "Anti/Syn Flip" buttons)
        self.__customPhosphate = None
        #self.__customC1        = None
        self.__customBase      = None
        
        #drawing objects for manually adjusting coordinates
        self.__customPhosObject = None
        
        #coordinates before the user begins manually adjusting them
        self.__adjustPhosOrig = None
        self.__adjustPhosWindow = None
        
        self.__adjustBaseOrig = None
        self.__adjustBaseWindow = None
        
        self.__translatedPhosCoords = None
        self.__adjustedBaseCoords = None
        
        self.__previousPhosTranslateValue = [0,0,0]
        self.__previousBaseTranslateValue = [0,0,0]
        self.__previousBaseRotateValue = [0,0,0]
        self.__previousBaseRotateValueChi = 0
        
        self.__deleteHandlerID = None
        
        #create GUI for selecting the initial atom (using nextPhos.firstPhos)
        self.selectInitialPhos(curCoords)
    
    
    def selectInitialPhos(self, curCoords):
        """Create a GUI for building the initial phosphate
        
        ARGUMENTS:
           curCoords - the current (x,y,z) coordinates for where to start building
        RETURNS:
           none
        EFFECTS:
           creates a GUI for building the initial phosphate for a chain trace
        REQUIREMENTS:
           a refinement map must be selected (i.e. imol_refinement_map must return a valid map ID)
               otherwise this function will display a map selection dialog and then return
        """
        
        #get the map that Coot is currently using
        mapNum = imol_refinement_map()
        
        #if there's no selected map, then show the select map dialog and return
        if mapNum == -1:
            print "No refinement map set for RCrane building"
            #add_status_bar_text("No refinement map set for RCrane building")
            #show_select_map_dialog()
            selectMapDialog(lambda: TraceGui(curCoords))
            return
        
        #create a window for selecting peaks
        window = createRCraneWindowObject()
        self.__deleteHandlerID = window.connect("delete_event", self.__initialPhosCloseWin)
        
        #initially, simply use the window to inform the user that peak selection is currently being done (since it will probably take a while)
        waitLabel = gtk.Label("\n  Calculating potential phosphate coordinates...  \n\n  Please wait...  \n")
        waitLabel.set_justify(gtk.JUSTIFY_CENTER)
        window.add(waitLabel)
        window.show_all()
        
        #make sure the window gets drawn before we start searching for phosphates
        while gtk.events_pending():
            gtk.main_iteration(False)
        
        
        #find all nearby electron density peaks
        self.__candidatePhosphates = self.__nextPhos.firstPhos(mapNum, curCoords)
        #self.__candidatePhosphates.append(curCoords)
        
        #draw crosses for each candidate peak
        self.__drawPhosPeaks()
        
        #highlight the first peak
        self.__drawBox(self.__candidatePhosphates[0])
        
        #create all buttons and labels
        topLabel        = gtk.Label     ("  Select starting phophate location:  ")
        nextButton      = buttonWithIcon("  Next  ",        "gtk-media-next")
        previousButton  = buttonWithIcon("  Previous  ",    "gtk-media-previous")
        adjustButton    = gtk.Button    ("  Manually Adjust  ")
        acceptButton    = buttonWithIcon("  Accept  ",      "gtk-ok")
        cancelButton    = buttonWithIcon("  Cancel  ",      "gtk-cancel")
        
        #connect the buttons to the appropriate function
        nextButton.connect      ("clicked", self.__initialPhosNext)
        previousButton.connect  ("clicked", self.__initialPhosPrevious)
        adjustButton.connect    ("clicked", self.__createAdjustInitialPhosDialog, window)
        acceptButton.connect    ("clicked", self.__initialPhosAccept, window)
        cancelButton.connect    ("clicked", self.__initialPhosCancel, window)
        
        #pack the buttons and lables into the window
        windowBox = gtk.VBox(False, 5)
        window.remove(waitLabel)
        window.add(windowBox)
        window.resize(1,1) #remove any size constraints that were caused by the waiting message
        
        frame = gtk.Frame()
        windowBox.pack_start(frame, False, False, 2)
        
        frameBox = gtk.VBox(False, 0)
        frame.add(frameBox)
        
        frameBox.pack_start(topLabel, False, False, HBOX_SPACING)
        nextPrevBox = gtk.HBox(True, HBOX_SPACING)
        nextPrevBox.pack_start(previousButton, True, True, BUTTON_SPACING)
        nextPrevBox.pack_start(nextButton, True, True, BUTTON_SPACING)
        frameBox.pack_start(nextPrevBox, False, False, BUTTON_SPACING)
        
        adjustBox = gtk.HBox(False, 0)
        adjustBox.pack_start(adjustButton, True, False, BUTTON_SPACING)
        frameBox.pack_start(adjustBox, False, False, BUTTON_SPACING)
        
        separator = gtk.HSeparator()
        windowBox.pack_start(separator)
        
        acceptCancelBox = gtk.HBox(False, HBOX_SPACING)
        acceptCancelBox.pack_start(acceptButton, True, False, BUTTON_SPACING)
        acceptCancelBox.pack_start(cancelButton, True, False, BUTTON_SPACING)
        windowBox.pack_start(acceptCancelBox, False, False, BUTTON_SPACING)
        
        #make the accept button the default focus
        acceptButton.grab_focus()
        
        #display the new window
        window.show_all()
    
    
    def selectNextNt(self, prevPhos, window = None, initialNtCoords = None, initialNtType = None, initialPhos = None):
        """Create a GUI for building the next nucleotide
        
        ARGUMENTS:
           prevPhos - the (x,y,z) coordinates of the 5' phosphate
        OPTIONAL ARGUMENTS:
           window          - a pyGTK+ window object (from the previous iteration of selectNextNt, if there was one)
           initialNtCoords - a dictionary of atom coordinates indicating what the initially selected nucleotide should be
                             This is typically used if the user has just clicked "Previous Nt"
                             If not given, the the most likely predicted base location is used
           initialNtType   - what type of base (i.e. A, C, G, or T) should be initially selected
                             this is typically used if the user has just clicked "Previous Nt"
                             If not given, then self.__baseType is used.
                             If initialNtCoords is given, then this argument MUST match the base type described by initialNtCoords.
           initialPhos     - the coordinates for the initially selected phosphate should be
                             If this is given, then initialNtCoords and must *NOT* be given
        RETURNS:
           none
        EFFECTS:
           creates a GUI for building the next nucleotide of a chain trace
        REQUIREMENTS:
           a refinement map must be selected (i.e. imol_refinement_map must return a valid map ID)
        """
        
        #get the map that Coot is currently using
        mapNum = imol_refinement_map()
        
        #create a window if we don't already have one
        if window is None:
            window = createRCraneWindowObject()
        else:
            #if we do have one, them empty it
            windowChild = window.get_child()
            window.remove(windowChild)
            window.disconnect(self.__deleteHandlerID)
        
        #make sure that closing the window will invoke the proper function
        self.__deleteHandlerID = window.connect("delete_event", self.__nextNtCloseWin)
        
        #initially, simply use the window to inform the user that peak selection is currently being done (since it will probably take a while)
        waitLabel = gtk.Label("\n  Calculating potential phosphate coordinates...  \n\n  Please wait...  \n")
        waitLabel.set_justify(gtk.JUSTIFY_CENTER)
        window.add(waitLabel)
        window.show_all()
        
        #make sure the window gets drawn before we start searching for phosphates
        while gtk.events_pending():
            gtk.main_iteration(False)
        
        
        #find candidate phosphate positions
        if self.__pseudoMolecule.getNumNts() == 1:
            #if this is the second nucleotide, then we can't factor in phosphate angles
            (self.__candidatePhosphates, self.__candidateC1s) = self.__nextPhos.secondPhos(mapNum, prevPhos, self.__direction)
        elif self.__direction == 3:
            (self.__candidatePhosphates, self.__candidateC1s) = self.__nextPhos.nextPhos  (mapNum, prevPhos, self.__pseudoMolecule.getPhosCoords(-2), self.__pseudoMolecule.getSugarCoords(-2), self.__direction)
        else: #self.__direction == 5
            (self.__candidatePhosphates, self.__candidateC1s) = self.__nextPhos.nextPhos  (mapNum, prevPhos, self.__pseudoMolecule.getPhosCoords(2), self.__pseudoMolecule.getSugarCoords(0), self.__direction)
        
        #initialize the __currentC1Index array
        self.__currentC1Index = [0] * len(self.__candidatePhosphates)
        
        #initialize the __candidateBases array
        self.__candidateBases = []
        for curC1s in self.__candidateC1s:
            self.__candidateBases.append([None] * len(curC1s))
        
        #draw the candidate phosphate peaks
        self.__drawPhosPeaks()
        
        self.__customPhosphate = None
        self.__customBase = None
        if initialNtType is not None:
            self.__baseType = initialNtType
        
        self.__currentPhosphateIndex = 0
        
        if initialNtCoords is not None or initialPhos is not None:
            if initialPhos is not None:
                initialPhosCoords = initialPhos
            else:
                initialPhosCoords = initialNtCoords["P"]
            
            #determine if it was a custom phosphate position or not
            for (index, coords) in enumerate(self.__candidatePhosphates):
                if coords[0] == initialPhosCoords[0] and coords[1] == initialPhosCoords[1] and coords[2] == initialPhosCoords[2]:
                    #Note that here (and below in the C1' loop) we're comparing floating point numbers using == and expecting sane results.
                    #This is normally a terrible idea due to rounding errors; however, the numbers we're comparing here have all
                    #come out of the exact same deterministic Coot function, so any values that are supposed to be equal
                    #should be exactly equal, not just equal in their relevant significant figures.
                    self.__currentPhosphateIndex = index
                    #print "Using phosphate index " + str(index)
                    break
            else:
                #if we haven't found a matching phosphate, then the user manually set these phosphate coordinates
                self.__customPhosphate = initialPhosCoords
                
                #which means that we don't know which C1' coordinates were used
                self.__customBase = [self.__baseType, deepcopy(initialNtCoords)]
                del self.__customBase[1]["P"] #we don't want the phoshate coordinate in the __customBase dictionary
                
                #print "Using custom phosphate"
            
            
            #if this function was called with initialPhos (rather than initialNtCoords) then we don't have to worry about
            #initial base coordinates
            #otherwise, if the phosphate was a standard one (i.e. non-custom coordaintes), then we should 
            #determine if there was a standard or custom base position
            if initialPhos is None and self.__customPhosphate is None:
                initialC1Coords = initialNtCoords["C1'"]
                
                for (index, coords) in enumerate(self.__candidateC1s[self.__currentPhosphateIndex]):
                    if coords[0] == initialC1Coords[0] and coords[1] == initialC1Coords[1] and coords[2] == initialC1Coords[2]:
                        self.__currentC1Index[self.__currentPhosphateIndex] = index
                        #print "Using C1' index " + str(index)
                        break
                    #else:
                    #    print "C1' doesn't match:"
                    #    pprint(coords)
                else:
                    #if we haven't found a matching phosphate, then the user manually set the base coordinates
                    self.__customBase = [self.__baseType, deepcopy(initialNtCoords)]
                    del self.__customBase[1]["P"] #we don't want the phoshate coordinate in the __customBase dictionary
                    #print "Using custom C1'"
                    #pprint(self.__customBase)
                
                #draw the C1' peaks
                self.__drawC1Peaks()
        else:
            #draw the C1' peaks
            self.__drawC1Peaks()
        
        
        #these variables are just to save typing
        phosIndex = self.__currentPhosphateIndex
        c1index = self.__currentC1Index[self.__currentPhosphateIndex]
        
        #calculate a base location for the current phosphate and C1' combination
        self.__candidateBases[phosIndex][c1index] = self.__nextPhos.findBase(mapNum, self.__candidateC1s[phosIndex][c1index], prevPhos, self.__candidatePhosphates[phosIndex], self.__baseType, self.__direction)
        
        #if we've been given an initialNtCoords, find out if the user set a custom base position but didn't move C1'
        if initialNtCoords is not None and self.__customBase is None:
            for (atomName, coords) in self.__candidateBases[phosIndex][c1index][1].iteritems():
                if initialNtCoords[atomName][0] != coords[0] or initialNtCoords[atomName][1] != coords[1] or initialNtCoords[atomName][2] != coords[2]:
                    #there is a custom base position
                    self.__customBase = [self.__baseType, deepcopy(initialNtCoords)]
                    del self.__customBase[1]["P"] #we don't want the phoshate coordinate in the __customBase dictionary
                    #print "Using custom base because of atom " + atomName
                    break
        
        
        #determine what coordinates we should draw
        curPhosCoords = None
        curC1Coords = None
        curBaseCoords = None
        
        #determine which phosphate coordinates to display initially
        if self.__customPhosphate is None:
            curPhosCoords = self.__candidatePhosphates[phosIndex]
        else:
            curPhosCoords = self.__customPhosphate
            #if we're starting with custom phosphate coordinates, then we need to draw a custom phosphate cross
            self.__drawAdjustPhosCross(curPhosCoords)
        
        #determine base and C1' coordinates
        if self.__customBase is None:
            curC1Coords = self.__candidateC1s[phosIndex][c1index]
            curBaseCoords = self.__candidateBases[phosIndex][c1index]
        else:
            curC1Coords = self.__customBase[1]["C1'"]
            curBaseCoords = self.__customBase
        
        #draw the selected candidate nucleotide
        self.__drawNt(prevPhos, curC1Coords, curPhosCoords, curBaseCoords)
        
        #create all buttons and labels
        topLabel            = gtk.Label     ("Trace the next nucleotide:")
        nextPhosButton      = buttonWithIcon("  Next Phos  ",        "gtk-media-next")
        prevPhosButton      = buttonWithIcon("  Previous Phos ",    "gtk-media-previous")
        adjustPhosButton    = gtk.Button    ("  Manually Adjust  ")
        
        nextC1Button      = buttonWithIcon(u"  Next C1\N{prime}  ",        "gtk-media-next")
        prevC1Button  = buttonWithIcon(u"  Previous C1\N{prime}  ",    "gtk-media-previous")
        adjustBaseButton    = gtk.Button    ("  Manually Adjust  ")
        
        baseSelect = gtk.combo_box_new_text()
        baseSelect.append_text("Adenine (A)")
        baseSelect.append_text("Cytosine (C)")
        baseSelect.append_text("Guanine (G)")
        baseSelect.append_text("Uracil (U)")
        
        #set the current base selection based on self.__baseType
        baseTypeNum = {"A":0, "C":1, "G":2, "U":3}[self.__baseType]
        baseSelect.set_active(baseTypeNum)
        
        flipBaseButton = gtk.Button("  Anti/Syn Flip  ")
        
        acceptButton        = buttonWithIcon("  Accept Nt  ",      "gtk-ok")
        prevNtButton        = buttonWithIcon("  Previous Nt  ",    "gtk-cancel")
        buildBackboneButton = buttonWithIcon("  Build Backbone  ", "gtk-apply")
        closeButton         = buttonWithIcon("  Delete Trace  ",   "gtk-close")
        
        #connect the buttons to the appropriate functions
        nextPhosButton.connect("clicked", self.__nextNtNextPhos)
        prevPhosButton.connect("clicked", self.__nextNtPrevPhos)
        adjustPhosButton.connect("clicked", self.__createAdjustPhosDialog, window)
        
        nextC1Button.connect("clicked", self.__nextNtNextC1)
        prevC1Button.connect("clicked", self.__nextNtPrevC1)
        adjustBaseButton.connect("clicked", self.__createAdjustBaseDialog, window)
        
        acceptButton.connect("clicked", self.__nextNtAccept, window)
        prevNtButton.connect("clicked", self.__nextNtPrevNt, window)
        buildBackboneButton.connect("clicked", self.__nextNtBuildBackbone, window)
        closeButton.connect("clicked", self.__nextNtClose, window)
        
        baseSelect.connect("changed", self.__nextNtBaseSelect)
        flipBaseButton.connect("clicked", self.__flipBase)
        
        
        if self.__direction == 3:
            directionString = "5' -> 3'"
        else:
            directionString = "3' -> 5'"
        
        numNts = self.__pseudoMolecule.getNumNts()
        if numNts == 1:
            directionLabel = gtk.Label("Tracing " + directionString)
            switchDirectionButton = buttonWithIcon("  Switch  ", "gtk-refresh")
            switchDirectionButton.connect("clicked", self.__switchDirection, prevPhos, window)
        else:
            directionLabel = gtk.Label("%i nucleotide%s traced (%s)" % ((numNts-1), "s" if numNts > 2 else "", directionString))
        
        
        #pack the buttons and label into the window
        windowBox = gtk.VBox(False, VBOX_SPACING)
        window.remove(waitLabel)
        window.add(windowBox)
        window.resize(1,1) #remove any size constraints that were caused by the waiting message
        
        topBox = gtk.VBox(False, 0)
        windowBox.pack_start(topBox, False, False, 6)
        topBox.pack_start(topLabel, False, False, 0)
        
        directionBox = gtk.HBox(False, HBOX_SPACING)
        topBox.pack_start(directionBox, False, False, VBOX_SPACING)
        
        if self.__pseudoMolecule.getNumNts() == 1:
            directionLabelAlign = gtk.Alignment(xalign=1, yalign=0.5)
            directionLabelAlign.add(directionLabel)
            
            switchDirectionAlign = gtk.Alignment(xalign=0, yalign=0.5)
            switchDirectionAlign.add(switchDirectionButton)
            
            directionBox.pack_start(directionLabelAlign, True, True, BUTTON_SPACING)
            directionBox.pack_start(switchDirectionAlign, True, True, BUTTON_SPACING)
        else:
            directionBox.pack_start(directionLabel, True, False, 0)
            #self.__directionBoxHeight was initially set when we created this window for tracing the first full nucleotide
            #i.e. when there was a switch direction button at the top of the window
            #We don't want to resize the window in between the first and second nucleotide, since that can
            #move buttons out from underneath the user's cursor.  To avoid this, we make sure that
            #directionBox is the same height regardless of whether or not the switch direction button is present.
            directionBox.set_size_request(-1, self.__directionBoxHeight)
        
        
        #the phosphate frame
        phosFrame = gtk.Frame("Phosphate")
        phosFrame.set_label_align(0.5, 0.5)
        windowBox.pack_start(phosFrame, False, False, 2)
        
        phosFrameBox = gtk.VBox(False, VBOX_SPACING)
        phosFrame.add(phosFrameBox)
        
        phosPrevNextBox = gtk.HBox(True, HBOX_SPACING)
        phosPrevNextBox.pack_start(prevPhosButton, True, True, BUTTON_SPACING)
        phosPrevNextBox.pack_start(nextPhosButton, True, True, BUTTON_SPACING)
        phosFrameBox.pack_start(phosPrevNextBox, False, False, BUTTON_SPACING)
        adjustPhosBox = gtk.HBox(False, 0)
        adjustPhosBox.pack_start(adjustPhosButton, True, False, BUTTON_SPACING)
        phosFrameBox.pack_start(adjustPhosBox, False, False, BUTTON_SPACING)
        
        #the C1' frame
        C1Frame = gtk.Frame("Base")
        C1Frame.set_label_align(0.5, 0.5)
        windowBox.pack_start(C1Frame, False, False, 2)
        
        C1FrameBox = gtk.VBox(False, VBOX_SPACING)
        C1Frame.add(C1FrameBox)
        
        C1PrevNextBox = gtk.HBox(True, HBOX_SPACING)
        C1PrevNextBox.pack_start(prevC1Button, True, True, BUTTON_SPACING)
        C1PrevNextBox.pack_start(nextC1Button, True, True, BUTTON_SPACING)
        C1FrameBox.pack_start(C1PrevNextBox, False, False, BUTTON_SPACING)
        baseBox = gtk.HBox(False, 0)
        baseBox.pack_start(baseSelect, True, True, BUTTON_SPACING)
        baseBox.pack_start(flipBaseButton, True, True, BUTTON_SPACING)
        C1FrameBox.pack_start(baseBox, False, False, BUTTON_SPACING)
        adjustC1Box = gtk.HBox(False, 0)
        adjustC1Box.pack_start(adjustBaseButton, True, False, BUTTON_SPACING)
        C1FrameBox.pack_start(adjustC1Box, False, False, BUTTON_SPACING)
        
        
        #the separator and accept, cancel, and close buttons
        separator = gtk.HSeparator()
        windowBox.pack_start(separator)
        
        acceptCancelBox = gtk.HBox(True, HBOX_SPACING)
        acceptCancelBox.pack_start(acceptButton, True, True, BUTTON_SPACING)
        acceptCancelBox.pack_start(prevNtButton, True, True, BUTTON_SPACING)
        windowBox.pack_start(acceptCancelBox, False, False, BUTTON_SPACING)
        buildCloseHbox = gtk.HBox(True, VBOX_SPACING)
        windowBox.pack_start(buildCloseHbox, False, False, BUTTON_SPACING)
        buildCloseBox = gtk.VBox(True, VBOX_SPACING)
        buildCloseHbox.pack_start(buildCloseBox, False, False, BUTTON_SPACING)
        buildCloseBox.pack_start(buildBackboneButton, True, False, BUTTON_SPACING)
        buildCloseBox.pack_start(closeButton, True, False, BUTTON_SPACING)
        
        acceptButton.grab_focus() #make the accept button the default focus
        
        #display the new window
        window.show_all()
        
        #record the height of the directionBox so that it always stays the same height
        #We don't want it to change between tracing the first nucleotide (when we have a switch
        #direction button) and the second nucleotide (when we don't have the button)
        self.__directionBoxHeight = directionBox.size_request()[1]
    
    
    def __mutateBase(self):
        """Mutate the currently selected base to self.__baseType
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            the current base is mutated to self.__baseType and redrawn
        """
        
        #if there is a Manually Adjust window open, click accept
        #Technically, we don't need to do this.  We could make it so the manual adjustment knows about the base mutation,
        #but just closing the window is far simpler, and should be good enough.
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosAccept(None, self.__adjustPhosWindow)
        
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseAccept(None, self.__adjustBaseWindow)
        
        #get the current base coordinates
        phosIndex = self.__currentPhosphateIndex
        c1Index = self.__currentC1Index[phosIndex]
        curBase = None
        c1loc = None
        if self.__customBase is not None:
            curBase = self.__customBase
            c1loc = self.__customBase[1]["C1'"]
        else:
            curBase = self.__candidateBases[phosIndex][c1Index]
            c1loc = self.__candidateC1s[phosIndex][c1Index]
        
        #if the base we're mutating to is the same as the base we are now, don't do anything
        if curBase[0] == self.__baseType:
            return
        
        newBase = self.__nextPhos.mutateBase(curBase, self.__baseType)
        if self.__customBase is not None:
            self.__customBase = newBase
        else:
            self.__candidateBases[phosIndex][c1Index] = newBase
        
        #determine the current phosphate location so we know where to redrawing it
        curPhosCoords = None
        if self.__customPhosphate is None:
            curPhosCoords = self.__candidatePhosphates[phosIndex]
        else:
            curPhosCoords = self.__customPhosphate
        
        self.__drawNt(self.__getPrevPhosCoords(), c1loc, curPhosCoords, newBase)
        
        
    def __flipBase(self, widget=None):
        """Perform an anti/syn rotation of the base.
        
        OPTIONAL ARGUMENTS:
            widget - the widget that called this function.  This is passed by pyGTK, but is ignored.
        RETURNS:
            None
        EFFECTS:
            Rotates the current base and redraws it
        NOTE:
            With pyrimidines, this function simply rotates the base by 180 about chi
            However, for purines, the base is moved slightly within the plane so that it occupies roughly the same
                location in space as it did before the rotation.
        """
        
        #if there is a Manually Adjust window open, click accept
        #Technically, we don't need to do this.  We could make it so the manual adjustment knows about the base flip,
        #but just closing the window is far simpler, and should be good enough.
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosAccept(None, self.__adjustPhosWindow)
        
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseAccept(None, self.__adjustBaseWindow)
        
        #get the current base coordinates
        phosIndex = self.__currentPhosphateIndex
        c1Index = self.__currentC1Index[phosIndex]
        
        c1loc = None
        if self.__customBase is not None:
            curBase = self.__customBase
            c1loc = self.__customBase[1]["C1'"]
        else:
            curBase = self.__candidateBases[phosIndex][c1Index]
            c1loc = self.__candidateC1s[phosIndex][c1Index]
        
        newBase = self.__nextPhos.flipBase(curBase)
        self.__customBase = newBase
        
        #determine the current phosphate location
        curPhosCoords = None
        if self.__customPhosphate is None:
            curPhosCoords = self.__candidatePhosphates[phosIndex]
        else:
            curPhosCoords = self.__customPhosphate
        
        self.__drawNt(self.__getPrevPhosCoords(), c1loc, curPhosCoords, newBase)
    
    
    ##########################################################################################
    #   Functions for the adjust phosphate dialog
    ##########################################################################################
    
    def __createAdjustPhosDialog(self, widget, parentWindow, sliderFunc = None, acceptFunc = None, cancelFunc = None, closeWinFunc = None):
        """Create a dialog box that allows the user to manually change the phosphate position.
        ARGUMENTS:
            widget - the widget that called this function.  Ignored.
            parentWindow - the RCrane window.  The newly created dialog box is set as a transient for parentWindow
        OPTIONAL ARGUMENTS:
            By default, this function assumes that it is being used to move a next phosphate, not an initial phosphate
            By provided all these arguments, a dialog appropriate for moving the initial phosphate can be created
            sliderFunc   - the function to call when moving the sliders
            acceptFunc   - the function to call when the user clicks OK
            cancelFunc   - the function to call when the user clicks Cancel
            closeWinFunc - the function to call when the user closes the window
        RETURNS:
            None
        EFFECTS:
            Creates and displays a dialog box for adjusting the phosphate location
        """
        
        #if there is still a Manually Adjust Base window open, click accept
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseAccept(None, self.__adjustBaseWindow)
        
        #by default, use fuctions for the nextNt dialog (not the initialNt)
        if sliderFunc is None:
            sliderFunc = self.__adjustPhos
        if acceptFunc is None:
            acceptFunc = self.__adjustPhosAccept
        if cancelFunc is None:
            cancelFunc = self.__adjustPhosCancel
        if closeWinFunc is None:
            closeWinFunc = self.__adjustPhosCloseWin
        
        #create the adjustment dialog
        dialogWindow = createRCraneWindowObject("Adjust Phosphate")
        dialogWindow.set_transient_for(parentWindow)
        
        dialogWindow.connect("delete_event", closeWinFunc)
        
        windowBox = gtk.VBox(False, VBOX_SPACING)
        dialogWindow.add(windowBox)
        
        xFrame = gtk.Frame("X translation")
        windowBox.pack_start(xFrame, False, False, 2)
        
        xAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        xScale = gtk.HScale(xAdjustment)
        xScale.set_draw_value(False)
        xFrame.add(xScale)
        
        
        yFrame = gtk.Frame("Y translation")
        windowBox.pack_start(yFrame, False, False, 2)
        
        yAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        yScale = gtk.HScale(yAdjustment)
        yScale.set_draw_value(False)
        yFrame.add(yScale)
        
        
        zFrame = gtk.Frame("Z translation")
        windowBox.pack_start(zFrame, False, False, 2)
        
        zAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        zScale = gtk.HScale(zAdjustment)
        zScale.set_draw_value(False)
        zFrame.add(zScale)
        
        
        separator = gtk.HSeparator()
        windowBox.pack_start(separator)
        
        okButton        = buttonWithIcon("  OK  ",      "gtk-ok")
        cancelButton    = buttonWithIcon("  Cancel  ",  "gtk-cancel")
        
        okCancelBox = gtk.HBox(True, HBOX_SPACING)
        okCancelBox.pack_start(okButton, True, True, BUTTON_SPACING)
        okCancelBox.pack_start(cancelButton, True, True, BUTTON_SPACING)
        windowBox.pack_start(okCancelBox, False, False, BUTTON_SPACING)
        
        #connect the sliders to the __adjustPhos function
        xAdjustment.connect("value_changed", sliderFunc, 0)
        yAdjustment.connect("value_changed", sliderFunc, 1)
        zAdjustment.connect("value_changed", sliderFunc, 2)
        
        okButton.connect("clicked", acceptFunc, dialogWindow)
        cancelButton.connect("clicked", cancelFunc, dialogWindow)
        
        #remember the original phopshate coordinates
        if self.__customPhosphate is not None:
            self.__adjustPhosOrig = deepcopy(self.__customPhosphate)
        else:
            self.__adjustPhosOrig = self.__candidatePhosphates[self.__currentPhosphateIndex]
        
        self.__translatedPhosCoords = deepcopy(self.__adjustPhosOrig)
        
        #reset the previous slider values
        self.__previousPhosTranslateValue = [0,0,0]
        
        #store the window so that we can destroy it if the user starts clicking buttons in the main RCrane window
        self.__adjustPhosWindow = dialogWindow
        
        #display the dialog
        dialogWindow.show_all()
    
    
    def __adjustPhos(self, adjustment, dim):
        """Move a next (i.e. not initial) phosphate atom in response to the user moving a slider in the adjust phosphate dialog
        
        ARGUMENTS:
            adjustment - the adjustment object
            dim        - which dimension to move the phosphate in
                         0 for the screen x-axis, 1 for the screen y-axis, and 2 for the screen z-axis
        RETURNS:
            None
        EFFECTS:
            Moves and redraws the phosphate
        """
        
        #figure out how much we moved since the last time this function was called
        adjDiff = adjustment.value - self.__previousPhosTranslateValue[dim]
        self.__previousPhosTranslateValue[dim] = adjustment.value #store the new value
        
        screenVectors = screen_vectors()
        
        #fetch and store the zoom factor, so we don't have to call the function three times
        #(since each call causes the creation and destruction of a C++ object, so it's faster to cache the result)
        zoom = zoom_factor()
        
        #determine how much to move the phosphate by
        #using the same formula as Coot's Rotate/Translate Zone function
        #(including multiplying y translations by -1)
        movement = [None, None, None]
        scaleFactor = 0.002
        if dim == 1:
            scaleFactor *= -1
        movement[0] = screenVectors[dim][0] * adjDiff * scaleFactor * zoom
        movement[1] = screenVectors[dim][1] * adjDiff * scaleFactor * zoom
        movement[2] = screenVectors[dim][2] * adjDiff * scaleFactor * zoom
        
        #calculate the new phosphate coordinates
        self.__translatedPhosCoords = plus(self.__translatedPhosCoords, movement)
        
        phosIndex = self.__currentPhosphateIndex
        c1Index = self.__currentC1Index[phosIndex]
        
        self.__drawAdjustPhosCross(self.__translatedPhosCoords)
        self.__drawNt(self.__getPrevPhosCoords(), self.__candidateC1s[phosIndex][c1Index], self.__translatedPhosCoords, self.__candidateBases[phosIndex][c1Index])
    
    
    def __adjustPhosCloseWin(self, window, widget):
        """Respond to the user closing the adjust phosphate dialog box
        
        ARUGMENTS:
            window - the window to close
            widget - the widget used to close the window
        NOTES:
            This function simply calls __adjustPhosCancel with a reversed argument list
        """
        self.__adjustPhosCancel(widget, window)
    
    
    def __adjustPhosCancel(self, widget, window):
        """Respond to the user clicking Cancel in the adjust phosphate dialog box
        
        ARUGMENTS:
            widget - the cancel button
            window - the dialog box window
        RETURNS:
            None
        EFFECTS:
            Closes the adjust phosphate dialog box and reverts the phosphate to its original coordinates
        """
        
        #destroy the AdjustPhosCross
        #(if we haven't moved the phosphate at all, then there's nothing to destroy)
        if self.__customPhosObject is not None:
            close_generic_object(self.__customPhosObject)
        
        #redraw the original version of the nucleotide
        phosIndex = self.__currentPhosphateIndex
        c1Index = self.__currentC1Index[phosIndex]
        
        origPhosCoords = None
        if self.__customPhosphate is not None:
            origPhosCoords = self.__customPhosphate
            self.__drawAdjustPhosCross(self.__customPhosphate)
        else:
            origPhosCoords = self.__candidatePhosphates[phosIndex]
        
        self.__drawNt(self.__getPrevPhosCoords(), self.__candidateC1s[phosIndex][c1Index], origPhosCoords, self.__candidateBases[phosIndex][c1Index])
        #self.__drawNt calls graphics_draw(), so we don't need to call that explicitly after closing self.__customPhosObject
        
        window.destroy()
        self.__adjustPhosWindow = None
        #if we don't explicitely do this, the __adjustPhosWindow stays as a GTK window object
        #which leads to problems in __nextNtAccept
    
    
    def __adjustPhosAccept(self, widget, window):
        """Respond to the user clicking OK in the adjust phosphate dialog box
        
        ARUGMENTS:
            widget - the OK button
            window - the dialog box window
        RETURNS:
            None
        EFFECTS:
            Closes the adjust phosphate dialog box and stores the coordinates that the user has selected
        """
        
        #store the current coordinates
        self.__customPhosphate = self.__translatedPhosCoords
        
        #destroy the C1' crosses
        close_generic_object(self.__candidateC1Object)
        graphics_draw()
        
        window.destroy()
        self.__adjustPhosWindow = None
        #if we don't explicitely do this, the __adjustPhosWindow stays as a GTK window object
        #which leads to problems in __nextNtAccept
    
    
    ##########################################################################################
    #   Functions for the adjust base dialog
    ##########################################################################################
    
    def __createAdjustBaseDialog(self, widget, parentWindow):
        """Creates a dialog that allows the user to manually change the base position.
        
        ARGUMENTS:
            widget 
        ARGUMENTS:
            widget - the widget that called this function.  Ignored.
            parentWindow - the RCrane window.  The newly created dialog box is set as a transient for parentWindow
        RETURNS:
            None
        EFFECTS:
            Creates and displays a dialog box for adjusting the base location
        """
        
        #if there is still a Manually Adjust Phosphate window open, click accept
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosAccept(None, self.__adjustPhosWindow)
        
        #create the adjustment dialog
        dialogWindow = createRCraneWindowObject("Adjust Base")
        dialogWindow.set_transient_for(parentWindow)
        
        dialogWindow.connect("delete_event", self.__adjustBaseCancel)
        
        windowBox = gtk.VBox(False, VBOX_SPACING)
        dialogWindow.add(windowBox)
        
        xFrame = gtk.Frame("X translation")
        windowBox.pack_start(xFrame, False, False, 2)
        xAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        xScale = gtk.HScale(xAdjustment)
        xScale.set_draw_value(False)
        xFrame.add(xScale)
        
        
        yFrame = gtk.Frame("Y translation")
        windowBox.pack_start(yFrame, False, False, 2)
        yAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        yScale = gtk.HScale(yAdjustment)
        yScale.set_draw_value(False)
        yFrame.add(yScale)
        
        
        zFrame = gtk.Frame("Z translation")
        windowBox.pack_start(zFrame, False, False, 2)
        zAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        zScale = gtk.HScale(zAdjustment)
        zScale.set_draw_value(False)
        zFrame.add(zScale)
        
        xRotFrame = gtk.Frame("X Rotation")
        windowBox.pack_start(xRotFrame, False, False, 2)
        xRotAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        xRotScale = gtk.HScale(xRotAdjustment)
        xRotScale.set_draw_value(False)
        xRotFrame.add(xRotScale)
        
        yRotFrame = gtk.Frame("Y Rotation")
        windowBox.pack_start(yRotFrame, False, False, 2)
        yRotAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        yRotScale = gtk.HScale(yRotAdjustment)
        yRotScale.set_draw_value(False)
        yRotFrame.add(yRotScale)
        
        zRotFrame = gtk.Frame("Z Rotation")
        windowBox.pack_start(zRotFrame, False, False, 2)
        zRotAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        zRotScale = gtk.HScale(zRotAdjustment)
        zRotScale.set_draw_value(False)
        zRotFrame.add(zRotScale)
        
        chiRotFrame = gtk.Frame("Chi Rotation")
        #chiRotFrame = gtk.Frame(u"\N{greek small letter chi} Rotation") #this works, but the chi looks too much like an x
        windowBox.pack_start(chiRotFrame, False, False, 2)
        chiRotAdjustment = gtk.Adjustment(0, -180, 360, 0.1, 1.0)
        chiRotScale = gtk.HScale(chiRotAdjustment)
        chiRotScale.set_draw_value(False)
        chiRotFrame.add(chiRotScale)
        
        separator = gtk.HSeparator()
        windowBox.pack_start(separator)
        
        okButton        = buttonWithIcon("  OK  ",      "gtk-ok")
        cancelButton    = buttonWithIcon("  Cancel  ",  "gtk-cancel")
        
        okCancelBox = gtk.HBox(True, HBOX_SPACING)
        okCancelBox.pack_start(okButton, True, True, BUTTON_SPACING)
        okCancelBox.pack_start(cancelButton, True, True, BUTTON_SPACING)
        windowBox.pack_start(okCancelBox, False, False, BUTTON_SPACING)
        
        #connect the sliders to the __adjustPhos function
        xAdjustment.connect("value_changed", self.__translateBase, 0)
        yAdjustment.connect("value_changed", self.__translateBase, 1)
        zAdjustment.connect("value_changed", self.__translateBase, 2)
        
        xRotAdjustment.connect("value_changed", self.__rotateBase, 0)
        yRotAdjustment.connect("value_changed", self.__rotateBase, 1)
        zRotAdjustment.connect("value_changed", self.__rotateBase, 2)
        chiRotAdjustment.connect("value_changed", self.__rotateBaseChi)
        
        okButton.connect("clicked", self.__adjustBaseAccept, dialogWindow)
        cancelButton.connect("clicked", self.__adjustBaseCancel, dialogWindow)
        
        #remember the original phopshate coordinates
        if self.__customBase is not None:
            self.__adjustBaseOrig = deepcopy(self.__customBase)
        else:
            phosIndex = self.__currentPhosphateIndex
            c1Index = self.__currentC1Index[phosIndex]
            self.__adjustBaseOrig = self.__candidateBases[phosIndex][c1Index]
        
        self.__adjustedBaseCoords = deepcopy(self.__adjustBaseOrig)
        
        #reset the previous slider values
        self.__previousBaseTranslateValue = [0,0,0]
        self.__previousBaseRotateValue = [0,0,0]
        self.__previousBaseRotateValueChi = 0
        
        #store the window so that we can destroy it if the user starts clicking buttons in the main RCrane window
        self.__adjustBaseWindow = dialogWindow
        
        #display the dialog
        dialogWindow.show_all()
    
    
    def __translateBase(self, adjustment, dim):
        """Move a base in response to the user moving a slider in the adjust base dialog.
        
        ARGUMENTS:
            adjustment - the adjustment object
            dim        - which dimension to move the base in
                         0 for the screen x-axis, 1 for the screen y-axis, and 2 for the screen z-axis
        RETURNS:
            None
        EFFECTS:
            Moves and redraws the base
        """
        
        #figure out how much we moved since the last time this function was called
        adjDiff = adjustment.value - self.__previousBaseTranslateValue[dim]
        self.__previousBaseTranslateValue[dim] = adjustment.value #store the new value
        
        #create a dummy screenVector object
        #screenVectors = [[1,0,0], [0,1,0], [0,0,1]]
        screenVectors = screen_vectors()
        
        #fetch and store the zoom factor, so we don't have to call the function three times
        #(since each call causes the creation and distruction of a C++ object, so it's faster to cache the result)
        zoom = zoom_factor()
        
        #determine how much to move the base by
        #using the same formula as Coot's Rotate/Translate Zone function
        #(including multiplying y translations by -1)
        movement = [None, None, None]
        scaleFactor = 0.002
        if dim == 1:
            scaleFactor *= -1
        movement[0] = screenVectors[dim][0] * adjDiff * scaleFactor * zoom
        movement[1] = screenVectors[dim][1] * adjDiff * scaleFactor * zoom
        movement[2] = screenVectors[dim][2] * adjDiff * scaleFactor * zoom
        
        #calculate the new base coordinates
        for curAtom, prevCoords in self.__adjustedBaseCoords[1].iteritems():
            self.__adjustedBaseCoords[1][curAtom] = plus(prevCoords, movement)
        
        #figure out the current phosphate coordinates so we know what to draw
        phosCoords = None
        if self.__customPhosphate is None:
            phosCoords  = self.__candidatePhosphates[self.__currentPhosphateIndex]
        else:
            phosCoords  = self.__customPhosphate
        
        self.__drawNt(self.__getPrevPhosCoords(), self.__adjustedBaseCoords[1]["C1'"], phosCoords, self.__adjustedBaseCoords)
    
    def __rotateBase(self, adjustment, dim):
        """Rotates a base in response to the user moving a slider in the adjust base dialog.
        
        ARGUMENTS:
            adjustment - the adjustment object
            dim        - which dimension to rotate the base about
                         0 for the screen x-axis, 1 for the screen y-axis, and 2 for the screen z-axis
        RETURNS:
            None
        EFFECTS:
            Rotates and redraws the base
        """
        
        #figure out how much we moved since the last time this function was called
        adjDiff = adjustment.value - self.__previousBaseRotateValue[dim]
        self.__previousBaseRotateValue[dim] = adjustment.value #store the new value
        
        #calculate the rotation axis
        screenVectors = screen_vectors()
        axis = screenVectors[dim] #the screen_vectors function returns unit vector, so we don't need to worry about that
        
        #translate the base to the origin
        c1loc = self.__adjustedBaseCoords[1]["C1'"]
        translatedCoords = dict()
        for curAtom, prevCoords in self.__adjustedBaseCoords[1].iteritems():
            translatedCoords[curAtom] = minus(prevCoords, c1loc)
        
        #do the rotation
        self.__adjustedBaseCoords[1] = rotateAtoms(translatedCoords, axis, adjDiff, c1loc)
        
        #figure out the current phosphate coordinates so we know what to draw
        phosCoords = None
        if self.__customPhosphate is None:
            phosCoords  = self.__candidatePhosphates[self.__currentPhosphateIndex]
        else:
            phosCoords  = self.__customPhosphate
        
        self.__drawNt(self.__getPrevPhosCoords(), c1loc, phosCoords, self.__adjustedBaseCoords)
    
    def __rotateBaseChi(self, adjustment):
        """Rotates a base about chi in response to the user moving a slider in the adjust base dialog.
        
        ARGUMENTS:
            adjustment - the adjustment object
        RETURNS:
            None
        EFFECTS:
            Rotates and redraws the base
        """
        
        #figure out how much we moved since the last time this function was called
        adjDiff = adjustment.value - self.__previousBaseRotateValueChi
        self.__previousBaseRotateValueChi = adjustment.value #store the new value
        
        #translate the base to the origin
        c1loc = self.__adjustedBaseCoords[1]["C1'"]
        translatedCoords = dict()
        for curAtom, prevCoords in self.__adjustedBaseCoords[1].iteritems():
            translatedCoords[curAtom] = minus(prevCoords, c1loc)
        
        axis = None
        if translatedCoords.has_key("N9"):
            axis = translatedCoords["N9"]
        else:
            axis = translatedCoords["N1"]
        
        #do the rotation
        self.__adjustedBaseCoords[1] = rotateAtoms(translatedCoords, axis, adjDiff, c1loc)
        
        #figure out the current phosphate coordinates so we know what to draw
        phosCoords = None
        if self.__customPhosphate is None:
            phosCoords  = self.__candidatePhosphates[self.__currentPhosphateIndex]
        else:
            phosCoords  = self.__customPhosphate
        
        self.__drawNt(self.__getPrevPhosCoords(), c1loc, phosCoords, self.__adjustedBaseCoords)
    
    def __adjustBaseCloseWin(self, window, widget):
        """Respond to the user closing the adjust base dialog box
        
        ARUGMENTS:
            window - the window to close
            widget - the widget used to close the window
        NOTES:
            This function simply calls __adjustBaseCancel with a reversed argument list
        """
        
        self.__adjustBaseCancel(widget, window)
    
    def __adjustBaseCancel(self, widget, window):
        """Respond to the user clicking Cancel in the adjust base dialog box
        
        ARUGMENTS:
            widget - the cancel button
            window - the dialog box window
        RETURNS:
            None
        EFFECTS:
            Closes the adjust base dialog box and reverts the base to its original coordinates
        """
        
        #redraw the original version of the nucleotide
        phosIndex = self.__currentPhosphateIndex
        c1Index = self.__currentC1Index[phosIndex]
        
        #figure out the current phosphate coordinates so we know what to draw
        phosCoords = None
        if self.__customPhosphate is None:
            phosCoords  = self.__candidatePhosphates[self.__currentPhosphateIndex]
        else:
            phosCoords  = self.__customPhosphate
        
        self.__drawNt(self.__getPrevPhosCoords(), self.__adjustBaseOrig[1]["C1'"], phosCoords, self.__adjustBaseOrig)
        #self.__drawNt calls graphics_draw(), so we don't need to call that explicitly after closing self.__customPhosObject
        
        window.destroy()
        self.__adjustBaseWindow = None
        #if we don't explicitely do this, the __adjustBaseWindow stays as a GTK window object
        #which leads to problems in __nextNtAccept
    
    def __adjustBaseAccept(self, widget, window):
        """Respond to the user clicking OK in the adjust base dialog box
        
        ARUGMENTS:
            widget - the OK button
            window - the dialog box window
        RETURNS:
            None
        EFFECTS:
            Closes the adjust base dialog box and stores the coordinates that the user has selected
        """
        
        #store the current coordinates
        self.__customBase = self.__adjustedBaseCoords
        
        window.destroy()
        self.__adjustBaseWindow = None
        #if we don't explicitely do this, the __adjustBaseWindow stays as a GTK window object
        #which leads to problems in __nextNtAccept
    
    
    ##########################################################################################
    #   Functions for buttons in the initial phosphate GUI (generated by selectInitialPhos)
    ##########################################################################################
    
    def __initialPhosNext(self, widget):
        """Select the next phosphate in the initial phosphate GUI.
        
        ARGUMENTS:
            widget - the button used to call this function.  Ignored.
        RETURNS:
            None
        EFFECTS:
            Highlights the next phosphate
        """
        
        #if we don't have a custom phophsate position, increment the phosphate index
        if self.__customPhosphate is None:
            self.__currentPhosphateIndex += 1
            if self.__currentPhosphateIndex >= len(self.__candidatePhosphates):
                self.__currentPhosphateIndex = 0
        else:
            self.__customPhosphate = None
            close_generic_object(self.__customPhosObject)
        
        #if there is still a Manually Adjust Phosphate window open, close it
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
        
        #forget about any coordinate adjustments the user has made via the "Manually Adjust" button
        self.__customPhosphate = None
        
        self.__drawBox(self.__candidatePhosphates[self.__currentPhosphateIndex])
    
    def __initialPhosPrevious(self, widget):
        """Select the previous phosphate in the initial phosphate GUI.
        
        ARGUMENTS:
            widget - the button used to call this function.  Ignored.
        RETURNS:
            None
        EFFECTS:
            Highlights the previous phosphate
        """
        
        #if we don't have a custom phophsate position, decrement the phosphate index
        if self.__customPhosphate is None:
            self.__currentPhosphateIndex -= 1
            if self.__currentPhosphateIndex < 0 :
                self.__currentPhosphateIndex = len(self.__candidatePhosphates) - 1
        else:
            self.__customPhosphate = None
            close_generic_object(self.__customPhosObject)
        
        #if there is still a Manually Adjust Phosphate window open, close it
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
        
        #forget about any coordinate adjustments the user has made via the "Manually Adjust" button
        self.__customPhosphate = None
        
        self.__drawBox(self.__candidatePhosphates[self.__currentPhosphateIndex])
    
    
    def __createAdjustInitialPhosDialog(self, widget, window):
        """Creates a dialog to adjusting the initial phosphate coordinates.
        
        ARGUMENTS:
            widget - the button used to call this function.  Ignored.
            window - the RCrane window
        RETURNS:
            None
        EFFECTS:
            Creates and displays the adjust phosphate dialog
        """
        
        self.__createAdjustPhosDialog(widget, window, self.__adjustInitialPhos, self.__adjustInitialPhosAccept, self.__adjustInitialPhosCancel, self.__adjustInitialPhosCloseWin)
    
    def __adjustInitialPhos(self, adjustment, dim):
        """Move an initial phosphate atom in response to the user moving a slider in the adjust phosphate dialog
        
        ARGUMENTS:
            adjustment - the adjustment object
            dim        - which dimension to move the phosphate in
                         0 for the screen x-axis, 1 for the screen y-axis, and 2 for the screen z-axis
        RETURNS:
            None
        EFFECTS:
            Moves and redraws the phosphate
        """
        
        #figure out how much we moved since the last time this function was called
        adjDiff = adjustment.value - self.__previousPhosTranslateValue[dim]
        self.__previousPhosTranslateValue[dim] = adjustment.value #store the new value
        
        screenVectors = screen_vectors()
        
        #fetch and store the zoom factor, so we don't have to call the function three times
        #(since each call causes the creation and distruction of a C++ object, so it's faster to cache the result)
        zoom = zoom_factor()
        
        #determine how much to move the phosphate by
        #using the same formula as Coot's Rotate/Translate Zone function
        #(including multiplying y translations by -1)
        movement = [None, None, None]
        scaleFactor = 0.002
        if dim == 1:
            scaleFactor *= -1
        movement[0] = screenVectors[dim][0] * adjDiff * scaleFactor * zoom
        movement[1] = screenVectors[dim][1] * adjDiff * scaleFactor * zoom
        movement[2] = screenVectors[dim][2] * adjDiff * scaleFactor * zoom
        
        #calculate the new phosphate coordinates
        self.__translatedPhosCoords = plus(self.__translatedPhosCoords, movement)
        
        self.__drawAdjustPhosCross(self.__translatedPhosCoords)
        self.__drawBox(self.__translatedPhosCoords)
    
    def __adjustInitialPhosCloseWin(self, window, widget):
        """Respond to the user closing the adjust initial phosphate dialog box
        
        ARUGMENTS:
            window - the window to close
            widget - the widget used to close the window
        NOTES:
            This function simply calls __adjustPhosCancel with a reversed argument list
        """
        
        self.__adjustInitialPhosCancel(widget, window)
    
    def __adjustInitialPhosCancel(self, widget, window):
        """Respond to the user clicking Cancel in the adjust initial phosphate dialog box
        
        ARUGMENTS:
            widget - the cancel button
            window - the dialog box window
        RETURNS:
            None
        EFFECTS:
            Closes the adjust phosphate dialog box and reverts the phosphate to its original coordinates
        """
        
        #destroy the AdjustPhosCross
        #(if we haven't moved the phosphate at all, then there's nothing to destroy)
        if self.__customPhosObject is not None:
            close_generic_object(self.__customPhosObject)
        
        #redraw the original version of the nucleotide
        phosIndex = self.__currentPhosphateIndex
        
        origPhosCoords = None
        if self.__customPhosphate is not None:
            origPhosCoords = self.__customPhosphate
            self.__drawAdjustPhosCross(self.__customPhosphate)
        else:
            origPhosCoords = self.__candidatePhosphates[phosIndex]
        
        self.__drawBox(origPhosCoords)
        
        window.destroy()
        self.__adjustPhosWindow = None
        #if we don't explicitely do this, the __adjustPhosWindow stays as a GTK window object
        #which leads to problems in __initialPhosAccept
    
    def __adjustInitialPhosAccept(self, widget, window):
        """Respond to the user clicking OK in the adjust initial phosphate dialog box
        
        ARUGMENTS:
            widget - the OK button
            window - the dialog box window
        RETURNS:
            None
        EFFECTS:
            Closes the adjust phosphate dialog box and stores the coordinates that the user has selected
        """
        
        #store the current coordinates
        self.__customPhosphate = self.__translatedPhosCoords
        
        window.destroy()
        self.__adjustPhosWindow = None
        #if we don't explicitely do this, the __adjustPhosWindow stays as a GTK window object
        #which leads to problems in __initialPhosAccept
    
    def __initialPhosCancel(self, widget, window):
        """Cancel the build and erase all objects
        
        ARUGMENTS:
            widget - the cancel button
            window - the RCrane window
        RETURNS:
            None
        EFFECTS:
            Closes the RCrane dialog and destroys all RCrane objects
        """
        #ARGUMENTS:
        #   window - the GUI window object
        
        close_generic_object(self.__candidatePhosObject)
        close_generic_object(self.__batonObject)
        self.__pseudoMolecule.closeBatonObject()
        
        if self.__customPhosObject is not None:
            close_generic_object(self.__customPhosObject)
        
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
        
        graphics_draw()
        window.destroy()
    
    def __initialPhosCloseWin(self, window, widget):
        """Cancel the build and erase all objects
        ARGUMENTS:
            window - the RCrane window
            widget - the cancel button
        RETURNS:
            None
        NOTE:
            This function simply calls __initialPhosCancel with a reversed argument list
            It is used when the window manager issues a close signal (rather than the user
            clicking on the cancel button).
        """
        self.__initialPhosCancel(widget, window)
    
    def __initialPhosAccept(self, widget, window):
        """Accept the coordinates of the currently selected phosphate
        ARGUMENTS:
            widget - the button used to invoke this function
            window - the RCrane window
        RETURNS:
            None
        EFFECTS:
            Creates a new Coot molecule using the accepting phospahate and creates a next nucleotide dialog
        """ 
        
        #if there is still a Manually Adjust window open, click accept in them
        if self.__adjustPhosWindow is not None:
            self.__adjustInitialPhosAccept(None, self.__adjustPhosWindow)
        
        #erase the selection objects
        close_generic_object(self.__candidatePhosObject)
        close_generic_object(self.__batonObject)
        if self.__customPhosObject is not None:
            close_generic_object(self.__customPhosObject)
        graphics_draw()
        #window.destroy()
        
        acceptedPoint = None
        if self.__customPhosphate is None:
            #if there are no manually set coordinates, use the appropriate phosphate index
            acceptedPoint = self.__candidatePhosphates[self.__currentPhosphateIndex][0:3]
        else:
            #if the user has used Manually Adjust, use that
            acceptedPoint = self.__customPhosphate
        
        #add the atom to the pseudoMolecule
        self.__pseudoMolecule.addPhos(acceptedPoint)
        
        set_rotation_centre(acceptedPoint[0], acceptedPoint[1], acceptedPoint[2])
        
        self.selectNextNt(acceptedPoint, window)
    
    ##########################################################################################
    #   Functions for buttons in the next nucleotide GUI (generated by selectNextNt)
    ##########################################################################################
    def __nextNtNextPhos(self, widget):
        """Select the next phosphate in the next nucleotide GUI.
        
        ARGUMENTS:
            widget - the button used to call this function.  Ignored.
        RETURNS:
            None
        EFFECTS:
            Highlights the next phosphate
        """
        
        #if we don't have a custom phosphate position, increment the phosphate index
        if self.__customPhosphate is None:
            self.__currentPhosphateIndex += 1
            if self.__currentPhosphateIndex >= len(self.__candidatePhosphates):
                self.__currentPhosphateIndex = 0
        else:
            self.__customPhosphate = None
            close_generic_object(self.__customPhosObject)
        
        #if there is still a Manually Adjust Phosphate window open, close it
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
        
        #forget about any coordinate adjustments the user has made via the "Manually Adjust" or "Anti/Syn Flip" buttons
        #we must do this before (potentially) calling self.__mutateBase(), otherwise __mutateBase will mutate the customBase
        self.__customPhosphate = None
        self.__customBase      = None
        
        phosIndex = self.__currentPhosphateIndex
        
        prevPhos = self.__getPrevPhosCoords()
        
        #make sure we have C1' predictions for the current phosphate
        if self.__candidateC1s[phosIndex] is None:
            #if we get here, something has gone horribly wrong
            print "WARNING: CANDIDATE C1'S HAVE NOT BEEN PRE-CALCULATED"
            mapNum = imol_refinement_map()
            self.__candidateC1s[phosIndex] = self.__nextPhos.findSugar(mapNum, prevPhos, self.__candidatePhosphates[phosIndex])
        
        c1Index = self.__currentC1Index[phosIndex]
        
        #make sure a base location has been calculated
        if self.__candidateBases[phosIndex][c1Index] is None:
            mapNum = imol_refinement_map()
            self.__candidateBases[phosIndex][c1Index] = self.__nextPhos.findBase(mapNum, self.__candidateC1s[phosIndex][c1Index], prevPhos, self.__candidatePhosphates[phosIndex], self.__baseType, self.__direction)
        else:
            #if it has, make sure that the base is the same as self.__baseType
            self.__mutateBase()
        
        #draw candidate C1' peaks
        self.__drawC1Peaks()
        
        self.__drawNt(prevPhos, self.__candidateC1s[phosIndex][c1Index], self.__candidatePhosphates[phosIndex], self.__candidateBases[phosIndex][c1Index])
    
    def __nextNtPrevPhos(self, widget):
        """Select the previous phosphate in the next nucleotide GUI.
        
        ARGUMENTS:
            widget - the button used to call this function.  Ignored.
        RETURNS:
            None
        EFFECTS:
            Highlights the previous phosphate
        """
        
        #if we don't have a custom phophsate position, decrement the phosphate index
        if self.__customPhosphate is None:
            self.__currentPhosphateIndex -= 1
            if self.__currentPhosphateIndex < 0 :
                self.__currentPhosphateIndex = len(self.__candidatePhosphates) - 1
        else:
            self.__customPhosphate = None
            close_generic_object(self.__customPhosObject)
        
        #if there is still a Manually Adjust Phosphate window open, close it
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
        
        #forget about any coordinate adjustments the user has made via the "Manually Adjust" or "Anti/Syn Flip" buttons
        #we must do this before (potentially) calling self.__mutateBase(), otherwise __mutateBase will mutate the customBase
        self.__customPhosphate = None
        #self.__customC1        = None
        self.__customBase      = None
        
        phosIndex = self.__currentPhosphateIndex
        prevPhos = self.__getPrevPhosCoords()
        
        #make sure we have C1' predictions for the current phosphate
        if self.__candidateC1s[phosIndex] is None:
            #if we get here, something has gone horribly wrong
            print "WARNING: CANDIDATE C1'S HAVE NOT BEEN PRE-CALCULATED"
            mapNum = imol_refinement_map()
            self.__candidateC1s[phosIndex] = self.__nextPhos.findSugar(mapNum, prevPhos, self.__candidatePhosphates[phosIndex])
        
        c1Index = self.__currentC1Index[phosIndex]
        
        #make sure a base location has been calculated
        if self.__candidateBases[phosIndex][c1Index] is None:
            mapNum = imol_refinement_map()
            self.__candidateBases[phosIndex][c1Index] = self.__nextPhos.findBase(mapNum, self.__candidateC1s[phosIndex][c1Index], prevPhos, self.__candidatePhosphates[phosIndex], self.__baseType, self.__direction)
        else:
            #if it has, make sure that the base is the same as self.__baseType
            self.__mutateBase()
        
        #draw candidate C1' peaks
        self.__drawC1Peaks()
        
        self.__drawNt(prevPhos, self.__candidateC1s[phosIndex][c1Index], self.__candidatePhosphates[phosIndex], self.__candidateBases[phosIndex][c1Index])
    
    def __nextNtNextC1(self, widget):
        """Select the next C1' in the next nucleotide GUI.
        
        ARGUMENTS:
            widget - the button used to call this function.  Ignored.
        RETURNS:
            None
        EFFECTS:
            Highlights the next C1'
        """
        
        #if the user has set a custom phosphate, then there's no list of C1' candidates, so don't do anything
        if self.__customPhosphate is not None:
            return
        
        #if there is still a Manually Adjust Phosphate window open, close it
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
        
        #if there is still a Manually Adjust Base window open, close it
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseWindow.destroy()
        
        phosIndex = self.__currentPhosphateIndex
        
        if self.__customBase is None:
            self.__currentC1Index[phosIndex] += 1
            if self.__currentC1Index[phosIndex]>= len(self.__candidateC1s[phosIndex]):
                self.__currentC1Index[phosIndex] = 0
        else:
            self.__customBase = None
        
        prevPhos = self.__getPrevPhosCoords()
        
        c1Index = self.__currentC1Index[phosIndex]
        
        #make sure a base location has been calculated
        if self.__candidateBases[phosIndex][c1Index] is None:
            mapNum = imol_refinement_map()
            self.__candidateBases[phosIndex][c1Index] = self.__nextPhos.findBase(mapNum, self.__candidateC1s[phosIndex][c1Index], prevPhos, self.__candidatePhosphates[phosIndex], self.__baseType, self.__direction)
        else:
            #if it has, make sure that the base is the same as self.__baseType
            self.__mutateBase()
        
        self.__drawNt(prevPhos, self.__candidateC1s[phosIndex][c1Index], self.__candidatePhosphates[phosIndex], self.__candidateBases[phosIndex][c1Index])
    
    def __nextNtPrevC1(self, widget):
        """Select the previous C1' in the next nucleotide GUI.
        
        ARGUMENTS:
            widget - the button used to call this function.  Ignored.
        RETURNS:
            None
        EFFECTS:
            Highlights the previous C1'
        """
        
        #if the user has set a custom phosphate, then there's no list of C1' candidates, so don't do anything
        if self.__customPhosphate is not None:
            return
        
        #if there is still a Manually Adjust Phosphate window open, close it
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
        
        #if there is still a Manually Adjust Base window open, close it
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseWindow.destroy()
        
        phosIndex = self.__currentPhosphateIndex
        
        if self.__customBase is None:
            self.__currentC1Index[phosIndex] -= 1
            if self.__currentC1Index[phosIndex] < 0:
                self.__currentC1Index[phosIndex] = len(self.__candidateC1s[phosIndex]) - 1
        else:
            self.__customBase = None
        
        prevPhos = self.__getPrevPhosCoords()
        
        c1Index = self.__currentC1Index[phosIndex]
        
        #make sure a base location has been calculated
        if self.__candidateBases[phosIndex][c1Index] is None:
            mapNum = imol_refinement_map()
            self.__candidateBases[phosIndex][c1Index] = self.__nextPhos.findBase(mapNum, self.__candidateC1s[phosIndex][c1Index], prevPhos, self.__candidatePhosphates[phosIndex], self.__baseType, self.__direction)
        else:
            #if it has, make sure that the base is the same as self.__baseType
            self.__mutateBase()
        
        self.__drawNt(prevPhos, self.__candidateC1s[phosIndex][c1Index], self.__candidatePhosphates[phosIndex], self.__candidateBases[phosIndex][c1Index])
    
    def __nextNtAccept(self, widget, window):
        """Accept the coordinates of the currently selected nucleotide
        
        ARGUMENTS:
            widget - the button used to invoke this function
            window - the RCrane window
        RETURNS:
            None
        EFFECTS:
            Adds the accepted nucleotide coordinates to the Coot molecule and begins tracing the next nucleotide
        """ 
        
        #if there is still a Manually Adjust window open, click accept
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosAccept(None, self.__adjustPhosWindow)
        
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseAccept(None, self.__adjustBaseWindow)
        
        #erase the selection objects
        close_generic_object(self.__candidatePhosObject)
        close_generic_object(self.__candidateC1Object)
        close_generic_object(self.__batonObject)
        if self.__customPhosObject is not None:
            close_generic_object(self.__customPhosObject)
        
        graphics_draw()
        
        phosIndex = self.__currentPhosphateIndex
        c1Index = self.__currentC1Index[phosIndex]
        
        phosCoords = None
        if self.__customPhosphate is None:
            phosCoords  = self.__candidatePhosphates[phosIndex]
        else:
            phosCoords  = self.__customPhosphate
        #c1coords    = self.__candidateC1s[phosIndex][c1Index]
        
        baseType, baseCoords = (None, None)
        if self.__customBase:
            (baseType, baseCoords) = self.__customBase
        else:
            (baseType, baseCoords) = self.__candidateBases[phosIndex][c1Index]
        
        #add the atom to the pseudoMolecule
        if self.__direction == 3:
            self.__pseudoMolecule.addBaseAndPhos(baseType, baseCoords, phosCoords)
        else: #self.__direction == 5
            self.__pseudoMolecule.addBaseAndPhos5p(baseType, baseCoords, phosCoords)
        
        set_rotation_centre(phosCoords[0], phosCoords[1], phosCoords[2])
        
        #forget about any coordinate adjustments the user has made via the "Manually Adjust" or "Anti/Syn Flip" buttons
        self.__customPhosphate = None
        self.__customBase      = None
        
        self.selectNextNt(phosCoords[0:3], window)
    
    def __nextNtClose(self, widget, window):
        """Close the next nucleotide dialog and delete the molecule
        
        ARGUMENTS:
            widget - the button used to invoke this function
            window - the RCrane window
        RETURNS:
            None
        EFFECTS:
            Destroys the next nucleotide dialog and the Coot molecule
        """
        
        self.__pseudoMolecule.deleteMolecule()
        self.__nextNtFinished(window)
        

    def __nextNtFinished(self, window):
        """Close the next nucleotide window without deleting the molecule
        
        ARGUMENTS:
            window - the RCrane window
        RETURNS:
            None
        EFFECTS:
            Destroys the next nucleotide dialog
        NOTE:
            This function is called by __nextNtClose and __nextNtBuildBackbone.
        """
        
        close_generic_object(self.__candidatePhosObject)
        close_generic_object(self.__candidateC1Object)
        close_generic_object(self.__batonObject)
        if self.__customPhosObject is not None:
            close_generic_object(self.__customPhosObject)
        self.__pseudoMolecule.closeBatonObject()
        graphics_draw()
        
        #if there is still a Manually Adjust Phosphate window open, close it
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
            
        #if there is still a Manually Adjust Base window open, close it
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseWindow.destroy()
        
        
        window.destroy()
    
    def __nextNtBuildBackbone(self, widget, window):
        """Finish tracing the backbone and calculate all-atom backbone coordinates.
        
        ARGUMENTS:
            widget - the button used to invoke this function
            window - the RCrane window
        RETURNS:
            None
        EFFECTS:
            Calculates an initial set of backbone coordinates and creates a review suites GUI dialog
        """
        
        #if there is still a Manually Adjust Phosphate window open, close it
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosWindow.destroy()
            
        #if there is still a Manually Adjust Base window open, close it
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseWindow.destroy()
        
        
        if self.__pseudoMolecule.getNumNts() < 2:
            #if the user has only built a phosphate, then we can't build any backbone
            #so just close the window and leave the phosphate atom
            self.__nextNtFinished(window)
            
        elif self.__pseudoMolecule.getNumNts() == 2:
            #close the graphics for the nucleotide that we haven't built yet
            close_generic_object(self.__candidatePhosObject)
            close_generic_object(self.__candidateC1Object)
            close_generic_object(self.__batonObject)
            
            #if there are two nucleotides, then we have a single sugar with both a 5' and a 3' phosphate
            #we can't determine a conformer, but we can predict the sugar pucker
            builtChain = self.__pseudoMolecule.createChainObject()
            pucker = determinePucker(builtChain.nucleotides[0].pperp())
            #window.disconnect(self.__deleteHandlerID)
            calcCoords(builtChain, pucker, self.__pseudoMolecule, window)
            self.__nextNtFinished(window)
            
        else:
            #close the graphics for the nucleotide that we haven't built yet
            close_generic_object(self.__candidatePhosObject)
            close_generic_object(self.__candidateC1Object)
            close_generic_object(self.__batonObject)
            
            builtChain = self.__pseudoMolecule.createChainObject()
            builtChainOrig = deepcopy(builtChain)
            #for curNuc in builtChain.nucleotides: pprint(curNuc.atoms)
            (bestPath, predictedProbs) = determineRotamerSeq(builtChain)
            print "Suite string: " + "".join(bestPath)
            
            #FOR DEBUGGING
            #print "OPENING foo3.txt"
            #output = open("foo3.txt", 'wb')
            #from cPickle import dump
            #dump([builtChain, bestPath, predictedProbs, self.__pseudoMolecule._PseudoMolecule__molecule[0][0][1], builtChainOrig], output)
            #output.close()
            #sleep(3)
            
            (intermediateAtomLocs, minimizationScores) = calcCoords(builtChain, bestPath, self.__pseudoMolecule, window)
            #print "minimizationScores =", minimizationScores
            #for curNuc in builtChain.nucleotides: pprint(curNuc.atoms)
            self.__pseudoMolecule.closeBatonObject()
            window.disconnect(self.__deleteHandlerID)
            ReviewSuitesGui(predictedProbs, bestPath, self.__pseudoMolecule, builtChainOrig, intermediateAtomLocs, minimizationScores, window)
            #self.__nextNtFinished(window)
        
        
    
    def __nextNtCloseWin(self, window, widget):
        """Respond to the user closing the next nucleotide window 
        
        ARUGMENTS:
            window - the window to close
            widget - the widget used to close the window
        NOTES:
            This function simply calls __nextNtClose with a reversed argument list
        """
        
        #note the different order of arguments between this and __nextNtClose
        self.__nextNtClose(widget, window)
    
    def __nextNtBaseSelect(self, combobox):
        """Mutate the base in response to the user selecting a base type from the drop-down
        
        ARGUMENTS:
            combobox - the drop down object used to select the new base type
        RETURNS:
            None
        EFFECTS:
            the current base is mutated to self.__baseType and redrawn
        """
        
        #determine the base type
        selectionNum = combobox.get_active()
        baseType = ["A", "C", "G", "U"][selectionNum]
        
        self.__baseType = baseType
        self.__mutateBase()
    
    def __nextNtPrevNt(self, widget, window):
        """Cancels building the current nucleotide and moves back one nucleotide
        
        ARGUMENTS:
            widget - the button used to invoke this function
            window - the RCrane window
        RETURNS:
            None
        EFFECTS:
            Undraws the current nucleotide and prompts the user to rebuild the previous nucleotide
        """
        
        #remove the last nucleotide from the pseudoMolecule object
        if self.__direction == 3:
            (prevNtType, prevNtCoords) = self.__pseudoMolecule.removeLastBaseAndPhos()
        else:
            (prevNtType, prevNtCoords) = self.__pseudoMolecule.removeFirstBaseAndPhos()
        #pprint(prevNtType)
        #pprint(prevNtCoords)
        
        #make sure that we're not on the first nucleotide
        if prevNtCoords is None:
            #if we are, then don't do anything
            print "Can't go back any further."
            add_status_bar_text("Can't go back any further.")
            return
        
        #close the current drawing objects
        close_generic_object(self.__candidatePhosObject)
        close_generic_object(self.__candidateC1Object)
        close_generic_object(self.__batonObject)
        if self.__customPhosObject is not None:
            close_generic_object(self.__customPhosObject)
        
        #forget about any coordinate adjustments the user has made via the "Manually Adjust" or "Anti/Syn Flip" buttons
        self.__customPhosphate = None
        self.__customBase      = None
        
        #shift the zoom back
        prevPhos = self.__getPrevPhosCoords()
        set_rotation_centre(prevPhos[0], prevPhos[1], prevPhos[2])
        
        self.selectNextNt(prevPhos, window, prevNtCoords, prevNtType)
        
        
    
    
    ##########################################################################################
    #   Drawing functions
    ########################################################################################## 
    
    def __drawCross(self, coords, crossLength, crossColor, object):
        """Draw a cross of the specified size and color
        
        ARGUMENTS:
            coords      - the coordinates for the center of the cross
            crossLength - the size of one "arm" of the cross (i.e. the distance from the center to the end, not end-to-end)
            crossColor  - the color of the cross
            object      - the Coot drawing object to add this cross to
        RETURNS:
            None
        EFFECTS:
            Adds the specified cross to object
        """
        
        to_generic_object_add_line(object, crossColor, 6, coords[0]-crossLength, coords[1],             coords[2],             coords[0]+crossLength, coords[1],             coords[2])
        to_generic_object_add_line(object, crossColor, 6, coords[0],             coords[1]-crossLength, coords[2],             coords[0],             coords[1]+crossLength, coords[2])
        to_generic_object_add_line(object, crossColor, 6, coords[0],             coords[1],             coords[2]-crossLength, coords[0],             coords[1],             coords[2]+crossLength)
        
    def __drawPhosCross(self, coords):
        """Draw a cross for a phosphate candidate
        
        ARGUMENTS:
            coords - the coordinates for the center of the cross
        RETURNS:
            None
        EFFECTS:
            Adds the cross to self.__candidatePhosObject
            Note that this function does not issue a graphics_draw() command
        """
        self.__drawCross(coords, CROSS_LENGTH_PHOS, CROSS_COLOR_PHOS, self.__candidatePhosObject)
    
    def __drawC1Cross(self, coords):
        """Draw a cross for a C1' candidate
        
        ARGUMENTS:
            coords - the coordinates for the center of the cross
        RETURNS:
            None
        EFFECTS:
            Adds the cross to self.__candidateC1Object
            Note that this function does not issue a graphics_draw() command
        """
        self.__drawCross(coords, CROSS_LENGTH_C1, CROSS_COLOR_C1, self.__candidateC1Object)
    
    def __drawAdjustPhosCross(self, coords):
        """Draw a cross for the manually adjusted phosphate
        
        ARGUMENTS:
            coords - the coordinates for the center of the cross
        RETURNS:
            None
        EFFECTS:
            Deletes the current manually adjusted phosphate cross if one is present and creates a new cross
            Note that this function does not issue a graphics_draw() command
        """
        
        if self.__customPhosObject:
            close_generic_object(self.__customPhosObject)
        
        self.__customPhosObject = new_generic_object_number(ADJUST_PHOS_TITLE)
        set_display_generic_object(self.__customPhosObject, 1)
        
        self.__drawCross(coords, CROSS_LENGTH_PHOS, CROSS_COLOR_PHOS, self.__customPhosObject)
    
    def __drawBox(self, coords):
        """Draw a new box around the specified point
        
        ARGUMENTS:
            coords - the coordinates for the center of the box
        RETURNS:
            None
        EFFECTS:
            Deletes the current box if one is present and creates a new box
            Then issues a graphics_draw() command
        """
        
        #delete the previous box if there is one
        if self.__batonObject:
            close_generic_object(self.__batonObject)
        
        #create a new object for the box
        self.__batonObject = new_generic_object_number(BATON_TITLE)
        set_display_generic_object(self.__batonObject, 1)
        
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]-BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]-BOX_LENGTH,   coords[0]-BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]-BOX_LENGTH)
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]-BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]-BOX_LENGTH,   coords[0]-BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]+BOX_LENGTH)
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]-BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]-BOX_LENGTH,   coords[0]+BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]-BOX_LENGTH)
        
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]-BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]+BOX_LENGTH,   coords[0]-BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]-BOX_LENGTH)
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]-BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]+BOX_LENGTH,   coords[0]-BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]+BOX_LENGTH)
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]-BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]+BOX_LENGTH,   coords[0]+BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]+BOX_LENGTH)
        
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]+BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]-BOX_LENGTH,   coords[0]-BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]-BOX_LENGTH)
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]+BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]-BOX_LENGTH,   coords[0]+BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]-BOX_LENGTH)
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]+BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]-BOX_LENGTH,   coords[0]+BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]+BOX_LENGTH)
        
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]+BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]+BOX_LENGTH,   coords[0]+BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]-BOX_LENGTH)
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]+BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]+BOX_LENGTH,   coords[0]+BOX_LENGTH, coords[1]+BOX_LENGTH, coords[2]+BOX_LENGTH)
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, coords[0]+BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]+BOX_LENGTH,   coords[0]-BOX_LENGTH, coords[1]-BOX_LENGTH, coords[2]+BOX_LENGTH)
        
        
        graphics_draw()
    
    
    def __drawNt(self, phos5, c1, phos3, base):
        """Draw batons for a new nucleotide candidate
        
        ARGUMENTS:
            phos5 - the coordinates of the 5' phosphate
            c1    - the coordinates of the C1'
            phos3 - the coordinates of the 3' phosphate
            base  - a 2 item list consisting of (baseType, baseCoords)
                    baseType is either C, U, G, or A
                    baseCoords is a dictionary of atom type => coordinates
        EFFECTS:
            Draws the specified nucleotide
            Issues a graphics_draw() command
        """
        
        #delete the previous batonObject if there is one
        if self.__batonObject:
            close_generic_object(self.__batonObject)
        
        #create a new object for the nt
        self.__batonObject = new_generic_object_number(BATON_TITLE)
        set_display_generic_object(self.__batonObject, 1)
        
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, phos5[0], phos5[1], phos5[2], c1[0], c1[1], c1[2])
        to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, c1[0], c1[1], c1[2], phos3[0], phos3[1], phos3[2])
        
        #draw the base
        (baseType, baseCoords) = base
        
        bondList = None
        if baseType == "C":
            bondList = [["C1'","N1"], ["N1","C2"], ["C2","N3"], ["N3","C4"], ["C4","C5"], ["C5","C6"], ["C6","N1"], ["C2","O2"], ["C4","N4"]]
        elif baseType == "U":
            bondList = [["C1'","N1"], ["N1","C2"], ["C2","N3"], ["N3","C4"], ["C4","C5"], ["C5","C6"], ["C6","N1"], ["C2","O2"], ["C4","O4"]]
        elif baseType == "G":
            bondList = [["C1'","N9"], ["N9","C8"], ["C8","N7"], ["N7","C5"], ["C5","C4"], ["C4","N3"], ["N3","C2"], ["C2","N1"], ["N1","C6"], ["C6","C5"], ["C4","N9"], ["C2","N2"], ["C6","O6"]]
        elif baseType == "A":
            bondList = [["C1'","N9"], ["N9","C8"], ["C8","N7"], ["N7","C5"], ["C5","C4"], ["C4","N3"], ["N3","C2"], ["C2","N1"], ["N1","C6"], ["C6","C5"], ["C4","N9"], ["C6","N6"]]
        
        for bondAtoms in bondList:
            (atom1, atom2) = bondAtoms
            atom1Coords = baseCoords[atom1]
            atom2Coords = baseCoords[atom2]
            to_generic_object_add_line(self.__batonObject, BATON_COLOR, 6, atom1Coords[0], atom1Coords[1], atom1Coords[2], atom2Coords[0], atom2Coords[1], atom2Coords[2])
        
        graphics_draw()
    
    
    def __drawPhosPeaks(self):
        """Draw crosses for all phosphate candidates
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            Draws crosses for each point in self.__candidatePhosphates
            Note that this function does not issue a graphics_draw() command
        """
        #create an object to display the candidate peaks
        self.__candidatePhosObject = new_generic_object_number(CANDIDATE_PHOS_TITLE)
        set_display_generic_object(self.__candidatePhosObject, 1)
        
        #draw crosses for each candidate peak
        for curPeak in self.__candidatePhosphates:
            self.__drawPhosCross(curPeak)
    
    def __drawC1Peaks(self):
        """Draw crosses for all C1' candidates for the currently selected phosphate
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            Draws crosses for each point in self.__candidateC1s for the current phosphate
            Note that this function does not issue a graphics_draw() command
        """
        
        #close the previous C1' candidate object if it exists
        if self.__candidateC1Object is not None and not (is_closed_generic_object_p(self.__candidateC1Object)):
            close_generic_object(self.__candidateC1Object)
        
        #create an object to display the candidate peaks
        self.__candidateC1Object = new_generic_object_number(CANDIDATE_C1_TITLE)
        set_display_generic_object(self.__candidateC1Object, 1)
        
        #draw crosses for each candidate peak
        for curPeak in self.__candidateC1s[self.__currentPhosphateIndex]:
            self.__drawC1Cross(curPeak)
    
    
    def __switchDirection(self, widget, prevPhos, window):
        """Switch the direction of the trace ( from 3'->5' to 5'->3', or vice-versa)
        ARGUMENTS:
            widget   - the button used to call this function
            prevPhos - the coordinates of the previous phosphate atom
            window   - the RCrane window
        RETURNS:
            None
        """
        
        if self.__direction == 3:
            self.__direction = 5
        else:
            self.__direction = 3
        
        
        #get the currently selected phosphate coordinates so we can preserve them after the switch
        phosCoords = None
        if self.__customPhosphate is None:
            phosCoords  = self.__candidatePhosphates[self.__currentPhosphateIndex]
        else:
            phosCoords  = self.__customPhosphate
        #note that we can't pass the phosphate index itself, since the phosphates will likely be ranked differently
        #after the direction switch.
        #We could pass information about whether it was a custom phosphate position or not, but that wouldn't likely
        #save any significant amount of time, and it could lead to problems later if I ever implement a more clever
        #phosphate finding algorithm where potential phosphate locations may be filtered before being presented to the user
        
        
        #if there is still a Manually Adjust window open, click accept
        if self.__adjustPhosWindow is not None:
            self.__adjustPhosAccept(None, self.__adjustPhosWindow)
        
        if self.__adjustBaseWindow is not None:
            self.__adjustBaseAccept(None, self.__adjustBaseWindow)
        
        #erase any drawing objects
        close_generic_object(self.__candidatePhosObject)
        close_generic_object(self.__candidateC1Object)
        close_generic_object(self.__batonObject)
        if self.__customPhosObject is not None:
            close_generic_object(self.__customPhosObject)
        graphics_draw()
        
        
        self.selectNextNt(prevPhos, window, initialPhos = phosCoords)
    
    
    def __getPrevPhosCoords(self):
        """Get the coordinates of the previous phosphate
        
        ARGUMENTS:
            None
        RETURNS:
            the coordinates of the previous phosphate (phosphate -1 if tracing 5'->3', phosphate 1 if tracing 3'->5')
        """
        
        
        if self.__direction == 3:
            return self.__pseudoMolecule.getPhosCoords(-1)
        else:
            return self.__pseudoMolecule.getPhosCoords(1)
