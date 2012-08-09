#!/usr/bin/env python
"""A class for the graphical interface used to review the predicted conformers."""

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

import gtk
import pango
from copy import deepcopy
#from pprint import pprint     #for debugging
#from cPickle import dump      #for debugging
import os
#import sys
from time import sleep      #for debugging

from coot import set_rotation_centre, add_atom_label, remove_atom_label, graphics_draw, close_generic_object, set_display_generic_object

from puckerList import puckerListFullName, puckerList
from calcCoords import calcCoords, recalcCoords
from guiUtils import buttonWithIcon, HBOX_SPACING, VBOX_SPACING, BUTTON_SPACING, createRCraneWindowObject
from rotamerSeq import determineAlternateConf

DIFFERENT_PUCKER_COLOR = "pink"
BEST_PATH_COLOR = "palegreen"
BOLD_WEIGHT = pango.WEIGHT_HEAVY


class ReviewSuitesGui:
    """A class for the graphical interface used to review the predicted conformers."""
    
    def __init__(self, predictedProbs, bestPath, pseudoMol, origCoords, intermediateAtomLocs, minimizationScores, window = None, useResNums = True, oldCoordsObj = None, ignoreDensity = False, prevTorsions = None, nextTorsions = None):
        """Initialize a ReviewSuitesGui object and create a GUI for reviewing the predicted suites
        
        ARGUMENTS:
            predictedProbs       - the probabilities for each rotamer calculated for each suite
                                   formatted as a list of suites, where each suite is represented as a dictionary of rotamer => probability
            bestPath             - a list of the most likely rotamers for each suite
                                   This is typically calculed by rotamerSeq using an HMM, so it may differ slightly from simply taking all the maxima from predictedProbs
            pseudoMol            - a PseudoMolecule object for the built molecule
            origCoords           - a copy of the coordinates contained in pseudoMol as a Chain object
            intermediateAtomLocs - phosphate, O3', and C3' coordinates for each nucleotide immediately before minimization of that nucleotide was started
            minimizationScores   - the minimization scores for each nucleotide
        OPTIONAL ARGUMENTS:
            window               - a window to put the GUI into
                                   If not provided, a new window will be created.
            useResNums           - whether or not to display the rotamer numbers (as opposed to just the indices) in the GUI
                                   Defaults to True
            oldCoordsObj         - a Coot generic object for drawing the coordinates we are potentially replacing when running rotamerize
                                   If not provided, then the toggle checkbox will not be displayed and the generic object will not
                                   be destroyed when closing the window.
            ignoreDensity        - ignore the density when performing the minimization
                                   defaults to False
            prevTorsions         - the pre-rotamerize torsions for the suite before the start of origCoords.  These torsions will be used
                                   as restraints when rebuilding the first suite
                                   defaults to None
            nextTorsions         - the pre-rotamerize torsions for the suite after the end of origCoords.  These torsions will be used
                                   as restraints when rebuilding the first suite
                                   defaults to None
        RETURNS:
            an initialized ReviewSuitesGui object
        EFFECTS:
            Displays an RCrane dialog box for reviewing the predicted suites
        NOTE:
            If the pseudoMol object has savedCoordinates, then an an "Accept New Coordinates" and a "Restore Original Coordinates" will be displayed.
            Otherwise, only a "Done" button will be displayed
        """
        
        self.__predictedProbs = predictedProbs
        self.__bestPath = bestPath
        self.__pseudoMolecule = pseudoMol
        self.__origCoords = origCoords
        self.__intermediateAtomLocs = intermediateAtomLocs
        self.__window = window
        self.__useResNums = useResNums
        self.__oldCoordsObj = oldCoordsObj
        self.__minimizationScores = minimizationScores
        self.__ignoreDensity = ignoreDensity
        self.__prevTorsions = prevTorsions
        self.__nextTorsions = nextTorsions
        #print minimizationScores
        
        self.__curPath  = deepcopy(bestPath)
        self.__curSuite = 1
        self.__numNucs = len(self.__bestPath) + 2
        self.__treeViewChangedHandlerID = None
        
        #if the initial nucleotide doesn't have a C1' atom, then we can't rebuilt it
        #(This happens when we're rotamerizing.  The origCoords object contains a few 5' atoms.)
        #as a result we set __skipFirstNuc so we know to ignore the first nucleotide
        if origCoords.nucs[0].hasAtom("C1'"):
            self.__skipFirstNuc = 0
        else:
            self.__skipFirstNuc = 1
        
        #use clearCache() to initialize the cache variables
        self.__clearCache()
        
        #print origCoords
        
        #output = open("foo2.txt", 'wb')
        #dump([predictedProbs, bestPath, pseudoMol._PseudoMolecule__molecule[0][0][1], origCoords, intermediateAtomLocs, minimizationScores], output)
        #output.close()
        
        #create the review suite GUI
        self.__reviewSuites()
        
    
    def __reviewSuites(self):
        """Create the review suites GUI
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            Displays an RCrane dialog box for reviewing the predicted suites
        """
        
        if self.__window is None:
            self.__window = createRCraneWindowObject()
        elif self.__window is not None:
            #remove the current contents of the window
            windowChild = self.__window.get_child()
            if windowChild is not None:
                self.__window.remove(windowChild)
            self.__window.resize(1,1) #remove any previous size constraints
        
        #connect the close window signal
        self.__window.connect("delete_event", self.__closeWin)
        
        
        # Create the buttons
        prevConfButton      = buttonWithIcon("  Previous Conformer  ",    "gtk-go-up")
        nextConfButton      = buttonWithIcon("  Next Conformer  ",        "gtk-go-down")
        
        prevSuiteButton     = buttonWithIcon("  Previous Suite  ",    "gtk-media-previous")
        nextSuiteButton     = buttonWithIcon("  Next Suite  ",        "gtk-media-next")
        
        if self.__pseudoMolecule.hasSavedCoordinates():
            doneButton     = buttonWithIcon("  Accept New Coordinates  ",      "gtk-ok")
            cancelButton   = buttonWithIcon("  Restore Original Coordinates  ", "gtk-cancel")
        else:
            doneButton     = buttonWithIcon("  Done  ",      "gtk-ok")
        
        
        #create the suiteFrame, since that needs to be passed to the previous and next suite buttons
        suiteNum = self.__origCoords.nucs[1 + self.__skipFirstNuc].resNum
        if self.__useResNums:
            frameTitle = "Suite %s (1 of %i)" % (str(suiteNum), len(self.__predictedProbs))
        else:
            frameTitle = "Suite 1 of %i" % len(self.__predictedProbs)
        
        suiteFrame = gtk.Frame(frameTitle)
        
        
        #create the TreeView objects
        (confListStore, initialConfNum) = self.__createConfListStore(1)
        (confTreeView, confTreeSelection, rowHeight) = self.__createConfTreeView(confListStore, initialConfNum)
        if len(self.__bestPath) > 1:
            overviewTreeView = self.__createOverviewTreeView(confTreeView, suiteFrame)
        else:
            overviewTreeView = None
        self.__treeViewChangedHandlerID = confTreeSelection.connect("changed", self.__changed, overviewTreeView)
        
        #create the ScrolledWindow objects used to display the TreeViews
        (treeViewScroll, treeViewScrollWidth) = self.__createConfScrolledWindow(confTreeView, rowHeight)
        if len(self.__bestPath) > 1:
            (overviewScroll, wideOverview) = self.__createOverviewScrolledWindow(overviewTreeView, prevConfButton, treeViewScrollWidth)
        
        #connect the buttons
        prevConfButton.connect("clicked", self.__prevConf, confTreeSelection)
        nextConfButton.connect("clicked", self.__nextConf, confTreeSelection)
        
        #we don't next previous and next suite buttons if there's only a single suite
        if len(self.__bestPath) > 1:
            prevSuiteButton.connect("clicked", self.__prevSuite, confTreeView, suiteFrame, overviewTreeView)
            nextSuiteButton.connect("clicked", self.__nextSuite, confTreeView, suiteFrame, overviewTreeView)
        
        doneButton.connect("clicked", self.__accept)
        
        if self.__pseudoMolecule.hasSavedCoordinates():
            cancelButton.connect("clicked", self.__restoreOldCoords)
        
        
        #put the new contents into the window
        windowBox = gtk.VBox(False, VBOX_SPACING)
        
        #create the suite overview frame if there is more than one suite
        if len(self.__bestPath) > 1:
            overviewFrame = gtk.Frame("Overview")
            overviewBox = gtk.HBox(True, 0)
            overviewBox.pack_start(overviewScroll, True, wideOverview, 0)
            overviewFrame.add(overviewBox)
            windowBox.pack_start(overviewFrame, False, True, BUTTON_SPACING)

        
        suiteFrameHBox = gtk.HBox(False, HBOX_SPACING)
        suiteFrameVBox = gtk.VBox(False, VBOX_SPACING)
        suiteFrame.add(suiteFrameVBox)
        suiteFrameVBox.pack_start(suiteFrameHBox, True, True, BUTTON_SPACING)
        windowBox.pack_start(suiteFrame, True, True, BUTTON_SPACING)
        
        prevNextConfBox = gtk.VBox(True, VBOX_SPACING)
        prevNextConfAlignment = gtk.Alignment(xalign=0.5, yalign=0.5)
        prevNextConfAlignment.add(prevNextConfBox)
        suiteFrameHBox.pack_start(prevNextConfAlignment, False, False, BUTTON_SPACING)
        prevNextConfBox.pack_start(prevConfButton, True, True, BUTTON_SPACING)
        prevNextConfBox.pack_start(nextConfButton, True, True, BUTTON_SPACING)
        
        #create a checkbox for displaying the old coordinates if we have a generic object
        if self.__oldCoordsObj is not None:
            oldCoordsCheckButton = gtk.CheckButton("Display old coordinates")
            oldCoordsCheckButton.set_active(True)
            oldCoordsCheckButton.connect("toggled", self.__toggleOldCoords)
            oldCoordsAlignment = gtk.Alignment(xalign = 0.5)
            oldCoordsAlignment.add(oldCoordsCheckButton)
            prevNextConfBox.pack_start(oldCoordsAlignment, True, True, BUTTON_SPACING)
        
        suiteFrameHBox.pack_start(treeViewScroll, True, False, BUTTON_SPACING)
        
        #we don't next previous and next suite buttons if there's only a single suite
        if len(self.__bestPath) > 1:
            prevNextSuiteBox = gtk.HBox(True, HBOX_SPACING)
            prevNextSuiteAlignment = gtk.Alignment(xalign=0.5)
            prevNextSuiteAlignment.add(prevNextSuiteBox)
            prevNextSuiteBox.pack_start(prevSuiteButton, True, True, BUTTON_SPACING)
            prevNextSuiteBox.pack_start(nextSuiteButton, True, True, BUTTON_SPACING)
            #If this is not a rotamerize call, then we put the previous and next suite buttons outside of the "Suite X of X"
            #frame.  If it is a rotamerize call, then the buttons go inside of the frame.  This is inconsistant, but it makes
            #things look a little more balanced (and it's a minor enough thing that I doubt anyone will ever notice the
            #inconsistancy)
            if self.__pseudoMolecule.hasSavedCoordinates():
                suiteFrameVBox.pack_start(prevNextSuiteAlignment, False, False, BUTTON_SPACING)
            else:
                windowBox.pack_start(prevNextSuiteAlignment, False, False, BUTTON_SPACING)
        
        #if this is a rotamerize call, then we pack the "Accept New Coordinates" and the "Restore Original Coordinates"
        #buttons
        if self.__pseudoMolecule.hasSavedCoordinates():
            doneBox = gtk.HBox(True, HBOX_SPACING)
            doneAlignment = gtk.Alignment(xalign = 0.5)
            doneAlignment.add(doneBox)
            windowBox.pack_start(doneAlignment, False, False, BUTTON_SPACING)
            doneBox.pack_start(doneButton, True, True, BUTTON_SPACING)
            doneBox.pack_start(cancelButton, True, True, BUTTON_SPACING)
            
        #if this isn't a rotamerize call, then we just pack the "Done" button
        else:
            doneAlignment = gtk.Alignment(xalign = 0.5)
            doneAlignment.add(doneButton)
            windowBox.pack_start(doneAlignment, False, False, BUTTON_SPACING)
        
        #set the current rot
        self.__selectedConf = self.__curPath[0]
        
        #go to the center of the first suite
        firstPhos = self.__pseudoMolecule.getPhosCoords(suiteNum)
        set_rotation_centre(*firstPhos)
        add_atom_label(self.__pseudoMolecule.molNum(), self.__pseudoMolecule.chain, int(suiteNum), " P  ") #label the phosphate
            #Note that add_atom_label can't handle insertion codes, so this command will label the wrong atom if we are
            #trying to center on a suite with an insertion code
        graphics_draw()
        
        #display the new window
        self.__window.add(windowBox)
        self.__window.show_all()
        
        if len(self.__bestPath) > 1:
            nextSuiteButton.grab_focus() #make the Next Suite button the default focus if it exists
        else:
            nextConfButton.grab_focus()  #otherwise, make the Next Conformer button the default focus
        
        
    def __changed(self, confTreeSelection, overviewTreeView):
        """Respond to the user selecting a different rotamer
        
        ARGUMENTS:
            confTreeSeletion - the TreeSelection object for the conformer TreeView
        RETURNS:
            None
        EFFECTS:
            Recalculates and reminimizes coordinates using the newly selected conformer for the current suite
        """
        
        #figure out what nucleotides we need to cache
        startingNucIndex = self.__curSuite + self.__skipFirstNuc - 1
        if startingNucIndex <= 0:
            startingNucIndex = 0
            hasPrevSuite = False #there's no suite before this one, so we don't have to worry about changes to the starting pucker
        else:
            hasPrevSuite = True
        endingNucIndex = self.__curSuite + self.__skipFirstNuc + 1
        #if endingNuc > self.__numNucs: endingNuc = self.__numNucs
        
        startingNucNum = self.__origCoords.nucs[startingNucIndex].resNum
        endingNucNum   = self.__origCoords.nucs[endingNucIndex].resNum
        
        #print "\n\n\n*****************__changed called\n*************************\n**************************"
        #print "startingNucIndex =", startingNucIndex
        #print "endingNucIndex =", endingNucIndex
        #print "hasPrevSuite =", hasPrevSuite        
        #print "startingNucNum =", startingNucNum
        #print "endingNucNum   =", endingNucNum
        #from time import sleep; sleep(3)
        
        previousConf = self.__selectedConf
        
        #if need be, cache the existing structure
        if not(self.__structureCache.has_key(previousConf)):
            #print "Storing structure for", previousConf
            self.__structureCache[previousConf] = self.__pseudoMolecule.getCootNucs(startingNucNum, endingNucNum)
        
        
        #figure out what the newly selected conformer is
        (confListStore, selectedConfIter) = confTreeSelection.get_selected()
        selectedConf = confListStore.get_value(selectedConfIter, 0)
        print "Conformer %s was selected" % selectedConf
        self.__selectedConf = selectedConf
        
        #see if we're changing the starting pucker (and we're not at the first suite or the first suite of a rotamerize)
        if hasPrevSuite and puckerList[previousConf][0] != puckerList[selectedConf][0] and not (startingNucIndex == 1 and self.__prevTorsions):
            
            preStartingNucNum = self.__origCoords.nucs[startingNucIndex-1].resNum
            #print "preStartingNucNum =", startingNucNum
            
            #if we are, then see if we need to cache the current structure
            if not self.__startingPuckerStructureCache.has_key(puckerList[previousConf][0]):
                self.__startingPuckerStructureCache[puckerList[previousConf][0]] = self.__pseudoMolecule.getCootNucs(preStartingNucNum, startingNucNum)
            
            #then see if we have the structure we need in the cache
            if self.__startingPuckerStructureCache.has_key(puckerList[selectedConf][0]):
                self.__pseudoMolecule.setCootNucs(preStartingNucNum, startingNucNum, self.__startingPuckerStructureCache[puckerList[selectedConf][0]])
                changingStartingPucker = False
            else:
                changingStartingPucker = puckerList[selectedConf][0]
            
            #see if we've calculated the appropriate conformer for the previous suite
            if len(self.__previousSuiteConf) == 2:
                previousSuiteConf = self.__previousSuiteConf[puckerList[selectedConf][0]]
            else:
                #if not, then calculate it now
                self.__previousSuiteConf[puckerList[previousConf][0]] = self.__curPath[self.__curSuite-2]
                
                previousSuiteLeadingPucker = puckerList[self.__curPath[self.__curSuite-2]][0]
                previousSuiteConf = determineAlternateConf(previousSuiteLeadingPucker, puckerList[selectedConf][0], self.__curSuite-2, self.__predictedProbs)
                self.__previousSuiteConf[puckerList[selectedConf][0]] = previousSuiteConf
                
        else:
            changingStartingPucker = False
            if self.__curSuite > 1:
                previousSuiteConf = self.__curPath[self.__curSuite-2]
            elif self.__prevTorsions:
                previousSuiteConf = self.__prevTorsions
            else:
                previousSuiteConf = None
        
        #print "changingStartingPucker =", changingStartingPucker
        #print "previousSuiteConf =", previousSuiteConf
        #print "self.__previousSuiteConf = ", self.__previousSuiteConf
        
        
        #see if we're changing the ending pucker
        #print "endingNuc =", endingNuc
        #print "self.__numNucs =", self.__numNucs
        if endingNucIndex == (self.__numNucs - 1 + self.__skipFirstNuc):
            #if this is the last suite, then we don't really care what the ending pucker is when dealing with caching the next nucleotide
            #so we're only going to store one cached structure for both 2' and 3' ending pucker
            #to see if we've already stored a structure, we check to see if __endingPuckerStructureCache is a dictionary (since it gets reset to an empty dictionary in __nextSuite)
            #we also only need to cache a single nucleotide, since there's no endingNuc+1 nucleotide
            if isinstance(self.__endingPuckerStructureCache, dict):
                self.__endingPuckerStructureCache = self.__pseudoMolecule.getCootNucs(endingNucNum, endingNucNum)
            else:
                #if __endingPuckerStructureCache isn't a dictionary, then we've already stored a structure there and we should restore it
                self.__pseudoMolecule.setCootNucs(endingNucNum, endingNucNum, self.__endingPuckerStructureCache)
            
            changingEndingPucker = False
            
            if self.__nextTorsions:
                nextSuiteConf = self.__nextTorsions
            else:
                nextSuiteConf = None
            
        elif puckerList[previousConf][1] != puckerList[selectedConf][1]:
            
            postEndingNucNum = self.__origCoords.nucs[endingNucIndex+1].resNum
            #print "postEndingNucNum =", postEndingNucNum
            
            #if we're changing the ending pucker and this isn't the last suite, then see if we need to cache the current structure
            if not self.__endingPuckerStructureCache.has_key(puckerList[previousConf][1]):
                self.__endingPuckerStructureCache[puckerList[previousConf][1]] = self.__pseudoMolecule.getCootNucs(endingNucNum, postEndingNucNum)
            
            #then see if we have the structure we need in the cache
            if self.__endingPuckerStructureCache.has_key(puckerList[selectedConf][1]):
                self.__pseudoMolecule.setCootNucs(endingNucNum, postEndingNucNum, self.__endingPuckerStructureCache[puckerList[selectedConf][1]])
                changingEndingPucker = False
            else:
                changingEndingPucker = puckerList[selectedConf][1]
            
            #see if we've calculated the appropriate conformer for the next suite
            if len(self.__nextSuiteConf) == 2:
                nextSuiteConf = self.__nextSuiteConf[puckerList[selectedConf][1]]
            else:
                #if not, then calculate it now
                self.__nextSuiteConf[puckerList[previousConf][1]] = self.__curPath[self.__curSuite]
                
                nextSuiteEndingPucker = puckerList[self.__curPath[self.__curSuite]][1]
                nextSuiteConf = determineAlternateConf(puckerList[selectedConf][1], nextSuiteEndingPucker, self.__curSuite, self.__predictedProbs)
                self.__nextSuiteConf[puckerList[selectedConf][1]] = nextSuiteConf
        else:
            
            postEndingNucNum = self.__origCoords.nucs[endingNucIndex+1].resNum
            #print "postEndingNucNum =", postEndingNucNum
            
            #even if we're not changing the ending pucker, we should still reset the structure of the next nucleotide back to how it was before we started changing rotamers
            if self.__endingPuckerStructureCache.has_key(puckerList[previousConf][1]):
                self.__pseudoMolecule.setCootNucs(endingNucNum, postEndingNucNum, self.__endingPuckerStructureCache[puckerList[previousConf][1]])
            else:
                #if the structure hasn't been cached, then store the current structure.  That means that this is the first new conformer we've selected, so we the structure is already how we want it
                self.__endingPuckerStructureCache[puckerList[previousConf][1]] = self.__pseudoMolecule.getCootNucs(endingNucNum, postEndingNucNum)
            
            changingEndingPucker = False
            if self.__curSuite < (self.__numNucs - 2):
                nextSuiteConf = self.__curPath[self.__curSuite]
            elif self.__nextTorsions:
                nextSuiteConf = self.__nextTorsions
            else:
                nextSuiteConf = None
        
        #print ""
        #print "changingEndingPucker =", changingEndingPucker
        #print "nextSuiteConf =", nextSuiteConf
        #print "self.__nextSuiteConf = ", self.__nextSuiteConf
        #pprint(self.__previousSuiteConf)
        
        
        #see if we have the structure for this conformer in the cache
        if self.__structureCache.has_key(selectedConf):
            #if we do, then just use that structure
            #print "Restoring structure for", selectedConf
            self.__pseudoMolecule.setCootNucs(startingNucNum, endingNucNum, self.__structureCache[selectedConf])
            self.__pseudoMolecule.drawExtraBonds()
        else:
            #if we don't, then we have to calculate one
            
            rotList = [previousSuiteConf, selectedConf, nextSuiteConf]
            #note that startingNuc and endingNuc start numbering at 0 (whereas Coot starts at 1, so we have to add one when we send these variables to recalcCoords)
            minimizationStartingNucIndex = self.__curSuite + self.__skipFirstNuc - 1
            minimizationEndingNucIndex = self.__curSuite + self.__skipFirstNuc 
            
            #print "minimizationStartingNucIndex =   ", minimizationStartingNucIndex
            #print "minimizationEndingNucIndex =     ", minimizationEndingNucIndex
            #print "changingStartingPucker =         ", changingStartingPucker
            #print "changingEndingPucker =           ", changingEndingPucker
            
            #see if we need to change any suites other than the current one
            if changingStartingPucker:
                rotList.insert(0, self.__curPath[self.__curSuite-3])
                minimizationStartingNucIndex -= 1
            
            if changingEndingPucker:
                if endingNucIndex < (self.__numNucs - 2 + self.__skipFirstNuc):
                    rotList.append(self.__curPath[self.__curSuite+1])
                else:
                    rotList.append(None)
                minimizationEndingNucIndex += 1
                
            minimizationStartingNucNum = self.__origCoords.nucs[minimizationStartingNucIndex].resNum
            minimizationEndingNucNum   = self.__origCoords.nucs[minimizationEndingNucIndex].resNum
            
            #reset the necessary nucleotides so that we can restart the minimization
            self.__pseudoMolecule.resetNucs(minimizationStartingNucNum, minimizationEndingNucNum, self.__origCoords, self.__intermediateAtomLocs)
            self.__pseudoMolecule.clearExtraBonds()
            
            recalcCoords(minimizationStartingNucNum, minimizationEndingNucNum, rotList, self.__origCoords, self.__pseudoMolecule, self.__window, self.__ignoreDensity)
            self.__pseudoMolecule.drawExtraBonds()
            
            #if necessary, restore the B-factors for the rebuilt atoms
            if self.__pseudoMolecule.hasSavedBfacs():
                self.__pseudoMolecule.restoreSavedBfacs(minimizationStartingNucNum, minimizationEndingNucNum)
                
        
        #update the conformers on the suite overview table
        if overviewTreeView is not None:
            self.__changeOverviewConformer(overviewTreeView, self.__curSuite, selectedConf)
            if isinstance(previousSuiteConf, basestring):
                #update the overview table if we've changed the previous conformer
                #previousSuiteConf can also be a list of torsions if we're changing the first nt of a rotamerization
                #but that shouldn't be reflected in the overview table
                self.__changeOverviewConformer(overviewTreeView, self.__curSuite-1, previousSuiteConf)
            if isinstance(nextSuiteConf, basestring):
                #update the overview table if we've changed the next conformer (ignoring lists as above)
                self.__changeOverviewConformer(overviewTreeView, self.__curSuite+1, nextSuiteConf)
    
    
    def __prevConf(self, widget, confTreeSelection):
        """Respond to the user clicking on the previous conformer button
        
        ARGUMENTS:
            widget           - the button used to invoke this function
            confTreeSeletion - the TreeSelection object for the conformer TreeView
        RETURNS:
            None
        EFFECTS:
            Highlights the previous conformer (which will implicitly invoke __changed)
        """
        
        (confListStore, selectedConfIter) = confTreeSelection.get_selected()
        
        #there is no iter_prev, so we have to jump through some hoops to get the previous row
        curPos = confListStore.get_path(selectedConfIter)[0]
        if curPos != 0:
            newIter = confListStore.get_iter((curPos-1))
            confTreeSelection.select_iter(newIter)
            
            confTreeView = confTreeSelection.get_tree_view()
            confTreeView.scroll_to_cell((curPos-1))
    
    
    def __nextConf(self, widget, confTreeSelection):
        """Respond to the user clicking on the next conformer button
        
        ARGUMENTS:
            widget           - the button used to invoke this function
            confTreeSeletion - the TreeSelection object for the conformer TreeView
        RETURNS:
            None
        EFFECTS:
            Highlights the next conformer (which will implicitly invoke __changed)
        """
        
        (confListStore, selectedConfIter) = confTreeSelection.get_selected()
        newIter = confListStore.iter_next(selectedConfIter)
        if newIter is not None:
            confTreeSelection.select_iter(newIter)
            
            confTreeView = confTreeSelection.get_tree_view()
            confTreeView.scroll_to_cell(confListStore.get_path(newIter))
    
    
    def __createConfListStore(self, curSuite):
        """Creates a ListStore containing data about the specified suite
        
        ARGUMENTS:
            curSuite - the suite to create a ListStore for
        RETURNS:
            confListStore  - a ListStore containing data about curSuite
            initialConfNum - the most likely suite (i.e. the row of confListStore to initially highlight)
        """
        
        curSuite -= 1 #the arrays start at 0, not 1
        
        confListStore = gtk.ListStore(str, str, str, str, 'gboolean', 'gboolean', 'gboolean', int)
        
        conformerList = self.__predictedProbs[curSuite].items()
        conformerList.sort(key=lambda(x): x[1], reverse=True)
        (defaultLeadingPucker, defaultEndingPucker) = puckerList[self.__bestPath[curSuite]]
        i = 0
        for (conf, score) in conformerList:
            confListStore.append([conf, "%4.2f" % score] + puckerListFullName[conf] + [(puckerList[conf][0] != defaultLeadingPucker), (puckerList[conf][1] != defaultEndingPucker), (conf == self.__bestPath[curSuite]), i+1])
            if conf == self.__curPath[curSuite]: initialConfNum = i
            i += 1
        
        return (confListStore, initialConfNum)
    
    
    def __createConfTreeView(self, confListStore, initialConfNum):
        """Create the TreeView object for the list of conformers
        
        ARGUMENTS:
            confListStore     - the ListStore object to attach to the TreeView
            initialConfNum    - the most likely suite (i.e. the row of confListStore to initially highlight)
        RETURNS:
            confTreeView      - the newly created TreeView object
            confTreeSelection - the TreeSelection object associated with the TreeView
            rowHeight         - the height of a single row in the TreeView (used to determine how high to size the ScrolledWindow)
        """
        
        #then create the columns and pack them into the tree view object
        confTreeView  = gtk.TreeView(confListStore)
        confTreeView.set_enable_search(False)
        
        renderRank = gtk.CellRendererText()
        renderRank.set_property('xalign', 0.5)
        renderRank.set_property('cell-background-gdk', self.__window.get_style().bg[gtk.STATE_NORMAL])
        columnRank = gtk.TreeViewColumn("Rank", renderRank, text=7)
        confTreeView.append_column(columnRank)
        
        renderConf = gtk.CellRendererText()
        renderConf.set_property('xalign', 0.5)
        renderConf.set_property('cell-background', BEST_PATH_COLOR)
        columnConf = gtk.TreeViewColumn("Conformer", renderConf, text=0, cell_background_set=6)
        confTreeView.append_column(columnConf)
        
        renderScore = gtk.CellRendererText()
        renderScore.set_property('xalign', 0.5)
        renderScore.set_property('cell-background', BEST_PATH_COLOR)
        columnScore = gtk.TreeViewColumn("Score", renderScore, text=1, cell_background_set=6)
        confTreeView.append_column(columnScore)
        
        renderLeadingPucker = gtk.CellRendererText()
        renderLeadingPucker.set_property('xalign', 0.5)
        columnLeadingPucker = gtk.TreeViewColumn("Leading Pucker", renderLeadingPucker, text=2)
        columnLeadingPucker.set_cell_data_func(renderLeadingPucker, self.__setPuckerBg, (4, 6))
        confTreeView.append_column(columnLeadingPucker)
        
        renderEndingPucker = gtk.CellRendererText()
        renderEndingPucker.set_property('xalign', 0.5)
        columnEndingPucker = gtk.TreeViewColumn("Ending Pucker", renderEndingPucker, text=3)
        columnEndingPucker.set_cell_data_func(renderEndingPucker, self.__setPuckerBg, (5, 6))
        confTreeView.append_column(columnEndingPucker)
        
        #select the appropriate row of the TreeView
        confTreeSelection = confTreeView.get_selection()
        confTreeSelection.select_iter(confListStore.get_iter((initialConfNum)))
        
        #figure out the height of a row (including the line dividing two rows)
        rowHeight = columnLeadingPucker.cell_get_size()[4] + confTreeView.style_get_property("vertical-separator")
        
        return (confTreeView, confTreeSelection, rowHeight)
    
    
    def __createConfScrolledWindow(self, confTreeView, rowHeight):
        """Create a ScrolledWindow object to display the conformer TreeView
        ARGUMENTS:
            confTreeView        - the TreeView object used to display the list of conformers
            rowHeight           - the height of a row in confTreeView
        RETURNS:
            treeViewScroll      - the newly created ScrolledWindow object
            treeViewScrollWidth - the width of treeViewScroll (used to figure out how wide to make the overview ScrolledWindow)
        """
        
        treeViewScroll = gtk.ScrolledWindow()
        treeViewScroll.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
            #we make sure that the window is big enough that it should never need a horizontal scrollbar, but sometimes GTK
            #insists on adding one, so we disallow it here
        treeViewScroll.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        treeViewScroll.add(confTreeView)
        
        #figure out the size of the treeview + vertical scrollbar is so we can set the scrolled window to the appropriate size
        (tableWidthReq, tableHeightReq) = confTreeView.size_request()
        try:
            vScrollWidthReq = treeViewScroll.get_vscrollbar().size_request()[0]
        except AttributeError:
            #get_vscrollbar only exists in PyGTK 2.8 or above
            #if we're using an older version of PyGTK, then just check the default size of a vertical scrollbar
            vScrollWidthReq = gtk.VScrollbar().size_request()[0]
        scrollbarSpacing = treeViewScroll.style_get_property("scrollbar-spacing")
        #rowHeight = columnLeadingPucker.cell_get_size()[4] + confTreeView.style_get_property("vertical-separator")
        headerHeight = tableHeightReq - len(puckerList)*rowHeight
        treeViewScrollWidth = tableWidthReq + vScrollWidthReq + scrollbarSpacing + 4
        treeViewScroll.set_size_request(treeViewScrollWidth, headerHeight + rowHeight * 3 + 4)
            #the "+ 4"'s are for the shadow, since there doesn't seem to be a GTK property corresponding to shadow size
            
        return (treeViewScroll, treeViewScrollWidth)
    
    
    def __createOverviewTreeView(self, confTreeView, suiteFrame):
        """Create a TreeView object for the suite overview
        
        ARGUMENTS:
            confTreeView      - the TreeView object for the list of conformers
            suiteFrame        - the Frame object containing information about the current suite
        RETURNS:
            overviewTreeView  - the newly created TreeView object
        """
        
        #create the overview list store
        typeList = [str] * len(self.__curPath)
        overviewListStore = gtk.ListStore(*typeList)
        
        resNumList = self.__origCoords.resNumList()[1+self.__skipFirstNuc:-1]
        overviewListStore.append(resNumList)
        overviewListStore.append([conf.replace("&", "&amp;") for conf in self.__curPath])
            #we need to escape the ampersands since we're using Pango markup
        
        #decide how wide to make the columns
        #they should be as wide as the longest residue number or conformer name (conformer names are always 2 characters)
        columnWidthChars = max([len(resNum) for resNum in resNumList] + [2])
        
        overviewTreeView = gtk.TreeView(overviewListStore)
        overviewTreeView.set_enable_search(False)
        overviewTreeView.set_headers_visible(False)
        overviewTreeView.get_selection().set_mode(gtk.SELECTION_NONE)
        
        #set_grid_lines() only exists in Gtk 2.10 or newer
        try:
            overviewTreeView.set_grid_lines(gtk.TREE_VIEW_GRID_LINES_VERTICAL)
        except AttributeError:
            #in Gtk's older than 2.10, there doesn't seem to be a way to set the gridlines to be vertical only
            #without modifying the theme, which we probably don't want to do
            #so we just leave the gridlines out
            #(The "enable-grid-lines" property can only be set to True or False, not vertical-only)
            pass

        
        setColumn = len(self.__curPath) #this column is False for the suite number row and True for the conformer row
        #we use it to selectively apply formatting changes
        
        #create columns and renderers for each suite
        for (colNum, conf) in enumerate(self.__curPath):
            #print "(colNum, conf) = ", colNum, conf
            
            renderer = gtk.CellRendererText()
            renderer.set_property('xalign', 0.5)
            renderer.set_property('width-chars', columnWidthChars+2) #the +2 ensures that we have enough width for bolding
            #renderer.set_property('width-chars', 8) #for debugging
            
            if colNum == 0:
                renderer.set_property('weight', BOLD_WEIGHT)
            

            #figure out the ratio of the best conformer score to the second best conformer score
            (bestProb, secondBestProb) = sorted(self.__predictedProbs[colNum].values(), reverse=True)[0:2]
            try:
                ratio = bestProb / secondBestProb
            except ZeroDivisionError:
                ratio = float('inf')
            
            #retrieve the minimization score for the second nucleotide of the suite
            #for the first suite, use either the first or the second nucleotide, whichever is worse
            if colNum == 0:
                minimizationScore = max(self.__minimizationScores[0], self.__minimizationScores[1])
            else:
                minimizationScore = self.__minimizationScores[colNum+1]
            
            #set the background color based on how sure we are about the prediction or how well the minimization worked (whichever is worse)
            if ratio <= 1.2 or minimizationScore >= 30:
                renderer.set_property('cell-background', 'red')
            elif ratio <= 1.75 or minimizationScore >= 20:
                renderer.set_property('cell-background', 'orange')
            elif ratio <= 2.5 or minimizationScore >= 10:
                renderer.set_property('cell-background', 'yellow')
            #otherwise, leave the background as white
            
            
            #store the column number as the title so we can easily retrieve it later
            column = gtk.TreeViewColumn(str(colNum+1), renderer, markup=colNum)
            column.set_expand(True)
            overviewTreeView.append_column(column)
            
            #make all cells the same height 
            #otherwise the height of the second row will change when the font switches from bold to bold+italics
            if colNum == 0:
                #we want the cells to be the same height that the first column starts as, since that one starts bold
                colHeight = column.cell_get_size()[4] - confTreeView.style_get_property("vertical-separator")
            renderer.set_property("height", colHeight)
        
        overviewTreeView.connect("button-press-event", self.__overviewClick, confTreeView, suiteFrame)
        
        return overviewTreeView
    
    
    def __createOverviewScrolledWindow(self, overviewTreeView, prevConfButton, treeViewScrollWidth):
        """Create a ScrolledWindow object to display the suite overview TreeView
        
        ARGUMENTS:
            overviewTreeView    - the TreeView object used to display the suite overview
            prevConfButton      - the previous conformer button (so we can calculate its width)
            treeViewScrollWidth - the width of the TreeView used to display the conformer list
        RETURNS:
            overviewScroll      - the newly created ScrolledWindow object
            wideOverview        - whether or not the ScrolledWindow takes up the full window width
        """
        
        #If the overview table is wider than the suite frame (i.e. wider than the window is going to be), then we want to
        #limit the overview table width to the size of the suite frame and include a horizontal scrollbar.
        #To know if the overview table is wider or narrower than the suite frame , first we have to determine the width of the suite frame.
        #Since the suite frame isn't realized (i.e. drawn on the screen) yet, we can just call suite_frame.size_request(),
        #so we have to add up the size requests for all of the lowest level objects.
        
        overviewScroll = gtk.ScrolledWindow()
        overviewScroll.set_shadow_type(gtk.SHADOW_ETCHED_IN)
        overviewScroll.add(overviewTreeView)
        
        (overviewTableWidthReq, overviewTableHeightReq) = overviewTreeView.size_request()
        
        #determine the width of the previous conformer button
        prevConfButtonHBox = prevConfButton.get_child()
        #add up the width of the icon and the text
        prevConfButtonWidth = sum([child.size_request()[0] for child in prevConfButtonHBox.get_children()])
        #add on the width of the icon and text padding (the 2* is for left and right padding for each)
        prevConfButtonWidth += (2 * prevConfButtonHBox.query_child_packing(prevConfButtonHBox.get_children()[0])[2] + 2 * prevConfButtonHBox.query_child_packing(prevConfButtonHBox.get_children()[1])[2])
        
        suiteFrameWidth = treeViewScrollWidth + prevConfButtonWidth + BUTTON_SPACING + 3*HBOX_SPACING
        
        if overviewTableWidthReq <= suiteFrameWidth:
            overviewScroll.set_policy(gtk.POLICY_NEVER, gtk.POLICY_NEVER)
            overviewScroll.set_size_request(overviewTableWidthReq + 4, overviewTableHeightReq + 4)
            wideOverview = False
        else:
            
            #figure out how tall the scrollbar is
            try:
                hScrollHeightReq = overviewScroll.get_hscrollbar().size_request()[1]
            except AttributeError:
                #get_hscrollbar only exists in PyGTK 2.8 or above
                #if we're using an older version of PyGTK, then just check the default size of a horizontal scrollbar
                hScrollHeightReq = gtk.HScrollbar().size_request()[1]
            scrollbarSpacing = overviewScroll.style_get_property("scrollbar-spacing")
            
            overviewScroll.set_policy(gtk.POLICY_ALWAYS, gtk.POLICY_NEVER)
            overviewScroll.set_size_request(suiteFrameWidth, hScrollHeightReq + scrollbarSpacing + overviewTableHeightReq + 4)
            wideOverview = True
            
        return (overviewScroll, wideOverview)
        
    def __setPuckerBg(self, column, renderer, confListStore, iter, columnList):
        """Appropriately set the background of a cell listing sugar pucker
        
        ARGUMENTS:
            column          - the TreeViewColumn to be changed
            renderer        - the CellRenderer for the column
            confListStore   - the ListStore containing the displayed data
            iter            - a TreeIter pointing to the specified row of the TreeView
            columnList      - a list of (the number of the column to be changed, the number of the column containing the best conformer)
        RETURNS:
            None
        EFFECTS:
            Sets the background color for the specified cell to:
                DIFFERENT_PUCKER_COLOR if it's different than the pucker for the best path
                BEST_PATH_COLOR is it's the best conformer
                None otherwise
        """
        
        (puckerColumn, bestPathColumn) = columnList
        isDifferentPucker = confListStore.get_value(iter, puckerColumn)
        if isDifferentPucker:
            renderer.set_property('cell-background', DIFFERENT_PUCKER_COLOR)
        else:
            isBestPath = confListStore.get_value(iter, bestPathColumn)
            if isBestPath:
                renderer.set_property('cell-background', BEST_PATH_COLOR)
            else:
                renderer.set_property('cell-background', None)
    
    
    def __prevSuite(self, widget, confTreeView, suiteFrame, overviewTreeView):
        """Select the previous suite
        
        ARGUMENTS:
            widget           - the button used to invoke this function.  Ignored.
            confTreeView     - the TreeView used to display data about conformers for the current suite
            suiteFrame       - the GTK frame containing the "Suite X of Y" label
            overviewTreeView - the TreeView used to display the suite overview
        RETURNS:
            None
        EFFECTS:
            Recenters the screen on the previous suite and updates the RCrane GUI to reflect that suite.
        """
        if self.__curSuite > 1:
            self.__prevOrNextSuite(confTreeView, suiteFrame, overviewTreeView, self.__curSuite, self.__curSuite-1)
            
    
    def __nextSuite(self, widget, confTreeView, suiteFrame, overviewTreeView):
        """Select the next suite
        
        ARGUMENTS:
            widget           - the button used to invoke this function.  Ignored.
            confTreeView     - the TreeView used to display data about conformers for the current suite
            suiteFrame       - the GTK frame containing the "Suite X of Y" label
            overviewTreeView - the TreeView used to display the suite overview
        RETURNS:
            None
        EFFECTS:
            Recenters the screen on the next suite and updates the RCrane GUI to reflect that suite.
        """
        if self.__curSuite < len(self.__predictedProbs):
            self.__prevOrNextSuite(confTreeView, suiteFrame, overviewTreeView, self.__curSuite, self.__curSuite+1)
    
    
    def __prevOrNextSuite(self, confTreeView, suiteFrame, overviewTreeView, oldSuiteIndex, newSuiteIndex):
        """Select the specified suite
        
        ARGUMENTS:
            confTreeView     - the TreeView used to display data about conformers for the current suite
            suiteFrame       - the GTK frame containing the "Suite X of Y" label
            overviewTreeView - the TreeView used to display the suite overview
            oldSuiteIndex    - the suite we are moving from
            newSuiteIndex    - the suite we are moving to
        RETURNS:
            None
        EFFECTS:
            Recenters the screen on newSuite and updates the RCrane GUI to reflect that suite.
        """
        
        #store the newly selected conformer
        self.__curPath[oldSuiteIndex-1] = self.__selectedConf
        
        #if we've switched puckers, store the previous or next conformers
        if len(self.__previousSuiteConf) == 2:
            #print "Previous suite was", self.__curPath[oldSuiteIndex-2]
            self.__curPath[oldSuiteIndex-2] = self.__previousSuiteConf[puckerList[self.__selectedConf][0]]
            #print "Setting previous suite to", self.__curPath[oldSuiteIndex-2]
        if len(self.__nextSuiteConf) == 2:
            #print "Next suite was", self.__curPath[oldSuiteIndex]
            self.__curPath[oldSuiteIndex  ] = self.__nextSuiteConf    [puckerList[self.__selectedConf][1]]
            #print "Setting next suite to", self.__curPath[oldSuiteIndex]
        
        #go to the center of the next suite
        oldSuiteNum = self.__origCoords.nucs[oldSuiteIndex + self.__skipFirstNuc].resNum        
        newSuiteNum = self.__origCoords.nucs[newSuiteIndex + self.__skipFirstNuc].resNum
        newPhos = self.__pseudoMolecule.getPhosCoords(newSuiteNum)
        set_rotation_centre(*newPhos)
        remove_atom_label(self.__pseudoMolecule.molNum(), self.__pseudoMolecule.chain, int(oldSuiteNum), " P  ")
        add_atom_label(self.__pseudoMolecule.molNum(), self.__pseudoMolecule.chain, int(newSuiteNum), " P  ") #label the phosphate
        graphics_draw()
        
        #update the suite label
        if self.__useResNums:
            frameTitle = "Suite %s (%i of %i)" % (str(newSuiteNum), newSuiteIndex, len(self.__predictedProbs))
        else:
            frameTitle = "Suite %i of %i" % (newSuiteIndex, len(self.__predictedProbs))
        suiteFrame.set_label(frameTitle)
        
        (newListStore, initialConfNum) = self.__createConfListStore(newSuiteIndex)
        
        confTreeSelection = confTreeView.get_selection()
        confTreeSelection.disconnect(self.__treeViewChangedHandlerID)
        
        confTreeView.set_model(newListStore)
        confTreeSelection.select_iter(newListStore.get_iter((initialConfNum)))
        confTreeView.scroll_to_cell((initialConfNum), use_align=True, row_align=1.0)
        self.__treeViewChangedHandlerID = confTreeSelection.connect("changed", self.__changed, overviewTreeView)
        
        
        #update the suite overview table
        
        #move the bolding to the new column
        overviewTreeView.get_column(oldSuiteIndex-1).get_cell_renderers()[0].set_property('weight', pango.WEIGHT_NORMAL)
        overviewTreeView.get_column(newSuiteIndex-1).get_cell_renderers()[0].set_property('weight', BOLD_WEIGHT)
        
        #scroll the treeview if necessary
        overviewTreeView.scroll_to_cell((0,), overviewTreeView.get_column(newSuiteIndex-1))
        
        #force a redraw of both rows of the overview stable so that the boldness gets updated properly
        overviewModel = overviewTreeView.get_model()
        overviewModel.row_changed((0,), overviewModel.get_iter((0,)))
        overviewModel.row_changed((1,), overviewModel.get_iter((1,)))
        
        
        #set up the cache variables
        self.__selectedConf = self.__curPath[newSuiteIndex-1]
        self.__curSuite = newSuiteIndex
        self.__clearCache()
    
    def __clearCache(self):
        """Reset all data about the current suite (appropriate for when we've just switched suites)
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            Resets (or initializes) all cache variables
        """
        
        self.__structureCache = {}
        self.__startingPuckerStructureCache = {}
        self.__endingPuckerStructureCache = {}
        self.__previousSuiteConf = {}
        self.__nextSuiteConf = {}
        self.__intermediateAtomLocsStartingPuckerCache = {}
        self.__intermediateAtomLocsEndingPuckerCache = {}
        self.__intermediateAtomLocsCache = {}
        
    
    
    def __closeWin(self, window, widget):
        """Close the RCrane window
        
        ARUGMENTS:
            window - the window to close
            widget - the widget used to close the window
        NOTES:
            This function simply calls __accept with the appropriate argument list
        """
        
        #self.__close(widget)
        self.__accept(widget)
    
    def __close(self, widget = None):
        """Close the RCrane window
        
        ARUGMENTS:
            widget - the widget used to close the window
        RETURNS:
            None
        EFFECTS:
            Closes the RCrane window
        """
        
        self.__pseudoMolecule.clearExtraBonds()
        if self.__oldCoordsObj is not None:
            close_generic_object(self.__oldCoordsObj)
        
        graphics_draw()
        self.__window.destroy()
    
    
    def __accept(self, widget = None):
        """Accept the new molecule and close the RCrane window
        
        ARUGMENTS:
            widget - the widget used to close the window
        RETURNS:
            None
        NOTE:
            This function re-orders atoms in the pseudomolecule (that's the "accept" part) and then calls __close
        """
        
        #perform the final molecule cleanup (reorder atoms and, if we just did a rotamerize, fix segids)
        self.__pseudoMolecule.finalMoleculeCleanup(fixSegids = bool(self.__pseudoMolecule.hasSavedCoordinates()))
        self.__close(widget)
    
    
    def __toggleOldCoords(self, widget):
        """Toggle display of the old coordinates object
        
        ARGUMENTS:
            widget - the check button that was toggled
        RETURNS:
            None
        """
        
        set_display_generic_object(self.__oldCoordsObj, int(widget.get_active()))
    
    
    def __restoreOldCoords(self, widget):
        """Set the coordinates back to the way they were before we called rotamerize
        
        ARGUMENTS:
            widget - the cancel button that was clicked to call this function.  Ignored
        RETURNS:
            None
        """
        
        self.__pseudoMolecule.restoreSavedCoordinates()
        self.__close()
    
    
    def __overviewClick(self, overviewTreeView, event, confTreeView, suiteFrame):
        """Respond to the user clicking on the overview TreeView object
        
        ARGUMENTS:
            overviewTreeView - the TreeView object for the suite overview
            event            - the click event
            confTreeView     - the TreeView object for the list of conformers
            suiteFrame       - the Frame object containing information about the current suite
        RETURNS:
            None
        """
        
        #figure out which column was just clicked on
        clickedColumn = overviewTreeView.get_path_at_pos(int(event.x), int(event.y))[1]
        
        #The suite index for each column is stored as the column title
        #(since there seems to be no easy way to get that information from the column object itself)
        newSuite = int(clickedColumn.get_title())
        
        
        if newSuite == self.__curSuite:
            #if the user clicked on the current suite, then just re-center
            newSuiteNum = self.__origCoords.nucs[newSuite + self.__skipFirstNuc].resNum
            newPhos = self.__pseudoMolecule.getPhosCoords(newSuiteNum)
            set_rotation_centre(*newPhos)
        else:
            #if the user clicked on a new suite
            #print "Changing suites in __overviewClick"
            self.__prevOrNextSuite(confTreeView, suiteFrame, overviewTreeView, self.__curSuite, newSuite)
    
    
    def __changeOverviewConformer(self, overviewTreeView, suiteNum, newConf):
        """Change the conformer listed in the overview table for the specified suite
        
        ARGUMENTS:
            overviewTreeView - the TreeView object for the suite overview
            suiteNum         - the number of the suite to change
            newConf          - the new conformer to list in the overview table
        RETURNS:
            None
        """
        
        suiteNum -= 1 #we need the suite numbers to start at 0, not 1
        
        overviewListStore = overviewTreeView.get_model()
        if newConf == self.__bestPath[suiteNum]:
            overviewListStore.set_value(overviewListStore.get_iter((1,)), suiteNum, newConf.replace("&", "&amp;"))
        else:
            overviewListStore.set_value(overviewListStore.get_iter((1,)), suiteNum, "<i>" + newConf.replace("&", "&amp;") + "</i>")
        