#!/usr/bin/env python
"""Recalculate coordinates for an already-built stretch of nucleotides."""

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

from copy import deepcopy
import gtk
import os.path
from coot import add_status_bar_text, new_generic_object_number, set_display_generic_object, to_generic_object_add_line, graphics_draw, clear_pending_picks, imol_refinement_map
from coot import user_defined_click_py as user_defined_click

from pseudoMolecule import PseudoMolecule, PseudoMoleculeError
from rotamerSeq import determineRotamerSeq
from calcCoords import recalcCoords
from guiUtils import createRCraneWindowObject, createRCraneErrorDialog
from reviewSuitesGui import ReviewSuitesGui
from bondLength import BOND_LIST_FULL
from selectMapDialog import selectMapDialog
from citationPopup import createCitationPopup
import extendChain #we only need clearPendingExtendChain, but if we import that function directly into this namespace
    #we'll wind up with circular dependencies
from reportInputErrors import reportPDB2Error, reportInsCodeError, STANDARD_BASES, reportModifiedNucError

#from pprint import pprint

ORIG_COORDS_COLOR = "yellowtint" #the color to use when drawing the original coordinates


#These two variables allow the user to cancel a rotamerize before selecting atoms
#We could have turned this module into a class and made these class variables, but we'll never, ever
#need more than one instance of these variables per execution of Coot, so module-level variables
#work just as well

__menuItem = None #the "Rotamerize existing structure..." menuItem object
__menuItemNoDensity = None #the "Rotamerize without density..." menuItem object
    #these get set in launch.py via a call to storeMenuItem()

waitingForClicks = False #whether or not we are waiting for the user to pick two atoms

MAX_NUCLEOTIDES_TO_ROTAMERIZE = 20
    #the maximum number of nucleotides to rotamerize at once.  This can be changed via setMaxNucleotdiesToRotamerize()
    #if this value is <= 0, then there is no limit


def pickAtoms(ignoreDensity = False):
    """Prompt the user to pick two atoms to rotamerize a stretch of nucleotides.  Alternatively, if
    we're already waiting for the user to pick atoms, then cancel the pending picks.
    
    ARGUMENTS:
        None
    OPTIONAL ARGUMENTS:
        ignoreDensity - ignore the density when performing the minimization
                        defaults to False
    RETURNS:
        None
    """
    
    #create the citation pop-up if necessary
    createCitationPopup()
    
    global waitingForClicks
    
    #make sure that we're not waiting on a pending extendChain call
    extendChain.clearPendingExtendChain()
    
    if waitingForClicks:
        #if we're already in the middle of a rotamerize call, then cancel the pending rotamerize
        print "Rotamerize cancelled"
        add_status_bar_text("Rotamerize cancelled")
        clear_pending_picks() #tell Coot to stop waiting for the user to click on atoms
        __setMenuToRotamerize()
        waitingForClicks = False
    
    else:
        #if we're not in the middle of a rotamerize call, then start one
        
        #make sure that there is a refinement map set
        if not ignoreDensity and imol_refinement_map() == -1:
            print "No refinement map set for RCrane rotamerize"
            selectMapDialog(pickAtoms)
            return
        
        print "Click on 2 atoms (in the same molecule)"
        add_status_bar_text("Pick 2 atoms [Ctrl Left-mouse rotates the view]...")
        __setMenuToCancel()
        waitingForClicks = True
        user_defined_click(2, lambda atomSpec1, atomSpec2: rotamerizeFromAtomSpecs(atomSpec1, atomSpec2, ignoreDensity))


def rotamerizeFromAtomSpecs(atomSpec1, atomSpec2, ignoreDensity = False):
    """Rotamerize a stretch of nucleotides
    
    ARGUMENTS:
        atomSpec1 - an atom spec specifying the start of the stretch
        atomSpec2 - an atom spec specifying the end of the stretch
    OPTIONAL ARGUMENTS:
        ignoreDensity - ignore the density when performing the minimization
                        defaults to False
    RETURNS:
        False if there is a problem with the selected nucleotides
        None otherwise
    """
    
    #stop waiting for the user to cancel the rotamerize
    __setMenuToRotamerize()
    global waitingForClicks
    waitingForClicks = False
    
    #print "atom 1 = ", atomSpec1
    #print "atom 2 = ", atomSpec2
    
    #unpack the atom descriptions
    (mol1, chain1, startResNum, startInsCode, startAtom, alt1) = atomSpec1[1:]
    (mol2, chain2, endResNum,   endInsCode,   endAtom,   alt2) = atomSpec2[1:]
    
    #make sure that the picked atoms define a single stretch of RNA
    if mol1 != mol2:
        print "Rotamerize atoms are not in the same molecule"
        add_status_bar_text("Rotamerize atoms are not in the same molecule")
        return False
    
    if chain1 != chain2:
        print "Rotamerize atoms are not in the same chain"
        add_status_bar_text("Rotamerize atoms are not in the same chain")
        return False
    
    molNum       = mol1
    chain        = chain1
    startResFull = str(startResNum) + str(startInsCode)
    endResFull   = str(endResNum)   + str(endInsCode)
    
    try:
        pseudoMol = PseudoMolecule(createFromMolecule = molNum, chain = chain, batons = False)
    except PseudoMoleculeError, err:
        #if creating the PseudoMolecule raises an error, then it means that molecule we clicked on contains
        #anisotropic temperature records and we're using a version of Coot that's older than rev 3631
        #so clear_and_update_molecule can't handle anisotropic records
        
        print err
        #add_status_bar_text(str(err))
        
        #Since this isn't an intuitive problem (i.e. it's far less obvious for the user than "you clicked on
        #two different chains"), we create a dialog box instead of a status bar message
        createRCraneErrorDialog(str(err) + "\n\nRotamerization canceled.")
        
        return False
    
    startIndex = pseudoMol.resIndex(startResFull)
    endIndex   = pseudoMol.resIndex(endResFull)
    
    if (startIndex > endIndex):
        #make sure that the starting atom is before the ending atom
        #if it isn't, swap them
        (startResFull, endResFull) = (endResFull, startResFull)
        (startIndex,   endIndex)   = (endIndex,   startIndex)
        (startAtom,    endAtom)    = (startAtom,  endAtom)
    
    #make sure that the molecule using PDB3 naming (as oppsed to PDB2)
    #TODO: check all residues for PDB3 naming instead of just the first one in case the molecule uses a mix of PDB2 and PDB3 naming
    #   this seems like an unlikely problem, but it can't hurt to double check
    if pseudoMol.checkPDB2FromIndex(startIndex):
        reportPDB2Error(startResFull)
        return False
    
    if (startIndex > 0 
      and startAtom.strip() in frozenset(["P", "OP1", "OP2", "O5'", "C5'"])
      and pseudoMol.connectedToPrevFromIndex(startIndex)):
        #if we're not at the first nucleotide and the user clicked on an atom belonging exclusively to the
        #starting suite of this residue, then start the minimization at the previous nucleotide
        
        startIndex -= 1
        startResFull = pseudoMol.resNumFull(startIndex)
    
    if ((endIndex + 1) < pseudoMol.getNumNts()
      and endAtom.strip() == "O3'"
      and pseudoMol.connectedToNextFromIndex(endIndex)):
        #if we're not at the last nucleotide and the user clicked on an atom belonging exclusively to the
        #ending suite of this residue, then include the next residue in the minimization
        endIndex += 1
        endResFull = pseudoMol.resNumFull(endIndex)
    
    #warn the user if they haven't selected an entire suite
    if startIndex == endIndex:
        print "You must select at least one suite (two residues) to rotamerize."
        add_status_bar_text("You must select at least one suite (two residues) to rotamerize.")
        return False
    
    #make sure that the selected residues are all connected
    for curIndex in range(startIndex, endIndex):
        if not pseudoMol.connectedToNextFromIndex(curIndex):
            print "Nucleotides " + pseudoMol.resNumFull(curIndex) + " and " + pseudoMol.resNumFull(curIndex+1) + " are not connected.  Cannot rotamerize."
            add_status_bar_text("Nucleotides " + pseudoMol.resNumFull(curIndex) + " and " + pseudoMol.resNumFull(curIndex+1) + " are not connected.  Cannot rotamerize.")
            return False
    #note that this loop does not check to see if the last nucleotide is connectedToNext, since that's not required
    
    
    for curIndex in range(startIndex, endIndex+1):
        #make sure that the selected residues aren't modified nucleotides
        curResType = pseudoMol.resTypeFromIndex(curIndex)
        if curResType not in STANDARD_BASES:
            reportModifiedNucError(pseudoMol.resNumFull(curIndex), curResType)
            return False
        
        #make sure that all selected residues have phosphates and glycosidic bonds
        missingAtom = pseudoMol.checkPhosAndGlycosidicFromIndex(curIndex)
        if missingAtom is not None:
            print "Nucleotide " + pseudoMol.resNumFull(curIndex) + " is missing " + missingAtom + " atom.  Cannot rotamerize without phosphate and glycosidic bond coordinates."
            add_status_bar_text("Nucleotide " + pseudoMol.resNumFull(curIndex) + " is missing " + missingAtom + " atom.  Cannot rotamerize without phosphate and glycosidic bond coordinates.")
            return False
    
    #if the ending res is the last residue of the chain, then change it to the second to last residue
    #since we need a terminal 3' phosphate
    #(Note that createPartialChainObject will then add the extra 3' phosphate)
    if (endIndex + 1) >= pseudoMol.getNumNts() or not pseudoMol.connectedToNextFromIndex(endIndex):
        endIndex -= 1
        endResFull = pseudoMol.resNumFull(endIndex)
        
        #make sure that we're not now down to a single residue
        if startIndex >= endIndex:
            print "Cannot rotamerize this suite without a 3' phosphate."
            add_status_bar_text("Cannot rotamerize this suite without a 3' phosphate.")
            return False
    
    
    #make sure that there are no insertion codes (since Coot's refine_zone can't handle these)
    #ideally, refine_zone will eventually get fixed and I can remove this check
    for curIndex in range(startIndex, endIndex+2):
        #we need to check the last residue and the one after that, since createPartialChainObject will add on the 3' phosphate
        if pseudoMol.resNum(curIndex)[1] != "":
            reportInsCodeError(pseudoMol.resNumFull(curIndex))
            return False
    
    #make sure that the user hasn't tried to rotamerize too many nucleotides
    #If they have, it's likely due to misclicking, so we don't want to force them to sit through rotamerizing hundreds of nucleotides
    numNucsToRotamerize = endIndex - startIndex + 1
    if (MAX_NUCLEOTIDES_TO_ROTAMERIZE > 0) and (numNucsToRotamerize > MAX_NUCLEOTIDES_TO_ROTAMERIZE):
        
        print "Warning:  Too many nucleotides to rotamerize."
        print "   Trying to rotamerize", numNucsToRotamerize, "nucleotides and your current maximum is", MAX_NUCLEOTIDES_TO_ROTAMERIZE
        print "   To increase the maximum, go to Calculate -> Scripting -> Python"
        print "   and type rcrane.setRotamerizeMaxNucleotides(" + str(numNucsToRotamerize) + ")"
        print "   or type rcrane.setRotamerizeMaxNucleotides(-1) to remove the limit."
        add_status_bar_text("Warning:  Too many nts to rotamerize.  Trying to rotamerize " + str(numNucsToRotamerize) + " nts and your current max is " + str(MAX_NUCLEOTIDES_TO_ROTAMERIZE))
        return False
    
    
    #################################################################################################
    #                                                                                               #
    #  now that we've finished sanity checking the input, we can move on to actually rotamerizing   #
    #                                                                                               #
    #################################################################################################
    
    #run the conformer predictions
    (chainObj, resNum5p, resNum3p) = pseudoMol.createPartialChainObject(startResFull, endResFull, True)
    (bestPath, predictedProbs) = determineRotamerSeq(chainObj)
    print "Suite string: " + "".join(bestPath)
    
    prevTorsions = None
    if resNum5p is None:
        #if we're at the start of the chain, add in a None rotamer as the previous suite
        rebuildRots = [None]+bestPath
    else:
        #if we're not at the start of the chain, record the current 5' torsions so they can be used as restraints in minimization
        prevTorsions = pseudoMol.calcSuiteTorsions(startResFull)
        rebuildRots = [prevTorsions]+bestPath
    
    #if we're not at the end of the chain, record the current 3' torsions so that they can be used as restraints in minimization
    nextTorsions = None
    if resNum3p is None:
        rebuildRots = rebuildRots+[None]
    else:
        nextTorsions = pseudoMol.calcSuiteTorsions(resNum3p)
        rebuildRots = rebuildRots+[nextTorsions]
    
    #pprint(rebuildRots)
    #from time import sleep; sleep(3)
    
    #store the B-factors
    pseudoMol.saveBfacs(startResFull, endResFull)
    #pprint(bfacs)
    
    #store a backup of the nucleotides to be rebuilt so we can restore them later if needed
    pseudoMol.saveCoordinates(resNum5p or startResFull, resNum3p or endResFull)
    
    #display the original coordinates as a drawing object
    oldCoordsObj = __drawOrigCoords(chainObj)
    #Note that this is the opposite of Coot's minimization behavior.  When Coot does a minimization, the original
    #coordinates are displayed normally, and the new coordinates are displayed with white carbons and faded colors.
    #Here, the new coordinates will be displayed normally and old coordinates are displayed in yellow.
    #This was done for several reasons:
    #  There is currently no way to display a molecule using the "faded" coloring of the minimization results.  This
    #    would be required to really mimic Coot's behaviour.  I could set the color rotation to try to get white
    #    carbons, but it still looks very different.
    #  With the minimization results, we only have one choice: accept or reject (or try dragging atoms around).  With
    #    the rotamerization results, we can go through a try out alternate rotamers.  As a result, the user is likely to
    #    spend more time examining rotamerization results than they would when examining minimization results.  Because
    #    of this, we want the new coordinates to be displayed in a way that's easy to see and evaluate.
    
    #create a window object
    window = createRCraneWindowObject()
    
    #start coordinate calculation
    pseudoMol.resetNucs(startResFull, endResFull)
    pseudoMol.setExtraBondRange(resNum5p or startResFull, resNum3p or endResFull)
    (intermediateAtomLocs, minimizationScores) = recalcCoords(startResFull, endResFull, rebuildRots, chainObj, pseudoMol, window, ignoreDensity)
    
    #if we've included 3' atoms, then intermediateAtomLocs will be one residue off from chainObj,
    #so add a None to the start of intermediateAtomLocs to sync them up
    if resNum3p is not None:
        intermediateAtomLocs = [None] + intermediateAtomLocs
    
    #reset the B-factors
    pseudoMol.restoreSavedBfacs()
    
    #display review GUI
    ReviewSuitesGui(predictedProbs, bestPath, pseudoMol, chainObj, intermediateAtomLocs, minimizationScores, window, True, oldCoordsObj, ignoreDensity = ignoreDensity, prevTorsions = prevTorsions, nextTorsions = nextTorsions)
    

def __drawOrigCoords(origCoords):
    """Display the original pre-rotamerize coordinates as a generic object
    
    ARGUMENTS:
        origCoords - a chain object containing the atoms to be drawn
    RETURNS:
        drawNum - the generic object number
    EFFECTS:
        creates and displays a new Coot generic object
    """
    
    #create a new generic display object
    drawNum = new_generic_object_number("Pre-rotamerize coordinates")
    set_display_generic_object(drawNum, 1)
    
    lineList = [] #a list of all the bonds to be drawn
    
    #draw bonds for each nucleotide
    prevO3 = None
    for curNuc in origCoords.nucs:
        
        #draw each bond if both atoms are present
        for (atom1, atom2) in (BOND_LIST_FULL["backbone"] + BOND_LIST_FULL[curNuc.type]):
            if curNuc.hasAtom(atom1) and curNuc.hasAtom(atom2):
                atom1Coords = curNuc.atoms[atom1]
                atom2Coords = curNuc.atoms[atom2]
                lineList.append(atom1Coords + atom2Coords)
        
        #draw a bond between O5' and the previous phosphate if both atoms are present
        if prevO3 is not None and curNuc.hasAtom("P"):
            lineList.append(prevO3 + curNuc.atoms["P"])
        
        #update prevPhos
        if curNuc.hasAtom("O3'"):
            prevO3 = curNuc.atoms["O3'"]
        else:
            prevO3 = None
    
    #actually draw the bonds
    for curLine in lineList:
        to_generic_object_add_line(drawNum, ORIG_COORDS_COLOR, 6, *curLine)
    graphics_draw()
    
    return drawNum


def storeMenuItem(menuItem, resetLabel = False, ignoreDensity = False):
    """Store the "Rotamerize existing structure..." and "Rotamerize without density..." menu entries so that we can later change them to "Cancel current rotamerize..."
    
    ARGUMENTS:
        menuItem   - the gtk MenuItem to store
    OPTIONAL ARGUMENTS:
        resetLabel    - whether we should reset the menu label to "Rotamerize existing structure..."
                        (and also clear all pending atom picks in Coot)
                        This is intended for use when we're reloading RCrane
                        Defaults to False
        ignoreDensity - whether this is the menu entry for "Rotamerize without density..."
    RETURNS:
        None
    EFFECTS:
        stores menuItem in the module-level variable __menuItem
    """
    
    if ignoreDensity:
        global __menuItemNoDensity
        __menuItemNoDensity = menuItem
    else:
        global __menuItem
        __menuItem = menuItem
    
    #if we're being called during an RCrane reload, then reset everything back to the appropriate starting position
    global waitingForClicks
    if resetLabel:
        waitingForClicks = False #this is redundant with module initialization, but can't hurt to do
                                 #in case this function is called in an unexpected way
        __setMenuToRotamerize()
        clear_pending_picks()


def __setMenuToRotamerize():
    """Set the menu entry labels to "Rotamerize existing structure..." and "Rotamerize without density..."
    
    ARGUMENTS:
        None
    RETURNS:
        None
    """
    
    if __menuItem is not None:
        __menuItem.get_child().set_text("Rotamerize existing structure...")
        
    if __menuItemNoDensity is not None:
        __menuItemNoDensity.get_child().set_text("Rotamerize without density...")


def __setMenuToCancel():
    """Set the menu entry label to "Cancel current rotamerize..."
    
    ARGUMENTS:
        None
    RETURNS:
        None
    """
    
    #print "__menuItem =", __menuItem
    if __menuItem is not None:
        __menuItem.get_child().set_text("Cancel current rotamerize...")
        
    if __menuItemNoDensity is not None:
        __menuItemNoDensity.get_child().set_text("Cancel current rotamerize...")


def setRotamerizeMaxNucleotides(newVal):
    """Set the maximum number of nucleotides that can be rotamerized at once
    
    ARGUMENTS:
        newVal - the new maximum number of nucleotides to rotamerize at once
                 if this value is <= 0, then there will be no limit
    RETURNS:
        None
    """
    
    global MAX_NUCLEOTIDES_TO_ROTAMERIZE
    MAX_NUCLEOTIDES_TO_ROTAMERIZE = int(newVal)


def getRotamerizeMaxNucleotides():
    """Get the maximum number of nucleotides that can be rotamerized at once
    
    ARGUMENTS:
        None
    RETURNS
        MAX_NUCLEOTIDES_TO_ROTAMERIZE - the maximum number of nucleotides that can be rotamerized at once
    """
    
    return MAX_NUCLEOTIDES_TO_ROTAMERIZE


def clearPendingRotamerize():
    """Cancel any pending rotamerize call (intended to be called by extendChain)
    
    ARGUMENTS:
        None
    RETURNS:
        True if there was a pending rotamerize call
        False otherwise
    """
    
    global waitingForClicks
    
    if waitingForClicks:
        #if we're already in the middle of a rotamerize call, then cancel the pending rotamerize
        print "Rotamerize cancelled"
        clear_pending_picks() #tell Coot to stop waiting for the user to click on atoms
        __setMenuToRotamerize()
        waitingForClicks = False
        return True
    else:
        return False
