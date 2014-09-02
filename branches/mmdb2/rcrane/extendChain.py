#!/usr/bin/env python
"""Extend an existing chain of nucleotides."""

# Copyright 2012 Kevin Keating
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

from coot import add_status_bar_text, clear_pending_picks, imol_refinement_map
from coot import user_defined_click_py as user_defined_click

from selectMapDialog import selectMapDialog
from citationPopup import createCitationPopup
from pseudoMolecule import PseudoMolecule
import rotamerize #we only need clearPendingRotamerize, but if we import that function directly into this namespace
    #we'll wind up with circular dependencies
from traceGui import TraceGui
from reportInputErrors import reportPDB2Error, reportInsCodeError, STANDARD_BASES, reportModifiedNucError


__menuItem = None #the "Add to terminal nts..." menuItem object
waitingForClicks = False #whether or not we are waiting for the user to pick an atom


def pickAtoms():
    """Prompt the user to pick an atoms to extend a stretch of nucleotides.  Alternatively, if
    we're already waiting for the user to pick an atom, then cancel the pending pick.
    
    ARGUMENTS:
        None
    RETURNS:
        None
    """
    
    #create the citation pop-up if necessary
    createCitationPopup()
    
    #make sure that we're not waiting on a pending rotamerize call
    rotamerize.clearPendingRotamerize()
    
    global waitingForClicks
    
    if waitingForClicks:
        #if we're already in the middle of an extend chain call, then cancel the pending call
        print "Rotamerize cancelled"
        add_status_bar_text("Extend chain cancelled")
        clear_pending_picks() #tell Coot to stop waiting for the user to click on atoms
        __setMenuToExtendChain()
        waitingForClicks = False
    
    else:
        #if we're not in the middle of a extend chain call, then start one
        
        #make sure that there is a refinement map set
        if imol_refinement_map() == -1:
            print "No refinement map set for RCrane extend chain"
            selectMapDialog(pickAtoms)
            return
        
        print "Click on a nucleotide to extend"
        add_status_bar_text("Pick a nucleotide [Ctrl Left-mouse rotates the view]...")
        __setMenuToCancel()
        waitingForClicks = True
        user_defined_click(1, extendChainFromAtomSpec)
    

def extendChainFromAtomSpec(atomSpec):
    """Rotamerize a stretch of nucleotides
    
    ARGUMENTS:
        atomSpec - an atom spec specifying the start of the chain extension
    RETURNS:
        False if there is a problem with the selected nucleotide
        None otherwise
    """
    
    #stop waiting for the user to cancel the rotamerize
    __setMenuToExtendChain()
    global waitingForClicks
    waitingForClicks = False
    
    #print "In extendChainFromAtomSpec"
    #print "atomSpec =", atomSpec
    
    (mol, chain, resNum, insCode, atom, alt) = atomSpec[1:]
    atom = atom.strip()
    try:
        pseudoMol = PseudoMolecule(createFromMolecule = mol, chain = chain)
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
    
    resNumFull = str(resNum) + str(insCode)
    resIndex = pseudoMol.resIndex(resNum)
    
    #figure out what direction to extend in
    connectedToNext = pseudoMol.connectedToNextFromIndex(resIndex)
    connectedToPrev = pseudoMol.connectedToPrevFromIndex(resIndex)
    curResIsPhosOnly = pseudoMol.isOnlyPhosGroupFromIndex(resIndex)
    
    nextResIsPhosOnly = False
    if connectedToNext:
        nextResIsPhosOnly = pseudoMol.isOnlyPhosGroupFromIndex(resIndex+1)
    
    
    extendDir = None #the direction to extend the chain in
    
    if curResIsPhosOnly:
        if connectedToPrev:
            #we're building on a 3' phosphate
            extendDir = 3
        else:
            #the selected residue is an isolated phosphate group, so we can't do anything
            print "Nucleotide " + resNumFull + " is an isolated phosphate group.  Cannot extend."
            add_status_bar_text("Nucleotide " + resNumFull + " is an isolated phosphate group.  Cannot extend.")
            return False
    elif nextResIsPhosOnly:
        if connectedToPrev:
            #we're at the 3' end of the chain and there is a 3' phosphate
            extendDir = 3
            resIndex += 1 #increment the resIndex so that it points to the real last residue of the segment
            curResIsPhosOnly = True
            resNumFull = pseudoMol.resNumFull(resIndex) #update the residue name in case we need to use it for later error messages
            
        else:
            #we're dealing with a single nt containing a 3' phosphate
            #we need to decide direction based on the atom clicked
            
            extendDir = __decideDirectionFromAtom(atom, resNumFull)
            if not extendDir: return False #if we couldn't determine a direction from the atom, then cancel
    
    elif connectedToNext and connectedToPrev:
        print "Nucleotide " + resNumFull + " is not at the end of a chain"
        add_status_bar_text("Nucleotide " + resNumFull + " is not at the end of a chain")
        return False
    elif connectedToNext:
        extendDir = 5
    elif connectedToPrev:
        extendDir = 3
    else:
        #we're dealing with a single nt not containing a 3' phosphate
        #we need to decide direction based on the atom clicked
        extendDir = __decideDirectionFromAtom(atom, resNumFull)
        if not extendDir: return False #if we couldn't determine a direction from the atom, then cancel
    
    
    
    #make sure that the residue we're going to extend doesn't have an insertion code
    #   (also check the next and previous residues)
    if pseudoMol.resNum(resIndex)[1] != "":
        reportInsCodeError(resNumFull)
        return False
    if connectedToNext and pseudoMol.resNum(resIndex+1)[1] != "":
        reportInsCodeError(pseudoMol.resNumFull(resIndex+1))
        return False
    if connectedToPrev and pseudoMol.resNum(resIndex-1)[1] != "":
        reportInsCodeError(pseudoMol.resNumFull(resIndex-1))
        return False
    
    #figure out which residue we should be checking for PDB2 naming and missing glycosidic bond atoms
    if curResIsPhosOnly:
        resIndexToCheck = resIndex - 1
        resNumFullToCheck = pseudoMol.resNumFull(resIndex - 1)
    else:
        resIndexToCheck = resIndex
        resNumFullToCheck = resNumFull
    
    #make sure that this nucleotide isn't mdofied
    resType = pseudoMol.resTypeFromIndex(resIndexToCheck)
    if resType not in STANDARD_BASES:
        reportModifiedNucError(resNumFullToCheck, resType)
        return False
    
    #make sure that the molecule (or at least the end of it we're going to be building on) uses PDB3 atom naming
    if pseudoMol.checkPDB2FromIndex(resIndexToCheck):
        reportPDB2Error(resNumFullToCheck)
        return False
    
    #the current residue must have a glycosidic bond
    baseType = pseudoMol.resTypeFromIndex(resIndexToCheck)
    if baseType == "A" or baseType== "G":
        glyN = "N9"
    else:
        #we already know that the base isn't modified because we checked above
        glyN = "N1"
    
    if pseudoMol.getAtomCoordsFromIndex("C1'", resIndexToCheck) is None:
        print "Nucleotide " + resNumFullToCheck + " does not have a C1' atom.  Cannot extend chain."
        add_status_bar_text("Nucleotide " + resNumFullToCheck + " does not have a C1' atom.  Cannot extend chain.")
        return False
    if pseudoMol.getAtomCoordsFromIndex(glyN, resIndexToCheck) is None:
        print "Nucleotide " + resNumFullToCheck + " does not have an " + glyN + " atom.  Cannot extend chain."
        add_status_bar_text("Nucleotide " + resNumFullToCheck + " does not have an " + glyN + " atom.  Cannot extend chain.")
        return False
    
    #make sure that this isn't an isolated nt without any phosphates
    curPhosCoords = pseudoMol.getPhosCoordsFromIndex(resIndex) #only used for sanity checking
    if not connectedToNext and not connectedToPrev and curPhosCoords is None:
        print "Nucleotide " + resNumFull + "is an isolated nucleotides without any phosphates.  Cannot extend chain."
        add_status_bar_text("Nucleotide " + resNumFull + "is an isolated nucleotides without any phosphates.  Cannot extend chain.")
        return False
    
    #make sure that we're not building 3'->5' on a nucleotide without a 3' phosphate (because we can't determine pseudotorsions w/o a 3' phos)
    if extendDir == 5 and not connectedToNext:
        print "Cannot extend nucleotide " + resNumFull + " in the 5' direction without a 3' phosphate."
        print "  Please add a 3' phosphate by extending the chain in the 3' direction."
        add_status_bar_text("Cannot extend nucleotide " + resNumFull + " in the 5' direction without a 3' phosphate.")
        return False
    
    print "About to extend the chain in the " + str(extendDir) + "' direction"
    
    #tell the PseudoMolecule object where to insert new residues
    pseudoMol.setResInsertionPoint(resIndex)
    
    #save the molecule state so we can restore it if we cancel the extend
    pseudoMol.saveMoleculeState()
    
    TraceGui(direction = extendDir, existingMolecule = pseudoMol, resIndexToExtend = resIndex)
    


ATOMS3P = frozenset(["C3'", "O3'"])
    #atoms that count as a 3' selection
ATOMS5P = frozenset(["P", "OP1", "OP2", "OP3", "O5'", "C5'", "C4'"])
    #atoms that count as a 5' selection
def __decideDirectionFromAtom(atomName, resNum):
    """Determine the direction to extend the chain based on the atom name (for use with isolated nuclotides)
    
    ARGUMENTS:
        atomName - the atom name
        resNum   - the residue number of the current nucleotide (used in reporting errors)
    RETURNS:
        False if the direction cannot be determined from the atom name
        3 or 5 otherwise
    EFFECTS:
        if the direction cannot be determined from the atom name, an error will be reported to the user
        (and False will be returned)
    """
    
    if atomName in ATOMS3P:
        return 3
    elif atomName in ATOMS5P:
        return 5
    else:
        print "Nucleotide " + resNum + " is an isolated nucleotide.  Please select either the O5' or"
        print "O3' atom to extend in the 5' or 3' direction, respectively."
        add_status_bar_text("Must select either O5' or O3' of isolated nucleotides.")
        return False


def storeMenuItem(menuItem, resetLabel = False):
    """ARGUMENTS:
        menuItem   - the gtk MenuItem to store
    OPTIONAL ARGUMENTS:
        resetLabel    - whether we should reset the menu label to "Extend chain..."
                        (and also clear all pending atom picks in Coot)
                        This is intended for use when we're reloading RCrane
                        Defaults to False
    RETURNS:
        None
    EFFECTS:
        stores menuItem in the module-level variable __menuItem
    """
    
    global __menuItem
    __menuItem = menuItem
    
    #if we're being called during an RCrane reload, then reset everything back to the appropriate starting position
    global waitingForClicks
    if resetLabel:
        waitingForClicks = False #this is redundant with module initialization, but can't hurt to do
                                 #in case this function is called in an unexpected way
        __setMenuToExtendChain()
        clear_pending_picks()


def __setMenuToExtendChain():
    """Set the menu entry labels to "Extend chain..."
    
    ARGUMENTS:
        None
    RETURNS:
        None
    """
    
    if __menuItem is not None:
        __menuItem.get_child().set_text("Extend chain...")


def __setMenuToCancel():
    """Set the menu entry label to "Cancel extend chain..."
    
    ARGUMENTS:
        None
    RETURNS:
        None
    """
    
    if __menuItem is not None:
        __menuItem.get_child().set_text("Cancel extend chain...")


def clearPendingExtendChain():
    """Cancel any pending extend chain call (intended to be called by rotamerize)
    
    ARGUMENTS:
        None
    RETURNS:
        True if there was a pending extend chain call
        False otherwise
    """
    
    global waitingForClicks
    
    if waitingForClicks:
        #if we're already in the middle of a extend chain call, then cancel the pending extend chain
        print "Extend chain cancelled"
        clear_pending_picks() #tell Coot to stop waiting for the user to click on atoms
        __setMenuToExtendChain()
        waitingForClicks = False
        return True
    else:
        return False