#!/usr/bin/env python
"""The pseudoMolecule class for storing atomic coordinate data and updating Coot as appropriate."""

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

from copy import deepcopy
#from pprint import pprint

from coot import new_generic_object_number, set_display_generic_object, to_generic_object_add_line, graphics_draw, close_generic_object, close_molecule, get_molecule_bonds_colour_map_rotation, set_molecule_bonds_colour_map_rotation
#these renames are normally done in coot_utils.py, but that has to be exec'ed, not imported
#so we manually rename them here
from coot import add_molecule_py as add_molecule
from coot import clear_and_update_molecule_py as clear_and_update_molecule
from coot import residue_info_py as residue_info
from coot import svn_revision #this can be removed if we ever require Coot rev 3631 or newer to run RCrane

#use the new version of python_representation if it exists
try:
    from coot_utils_adapter import python_representation_kk as python_representation
    #print "Using python_representation_kk"
except ImportError:
    from coot_utils_adapter import python_representation

from chain import Chain
from nucleotide import Nucleotide
from strucCalc import dist
from bondLength import BOND_LIST, BOND_DIST_CUTOFF

BATONCOLOR = "orange"
MOLECULE_NAME = "RCrane-molecule-%i"
DEFAULT_B_FACTOR = 20.00

ATOMS_TO_ERASE_FOR_REBUILDING_MOLECULE = frozenset("C2' O2' C3' O3' C4' O4' C5' O5' OP1 OP2".split(" "))
BASE_ATOMS = frozenset("N6 O6 N2 N7 N3 C5 C4 N1 C8 C6 O4 O2 N9 N4 C2".split(" "))
PHOS_GROUP_ATOMS = frozenset(["P", "OP1", "OP2", "OP3"])

CONNECTED_DISTANCE_CUTOFF = 2.5 #how long can the O3'-P bond be before two nucleotides are considered not connected

class PseudoMolecule:
    """A class for storing atomic coordinate data and update a Coot molecule and drawing objects as appropriate"""
    
    def __init__(self, batons = True, initialResList = None, moleculeNumber = None, createFromMolecule = None, chain = None):
        """Initialize a new pseudoMolecule object.
        
        OPTIONAL ARGUMENTS:
            batons             - Whether or not to create a generic display object for drawing batons between phosphates and bases.  Defaults to True
            initialResList     - A list of nucleotide coordinates to store in the pseudoMolecule object.  Defaults to None.
                                 Formatted as a list of dictionaries, where each dictionary contains atomName: [x, y, z]
            moleculeNumber     - If given, this Coot molecule will be updated/overwritten.  Otherwise a new Coot molecule will be created.
            createFromMolecule - create the PseudoMolecule object using the specified existing Coot molecule
                                 This differs from moleculeNumber in that moleculeNumber will overwrite an existing molecule.  CreateFromMolecule will
                                 use the coordinates of the existing molecule.
            chain              - the chain contains nucleotides to be accessed/modified via this object.
                                 Only relevant when createFromMolecule is defined.
                                 If not given, the first chain in the molecule will be used
        RETURNS:
            an initialized pseudoMolecule object
        EXCEPTIONS:
            raises PseudoMoleculeError if the specified chain doesn't exist
            or if molecule createFromMolecule contains anisotropic temperature records and we're using Coot older than rev 3631
        """
        
        if batons:
            #create an object for the batons and set it to display
            self.__batonObject = new_generic_object_number("Accepted Batons")
            set_display_generic_object(self.__batonObject, 1)
        
        self.__molecule = [ [ ["A", [] ] ] ] #a data structure containing all atoms of the molecule
                                             #in the format used by Coot
        self.__moleculeNumber = moleculeNumber #the number of this molecule
                                               #this number is used to interact with Coot's representation of the molecule
        self.__extraBondObject        = None
        self.__resNumDict             = {}
        self.__extraBondStartingIndex = None
        self.__extraBondEndingIndex   = None
        self.__savedCoordinates       = None
        self.__savedBfacs             = None
        self.chain                    = "A"
        self.__chainIndex             = 0
        self.resInsertionPoint      = None
        self.origResInsertionPoint  = None
        self.__savedMoleculeState     = None
        
        if createFromMolecule is not None:
            self.__moleculeNumber = createFromMolecule
            self.__molecule = python_representation(createFromMolecule)
            self.__checkAniso() #if the version of Coot isn't new enough, then we can't handle structures with anisotropic temperature records
            self.__fixTerAtoms() #workaround a Coot bug with TER records in the PDB file
            
            #figure out which chain we want
            if chain is not None:
                self.chain = chain
                try:
                    self.__chainIndex = [chainList[0] for chainList in self.__molecule[0]].index(chain)
                    #print "__chainIndex =", self.__chainIndex
                    #print "chains:", [chainList[0] for chainList in self.__molecule[0]]
                except ValueError:
                    raise PseudoMoleculeError("Chain " + str(chain) + " not found")
            
            #initialize the residue index lookup
            self.__initializeResNumDict()
            
            #we don't need to clear_and_update_molecule since we haven't changed anything about it
        
        #if there's an initial residue list, then add those residues and create the molecule in Coot
        elif initialResList:
            self.__molecule[0][0][1] = initialResList
            self.__initializeResNumDict()
            
            #create or update Coot's representation of the molecule
            if self.__moleculeNumber is None:
                self.__moleculeNumber = add_molecule(self.__molecule, MOLECULE_NAME % moleculeNameCount.next())
            else:
                clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        

    
    
    def addPhos (self, coords):
        """Add a new phosphate atom to the pseudoMolecule object.
        
        ARGUMENTS:
            coords - the coordinates of the phosphate to add
        RETURNS:
            None
        EFFECTS:
            The phosphate will be added to the Coot molecule.
            If this is the first phosphate to be added, a new Coot molecule object will be
                created and displayed.
        """
        
        resList = self.__molecule[0][self.__chainIndex][1]
        
        if self.resInsertionPoint is None:
            resInsertionPoint = len(resList)
        else:
            resInsertionPoint = self.resInsertionPoint + 1
            self.resInsertionPoint += 1
        
        #determine the residue number for the next residue
        nextResNum = None
        if len(resList) == 0:
            #if this is the first residue, then use residue number 1
            nextResNum = 1
        else:
            #if this isn't the first residue, then increment the previous residue's number
            nextResNum = resList[resInsertionPoint-1][0] + 1
        
        #update the residue number lookup dictionary
        if self.resInsertionPoint is None:
            self.__resNumDict[str(nextResNum)] = len(resList)
        else:
            self.__initializeResNumDict()
        
        #create a new residue containing the specified atom and add it to self.__molecule
        newRes = [nextResNum, "", "G", [[[" P  ",""], [1.0, DEFAULT_B_FACTOR, " P"], coords[0:3]]]]
        resList.insert(resInsertionPoint, newRes)
        self.__molecule[0][self.__chainIndex][1] = resList
        
        #create or update Coot's representation of the molecule
        if self.__moleculeNumber is None:
            self.__moleculeNumber = add_molecule(self.__molecule, MOLECULE_NAME % moleculeNameCount.next())
        else:
            clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        
        #if this isn't the first atom, then draw a baton to it
        #if len(resList) > 1:
        #    prevPhosCoords = self.getPhosCoords(-1)
        #    to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, prevPhosCoords[0], prevPhosCoords[1], prevPhosCoords[2], coords[0], coords[1], coords[2])
        #    graphics_draw()
    
    
    def addPhos5p(self, coords):
        """Add a new 5' phosphate atom to the pseudoMolecule object.
        
        ARGUMENTS:
            coords - the coordinates of the phosphate to add
        RETURNS:
            None
        NOTE:
            This function should only be used during an extendChain, not when tracing a new molecule
        """
        
        self.__molecule[0][self.__chainIndex][1][self.resInsertionPoint][3].insert(0, [[" P  ",""], [1.0, DEFAULT_B_FACTOR, " P"], coords[0:3]])
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        
    
    def addBaseAndPhos(self, baseType, baseCoords, phosCoords):
        """Add a base and the 3' phosphate (i.e. the next phosphate) to the pseudoMolecule object
        
        ARGUMENTS:
            baseType    - the type of base (i.e. "A", "G", "C", or "U")
            baseCoords  - coordinates of the base in the format atomName: [x, y, z]
                          Note that the C1' coordinates MUST be included
            phosCoords  - the coordinates of the 3' phosphate
        RETURNS:
            None
        EFFECTS:
            The base and phosphate will be added to the Coot molecule and batons will be drawn connecting
            the phosphates and the C1' atom.
        """
        
        resList = self.__molecule[0][self.__chainIndex][1]
        
        if self.resInsertionPoint is None:
            resInsertionPoint = len(resList)
        else:
            self.resInsertionPoint += 1
            resInsertionPoint = self.resInsertionPoint
        
        #pprint(resList)
        #print "Adding %s base" % baseType
        
        #modify the base type of the last residue (since we didn't know what type of base it was when we were adding the phosphate)
        resList[resInsertionPoint-1][2] = baseType
        
        #start the list of base atoms with the C1' atom so we can easily find it later (for getSugarCoords)
        baseAtomList = [[[" C1'",""], [1.0, DEFAULT_B_FACTOR, " C"], baseCoords["C1'"][0:3]]]
        
        for curAtom in baseCoords:
            #we've already added the C1' to atomList
            if curAtom == "C1'": continue
            
            baseAtomList.append([[curAtom.center(4),""], [1.0, DEFAULT_B_FACTOR, " %s" % curAtom[0:1]], baseCoords[curAtom][0:3]])
        
        resList[resInsertionPoint-1][3].extend(baseAtomList)
        
        
        #determine the residue number for the next residue
        nextResNum = None
        if len(resList) == 0:
            #if this is the first residue, then use residue number 1
            nextResNum = 1
        else:
            #if this isn't the first residue, then increment the previous residue's number
            nextResNum = resList[resInsertionPoint-1][0] + 1
            if nextResNum == 0: nextResNum = 1 #go from -1 to 1
        
        #create a new residue containing the specified atom and add it to self.__molecule
        newRes = [nextResNum, "", baseType, [[[" P  ",""], [1.0, DEFAULT_B_FACTOR, " P"], phosCoords[0:3]]]]
        resList.insert(resInsertionPoint, newRes)
        self.__molecule[0][self.__chainIndex][1] = resList
        
        #update the residue number lookup dictionary
        if self.resInsertionPoint is None:
            #if we're building a new chain, then we just have to add the new residue number
            self.__resNumDict[str(nextResNum)] = resInsertionPoint
        else:
            #if we're extending a chain, then we may be inserting a residue into the middle of the molecule
            #so we just regenerate the entire residue number dictionary
            self.__initializeResNumDict()
            #we could try to be clever and update only the dictionary elements that actually change
            #but I suspect that the cleverness would take longer than just regenerating the entire dictionary
            #(not to mention be far more bug-prone)
        
        #create or update Coot's representation of the molecule
        if self.__moleculeNumber is None:
            self.__moleculeNumber = add_molecule(self.__molecule, MOLECULE_NAME % moleculeNameCount.next())
        else:
            clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        
        #if this isn't the first atom, then draw a baton to it
        #if len(resList) > 1:
        prevPhosCoords = self.getPhosCoordsFromIndex(resInsertionPoint - 1)
        c1coords = baseCoords["C1'"]
        #pprint(prevPhosCoords)
        to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, prevPhosCoords[0], prevPhosCoords[1], prevPhosCoords[2], c1coords[0], c1coords[1], c1coords[2])
        to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, c1coords[0], c1coords[1], c1coords[2], phosCoords[0], phosCoords[1], phosCoords[2])
        graphics_draw()
        
    
    def addBaseAndPhos5p(self, baseType, baseCoords, phosCoords):
        """Add a 5' base and phosphate to the pseudoMolecule object
        
        ARGUMENTS:
            baseType    - the type of base (i.e. "A", "G", "C", or "U")
            baseCoords  - coordinates of the base in the format atomName: [x, y, z]
                          Note that the C1' coordinates MUST be included
            phosCoords  - the coordinates of the phosphate
        RETURNS:
            None
        EFFECTS:
            The base and phosphate will be added to the Coot molecule and batons will be drawn connecting
            the phosphates and the C1' atom.
        """
        
        resList = self.__molecule[0][self.__chainIndex][1]
        
        if self.resInsertionPoint is None:
            resInsertionPoint = 0
            newResNum = 1
        else:
            newResNum = resList[self.resInsertionPoint][0] - 1
            if newResNum == 0:
                newResNum = -1
            self.origResInsertionPoint += 1
            resInsertionPoint = self.resInsertionPoint
        
        #create the atom list for the new residue, starting with the phosphate and the C1'
        atomList = [[[" P  ",""], [1.0, DEFAULT_B_FACTOR, " P"], phosCoords[0:3]],
                    [[" C1'",""], [1.0, DEFAULT_B_FACTOR, " C"], baseCoords["C1'"][0:3]]]
        
        for curAtom in baseCoords:
            #we've already added the C1' to atomList
            if curAtom == "C1'": continue
            
            atomList.append([[curAtom.center(4),""], [1.0, DEFAULT_B_FACTOR, " %s" % curAtom[0:1]], baseCoords[curAtom][0:3]])
        
        #create a new residue (numbered 1) containing the atomList
        newRes = [newResNum, "", baseType, atomList]
        
        #if we're building a new chain (as opposed to extending an existing chain), increment the residue number of all other residues in the chain
        if self.resInsertionPoint is None:
            for curRes in resList:
                curRes[0] += 1
        
        #prepend the new residue at the beginning of the residue list (or the appropriate spot in the residue list if we're doing an extend chain)
        resList.insert(resInsertionPoint, newRes)
        
        self.__molecule[0][self.__chainIndex][1] = resList
        
        #update the resNumDict
        if self.resInsertionPoint is None:
            new5pResNum = resList[-1][0] #this is the new 5' residue number, since we're incrementing all the existing residues
            self.__resNumDict[str(new5pResNum)] = len(resList) - 1
        else:
            #self.__resNumDict[str(newResNum)] = resInsertionPoint
            
            #if we're extending a chain, then we may be inserting a residue into the middle of the molecule
            #so we just regenerate the entire residue number dictionary
            self.__initializeResNumDict()
            #we could try to be clever and update only the dictionary elements that actually change
            #but I suspect that the cleverness would take longer than just regenerating the entire dictionary
            #(not to mention be far more bug-prone)
        
        #create or update Coot's representation of the molecule
        if self.__moleculeNumber is None:
            self.__moleculeNumber = add_molecule(self.__molecule, MOLECULE_NAME % moleculeNameCount.next())
        else:
            clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        
        #draw batons to it
        prevPhosCoords = self.getPhosCoordsFromIndex(resInsertionPoint+1)
        c1coords = baseCoords["C1'"]
        #pprint(prevPhosCoords)
        to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, prevPhosCoords[0], prevPhosCoords[1], prevPhosCoords[2], c1coords[0], c1coords[1], c1coords[2])
        to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, c1coords[0], c1coords[1], c1coords[2], phosCoords[0], phosCoords[1], phosCoords[2])
        graphics_draw()
    
    
    def getPhosCoords(self, resNum):
        """Get the coordinates of a specified phosphate.
        
        ARGUMENTS:
            resNum - the residue number
        RETURNS:
            the coordinates of the phoshate from residue resNum 
        NOTE:
            Negative integers can be used to refer to phosphates at the end of the chain (i.e. resNum = -1 will return
            the coordinates of the last phosphate).  If you want to refer to a nucleotide with an actual negative residue
            number (very uncommon, but they do exist), then pass the number as a string (i.e. "-1")
        """
        
        if isinstance(resNum, int) and resNum < 0:
            resIndex = resNum
        else:
            resIndex = self.resIndex(resNum)
        
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            if curAtom[0][0].strip() == "P":
                return curAtom[2]
        
        #if we haven't found a phosphate atom
        return None
    
    
    def getPhosCoordsFromIndex(self, resNum):
        """Get the coordinates of a specified phosphate.
        
        ARGUMENTS:
            resIndex - the residue index
        RETURNS:
            the coordinates of the phoshate from residue resIndex
        """
        
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            if curAtom[0][0].strip() == "P":
                return curAtom[2]
        
        #if we haven't found a phosphate atom
        return None
    
    
    def getSugarCoordsFromIndex(self, resIndex):
        """Get the coordinates of a specified C1' atom (NOT the entire sugar).
        
        ARGUMENTS:
            resIndex - the residue index
        RETURNS:
            the coordinates of the C1' from the specified residue
        NOTE:
            Negative numbers can be used to refer to C1' atoms at the end of the chain (i.e. resNum = -1 will return the coordinates of the last C1')
        """
        
        #if the specified residue was created/processed by RCrane, then the C1' atom will be the 2nd atom in the resideu
        #but this isn't necessarily true (and probably isn't) for other RNA nts
        #so we search through all the atoms in the nucleotide
        
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            if curAtom[0][0].strip() == "C1'":
                return curAtom[2]
    
    #def getCootRes(self, resNum):
    #    return deepcopy(self.__molecule[0][self.__chainIndex[1][resNum])
    
    def getAtomCoordsFromIndex(self, atomName, resIndex):
        """Get the coordinates of a specified atom
        
        ARGUMENTS:
            atomName - the name of the atom to get coordinates for
            resIndex - the residue index
        RETURNS:
            the coordinates of the specified atom
        """
        
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            if curAtom[0][0].strip() == atomName:
                return curAtom[2]
        
        #if we haven't found the desired atom
        return None
    
    def closeBatonObject(self):
        """close the Coot generic object for the batons connecting the phosphate and C1' atoms.
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            The batons generic object will be closed
        """
        
        close_generic_object(self.__batonObject)
        
    def getNumNts(self):
        """Get the number of nucleotides in the molecule.
        
        ARGUMENTS:
            None
        RETURNS:
            the number of nucleotides in the molecule
        """
        
        return len(self.__molecule[0][self.__chainIndex][1])
    
    def removeLastBaseAndPhos(self):
        """Remove and return the last nucleotide of the chain
        
        ARGUMENTS:
            None
        RETURNS:
            nucType - the type of base (i.e. "A", "G", "C", or "U")
            coords  - the coordinates of the removed nucleotide in the format atomName: [x, y, z]
        EFFECTS:
            The last nucleotide of the chain will be removed from the Coot molecule, as will the batons
            connecting the phosphates and the C1' atom of that nucleotide
        NOTE:
            This function assumes that the phosphate is the first atom of the penultimate nucleotide.  This
            is always true for molecules built using the addPhos() and addBaseAndPhos() functions, but may
            not be true for arbitrary molecules.  As a result, this function should *not* be called for
            pseudoMolecule objects created via the createFromMolecule argument.
        """
        
        if (self.getNumNts() == 1):
            #Coot doesn't like displaying empty molecules, so never shrink beyond the initial phosphate
            return (None, None)
        
        if self.resInsertionPoint is None:
            lastPhosNtIndex = len(self.__molecule[0][self.__chainIndex][1])
        else:
            lastPhosNtIndex = self.resInsertionPoint
        lastFullNtIndex =  lastPhosNtIndex - 1
        
        #retrieve the 5' phosphate and the last full nucleotide from the __molecule object
        lastPhosNt = self.__molecule[0][self.__chainIndex][1][lastPhosNtIndex]
        lastFullNt = self.__molecule[0][self.__chainIndex][1][lastFullNtIndex]
        
        #get the nucleotide type
        nucType = lastFullNt[2]
        
        #convert the nucleotide into a dictionary
        coords = dict()
        for curAtomData in lastFullNt[3][1:]:
            atomName = curAtomData[0][0].strip()
            curCoords = curAtomData[2]
            coords[atomName] = curCoords
        coords["P"] = lastPhosNt[3][0][2]
        #Note that coords has a base and it's 5' phosphate, so it's not a proper nucleotide
        
        #remove the last phosphate nt
        del self.__molecule[0][self.__chainIndex][1][lastPhosNtIndex]
        
        #truncate the atom list of last nucleotide to just the phosphate
        self.__molecule[0][self.__chainIndex][1][lastFullNtIndex][3] = self.__molecule[0][self.__chainIndex][1][lastFullNtIndex][3][0:1]
        
        #update the resNumDict
        if self.resInsertionPoint is None:
            del self.__resNumDict[str(self.__molecule[0][self.__chainIndex][1][-1][0])]
        else:
            self.__initializeResNumDict()
            
            #also, decrement resInsertionPoint
            self.resInsertionPoint -= 1
        
        #there's no way to selectively remove lines from a drawing object, so we have to destroy the Accepted Batons object
        close_generic_object(self.__batonObject)
        self.__batonObject = new_generic_object_number("Accepted Batons")
        set_display_generic_object(self.__batonObject, 1)
        
        #redraw the entire molecule (or the new portion of the molecule if we're doing an extendChain)
        if self.resInsertionPoint is None:
            nucsToDraw = xrange(len(self.__molecule[0][self.__chainIndex][1])-1)
        else:
            nucsToDraw = xrange(self.origResInsertionPoint, self.resInsertionPoint)
        #print list(nucsToDraw)
        
        prevPhosCoords = self.__molecule[0][self.__chainIndex][1][0][3][0][2]
        for i in nucsToDraw:
            
            phosCoords     = self.getPhosCoordsFromIndex(i)
            c1Coords       = self.getSugarCoordsFromIndex(i)
            nextPhosCoords = self.getPhosCoordsFromIndex(i+1)
            
            to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, phosCoords[0], phosCoords[1], phosCoords[2], c1Coords[0], c1Coords[1], c1Coords[2])
            to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, c1Coords[0], c1Coords[1], c1Coords[2], nextPhosCoords[0], nextPhosCoords[1], nextPhosCoords[2])
        
        #redraw the graphics
        graphics_draw()
        
        #update Coot's representation of the molecule
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        
        return (nucType, coords)
    
    
    def removeFirstBaseAndPhos(self):
        """Remove and return the first nucleotide of the chain
        
        ARGUMENTS:
            None
        RETURNS:
            nucType - the type of base (i.e. "A", "G", "C", or "U")
            coords  - the coordinates of the removed nucleotide in the format atomName: [x, y, z]
        EFFECTS:
            The first nucleotide of the chain will be removed from the Coot molecule, as will the batons
            connecting the phosphates and the C1' atom of that nucleotide
        NOTE:
            removeLastBaseAndPhos is used for Previous Nt when tracing 5'->3'.  This function provides
            the equivalent functionality when tracing 3'->5'.
            This function should *not* be called for pseudoMolecule objects created via the createFromMolecule argument.
        """
        
        if (self.getNumNts() == 1):
            #Coot doesn't like displaying empty molecules, so never shrink beyond the initial phosphate
            return (None, None)
        
        #if we're doing an extendChain, update resInsertionPoint
        if self.origResInsertionPoint is not None:
            self.origResInsertionPoint -= 1
        
        #remove the 3' nucleotide from __molecule
        firstNt = self.__molecule[0][self.__chainIndex][1].pop(0)
        
        #get the nucleotide type
        nucType = firstNt[2]
        
        #convert the nucleotide into a dictionary
        coords = dict()
        for curAtomData in firstNt[3]:
            atomName = curAtomData[0][0].strip()
            curCoords = curAtomData[2]
            coords[atomName] = curCoords
        
        if self.resInsertionPoint is None:
            #update resNumDict
            del self.__resNumDict[str(self.__molecule[0][self.__chainIndex][1][-1][0])]
            
            #decrement the residue number of all other residues in the chain
            for curRes in self.__molecule[0][self.__chainIndex][1]:
                curRes[0] -= 1
        else:
            #if we're doing an extend chain, then regenerate the entire resNumDict
            self.__initializeResNumDict()
            
            #and don't renumber any residues
        
        #there's no way to selectively remove lines from a drawing object, so we have to destroy the Accepted Batons object
        close_generic_object(self.__batonObject)
        self.__batonObject = new_generic_object_number("Accepted Batons")
        set_display_generic_object(self.__batonObject, 1)
        
        #redraw the entire molecule (or the new portion of the molecule if we're doing an extendChain)
        if self.resInsertionPoint is None:
            nucsToDraw = xrange(len(self.__molecule[0][self.__chainIndex][1])-1)
        else:
            nucsToDraw = xrange(self.resInsertionPoint, self.origResInsertionPoint)
        prevPhosCoords = self.__molecule[0][self.__chainIndex][1][0][3][0][2]
        for i in nucsToDraw:
            
            phosCoords     = self.getPhosCoordsFromIndex(i)
            c1Coords       = self.getSugarCoordsFromIndex(i)
            nextPhosCoords = self.getPhosCoordsFromIndex(i+1)
            
            to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, phosCoords[0], phosCoords[1], phosCoords[2], c1Coords[0], c1Coords[1], c1Coords[2])
            to_generic_object_add_line(self.__batonObject, BATONCOLOR, 6, c1Coords[0], c1Coords[1], c1Coords[2], nextPhosCoords[0], nextPhosCoords[1], nextPhosCoords[2])
        
        #redraw the graphics
        graphics_draw()
        
        #update Coot's representation of the molecule
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        
        return (nucType, coords)
    
    
    def createChainObject(self):
        """Return a chain object containing all coordinates
        
        ARGUMENTS:
            None
        RETURNS:
            a list of coordinates for each nucleotide, where the coordinates are formatted as atomName: [x, y, z]
        """
        
        chainObj = Chain()
        
        for curRes in self.__molecule[0][self.__chainIndex][1]:
            resName = curRes[2]
            resNum  = str(curRes[0]) + str(curRes[1]).strip()
            
            #create a dictionary at atomic coordinates
            atomDict = {}
            for curAtom in curRes[3]:
                atomName = curAtom[0][0].strip()
                atomCoords = curAtom[2]
                atomDict[atomName] = atomCoords
            
            nucObj = Nucleotide(resName, atomDict, resNum)
            chainObj.addNuc(nucObj)
        
        return chainObj
    
    
    def createPartialChainObject(self, startRes, endRes, addFlankingAtoms = False):
        """Create a chain object containing coordinates for some nucleotides
        
        ARGUMENTS:
            startRes - the starting residue
            endRes   - the ending residue
        OPTIONAL ARGUMENTS:
            addFlankingAtoms - if True, add the 5' and 3' atoms that will move during minimization
                               Note that these atoms are only added if the 5' and 3' nucleotides are connected
                               to startRes and endRes, respectively.
                               In particular, the 3' phosphate, non-bridging oxygens, O5', and C5' will be added
                               and the 5' O3', and C3'.
                               (The 3' C5' and 5' C3' don't actually move during the minimization, but the C5'-O5' and
                               C3'-O3' bond can move due to movement in the O5' and O3' atoms, so we include the carbons
                               in the chain object so that the bonds can be drawn.)
                               Defaults to False
        RETURNS:
            a list of coordinates for each nucleotide, where the coordinates are formatted as atomName: [x, y, z]
            if addFlankingAtoms is True, then the following are also returned:
                resNum5p - the number of the 5' nucleotide, or None if the 5' nucleotide is not connected (or non-existant)
                resNum3p - the number of the 3' nucleotide, or None if the 3' nucleotide is not connected (or non-existant)
        """
        
        #convert the residue numbers to indices
        startIndex = self.resIndex(startRes)
        endIndex = self.resIndex(endRes)
        
        return self.createPartialChainObjectFromIndex(startIndex, endIndex, addFlankingAtoms)
    
    
    def createPartialChainObjectFromIndex(self, startIndex, endIndex, addFlankingAtoms = False):
        """Create a chain object containing coordinates for some nucleotides
        
        ARGUMENTS:
            startIndex - the starting residue index
            endIndex   - the ending residue index
        OPTIONAL ARGUMENTS:
            addFlankingAtoms - if True, add the 5' and 3' atoms that will move during minimization
                               Note that these atoms are only added if the 5' and 3' nucleotides are connected
                               to startRes and endRes, respectively.
                               In particular, the 3' phosphate, non-bridging oxygens, O5', and C5' will be added
                               and the 5' O3', and C3'.
                               (The 3' C5' and 5' C3' don't actually move during the minimization, but the C5'-O5' and
                               C3'-O3' bond can move due to movement in the O5' and O3' atoms, so we include the carbons
                               in the chain object so that the bonds can be drawn.)
                               Defaults to False
        RETURNS:
            a list of coordinates for each nucleotide, where the coordinates are formatted as atomName: [x, y, z]
            if addFlankingAtoms is True, then the following are also returned:
                resNum5p - the number of the 5' nucleotide, or None if the 5' nucleotide is not connected (or non-existant)
                resNum3p - the number of the 3' nucleotide, or None if the 3' nucleotide is not connected (or non-existant)
        """
        
        #make sure that the start is before the end, and swap if its not
        if startIndex > endIndex:
            (endIndex, startIndex) = (startIndex, endIndex)
        
        chainObj = Chain() #initialize the chain object
        resNum5p = None
        resNum3p = None
        
        #add the 5' flanking atoms if we're not at the very start of the chain
        if addFlankingAtoms and self.connectedToPrevFromIndex(startIndex):
            curRes = self.__molecule[0][self.__chainIndex][1][startIndex-1]
            resName = curRes[2]
            resNum5p  = str(curRes[0]) + str(curRes[1]).strip()
            
            atomDict = {}
            for curAtom in curRes[3]:
                atomName = curAtom[0][0].strip()
                #print "examining atomName =", atomName
                if atomName in frozenset(["O3'", "C3'"]):
                    atomDict[atomName] = curAtom[2]
            
            nucObj = Nucleotide(resName, atomDict, resNum5p)
            chainObj.addNuc(nucObj)
        
        #add the specified nucleotides
        for curRes in self.__molecule[0][self.__chainIndex][1][startIndex:endIndex+1]:
            resName = curRes[2]
            resNum  = str(curRes[0]) + str(curRes[1]).strip()
            
            #create a dictionary at atomic coordinates
            atomDict = {}
            for curAtom in curRes[3]:
                atomName = curAtom[0][0].strip()
                atomCoords = curAtom[2]
                atomDict[atomName] = atomCoords
            
            nucObj = Nucleotide(resName, atomDict, resNum)
            chainObj.addNuc(nucObj)
        
        #add the 3' flanking atoms if we're not at the very end of the chain
        if addFlankingAtoms and self.connectedToNextFromIndex(endIndex):
            #add on the phosphate from the next nucleotide
            curRes = self.__molecule[0][self.__chainIndex][1][endIndex+1]
            resName = curRes[2]
            resNum3p  = str(curRes[0]) + str(curRes[1]).strip()
            
            atomDict = {}
            for curAtom in curRes[3]:
                atomName = curAtom[0][0].strip()
                #print "examining atomName =", atomName
                if atomName in frozenset(["P", "O5'", "C5'", "OP1", "OP2"]):
                    atomDict[atomName] = curAtom[2]
            
            nucObj = Nucleotide(resName, atomDict, resNum3p)
            chainObj.addNuc(nucObj)
            
        if addFlankingAtoms:
            return (chainObj, resNum5p, resNum3p)
        else:
            return chainObj
    
    def addSugar(self, resNum, sugarCoords):
        """Add sugar atoms to an already existing residue
        
        ARGUMENTS:
            resNum      - the residue to add the coordinates to (using Coot numbering)
            sugarCoords - the sugar coordinates in the form atomName: [x, y, z]
        RETURNS:
            None
        EFFECTS:
            The new sugar atoms will be drawn, but the batons connecting the phosphate and C1' atom will NOT be removed
        """
        
        atomList = []
        
        for (atomName, coords) in sugarCoords.iteritems():
            #convert PDB3 format atoms names to PDB2
            atomList.append([["%4s" % atomName,""], [1.0, DEFAULT_B_FACTOR, " %s" % atomName[0:1]], coords])
        
        resIndex = self.resIndex(resNum) #convert resNum (Coot numbering) to an index (internal numbering)
        self.__molecule[0][self.__chainIndex][1][resIndex][3].extend(atomList)
        
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
    
    
    def replaceSugar(self, resNum, sugarCoords):
        """Replace the sugar atoms in an already existing residue
        
        ARGUMENTS:
            resNum      - the residue to change the coordinates to (using Coot numbering)
            sugarCoords - the sugar coordinates in the form atomName: [x, y, z]
        RETURNS:
            None
        EFFECTS:
            The new sugar atoms will be drawn
        """
        
        resIndex = self.resIndex(resNum) #convert resNum (Coot numbering) to an index (internal numbering)
        
        for atom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            atomName = atom[0][0].strip()
            if sugarCoords.has_key(atomName):
                atom[2] = sugarCoords[atomName]
        
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
    
    
    def addPhosOxy(self, resNum, phosOxyCoords):
        """Add non-bridging phosphoryl oxygens to an already existing residue
        
        ARGUMENTS:
            resNum        - the residue to add the coordinates to using Coot numbering
            phosOxyCoords - the non-bridging oxygen coordinates in the form atomName: [x, y, z]
        RETURNS:
            None
        EFFECTS:
            The new oxygen atoms will be added to the Coot molecule
        """
        
        resIndex = self.resIndex(resNum)
        self.addPhosOxyFromIndex(resIndex, phosOxyCoords)
    
    
    def addPhosOxyFromIndex(self, resIndex, phosOxyCoords):
        """Add non-bridging phosphoryl oxygens to an already existing residue
        
        ARGUMENTS:
            resIndex      - the residue to add the coordinates to
            phosOxyCoords - the non-bridging oxygen coordinates in the form atomName: [x, y, z]
        RETURNS:
            None
        EFFECTS:
            The new oxygen atoms will be added to the Coot molecule
        """
        
        atomList = []
        for (atomName, coords) in phosOxyCoords.iteritems():
            #convert PDB3 format atoms names to PDB2
            if atomName == "OP1" or atomName == "OP2":
                atomName = " " + atomName
            else:
                continue
            
            atomList.append([[atomName,""], [1.0, DEFAULT_B_FACTOR, " O", ""], coords])
        
        self.__molecule[0][self.__chainIndex][1][resIndex][3].extend(atomList)
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
    
    
    def getPhosOxy(self, resNum):
        """Get the non-bridging phosphoryl oxygen coordinates for the specified residue.
        
        ARGUMENTS:
            resNum - the residue number
        RETURNS:
            the non-bridging oxygen coordinates of the specified residue as a dictionary in the form atomName: [x, y, z]
        NOTE:
            Negative numbers can be used to refer to residues at the end of the chain
        """
        
        phosOxyCoords = {}
        for curAtom in self.__molecule[0][self.__chainIndex][1][resNum][3][-2:]:
            atomName = curAtom[0][0]
            
            if atomName == " OP1" or atomName == " OP2":
                atomName = atomName.strip()
            else:
                continue
            
            phosOxyCoords[atomName] = curAtom[2]
        
        return phosOxyCoords
    
    
    def molNum(self):
        """Get the Coot molecule number used to display atoms from this object.
        
        ARGUMENTS:
            None
        RETURNS:
            the molecule number
        """
        
        return self.__moleculeNumber;
    
    
    def updateRes(self, startingResNum, endingResNum):
        """Update the Python representation of the specified residues.
        ARGUMENTS:
            startingResNum - the first residue to update, using Coot numbering
            endingResNum   - the last residue to update, using Coot numbering
        RETURNS:
            newCoordDicts - a list of dictionaries containing the updated atom coordinates
        EFFECTS:
            updates the Python representation of the specified residues
        NOTE:
            This function must be run whenever atom coordinates are changed in Coot (such as after a minimization run)
        """
        
        #note that minimization in Coot should not change the order of the atom list
        # if it does, this function needs to be more complicated
        
        startingResIndex = self.resIndex(startingResNum)
        endingResIndex   = self.resIndex(endingResNum)
        
        #print "in pseudoMol.updateRes"
        #print "  startingResIndex =", startingResIndex
        #print "  endingResIndex =", endingResIndex
        #print "  startingResNum =", startingResNum
        #print "  endingResNum =", endingResNum
        
        newCoordDicts = []
        for curResIndex in xrange(startingResIndex, endingResIndex+1):
            (curResNum, curInsCode) = self.resNum(curResIndex)
            
            #print "Updating residue %i" % curRes
            # 201805115-PE
            # the atoms from residue_info contain an atom index now - let's
            # strip it off to restore old functionality
            newCoords_4 = residue_info(self.__moleculeNumber, self.chain, curResNum, curInsCode)
            newCoords = [atom_bits[:3] for atom_bits in newCoords_4]
            self.__molecule[0][self.__chainIndex][1][curResIndex][3] = newCoords
            
            
            #convert the new atoms into a dictionary
            curResDict = {}
            for curAtom in newCoords:
                atomName = curAtom[0][0].strip()
                curResDict[atomName] = curAtom[2]
            newCoordDicts.append(curResDict)
            
        return newCoordDicts
    
    
    def deleteMolecule(self):
        """Delete the Coot molecule associated with this object.
        
        ARGUMENTS:
            None
        RETURNS:
            None
        """
        
        close_molecule(self.__moleculeNumber)
    

    def getCootNucs(self, startingRes, endingRes):
        """Retrieve the Coot structure for the specified nucleotides.
        
        ARGUMENTS:
            startingRes - the first residue to retrieve the structure for
            endingRes   - the last residue to retrieve the structure for
        RETURNS:
            a copy of the Coot representation of residues startingRes-endingRes
        NOTE
            This function is used to cache structure so that it can easily be restored later.
        """
        
        startingIndex = self.resIndex(startingRes) #convert startingRes (Coot numbering) to an index (internal numbering)
        endingIndex   = self.resIndex(endingRes)   #convert endingRes (Coot numbering) to an index (internal numbering)
        
        return deepcopy(self.__molecule[0][self.__chainIndex][1][startingIndex:endingIndex+1])
    
    def setCootNucs(self, startingRes, endingRes, cootNucs, updateMol = True):
        """Set the structure for the given nucleotides given a Coot-formatted list (typically something that had previously been returned by getCootNucs)
        
        ARGUMENTS:
            startingRes - the first residue to set the structure for
            endingRes   - the last residue to set the structure for
            cootNucs    - a list of the nucleotides to restore, formatted in the Coot internal representation format (typically returned from the getCootNucs function)
        OPTIONAL ARGUMENTS:
            updateMol   - whether to update the Coot molecule via clear_and_update_molecule
                          defaults to True
                          ONLY use this if you know the Coot molecule will be updated at a later time
        RETURNS:
            None
        EFFECTS:
            sets the coordinates of the specified nucleotides to cootNucs
        """
        
        startingIndex = self.resIndex(startingRes) #convert startingRes (Coot numbering) to an index (internal numbering)
        endingIndex   = self.resIndex(endingRes)   #convert endingRes (Coot numbering) to an index (internal numbering)
        
        self.__molecule[0][self.__chainIndex][1][startingIndex:endingIndex+1] = deepcopy(cootNucs)
        if updateMol:
            clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
    
    
    def resetNucs(self, startingRes, endingRes, origCoords = None, intermediateAtomLocs = None):
        """Reset the specified nucleotides to the condition they were in immediately before their minimization
        
        ARGUMENTS:
            startingRes          - the first residue to reset
            endingRes            - the last residue to reset
        OPTIONAL ARGUMENTS
            origCoords           - a list of the orginal phosphate and base coordinates for all nucleotides
                                   if not provided, then the phosphate and base coordinates will not be reset
                                   but all other atoms will still be deleted
                                   must be given if intermediateAtomLocs is given
            intermediateAtomLocs - a list of P, O5', and C5' coodinates as they were immediately after the
                                   minimization of the previous nucleotide, but before the minimization of
                                   that nucleotide
                                   if not provided, then these atom coordinates will not be reset
        RETURNS:
            None
        EFFECTS:
            Resets the coordinates of the specified nucleotides
        """
        
        #figure out the starting indices for the self.__molecule residue list
        startingIndex = self.resIndex(startingRes) #convert startingRes (Coot numbering) to an index (internal numbering)
        endingIndex   = self.resIndex(endingRes)   #convert endingRes (Coot numbering) to an index (internal numbering)
        
        #origCoords and intermediateAtomLocs may not span the full molecule, so their indices may be offset
        if origCoords is not None:
            origCoordsOffset = startingIndex - origCoords.resIndex(startingRes)
        else:
            origCoordsOffset = 0
        
        resList = self.__molecule[0][self.__chainIndex][1]
        
        for curResIndex in xrange(startingIndex, endingIndex+1):
            
            origCoordsCurResIndex = curResIndex - origCoordsOffset
            
            #print "In resetNucs loop:"
            #print "\tcurResIndex           =", curResIndex
            #print "\torigCoordsCurResIndex =", origCoordsCurResIndex
            
            #erase the sugar atoms, as we're going to reset these to their default locations
            #and erase hydrogens, since Coot doesn't have restraints for them (and RCrane doesn't know how to build them correctly anyway)
            for curAtomIndex in xrange(len(resList[curResIndex][3])-1, -1, -1):
                curAtomName = resList[curResIndex][3][curAtomIndex][0][0].strip()
                curAtomType = resList[curResIndex][3][curAtomIndex][1][2].strip()
                
                #print "Examining residue", curResIndex, "atom", curAtomName, "type", curAtomType
                if curAtomType == "H" or curAtomName in ATOMS_TO_ERASE_FOR_REBUILDING_MOLECULE:
                    del resList[curResIndex][3][curAtomIndex]
            
            #if this is the first residue to be reset but NOT the first residue of the chain,
            #then replace the phosphate with the intermediate phosphate coords (and the previous O3' as well)
            if intermediateAtomLocs is not None and curResIndex == startingIndex and curResIndex != 0 and origCoordsCurResIndex != 0 and intermediateAtomLocs[origCoordsCurResIndex] is not None:
                
                #replace the phosphate of the current nucleotide
                for curAtomIndex in xrange(len(resList[curResIndex][3])):
                    if resList[curResIndex][3][curAtomIndex][0][0].strip() == "P":
                        resList[curResIndex][3][curAtomIndex][2] = intermediateAtomLocs[origCoordsCurResIndex][0]
                        break
                
                #replace the O3' of the previous residue
                for curAtomIndex in xrange(len(resList[curResIndex-1][3])):
                    if resList[curResIndex-1][3][curAtomIndex][0][0].strip() == "O3'":
                        resList[curResIndex-1][3][curAtomIndex][2] = intermediateAtomLocs[origCoordsCurResIndex][1]
            elif origCoords is not None:
                for curAtomIndex in xrange(len(resList[curResIndex][3])):
                    if resList[curResIndex][3][curAtomIndex][0][0].strip() == "P":
                        resList[curResIndex][3][curAtomIndex][2] = origCoords.nucleotides[origCoordsCurResIndex].atoms["P"]
                        break
                
            #reset base location, since it can move slightly due to the harmonic restraints
            #NOTE: if I ever turn off harmonic restraints on the base atoms and go back to fixing them,
            #then I can safely comment out this loop
            if origCoords is not None:
                for curAtomIndex in xrange(len(resList[curResIndex][3])):
                    curAtomName = resList[curResIndex][3][curAtomIndex][0][0].strip()
                    if curAtomName in BASE_ATOMS or curAtomName == "C1'":
                        resList[curResIndex][3][curAtomIndex][2] = origCoords.nucleotides[origCoordsCurResIndex].atoms[curAtomName]
                        #print "Resetting base coordinates for", curAtomName, "in residue", curResIndex, "to", origCoords.nucleotides[curResIndex].atoms[curAtomName]
        
        #if the endingRes is the last residue of a segment, then get rid of the 3' non-bridging oxygens
        curResIndex = endingIndex+1
        
        if self.connectedToNextFromIndex(endingIndex) and (len(resList[curResIndex][3]) <= 4) and not self.connectedToNextFromIndex(curResIndex):
            #if the next residue has more than four atoms, then it must be more than just a phosphate
            #(remember that terminal phosphates can have a third non-bridging oxygen)
            for curAtomIndex in xrange(len(resList[curResIndex][3])-1, -1, -1):
                curAtomName = resList[curResIndex][3][curAtomIndex][0][0].strip()
                #print "Examining residue", curIndex, "atom", curAtomName
                if curAtomName in frozenset(["OP1", "OP2", "OP3"]):
                    del resList[curResIndex][3][curAtomIndex]
            
        
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        #graphics_draw() #redrawing the screen isn't necessary, although it can be useful when debugging
        
    
    def drawExtraBonds(self):
        """Draw any bonds that are too long for Coot to render
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            Draws bonds into drawing object self.__extraBondObject
        NOTE:
            If the extra bond starting and ending residues have been set via a setExtraBondRange() call, then this
            function will only draw extra bonds for the specified residues.  Otherwise extra bonds will be drawn for
            the entire molecule
        """
        
        self.clearExtraBonds()
        
        resList = self.__molecule[0][self.__chainIndex][1][self.__extraBondStartingIndex:self.__extraBondEndingIndex]
        lineList = []
        
        prevO3 = None
        #examine all residues except the final one (which is just the terminal 3' P and OP1+OP2)
        for curRes in resList[:-1]:
            #create a dictionary with all of the atomic coordinates
            curResCoords = dict([(curAtom[0][0], curAtom[2]) for curAtom in curRes[3]])
            #pprint(curRes)
            #pprint(curResCoords)
            
            for (atom1, atom2, distCutoff) in BOND_LIST:
                if curResCoords.has_key(atom1) and curResCoords.has_key(atom2) and dist(curResCoords[atom1], curResCoords[atom2]) > distCutoff:
                    lineList.append([BATONCOLOR, 6, curResCoords[atom1][0], curResCoords[atom1][1], curResCoords[atom1][2], curResCoords[atom2][0], curResCoords[atom2][1], curResCoords[atom2][2]])
            
            #draw the O3'-1 - P bond if necessary
            if prevO3 is not None and (dist(prevO3, curResCoords[" P  "]) > BOND_DIST_CUTOFF):
                lineList.append([BATONCOLOR, 6, prevO3[0], prevO3[1], prevO3[2], curResCoords[" P  "][0], curResCoords[" P  "][1], curResCoords[" P  "][2]])
            prevO3 = curResCoords[" O3'"]
            
            #draw the glycosidic bond if necessary
            resName = curRes[2].strip()
            if resName == "A" or resName == "G":
                if dist(curResCoords[" C1'"], curResCoords[" N9 "]) > BOND_DIST_CUTOFF:
                    lineList.append([BATONCOLOR, 6, curResCoords[" C1'"][0], curResCoords[" C1'"][1], curResCoords[" C1'"][2], curResCoords[" N9 "][0], curResCoords[" N9 "][1], curResCoords[" N9 "][2]])
            else:
                if dist(curResCoords[" C1'"], curResCoords[" N1 "]) > BOND_DIST_CUTOFF:
                    lineList.append([BATONCOLOR, 6, curResCoords[" C1'"][0], curResCoords[" C1'"][1], curResCoords[" C1'"][2], curResCoords[" N1 "][0], curResCoords[" N1 "][1], curResCoords[" N1 "][2]])
        
        
        #handle OP1 and OP2 of the last nucleotide
        #I doubt that these bonds will ever get overly stretched, but it can't hurt to double check
        curRes = resList[-1]
        curResCoords = dict([(curAtom[0][0], curAtom[2]) for curAtom in curRes[3]])
        for (atom1, atom2, distCutoff) in BOND_LIST[0:2]: #the first two elements of BOND_LIST are the OP1 and OP2 bonds
            if dist(curResCoords[atom1], curResCoords[atom2]) > distCutoff:
                lineList.append([BATONCOLOR, 6, curResCoords[atom1][0], curResCoords[atom1][1], curResCoords[atom1][2], curResCoords[atom2][0], curResCoords[atom2][1], curResCoords[atom2][2]])
                
        #draw the last O3'-1 - P bond if necessary
        if prevO3 is not None and (dist(prevO3, curResCoords[" P  "]) > BOND_DIST_CUTOFF):
            lineList.append([BATONCOLOR, 6, prevO3[0], prevO3[1], prevO3[2], curResCoords[" P  "][0], curResCoords[" P  "][1], curResCoords[" P  "][2]])
        
        #actually draw the lines
        if len(lineList):
            self.__extraBondObject = new_generic_object_number("Extra bonds")
            set_display_generic_object(self.__extraBondObject, 1)
            for curLine in lineList:
                to_generic_object_add_line(self.__extraBondObject, *curLine)
        graphics_draw()
        
    
    def clearExtraBonds(self):
        """Clears any extra bonds drawn by drawExtraBonds
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            Clears all bonds in drawing object self.__extraBondObject
        """
        
        if self.__extraBondObject is not None:
            close_generic_object(self.__extraBondObject)
        self.__extraBondObject = None
        
    
    def setExtraBondRange(self, startingResNum, endingResNum):
        """Set the nucleotides that will have extra bonds drawn during any subsequent drawExtraBonds() call
        
        ARGUMENTS
            startingResNum - the residue number of the first residue to draw bonds for
            endingResNum   - the residue number of the lsat residue to draw bonds for
        RETURNS:
            None
        """
        
        #convert the residue numbers to indices
        startingIndex = self.resIndex(startingResNum)
        endingIndex   = self.resIndex(endingResNum) + 1
            #we add 1 to the ending index because the Python slice operator is exclusive of the final index, so we need a +1 to include it
        
        self.__extraBondStartingIndex = startingIndex
        self.__extraBondEndingIndex   = endingIndex
    
    def clearExtraBondRange(self):
        """Make any subsequent drawExtraBonds() call draw extra bonds for the entire molecule
        
        ARGUMENTS
            None
        RETURNS:
            None
        """
        
        self.__extraBondStartingIndex = None
        self.__extraBondEndingIndex   = None
    
    
    def resIndex(self, resNum):
        """Find the index (i.e. position in the chain) of a given residue number
        
        ARGUMENTS:
            resNum - the residue number to find (which may include an insertion code)
        RETURNS:
            the index of the specified residue number in this chain
        """
        
        return self.__resNumDict[str(resNum)]
        
    
    def resNum(self, resIndex):
        """Given a residue index, retrieves the residue number and insertion code
        
        ARGUMENTS:
            resIndex - the residue index
        RETURNS:
            the residue number (without the insertion code)
            the insertion code
        """
        
        res = self.__molecule[0][self.__chainIndex][1][resIndex]
        return (res[0], res[1])
    
    
    def resNumFull(self, resIndex):
        """Given a residue index, retrieves the full residue number (i.e. including insertion code)
        
        ARGUMENTS:
            resIndex - the residue index
        RETURNS:
            the residue number as a string including the insertion code
        """
        
        res = self.__molecule[0][self.__chainIndex][1][resIndex]
        return (str(res[0]) + str(res[1])).strip()
        
    
    def __initializeResNumDict(self):
        """initialize the resNum lookup dictionary if this object was created with an initial residue list
        
        ARGUMENTS:
            None
        RETURNS:
            None
        """
        
        resNumList = [str(res[0])+str(res[1]) for res in self.__molecule[0][self.__chainIndex][1]]
        self.__resNumDict = dict([(resNum, resIndex) for (resIndex, resNum) in enumerate(resNumList)])
        #print "resNumDict = ", self.__resNumDict
    
    
    def saveBfacs(self, startingRes, endingRes):
        """Save the B-factors for the specified nucleotides.
        
        ARGUMENTS:
            startingRes - the first residue to save the B-factors for
            endingRes   - the last residue to save the B-factors for
        RETURNS:
            None
        """
        
        startingIndex = self.resIndex(startingRes) #convert startingRes (Coot numbering) to an index (internal numbering)
        endingIndex   = self.resIndex(endingRes)   #convert endingRes (Coot numbering) to an index (internal numbering)
        
        
        bfacs = {}
        for curResIndex in xrange(startingIndex, endingIndex+1):
            resNum = self.resNumFull(curResIndex)
            bfacs[resNum] = {}
            for curAtom in self.__molecule[0][self.__chainIndex][1][curResIndex][3]:
                atomName = curAtom[0][0].strip()
                bfacs[resNum][atomName] = curAtom[1][1]
        
        self.__savedBfacs = [startingIndex, endingIndex, bfacs]
    
    
    def restoreSavedBfacs(self, startingRes = None, endingRes = None):
        """Restore the B-factors saved using saveBfacs()
        
        OPTIONAL ARGUMENTS:
            startingRes - the first residue to restore the B-factors for
                          if not given, B-factors will be restored starting from the first residue saved
            endingRes   - the last residue to restore the B-factors for
                          if not given, B-factors will be restored to the last residue saved
        RETURNS:
            None
        NOTE:
            This function will *not* clear the saved B-factors.
        """
        
        #figure out the first residue to restore
        if startingRes:
            startingIndex = self.resIndex(startingRes)
        else:
            startingIndex = self.__savedBfacs[0]
        
        #figure out the last residue to restore
        if endingRes:
            endingIndex = self.resIndex(endingRes)
        else:
            endingIndex = self.__savedBfacs[1]
            
        bfacs = self.__savedBfacs[2]
        
        #go through each residue and restore B-factors
        for curResIndex in xrange(startingIndex, endingIndex+1):
            resNum = self.resNumFull(curResIndex)
            for curAtom in self.__molecule[0][self.__chainIndex][1][curResIndex][3]:
                atomName = curAtom[0][0].strip()
                
                #it's possible that the molecule we're rotamerizing was missing atoms, so we need to make sure
                #that a given atom had a B-factor to begin with
                if bfacs[resNum].has_key(atomName):
                    curAtom[1][1] = bfacs[resNum][atomName]
                #if an atoms not in bfacs, then just leave it as the default
                #another option would be to set its B-factor to the average of all the other atoms in the nucleotide
                #but we don't do that here
        
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
    
    
    def hasSavedBfacs(self):
        """Does this object have B-factors saved using saveBfacs()?
        
        ARGUMENTS:
            None
        RETURNS:
            True if there are saved B-factors, False otherwise
        """
        
        return (self.__savedBfacs is not None)
    
    
    def clearSavedBfacs(self):
        """Clear coordinates saved using saveBfacs()
        
        ARGUMENTS:
            None
        RETURNS:
            None
        """
        
        self.__savedBfacs = None
    
    
    def saveCoordinates(self, startingResNum, endingResNum):
        """Save coordinates of the specified residues so they can be restored later
        
        ARGUMENTS:
            startingRes - the first residue to save the structure of
            endingRes   - the last residue to save the structure of
        RESTURNS:
            None
        """
        
        self.__savedCoordinates = [startingResNum, endingResNum, self.getCootNucs(startingResNum, endingResNum)]
    
    
    def restoreSavedCoordinates(self):
        """Restore coordinates saved using saveCoordinates
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            Will clear the saved coordinates after restoring them (i.e. subsequent calls to hasSavedCoordinates will return False)
        """
        
        (startingResNum, endingResNum, cootNucs) = self.__savedCoordinates
        self.setCootNucs(startingResNum, endingResNum, cootNucs)
        self.__savedCoordinates = None
        
    
    def hasSavedCoordinates(self):
        """Does this object have coordinates saved using saveCoordinates()?
        
        ARGUMENTS:
            None
        RETURNS:
            True if there are saved coordinates, False otherwise
        """
        
        return (self.__savedCoordinates is not None)
    
    
    def clearSavedCoordinates(self):
        """Clear coordinates saved using saveCoordinates()
        
        ARGUMENTS:
            None
        RETURNS:
            None
        """
        
        self.__savedCoordinates = None
    
    
    def connectedToNextFromIndex(self, resIndex, distanceCutoff = CONNECTED_DISTANCE_CUTOFF):
        """Determine if a nucleotide is connected to the next nucleotide
        
        ARGUMENTS:
            resIndex - the residue index
        OPTIONAL ARGUMENTS:
            distanceCutoff - how long can the O3'-P bond be before the nucleotides are considered not connected
        RETURNS:
            True if the residue is connected to the next one
            False otherwise
        """
        
        #make sure this isn't the last residue
        if resIndex + 1 == self.getNumNts():
            return False
        
        #get the O3' coordinates of the current nucleotide
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            if curAtom[0][0].strip() == "O3'":
                o3Coords = curAtom[2]
                break
        else:
            #if we didn't find the O3' coordinates, then assume that the nucleotides aren't connected
            return False
        
        #get the phosphate coordinates of the next nucleotide
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex+1][3]:
            if curAtom[0][0].strip() == "P":
                phosCoords = curAtom[2]
                break
        else:
            #if we didn't find the phosphate coordinates, then assume that the nucleotides aren't connected
            return False
        
        #make sure that the next residue doesn't have an insertion code (since we can't pass insertion codes to the minimizer)
        if self.__molecule[0][self.__chainIndex][1][resIndex+1][1] != "":
            return False
        
        if dist(o3Coords, phosCoords) > distanceCutoff:
            return False
        else:
            return True
    
    
    def connectedToPrevFromIndex(self, resIndex, distanceCutoff = CONNECTED_DISTANCE_CUTOFF):
        """Determine if a nucleotide is connected to the previous nucleotide
        
        ARGUMENTS:
            resIndex - the residue index
        OPTIONAL ARGUMENTS:
            distanceCutoff - how long can the O3'-P bond be before the nucleotides are considered not connected
        RETURNS:
            True if the residue is connected to the previous one
            False otherwise
        """
        
        if resIndex == 0:
            return False
        else:
            return self.connectedToNextFromIndex(resIndex-1, distanceCutoff)
    
    
    def connectedToNext(self, resNum, distanceCutoff = 2.5):
        """Determine if a nucleotide is connected to the next nucleotide
        
        ARGUMENTS:
            resNum - the residue number
        OPTIONAL ARGUMENTS:
            distanceCutoff - how long can the O3'-P bond be before the nucleotides are considered not connected
        RETURNS:
            True if the residue is connected to the next one
            False otherwise
        """
        
        resIndex = self.resIndex(resNum)
        self.connectedToNextFromIndex(resIndex, distanceCutoff)
        
        
    def connectedToPrev(self, resNum, distanceCutoff = 2.5):
        """Determine if a nucleotide is connected to the next nucleotide
        
        ARGUMENTS:
            resNum - the residue number
        OPTIONAL ARGUMENTS:
            distanceCutoff - how long can the O3'-P bond be before the nucleotides are considered not connected
        RETURNS:
            True if the residue is connected to the previous one
            False otherwise
        """
        
        resIndex = self.resIndex(resNum)
        self.connectedToPrevFromIndex(resIndex, distanceCutoff)
    
    
    def checkPhosAndGlycosidicFromIndex(self, resIndex):
        """Make sure that the specified nucleotide contains phosphate and glycosidic bond coordinates
        
        ARGUMENTS:
            resIndex - the residue index
        RETURNS:
            None if the nucleotide contains the required atoms
            atom name of a missing atom otherwise
        """
        
        resType = self.__molecule[0][self.__chainIndex][1][resIndex][2]
        if resType == "C" or resType == "U":
            baseN = "N1"
        elif resType == "G" or resType == "A":
            baseN = "N9"
        else:
            return "N1/9"
        
        #make a set containing all of the atom names
        atomList = [atom[0][0].strip() for atom in self.__molecule[0][self.__chainIndex][1][resIndex][3]]
        atomSet = frozenset(atomList)
        
        #make sure all of the required atoms are in the set
        for curAtom in ("P", "C1'", baseN):
            if curAtom not in atomSet:
                return curAtom
        else:
            return None
    
    
    def checkPDB2FromIndex(self, resIndex):
        """Find out if a given residue using PDB2-style naming
        
        ARGUMENTS:
            resIndex - the residue index to check
        RETURNS:
            True if this residue uses PDB2-style naming
            False otherwise
        NOTES:
            Currently, this function simply checks for a C1* atom.  If it exists, then the residue
            is assumed to use PDB2 naming.  Otherwise, the residue is assumed to be okay
        """
        
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            if curAtom[0][0].strip() == "C1*":
                return True
        else:
            return False
    
    
    def numAtomsFromIndex(self, resIndex):
        """Return the number of atoms within a given residue (useful for deciding if a residue is just a phosphate)
        ARGUMENTS
            resIndex - the residue index to check
        RETURNS:
            the number of atoms in residue resIndex
        """
        
        return len(self.__molecule[0][self.__chainIndex][1][resIndex][3])
    
    
    def calcSuiteTorsions(self, resNum):
        """Calculate all torsions for the specified suite.
        
        ARGUMENTS:
            resNum - the suite number to calculate torsions for
        RETURNS
            a list containing the seven suite torsions
            note that missing suite atoms will result in None entries for the undefined torsions
        """
        
        resIndex = self.resIndex(resNum)
        return self.calcSuiteTorsionsFromIndex(resIndex)
    
    
    def calcSuiteTorsionsFromIndex(self, resIndex):
        """Calculate all torsions for the specified suite.
        
        ARGUMENTS:
            resIndex - the suite index to calculate torsions for
        RETURNS
            a list containing the seven suite torsions
            note that missing suite atoms will result in None entries for the undefined torsions
        """
        
        #make sure that there actually is a previous residue to measure torsions in
        if resIndex == 0 or not self.connectedToPrevFromIndex(resIndex):
            return None
        
        chainObj = self.createPartialChainObjectFromIndex(resIndex-1, resIndex)
        
        return [chainObj.nucs[0].delta(),
                chainObj.nucs[0].epsilon(),
                chainObj.nucs[0].zeta(),
                chainObj.nucs[1].alpha(),
                chainObj.nucs[1].beta(),
                chainObj.nucs[1].gamma(),
                chainObj.nucs[1].delta()]
    
    
    def __fixTerAtoms(self):
        """Work around a bug in Coot's python_representation function where any TER lines in the PDB file get converted to an unnamed atom with no occupancy at coordinates (0,0,0)
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            Removes any "blank" atoms from self.__molecule
        """
        
        #Check the last atom of each residue of each chain to see if it's a blank atom
        #This assumes that any TER records are at the end of a residue, which is probably
        #a safe assumption (although it's not guaranteed).
        #We can't assume that a TER record will only be in the last residue of a chain, since there may
        #be solvent atoms listed at the end of a chain
        
        #Note that this function removes the blank atom, rather than replacing it with a TER record.
        #Replacing the atom with a TER record would be preferable, but I don't think that there's any way to do that.
        
        
        for curChain in self.__molecule[0]:
            for curRes in curChain[1]:
                finalAtom = curRes[3][-1]
                if finalAtom[0][0] == "" and finalAtom[1][0] == 0.0 and finalAtom[2][0] == 0.0 and finalAtom[2][1] == 0.0 and finalAtom[2][1] == 0.0:
                    #if the final atom of the residue has no name, an occupancy of 0, and coordinates of all 0's, then remove it
                    curRes[3].pop()
    
    
    def __checkAniso(self):
        """Make sure that we're not going to have any problems with anisotropic temperature records
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EXCEPTIONS:
            raises a PseudoMoleculeError if Coot is older than revision 3631 and self.__molecule contains
            anisotropic temperature records, since clear_and_update_molecule() won't be able to handle
            reading the molecule back in
        """
        
        #this function can be removed if we ever require Coot rev 3631 or newer to run RCrane
        if svn_revision() >= 3631:
            return
        
        for curChain in self.__molecule[0]:
            for curRes in curChain[1]:
                for curAtom in curRes[3]:
                    if isinstance(curAtom[1][1], list):
                        raise PseudoMoleculeError("This molecule contains anisotropic temperature records.  RCrane requires Coot 0.7-pre r3631 or newer to rotamerize molecules with anisotropic temperature records.")
    
    
    def getAtomNames(self, resNum, strip = True):
        """Get a list of all the atom names in a specific residue
        
        ARGUMENTS:
            resNum    - the residue number to get atom names from
        OPTIONAL ARGUMENTS:
            strip     - whether to remove leading and trailing spaces from the atom names
                        defaults to True
        RETURNS:
            atomNames - a list of atom names
        """
        
        #convert the residue number to an index
        resIndex = self.resIndex(resNum)
        
        #get a list of all the atom names
        atomNames = [curAtom[0][0] for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]]
        
        #strip spaces from the atom names if desired
        if strip:
            return [curAtom.strip() for curAtom in atomNames]
        else:
            return atomNames
        
    def resTypeFromIndex(self, resIndex):
        """Get the residue type (i.e. A, G, C, or U) of the specified residue
        
        ARGUMENTS:
            resIndex - the residue index
        RETURNS:
            the residue type of the specified residue
        """
        
        return str(self.__molecule[0][self.__chainIndex][1][resIndex][2])
        #the residue type should always be a string.  The explicit casting here should only matter when something has already gone wrong.
    
    
    #initialize a dictionary with the desired atom order for __reorderAtoms
    __reorderAtomsOrder = {}
    __reorderAtomsOrder["A"] = dict(zip("P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 N6 N1 C2 N3 C4".split(" "), range(40)))
    __reorderAtomsOrder["G"] = dict(zip("P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4".split(" "), range(40)))
    __reorderAtomsOrder["C"] = dict(zip("P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 N4 C5 C6".split(" "), range(40)))
    __reorderAtomsOrder["U"] = dict(zip("P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' N1 C2 O2 N3 C4 O4 C5 C6".split(" "), range(40)))
    
    def finalMoleculeCleanup(self, fixSegids = False):
        """Perform a final cleanup on the coot molecule, consisting of:
            Reordering the newly created atoms so that they follow the standard ordering for RNA residues
            Restoring the segid for all newly created atoms (if requested)
        
        ARGUMENTS:
            None
        OPTIONAL ARGUMENTS:
            fixSegids - whether to properly set segid (segment ID, such as used in CNS) for all newly created atoms
                Defaults to False
                Note that this is only required when doing a rotamerize, not when creating a new trace,
                since new molecules created by RCrane always have a blank segid
        RETURNS:
            None
        NOTE:
            If the extra bond range has been set via setExtraBondRange, then only nucleotides in that range will
            be cleaned up.
        """
        
        startIndex = None
        endIndex = None
        if self.__extraBondStartingIndex is not None and self.__extraBondEndingIndex is not None:
            startIndex = self.__extraBondStartingIndex
            endIndex   = self.__extraBondEndingIndex
        
        for curRes in self.__molecule[0][self.__chainIndex][1][startIndex:endIndex]:
            #sort all the atoms by name, using the ordering in __reorderAtomsOrder for the current residue type
            curRes[3].sort(key = lambda x: self.__reorderAtomsOrder[str(curRes[2])][x[0][0].strip()])
        
        #if this residue previously had a segid, apply it to all the new atoms
        for curRes in self.__molecule[0][self.__chainIndex][1][startIndex:endIndex]:
            #figure out the segid of the phosphate
            #if we've done anything to this nucleotide, it must have had a phosphate.
            #We just sorted the atoms, so there must be a phosphate as the first atom
            desiredSegid = curRes[3][0][1][3]
            
            #if there is a segid, apply it to all the atoms
            if desiredSegid != '':
                for curAtom in curRes[3]:
                    curAtom[1][3] = desiredSegid
        
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
        
    
    def isOnlyPhosGroupFromIndex(self, resIndex):
        """Determine if the specified residue contains no more than a phosphate group (i.e. P, OP1, OP2, and/or OP3)
        
        ARGUMENTS:
            resIndex - the residue index
        RETURNS:
            True if the specified residue is only a phosphate group (or a subset of a phosphate group)
            False otherwise
        """
        
        atomList = self.__molecule[0][self.__chainIndex][1][resIndex][3]
        #from pprint import pprint; pprint(atomList)
        
        #if we have more than four atoms, then it must not be a phosphate group
        if len(atomList) > 4:
            #print "not phosphate group because atomList too long"
            return False
        
        #if we have four or fewer atoms, check their names
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            atomName = curAtom[0][0].strip()
            if atomName not in PHOS_GROUP_ATOMS:
                #print "not phosphate group because of", atomName
                return False
        else:
            #print "is phosphate group"
            return True
    
    
    def isOnlyPhosGroup(self, resNum):
        """Determine if the specified residue contains no more than a phosphate group (i.e. P, OP1, OP2, and/or OP3)
        
        ARGUMENTS:
            resNum - the residue number
        RETURNS:
            True if the specified residue is only a phosphate group (or a subset of a phosphate group)
            False otherwise
        """
        
        resIndex = self.resIndex(resNum)
        return self.isOnlyPhosGroupFromIndex(resIndex)
    
    
    def getPhosCoordsFromIndex(self, resIndex):
        """Get the coordinates of a specified phosphate.
        
        ARGUMENTS:
            resIndex - the residue index
        RETURNS:
            the coordinates of the phoshate from the specified residue
        """
        
        for curAtom in self.__molecule[0][self.__chainIndex][1][resIndex][3]:
            if curAtom[0][0].strip() == "P":
                return curAtom[2]
        
        #if we haven't found a phosphate atom
        return None
    
    
    def setResInsertionPoint(self, resIndex):
        """Set the residue insertion point (i.e. where in the chain new nucleotides should be added).
        This is only necessary when doing an extend chain.
        
        ARGUMENTS:
            resIndex - the residue index to use as the residue insertion point
        RETURNS:
            None
        NOTE:
            This sets both resInsertionPoint and origResInsertionPoint
        """
        
        self.resInsertionPoint = resIndex
        self.origResInsertionPoint = resIndex
    
    
    def saveMoleculeState(self):
        """Store the state of the entire molecule (can be restored via restoreMoleculeState)
        
        ARGUMENTS:
            None
        RETURNS:
            None
        """
        
        #TODO: replace this with saveCoordinates and removing elements from the residue list so that we don't waste memory storing an entire extra copy of
        #the molecule (which could be large if it's a ribosome)
        
        self.__savedMoleculeState = deepcopy(self.__molecule[0][self.__chainIndex][1])
    
    
    def restoreMoleculeState(self):
        """Restore the molecule state that was previously saved using saveMoleculeState()
        
        ARGUMENTS:
            None
        RETURNS:
            None
        """
        
        if self.__savedMoleculeState is not None:
            self.__molecule[0][self.__chainIndex][1] = self.__savedMoleculeState
            clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
            return True
        else:
            return False
    
    
    def hasSavedMoleculeState(self):
        """Determine if this object has a saved molecule state (stored via saveMoleculeState())
        
        ARGUMENTS:
            None
        RETURNS:
            True if there is a saved molecule state
            False otherwise
        """
        
        if self.__savedMoleculeState is None:
            return False
        else:
            return True
    
    
    def clearNonBridgingOxysFromIndex(self, resIndex):
        """Delete non-bridging oxygens from the specified residue
        
        ARGUMENTS:
            resIndex - the residue index to clear non-bridging oxygens from
        RETURNS:
            None
        """
        
        atomList = self.__molecule[0][self.__chainIndex][1][resIndex][3]
        newAtomList = [curAtom for curAtom in atomList if curAtom[0][0].strip() == "P"]
        self.__molecule[0][self.__chainIndex][1][resIndex][3] = newAtomList
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
     
    
    def deleteResFromIndex(self, resIndex):
        """Delete an entire residue
        
        ARGUMENTS:
            resIndex - the index number of the residue to delete
        RETURNS:
            None
        NOTES:
            No redraw will be done
        """
        
        del self.__molecule[0][self.__chainIndex][1][resIndex]
        self.__initializeResNumDict()
        clear_and_update_molecule(self.__moleculeNumber, self.__molecule)
    
    
    def mergeRes(self, phosGroupIndex, res2Index):
        """Merge the atom lists of two residues
        
        ARGUMENTS:
            phosGroupIndex - the first residue to merge.  Presumably, this contains only a phosphate group
            res2Index      - the second residue to merge
        RETURNS:
            None
        NOTES:
            The atoms from phosGroup will be put before the atoms of res2, but the residue information
            (i.e. residue type) from res2 will be kept.
        """
        
        #prepend the phosGroup atoms onto res2
        phosGroupAtoms = self.__molecule[0][self.__chainIndex][1][phosGroupIndex][3]
        #res2Atoms      = self.__molecule[0][self.__chainIndex][1][res2Index][3]
        self.__molecule[0][self.__chainIndex][1][res2Index][3][0:0] = phosGroupAtoms
        
        self.deleteResFromIndex(phosGroupIndex)


#moleculeNameCount stores the number that we should append to the molecule name, i.e., RCrane Molecule 1, RCrane Molecule 2, RCrane Molecule 3...
#it's implemented as a simple generator function
def incGenerator():
    """Create a generator that starts at 1 and increments by 1"""
    num = 1
    while True:
        yield num
        num += 1
moleculeNameCount = incGenerator()


class PseudoMoleculeError(Exception):
    pass
