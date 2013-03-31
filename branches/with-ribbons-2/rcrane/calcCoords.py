#!/usr/bin/env python
"""Functions for calculating nucleotide coordinates"""

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

import os.path
import re
import gtk
from copy import deepcopy
#from time import sleep      #for debugging
#from pprint import pprint   #for debugging

from coot import refinement_immediate_replacement_state, set_refinement_immediate_replacement, accept_regularizement, clear_all_fixed_atoms, add_extra_torsion_restraint, set_refine_with_torsion_restraints, refine_with_torsion_restraints_state, matrix_state, set_matrix, delete_all_extra_restraints, add_extra_start_pos_restraint, refine_zone, set_use_only_extra_torsion_restraints_for_torsions
from coot import svn_revision #needed to see if Coot is new enough for Phenix restraints
from coot import monomer_restraints_py           as monomer_restraints
from coot import set_monomer_restraints_py       as set_monomer_restraints
from coot import refine_residues_py              as refine_residues
from coot import refine_zone_with_score_py       as refine_zone_with_score
from coot import mark_multiple_atoms_as_fixed_py as mark_multiple_atoms_as_fixed

#regularize_zone_with_score_py only exists in Coot newer than 3728 (0.7-pre)
#so we won't be able to find it in Coot 0.6.2
#this function is only used in Rotamerize without density, and the menu option for that won't be created
#unless Coot is newer than 3728
#so we can safely ignore the ImportError
try:
    from coot import regularize_zone_with_score_py   as regularize_zone_with_score
except ImportError:
    pass


#use_only_extra_torsion_restraints_for_torsions_state() only exists in Coot newer than ~3902
#if the function doesn't exist, then just assume that this variable was turned off
#(since RCrane is probably the only thing that uses it)
try:
    from coot import use_only_extra_torsion_restraints_for_torsions_state
except ImportError:
    def use_only_extra_torsion_restraints_for_torsions_state(): return 0

from buildInitSugar import BuildInitSugar, rotateSugar
from puckerList import puckerList
from buildPhosOxy import buildPhosOxy, buildInitOrTerminalPhosOxy
from rotData import RotDataLoader
from puckerList import puckerList
from guiUtils import HBOX_SPACING, VBOX_SPACING, createRCraneWindowObject
import phenixRestraints


#initialize a BuildInitSugar object when this module is loaded
rcranePath = os.path.dirname(os.path.abspath(__file__))
dataPath = os.path.join(rcranePath, "data")
sugarBuilder = BuildInitSugar(c3pStruc = os.path.join(dataPath, "c3p.pdb"),
                              c2pStruc = os.path.join(dataPath, "c2p.pdb"))

#initialize a rotData object when this module is loaded
rotData = RotDataLoader(os.path.join(dataPath, "dihedData.csv"))

TORSION_STD_DEV_MOD = 0.1  #for minimization restraints, all the standard deviations for all rotamer torsions are multiplied by this number
MANUAL_TORSION_STD_DEV = 2 #when rotamerizing an already existing structure, this will be used as the standard deviation for the
                           #non-predicted torsions (i.e. the first and last half-nucleotides)

#the standard deviations for the harmonic (start position) restraints
#the larger the number, the more freedom the restrained atoms have to move
HARMONIC_STD_DEV_PHOSPHATE = 0.25 #for phosphates
HARMONIC_STD_DEV_BASE      = 0.1  #for base atoms

#sugar torsions for C3' and C2' sugars (taken from Phenix (phenix/1.6.2-432/phenix-1.6.2-432/chem_data/geostd/rna_dna/mod_rna2p.cif and mod_rna3p.cif))
NU0_C3 = 3.0
NU1_C3 = 335.0
NU4_C3 = 145.0      #Coot refers to C5'-C4'-O4'-C1' as nu4.  I think nu4 is actually C3'-C4'-O4'-C1'.
                    #We're using the C5' version of the torsion here

NU0_C2 = 339.0
NU1_C2 = 35.0
NU4_C2 = 123.0

#TODO: use realnu4 and wasnu4? that's what Phenix uses

#torsions for chi (also taken from Phenix)
CHI_MEAN    = -123.0
CHI_STD_DEV = 24.3

NU_STD_DEV = 4 #Coot default is 40, but we want to make sure that the sugars don't get flattened
#TODO: try a value of 8 here? that's what Phenix uses


REFINE_MAP_WEIGHT = 10  #the weight of the map term during the minimization (Coot default is 60)
    #if this value it too high, Coot will distort the sugar to try to fit the O2' into density
    #TODO: try to balance this value against NU_STD_DEV

SYN_REBUILDING_CUTOFF = 3 #any minimization scores above this will cause the minimiztion to be restarted using a syn sugar
HIGH_ANTI_REBUILDING_CUTOFF = 8 #if both the anti and the syn minimizations are above this score, then the minimization
    #will be restarted using a high-anti sugar
REFINEMENT_FAIL_SCORE = 99999 #what score to assign a minimization that didn't return a score (which means that the refinement failed)
PRINT_SUMMARY_TABLE = False #whether to print a summary table of minimization scores
    #note that setting rcrane_debug to True before launching RCrane will set PRINT_SUMMARY_TABLE to True
    #(launch.py will set PRINT_SUMMARY_TABLE after loading this module)

#atoms to not fix during minimization
PREV_RES_MOBILE_ATOMS = frozenset(["O3'"]) #OP1 and OP2 may be added at runtime
CUR_RES_MOBILE_ATOMS  = frozenset("P C2' O2' C3' O3' C4' O4' C5' O5'".split(" ")) #C1' should be the only restrained backbone atom
NEXT_RES_MOBILE_ATOMS  = frozenset("P O5' OP1 OP2".split(" ")) #The O5' will only be present if the next residue is already built

#default to using the Coot/CCP4 restraints, as opposed to Phenix's pucker-specific restraints
#This default may change in the future
USE_PHENIX_RESTRAINTS = False


PHENIX_NEW_RESTRAINT_ATOMS = frozenset("N1 N9 C1' C2' O2' C3' O3' C4' O4' C5' O5' P OP1 OP2".split(" "))
    #bond and angle restraints where all atoms are on this list will be rewritten if USE_PHENIX_RESTRAINTS is True

def calcCoords(builtChain, bestPath, pseudoMol, window):
    """Calculate coordinates for a chain of nucleotides
    
    ARGUMENTS:
        builtChain - a chain object containing phosphate and base coordinates for all nucleotides to be built
        bestPath   - a list of the desired conformers to be built
        pseudoMol  - a pseudoMolecule object currently being used to display the chain
        window     - the window contianing the GUI
    RETURNS:
        intermediateAtomLocs - phosphate and O3' coordinates for each nucleotide immediately before minimization of that nucleotide was started
        minimizationScores   - a list of the minimization scores for each nucleotide
    NOTE:
        If only a single nucleotide is being built, then we don't have a full suite, so there are no confmers
        In this case, bestPath should be a scalar containing the intended pucker of the nucleotide to be built
    """
    
    #put a progress bar in window
    progressDialog = ProgressDialogObject(window, builtChain.numNucs()-1)
    
    #we're going to have to change some Coot settings during the coordinate calculation
    #so check what the existing values are so we can set them back
    origCootSettings = __changeCootSettings()
    
    intermediateAtomLocs = None
    minimizationScores = None
    
    #enclose the rest of the function in a try clause so that we can still reset the Coot variables even if something goes wrong
    try:
    #for x in (1,):
        if builtChain.numNucs() == 1:
            #if there's only one nucleotide, then there's nothing we can do (the nucleotide should contain only a phosphate)
            return
        
        elif builtChain.numNucs() == 2:
            #if there are two nucleotides, then we have a single sugar with both a 5' and a 3' phosphate
            #we can't determine a conformer, but we can predict the sugar pucker and minimize things without any torsions
            
            __minCoords(pseudoMol,                                                                          #pseudoMolecule object
                        None, None, 1, 2,                                                                   #residue numbers
                        None, None,                                                                         #rotamers
                        None, builtChain.nucleotides[0].type, None,                                         #residue types
                        builtChain.nucleotides[0].atoms,                                                    #atomic coordinates
                        pucker = bestPath, nextResAtoms = builtChain.nucleotides[1].atoms)
        
        else:
            
            #if we recalculate coordinates later, we want to mimic the conditions of this minimization as closely as possible
            #this means we need to store intermediate locations for the phophate and O3'
            #(i.e. for phosphate i, we need to store it's location after minimizing nucleotide i-1 but before minimizing nucleotide i)
            #the first phosphate doesn't have an intermediate location, so just store None for that nucleotide
            intermediateAtomLocs = [None]
            minimizationScores = []
            
            #build the first nucleotide (don't put torsional constraints on the initial alpha, beta, and gamma)
            
            #minimize the structure
            (newCoords, score) = __minCoords(pseudoMol,                                                  #pseudoMolecule object
                                             None, None, 1, 2,                                           #residue numbers
                                             None, bestPath[0],                                          #rotamers
                                             None, builtChain.nucs[0].type, builtChain.nucs[1].type,     #residue types
                                             builtChain.nucs[0].atoms)                                   #atomic coordinates
            #update the builtChain object with the new coordinates
            (builtChain.nucs[0].atoms, builtChain.nucs[1].atoms) = newCoords
            #store the phosphate and previous O3' location
            intermediateAtomLocs.append([builtChain.nucs[1].atoms["P"], builtChain.nucs[0].atoms["O3'"], builtChain.nucs[0].atoms["C3'"]])
            minimizationScores.append(score)
            #increment the progress bar
            progressDialog.progress()
            
            #return
            
            fixPrevPhosOxy = True #don't minimize the first previous phosphoryl oxygens since it's going to take a long time and isn't
                #going to make them any more accurate since there's no O3' to use to position them
            
            #built the middle nucleotides
            for resNum in xrange(1, len(bestPath)):
                
                #minimize the structure
                (newCoords, score) = __minCoords(pseudoMol,                                                                             #pseudoMolecule object
                                        resNum-1 or None, resNum, resNum+1, resNum+2,                                                   #residue numbers
                                        bestPath[resNum-1], bestPath[resNum],                                                           #rotamers
                                        builtChain.nucs[resNum-1].type, builtChain.nucs[resNum].type, builtChain.nucs[resNum+1].type,   #residue types
                                        builtChain.nucs[resNum].atoms,                                                                  #atomic coordinates
                                        fixPrevPhosOxy = fixPrevPhosOxy)
                #update the builtChain object with the new coordinates
                (builtChain.nucs[resNum-1].atoms, builtChain.nucs[resNum].atoms, builtChain.nucs[resNum+1].atoms) = newCoords
                #store the phosphate location
                intermediateAtomLocs.append([builtChain.nucs[resNum+1].atoms["P"], builtChain.nucs[resNum].atoms["O3'"], builtChain.nucs[resNum].atoms["C3'"]])
                minimizationScores.append(score)
                progressDialog.progress() #increment the progress bar
                fixPrevPhosOxy = False    #minimize all of the non-bridging oxygens from here on
                
            #build the last nucleotide (don't put torsional constraints on the final epsilon and zeta)
            resNum = len(bestPath)
            
            #minimize the structure
            (newCoords, score) = __minCoords(pseudoMol,                                                 #pseudoMolecule object
                                    resNum-1 or None, resNum, resNum+1, resNum+2,                       #residue numbers
                                    bestPath[resNum-1], None,                                           #rotamers
                                    builtChain.nucs[resNum-1].type, builtChain.nucs[resNum].type, None, #residue types
                                    builtChain.nucs[resNum].atoms,                                      #atomic coordinates
                                    fixPrevPhosOxy = fixPrevPhosOxy)
            #update the builtChain object with the new coordinates
            (builtChain.nucs[resNum-1].atoms, builtChain.nucs[resNum].atoms, builtChain.nucs[resNum+1].atoms) = newCoords
            minimizationScores.append(score)
            #increment the progress bar before we do the last thing so that the user sees it at 100% for a few seconds
            progressDialog.progress() 
            
            #We no longer minimize the terminal phosphoryl oxygens, since it occasionally takes a long time and doesn't seem
            #to improve their placement much.
    
    finally:
        #restore the original Coot settings even if something went wrong during the minimization
        __restoreCootSettings(origCootSettings)
    
    #only draw extra bonds if the user is going to be presented with a GUI
    #otherwise, there won't be any way to delete the extra bonds
    if builtChain.numNucs() > 2:
        pseudoMol.drawExtraBonds()
    
    return (intermediateAtomLocs, minimizationScores)


def recalcCoords(startingRes, endingRes, rots, origCoords, pseudoMol, window, ignoreDensity =  False):
    """Recalculate coordinates (using different rotamers) for part of a chain of nucleotides starting with a chain built by calcCoords
    
    ARGUMENTS:
        startingRes             - the first residue to rebuild
        endingRes               - the last residue to rebuilt
        rots                    - what rotamers to use for the rebuild
        origCoords              - a chain object containing the current coordinates
        pseudoMol               - a pseudoMolecule object currently being used to display the chain
        window                  - the window contianing the GUI
    OPTIONAL ARGUMENTS:
        ignoreDensity           - ignore the density when performing the minimization
                                  defaults to False
    RETURNS:
        intermediateAtomLocs - phosphate and O3' coordinates for each nucleotide immediately before minimization of that nucleotide was started
        minimizationScores   - a list of the minimization scores for each nucleotide
    NOTE:
        Currently, the return values from this function are only used for initial rotamerize minimization.  They are ignored for all
        other calls to this function
    """
    
    #convert the starting and ending residues to indices
    startingResIndex = origCoords.resIndex(startingRes)
    endingResIndex   = origCoords.resIndex(endingRes)
    #print "Indices:", startingResIndex, ",", endingResIndex
    
    progressDialog = ProgressDialogObject(window, endingResIndex - startingResIndex + 1)
    
    while gtk.events_pending(): gtk.main_iteration(False)
        #at least on Windows, the progressDialog isn't showing until after the atoms are frozen without this line
    
    #we're going to have to change some Coot settings during the coordinate calculation
    #so check what the existing values are so we can set them back
    origCootSettings = __changeCootSettings()
    builtChain = deepcopy(origCoords) #we don't want to modify this object
    
    print "Recalculating ", startingRes, "-", endingRes
    #from time import sleep; sleep(3)
    
    #initialize intermediateAtomLocs (the first nucleotide doesn't have an intermediate location)
    intermediateAtomLocs = [None]
    minimizationScores = []
    
    try:
        
        if startingResIndex == 0:
            #if we're at the first nucleotide of the chain (or we're right after a chain break)
            
            print "Rebuilding initial nucleotide"
            #if we have to build the first nucleotide of the chain
            
            #minimize the structure
            (newCoords, score) = __minCoords(pseudoMol,                                                     #pseudoMolecule object
                                    None, None, startingRes, builtChain.nucs[1].resNum,                     #residue numbers
                                    None, rots[1],                                                          #rotamers
                                    None, builtChain.nucleotides[0].type, builtChain.nucleotides[1].type,   #residue types
                                    builtChain.nucleotides[0].atoms,                                        #atomic coordinates
                                    ignoreDensity = ignoreDensity)
            #update the builtChain object with the new coordinates
            (builtChain.nucleotides[0].atoms, builtChain.nucleotides[1].atoms) = newCoords
            intermediateAtomLocs.append([builtChain.nucs[1].atoms["P"], builtChain.nucs[0].atoms["O3'"], builtChain.nucs[0].atoms["C3'"]])
            minimizationScores.append(score)
            
            startingResIndex += 1
            rotsIndex         = 2
            
            #increment the progress bar
            progressDialog.progress()
        else:
            rotsIndex  = 1
        
        fixPrevPhosOxy = True #don't minimize the first previous phosphoryl oxygens
            #if we just built the first nucleotide, then minimizing them is going to take a long time and isn't going to make them 
            #any more accurate since there's no O3' to use to position them
            #if we didn't just build the first nucleotide, then we don't want to adjust the non-bridging oxygens that are outside
            #of our minimization range
        for resIndex in xrange(startingResIndex, endingResIndex):
            
            
            
            #minimize the structure
            (newCoords, score) = __minCoords(pseudoMol,                                                                                         #pseudoMolecule object
                                    builtChain.nucs[resIndex-2].resNum if resIndex > 1 else None, builtChain.nucs[resIndex-1].resNum, builtChain.nucs[resIndex].resNum, builtChain.nucs[resIndex+1].resNum,   #residue numbers
                                    rots[rotsIndex-1], rots[rotsIndex],                                                                         #rotamers
                                    builtChain.nucs[resIndex-1].type, builtChain.nucs[resIndex].type, builtChain.nucs[resIndex+1].type,         #residue types
                                    builtChain.nucs[resIndex].atoms,                                                                            #atomic coordinates
                                    fixPrevPhosOxy = fixPrevPhosOxy,                                                                            #if this is the first residue that we're rebuilding but not the first residue
                                                                                                                                                #of the chain, then don't move the non-bridging oxygens of the previous residue
                                    ignoreDensity = ignoreDensity)
            #update the builtChain object with the new coordinates
            (builtChain.nucleotides[resIndex-1].atoms, builtChain.nucleotides[resIndex].atoms, builtChain.nucleotides[resIndex+1].atoms) = newCoords
            intermediateAtomLocs.append([builtChain.nucs[resIndex+1].atoms["P"], builtChain.nucs[resIndex].atoms["O3'"], builtChain.nucs[resIndex].atoms["C3'"]])
            minimizationScores.append(score)
            #increment the progress bar
            progressDialog.progress()
            
            rotsIndex += 1
            firstRes  = False
            fixPrevPhosOxy = False #minimize all of the non-bridging oxygens from here on
        
        resIndex = endingResIndex
        
        #for the last rebuilt nucleotide, we have to call __minCoords with nextResAlreadyBuilt=True
        #we don't need to worry about doing anything special if this is the last suite of the chain, though, since the 3' phosphoryl oxygens are already built
        #only do this if we have a final 3' phosphate
        
        #minimize the structure
        if (resIndex + 1) < len(builtChain.nucs):
            #print "*****Minimizing last nt*****"
            #from time import sleep; sleep(3)
            (newCoords, score) = __minCoords(pseudoMol,                                                                                         #pseudoMolecule obj
                                    builtChain.nucs[resIndex-2].resNum if resIndex > 1 else None, builtChain.nucs[resIndex-1].resNum, builtChain.nucs[resIndex].resNum, builtChain.nucs[resIndex+1].resNum,   #residue numbers
                                    rots[rotsIndex-1], rots[rotsIndex],                                                                         #rotamers
                                    builtChain.nucs[resIndex-1].type, builtChain.nucs[resIndex].type, builtChain.nucs[resIndex+1].type,         #residue types
                                    builtChain.nucs[resIndex].atoms,                                                                            #atomic coordinates
                                    nextResAlreadyBuilt=True, fixPrevPhosOxy = fixPrevPhosOxy, ignoreDensity = ignoreDensity)
            #update the builtChain object with the new coordinates
            (builtChain.nucleotides[resIndex-1].atoms, builtChain.nucleotides[resIndex].atoms, builtChain.nucleotides[resIndex+1].atoms) = newCoords
            minimizationScores.append(score)
            #increment the progress bar (even though the user probably won't actually see this, since it will be replaced almost immediately by the review suites GUI))
        progressDialog.progress()
        
        #don't bother to do a separate minimization run for the phosphoryl oxygens of the final nucleotide.  They should be good enough
        
    finally:
        __restoreCootSettings(origCootSettings) #restore the original Coot settings even if something went wrong during the minimization
        progressDialog.restoreWindow() #never return from this function without running restoreWindow()
            #GTK can crash Coot if it tries to update a GUI element that had been removed from the window and not restored (and possibly gone out of scope)
    
    pseudoMol.drawExtraBonds() #redraw overly-long bonds for the entire molecule
    
    return (intermediateAtomLocs, minimizationScores)


def __minCoords(pseudoMol, prevPrevResNum, prevResNum, curResNum, nextResNum, curRot, nextRot, prevResType, curResType, nextResType, curResAtoms, nextResAlreadyBuilt = False, fixPrevPhosOxy = False, pucker = None, nextResAtoms = None, ignoreDensity = False):
    """Build a single nucleotide and minimize its coordinates.
    
    ARGUMENTS:
        pseudoMol           - the pseudoMolecule object currently being used to display the chain
        prevPrevResNum      - the residue number of residue i-2
                              only used if we're using Phenix's pucker-specific restraints, in which case it's needed
                              to restrain the non-bridging oxygens of residue prevResNum
        prevResNum          - the residue number of the previous residue using Coot numbering (i.e. starting at 1, not 0)
        curResNum           - the residue number to build using Coot numbering (i.e. starting at 1, not 0)
        nextResNum          - the residue number of the next residue using Coot numbering (i.e. starting at 1, not 0)
        curRot              - the rotamer to use when building the current residue
        nextRot             - the rotamer to use when building the next residue
        prevResType         - the residue type (i.e. A, G, C, or U) of the previous residue
        curResType          - the residue type (i.e. A, G, C, or U) of the current residue
        nextResType         - the residue type (i.e. A, G, C, or U) of the next residue
        curResAtoms         - a hash of atom coordinates for the current residue
    OPTIONAL ARGUMENTS:
        nextResAlreadyBuilt - should be True if the next residue is already built (i.e. if we are being
                              called from recalcCoords instead of calcCoords)
                              Defaults to False.
        fixPrevPhosOxy      - whether to fix the phosphoryl oxygens of the previous nucleotide during minimization
        pucker              - the sugar pucker for this nucleotide.  Only necessary when both curRot
                              and nextRot are None (i.e. when we're building the only nucleotide of the chain)
                              Defaults to None
        nextResAtoms        - a hash of atoms containing coordinates for the next phosphate.
                              Should only be provided if we are building the only nucleotide of the chain
                              Defaults to None
        ignoreDensity       - ignore the density when performing the minimization
                              defaults to False
    RETURNS:
        newCoords           - new atomic coordinates calculated by the minimization procedure, formatted as a list of dictionaries
        score               - the score received from Coot's minimization
    """
    
    #print "*** Running __minCoords with:"
    #print "*** prevResNum =", prevResNum, ",\tcurResNum =", curResNum, ",\tnextResNum =", nextResNum
    #print "*** curRot =", curRot
    #print "*** nextRot =", nextRot
    #print "*** prevResType =", prevResType, "curResType =", curResType, "nextResType =", nextResType
    #print "*** fixPrevPhosOxy =", fixPrevPhosOxy
    #sleep(2)
    
    chain = pseudoMol.chain
    
    #initialize the refinement score variables so that we can easily print them later without worrying about them being undefined
    antiRefinementScore     = None
    synRefinementScore      = None
    highAntiRefinementScore = None
    selectedStartingStruc   = None
    
    #if we're not using Phenix restraints, then we don't need the i-2 residue number
    if not USE_PHENIX_RESTRAINTS:
        #By setting prevPrevResNum to None, that residue won't be included in the minimization
        #(which is what we want if we're not using Phenix restraints)
        prevPrevResNum = None
    
    #build an anti sugar
    if isinstance(curRot, str):
        #if curRot is None or a list (indicating that we're manually setting torsions), then don't use it to determine the sugar pucker
        curPucker = puckerList[curRot][1]
    elif nextRot is not None:
        curPucker = puckerList[nextRot][0]
    else:
        curPucker = pucker
    initSugarCoords = sugarBuilder.buildSugar(curResAtoms, curPucker)
    pseudoMol.addSugar(curResNum, initSugarCoords)
    
    #if we're building the only nucleotide of the chain, then build the phosphates before minimization
    if nextResAtoms is not None:
        sugarAndP = dict(curResAtoms.items()+initSugarCoords.items())
        phosOxyCoords5 = buildInitOrTerminalPhosOxy(sugarAndP)
        pseudoMol.addPhosOxy(curResNum, phosOxyCoords5)
        phosOxyCoords3 = buildInitOrTerminalPhosOxy(nextResAtoms, sugarAndP)
        pseudoMol.addPhosOxy(nextResNum, phosOxyCoords3)
        
    #store the structure of the previous and next nucleotides, so we can restore them if we restart minimization with a syn sugar
    preMinimizationStruc = pseudoMol.getCootNucs(prevResNum or curResNum, nextResNum)
    
    #uncomment these lines to produce a Coot script for the minimization (i.e a script that will replicate the minimization we're about to do)
    #if curResNum == 8:
    #    import sys
    #    if sys.modules.has_key("genMinimizationScript"): del sys.modules["genMinimizationScript"]
    #    from genMinimizationScript import genMinimizationScript
    #    genMinimizationScript("lastNucMin.py", pseudoMol, prevResNum, curResNum, nextResNum, curRot, nextRot, prevResType, curResType, nextResType, curResAtoms, nextResAlreadyBuilt, fixPrevPhosOxy)
    #    raise Exception #don't actually run the minimization
    
    #set the appropriate restraints
    molNum = pseudoMol.molNum()
    __fixAtoms(pseudoMol, prevResNum, curResNum, nextResNum, fixPrevPhosOxy)
    __setTorsionRestraints(molNum, chain, prevResNum, curResNum, nextResNum, curRot, nextRot, curResType, nextResAlreadyBuilt)
    
    if USE_PHENIX_RESTRAINTS:
        __addPhenixRestraints(molNum, chain, prevPrevResNum, prevResNum, curResNum, nextResNum, curRot, nextRot, curResType, nextResAlreadyBuilt)
        if prevPrevResNum is not None:
            #when USE_PHENIX_RESTRAINTS is true, we remove the polymer-type from RNA nucleotides so that
            #we can overwrite the link restraints (i.e. restraints spanning more than one nucleotide)
            #as a side-effect of this, Coot no longer implicitely includes flanking nucleotides in the minimization
            #as a result, the non-bridging oxygens of prevResNum will be improperly positioned unless we manually
            #include the fully immobilized prevPrevResNum in the minimization
            __fixEntireResidue(molNum, chain, prevPrevResNum, pseudoMol.getAtomNames(prevPrevResNum, strip=False))
    
    #minimize the structure
    antiRefinementScore = __runRefinement(molNum, chain, prevPrevResNum or prevResNum or curResNum, nextResNum, ignoreDensity)
    score = antiRefinementScore
    selectedStartingStruc = "anti"
    
    #check how good the results of the minimization are to decide if we should restart with a syn chi
    if antiRefinementScore > SYN_REBUILDING_CUTOFF:
    #if False:
        #if the minimization did not end well, then we should restart it with a syn sugar
        
        #first, store the results of this minimization
        newCoords = pseudoMol.updateRes(prevResNum or curResNum, nextResNum)
        #antiRefinementStruc = pseudoMol.getCootNucs((resNum-1 or resNum) - 1, resNum)
        antiRefinementStruc = pseudoMol.getCootNucs(prevResNum or curResNum, nextResNum)
        
        #restore the structure as it was before minimization
        #pseudoMol.setCootNucs((resNum-1 or resNum) - 1, resNum, preMinimizationStruc, False)
        pseudoMol.setCootNucs(prevResNum or curResNum, nextResNum, preMinimizationStruc, False)
            #Note: the False fourth argument skips running clear_and_update_molecule, since that will be run by
            #replaceSugar() below
        
        #calculate syn sugar coordinates
        synSugarCoords = rotateSugar(initSugarCoords, curResAtoms, "syn")
        pseudoMol.replaceSugar(curResNum, synSugarCoords)
        
        #uncomment these lines to produce a Coot script for the minimization (i.e a script that will replicate the minimization we're about to do)
        #if curResNum == 8:
        #    import sys
        #    if sys.modules.has_key("genMinimizationScript"): del sys.modules["genMinimizationScript"]
        #    from genMinimizationScript import genMinimizationScript
        #    genMinimizationScript("lastNucMin.py", pseudoMol, prevResNum, curResNum, nextResNum, curRot, nextRot, prevResType, curResType, nextResType, curResAtoms, nextResAlreadyBuilt, fixPrevPhosOxy)
        #    raise Exception
        
        synRefinementScore = __runRefinement(molNum, chain, prevPrevResNum or prevResNum or curResNum, nextResNum, ignoreDensity)
        
        #if the syn refinement wasn't better than the anti refinement, then go back to the anti structure
        if synRefinementScore >= antiRefinementScore:
            pseudoMol.setCootNucs(prevResNum or curResNum, nextResNum, antiRefinementStruc)
            #we don't need to run pseudoMol.updateRes() here since it was run above
        else:
            newCoords = pseudoMol.updateRes(prevResNum or curResNum, nextResNum)
            score = synRefinementScore
            selectedStartingStruc = "syn"
    else:
        #update the Pseudomolecule object so that it knows about the moved atoms
        newCoords = pseudoMol.updateRes(prevResNum or curResNum, nextResNum)
    
    
    #try a high-anti structure if the best score is still too high
    if score > HIGH_ANTI_REBUILDING_CUTOFF:
        
        prevRefinementStruc = pseudoMol.getCootNucs(prevResNum or curResNum, nextResNum)
        
        pseudoMol.setCootNucs(prevResNum or curResNum, nextResNum, preMinimizationStruc, False)
        highAntiSugarCoords = rotateSugar(initSugarCoords, curResAtoms, "high-anti")
        pseudoMol.replaceSugar(curResNum, highAntiSugarCoords)
        #raise Exception
        
        highAntiRefinementScore = __runRefinement(molNum, chain, prevPrevResNum or prevResNum or curResNum, nextResNum, ignoreDensity)
        
        if highAntiRefinementScore >= score:
            pseudoMol.setCootNucs(prevResNum or curResNum, nextResNum, prevRefinementStruc)
        else:
            newCoords = pseudoMol.updateRes(prevResNum or curResNum, nextResNum)
            score = highAntiRefinementScore
            selectedStartingStruc = "high-anti"
    
    #if desired, print a summary table listing the refinement scores and the selected structure
    if PRINT_SUMMARY_TABLE:
        print "********************************"
        print "*  Starting struc  *   Score   *"
        print "********************************"
        print "* Anti             *  %s  *" % __stringifyScore(antiRefinementScore)
        print "* Syn              *  %s  *" % __stringifyScore(synRefinementScore)
        print "* High-anti        *  %s  *" % __stringifyScore(highAntiRefinementScore)
        print "********************************"
        print "* Using %s minimization" % selectedStartingStruc + " " * (9-len(selectedStartingStruc)) + " *"
        print "********************************"
        #from time import sleep; sleep(2)
    
    
    #pprint(newCoords)
    
    #if we haven't already, build the non-bridging oxygens
    if nextResAtoms is None:
        curResIndex = pseudoMol.resIndex(curResNum)
        #print "curResIndex =", curResIndex
        #print "curResNum =", curResNum
        
        #if curResIndex == 0:
        if prevResNum is None:
            #if this is the first (but not only) nucleotide
            (newCurResAtoms, newNextResAtoms) = newCoords
            phosOxyCoords = buildInitOrTerminalPhosOxy(newCurResAtoms)
            if phosOxyCoords is not None:
                pseudoMol.addPhosOxy(curResNum, phosOxyCoords)
            
        elif pseudoMol.numAtomsFromIndex(curResIndex+1) == 1 and not pseudoMol.connectedToNextFromIndex(curResIndex+1):
            #if the next residue is just a terminal 3' phosphate, then we need to add non-bridging oxygens
            (newPrevResAtoms, newCurResAtoms, newNextResAtoms) = newCoords
            phosOxyCoords5 = buildPhosOxy(newCurResAtoms, newPrevResAtoms)
            if phosOxyCoords5 is not None:
                pseudoMol.addPhosOxy(curResNum, phosOxyCoords5)
            phosOxyCoords3 = buildInitOrTerminalPhosOxy(newNextResAtoms, newCurResAtoms)
            if phosOxyCoords3 is not None:
                pseudoMol.addPhosOxy(nextResNum, phosOxyCoords3)
            
        else:
            #if this is a middle nucleotide of a chain
            #from pprint import pprint; pprint(newCoords)
            (newPrevResAtoms, newCurResAtoms, newNextResAtoms) = newCoords
            phosOxyCoords = buildPhosOxy(newCurResAtoms, newPrevResAtoms)
            if phosOxyCoords is not None:
                pseudoMol.addPhosOxy(curResNum, phosOxyCoords)
    
    #clear the fixed atoms and torsional restraints
    __clearResRestraints(molNum)
    
    return (newCoords, score)


def __stringifyScore(score):
    """Format a refinement score for printing in the summary table
    ARGUMENTS:
        score - the refinement
    RETURNS:
        score formatted as a fixed length string or "Not run" if score is None
    """
    
    if score == REFINEMENT_FAIL_SCORE:
        return " Failed"
    elif score is not None:
        return "%7.3f" % score
    else:
        return "Not run"


def __fixAtoms(pseudoMol, prevResNumFull, curResNumFull, nextResNumFull, fixPrevPhosOxy):
    """Fix and restrain the appropriate atoms for minimization
    
    ARGUMENTS:
        pseudoMol           - the PseudoMolecule object representing the molecule to fix atoms in
        prevResNumFull      - the residue number of the previous residue using Coot numbering (i.e. starting at 1, not 0) including insertion code
        curResNumFull       - the residue number using Coot numbering (i.e. starting at 1, not 0) including insertion code
        nextResNumFull      - the residue number of the next residue using Coot numbering (i.e. starting at 1, not 0) including insertion code
    OPTIONAL ARGUMENTS:
        fixPrevPhosOxy      - whether to fix the phosphoryl oxygens of the previous nucleotide during minimization
                              Defaults to False
    RETURNS:
        None
    """
    
    molNum = pseudoMol.molNum()
    chain = pseudoMol.chain
    
    (prevResNum, prevResInsCode) = __splitResNum(prevResNumFull)
    (curResNum,  curResInsCode)  = __splitResNum(curResNumFull)
    (nextResNum, nextResInsCode) = __splitResNum(nextResNumFull)
    
    #debugOut = open("fixAtoms.txt", "a")
    #debugOut.write("Minimizing "+ ", ".join(map(str, [prevResNumFull, curResNumFull, nextResNumFull])) + "\n")
    
    fixList = [] #the list of atoms to fix
    
    #fix atoms from the previous residue
    if prevResNum is not None:
        
        #fix the previous non-bridging oxygens if this is the first minimization of a recalcSuites
        if fixPrevPhosOxy:
            prevResMobileAtoms = PREV_RES_MOBILE_ATOMS
        else:
            prevResMobileAtoms = PREV_RES_MOBILE_ATOMS.union(["OP1", "OP2"])
        
        
        for curAtom in  pseudoMol.getAtomNames(prevResNum, strip=False):
            if curAtom.strip() not in prevResMobileAtoms:
                fixList.append([chain, prevResNum, prevResInsCode, curAtom, ""])
                #debugOut.write(",".join(map(str, [chain, prevResNum, prevResInsCode, curAtom, ""])) + "\n")
    
    
    #restrain the starting and ending phosphates
    add_extra_start_pos_restraint(molNum, chain, curResNum , curResInsCode, " P  ", "", HARMONIC_STD_DEV_PHOSPHATE)
    add_extra_start_pos_restraint(molNum, chain, nextResNum, nextResInsCode, " P  ", "", HARMONIC_STD_DEV_PHOSPHATE)
    #debugOut.write(",".join(map(str, [molNum, chain, curResNum , curResInsCode, " P  ", "", HARMONIC_STD_DEV_PHOSPHATE])) + "\n")
    #debugOut.write(",".join(map(str, [molNum, chain, nextResNum, nextResInsCode, " P  ", "", HARMONIC_STD_DEV_PHOSPHATE])) + "\n")
    
    #restrain atoms from the current residue base + C1'
    for curAtom in  pseudoMol.getAtomNames(curResNum, strip=False):
        if curAtom.strip() not in CUR_RES_MOBILE_ATOMS:
            #fixList.append([chain, resNum, "", curAtom, ""])
            add_extra_start_pos_restraint(molNum, chain, curResNum  , curResInsCode, curAtom, "", HARMONIC_STD_DEV_BASE)
            #debugOut.write(",".join(map(str, [molNum, chain, curResNum  , curResInsCode, curAtom, "", HARMONIC_STD_DEV_BASE])) + "\n")
    
    #fix atoms from the next residue
    if nextResNum is not None:
        for curAtom in  pseudoMol.getAtomNames(nextResNum, strip=False):
            if curAtom.strip() not in NEXT_RES_MOBILE_ATOMS:
                fixList.append([chain, nextResNum, nextResInsCode, curAtom, ""])
                #debugOut.write(",".join(map(str, [chain, nextResNum, nextResInsCode, curAtom, ""])) + "\n")
    
    mark_multiple_atoms_as_fixed(molNum, fixList, 1)
    
    #debugOut.close()


def __setTorsionRestraints(molNum, chain, prevResNumFull, curResNumFull, nextResNumFull, curRot, nextRot, curResType, nextResAlreadyBuilt = False):
    """Set the appropriate torsional restraints for the minimization
    
    ARGUMENTS:
        molNum              - the Coot molecule number
        chain               - the chain name
        prevResNumFull      - the residue number of the previous residue including insertion codes
        curResNumFull       - the residue number to build including insertion codes
        nextResNumFull      - the residue number of the next residue including insertion codes
        curRot              - the rotamer to use when building the current residue
        nextRot             - the rotamer to use when building the next residue
        curResType          - the residue type (i.e. A, G, C, or U) of the current residue
    OPTIONAL ARGUMENTS:
        nextResAlreadyBuilt - should be True if the next residue is already built (i.e. if we are being
                              called from recalcCoords instead of calcCoords)
                              Defaults to False.
    RETURNS:
        None
    """
    
    (prevResNum, prevResInsCode) = __splitResNum(prevResNumFull)
    (curResNum,  curResInsCode)  = __splitResNum(curResNumFull)
    (nextResNum, nextResInsCode) = __splitResNum(nextResNumFull)
    
    if curRot is not None:
        #print "Adding custom torsional restraints for residue", curResNum, "currot %s" % curRot
        
        if isinstance(curRot, str):
            
            prevDelta = rotData.prevDeltaMean(curRot)
            ep        = rotData.epMean   (curRot)
            zeta      = rotData.zetaMean (curRot)
            alpha     = rotData.alphaMean(curRot)
            beta      = rotData.betaMean (curRot)
            gamma     = rotData.gammaMean(curRot)
            
            prevDeltaSD = rotData.prevDeltaSD(curRot) * TORSION_STD_DEV_MOD
            epSD        = rotData.epSD(curRot)        * TORSION_STD_DEV_MOD
            zetaSD      = rotData.zetaSD (curRot)     * TORSION_STD_DEV_MOD
            alphaSD     = rotData.alphaSD(curRot)     * TORSION_STD_DEV_MOD
            betaSD      = rotData.betaSD (curRot)     * TORSION_STD_DEV_MOD
            gammaSD     = rotData.gammaSD(curRot)     * TORSION_STD_DEV_MOD
            
        else:
            prevDelta = curRot[0]
            ep        = curRot[1]
            zeta      = curRot[2]
            alpha     = curRot[3]
            beta      = curRot[4]
            gamma     = curRot[5]
            
            prevDeltaSD = MANUAL_TORSION_STD_DEV
            epSD        = MANUAL_TORSION_STD_DEV
            zetaSD      = MANUAL_TORSION_STD_DEV
            alphaSD     = MANUAL_TORSION_STD_DEV
            betaSD      = MANUAL_TORSION_STD_DEV
            gammaSD     = MANUAL_TORSION_STD_DEV
        
        #previous delta
        if prevDelta is not None:
            add_extra_torsion_restraint(molNum, chain, prevResNum, prevResInsCode, " C5'", "",
                                                chain, prevResNum, prevResInsCode, " C4'", "",
                                                chain, prevResNum, prevResInsCode, " C3'", "",
                                                chain, prevResNum, prevResInsCode, " O3'", "",
                                                prevDelta, prevDeltaSD, 1)
        
        #epsilon
        if ep is not None:
            add_extra_torsion_restraint(molNum, chain, prevResNum, prevResInsCode, " C4'", "",
                                                chain, prevResNum, prevResInsCode, " C3'", "",
                                                chain, prevResNum, prevResInsCode, " O3'", "",
                                                chain, curResNum,  curResInsCode,  " P  ", "",
                                                ep, epSD, 1)
        
        #zeta
        if zeta is not None:
            add_extra_torsion_restraint(molNum, chain, prevResNum, prevResInsCode, " C3'", "",
                                                chain, prevResNum, prevResInsCode, " O3'", "",
                                                chain, curResNum,  prevResInsCode, " P  ", "",
                                                chain, curResNum,  curResInsCode,  " O5'", "",
                                                zeta, zetaSD, 1)
        
        #alpha
        if alpha is not None:
            add_extra_torsion_restraint(molNum, chain, curResNum,  curResInsCode,  " C5'", "",
                                                chain, curResNum,  curResInsCode,  " O5'", "",
                                                chain, curResNum,  curResInsCode,  " P  ", "",
                                                chain, prevResNum, prevResInsCode, " O3'", "",
                                                alpha, alphaSD, 1)
        
        #beta
        if beta is not None:
            add_extra_torsion_restraint(molNum, chain, curResNum,  curResInsCode, " P  ", "",
                                                chain, curResNum,  curResInsCode, " O5'", "",
                                                chain, curResNum,  curResInsCode, " C5'", "",
                                                chain, curResNum,  curResInsCode, " C4'", "",
                                                beta, betaSD, 1)
        
        #gamma
        if gamma is not None:
            add_extra_torsion_restraint(molNum, chain, curResNum  , curResInsCode, " O5'", "",
                                                chain, curResNum  , curResInsCode, " C5'", "",
                                                chain, curResNum  , curResInsCode, " C4'", "",
                                                chain, curResNum  , curResInsCode, " C3'", "",
                                                gamma, gammaSD, 1)
        
    
    #add constraints on delta and the ring torsions
    curPucker = None
    delta = None
    if isinstance(curRot, str):
        curPucker = puckerList[curRot][1]
        delta   = rotData.curDeltaMean(curRot)
        deltaSD = rotData.curDeltaSD(curRot)*TORSION_STD_DEV_MOD
        
    elif isinstance(nextRot, str):
        curPucker = puckerList[nextRot][0]
        delta   = rotData.prevDeltaMean(nextRot)
        deltaSD = rotData.prevDeltaSD(nextRot)*TORSION_STD_DEV_MOD
    
    if delta is not None:
        add_extra_torsion_restraint(molNum, chain, curResNum  , curResInsCode, " C5'", "",
                                            chain, curResNum  , curResInsCode, " C4'", "",
                                            chain, curResNum  , curResInsCode, " C3'", "",
                                            chain, curResNum  , curResInsCode, " O3'", "",
                                            delta, deltaSD, 1)        
        
        
    if curPucker is not None:
        #print "adding ring constraints"
        
        if curPucker == 2:
            curNu0 = NU0_C2
            curNu1 = NU1_C2
            curNu4 = NU4_C2
        else:
            curNu0 = NU0_C3
            curNu1 = NU1_C3
            curNu4 = NU4_C3
        
        #current nu0
        add_extra_torsion_restraint(molNum, chain, curResNum  , curResInsCode, " C4'", "",
                                            chain, curResNum  , curResInsCode, " O4'", "",
                                            chain, curResNum  , curResInsCode, " C1'", "",
                                            chain, curResNum  , curResInsCode, " C2'", "",
                                            curNu0, NU_STD_DEV, 1)
        
        #current nu1
        add_extra_torsion_restraint(molNum, chain, curResNum  , curResInsCode, " O4'", "",
                                            chain, curResNum  , curResInsCode, " C1'", "",
                                            chain, curResNum  , curResInsCode, " C2'", "",
                                            chain, curResNum  , curResInsCode, " C3'", "",
                                            curNu1, NU_STD_DEV, 1)
        
        #current nu4
        add_extra_torsion_restraint(molNum, chain, curResNum  , curResInsCode, " C5'", "",
                                            chain, curResNum  , curResInsCode, " C4'", "",
                                            chain, curResNum  , curResInsCode, " O4'", "",
                                            chain, curResNum  , curResInsCode, " C1'", "",
                                            curNu4, NU_STD_DEV, 1)
        
        
    #print "adding chi restraint"
    #redifine the chi restraints, since Coot uses a very strange default value
    #(Note that even though the base is fixed, this chi restraint can affect the sugar positioning)
    if curResType == "G" or curResType == "A":
        baseAtom1 = " N9 "
        baseAtom2 = " C4 "
    else:
        baseAtom1 = " N1 "
        baseAtom2 = " C2 "
    add_extra_torsion_restraint(molNum, chain, curResNum  , curResInsCode, " O4'",    "",
                                        chain, curResNum  , curResInsCode, " C1'",    "",
                                        chain, curResNum  , curResInsCode, baseAtom1, "",
                                        chain, curResNum  , curResInsCode, baseAtom2, "",
                                        CHI_MEAN, CHI_STD_DEV, 2)
    #Note that this torsion has a period of 2 (which means a period of 360/2=180 degrees) to account for anti and syn rotations
        
    if nextRot is not None:
        
        if isinstance(nextRot, str):
            ep    = rotData.epMean(nextRot)
            epSD  = rotData.epSD(nextRot) * TORSION_STD_DEV_MOD
            
            if nextResAlreadyBuilt:
                zeta  = rotData.zetaMean(nextRot)
                alpha = rotData.alphaMean(nextRot)
                beta  = rotData.betaMean(nextRot)
                gamma = rotData.gammaMean(nextRot)
                
                zetaSD  = rotData.zetaSD(nextRot)  * TORSION_STD_DEV_MOD
                alphaSD = rotData.alphaSD(nextRot) * TORSION_STD_DEV_MOD
                betaSD  = rotData.betaSD(nextRot)  * TORSION_STD_DEV_MOD
                gammaSD = rotData.gammaSD(nextRot) * TORSION_STD_DEV_MOD
            
        else:
            ep = nextRot[1]
            epSD = MANUAL_TORSION_STD_DEV
            
            if nextResAlreadyBuilt:
                zeta  = nextRot[2]
                alpha = nextRot[3]
                beta  = nextRot[4]
                gamma = nextRot[5]
                
                zetaSD  = MANUAL_TORSION_STD_DEV
                alphaSD = MANUAL_TORSION_STD_DEV
                betaSD  = MANUAL_TORSION_STD_DEV
                gammaSD = MANUAL_TORSION_STD_DEV
        
        
        #print "Adding custom torsional restraints for nextrot %s" % nextRot
        if ep is not None:
            add_extra_torsion_restraint(molNum, chain, curResNum  , curResInsCode,  " C4'", "",
                                                chain, curResNum  , curResInsCode,  " C3'", "",
                                                chain, curResNum  , curResInsCode,  " O3'", "",
                                                chain, nextResNum , nextResInsCode, " P  ", "",
                                                ep, epSD, 1)
        
        if nextResAlreadyBuilt:
            #if we're rebuilding a section and about to meet up with the existing structure, we need more constraints
            #zeta
            if zeta is not None:
                #torsions can be None when we're rotamerizing and the last nucleotide is missing atoms
                add_extra_torsion_restraint(molNum, chain, curResNum , curResInsCode,  " C3'", "",
                                                    chain, curResNum , curResInsCode,  " O3'", "",
                                                    chain, nextResNum, nextResInsCode, " P  ", "",
                                                    chain, nextResNum, nextResInsCode, " O5'", "",
                                                    zeta, zetaSD, 1)
            
            #alpha
            if alpha is not None:
                add_extra_torsion_restraint(molNum, chain, nextResNum, nextResInsCode, " C5'", "",
                                                    chain, nextResNum, nextResInsCode, " O5'", "",
                                                    chain, nextResNum, nextResInsCode, " P  ", "",
                                                    chain, curResNum , curResInsCode,  " O3'", "",
                                                    alpha, alphaSD, 1)
            
            #beta
            if beta is not None:
                add_extra_torsion_restraint(molNum, chain, nextResNum, nextResInsCode, " P  ", "",
                                                    chain, nextResNum, nextResInsCode, " O5'", "",
                                                    chain, nextResNum, nextResInsCode, " C5'", "",
                                                    chain, nextResNum, nextResInsCode, " C4'", "",
                                                    beta, betaSD, 1)
            
            #gamma
            if gamma is not None:
                add_extra_torsion_restraint(molNum, chain, nextResNum, nextResInsCode, " O5'", "",
                                                    chain, nextResNum, nextResInsCode, " C5'", "",
                                                    chain, nextResNum, nextResInsCode, " C4'", "",
                                                    chain, nextResNum, nextResInsCode, " C3'", "",
                                                    gamma, gammaSD, 1)


def __runRefinement(molNum, chain, startingResNumFull, endingResNumFull, ignoreDensity = False):
    """Run the minimization using Coot's built-in refinement functions
    
    ARGUMENTS:
        molNum              - the Coot molecule number
        chain               - the chain name
        startingResNumFull  - the residue number (including insertion code) of the first nucleotide to refine
        endingResNumFull    - the residue number (including insertion code) of the last nucleotide to refine
    OPTIONAL ARGUMENTS:
        ignoreDensity - ignore the density when performing the minimization
                        defaults to False
    RETURNS:
        the refinement score
    """
    
    print "*************************"
    print "About to refine residues", startingResNumFull, "-", endingResNumFull
    print "*************************"
    
    (startingResNum, startingResInsCode) = __splitResNum(startingResNumFull)
    (endingResNum,   endingResInsCode)   = __splitResNum(endingResNumFull)
    #refine_zone and refine_zone_with_score don't support insertion codes, so for now the insertion codes gets ignored
    #this will cause problems if trying to rotamerize a chain that uses insertion codes
    
    if ignoreDensity:
        refinementResults = regularize_zone_with_score(molNum, chain, startingResNum, endingResNum, "")
    else:
        #note that we need to use refine_zone_with_score, not refine_residues, since refine_residues assumes that things aren't bonded if they're more than three Angstroms apart
        refinementResults = refine_zone_with_score(molNum, chain, startingResNum, endingResNum, "")
    
    accept_regularizement()
    
    return __calcRefinementScore(refinementResults)


def __clearResRestraints(molNum):
    """Clear the torsional restraints and fixed atoms that applied to the residues we were just minimizing
    
    ARGUMENTS:
        molNum  - the Coot molecule number
    RETURNS:
        None
    """
    clear_all_fixed_atoms(molNum)
    delete_all_extra_restraints(molNum)


#a regex to separate a residue name and insertion code
__splitResNumRE = re.compile("(-?\d+)([A-Za-z]?)$")

def __splitResNum(resNumFull):
    """Split a residue number into number and insertion code
    
    ARGUMENTS:
        resNumFull - a residue number potentially contianing an insertion code
    RETURNS:
        resNum     - the residue number without the insertion code (as an int)
        insCode    - the insertion code
    """
    
    if resNumFull is None:
        return (None, None)
    else:
        (resNum, insCode) = __splitResNumRE.match(str(resNumFull)).groups()
        resNum = int(resNum) #the Coot SWIG functions require an integer argument, so we do explicit cast here
        return (resNum, insCode)


def __changeMonomerRestraints():
    """Change the monomer restraints for RNA residues
    
    ARGUMENTS:
        None
    RETURNS:
        origRestraints - the original monomer restraints (so we can restore them later)
    """
    
    origRestraints = {}
    
    for curNuc in ("A", "G", "C", "U"):
        origRestraints[curNuc] = monomer_restraints(curNuc)
        newRestraints = monomer_restraints(curNuc) #we need a new copy of the dictionary so we don't modify the originalRestraints dict
        
        #this next line is redundant with use_only_extra_torsion_restraints_for_torsions(1), but it certainly won't hurt anything
        #use_only_extra_torsion_restraints_for_torsions(1) also has the advantage of turning off the link torsion restraints
        #(i.e. restraints that span two nucleotides), but the link torsion restraints don't seem to be used during minimization
        #regardless
        newRestraints['_chem_comp_tor'] = []
            #we don't need any of the existing restraints since we're rewriting all of the named restraints
            #and the CONST restraints are redundant with the planar restraints (and seem to be ignored by Coot anyway)
        
        #make the plane restraints tighter (change standard deviation from 0.020 to 0.001)
        #since we really don't want non-planar bases (they can be a pain to fix later)
        for curPlane in newRestraints['_chem_comp_plane_atom']:
            curPlane[2] = 0.001
        
        #if we're going to be using Phenix's restraints, then remove the default backbone bond and angle restraints
        if USE_PHENIX_RESTRAINTS == True:
            newRestraints['_chem_comp_angle'] = [curAngle for curAngle in newRestraints['_chem_comp_angle'] if not(curAngle[0].strip() in PHENIX_NEW_RESTRAINT_ATOMS and curAngle[1].strip() in PHENIX_NEW_RESTRAINT_ATOMS and curAngle[2].strip() in PHENIX_NEW_RESTRAINT_ATOMS)]
            
            #add the Phenix restraints that aren't pucker specific
            newRestraints['_chem_comp_angle'].extend([[" C2'", " C3'", " C4'", 102.6, 1.0],
                                                      [" C5'", " O5'", " P  ", 120.9, 1.0],
                                                      [" O5'", " P  ", " OP1", 110.7, 1.000],
                                                      [" OP1", " P  ", " OP2", 119.6, 1.000],
                                                      [" O5'", " P  ", " OP2", 110.7, 1.000]])
            
            
            newRestraints['_chem_comp_bond'] = [curBond for curBond in newRestraints['_chem_comp_bond'] if not (curBond[0].strip() in PHENIX_NEW_RESTRAINT_ATOMS and curBond[1].strip() in PHENIX_NEW_RESTRAINT_ATOMS)]
            
            #add in the bond restraints that aren't pucker specific
            newRestraints['_chem_comp_bond'].extend([[" P  ", " O5'", 'single', 1.593, 0.015],
                                                     [" P  ", " OP1", 'deloc',  1.485, 0.020],
                                                     [" P  ", " OP2", 'deloc',  1.485, 0.020]])
            
            #remove the polymer type so that we can replace the link restraints (restraints that span two nucleotides)
            #there's no way to set proper link restraints, so all link restraints are set as extra restraints in phenixRestraints.py
            newRestraints['_chem_comp'][3] = "" #previously, this value was "RNA" 
        
        set_monomer_restraints(curNuc, newRestraints)
    
    return origRestraints


def __changeCootSettings():
    """Change Coot's minimization settings
    
    ARUGMENTS:
        None
    RETURNS:
        these values are returned so the minimization settings can be restored later
        origImmediateReplaceValue    - the original refinement_immediate_replacement_state() value
        origRefineWithTorsionsValue  - the original refine_with_torsion_restraints_state()
        origWeightValue              - the original matrix_state() value
        origMonomerRestraints        - the original monomer restraints (i.e. torsion restraints)
        origUseOnlyExtraTorsionValue - the original settings for use_only_extra_torsion_restraints_for_torsions
    """
    
    #we're going to have to change some Coot settings during the coordinate calculation
    #so check what the existing values are so we can set them back
    origImmediateReplaceValue    = refinement_immediate_replacement_state()         #don't ask the user to accept the refinement
    origRefineWithTorsionsValue  = refine_with_torsion_restraints_state()           #whether or not we use torsional contraints
    origWeightValue              = matrix_state()                                   #how much weight is placed on geometry vs. map constraints
    #origNumStepsPerFrameValue   = refinement_refine_per_frame_state()              #how many steps of refinement occurr in between updating graphics
    origUseOnlyExtraTorsionValue = use_only_extra_torsion_restraints_for_torsions_state() #should we ignore the standard torsion restraints and only use the user-defined ones
    
    #set the values that we want
    set_refinement_immediate_replacement(1)
    set_refine_with_torsion_restraints(1)
    set_matrix(REFINE_MAP_WEIGHT)
    #dragged_refinement_steps_per_frame()
        #in theory, a large value here should cause Coot to not update the graphics as it refines (which should speed up refinement)
        #but Coot doesn't seem to like it when this value is too high
        #and Windows thinks that Coot has frozen if it takes too long
        #so we just leave this at the default
    #set_use_only_extra_torsion_restraints_for_torsions(1)
    set_use_only_extra_torsion_restraints_for_torsions(0)
        #this function seems to entirely disable torsion restraints, which is the opposite of what it's supposed to do
        #as a result, we make sure it's turned off rather than making sure it's turned on
        #I should probably fix this function/flag in Coot
        
    #remove all built-in torsional restraints from RNA nucleotides
    origMonomerRestraints = __changeMonomerRestraints()
    
    return (origImmediateReplaceValue, origRefineWithTorsionsValue, origWeightValue, origMonomerRestraints, origUseOnlyExtraTorsionValue)


def __restoreCootSettings(origCootSettings):
    """Restore Coot's minimization settings
    
    ARGUMENTS:
        origCootSettings - Coot settings to be restored (as returned by __changeCootSettings)
    RETURNS:
        None
    """
    
    (origImmediateReplaceValue, origRefineWithTorsionsValue, origWeightValue, origMonomerRestraints, origUseOnlyExtraTorsionValue) = origCootSettings
    
    set_refinement_immediate_replacement(origImmediateReplaceValue)
    set_refine_with_torsion_restraints(origRefineWithTorsionsValue)
    set_matrix(origWeightValue)
    set_use_only_extra_torsion_restraints_for_torsions(origUseOnlyExtraTorsionValue)
    for curNuc in ("A", "G", "C", "U"):
        set_monomer_restraints(curNuc, origMonomerRestraints[curNuc])


def __calcRefinementScore(refinementResults):
    """Calculate an overall score for the quality of the refinement results
    
    ARGUMENTS:
        refinementResults - the output of refine_residues()
    RETURNS:
        overallScore      - an overall refinement score based on the bond length, angles, and torsions
    """
    
    #from pprint import pprint; pprint(refinementResults)
    
    scoreList = refinementResults[2]
    
    if not scoreList:
        #if there is no scorelist, then Coot gave up on the minimization
        #(Modifying Coot to return a score when the minimization times out is do-able, but it would slightly slow
        #down all minimizations unless more significant changes to the minimizer code are made.  Typically, if we've
        #timed out, then something has gone wrong anyway, so simply using REFINEMENT_FAIL_SCORE is a workable solution.)
        
        return REFINEMENT_FAIL_SCORE 
    
    overallScore = 0
    for curScore in scoreList:
        if curScore[0] in ["Bonds", "Angles", "Torsions"]:
            #TODO: add in start_pos scores?
            overallScore += curScore[2]

    return overallScore
    
    
def __addPhenixRestraints(molNum, chain, prevPrevResNumFull, prevResNumFull, curResNumFull, nextResNumFull, curRot, nextRot, curResType, nextResAlreadyBuilt = False, onlyPucker = None):
    """Add angle and bond restraints for the specified nucleotides using the Phenix pucker-specific restraint values
    
    ARGUMENTS
        molNum              - the Coot molecule number
        chain               - the chain name
        prevPrevResNumFull  - the residue number of residue i-2 (needed to restrain the non-bridging oxygens of residue prevResNum)
        prevResNumFull      - the residue number of the previous residue using Coot numbering (i.e. starting at 1, not 0)
        curResNumFull       - the residue number to build using Coot numbering (i.e. starting at 1, not 0)
        nextResNumFull      - the residue number of the next residue using Coot numbering (i.e. starting at 1, not 0)
        curRot              - the rotamer to use when building the current residue
        nextRot             - the rotamer to use when building the next residue
        curResType          - the residue type (i.e. A, G, C, or U) of the current residue
        OPTIONAL ARGUMENTS:
        nextResAlreadyBuilt - should be True if the next residue is already built (i.e. if we are being
                              called from recalcCoords instead of calcCoords)
                              Defaults to False.
        onlyPucker          - should only be provided if curRot and nextRot are None.  This is used when we're only building
                              a single nucleotide, so we don't have a full suite and therefore can't predict a conformer
    RETURNS:
        None
    """
    
    if curResType == "A" or curResType == "G":
        glycosidicBaseAtom = " N9 "
    else:
        glycosidicBaseAtom = " N1 "
    
    #figure out the sugar puckers
    prevPucker = None
    curPucker = None
    nextPucker = None
    
    
    if isinstance(curRot, basestring): #basestring includes unicode and non-unicode strings.  It's overkill here, but it can't hurt:
        (prevPucker, curPucker) = puckerList[curRot]
    elif isinstance(curRot, list):
        #if curRot is a list, then it's not really a rotamer, it's just a list of the current torsion values
        #in this case, we can use delta to figure out the appropriate pucker
        
        if curRot[0] < 114.5:
            prevPucker = 3
        else:
            prevPucker = 2
        
        if not isinstance(nextRot, basestring):
            #if nextRot is a real rotamer, then don't bother with delta-based guesses
            
            if curRot[6] < 114.5:
                curPucker = 3
            else:
                curPucker = 2
    
    
    if isinstance(nextRot, basestring): #basestring includes unicode and non-unicode strings.  It's overkill here, but it can't hurt
        #if nextRot is a rotamer
        (curPucker, nextPucker) = puckerList[nextRot]
    elif isinstance(nextRot, list):
        #if nextRot is a list, then it's not really a rotamer, it's just a list of the current torsion values
        #in this case, we can use delta to figure out the appropriate pucker
        
        if curPucker is None:
            #if we had a real rotamer for curRot, then don't overwrite it with our delta-based guess
            if nextRot[0] < 114.5:
                curPucker = 3
            else:
                curPucker = 2
        
        if nextRot[6] < 114.5:
            nextPucker = 3
        else:
            nextPucker = 2
    
    #curPucker may get set twice here, but that doesn't matter, since curRot and nextRot *must* agree on their overlapping pucker
    
    #if we don't have any conformers, then set curPucker to onlyPucker
    if onlyPucker is not None and curPucker is None:
        curPucker = onlyPucker
    
    #split the residue numbers
    (prevPrevResNum,  prevPrevResInsCode) = __splitResNum(prevPrevResNumFull)
    (prevResNum,      prevResInsCode)     = __splitResNum(prevResNumFull)
    (curResNum,       curResInsCode)      = __splitResNum(curResNumFull)
    (nextResNum,      nextResInsCode)     = __splitResNum(nextResNumFull)
    
    #figure out the number of the residue before prevResNum
    phenixRestraints.setAngleRestraints(molNum, chain, prevPrevResNum, prevResNum, curResNum, nextResNum, glycosidicBaseAtom, prevPucker, curPucker, nextPucker, nextResAlreadyBuilt)
    phenixRestraints.setBondRestraints(molNum, chain, prevResNum, curResNum, nextResNum, glycosidicBaseAtom, prevPucker, curPucker, nextPucker, nextResAlreadyBuilt)


def enablePhenixRestraints():
    """Enable Phenix's pucker-specific restraints for minimization
    
    ARGUMENTS:
        None
    RETURNS:
        True if enabling the restraints succceeded, False otherwise
    NOTE:
        Phenix restraints are only available in versions of Coot newer than 3926
    """
    
    global USE_PHENIX_RESTRAINTS
    
    if svn_revision() >= 3926:
        #Coot revision 3926 has a bug fix to fix non-bonded constraints when using extra bond and angle restraints
        USE_PHENIX_RESTRAINTS = True
        print "Using Phenix restraints"
        return True
    else:
        print "Coot must be newer than 0.7-pre r3926 to use Phenix restraints.  Using Coot/CCP4 restraints."
        return False


def disablePhenixRestraints():
    """Disable Phenix's pucker-specific restraints for minimization and use the Coot/CCP4 restraints instead
    
    ARGUMENTS:
        None
    RETURNS:
        True (for consistancy with enablePhenixRestraints())
    """
    
    global USE_PHENIX_RESTRAINTS
    USE_PHENIX_RESTRAINTS = False
    print "Using Coot/CCP4 restraints"
    return True


def usingPhenixRestraints():
    """Check if we are using Phenix's pucker-specific restraints for minimization
    
    ARGUMENTS:
        None
    RETURNS:
        True if we are using Phenix restraints
        False if we are using Coot/CCP4 restraints
    """
    
    return USE_PHENIX_RESTRAINTS


def __fixEntireResidue(molNum, chain, resNumFull, atomList):
    """Fix all atoms in a residue
    
    ARGUMENTS:
        molNum      - the molecule number of the residue to fix
        chain       - the chain of the residue to fix
        resNumFull  - the residue number of the residue to fix including insertion code
        atomList    - a list of all atoms in the specified residue
    RETURNS:
        None
    """
    
    (resNum,  insCode) = __splitResNum(resNumFull)
    fixList = [[chain, resNum, "", curAtom, ""] for curAtom in atomList]
    mark_multiple_atoms_as_fixed(molNum, fixList, 1)


def buildOnlyPhosOxy(pseudoMol, resIndex, direction = 3):
    """Build the non-bridging oxygens onto a terminal nucleotide
    
    ARGUMENTS:
        resIndex - the index of the residue to build the phosphates onto
    OPTIONAL ARGUMENTS:
        direction - which side of the residue to add the phosphoryl oxygens to
    RETURNS:
        NONE
    """
    #print "In buildOnlyPhosOxy with resIndex =", resIndex
    
    if direction == 3:
        chain = pseudoMol.createPartialChainObjectFromIndex(resIndex-1, resIndex)
        
        prevResAtoms = chain.nucleotides[0].atoms
        curResAtoms  = chain.nucleotides[1].atoms
        
    else:
        chain = pseudoMol.createPartialChainObjectFromIndex(resIndex, resIndex)
        
        curResAtoms = chain.nucleotides[0].atoms
        prevResAtoms = None
    
    #from pprint import pprint
    #pprint(prevResAtoms)
    #pprint(curResAtoms)
    
    phosOxyCoords = buildInitOrTerminalPhosOxy(curResAtoms, prevResAtoms)
    pseudoMol.addPhosOxyFromIndex(resIndex, phosOxyCoords)


class ProgressDialogObject:
    """A class for displaying a progress dialog while calculating/minimizing backbone coordinates"""
    
    def __init__(self, window, totalNum):
        """Initialize a ProgressDialogObject object
        
        ARGUMENTS:
            window   - the window object to display the progress dialog in
                       if None, a new window will be created
            totalNum - the total number of steps for the progress bar
                       akin to the deprecated set_discrete_blocks functionality
        RETURNS:
            an initialized ProgressDialogObject object
        EFFECTS:
            If provided, the contents of window are replaced by a progress bar
            The contents of window can later be restored using restoreWindow()
        """
        
        self.__window = window
        self.__oldWindowSize = None
        self.__oldWindowContents = None
        self.__progressFracIncrement = None
        self.__curProgressFrac = 0
        self.__curProgressInt = 0
        self.__totalNum = totalNum
        self.__progressbar = None
        
        if self.__window is not None:
            windowChild = self.__window.get_child()
            if windowChild is not None:
                self.__oldWindowSize = window.get_size()
                self.__oldWindowContents = windowChild
                self.__window.remove(windowChild)
        else:
            self.__window = createRCraneWindowObject()
        
        self.__createProgressDialog()
    
    
    def __createProgressDialog(self):
        """Add the ProgressBar to self.__window
        
        ARGUMENTS:
            None
        RETURNS:
            None
        """
        
        windowBox = gtk.VBox(False, VBOX_SPACING)
        self.__window.add(windowBox)
        self.__window.resize(1,1) #remove any size constraints
        
        calculatingLabel = gtk.Label("\n  Calculating backbone coordinates...  \n\n  Please wait...  \n")
        calculatingLabel.set_justify(gtk.JUSTIFY_CENTER)
        windowBox.pack_start(calculatingLabel, False, False, HBOX_SPACING)
        
        self.__progressbar = gtk.ProgressBar()
        self.__progressbar.set_text("Built 0 of " + str(self.__totalNum) + " nucleotides")
        progressAlign = gtk.Alignment(xscale=0.9, xalign=0.5)
        progressAlign.add(self.__progressbar)
        windowBox.pack_start(progressAlign, False, False, HBOX_SPACING)
        
        self.__window.show_all()
        
        self.__progressFracIncrement = 1.0 / (self.__totalNum)
    
    
    def progress(self, step = 1):
        """Increment the progress bar
        
        OPTIONAL ARGUMENTS:
            step - how many steps to increment the progress bar by, defaults to 1
        RETURNS:
            None
        """
        
        self.__curProgressFrac += (step * self.__progressFracIncrement)
        self.__curProgressInt += step
        
        #make sure __curProgressFrac isn't > 1 due to a rounding error
        if self.__curProgressFrac > 1:
            self.__curProgressFrac = 1
        
        self.__progressbar.set_fraction(self.__curProgressFrac)
        self.__progressbar.set_text("Built " + str(self.__curProgressInt) + " of " + str(self.__totalNum) + " nucleotides")
    
    
    def restoreWindow(self):
        """Restore the original contents of the window and destroy the ProgressDialogObject object
        
        ARGUMENTS:
            None
        RETURNS:
            None
        EFFECTS:
            the original contents of the window are destroyed
            the ProgressDialogObject is set to None, since it is no longer useful
        """
        
        windowChild = self.__window.get_child()
        self.__window.remove(windowChild)
        
        if self.__oldWindowContents is not None:
            self.__window.add(self.__oldWindowContents)
        
        #set the window back to the size it was before
        if self.__oldWindowSize is not None:
            self.__window.resize(*self.__oldWindowSize)
        
        #set the object to none, since it's useless now
        #but don't actually call a destructor
        self = None
