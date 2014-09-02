#!/usr/bin/env python
"""Add minimization restraint using Phenix pucker-specific restraint values"""

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

from coot import add_extra_bond_restraint
#add_extra angle_restraint only exists in versions of Coot newer than 3901
#we can't use this module in version of Coot older than that
try:
    from coot import add_extra_angle_restraint
except ImportError:
    pass

PHENIX_ANGLE_SD = 1.0

#Note: If I ever want to use the non-bonded score that's output from the minimization, I'll need to add extra bond and angle restraints
#even for the fixed atoms.  Currently, restraints are only added if they can impact a mobile atom.  The restraints for fixed
#atoms won't affect the minimization itself, but a non-bonded constraint will still be calculated between fixed atoms, which
#will result in a significantly higher non-bonded score when using Phenix restraints as compared to regular Coot/CCP4 restraints
#Currently, __calcRefinementScore() in calcCoords.py only looks at bond, angle, and torsion restraints, so this increased
#score doesn't matter.


def setAngleRestraints(molNum, chain, prevPrevResNum, prevResNum, curResNum, nextResNum, glycosidicBaseAtom, prevPucker, curPucker, nextPucker, restrainNextRes = False):
    """Set angle restraints on the specified nucleotides to match the Phenix pucker-specific restraints
    
    ARGUMENTS:
        molNum              - the molecule number
        chain               - the chain
        prevPrevResNum      - the residue number of residue i-2
        prevResNum          - the residue number of the previous residue (using Coot numbering)
        curResNum           - the residue number of the current residue (using Coot numbering)
        nextResNum          - the residue number of the next residue (using Coot numbering)
        glycosidicBaseAtom  - the name of the base nitrogen atom (N1 or N9) for the current residue
        prevPucker          - the pucker of the previous residue (as an int: 2 or 3)
        curPucker           - the pucker of the current residue (as an int: 2 or 3)
        nextPucker          - the pucker of the previous residue (as an int: 2 or 3)
    OPTIONAL ARGUMENTS:
        restrainNextRes     - whether we should add restraints for atoms of the next residue (other than the phosphate, which is always restrained)
                              Defaults to False
    RETURNS:
        None
    NOTES:
        Note that restraints that are the same for C2'-endo and C3'-endo nucleotides are added in calcCoords.__setTorsionRestraints
    """
    
    restraintList = []
    linkRestraintList = []
    
    if prevResNum is not None:
        if prevPucker == 2:
            restraintList.extend([[prevResNum, " C2'", " C3'", " O3'", 109.5, 1.0],
                                  [prevResNum, " O3'", " C3'", " C4'", 109.4, 1.0]])
        elif prevPucker == 3:
            restraintList.extend([[prevResNum, " C2'", " C3'", " O3'", 113.7, 1.0],
                                  [prevResNum, " O3'", " C3'", " C4'", 113.0, 1.0]])
        
        #add the link restraints
        linkRestraintList.extend([[prevResNum, " O3'", curResNum,  " P  ", curResNum, " O5'", 103.579, 1.0],
                                  [prevResNum, " O3'", curResNum,  " P  ", curResNum, " OP1", 108.000, 1.0],
                                  [prevResNum, " O3'", curResNum,  " P  ", curResNum, " OP2", 108.000, 1.0],
                                  [prevResNum, " C3'", prevResNum, " O3'", curResNum, " P  ", 120.200, 1.0]])
        
        #restrain the previous phosphate
        if prevPrevResNum is not None:
            linkRestraintList.extend([[prevPrevResNum, " O3'", prevResNum,  " P  ", prevResNum, " OP1", 108.000, 1.0],
                                      [prevPrevResNum, " O3'", prevResNum,  " P  ", prevResNum, " OP2", 108.000, 1.0]])
            #print "Restraining previous non-bridging oxygens:", prevPrevResNum, "O3' - ",prevResNum, "P - ", prevResNum, "OP1/2"
            #from time import sleep; sleep(7)
        
    if curPucker == 2:
        restraintList.extend([[curResNum, glycosidicBaseAtom,  " C1'", " O4'", 108.2, 1.0],
                              [curResNum, glycosidicBaseAtom,  " C1'", " C2'", 114.0, 1.0],
                              [curResNum, " O4'",  " C1'", " C2'", 105.8, 1.0],
                              [curResNum, " C1'",  " C2'", " C3'", 101.5, 1.0],
                              [curResNum, " C2'",  " C3'", " O3'", 109.5, 1.0],
                              [curResNum, " O3'",  " C3'", " C4'", 109.4, 1.0],
                              [curResNum, " C3'",  " C4'", " C5'", 115.2, 1.0],
                              [curResNum, " C3'",  " C4'", " O4'", 106.1, 1.0],
                              [curResNum, " O4'",  " C4'", " C5'", 109.1, 1.0],
                              [curResNum, " C4'",  " O4'", " C1'", 109.7, 1.0],
                              [curResNum, " C4'",  " C5'", " O5'", 111.7, 1.0],
                              [curResNum, " C1'",  " C2'", " O2'", 111.8, 1.0],
                              [curResNum, " O2'",  " C2'", " C3'", 114.6, 1.0]])
    elif curPucker == 3:
        restraintList.extend([[curResNum, glycosidicBaseAtom,  " C1'", " O4'", 108.5, 1.0],
                              [curResNum, glycosidicBaseAtom,  " C1'", " C2'", 112.0, 1.0],
                              [curResNum, " O4'",  " C1'", " C2'", 107.6, 1.0],
                              [curResNum, " C1'",  " C2'", " C3'", 101.3, 1.0],
                              [curResNum, " C2'",  " C3'", " O3'", 113.7, 1.0],
                              [curResNum, " O3'",  " C3'", " C4'", 113.0, 1.0],
                              [curResNum, " C3'",  " C4'", " C5'", 116.0, 1.0],
                              [curResNum, " C3'",  " C4'", " O4'", 104.0, 1.0],
                              [curResNum, " O4'",  " C4'", " C5'", 109.8, 1.0],
                              [curResNum, " C4'",  " O4'", " C1'", 109.9, 1.0],
                              [curResNum, " C4'",  " C5'", " O5'", 111.5, 1.0],
                              [curResNum, " C1'",  " C2'", " O2'", 108.4, 1.0],
                              [curResNum, " O2'",  " C2'", " C3'", 110.7, 1.0]])
    
    if nextResNum is not None:
        
        #add the link restraints
        linkRestraintList.extend([[curResNum, " O3'", nextResNum, " P  ", nextResNum, " O5'", 103.579, 1.0],
                                  [curResNum, " O3'", nextResNum, " P  ", nextResNum, " OP1", 108.000, 1.0],
                                  [curResNum, " O3'", nextResNum, " P  ", nextResNum, " OP2", 108.000, 1.0],
                                  [curResNum, " C3'", curResNum,  " O3'", nextResNum, " P  ", 120.200, 1.0]])
        
        if restrainNextRes:
            if nextPucker == 2:
                restraintList.append([nextResNum, " C4'",  " C5'", " O5'", 111.7, 1.0])
            elif nextPucker == 3:
                restraintList.append([nextResNum, " C4'",  " C5'", " O5'", 111.5, 1.0])
            #all other restraints that would apply here are identical between C2' and C3' puckers, so they're inlucded as monomer restraints
        
    
    #apply the restraints
    for curRestraint in restraintList:
        (resNum, atom1, atom2, atom3, mean, sd) = curRestraint
        add_extra_angle_restraint(molNum, chain, resNum, "", atom1, "",
                                          chain, resNum, "", atom2, "",
                                          chain, resNum, "", atom3, "",
                                          mean, PHENIX_ANGLE_SD)
        
    for curRestraint in linkRestraintList:
        (resNum1, atom1, resNum2, atom2, resNum3, atom3, mean, sd) = curRestraint
        add_extra_angle_restraint(molNum, chain, resNum1, "", atom1, "",
                                          chain, resNum2, "", atom2, "",
                                          chain, resNum3, "", atom3, "",
                                          mean, PHENIX_ANGLE_SD)


def setBondRestraints(molNum, chain, prevResNum, curResNum, nextResNum, glycosidicBaseAtom, prevPucker, curPucker, nextPucker, restrainNextRes = False):
    """Set angle restraints on the specified nucleotides to match the Phenix pucker-specific restraints
    
    ARGUMENTS:
        molNum              - the molecule number
        chain               - the chain
        prevResNum          - the residue number of the previous residue (using Coot numbering)
        curResNum           - the residue number of the current residue (using Coot numbering)
        nextResNum          - the residue number of the next residue (using Coot numbering)
        glycosidicBaseAtom  - the name of the base nitrogen atom (N1 or N9) for the current residue
        prevPucker          - the pucker of the previous residue (as an int: 2 or 3)
        curPucker           - the pucker of the current residue (as an int: 2 or 3)
        nextPucker          - the pucker of the previous residue (as an int: 2 or 3)
    OPTIONAL ARGUMENTS:
        restrainNextRes     - whether we should add restraints for atoms of the next residue (other than the phosphate, which is always restrained)
                              Defaults to False
    RETURNS:
        None
    """
    
    restraintList = []
    linkRestraintList = []
    
    if prevResNum is not None:
        if prevPucker == 2:
            restraintList.append([prevResNum, " C3'", " O3'", 1.427, 0.015])
        elif prevPucker == 3:
            restraintList.append([prevResNum, " C3'", " O3'", 1.417, 0.015])
        
        linkRestraintList.append([prevResNum, " O3'", curResNum, " P  ", 1.607, 0.015])
    
    
    if curPucker == 2:
        restraintList.extend([[curResNum, " C1'", " C2'", 1.526, 0.015],
                              [curResNum, " C2'", " C3'", 1.525, 0.015],
                              [curResNum, " C3'", " C4'", 1.527, 0.015],
                              [curResNum, " C3'", " O3'", 1.427, 0.015],
                              [curResNum, " C4'", " C5'", 1.509, 0.015],
                              [curResNum, " C4'", " O4'", 1.454, 0.015],
                              [curResNum, " O4'", " C1'", 1.415, 0.015],
                              [curResNum, " C5'", " O5'", 1.424, 0.015],
                              [curResNum, " C2'", " O2'", 1.412, 0.015],
                              [curResNum, glycosidicBaseAtom, " C1'", 1.465, 0.015]])
                              #the P-O5' restraint is included as a monomer restraint
                              #In the actual Phenix restraint files, the glycosidic bond restraint appears to be
                              #pucker-specific only for purines.  With pyrimidines, the ideal value is 1.478 A
                              #regardless of pucker.  I'm assuming that this is a mistake.
    elif curPucker == 3:
        restraintList.extend([[curResNum, " C1'", " C2'", 1.529, 0.015],
                              [curResNum, " C2'", " C3'", 1.523, 0.015],
                              [curResNum, " C3'", " C4'", 1.521, 0.015],
                              [curResNum, " C3'", " O3'", 1.417, 0.015],
                              [curResNum, " C4'", " C5'", 1.508, 0.015],
                              [curResNum, " C4'", " O4'", 1.451, 0.015],
                              [curResNum, " O4'", " C1'", 1.412, 0.015],
                              [curResNum, " C5'", " O5'", 1.420, 0.015],
                              [curResNum, " C2'", " O2'", 1.420, 0.015],
                              [curResNum, glycosidicBaseAtom, " C1'", 1.475, 0.015]])
                              #the P-O5' restraint is included as a monomer restraint
    
    if nextResNum is not None:
        linkRestraintList.append([curResNum, " O3'", nextResNum, " P  ", 1.607, 0.015])
        
        if restrainNextRes:
            #print "Adding C5'-O5' restraint for residue", nextResNum, "pucker", nextPucker
            #from time import sleep; sleep(2)
            if nextPucker == 2:
                restraintList.append([nextResNum, " C5'", " O5'", 1.424, 0.015])
            elif nextPucker == 3:
                restraintList.append([nextResNum, " C5'", " O5'", 1.420, 0.015])
        #else:
        #    print "NOT adding C5'-O5' restraint for residue", nextResNum, "pucker", nextPucker
        #    from time import sleep; sleep(2)
        
    
    #apply the restraints
    for curRestraint in restraintList:
        (resNum, atom1, atom2, mean, sd) = curRestraint
        add_extra_bond_restraint(molNum, chain, resNum, "", atom1, "",
                                         chain, resNum, "", atom2, "",
                                         mean, sd)
        #print "add_extra_bond_restraint(", ", ".join(map(str, [molNum, chain, resNum, "", atom1, "", chain, resNum, "", atom2, "", mean, sd])), ")"
        
    for curRestraint in linkRestraintList:
        (resNum1, atom1, resNum2, atom2, mean, sd) = curRestraint
        add_extra_bond_restraint(molNum, chain, resNum1, "", atom1, "",
                                         chain, resNum2, "", atom2, "",
                                         mean, sd)
        #print "add_extra_bond_restraint(", ", ".join(map(str, [molNum, chain, resNum1, "", atom1, "", chain, resNum2, "", atom2, "", mean, sd])), ")"

    
    