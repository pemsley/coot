#!/usr/bin/env python
"""A class for storing information about individual nucleotides."""

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

from strucCalc import torsion, dist, distToLine

SUGARATOM = "C1'"

class Nucleotide:
    """A class for storing information about individual nucleotides (which may be part of a Chain object)."""
    
    def __init__(self, type, atoms, resNum = None):
        """Initialize a Nucleotide object
        
        ARGUMENTS:
            type  - the nucleotide type (i.e. A, C, G, or U)
            atoms - a hash of atom coordinates (in the form atomName: [x, y, z])
        OPTIONAL ARGUMENTS:
            resNum - the residue number (possibly including an insertion code)
        RETURNS:
            an initialized Nucleotide object
        """
        
        self.type = type
        self.atoms = atoms
        self.resNum = resNum
        self.chain = None
        self.seqNum = None
    
    
    def prevNuc(self):
        """If this nucleotide is part of a chain, return the previous nucleotide in the chain
        
        ARGUMENTS:
            None
        RETURNS:
            The Nucleotide object of the previous nucleotide if there is one
            None if this is the first nucleotide of the chain,or if this nucleotide is not part of a chain.
        """
        
        if self.seqNum > 0:
            return self.chain.nucleotides[self.seqNum-1]
        else:
            return None
    
    
    def nextNuc(self):
        """If this nucleotide is part of a chain, return the next nucleotide in the chain
        
        ARGUMENTS:
            None
        RETURNS:
            The Nucleotide object of the next nucleotide if there is one
            None if this is the last nucleotide of the chain or if this nucleotide is not part of a chain.
        """
        
        if self.seqNum < self.chain.len():
            return self.chain.nucleotides[self.seqNum+1]
        else:
            return None
    
    
    def connectedToPrev(self):
        """Is this nucleotide connected to the previous nucleotide of the chain?
        
        ARGUMENTS:
            None
        RETURNS:
            True if this nucleotide is connected to the previous nucleotide of the chain
            False if it isn't or if this nucleotide is not part of a chain.
        NOTE:
            Because the chain module does not currently support unconnected nucleotides, this function only
            tells you if this is the first nucleotide of the chain or not.
        """
        
        return (self.seqNum > 0)
    
    
    def connectedToNext(self):
        """Is this nucleotide connected to the next nucleotide of the chain?
        
        ARGUMENTS:
            None
        RETURNS:
            True if this nucleotide is connected to the previous nucleotide of the chain
            False if it isn't
            Note that the behavior of this function (unlike that of connectedToPrev()) is not defined
                if this nucleotide is not part of a chain.  (Currently, it raises an error.)
        NOTE:
            Because the chain module does not currently support unconnected nucleotides, this function only
            tells you if this is the last nucleotide of the chain or not.
        """
        
        return (self.seqNum < (self.chain.len()-1))
    
    
    def hasAtom(self, atomName):
        """Does this nucleotide contain the specified atom?
        
        ARGUMENTS:
            atomName - the name of the atom to check for
        RETURNS:
            True if this nucleotide contains atomName, False otherwise
        """
        
        return self.atoms.has_key(atomName)
    
    def eta(self):
        """Calculate the eta pseudotorsion for the current nucleotide
        
        ARGUMENTS:
            None
        RETURNS:
            The value of eta if it is defined, None otherwise.
        NOTE:
            Calculating eta requires the sugar atom of the previous nucleotide and the phosphate atom
            of the next nucleotide.  As a result, behaviour of this function is undefined if this nucleotide
            is not part of a chain.  (Currently, it raises an error.)
        """
        
        if (self.connectedToPrev() and self.connectedToNext()
              and self.prevNuc().hasAtom(SUGARATOM)
              and self.hasAtom("P")
              and self.hasAtom(SUGARATOM)
              and self.nextNuc().hasAtom("P")):
            
            return torsion(self.prevNuc().atoms[SUGARATOM],
                               self.atoms["P"],
                               self.atoms[SUGARATOM],
                               self.nextNuc().atoms["P"])
        else:
            return None
    
    
    def theta(self):
        """Calculate the theta pseudotorsion for the current nucleotide
        
        ARGUMENTS:
            None
        RETURNS:
            The value of theta if it is defined, None otherwise.
        NOTE:
            Calculating theta requires the phosphate and sugar atom of the next nucleotide.
            As a result, behaviour of this function is undefined if this nucleotide
            is not part of a chain.  (Currently, it raises an error.)
        """
        
        if (self.connectedToNext()
              and self.hasAtom("P")
              and self.hasAtom(SUGARATOM)
              and self.nextNuc().hasAtom("P")
              and self.nextNuc().hasAtom(SUGARATOM)):
            
            return torsion(self.atoms["P"],
                               self.atoms[SUGARATOM],
                               self.nextNuc().atoms["P"],
                               self.nextNuc().atoms[SUGARATOM])
        else:
            return None
    
    
    def pperp(self):
        """Calculate the base-phosphate perpendicular distance (P-perp distance) for the current nucleotide
        
        ARGUMENTS:
            None
        RETURNS:
            The P-perp distance if it is defined, None otherwise
        NOTE:
            Calculating P-perp requires the phosphate of the next nucleotide.  As a result, behaviour of this
            function is undefined if this nucleotide is not part of a chain.  (Currently, it raises an error.)
        """
        
        
        #determine the coordinates of the base nitrogen
        if self.hasAtom("N9"):
            baseN = self.atoms["N9"]
        elif self.hasAtom("N1"):
            baseN = self.atoms["N1"]
        else:
            baseN = None
        
        if baseN is not None and self.hasAtom(SUGARATOM) and self.connectedToNext() and self.nextNuc().hasAtom("P"):
            return distToLine(self.nextNuc().atoms["P"],
                              self.atoms[SUGARATOM],
                              baseN)
        else:
            return None
    
    
    def phosDist(self):
        """Calculate the phosphate-phosphate distance
        
        ARGUMENTS:
            None
        RETURNS:
            The distance from the phosphate of this nucleotide to the phosphate of the next nucleotide
            None if no next phosphate if present
        NOTE:
            Behavior of this function is undefined if this nucleotide is not part of a chain.
            (Currently, it raises an error.)
        """
        
        if self.connectedToNext() and self.hasAtom("P") and self.nextNuc().hasAtom("P"):
            return dist(self.atoms["P"], self.nextNuc().atoms["P"])
        else:
            return None
    
    
    def startingSugarDist(self):
        """Calculate the sugar-sugar distance for the starting suite of this nucleotide
        
        ARGUMENTS:
            None
        RETURNS:
            the distance from the sugar atom (i.e. C1') of the previous nucleotide to the sugar atom (C1') of this nucleotide
            None if there is no previous nucleotide
        NOTE:
            Behavior of this function is undefined if this nucleotide is not part of a chain.
            (Currently, it returns None.)
        """
        
        if self.connectedToPrev() and self.hasAtom(SUGARATOM) and self.prevNuc().hasAtom(SUGARATOM):
            return dist(self.atoms[SUGARATOM], self.prevNuc().atoms[SUGARATOM])
        else:
            return None
    
    
    def endingSugarDist(self):
        """Calculate the sugar-sugar distance for the ending suite of this nucleotide
        
        ARGUMENTS:
            None
        RETURNS:
            the distance from the sugar atom (i.e. C1') of this nucleotide to the sugar atom (C1') of the next nucleotide
            None if there is no previous nucleotide
        NOTE:
            Behavior of this function is undefined if this nucleotide is not part of a chain.
            (Currently, it raises an error.)
        """
        
        if self.connectedToNext():
            return self.nextNuc().startingSugarDist()
        else:
            return None
    
    def __str__(self):
        """Create a string representation of the object"""
        
        retVal  = ("Nuc Num:  " + str(self.resNum) + "\n")
        retVal += ("Seq Num:  " + str(self.seqNum) + "\n")        
        retVal += ("Nuc type: " + str(self.type) + "\n")        
        retVal += "Atoms: \n"
        for (name, coords) in self.atoms.iteritems():
            retVal += ("\t" + name + ":  " + str(coords) + "\n")
        return retVal
    
    
    def alpha(self):
        """Calculate the alpha torsion
        
        ARGUMENTS:
            None
        RETURNS:
            The value of alpha if it is defined, None otherwise.
        """
        
        if (self.connectedToPrev()
              and self.prevNuc().hasAtom("O3'")
              and self.hasAtom("P")
              and self.hasAtom("O5'")
              and self.hasAtom("C5'")):
            
            return torsion(self.atoms["C5'"],
                           self.atoms["O5'"],
                           self.atoms["P"],
                           self.prevNuc().atoms["O3'"])
                
        else:
            return None

    def beta(self):
        """Calculate the beta torsion
        
        ARGUMENTS:
            None
        RETURNS:
            The value of beta if it is defined, None otherwise.
        """
        
        if   (self.hasAtom("P")
          and self.hasAtom("O5'")
          and self.hasAtom("C5'")
          and self.hasAtom("C4'")):
            
            return torsion(self.atoms["P"],
                           self.atoms["O5'"],
                           self.atoms["C5'"],
                           self.atoms["C4'"])
                
        else:
            return None
        
    def gamma(self):
        """Calculate the gamma torsion
        
        ARGUMENTS:
            None
        RETURNS:
            The value of gamma if it is defined, None otherwise.
        """
        
        if   (self.hasAtom("O5'")
          and self.hasAtom("C5'")
          and self.hasAtom("C4'")
          and self.hasAtom("C3'")):
            
            return torsion(self.atoms["O5'"],
                           self.atoms["C5'"],
                           self.atoms["C4'"],
                           self.atoms["C3'"])
                
        else:
            return None
        
    def delta(self):
        """Calculate the delta torsion
        
        ARGUMENTS:
            None
        RETURNS:
            The value of delta if it is defined, None otherwise.
        """
        
        if   (self.hasAtom("C5'")
          and self.hasAtom("C4'")
          and self.hasAtom("C3'")
          and self.hasAtom("O3'")):
            
            return torsion(self.atoms["C5'"],
                           self.atoms["C4'"],
                           self.atoms["C3'"],
                           self.atoms["O3'"])
                
        else:
            return None
    
    def epsilon(self):
        """Calculate the epsilon torsion
        
        ARGUMENTS:
            None
        RETURNS:
            The value of epsilon if it is defined, None otherwise.
        """
        
        if (self.connectedToNext()
              and self.hasAtom("C4'")
              and self.hasAtom("C3'")
              and self.hasAtom("O3'")
              and self.nextNuc().hasAtom("P")):
            
            return torsion(self.atoms["C4'"],
                           self.atoms["C3'"],
                           self.atoms["O3'"],
                           self.nextNuc().atoms["P"])
                
        else:
            return None
        
    
    def zeta(self):
        """Calculate the zeta torsion
        
        ARGUMENTS:
            None
        RETURNS:
            The value of zeta if it is defined, None otherwise.
        """
        
        if (self.connectedToNext()
              and self.hasAtom("C3'")
              and self.hasAtom("O3'")
              and self.nextNuc().hasAtom("P")
              and self.nextNuc().hasAtom("O5'")):
            
            return torsion(self.atoms["C3'"],
                           self.atoms["O3'"],
                           self.nextNuc().atoms["P"],
                           self.nextNuc().atoms["O5'"])
                
        else:
            return None
    