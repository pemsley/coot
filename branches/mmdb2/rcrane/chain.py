#!/usr/bin/env python
"""A class for storing chains of nucleotides."""

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

from suite import Suite

class Chain:
    """A class for storing chains of nucleotides (which must all be instances of the nucleotide class)."""
    
    def __init__(self, nucs = None):
        """Initialize a Chain object.
        
        OPTIONAL ARGUMENTS:
            nucs - a list of nucleotide objects to add to the chain
                   If not given, then an empty chain is created
        RETURNS:
            an initialized Chain object
        """
        
        #self.__nucNums = []
        self.__resNumDict = {}
        
        if nucs is not None:
            for (curSeqNum, curNuc) in enumerate(nucs):
                curNuc.chain = self
                curNuc.seqNum = curSeqNum
                self.__resNumDict[str(curNuc.resNum)] = curSeqNum
                #self.__nucNums.append(str(curNuc.resNum))
        else:
            nucs = []
        
        self.nucleotides = nucs
    
    @property
    def nucs(self):
        """A synonym for Chain.nucleotides"""
        return self.nucleotides
    
    
    def addNuc(self, nuc):
        """Add a nucleotide to the end of the chain
        
        ARGUMENTS:
            nuc - a nucleotide object to add to the chain
        RETURNS:
            None
        EFFECTS:
            Adds nuc to the chain
        """
        
        nuc.chain = self
        nuc.seqNum = self.numNucs()
        self.nucleotides.append(nuc)
        self.__resNumDict[str(nuc.resNum)] = nuc.seqNum
        #self.__nucNums.append(str(nuc.resNum))
    
    
    def addNuc5p(self, nuc):
        """Add a nucleotide to the start (i.e. 5' end) of the chain
        
        ARGUMENTS:
            nuc - a nucleotide object to add to the chain
        RETURNS:
            None
        EFFECTS:
            Adds nuc to the chain
        """
        
        nuc.chain = self
        nuc.seqNum = 0
        
        #update the sequence number of the nucleotides
        for curNuc in self.nucleotides:
            curNuc.seqNum += 1
        
        #update the residue number lookup dictionary
        for (curKey, curVal) in self.__resNumDict.iteritems():
            self.__resNumDict[curKey] = curVal + 1
        
        self.nucleotides.insert(0, nuc)
    
    
    def suites(self):
        """A generator for suites in the chain
        
        ARGUMENTS:
            None
        YIELDS:
            The next suite of the chain
        """
        
        #if the first nucleotide doesn't have a C1', then assume it's not a complete nucleotide and skip it
        if self.nucleotides[0].hasAtom("C1'"):
            startingNuc = 0
        else:
            startingNuc = 1
        
        for i in xrange(startingNuc, len(self.nucleotides)-2):
            yield Suite(self.nucleotides[i], self.nucleotides[i+1])
    
    
    def len(self):
        """The length of the chain in nucleotides
        
        ARGUMENTS:
            None
        RETURNS:
            The length of the chain in nucleotides
        """
        
        return len(self.nucleotides)
    
    
    def numNucs(self):
        """The length of the chain in nucleotides
        
        ARGUMENTS:
            None
        RETURNS:
            The length of the chain in nucleotides
        """
        
        return len(self.nucleotides)
    
    
    def numSuites(self):
        """The length of the chain in suites
        
        ARGUMENTS:
            None
        RETURNS:
            The length of the chain in suites
        """
        
        return len(self.nucleotides) - 1
    
    
    def resIndex(self, resNum):
        """Find the index (i.e. position in the chain) of a given residue number
        
        ARGUMENTS:
            resNum - the residue number to find (which may include an insertion code)
        RETURNS:
            the index of the specified residue number in this chain
        """
        
        return self.__resNumDict[str(resNum)]
    
    
    def resNumList(self):
        return [nuc.resNum for nuc in self.nucleotides]
        
    def resTypeList(self):
        return [nuc.type for nuc in self.nucleotides]
    
    def __str__(self):
        """Create a string representation of the object"""
        return "\n\n".join(map(str, self.nucleotides))