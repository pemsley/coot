#!/usr/bin/env python
"""A class for storing data about a suite of RNA"""

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

class Suite:
    """A class for storing data about a suite of RNA
    
    NOTE:
        No actual data is stored in the Suite object other than a link to the starting and
        ending nucleotides.  All functions defined here simply call the appropriate
        Nucleotide object function.
    """
    
    
    def __init__(self, startingNuc, endingNuc):
        """Initialize a Suite object
        
        ARGUMENTS:
            startingNuc - a Nucleotide object representing the starting nucleotide of the suite
            endingNuc   - a Nucleotide object representing the ending nucleotide of the suite
        RETURNS:
            an initialized Suite object
        """
        
        self.startingNuc = startingNuc
        self.endingNuc   = endingNuc
    
    def sugarDist(self):
        """Calculate the sugar distance for the suite"""
        return self.endingNuc.startingSugarDist()
    
    def theta(self):
        """Calculate theta for the suite"""
        return self.startingNuc.theta()
    
    def eta(self):
        """Calculate eta for the suite"""
        return self.endingNuc.eta()
    
    def startingPperp(self):
        """Calculate the starting base-phosphate perpendicular (P-perp) distance for the suite"""
        return self.startingNuc.pperp()
        
    def endingPperp(self):
        """Calculate the ending base-phosphate perpendicular (P-perp) distance for the suite"""
        return self.endingNuc.pperp()
    
    def startingPhosDist(self):
        """Calculate the starting phosphate-phosphate distance for the suite"""
        return self.startingNuc.phosDist()
        
    def endingPhosDist(self):
        """Calculate the ending phosphate-phosphate distance for the suite"""
        return self.endingNuc.phosDist()
    
    def suiteNum(self):
        """The suite number"""
        return self.endingNuc.resNum
