#!/usr/bin/env python
"""A class to access torsion mean and standard deviation information for the backbone rotamers"""

# Copyright 2010 Kevin Keating
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

class RotDataLoader:
    """A class to access torsion mean and standard deviation information for the backbone rotamers"""
    
    def __init__(self, filename):
        """Initialize a RotDataLoader object by reading data from a csv file
        
        ARGUMENTS:
            filename - a csv file containing data on the rotamers
        RETURNS:
            an initialized RotDataLoader object
        """
        
        input = open(filename, 'r')
        input.readline() #skip the header
        
        self.rotData = {}
        for curline in input:
            curdata = curline.split(",")
            rot = curdata[0].strip()
            if rot == "": continue #skip blank lines
            self.rotData[rot] = map(float, curdata[3:10]+curdata[13:20])
        input.close()
            
    
    
    #functions to access each mean and standard deviation
    def prevDeltaMean(self, rot): return self.rotData[rot][0]
    def epMean       (self, rot): return self.rotData[rot][1]
    def zetaMean     (self, rot): return self.rotData[rot][2]
    def alphaMean    (self, rot): return self.rotData[rot][3]
    def betaMean     (self, rot): return self.rotData[rot][4]
    def gammaMean    (self, rot): return self.rotData[rot][5]
    def curDeltaMean (self, rot): return self.rotData[rot][6]
    
    def prevDeltaSD  (self, rot): return self.rotData[rot][7]
    def epSD         (self, rot): return self.rotData[rot][8]
    def zetaSD       (self, rot): return self.rotData[rot][9]
    def alphaSD      (self, rot): return self.rotData[rot][10]
    def betaSD       (self, rot): return self.rotData[rot][11]
    def gammaSD      (self, rot): return self.rotData[rot][12]
    def curDeltaSD   (self, rot): return self.rotData[rot][13]
