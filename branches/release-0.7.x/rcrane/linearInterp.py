#!/usr/bin/env python
"""A class to perform linear interpolation on data read from csv file."""

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

from math import floor

class LinearInterp:
    """A class to perform linear interpolation on data read from csv file."""
    
    def __init__(self, fileName):
        """Initialize the class and read data from the specified a file name.
        
        ARGUMENTS:
            fileName - the name of the CSV data file.
                       This CSV file should have two columns: x and y, as well as a header row (which will be ignored).
        RETURNS:
            an initialized LinearInterp object
        """
        
        self.__xStart   = None #the starting value for the x-axis
        self.__xStep    = None #the size of the increment along the x-axis
        self.__fileName = None #the input filename
        self.__yData    = []   #the y values read from fileName
        self.__fileName = fileName
        
        input = open(fileName, 'r')
        input.readline() #skip the header
        
        #read in the first line, which determines xStart
        curdata = input.readline().split(",")
        self.__xStart = float(curdata[0])
        self.__yData.append(float(curdata[1]))
        
        #read in the second line, which determines xStep
        curdata = input.readline().split(",")
        self.__xStep = float(curdata[0]) - float(self.__xStart)
        self.__yData.append(float(curdata[1]))
        
        ep = self.__xStep / 10000 #since the step size is a float, we can't test that each step size is exactly right
                                #instead, we make sure that the step size is within .01% of the correct size
                                #using this epsilon value
        
        #read in the rest of the data file
        for curline in input:
            
            #skip blank lines
            if len(curline.rstrip()) == 0:
                continue
            curdata = curline.split(",")
            
            #make sure that step size is correct
            xVal = self.__xStart + self.__xStep * len(self.__yData)
            if (xVal - ep <= curdata[0]) and (curdata[0] <= xVal + ep):
                print "Incorrect step size in file " + self.__fileName
            
            self.__yData.append(float(curdata[1]))
        input.close()
    
    
    def interp(self, xVal):
        """Determine the interpolated (or extrapolated) y value for the given x value
        
        ARGUMENTS:
            xVal - the x value to determine the y value for
        RETURNS:
            yVal - the interpolated (or extrapolated) y value
        """
        
        xPos = int(floor( (xVal-self.__xStart) / self.__xStep ))
        
        #determine if this is an interpolation or an extrapolation
        yVal = None
        if xPos < 0:
            #extrapolation to lower x values
            yStep = self.__yData[0] - self.__yData[1]
            numSteps = (xVal-self.__xStart) / self.__xStep 
            yVal = self.__yData[0] - numSteps * yStep
            #print "Extrapolating to lower x value (" + self.__fileName + "): (%(xVal)f, %(yVal)f)" %vars()
        elif xPos > len(self.__yData)-1:
            #extrapolation to higher x values
            yStep = self.__yData[-1] - self.__yData[-2]
            numSteps = (xVal-self.__xStart) / self.__xStep - (len(self.__yData) - 1)
            yVal = self.__yData[-1] + numSteps * yStep
            #print "Extrapolating to higher x value (" + self.__fileName + "): (%(xVal)f, %(yVal)f)" %vars()
        else:
            #interpolation
            remainder = xVal - self.__xStart - self.__xStep * xPos
            yVal = self.__yData[xPos] * ( 1 - remainder/self.__xStep) + self.__yData[xPos+1] * (remainder/self.__xStep)
            #print "Interpolating (" + self.__fileName + "): (%(xVal)f, %(yVal)g)" %vars()
        return yVal
