#!/usr/bin/env python
"""Calculate the probabilities of various models for a given distance measure using linear interpolation or extrapolation"""

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

from math import floor


import sys
#from pprint import pprint

#from linearInterp import LinearInterp

DELIM = ","

class SmoothProb:
    """Calculate the probabilities of various models for a given distance measure using linear interpolation or extrapolation"""
    
    def __init__(self, fileName):
        """Creates and initialize a smoothProb object
        
        ARGUMENTS:
            inputFile - a file containing smoothed training data
        RETURNS:
            an initialized SmoothProb object
        """
        
        self.__xStart   = None #the starting value for the x-axis
        self.__xStep    = None #the size of the increment along the x-axis
        self.__fileName = None #the input filename
        self.__yData    = {}   #the y values read from fileName
        self.__yLen     = None #the number of Y values 
        self.__labels   = None #the names of the columns
        self.__fileName = fileName
        
        input = open(fileName, 'r')
        
        #read in the header line, which contains the labels for our probabilities
        labels = input.readline().split(DELIM)
        labels.pop(0) #we don't care about the header for the first column, since that's just distance
        labels = [x.strip() for x in labels] #get rid of any leading and trailing spaces (at the very least, get rid of the \n at the end of the line)
        self.__labels = labels
        
        #read in the first line, which determines xStart
        curdata = input.readline().split(DELIM)
        self.__xStart = float(curdata.pop(0))
        #create arrays for all y values
        for curindex, curvalue in enumerate(curdata):
            self.__yData[labels[curindex]] = [float(curvalue)]
        
        #read in the second line, which determines xStep
        curdata = input.readline().split(DELIM)
        self.__xStep = float(curdata.pop(0)) - float(self.__xStart)
        for curindex, curvalue in enumerate(curdata):
            self.__yData[labels[curindex]].append(float(curvalue))
        
        ep = self.__xStep / 10000 #since the step size is a float, we can't test that each step size is exactly
                                  #right without running into rounding issues.
                                  #Instead, we make sure that the step size is within .01% of the correct size
                                  #using this epsilon value
        
        #go through the rest of the input file
        for curline in input:
            #skip blank lines
            if len(curline.rstrip()) == 0:
                continue
            curdata = curline.split(DELIM)
            
            #make sure that step size is correct
            xVal = self.__xStart + self.__xStep * len(self.__yData)
            curdist = curdata.pop(0)
            if (xVal - ep <= curdist) and (curdist <= xVal + ep):
                raise "Incorrect smoothProb step size in file " + self.__fileName
            
            for curindex, curvalue in enumerate(curdata):
                self.__yData[labels[curindex]].append(float(curvalue))
        
        self.__yLen = len(self.__yData[labels[0]])
        
        input.close()
        
    
    def calcProb(self, xVal, label = None):
        """calculate the probabilities for a given distance
        
        ARGUMENTS:
            dist   - the distance
        OPTIONAL ARGUMENTS:
            label  - what model to calculate the probability for
                     must match a label in the first line of the input file
        RETURNS:
            if label is given, the probability that the model label is correct
            if label is not given, a dictionary containing probabilities for all models
        """
        
        xPos = int(floor( (xVal-self.__xStart) / self.__xStep ))
        
        #go through each label and determine a raw probability
        yDict = {}
        if xPos < 0:
            #if we're extrapolating down
            numSteps = (xVal-self.__xStart) / self.__xStep
            for curlabel in self.__labels:
                yStep = self.__yData[curlabel][0] - self.__yData[curlabel][1]
                yDict[curlabel] = self.__yData[curlabel][0] - numSteps * yStep
        elif xPos > self.__yLen-1:
            #if we're extrapolating up
            numSteps = (xVal-self.__xStart) / self.__xStep - (self.__yLen - 1)
            for curlabel in self.__labels:
                yStep = self.__yData[curlabel][-1] - self.__yData[curlabel][-2]
                yDict[curlabel] = self.__yData[curlabel][-1] + numSteps * yStep
        else:
            #if we're interpolating
            remainder = xVal - self.__xStart - self.__xStep * xPos
            for curlabel in self.__labels:
                yDict[curlabel] = self.__yData[curlabel][xPos] * ( 1 - remainder/self.__xStep) + self.__yData[curlabel][xPos+1] * (remainder/self.__xStep)
        
        #scale the raw probabilities so they add up to 1
        totalProb = sum(yDict.values())
        
        try:
            yDict = dict([(curlabel, curval/totalProb) for (curlabel, curval) in yDict.iteritems()])
        except ZeroDivisionError:
            #if totalProb = 0, then we're trying to make predictions with *very* bad data
            raise SmoothProbError("Underflow error: all probabilities equal to zero.")
        
        #return either the desired probability or the hash containing all probabilities
        if label is None:
            return yDict
        else:
            return yDict[label]


class SmoothProbError(Exception):
    #this error gets raised when all probabilities are 0
    def __init__(self, description):
        self.value = description
    def __str__(self):
        return repr(self.value)