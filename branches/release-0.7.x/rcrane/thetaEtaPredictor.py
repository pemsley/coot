#!/usr/bin/env python
"""Predict rotamers from (theta, eta) values"""

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

from math import radians, sin, cos, pi, sqrt, exp

DELIM = "," #the delimeter used in the input file (i.e. csv)

#constants used in the cluster data list
STARTPUCKER = 0
ENDPUCKER   = 1
CENTER      = 2
VARIANCE    = 3
ROTATION    = 4

class ThetaEtaPredictor:
    """Predict rotamers from (theta, eta) values
    
    NOTE:
        This prediction is based in-part on the soft K-means clustering algorithm from:
        MacKay, D.J.C. 2002. Maximum Likelihood and Clustering. In Information Theory, Inference & Learning
            Algorithms, pp. 300-306,  Cambridge University Press, New York, NY.
    """
    
    def __init__(self, fileName):
        """Create and initialize a ThetaEtaPredictor object
        
        ARGUMENTS:
            inputFile - a csv file containing training data for theta-eta Gaussian clusters
        RETURNS:
            a ThetaEtaPredictor object initialized with the cluster data from inputFile
        """
        
        input = open(fileName, 'r')
        input.readline() #ignore the header
        
        clusters = {}
        
        for curline in input:
            data = curline.split(DELIM)
            
            rotName     = data[0]
            startPucker = data[1]
            endPucker   = data[2]
            center      = map(float, data[3:5])
            variance    = map(float, data[5:7])
            rotation    = float(data[7])
           
            clusters[rotName] = [startPucker, endPucker, center, variance, rotation]
        
        input.close()
        
        self.__clusters = clusters
    
    
    def calcProb(self, theta, eta, rot = None):
        """Calculate the probability of each rotamer for a given (theta, eta) value
        
        ARGUMENTS:
            theta  - the theta value of the point
            eta    - the eta value of the point
        OPTIONAL ARGUMENTS
            rot    - a specific rotamer
        RETURNS:
            if rot is not given, a dictionary of {rotamer => probability}
            if rot is given, the probability that the rotamer is rot
        """
        
        probs = {}
        
        #calculate the raw (unscaled) probability for each rotamer
        for currot in self.__clusters.keys():
            probs[currot] = self.__calcRawProb(theta, eta, currot)
        
        #scale the raw probabilities so they add up to 1
        totalProb = sum(probs.values())
        probs = dict([(currot, curprob/totalProb) for (currot, curprob) in probs.iteritems()])
        
        #return either the desired probability or the dictionary containing all probabilities
        if rot is None:
            return probs
        else:
            return probs[rot]
    
    
    def __calcRawProb(self, theta, eta, rot):
        """Calculate the responsibility of a given rotamer cluster for a given point based only on (theta, eta) values
        
        ARGUMENTS:
            theta  - the theta value of the point
            eta    - the eta value of the point
            rot    - the rotamer
        RETURNS:
            the unscaled responsibility of rot for (theta, eta)
        """
        
        center   = self.__clusters[rot][CENTER]
        variance = self.__clusters[rot][VARIANCE]
        angle    = self.__clusters[rot][ROTATION]
        point = [theta, eta]
        
        #first, translate the point
        x = subtractCoords(point[0], center[0])
        y = subtractCoords(point[1], center[1])
        
        #rotate the point by the opposite of the angle
        angle = radians(angle)
        point[0] = x*cos(angle) - y*sin(angle)
        point[1] = x*sin(angle) + y*cos(angle)
        
        #the following caculation is based on the soft K-means clustering algorithm from
        #MacKay, D.J.C. 2002. Maximum Likelihood and Clustering. In Information Theory, Inference & Learning
        #    Algorithms, pp. 300-306,  Cambridge University Press, New York, NY.
        
        product = 1
        sum = 0
        for curaxis in (0, 1):
            product *= sqrt(2*pi * variance[curaxis])
            sum     += (point[curaxis]**2 / (2*variance[curaxis]))
        
        return 1.0/product * exp(-sum)


def subtractCoords(num1, num2):
    """Calculate the minimum distance between two measures, keeping the appropriate sign
    
    ARGUMENTS:
        num1 - the first measurement
        num2 - the second measurement
    RETURNS:
        the distance between these two angles of no more than 180 degrees
    """
    
    return (num1 - num2 + 180) % 360 - 180
