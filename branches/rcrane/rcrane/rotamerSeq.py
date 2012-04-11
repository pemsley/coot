#!/usr/bin/env python
"""Determine the sequence of rotamers for a structure using a hidden Markov model"""

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

import os.path
from math import log

from pseudoPredictor import PseudoPredictor
from puckerList import puckerList, rotList, rotsByPucker
from safelog import ln, negInf

#initialize a PseudoPredictor object when this module is loaded
dataPath = os.path.dirname(os.path.abspath(__file__))
dataPath = os.path.join(dataPath, "data")
pseudoPredic = PseudoPredictor( thetaEta      = os.path.join(dataPath, "thetaEtaClusts.csv"),
                                pucker        = os.path.join(dataPath, "smoothedPuckerDist.csv"),
                                sugarDist     = os.path.join(dataPath, "sugarDists.csv"),
                                startPhosDist = os.path.join(dataPath, "startingPDists.csv"),
                                endPhosDist   = os.path.join(dataPath, "endingPDists.csv"))



def determineRotamerSeq(builtChain):
    """Determine the best sequence of rotamers for a structure.
    
    ARGUMENTS:
        builtChain      - a Chain object containing phosphate and base coordinates
    RETURNS:
        bestPath        - a list of the most likely rotamer for each suite
        predictedProbs  - the probability for each suite for each rotamer formatted as a list of dictionaries
                          predictedProbs[suiteNum][rotamer] = probability
    """
    
    #calculate the probabilities for each suite
    predictedProbs = []
    for curSuite in builtChain.suites():
        curProbs = pseudoPredic.calcProb( theta         = curSuite.theta(),
                                          eta           = curSuite.eta(),
                                          startPperp    = curSuite.startingPperp(),
                                          endPperp      = curSuite.endingPperp(),
                                          sugarDist     = curSuite.sugarDist(),
                                          startPhosDist = curSuite.startingPhosDist(),
                                          endPhosDist   = curSuite.endingPhosDist())
        
        predictedProbs.append(curProbs)
    
    #determine the best path using an HMM
    bestPath = rotamerHMM(predictedProbs)
    
    return (bestPath, predictedProbs)

def rotamerHMM(predictedProbs):
    """Given rotamer likelihoods for all suites, predict the most likely rotamer string using a Hidden Markov Model
    
    ARGUMENTS:
        predictedProbs  - the probability for each suite for each rotamer formatted as a list of dictionaries
                          predictedProbs[suiteNum][rotamer] = probability
    RETURNS:
        bestPath        - a list of the most likely rotamer for each suite
    """
    
    numSuites = len(predictedProbs)
    pathProbs = [{} for i in xrange(numSuites)]  #the log probability of having followed a given path (the delta or phi array)
    path      = [{} for i in xrange(numSuites)]  #the path followed for $pathProbs (the psi array)
    
    #initialize the pathProbs list
    for curRot in rotList:
        pathProbs[0][curRot] = ln(predictedProbs[0][curRot])
    
    for curPos in xrange(1, numSuites):  #for each suite
        for curRot in rotList:                     #for each rotamer
            
            #figure out what the best previous rotamer is for ending up at the current rotamer
            bestPrevRot = max(pathProbs[curPos-1], key = lambda prevRot: pathProbs[curPos-1][prevRot] + __transitionProb(prevRot, curRot))
            path[curPos][curRot]      = bestPrevRot
            pathProbs[curPos][curRot] = pathProbs[curPos-1][bestPrevRot] + __transitionProb(bestPrevRot, curRot) + ln(predictedProbs[curPos][curRot])

    #initialize bestPath to the appropriate length
    bestPath = [None] * numSuites 
    
    #figure out the best last position
    curPos = numSuites - 1
    bestPath[curPos] = max(pathProbs[curPos], key = lambda curRot: pathProbs[curPos][curRot])
    
    #follow the path back to figure out what the best path was
    for curPos in xrange(numSuites-1, 0, -1):
        bestPath[curPos-1] = path[curPos][bestPath[curPos]]
    
    return bestPath
    
    

def __transitionProb(startingRot, endingRot):
    """Calculate the log of the transition probability between two rotamers
    
    ARGUMENTS:
        startingRot - the rotamer to transition from
        endingRot   - the rotamer to transition to
    RETURNS:
        0 if the ending pucker of the starting rotamer is the same as the starting pucker of the ending rotamer
        negative infinity otherwise
    """
    
    if puckerList[startingRot][1] == puckerList[endingRot][0]:
        return 0
    else:
        return negInf


def determinePucker(pperp):
    """Predict only the pucker (used when the user has only built a single nucleotide)
    
    ARGUMENTS:
        pperp - the base-phosphate perpendicular (P-perp) distance
    RETURNS:
        3 if a C3'-endo sugar pucker is more likely, 2 otherwise
    NOTE:
        This function simply calls the calcPucker() function from the pseudoPredictor module.  The
        determinePucker() function exists only so the traceGui module doesn't have to separately import
        and initialize the pseudoPredictor module.
    """
    
    return pseudoPredic.calcPucker(pperp)

def determineAlternateConf(leadingPucker, endingPucker, suiteNum, predictedProbs):
    """Determine the best conformer for a suite given required starting and ending sugar puckers
    
    ARGUMENTS:
        leadingPucker   - the pucker of the starting sugar of the suite (either 2 or 3)
        endingPucker    - the pucker of the ending sugar of the suite (either 2 or 3)
        suiteNum        - the number of the suite
        predictedProbs  - the probability for each suite for each rotamer, as returned by determineRotamerSeq
    RETURNS:
        the most likely conformer
    """
    
    return max(rotsByPucker[leadingPucker][endingPucker], key = lambda rot: predictedProbs[suiteNum][rot])
