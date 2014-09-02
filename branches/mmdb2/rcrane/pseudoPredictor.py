#!/usr/bin/env python
"""Predict rotamers from low resolution information"""

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

from thetaEtaPredictor import ThetaEtaPredictor
from smoothProb import SmoothProb, SmoothProbError
from puckerList import rotList, puckerListFullName as puckerList

class PseudoPredictor:
    """Predict rotamers from low resolution information"""
    
    def __init__(self, thetaEta = None, pucker = None, sugarDist = None, startPhosDist = None, endPhosDist = None):
        """Create and initialize a PseudoPredictor object
        OPTIONAL ARGUMENTS:
            Note that at least one of these arugments must be provided
            thetaEta        - A file contianing training data for theta eta Gaussian clusters
            pucker          - A file contianing training data for P-perp distances
            sugarDist       - A file contianing training data for sugar-sugar (C1'-C1') distances
            startPhosDist   - A file contianing training data for starting phosphate-phosphate distances
            endPhosDist     - A file contianing training data for ending phosphate-phosphate distances
        RETURNS:
            an initialized PseudoPredictor object
        """           
        
        #create a dictionary of the input filenames
        inputFiles = locals()
        del inputFiles["self"]
        
        self.__inputFiles = inputFiles
        self.__thetaEtaPredic      = None
        self.__puckerPredic        = None
        self.__sugarDistPredic     = None
        self.__startPhosDistPredic = None
        self.__endPhosDistPredic   = None
        
        if thetaEta is not None:
            self.__thetaEtaPredic = ThetaEtaPredictor(thetaEta)
        
        if pucker is not None:
            self.__puckerPredic = SmoothProb(pucker)
        
        if sugarDist is not None:
            self.__sugarDistPredic = SmoothProb(sugarDist)
        
        if startPhosDist is not None:
            self.__startPhosDistPredic = SmoothProb(startPhosDist)
        
        if endPhosDist is not None:
            self.__endPhosDistPredic = SmoothProb(endPhosDist)
    
    def calcProb(self, rot = None, thetaEta = None, theta = None, eta = None, startPperp = None, endPperp = None, sugarDist = None, startPhosDist = None, endPhosDist = None):
        """Calculate the probability of each rotamer for a suite from all given information
        
        OPTIONAL ARGUMENTS:
            rot - which rotamer to return the probability for
                if not provided, a hash of all probabilities is returned
            Note that at least one of the following arugments must be provided:
            thetaEta      - a list of the (theta, eta) coordinates of the suite
            theta         - the theta coordinates of the suite (ignored if thetaEta is given)
            eta           - the eta coordinates of the suite (ignored if thetaEta is given)
            startPperp    - the starting base-phosphate perpendicular (P-perp) distance
            endPperp      - the ending base-phosphate perpendicular (P-perp) distance
            sugarDist     - the C1'-C1' distance for the suite
            startPhosDist - the phosphate-phosphate distance for the starting nucleotide of the suite
            endPhosDist   - the phosphate-phosphate distance for the ending nucleotide of the suite
        RETURNS:
            if rot is given, returns the probability that the described suite is of rotamer rot
            if rot is not given, returns a hash in the form rotamer => probability
        """
        
        #if thetaEta is given, use those values instead of theta and eta
        if thetaEta is not None:
            (theta, eta) = thetaEta
        
        #initialize variables to None
        (thetaEtaProbs, startPuckerProbs, endPuckerProbs, sugarDistProbs, startPhosDistProbs, endPhosDistProbs) = [None] * 6
        
        #calculate individual probabilites using the given information
        if self.__thetaEtaPredic is not None and theta is not None and eta is not None:
            thetaEtaProbs = self.__thetaEtaPredic.calcProb(theta, eta)
        
        if self.__puckerPredic is not None and startPperp is not None:
            try:
                startPuckerProbs = self.__puckerPredic.calcProb(startPperp)
            except SmoothProbError:
                startPuckerProbs = None
        
        if self.__puckerPredic is not None and endPperp is not None:
            try:
                endPuckerProbs = self.__puckerPredic.calcProb(endPperp)
            except SmoothProbError:
                endPuckerProbs = None
        
        if self.__sugarDistPredic is not None and sugarDist is not None:
            try:
                sugarDistProbs = self.__sugarDistPredic.calcProb(sugarDist)
            except SmoothProbError:
                sugarDistProbs = None
        
        if self.__startPhosDistPredic is not None and startPhosDist is not None:
            try:
                startPhosDistProbs = self.__startPhosDistPredic.calcProb(startPhosDist)
            except SmoothProbError:
                startPhosDistProbs = None
        
        if self.__endPhosDistPredic is not None and endPhosDist is not None:
            try:
                endPhosDistProbs = self.__endPhosDistPredic.calcProb(endPhosDist)
            except SmoothProbError:
                endPhosDistProbs = None
        
        
        #calculate the overall probability for each rotamer
        probs = {}
        for currot in rotList:
            curprob = 1
            
            if thetaEtaProbs      is not None: curprob *= thetaEtaProbs[currot]
            if sugarDistProbs     is not None: curprob *= sugarDistProbs[currot]
            if startPhosDistProbs is not None: curprob *= startPhosDistProbs[currot]
            if endPhosDistProbs   is not None: curprob *= endPhosDistProbs[currot]
            
            if startPuckerProbs   is not None: curprob *= startPuckerProbs[puckerList[currot][0]]
            if endPuckerProbs     is not None: curprob *= endPuckerProbs[puckerList[currot][1]]
            
            probs[currot] = curprob
        
        #scale the raw probabilities so they add up to 1
        totalProb = sum(probs.values())
        
        try:
            probs = dict([(currot, curprob/totalProb) for (currot, curprob) in probs.iteritems()])
        except ZeroDivisionError:
            #if totalProb is 0, then the structure is very, very bad
            #our only reasonable options are to give up, or to give all rotamers equal likelihood
            #we give all rotamers equal likelihood so that the user at least has a chance of improving the structure
            equalProb = 1.0/len(rotList)
            probs = dict([(currot, equalProb) for currot in rotList])
        
        #return either the desired probability or the dictionary containing all probabilities
        if rot is None:
            return probs
        else:
            return probs[rot]
    
    def calcPucker(self, pperp):
        """Calculate the most likely pucker of a sugar given only the base-phosphate perpendicular distance
        
        ARGUMENTS:
            pperp - the base-phosphate perpendicular (P-perp) distance
        RETURNS:
            3 if a C3'-endo sugar pucker is more likely, 2 otherwise
        """
        
        if self.__puckerPredic is None:
            return None
        else:
            puckerProbs = self.__puckerPredic.calcProb(pperp)
            if puckerProbs["C3'"] > puckerProbs["C2'"]:
                return 3
            else:
                return 2