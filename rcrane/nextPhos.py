#!/usr/bin/env python
"""A class to determine the coordinates of the next phosphate atom in electron density"""

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

#from pprint import pprint
from math import radians, pi, cos, sin, sqrt

from linearInterp import LinearInterp
from strucCalc import dist, angle, crossProd, scalarProd, minus, plus, magnitude, torsion, rotateAtoms
from safelog import ln
from stats import median, lowerQuartile
from coot import map_peaks_py, map_peaks_near_point_from_list_py, density_at_point, add_molecule_py
add_molecule = add_molecule_py

RADIUS = 10                             #the maximum distance for the next phosphate
SUGAR_ROTATION_INTERVAL = radians(5)    #how much to rotate the sugar center by when trying to place it in density
DENSITY_CHECK_POINTS = 10               #when trying to place the sugar center into density, how many points should we check along each phosphate-sugar baton
GET_PEAKS_STARTING_SIGMA = 7            #the map sigma level to start the peak search at when searching for potential phosphate locations
GET_PEAKS_ENDING_SIGMA = 1              #the map sigma level to end the peak search at
GET_PEAKS_SIGMA_STEP = 0.25             #the interval for map sigma level for the peak search
#GET_PEAKS_MAX_PEAKS = 15        
DENSITY_WEIGHT = 15                     #when scoring peaks for the next phosphate location, the weight for the map density term
PHOS_DIST_WEIGHT = 5                    #when scoring peaks for the next phosphate location, the weight for the phosphate-phosphate distance term
                                        #   Note that when scoring peaks for the next phosphate location, the P-P-P angle and sugar-P-sugar angle terms both have weight 1
SECOND_PHOS_DENSITY_WEIGHT = 10         #when scoring peaks for the second phosphate in a chain, the weight to place on the map density term
                                        #   Here, the phosphate distance term has weight 1 (and is the only other term considered)


class NextPhos:
    """A class to determine the coordinates of the next phosphate atom in electron density"""
    
    
    __peakList = {} #a cached list of map peaks
        #by declaring this variable outside of the __init__ function, it's equivalent to a static variable in C++
        #Note, this caching assumes that imol numbers are never reused (which was confirmed to be true by Paul)
    
    def __init__(self, phosDistFilename, phosAngleFilename, sugarPhosSugarAngleFilename, basesFilename, pseudoChiFilename):
        """Initialize the NextPhos object with distance and angle data
        
        ARGUMENTS:
           phosDistFilename             - the name of a file containing smoothed data about inter-phosphate distances
           phosAngleFilename            - the name of a file containing smoothed data about inter-phosphate angles
           sugarPhosSugarAngleFilename  - the name of a file containing smoothed data about sugar-phosphate-sugar angles
           basesFilename                - the name of a PDB file containing structures of all four bases
           pseudoChiFilename            - the name of a file containing smoothed data about the pseudo-chi torsion (i.e. about base rotation relative to a simple backbone trace)
        RETURNS:
            an initialized NextPhos object
        """
        
        self.__phosDistInterp = LinearInterp(phosDistFilename)
        self.__phosAngleInterp = LinearInterp(phosAngleFilename)
        self.__sugarPhosSugarAngleInterp = LinearInterp(sugarPhosSugarAngleFilename)
        self.__baseStrucs = self.__readBases(basesFilename)
        self.__pseudoChiInterp = LinearInterp(pseudoChiFilename)
    
    
    def __readBases(self, filename):
        """Read nucleotide structures from the specified PDB file
        ARGUMENTS:
            filename - the name of a PDB file containing structures of all four bases
        RETURNS:
            struc    - a dictionary of file strucutres in the form base type => atom name => coordinates
        """
        
        input = open(filename, 'r')
        struc = dict(A = dict(), C = dict(), G = dict(), U = dict())
        
        #read in the PDB file
        for curline in input.readlines():
            if curline[0:6] != "ATOM  ":
                continue
            atomName = curline[12:16].strip()
            resName  = curline[17:20].strip()
            x = float(curline[30:38])
            y = float(curline[38:46])
            z = float(curline[46:54])
            
            #convert the atom name to PDB3 format
            #we'll need to convert it back to PDB2 later, since that's what Coot uses, but internally, RCrane uses PDB3
            atomName = atomName.replace("*", "'")
            
            struc[resName][atomName] = [x,y,z]
        input.close()
            
        #translate each base so that the C1' atom is at [0,0,0]
        for curRes in struc.keys():
            c1loc = struc[curRes]["C1'"]
            for curAtom in struc[curRes].keys():
                struc[curRes][curAtom] = minus(struc[curRes][curAtom], c1loc)
        
        return struc
        
    
    def nextPhos(self, mapNum, curPhos, prevPhos, prevSugar, direction = 3):
        """Find the next phosphate in a chain using phosphate distance and angle data
        
        ARGUMENTS:
            mapNum    - the molecule number of the Coot map to use
            curPhos   - the coordinates of the current phosphate
            prevPhos  - the coordinates of the previous phosphate
            prevSugar - the coordinates of the previous C1' atom
        OPTIONAL ARGUMENTS:
            direction - which direction to trace the chain: 3 implies 5'->3'
                                                            5 implies 3'->5'
                        defaults to 3 (5'->3')
        RETURNS:
            peakList  - a list of potential phosphate peaks
            sugarList - a list of potential C1' locations for each phosphate peak
        """
        
        peaks = self.getPeaks(mapNum, curPhos)
        
        #calculate the score for each peak
        peakScores = []
        for curPeak in peaks:
            #exclude any peaks that are within an Angstom of the current phosphate position
            #this will exclude the peak corresponding to the current phosphate position
            if dist(curPhos, curPeak) >= 1:
                
                #do a sugar search for the current phosphate
                if direction == 3:
                    sugarLocs = self.findSugar(mapNum, curPhos, curPeak)
                else:
                    sugarLocs = self.findSugar(mapNum, curPeak, curPhos)
                
                densityScore = sugarLocs[0][3] + curPeak[3]
                
                score = ((PHOS_DIST_WEIGHT * ln(self.__phosDistInterp.interp(dist(curPhos, curPeak))))
                         + ln(self.__phosAngleInterp.interp(angle(prevPhos, curPhos, curPeak)))
                         + ln(self.__sugarPhosSugarAngleInterp.interp(angle(prevSugar, curPhos, sugarLocs[0][0:3])))
                         + (DENSITY_WEIGHT * ln(densityScore))
                        )
                peakScores.append((curPeak, score, sugarLocs))
                
                #FOR DEBUGGING OUTPUT
                #peakScores.append((curPeak, score, sugarLocs, ln(self.__phosDistInterp.interp(dist(curPhos, curPeak))), ln(self.__phosAngleInterp.interp(angle(prevPhos, curPhos, curPeak))), ln(self.__sugarPhosSugarAngleInterp.interp(angle(prevSugar, prevPhos, sugarLocs[0][0:3]))), sugarLocs[0][3], sugarLocs[0][4], sugarLocs[0][5], sugarLocs[0][6], curPeak[3], sugarLocs[0][7]))
                #peakScores.append((curPeak, score, sugarLocs, ln(self.__phosDistInterp.interp(dist(curPhos, curPeak))), ln(self.__phosAngleInterp.interp(angle(prevPhos, curPhos, curPeak))), ln(self.__sugarPhosSugarAngleInterp.interp(angle(prevSugar, curPhos, sugarLocs[0][0:3]))), ln(densityScore), DENSITY_WEIGHT * ln(densityScore), sugarLocs[0][3], curPeak[3]))
                #print "Score for peak %(curPeak)s is %(score)f" % vars()
                #print "\tDist: " + str(dist(curPhos, curPeak)) + "\t" + str(self.__phosDistInterp.interp(dist(curPhos, curPeak)))
                #print "\tAngle: " + str(angle(prevPhos, curPhos, curPeak)) + "\t" + str(self.__phosAngleInterp.interp(angle(prevPhos, curPhos, curPeak)))
                
        
        #sort the peaks according to score
        peakScores.sort(key = lambda x: x[1], reverse = 1)
        
        #return only the coordinates, not the scores
        peakList  = [x[0] for x in peakScores]
        sugarList = [x[2] for x in peakScores]
        
        #FOR DEBUGGING OUTPUT
        #print "Current scores:"
        #print "overall\t\tphosDist\tphosAngle\tsPsAngle\tdensity score\tw. den score\tsugar score\tphos score"
        #for x in peakScores:
        #    #print "SCORE: %f\tSUGAR SCORE: %f" % (x[1], x[2][0][3])
        #    print "\t".join(map(str, (x[1],) + x[3:]))
        
        return (peakList, sugarList)
    
    
    def secondPhos(self, mapNum, curPhos, direction = 3):
        """find the second phosphate in a chain using phosphate distance data
        
        ARGUMENTS:
            mapNum    - the molecule number of the Coot map to use
            curPhos   - the coordinates of the first phosphate
        OPTIONAL ARGUMENTS:
            direction - which direction to trace the chain: 3 implies 5'->3'
                                                            5 implies 3'->5'
                        defaults to 3 (5'->3')
        RETURNS:
            peakList  - a list of potential phosphate peaks
            sugarList - a list of potential C1' locations for each phosphate peak
        """
        
        peaks = self.getPeaks(mapNum, curPhos)
        
        #calculate the score for each peak
        peakScores = []
        for curPeak in peaks:
            #exclude any peaks that are within an Angstom of the current phosphate position
            #this will exclude the peak corresponding to the current phosphate position
            if dist(curPhos, curPeak) >= 1:
                
                #do a sugar search for the current phosphate
                if direction == 3:
                    sugarLocs = self.findSugar(mapNum, curPhos, curPeak)
                else: #direction == 5:
                    sugarLocs = self.findSugar(mapNum, curPeak, curPhos)
                densityScore = sugarLocs[0][3] + curPeak[3]
                
                score = ln(self.__phosDistInterp.interp(dist(curPhos, curPeak))) + (SECOND_PHOS_DENSITY_WEIGHT * ln(densityScore))
                peakScores.append((curPeak, score, sugarLocs))
        
        #sort the peaks according to score
        peakScores.sort(key = lambda x: x[1], reverse = 1)
        
        #return only the coordinates, not the scores
        peakList  = [x[0] for x in peakScores]
        sugarList = [x[2] for x in peakScores]
        
        return (peakList, sugarList)
    
    
    def firstPhos(self, mapNum, curPos):
        """find the first phosphate in a chain
        
        ARGUMENTS:
            mapNum    - the molecule number of the Coot map to use
            curPos    - where to start the search (i.e. the current screen center coordinates)
        RETURNS:
            peaks     - a list of potential phosphate peaks
        """
        
        peaks = self.getPeaks(mapNum, curPos)
        
        #sort the peaks based on how close they are to curPos
        peaks.sort(key=lambda x: dist(curPos, x))
        
        return peaks
    
    
    def findSugar(self, mapNum, phos5, phos3):
        """find potential C1' locations between the given 5' and 3' phosphate coordinates
        
        ARGUMENTS:
            mapNum      - the molecule number of the Coot map to use
            phos5       - the coordinates of the 5' phosphate
            phos3       - the coordinates of the 3' phosphate
        RETURNS:
            sugarMaxima - a list of potential C1' locations, each listed as [x, y, z, score]
        """
        
        #calculate the distance between the two phosphatse
        phosPhosDist = dist(phos5, phos3)
        
        #calculate a potential spot for the sugar based on the phosphate-phosphate distance
        #projDist is how far along the 3'P-5'P vector the sugar center should be (measured from the 3'P)
        #perpDist is how far off from the 3'P-5'P vector the sugar center should be
        
        #these functions are for the sugar center, which is what we try to find here
        #since it will be more in the center of the density blob
        perpDist = -0.185842*phosPhosDist**2 + 1.62296*phosPhosDist - 0.124146
        projDist = 0.440092*phosPhosDist + 0.909732
        
        #if we wanted to find the C1' instead of the sugar center, we'd use these functions
        #however, finding the C1' directly causes our density scores to be less accurate
        #so we instead use the functions above to find the sugar center and later adjust our
        #coordinates to get the C1' location
        #perpDist = -0.124615*phosPhosDist**2 + 0.955624*phosPhosDist + 2.772573
        #projDist = 0.466938*phosPhosDist + 0.649833
        
        
        #calculate the normal to the plane defined by 3'P, 5'P, and a dummy point
        normal = crossProd([10,0,0], minus(phos3, phos5))
        
        #make sure the magnitude of the normal is not zero (or almost zero)
        #if it is zero, that means that our dummy point was co-linear with the 3'P-5'P vector
        #and we haven't calculated a normal
        #if the magnitude is almost zero, then the dummy point was almost co-linear and we have to worry about rounding error
        #in either of those cases, just use a different dummy point
        #they should both be incredibly rare cases, but it doesn't hurt to be safe
        if magnitude(normal) < 0.001:
            #print "Recalculating normal"
            normal = crossProd([0,10,0], minus(phos5, phos3))
        
        
        #scale the normal to the length of perpDist
        perpVector = scalarProd(normal, perpDist/magnitude(normal))
        
        #calculate the 3'P-5'P vector and scale it to the length of projDist
        projVector = minus(phos3, phos5)
        projVector = scalarProd(projVector, projDist/magnitude(projVector))
        
        #calculate a possible sugar location
        sugarLoc = plus(phos5, projVector)
        sugarLoc = plus(sugarLoc, perpVector)
        
        #rotate the potential sugar location around the 3'P-5'P vector to generate a list of potential sugar locations
        sugarRotationPoints = self.__rotateSugarCenter(phos5, phos3, sugarLoc)
        
        #test each potential sugar locations to find the one with the best electron density
        for curSugarLocFull in sugarRotationPoints:
            curSugarLoc = curSugarLocFull[0:3] #the rotation angle is stored as curSugarLocFull[4], so we trim that off for curSugarLoc
            curDensityTotal = 0
            #densityList = []   #if desired, this could be used to generate additional statistics on the density (such as the median or quartiles)
            
            #check density along the 5'P-sugar vector
            phosSugarVector = minus(curSugarLoc, phos5)
            phosSugarVector = scalarProd(phosSugarVector, 1.0/(DENSITY_CHECK_POINTS+1))
            for i in range(1, DENSITY_CHECK_POINTS+1):
                (x, y, z) = plus(phos5, scalarProd(i, phosSugarVector))
                curPointDensity = density_at_point(mapNum, x, y, z)
                curDensityTotal += curPointDensity
                #densityList.append(curPointDensity)
                
            
            #check at the sugar center
            (x, y, z) = curSugarLoc
            curPointDensity = density_at_point(mapNum, x, y, z)
            curDensityTotal += curPointDensity
            #densityList.append(curPointDensity)
            
            #check along the sugar-3'P vector
            sugarPhosVector = minus(phos3, curSugarLoc)
            sugarPhosVector = scalarProd(sugarPhosVector, 1.0/(DENSITY_CHECK_POINTS+1))
            for i in range(1, DENSITY_CHECK_POINTS+1):
                (x, y, z) = plus(curSugarLoc, scalarProd(i, sugarPhosVector))
                curPointDensity = density_at_point(mapNum, x, y, z)
                curDensityTotal += curPointDensity
                #densityList.append(curPointDensity)
            
            curSugarLocFull.append(curDensityTotal)
            #curSugarLocFull.extend([curDensityTotal, median(densityList), lowerQuartile(densityList), min(densityList)])#, pointList])
        
        #find all the local maxima
        sugarMaxima = []
        curPeakHeight = sugarRotationPoints[-1][4]
        nextPeakHeight = sugarRotationPoints[0][4]
        sugarRotationPoints.append(sugarRotationPoints[0]) #copy the first point to the end so that we can properly check the last point
        for i in range(0, len(sugarRotationPoints)-1):
            prevPeakHeight = curPeakHeight
            curPeakHeight  = nextPeakHeight
            nextPeakHeight = sugarRotationPoints[i+1][4]
            if prevPeakHeight < curPeakHeight and curPeakHeight >=  nextPeakHeight:
                sugarMaxima.append(sugarRotationPoints[i])
        
        #sort the local maxima by their density score
        sugarMaxima.sort(key = lambda x: x[4], reverse = True)
        
        #adjust all the sugar center coordinates so that they represent the corresponding C1' coordinates
        for i in range(0, len(sugarMaxima)):
            curSugar = sugarMaxima[i][0:3]
            #rotate a vector 148 degrees from the phosphate bisector
            phosAngle = angle(phos5, curSugar, phos3)
            phos5vector = minus(phos5, curSugar)
            axis = crossProd(minus(phos3, curSugar), phos5vector)
            axis = scalarProd(axis, 1/magnitude(axis))
            c1vec = rotate(phos5vector, axis, 148.539123-phosAngle/2)
            
            #scale the vector to the appropriate length
            c1vec = scalarProd(c1vec, 1.235367/magnitude(c1vec))
            
            #rotate the vector about the phosphate bisector
            phosBisectorAxis = rotate(phos5vector, axis, -phosAngle/2)
            phosBisectorAxis = scalarProd(phosBisectorAxis, 1/magnitude(phosBisectorAxis)) 
            c1vec = rotate(c1vec, phosBisectorAxis, -71.409162)
            
            sugarMaxima[i][0:3] = plus(c1vec, curSugar)
        
        return sugarMaxima
        
    
    def __rotateSugarCenter (self, phos5, phos3, sugarCenter):
        """rotate the sugar center by 360 degrees in ROTATE_SUGAR_INTERVAL increments
        
        ARGUMENTS:
            phos5         - the coordinates of the 5' phosphate
            phos3         - the coordinates of the 3' phosphate
            sugarCenter   - the coordinates of the sugar center to be rotated
        RETURNS:
            rotatedPoints - a list of the rotated points, each listed as [x, y, z, rotation angle]
        """
        
        #calculate a unit vector along the rotation axis
        axis = minus(phos3, phos5)
        axis = scalarProd(axis, 1/magnitude(axis))
        
        #perform the rotation
        (u, v, w) = axis
        (x, y, z) = minus(sugarCenter, phos5)
        
        #make sure that the original location appears on the list with a rotation value of 0
        sugarCenterRot = sugarCenter + [0]
        rotatedPoints = [sugarCenterRot]
        
        curAngle = SUGAR_ROTATION_INTERVAL
        while curAngle < (2*pi - 0.5*SUGAR_ROTATION_INTERVAL):
            cosTheta = cos(curAngle)
            sinTheta = sin(curAngle)
            
            a = u*x + v*y + w*z;
            newX = a*u + (x-a*u)*cosTheta + (v*z-w*y)*sinTheta + phos5[0];
            newY = a*v + (y-a*v)*cosTheta + (w*x-u*z)*sinTheta + phos5[1];
            newZ = a*w + (z-a*w)*cosTheta + (u*y-v*x)*sinTheta + phos5[2];
            
            rotatedPoints.append([newX, newY, newZ, curAngle])
            curAngle += SUGAR_ROTATION_INTERVAL
            
        return rotatedPoints
    
    
    def findBase(self, mapNum, sugar, phos5, phos3, baseType, direction = 3):
        """Rotate the sugar center by 360 degrees in ROTATE_SUGAR_INTERVAL increments
        
        ARGUMENTS:
            mapNum   - the molecule number of the Coot map to use
            sugar    - the coordinates of the C1' atom
            phos5    - the coordinates of the 5' phosphate
            phos3    - the coordinates of the 3' phosphate
            baseType - the base type (A, C, G, or U)
        OPTIONAL ARGUMENTS:
            direction - which direction are we tracing the chain
                        if it is 5 (i.e. 3'->5'), then phos5 and phos3 will be flipped
                        all other values will be ignored
                        defaults to 3 (i.e. 5'->3')
        RETURNS:
            baseObj  - a list of [baseType, baseCoordinates]
        """
        
        if direction == 5:
            (phos5, phos3) = (phos3, phos5)
        
        #calculate the bisector of the phos-sugar-phos angle
        #first, calculate a normal to the phos-sugar-phos plane
        sugarPhos5Vec = minus(phos5, sugar)
        sugarPhos3Vec = minus(phos3, sugar)
        normal = crossProd(sugarPhos5Vec, sugarPhos3Vec)
        normal = scalarProd(normal, 1.0/magnitude(normal))
        
        phosSugarPhosAngle = angle(phos5, sugar, phos3)
        
        bisector = rotate(sugarPhos5Vec, normal, phosSugarPhosAngle/2.0)
        
        
        #flip the bisector around (so it points away from the phosphates) and scale its length to 5 A
        startingBasePos = scalarProd(bisector, -1/magnitude(bisector))
        
        #rotate the base baton by 10 degree increments about half of a sphere
        rotations = [startingBasePos] #a list of coordinates for all of the rotations
        for curTheta in range(-90, -1, 10) + range(10, 91, 10):
            curRotation = rotate(startingBasePos, normal, curTheta)
            rotations.append(curRotation) #here's where the phi=0 rotation is accounted for
            
            for curPhi in range(-90, -1, 10) + range(10, 91, 10):
                rotations.append(rotate(curRotation, startingBasePos, curPhi))
                
        #test electron density along all base batons
        for curBaton in rotations:
            curDensityTotal = 0
            densityList = []
            for i in range(1, 9):
                (x, y, z) = plus(sugar, scalarProd(i/2.0, curBaton))
                curPointDensity = density_at_point(mapNum, x, y, z)
                curDensityTotal += curPointDensity
                densityList.append(curPointDensity)
            curBaton.append(curDensityTotal)        #the sum of the density (equivalent to the mean for ordering purposes)
            curBaton.append(median(densityList))    #the median of the density
            curBaton.append(min(densityList))       #the minimum of the density
        
        #find the baton with the max density (as measured using the median)
        #Note that we ignore the sum and minimum of the density.  Those calculations could be commented out,
        #   but they may be useful at some point in the future.  When we look at higher resolutions maybe?
        #   Besides, they're fast calculations.)
        baseDir = max(rotations, key = lambda x: x[4])
        
        #rotate the stock base+sugar structure to align with the base baton
        rotationAngle = angle(self.__baseStrucs["C"]["C4"], [0,0,0], baseDir)
        axis = crossProd(self.__baseStrucs["C"]["C4"], baseDir[0:3])
        
        orientedBase = rotateAtoms(self.__baseStrucs["C"], axis, rotationAngle)
        
        #rotate the base about chi to find the best fit to density
        bestFitBase = None
        maxDensity = -999999
        for curAngle in range(0,360,5):
            rotatedBase = rotateAtoms(orientedBase, orientedBase["C4"], curAngle, sugar)
            curDensity = 0
            for curAtom in ["N1", "C2", "N3", "C4", "C5", "C6"]:
                curDensity += density_at_point(mapNum, rotatedBase[curAtom][0], rotatedBase[curAtom][1], rotatedBase[curAtom][2])
            
            #this is "pseudoChi" because it uses the 5' phosphate in place of the O4' atom
            pseudoChi = torsion(phos5, sugar, rotatedBase["N1"], rotatedBase["N3"])
            curDensity *= self.__pseudoChiInterp.interp(pseudoChi)
            
            if curDensity > maxDensity:
                maxDensity = curDensity
                bestFitBase = rotatedBase
        
        baseObj = ["C", bestFitBase]
        
        #mutate the base to the appropriate type
        if baseType != "C":
            baseObj = self.mutateBase(baseObj, baseType)
        
        return baseObj
    
    
    def mutateBase(self, curBase, newBaseType):
        """Change the base type.
        
        ARGUMENTS:
            curBase - the current base object in the form [baseType, baseCoordinates]
            newBaseType - the base type to mutate to
        RETURNS:
            baseObj  - a list of [baseType, baseCoordinates]
        """
        
        (curBaseType, curBaseCoords) = curBase
        
        #calculate the vectors used to align the old and new bases
        #for pyrimidines, the vector is C1'-C4
        #for purines, the vector is from C1' to the center of the C4-C5 bond
        curAlignmentVector = None
        if curBaseType == "C" or curBaseType == "U":
            curAlignmentVector = minus(curBaseCoords["C4"], curBaseCoords["C1'"])
        else:
            curBaseCenter = plus(curBaseCoords["C4"], curBaseCoords["C5"])
            curBaseCenter = scalarProd(1.0/2.0, curBaseCenter)
            curAlignmentVector = minus(curBaseCenter, curBaseCoords["C1'"])
        
        #calculate the alignment vector for the new base
        newBaseCoords = self.__baseStrucs[newBaseType]
        newAlignmentVector = None
        if newBaseType == "C" or newBaseType == "U":
            newAlignmentVector = newBaseCoords["C4"]
        else:
            newAlignmentVector = plus(newBaseCoords["C4"], newBaseCoords["C5"])
            newAlignmentVector = scalarProd(1.0/2.0, newAlignmentVector)
        
        #calculate the angle between the alignment vectors
        rotationAngle = -angle(curAlignmentVector, [0,0,0], newAlignmentVector)
        axis = crossProd(curAlignmentVector, newAlignmentVector)
        
        #rotate the new base coordinates
        newBaseCoords = rotateAtoms(newBaseCoords, axis, rotationAngle)
        
        #calculate the normals of the base planes
        curNormal = None
        if curBaseType == "C" or curBaseType == "U":
            curNormal = crossProd(minus(curBaseCoords["N3"], curBaseCoords["N1"]), minus(curBaseCoords["C6"], curBaseCoords["N1"]))
        else:
            curNormal = crossProd(minus(curBaseCoords["N3"], curBaseCoords["N9"]), minus(curBaseCoords["N7"], curBaseCoords["N9"]))
        
        newNormal = None
        if newBaseType == "C" or newBaseType == "U":
            newNormal = crossProd(minus(newBaseCoords["N3"], newBaseCoords["N1"]), minus(newBaseCoords["C6"], newBaseCoords["N1"]))
        else:
            newNormal = crossProd(minus(newBaseCoords["N3"], newBaseCoords["N9"]), minus(newBaseCoords["N7"], newBaseCoords["N9"]))
        
        #calculate the angle between the normals
        normalAngle = -angle(curNormal, [0,0,0], newNormal);
        normalAxis = crossProd(curNormal, newNormal)
        
        #rotate the new base coordinates so that it falls in the same plane as the current base
        #and translate the base to the appropriate location
        newBaseCoords = rotateAtoms(newBaseCoords, normalAxis, normalAngle, curBaseCoords["C1'"])
        
        return [newBaseType, newBaseCoords]
    
    
    def flipBase(self, curBase):
        """Flip the base anti/syn
        
        ARGUMENTS:
            curBase - the current base object in the form [baseType, baseCoordinates]
        RETURNS:
            the base object with rotated coordinates
        NOTE:
            This function does not necessarily simply rotate about chi.  When flipping a purine, it also adjusts
            the glycosidic bond position so that it the base a chance of staying in the density
        """
        
        (baseType, curBaseCoords) = curBase
        
        curC1coords = curBaseCoords["C1'"]
        
        newBaseCoords = {}
        #translate the base to the origin
        for atomName, curAtomCoords in curBaseCoords.items():
            newBaseCoords[atomName] = minus(curBaseCoords[atomName], curC1coords)
        
        #the rotation axis is the same as the alignment vector for mutating the base
        #for pyrimidines, the axis is C1'-C4
        #for purines, the axis is from C1' to the center of the C4-C5 bond
        axis = None
        if baseType == "C" or baseType == "U":
            axis = newBaseCoords["C4"]
        else:
            axis = plus(newBaseCoords["C4"], newBaseCoords["C5"])
            axis = scalarProd(1.0/2.0, axis)
        
        #rotate the base 180 degrees about the axis
        newBaseCoords = rotateAtoms(newBaseCoords, axis, 180, curC1coords)
        
        return [baseType, newBaseCoords]
    
    
    def getPeaks(self, mapNum, pos, radius=RADIUS):
        """Get a list of peaks in the map near a given position
        
        ARGUMENTS:
            mapNum      - the molecule number of the Coot map to use
            pos         - the coordinates to search near
        OPTIONAL ARGUMENTS:
            radius      - how far from pos should we look for peaks?
                          defaults to RADIUS (defined at the top of this file)
        RETURNS:
            closePeaks  - a list of peaks in the format [x, y, z, density]
        """
        
        (x, y, z) = pos
        
        #determine if we've already done a peak search for this map
        if not mapNum in self.__peakList:
            
            peakList = []
            peakHash = dict()
            
            #go through all sigma levels and do a peak search
            curSigma = GET_PEAKS_STARTING_SIGMA
            while curSigma >= GET_PEAKS_ENDING_SIGMA:
                
                #do a peak search at the current sigma level
                curPeaks = map_peaks_py(mapNum, curSigma)
                
                #go through each peak we found and make sure we hadn't found it already
                for peak in curPeaks:
                    #check peakHash for this peak
                    try:
                        peakHash[peak[0]][peak[1]][peak[2]]
                    except KeyError:
                        #if it's not there, add it to peakList and peakHash
                        peakList.append(peak)
                        if not(peakHash.has_key(peak[0])):
                            peakHash[peak[0]] = dict()
                        if not(peakHash[peak[0]].has_key(peak[1])):
                            peakHash[peak[0]][peak[1]] = dict()
                        if not(peakHash[peak[0]][peak[1]].has_key(peak[2])):
                            peakHash[peak[0]][peak[1]][peak[2]] = True
                
                curSigma -= GET_PEAKS_SIGMA_STEP
            
            #once we've found all the peaks in this map, store them
            self.__peakList[mapNum] = peakList
        
        
        closePeaks = map_peaks_near_point_from_list_py(mapNum, self.__peakList[mapNum], x, y, z, radius)
        
        for curPeak in closePeaks:
            (x, y, z) = curPeak
            density = density_at_point(mapNum, x, y, z)
            curPeak.append(density)
        
        return closePeaks


def rotate(vector, axis, angle):
    """Rotate a vector.
    
    ARGUMENTS:
        vector - the vector to rotate
        axis   - the axis to rotate about (note that this MUST be a unit vector)
        angle  - the angle to rotate by
    RETURNS:
        the rotated vector
    NOTE:
        This function is very similar to the rotateAtoms() function defined in strucCalc.py, but rotate()
        rotates only a single vector, not a list of atoms; rotate() doesn't normalize the axis vector; and
        rotate() has no optional translate argument.  This module could easily be rewritten to use
        rotateAtoms() instead, but rotate() may be slightly faster, so I've kept it.
    """
    
    ##make sure the axis is a unit vector
    #axis = scalarProd(axis, 1.0/axisMag)
    
    (u, v, w) = axis
    (x, y, z) = vector
    
    angle = radians(angle)
    cosTheta = cos(angle)
    sinTheta = sin(angle)
    
    a = u*x + v*y + w*z;
    newX = a*u + (x-a*u)*cosTheta + (v*z-w*y)*sinTheta
    newY = a*v + (y-a*v)*cosTheta + (w*x-u*z)*sinTheta
    newZ = a*w + (z-a*w)*cosTheta + (u*y-v*x)*sinTheta
    
    return [newX, newY, newZ]
