#!/usr/bin/env python
"""Functions for building non-bridging phosphoryl oxygens"""

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

from math import radians, cos, sin

from strucCalc import minus, plus, dotProd, crossProd, scalarProd, magnitude, normalize

#constants
PHOSBONDLENGTH = 1.485
PHOSBONDANGLE  = 119.6
INITPHOSANGLE  = 108     #the angle between the O5' and the phosphoryl oxygens in the first phosphate of a segment

def buildPhosOxy(curResAtoms, prevResAtoms):
    """Calculate non-bridging phosphoryl oxygen coordinates using the coordinates of the current and previous nucleotides
    
    ARGUMENTS:
        curResAtoms - a dictionary of the current nucleotide coordinates (i.e. the nucleotide to build the phosphoryl oxygens on) in the format atomName: [x, y, z]
        prevResAtoms - a dictionary of the previous nucleotide coordinates in the format atomName: [x, y, z]
    RETURNS:
        phosOxyCoords - a dictionary of the phosphoryl oxygen coordinates in the format atomName: [x, y, z]
    """
    
    P  = curResAtoms ["P"]
    O5 = curResAtoms ["O5'"]
    O3 = prevResAtoms["O3'"]
    
    #calculate a line from O5' to O3'
    norm = minus(O5, O3)
    
    #calculate the intersection of a plane (with normal $norm and point $P) and a line (from O5' to O3')
    #using formula from http://local.wasp.uwa.edu.au/~pbourke/geometry/planeline/
    #       norm dot (P - O3)
    # i = ------------------------
    #       norm dot (O5 - O3)
    #intersecting point = O3 + i(O5 - O3)
    i = dotProd(norm, minus(P, O3)) / dotProd(norm, minus(O5, O3))
    interPoint = plus(O3, scalarProd(i, minus(O5, O3)))
    
    #move $interPoint so that the distance from $P to $interPoint is 1.485 (the length of the P-OP1 bond)
    #we also reflect the point about P
    PIline = minus(P, interPoint) #here's is where the reflection occurs, because we do $P-$interPoint instead of $interPoint-$P
    #scaledPoint = scalarProd(1/magnitude(PIline) * PHOSBONDLENGTH, PIline)
    scaledPoint = normalize(PIline, PHOSBONDLENGTH)
    #to get the new point location, we would do P + scaledPoint
    #but we need to rotate the point first before translating it back
    
    #rotate this new point by 59.8 and -59.8 degrees to determine the phosphoryl oxygen locations
    #we rotate about the axis defined by $norm
    angle = radians(PHOSBONDANGLE / 2)
    (x, y, z) = scaledPoint
    
    #unitnorm = scalarProd( 1/magnitude(norm), norm)
    unitnorm = normalize(norm)
    (u, v, w) = unitnorm
    
    a = u*x + v*y + w*z
    
    phosOxyCoords = {}
    for (atomName, theta) in (("OP1", angle), ("OP2", -angle)):
        cosTheta = cos(theta)
        sinTheta = sin(theta)
        
        #perform the rotation, and then add $P to the coordinates
        newX = a*u + (x-a*u)*cosTheta + (v*z-w*y)*sinTheta + P[0]
        newY = a*v + (y-a*v)*cosTheta + (w*x-u*z)*sinTheta + P[1]
        newZ = a*w + (z-a*w)*cosTheta + (u*y-v*x)*sinTheta + P[2]
        
        phosOxyCoords[atomName] = [newX, newY, newZ]
    
    return phosOxyCoords



def buildInitOrTerminalPhosOxy(curResAtoms, prevResAtoms = None):
    """build phosphoryl oxygens using only the O3' or O5'atom (intended for the first or last nucleotide of a chain/segment)
    
    ARGUMENTS:
        curResAtoms - a dictionary of the current nucleotide coordinates (i.e. the nucleotide to build the phosphoryl oxygens on) in the format atomName: [x, y, z]
                      this dictionary must contain at least the phosphate
    OPTIONAL ARGUMENTS:
        prevResAtoms - a dictionary of the previous nucleotide coordinates in the format atomName: [x, y, z]
                       if provided, the curResAtoms["P"] and prevResAtoms["O3'"] will be used to place the phosphoryl oxygens
                       if not provided, the curResAtoms["P"] and curResAtoms["O5'"]
    RETURNS:
        phosOxyCoords - a dictionary of the phosphoryl oxygen coordinates in the format atomName: [x, y, z]
    """
    
    P  = curResAtoms ["P"]
    if prevResAtoms is not None:
        O = prevResAtoms["O3'"]
        C = prevResAtoms["C3'"]
    else:
        O = curResAtoms["O5'"]
        C = curResAtoms["C5'"]
        
    #place atom along the P-O5' bond that is the appropriate distance from P
    phosOxy = minus(O, P)
    phosOxy = normalize(phosOxy, PHOSBONDLENGTH)
    
    #define plane with C5'-O5'-P
    norm = crossProd( minus(C, O), minus(P, O))
    norm = normalize(norm)
    
    #rotate dummy atom in plane about P by INITPHOSANGLE (which is the appropriate O5'-OP1 angle)
    (x, y, z) = phosOxy
    (u, v, w) = norm
    theta = -radians(INITPHOSANGLE) #use the negative angle so that we rotate away from the C5'
    
    cosTheta = cos(theta)
    sinTheta = sin(theta)
    
    #perform the rotation
    a = u*x + v*y + w*z
    phosOxy[0] = a*u + (x-a*u)*cosTheta + (v*z-w*y)*sinTheta
    phosOxy[1] = a*v + (y-a*v)*cosTheta + (w*x-u*z)*sinTheta
    phosOxy[2] = a*w + (z-a*w)*cosTheta + (u*y-v*x)*sinTheta
    
    #rotate dummy atom about O5'-P axis by 59.8 and -59.8 degrees
    norm = minus(O, P)
    norm = normalize(norm)
    
    (x, y, z) = phosOxy
    (u, v, w) = norm
    angle = radians(PHOSBONDANGLE / 2)
    
    phosOxyCoords = {}
    for (atomName, theta) in (("OP1", angle), ("OP2", -angle)):
        cosTheta = cos(theta)
        sinTheta = sin(theta)
        
        #perform the rotation, and then add $P to the coordinates
        newX = a*u + (x-a*u)*cosTheta + (v*z-w*y)*sinTheta + P[0]
        newY = a*v + (y-a*v)*cosTheta + (w*x-u*z)*sinTheta + P[1]
        newZ = a*w + (z-a*w)*cosTheta + (u*y-v*x)*sinTheta + P[2]
        
        phosOxyCoords[atomName] = [newX, newY, newZ]
    
    return phosOxyCoords
    
    