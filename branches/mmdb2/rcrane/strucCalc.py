#!/usr/bin/env python
"""functions for working with 3D "vectors" (really just arrays of length three)"""

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

from math import sqrt, acos, degrees, atan2, radians, cos, sin

def plus (a, b):
    """Add two 3D vectors"""
    
    return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]

def minus (a, b):
    """Subtract two 3D vectors"""
    
    return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]

def magnitude(a):
    """Calculate the magnitude of a 3D vector"""
    
    return sqrt(a[0]**2 + a[1]**2 + a[2]**2)

def normalize(a, length = 1):
    """Set a vector to length 1 (or a custom length)
    
    ARGUMENTS:
        a - the vector to normalize
    OPTIONAL ARGUMENTS:
        length - the length to set the vector to (default 1)
    RETURNS:
        vector a normalized to length length
    """
    
    return scalarProd(length/magnitude(a), a)

def dist(a, b):
    """calculate the distance between two 3D vectors"""
    
    return magnitude(minus(a, b))

def scalarProd(scalar, vector):
    """Calculate the scalar product between a scalar and a 3D vector
    
    NOTE:
        this function can be called with the arguments reversed if desired
    """
    
    #flip the arguments if they were called in the wrong order
    if isinstance(scalar, list) or isinstance(scalar, tuple):
        (scalar, vector) = (vector, scalar)
    
    return [scalar * vector[0], scalar * vector[1], scalar * vector[2]]

def dotProd(a, b):
    """Calculate the dot product of two 3D vectors"""
    
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def crossProd(a,b):
    """Calculate the cross product of two vectors"""
    
    (ax, ay, az) = a
    (bx, by, bz) = b
    
    return [ ay*bz - az*by,
             az*bx - ax*bz,
             ax*by - ay*bx ]


def angle(a, b, c):
    """Calculate the angle between three 3D vectors"""
    
    #calculate the two vectors
    vec1 = minus(a,b)
    vec2 = minus(c,b)
    
    #normalize the two vectors
    vec1 = scalarProd(1/magnitude(vec1), vec1)
    vec2 = scalarProd(1/magnitude(vec2), vec2)
    
    #make sure that the dot product of the vectors is between -1 and 1
    #since if it's slightly over 1 (or less than -1) due to rounding error, the acos will fail
    dotprod = dotProd(vec1, vec2)
    if dotprod < -1:
        dotprod = -1
    elif dotprod > 1:
        dotprod = 1
        
    return degrees(acos(dotprod))

def torsion(atom1, atom2, atom3, atom4):
    """Calculate to torsion between four atoms"""
    
    vector1 = minus(atom2, atom1)
    vector2 = minus(atom3, atom2)
    vector3 = minus(atom4, atom3)
    
    #torsion angle = atan2 (|b2|b1 . (b2 x b3), (b1 x b2) . (b2 x b3)
    #source: http://en.wikipedia.org/wiki/Dihedral_angle
    
    y1 = scalarProd( magnitude(vector2), vector1)
    y2 = crossProd(vector2, vector3)
    y = dotProd(y1, y2)
    
    x1 = crossProd(vector1, vector2)
    x2 = y2
    x = dotProd(x1, x2)
    
    torsion = atan2(y, x)
    torsion = degrees(torsion)
    if torsion < 0:
        torsion += 360
    
    return torsion

def distToLine(pointP, pointC, pointN):
    """Calculate the distance from a point to a line, where the line is defined by two points on the line
    
    ARGUMENTS:
        pointP - the point not on the line
        pointC - point 1 that defines the line
        pointN - point 2 that defines the line
    RETURNS:
        the distance from pointP to a line through pointC and pointN
    NOTE:
        the distance from point P to line CN is calculated as:
                        |(C-N) x (N-P)|
            distance = -----------------
                             |C-N|
    """
    
    cMinusN = minus(pointC, pointN)
    nMinusP = minus(pointN, pointP)
    
    crossProduct = crossProd(cMinusN, nMinusP)
    
    return magnitude(crossProduct)/magnitude(cMinusN)
    

def rotateAtoms(atoms, axis, angle, translate = None):
    """Rotate a group of atoms about a specified axis
    
    ARGUMENTS:
        atoms        - the list of atoms to rotate
        axis         - the vector to rotate the atoms about
        angle        - the angle to rotate the atoms by
    OPTIONAL ARGUMENTS:
        translate    - if given, all atoms are tranlated by this vector *after* the rotation is performed
    RETURNS:
        rotatedAtoms - the rotated (and possibly translated) atoms
    """
    
    #make sure the axis is a unit vector
    axis = scalarProd(axis, 1.0/magnitude(axis))
    
    (u, v, w) = axis
    
    angle = radians(angle)
    cosTheta = cos(angle)
    sinTheta = sin(angle)
    
    rotatedAtoms = dict()
    for curAtom in atoms:
        (x, y, z) = atoms[curAtom]
        a = u*x + v*y + w*z;
        newX = a*u + (x-a*u)*cosTheta + (v*z-w*y)*sinTheta
        newY = a*v + (y-a*v)*cosTheta + (w*x-u*z)*sinTheta
        newZ = a*w + (z-a*w)*cosTheta + (u*y-v*x)*sinTheta
        rotatedAtoms[curAtom] = [newX, newY, newZ]
    
    if translate is not None:
        for curAtom in rotatedAtoms:
            rotatedAtoms[curAtom] = plus(translate, rotatedAtoms[curAtom])
    
    return rotatedAtoms
