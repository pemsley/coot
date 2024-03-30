/*
 * MoleculesToTriangles/CXXClasses/RotatedTranslatedScaledEntity.cpp
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "RotatedTranslatedScaledEntity.h"
#include "Quaternion.hpp"

void RotatedTranslatedScaledEntity::rotateBy(const FCXXCoord  &vector)
{
    Quaternion<CXXCoord_ftype>rotateBy;
    CXXCoord_ftype rotateByArray[4];
    for (auto i=0; i<4; i++) rotateByArray[i] = vector.xyzr[i];
    rotateBy.RotationToQuaternion(rotateByArray);
    
    Quaternion<CXXCoord_ftype>beforeRotation;
    CXXCoord_ftype beforeRotationArray[4];
    for (auto i=0; i<4; i++) beforeRotationArray[i] = getRotation().xyzr[i];
    beforeRotation.RotationToQuaternion(beforeRotationArray);
    
    Quaternion<CXXCoord_ftype>afterRotation(beforeRotation^rotateBy);
    afterRotation.normalize();
    CXXCoord_ftype afterRotationArray[4];
    if (!afterRotation.IsIdentityRotation(1e-9, afterRotationArray)){
        afterRotation.QuaternionToRotation(afterRotationArray);
    }
    FCXXCoord  afterRotationCoord(afterRotationArray[0], afterRotationArray[1], afterRotationArray[2], afterRotationArray[3]);
    setRotation(afterRotationCoord);
}
