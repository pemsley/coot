/*
 *  RotatedTranslatedScaledEntity.mm
 *  Aesop
 *
 *  Created by Martin Noble on 14/06/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
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
