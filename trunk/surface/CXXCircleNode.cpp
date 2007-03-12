/*
 *  CXXCircleNode.cpp
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#if defined _MSC_VER
#include <string>
#endif

#include "CXXCircleNode.h"
#include <CXXCoord.h>
#include <CXXCircle.h>
#include <CXXNewHood.h>


CXXCircleNode::CXXCircleNode():
theParent(0),
theAngle(0),
theFlag(-1),
theAtomI(0),
theAtomJ(0),
theAtomK(0),
thisIsDeleted(1)
{
}

CXXCircleNode::CXXCircleNode ( const CXXCircle *aParent, PCAtom atomK, const CXXCoord &crd, int aFlag) :
theParent(aParent),
theCoord(crd),
theFlag(aFlag),
theAtomJ(aParent->getAtomJ()),
theAtomK(atomK), 
thisIsDeleted(1)
{
	if (aParent->getParent()) theAtomI = aParent->getParent()->getAtomI();
	unitRadius = crd - aParent->getCentreOfCircle();
	unitRadius /= aParent->getRadiusOfCircle();
}

int CXXCircleNode::setReference(const CXXCoord &referenceVector){
	theAngle = theParent->getNormal().angleBetween(unitRadius, referenceVector);
	return 0; 
};
