/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
/*
 *  CXXCircleNode.cpp
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  
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
