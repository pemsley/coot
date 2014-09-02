/*
 *  CXXSphereNode.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXSphereNode.h"
#include "CXXCoord.h"
#include "CXXTriangle.h"
#include <list>

CXXSphereNode::CXXSphereNode() : 
theVertex(CXXCoord(0.,0.,0.)),
theIntersector(0),
shouldBeDrawn(1),
theAtom(0)
{
//	references.resize(0);
}

CXXSphereNode::CXXSphereNode(const CXXCoord &aCoord) : 
theVertex(aCoord), 
theIntersector(0), 
shouldBeDrawn(1) ,
theAtom(0)
{
//	references.resize(0);
}

const CXXCoord &CXXSphereNode::vertex() const{
	return theVertex;
}

int CXXSphereNode::doDraw() const{
	return shouldBeDrawn;
}

int CXXSphereNode::setDoDraw(const int yesNo){
	shouldBeDrawn = yesNo;
	return shouldBeDrawn;
}

int CXXSphereNode::setVertex(const CXXCoord &crd){
	theVertex = crd;
	return 0;
}

void CXXSphereNode::setIntersector (const CXXCircle *aCircle){
	theIntersector = aCircle;
}

const CXXCircle *CXXSphereNode::getIntersector() const{
	return theIntersector;
}

