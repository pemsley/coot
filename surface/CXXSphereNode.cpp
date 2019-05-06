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

CXX_mot::CXXSphereNode::CXXSphereNode() : 
theVertex(CXXCoord(0.,0.,0.)),
theIntersector(0),
shouldBeDrawn(1),
theAtom(0)
{
//	references.resize(0);
}

CXX_mot::CXXSphereNode::CXXSphereNode(const CXXCoord &aCoord) : 
theVertex(aCoord), 
theIntersector(0), 
shouldBeDrawn(1) ,
theAtom(0)
{
//	references.resize(0);
}

const CXX_mot::CXXCoord &CXX_mot::CXXSphereNode::vertex() const{
	return theVertex;
}

int CXX_mot::CXXSphereNode::doDraw() const{
	return shouldBeDrawn;
}

int CXX_mot::CXXSphereNode::setDoDraw(const int yesNo){
	shouldBeDrawn = yesNo;
	return shouldBeDrawn;
}

int CXX_mot::CXXSphereNode::setVertex(const CXXCoord &crd){
	theVertex = crd;
	return 0;
}

void CXX_mot::CXXSphereNode::setIntersector (const CXXCircle *aCircle){
	theIntersector = aCircle;
}

const CXX_mot::CXXCircle *CXX_mot::CXXSphereNode::getIntersector() const{
	return theIntersector;
}

