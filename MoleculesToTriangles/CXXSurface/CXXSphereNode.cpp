/*
 * MoleculesToTriangles/CXXSurface/CXXSphereNode.cpp
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

#include "CXXSphereNode.h"
#include "CXXCoord.h"
#include "CXXTriangle.h"
#include <list>

CXXSphereNode::CXXSphereNode() : 
theVertex(CXXCoord<CXXCoord_ftype>(0.,0.,0.)),
theIntersector(0),
shouldBeDrawn(1),
theAtom(0)
{
//	references.resize(0);
}

CXXSphereNode::CXXSphereNode(const CXXCoord<CXXCoord_ftype>&aCoord) : 
theVertex(aCoord), 
theIntersector(0), 
shouldBeDrawn(1) ,
theAtom(0)
{
//	references.resize(0);
}

const CXXCoord<CXXCoord_ftype>&CXXSphereNode::vertex() const{
	return theVertex;
}

int CXXSphereNode::doDraw() const{
	return shouldBeDrawn;
}

int CXXSphereNode::setDoDraw(const int yesNo){
	shouldBeDrawn = yesNo;
	return shouldBeDrawn;
}

int CXXSphereNode::setVertex(const CXXCoord<CXXCoord_ftype>&crd){
	theVertex = crd;
	return 0;
}

void CXXSphereNode::setIntersector (const CXXCircle *aCircle){
	theIntersector = aCircle;
}

const CXXCircle *CXXSphereNode::getIntersector() const{
	return theIntersector;
}

