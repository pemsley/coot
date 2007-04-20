/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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
 *  CXXSphereNode.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  
 *
 */

#include "CXXSphereNode.h"
#include <CXXCoord.h>
#include <CXXTriangle.h>
#include <list>

CXXSphereNode::CXXSphereNode() : 
theVertex(CXXCoord(0.,0.,0.)),
theIntersector(0),
shouldBeDrawn(1),
theAtom(0)
{
//	references.resize(0);
}

CXXSphereNode::CXXSphereNode(CXXCoord &aCoord) : 
theVertex(aCoord), 
theIntersector(0), 
shouldBeDrawn(1) ,
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

CXXSphereNode::~CXXSphereNode(){}

const CXXCoord &CXXSphereNode::vertex() const{
	return theVertex;
}
/*
int CXXSphereNode::nReferences() const{
	return references.size();
}

vector<int>::const_iterator CXXSphereNode::beginReference() const {
	return references.begin();
}

vector<int>::const_iterator CXXSphereNode::endReference() const {
	return references.end();
}

int CXXSphereNode::reference(const int i) const {
	return references[i];
}

int CXXSphereNode::addReference (int iReferringTriangle){
//	cout << "Count was " << references.size() ;
	references.push_back(iReferringTriangle);
//	cout << " now " << references.size() << endl;
	return references.size();
}

int CXXSphereNode::removeReference (int iReferringTriangle){
//	cout << "Count was " << references.size() ;
	//references.remove(iReferringTriangle);
//	cout << " now " << references.size() << endl;
	return references.size();
}
*/
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
/*
int CXXSphereNode::eraseReferences(){
	references.resize(0);
	return 0;
}
*/
void CXXSphereNode::setIntersector (const CXXCircle *aCircle){
	theIntersector = aCircle;
}

const CXXCircle *CXXSphereNode::getIntersector() const{
	return theIntersector;
}

