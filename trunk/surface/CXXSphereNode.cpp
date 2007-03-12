/*
 *  CXXSphereNode.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
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

