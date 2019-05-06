/*
 *  CXXSphereTriangle.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXSphereTriangle.h"
#include "CXXCoord.h"
#include "CXXSphereTriangleEdge.h"
#include "CXXSphereElement.h"
#include "CXXSphereNode.h"

CXX_mot::CXXSphereTriangle::CXXSphereTriangle() :
theRadius(0),
theCentre(CXXCoord(0.,0.,0.)){
	
}

CXX_mot::CXXSphereTriangle::CXXSphereTriangle(CXXSphereElement *se, int *inputVertices, int *inputEdges, 
									 double inRad, CXXCoord &incent){
	CXXSphereTriangle newTriangle(se, inputVertices, inputEdges, inRad, incent, 0);
	*this = newTriangle;
}

CXX_mot::CXXSphereTriangle::CXXSphereTriangle(CXXSphereElement *se, int *inputVertices, int *inputEdges, 
									 double inRad, CXXCoord &incent, mmdb::Atom *anAtom){
	theAtom = anAtom;
	theSphereElement = se;
	for (int i=0; i<3; i++){
		triangleVertices[i] = inputVertices[i];
		triangleEdges[i] = inputEdges[i];
	}
	theRadius = inRad;
	theCentre = incent;
}
//CXX_mot::CXXSphereTriangle::~CXXSphereTriangle(){}
int CXX_mot::CXXSphereTriangle::vertex(int i) const{
	return triangleVertices[i];
}
int CXX_mot::CXXSphereTriangle::edge(int i) const{
	return triangleEdges[i];
}
double CXX_mot::CXXSphereTriangle::radius() const {
	return theRadius;
}
const CXX_mot::CXXCoord &CXX_mot::CXXSphereTriangle::centre() const {
	return theCentre;
}

int CXX_mot::CXXSphereTriangle::bisect(double radians){
	int iLongest = 0;
	double longestLength = -1e30;
	for (int i=0; i<3; i++){
		if (theSphereElement->edge(triangleEdges[i]).length()>longestLength){
			longestLength = theSphereElement->edge(triangleEdges[i]).length();
			iLongest = i;
		}
	}
	if (longestLength>radians){
		const CXXSphereTriangleEdge &theEdge(theSphereElement->edge(triangleEdges[iLongest]));
		const CXXSphereNode &vertex0 = theSphereElement->vertex(theEdge.vertex(0));
		CXXSphereNode newVertex(theEdge.midpoint());
		newVertex.setAtom(vertex0.getAtom());
		int iNewVertex = theSphereElement->addVertex(newVertex);
		
		double oldLength = theSphereElement->edge(triangleEdges[iLongest]).length();
		CXXSphereTriangleEdge newEdge0(theSphereElement->edge(triangleEdges[iLongest]));
		newEdge0.setVertex(1, iNewVertex);
		newEdge0.setLength(oldLength/2.);
		int iNewEdge0 = theSphereElement->addEdge(newEdge0);
		
		CXXSphereTriangleEdge newEdge1(theSphereElement->edge(triangleEdges[iLongest]));
		newEdge1.setVertex(0, iNewVertex);
		newEdge1.setLength(oldLength/2.);
		int iNewEdge1 = theSphereElement->addEdge(newEdge1);
		
		CXXCoord v1 = newVertex.vertex() - theCentre;
		CXXCoord v2 = theSphereElement->vertex(triangleVertices[(iLongest+2)%3]).vertex() - theCentre;
		CXXCoord newNormal = v1^v2;
		newNormal.normalise();
		if (newNormal*(theSphereElement->vertex(triangleVertices[iLongest]).vertex()-theCentre)<0.) 
			newNormal *= -1.;
		CXXSphereTriangleEdge newEdge2(newNormal, 
									   iNewVertex, 
									   vertex((iLongest+2)%3), 
									   theCentre, theCentre, theRadius,
									   theSphereElement);
		newEdge2.setCircle(0);	
		int iNewEdge2 = theSphereElement->addEdge(newEdge2);
		
		newNormal.scale(-1);
		CXXSphereTriangleEdge newEdge3(newNormal, 
									   vertex((iLongest+2)%3), 
									   iNewVertex, 
									   theCentre, theCentre, theRadius,
									   theSphereElement);
		newEdge3.setCircle(0);
		int iNewEdge3 = theSphereElement->addEdge(newEdge3);
		
		//Copy this triangle
		CXXSphereTriangle newTriangle(*this);

		//Update this triangle
		setVertex((iLongest+1)%3, iNewVertex);
		setEdge(iLongest, iNewEdge0);
		setEdge((iLongest+1)%3, iNewEdge2);

		//Update the new triangle
		newTriangle.setVertex(iLongest, iNewVertex);
		newTriangle.setEdge(iLongest, iNewEdge1);
		newTriangle.setEdge((iLongest+2)%3, iNewEdge3);
		
		theSphereElement->addTriangle(newTriangle);
		
		return 1;
	}
	else return 0;
}

int CXX_mot::CXXSphereTriangle::setVertex(int i, int newOne){
	triangleVertices[i] = newOne;
	return 0;
}

CXX_mot::CXXSphereElement *CXX_mot::CXXSphereTriangle::sphereElement() const{
	return theSphereElement;
}

int CXX_mot::CXXSphereTriangle::setEdge(int i, int newOne){
	triangleEdges[i] = newOne;
	return 0;
}

int CXX_mot::CXXSphereTriangle::setSphereElement(CXXSphereElement *se){
	theSphereElement = se;
	return 0;
}

int CXX_mot::CXXSphereTriangle::setCentre(const CXXCoord &crd){
	theCentre = crd;
	return 0;
}

int CXX_mot::CXXSphereTriangle::setRadius(const double rad){
	theRadius = rad;
	return 0;
}

mmdb::Atom *CXX_mot::CXXSphereTriangle::getAtom() const {
	return theAtom;
}

int CXX_mot::CXXSphereTriangle::setAtom(mmdb::Atom *anAtom) {
	theAtom = anAtom;
	return 0;
}
