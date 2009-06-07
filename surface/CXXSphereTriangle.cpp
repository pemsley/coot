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
 *  CXXSphereTriangle.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  
 *
 */

#include "CXXSphereTriangle.h"
#include <CXXCoord.h>
#include <CXXSphereTriangleEdge.h>
#include <CXXSphereElement.h>
#include <CXXSphereNode.h>

CXXSphereTriangle::CXXSphereTriangle() :
theRadius(0),
theCentre(CXXCoord(0.,0.,0.)){
	
}

CXXSphereTriangle::CXXSphereTriangle(CXXSphereElement *se, int *inputVertices, int *inputEdges, 
									 double inRad, CXXCoord &incent){
	CXXSphereTriangle newTriangle(se, inputVertices, inputEdges, inRad, incent, 0);
	*this = newTriangle;
}

CXXSphereTriangle::CXXSphereTriangle(CXXSphereElement *se, int *inputVertices, int *inputEdges, 
									 double inRad, CXXCoord &incent, CAtom *anAtom){
	theAtom = anAtom;
	theSphereElement = se;
	for (int i=0; i<3; i++){
		triangleVertices[i] = inputVertices[i];
		triangleEdges[i] = inputEdges[i];
	}
	theRadius = inRad;
	theCentre = incent;
}
CXXSphereTriangle::~CXXSphereTriangle(){}
int CXXSphereTriangle::vertex(int i) const{
	return triangleVertices[i];
}
int CXXSphereTriangle::edge(int i) const{
	return triangleEdges[i];
}
double CXXSphereTriangle::radius() const {
	return theRadius;
}
const CXXCoord &CXXSphereTriangle::centre() const {
	return theCentre;
}

int CXXSphereTriangle::bisect(double radians){
	int iLongest = 0;
	double longestLength = -1e30;
	for (int i=0; i<3; i++){
		if (theSphereElement->edge(triangleEdges[i]).length()>longestLength){
			longestLength = theSphereElement->edge(triangleEdges[i]).length();
			iLongest = i;
		}
	}
	if (longestLength>radians){
		CXXSphereNode newVertex(theSphereElement->edge(triangleEdges[iLongest]).midpoint());
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
		if (newNormal.dot(theSphereElement->vertex(triangleVertices[iLongest]).vertex()-theCentre)<0.) 
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

int CXXSphereTriangle::setVertex(int i, int newOne){
	triangleVertices[i] = newOne;
	return 0;
}

CXXSphereElement *CXXSphereTriangle::sphereElement() const{
	return theSphereElement;
}

int CXXSphereTriangle::setEdge(int i, int newOne){
	triangleEdges[i] = newOne;
	return 0;
}

int CXXSphereTriangle::setSphereElement(CXXSphereElement *se){
	theSphereElement = se;
	return 0;
}

int CXXSphereTriangle::setCentre(const CXXCoord &crd){
	theCentre = crd;
	return 0;
}

int CXXSphereTriangle::setRadius(const double rad){
	theRadius = rad;
	return 0;
}

CAtom *CXXSphereTriangle::getAtom() const {
	return theAtom;
}

int CXXSphereTriangle::setAtom(CAtom *anAtom) {
	theAtom = anAtom;
	return 0;
}
