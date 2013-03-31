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
 *  CXXTorusTriangle.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  
 *
 */

#include "CXXTorusTriangle.h"
#include "CXXTorusElement.h"
#include "CXXTorusNode.h"
#include <math.h>
#include <iostream>


CXXTorusTriangle::CXXTorusTriangle()
{
}

CXXTorusTriangle::CXXTorusTriangle(CXXTorusElement *aParent,  int n0, 
								    int n1,  int n2,
								   double dt0, double dt1, double dt2, 
								   double dw0, double dw1, double dw2)
{
	theTorusElement = aParent;
	nodes[0] = n0;
	nodes[1] = n1;
	nodes[2] = n2;
	deltaThetas[0] = dt0;
	deltaThetas[1] = dt1;
	deltaThetas[2] = dt2;
	deltaOmegas[0] = dw0;
	deltaOmegas[1] = dw1;
	deltaOmegas[2] = dw2;
}

CXXTorusTriangle::CXXTorusTriangle(const CXXTorusTriangle &oldTriangle)
{
	theTorusElement = oldTriangle.getTorusElement();
	for (int i=0; i<3; i++){ 
		nodes[i] = oldTriangle.getNode(i);
		deltaThetas[i] = oldTriangle.getDeltaTheta(i);
		deltaOmegas[i] = oldTriangle.getDeltaOmega(i);
	}
}

int CXXTorusTriangle::getNode(const int i) const
{
	return nodes[i];
}

CXXTorusElement *CXXTorusTriangle::getTorusElement(void) const
{
	return theTorusElement;
}

int CXXTorusTriangle::bisect(double degrees){
	int theEdge;
	int debug = 0;
	if ((theEdge = tooLongEdge(degrees)) >= 0){
		if (debug) cout << "Bisecting edge " << theEdge << " of triangle:" << endl;
		if (debug) selfReport();
		double newTheta = theTorusElement->getNode(nodes[theEdge]).getTheta() + (deltaThetas[theEdge]/2.);
		double newOmega = theTorusElement->getNode(nodes[theEdge]).getOmega() + (deltaOmegas[theEdge]/2.);
		
		CXXTorusNode newNode(newTheta, newOmega);
		int iNewNode = theTorusElement->addNode(newNode);
		
//Copy the current triangle, and then modify it before uploading it into the triangle array
		
		CXXTorusTriangle newTriangle1(*this);
		
		newTriangle1.setNode(theEdge, iNewNode);
		newTriangle1.setDeltaTheta(theEdge, deltaThetas[theEdge]/2.);
		newTriangle1.setDeltaOmega(theEdge, deltaOmegas[theEdge]/2.);
		newTriangle1.setDeltaTheta((theEdge+2)%3, -1. * 
								   (newTriangle1.getDeltaTheta(theEdge) + 
									newTriangle1.getDeltaTheta((theEdge+1)%3)));
		newTriangle1.setDeltaOmega((theEdge+2)%3, -1. * 
								   (newTriangle1.getDeltaOmega(theEdge) + 
									newTriangle1.getDeltaOmega((theEdge+1)%3)));
		theTorusElement->addTriangle(newTriangle1);

		// Modify this triangle to be the "other half" ofitself
		
		setNode((theEdge+1)%3, iNewNode);
		setDeltaTheta(theEdge, deltaThetas[theEdge]/2.);
		setDeltaOmega(theEdge, deltaOmegas[theEdge]/2.);
		setDeltaTheta((theEdge+1)%3, -1. * 
								   (getDeltaTheta(theEdge) + 
									getDeltaTheta((theEdge+2)%3)));
		setDeltaOmega((theEdge+1)%3, -1. * 
								   (getDeltaOmega(theEdge) + 
									getDeltaOmega((theEdge+2)%3)));
		return 1;
	}
	else {
		if (debug) cout << "Finished with triangle : " << endl;
		if (debug) selfReport();
		if (debug) cout << "*******" << endl;
		return 0;
	}
}

int CXXTorusTriangle::tooLongEdge(double radians){
	int longestEdge = 0;
	double longestEdgeLength  = -1e30;
	double edgeLength;
		
	for (int i=0; i<3; i++){
		edgeLength = max(fabs(deltaThetas[i]),fabs(deltaOmegas[i]));
		if (edgeLength > longestEdgeLength){
			longestEdge = i;
			longestEdgeLength = edgeLength;
		}
	}
	
	if (longestEdgeLength > radians) return longestEdge;
	else return -1;
}

void CXXTorusTriangle::setNode(int iNode, int aNode){
	nodes[iNode] = aNode;
}

void CXXTorusTriangle::selfReport() const
{
	double R2D = 360.0/(2*M_PI);
	for (int i=0; i<3; i++){
		fprintf(stderr," %8.3f%10.3f\n", 
				theTorusElement->getNode(i).getOmega() * R2D, 
				theTorusElement->getNode(i).getTheta() * R2D);
		fprintf(stderr,"+%8.3f%10.3f\n\n", deltaOmegas[i]*R2D,deltaThetas[i]*R2D);
	}	
}

double CXXTorusTriangle::getDeltaTheta(const int iEdge) const{
	return deltaThetas[iEdge];
}

double CXXTorusTriangle::getDeltaOmega(const int iEdge) const{
	return deltaOmegas[iEdge];
}

void CXXTorusTriangle::setDeltaTheta(const int iEdge, const double newDeltaTheta){
	deltaThetas[iEdge] = newDeltaTheta;
}

void CXXTorusTriangle::setDeltaOmega(const int iEdge, const double newDeltaOmega){
	deltaOmegas[iEdge] = newDeltaOmega;
}

void CXXTorusTriangle::setTorusElement(CXXTorusElement *anElement)
{
	theTorusElement = anElement;
}


