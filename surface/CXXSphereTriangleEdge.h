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
 *  CXXSphereTriangleEdge.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  
 *
 */


#ifndef CXXSphereTriangleEdge_included
#define CXXSphereTriangleEdge_included
#include <vector>
#include "CXXCoord.h"

using namespace std;
class CXXSphereNode;
class CXXSphereElement;
class CAtom;
class CXXCircle;
 
class CXXSphereTriangleEdge{
private:
	double edgeLength;
	CXXCoord edgeNormal;
	int edgeVertices[2];
	CXXCoord theEdgeCentre;
	CXXCoord theSphereCentre;
//Methods which set or update edgeVertices or edgeNormal *MUST* call this once the final
	//values for vertices and normal are reached
	double calculateLength();
	double edgeRadius;
	CXXSphereElement *theSphereElement;
	CXXCircle *theCircle;
public:
		CXXSphereTriangleEdge();
	~CXXSphereTriangleEdge();
	CXXSphereTriangleEdge(const CXXCoord &aNormal, int iV1,int iV2, 
						  const CXXCoord &anEdgeCentre, const CXXCoord &aSphereCentre,
				 		  double anEdgeRadius,
						  CXXSphereElement *aSphereElement);
	double length() const{
		return edgeLength;
	};
	
	const CXXCoord &normal() const;
	int vertex(int) const;
	const CXXCoord &edgeCentre() const;
	const CXXCoord &sphereCentre() const;
	CXXCoord midpoint() const;
	double radius() const;
	CXXSphereElement *sphereElement() const ;
	int setSphereElement(CXXSphereElement *se);
	int setSphereCentre (const CXXCoord &crd);
	int setEdgeCentre (const CXXCoord &crd);
	void setLength(double aLength) {
		edgeLength = aLength;
	};
	void setCircle(CXXCircle *aCircle){
		theCircle = aCircle;
	};
	CXXCircle *getCircle() const {
		return theCircle;
	};
	void setVertex(int index, int value){
		edgeVertices[index] = value;
	};
};

#endif
