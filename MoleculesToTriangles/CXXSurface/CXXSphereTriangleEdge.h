/*
 * MoleculesToTriangles/CXXSurface/CXXSphereTriangleEdge.h
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


#ifndef CXXSphereTriangleEdge_included
#define CXXSphereTriangleEdge_included
#include <vector>
#include "CXXCoord.h"
#include "mmdb2/mmdb_manager.h"

using namespace std;
class CXXSphereNode;
class CXXSphereElement;
class CXXCircle;
 
class CXXSphereTriangleEdge{
private:
	double edgeLength;
	CXXCoord<CXXCoord_ftype>edgeNormal;
	size_t edgeVertices[2];
	CXXCoord<CXXCoord_ftype>theEdgeCentre;
	CXXCoord<CXXCoord_ftype>theSphereCentre;
//Methods which set or update edgeVertices or edgeNormal *MUST* call this once the final
	//values for vertices and normal are reached
	double calculateLength();
	double edgeRadius;
	CXXSphereElement *theSphereElement;
	CXXCircle *theCircle;
public:
		CXXSphereTriangleEdge();
	~CXXSphereTriangleEdge();
	CXXSphereTriangleEdge(const CXXCoord<CXXCoord_ftype>&aNormal, size_t iV1,size_t iV2,
						  const CXXCoord<CXXCoord_ftype>&anEdgeCentre, const CXXCoord<CXXCoord_ftype>&aSphereCentre,
				 		  double anEdgeRadius,
						  CXXSphereElement *aSphereElement);
	double length() const{
		return edgeLength;
	};
	
	const CXXCoord<CXXCoord_ftype>&normal() const;
	size_t vertex(size_t) const;
	const CXXCoord<CXXCoord_ftype>&edgeCentre() const;
	const CXXCoord<CXXCoord_ftype>&sphereCentre() const;
	CXXCoord<CXXCoord_ftype>midpoint() const;
	double radius() const;
	CXXSphereElement *sphereElement() const ;
	int setSphereElement(CXXSphereElement *se);
	int setSphereCentre (const CXXCoord<CXXCoord_ftype>&crd);
	int setEdgeCentre (const CXXCoord<CXXCoord_ftype>&crd);
	void setLength(double aLength) {
		edgeLength = aLength;
	};
	void setCircle(CXXCircle *aCircle){
		theCircle = aCircle;
	};
	CXXCircle *getCircle() const {
		return theCircle;
	};
	void setVertex(size_t index, size_t value){
		edgeVertices[index] = value;
	};
};

#endif
