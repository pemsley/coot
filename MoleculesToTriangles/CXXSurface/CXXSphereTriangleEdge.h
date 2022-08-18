/*
 *  CXXSphereTriangleEdge.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
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
