/*
 *  CXXSphereTriangleEdge.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef CXX_mot_CXXSphereTriangleEdge_included
#define CXX_mot_CXXSphereTriangleEdge_included
#include <vector>
#include "CXXCoord.h"

using namespace std;

namespace CXX_mot {
class CXXSphereNode;
class CXXSphereElement;
// class CAtom;
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

}
#endif
