/*
 *  CXXSphereTriangleEdge.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
#include "CXXSphereTriangleEdge.h"
#include "CXXSphereElement.h"
#include "CXXSphereNode.h"
#include "CXXCoord.h"

CXXSphereTriangleEdge::CXXSphereTriangleEdge(){
}

CXXSphereTriangleEdge::~CXXSphereTriangleEdge(){
}

CXXSphereTriangleEdge::CXXSphereTriangleEdge(const CXXCoord &aNormal, int iV1, int iV2, 
											 const CXXCoord &anEdgeCentre, const CXXCoord &aSphereCentre,
											 double anEdgeRadius, CXXSphereElement *aSphereElement) :
edgeNormal ( aNormal),
theEdgeCentre ( anEdgeCentre),
theSphereCentre ( aSphereCentre),
edgeRadius ( anEdgeRadius),
theSphereElement ( aSphereElement )
{
	edgeVertices[0] = iV1;
	edgeVertices[1] = iV2;
	calculateLength();
}

const CXXCoord &CXXSphereTriangleEdge::normal() const{
	return edgeNormal;
}
int CXXSphereTriangleEdge::vertex(int i) const{
	return edgeVertices[i];
}
const CXXCoord &CXXSphereTriangleEdge::edgeCentre() const {
	return theEdgeCentre;
}
const CXXCoord &CXXSphereTriangleEdge::sphereCentre() const {
	return theSphereCentre;
}

double CXXSphereTriangleEdge::calculateLength()  {
	CXXCoord v1 = theSphereElement->vertex(edgeVertices[0]).vertex() - theEdgeCentre;
	CXXCoord v2 = theSphereElement->vertex(edgeVertices[1]).vertex() - theEdgeCentre;
	v1.normalise();
	v2.normalise();
	edgeLength = edgeNormal.angleBetween(v1, v2);
//	cout << edgeLength*360./(2.*M_PI) << endl;
	return edgeLength;
}

CXXCoord CXXSphereTriangleEdge::midpoint() const{
	CXXCoord result;
	CXXCoord v1 = theSphereElement->vertex(edgeVertices[0]).vertex() - theEdgeCentre;
	CXXCoord v2 = theSphereElement->vertex(edgeVertices[1]).vertex() - theEdgeCentre;
	
	// If the radii from the centre of the edge's circle to the vertices
	// are not parallel ( or anti parallel), then there is an easy way 
	// to identify the midpoint.  
	if (fabs(edgeLength-0.)> 1e-8 && fabs(edgeLength-M_PI)>1e-8){
		CXXCoord v3 = v1 + v2;
		v3.normalise();
		v3.scale(edgeRadius);
		// Answeer is plus if length > PI
		if (edgeLength>M_PI) v3.scale (-1.);
		result = theEdgeCentre + v3;
	}
	else if (fabs(edgeLength-0.)<= 1e-8){
		result = theSphereElement->vertex(edgeVertices[0]).vertex();
	}
	else {
		result = edgeNormal^v1;
		result = result + theEdgeCentre;
	}
	return CXXCoord(result);
}

double CXXSphereTriangleEdge::radius() const{
	return edgeRadius;
}

CXXSphereElement *CXXSphereTriangleEdge::sphereElement()const{
	return theSphereElement;
}

int CXXSphereTriangleEdge::setSphereElement(CXXSphereElement *se){
	theSphereElement = se;
	return 0;
}

int CXXSphereTriangleEdge::setSphereCentre (const CXXCoord &crd){
	theSphereCentre = crd;
	return 0;
}

int CXXSphereTriangleEdge::setEdgeCentre (const CXXCoord &crd){
	theEdgeCentre = crd;
	return 0;
}

