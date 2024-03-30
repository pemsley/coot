/*
 *  CXXSphereFlatTriangle.h
 *  CXXSurface
 *
 *  Created by Martin Noble on 23/05/2005.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXXSphereFlatTriangle_included
#define CXXSphereFlatTriangle_included

#include "CXXTriangle.h"
#include "CXXCircleNode.h"


class CXXSphereFlatTriangle : public CXXTriangle{
private:
	const CXXCircle *edgeCircles[3];
	CXXCircleNode circleNodes[3];
public:
	CXXSphereFlatTriangle() : CXXTriangle(){
		
		for (int i=0; i<3; i++){
			edgeCircles[i] = 0;
		}
	};
	
	CXXSphereFlatTriangle(size_t i, size_t j, size_t k, size_t l) : CXXTriangle (i, j, k, l){
		for (size_t i=0; i<3; i++){
			edgeCircles[i] = 0;
		}
	};
	CXXSphereFlatTriangle(int i, int j, int k) : CXXTriangle (i, j, k){
		for (int i=0; i<3; i++){
			edgeCircles[i] = 0;
		}
	};
	void setEdgeCircle(const int i, const CXXCircle *edgeCircle) {
		edgeCircles[i] = edgeCircle;
	};
	const CXXCircle *getEdgeCircle(const int i) const{
		return edgeCircles[i];
	};
	void setCircleNode(const int i, const CXXCircleNode &aNode) {
		circleNodes[i] = aNode;
	};
	const CXXCircleNode &getCircleNode(const int i) const {
		return circleNodes[i];
	};
};

#endif