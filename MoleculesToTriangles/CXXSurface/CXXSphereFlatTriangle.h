/*
 * MoleculesToTriangles/CXXSurface/CXXSphereFlatTriangle.h
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