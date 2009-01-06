/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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
 *  CXXSphereTriangle.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  
 *
 */

#ifndef CXXSphereTriangle_included
#define CXXSphereTriangle_included

#include "CXXCoord.h"
#include "CXXSphereTriangleEdge.h"
#include "CXXSphereElement.h"
#include "CXXSphereNode.h"
//class CXXSphereTriangleEdge;
//class CXXSphereElement;
//class CXXCoord;
//class CXXSphereNode;

class CXXSphereTriangle {
private:
	int triangleVertices[3];
	int triangleEdges[3];
	double theRadius;
	CXXCoord theCentre;
	CXXSphereElement *theSphereElement;
	class CAtom *theAtom;
public:
		CXXSphereTriangle();
	CXXSphereTriangle(CXXSphereElement *se, int *vertices, int *edges, 
					  double aRadius, CXXCoord &aCentre);
	CXXSphereTriangle(CXXSphereElement *se, int *vertices, int *edges, 
					  double aRadius, CXXCoord &aCentre, CAtom *anAtom);
	~CXXSphereTriangle();
	
	int vertex(int) const;
	int edge(int) const;
	double radius() const;
	const CXXCoord &centre() const;
	CXXSphereElement *sphereElement() const;
	int bisect(double);
	CAtom *getAtom() const;
	
	int setAtom(CAtom *anAtom);
	int setEdge(int i, int);
	int setVertex(int i, int);
	int setSphereElement(CXXSphereElement *se);
	int setCentre(const CXXCoord &);
	int setRadius(const double radius);
};
#endif
