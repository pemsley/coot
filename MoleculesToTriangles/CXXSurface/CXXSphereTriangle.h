/*
 * MoleculesToTriangles/CXXSurface/CXXSphereTriangle.h
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

#ifndef CXXSphereTriangle_included
#define CXXSphereTriangle_included

#include "CXXCoord.h"
#include "CXXSphereTriangleEdge.h"
//#include "CXXSphereElement.h"
#include "CXXSphereNode.h"
//class CXXSphereTriangleEdge;
class CXXSphereElement;
//class CXXCoord;
//class CXXSphereNode;
#include "mmdb2/mmdb_manager.h"

class CXXSphereTriangle {
private:
	size_t triangleVertices[3];
	size_t triangleEdges[3];
	double theRadius;
	CXXCoord<CXXCoord_ftype>theCentre;
	CXXSphereElement *theSphereElement;
    mmdb::Atom* theAtom;

public:
		CXXSphereTriangle();
	CXXSphereTriangle(CXXSphereElement *se, size_t *vertices, size_t *edges,
					  double aRadius, CXXCoord<CXXCoord_ftype>&aCentre);
	CXXSphereTriangle(CXXSphereElement *se, size_t *vertices, size_t *edges,
					  double aRadius, CXXCoord<CXXCoord_ftype>&aCentre, mmdb::Atom* anAtom);
	//~CXXSphereTriangle();
	
	size_t vertex(size_t) const;
	size_t edge(size_t) const;
	double radius() const;
	const CXXCoord<CXXCoord_ftype>&centre() const;
	CXXSphereElement *sphereElement() const;
	int bisect(double);
	mmdb::Atom* getAtom() const;
	
	int setAtom(mmdb::Atom* anAtom);
	int setEdge(size_t i, size_t);
	int setVertex(size_t, size_t);
	int setSphereElement(CXXSphereElement *se);
	int setCentre(const CXXCoord<CXXCoord_ftype>&);
	int setRadius(const double radius);
};
#endif
