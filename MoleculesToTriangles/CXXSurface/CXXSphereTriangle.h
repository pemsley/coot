/*
 *  CXXSphereTriangle.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
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
