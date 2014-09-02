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
	//~CXXSphereTriangle();
	
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
