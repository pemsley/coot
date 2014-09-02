/*
 *  CXXSphereNode.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXXSphereNode_included
#define CXXSphereNode_included
#include <vector>
#include <iostream>
#include <algorithm>
#include "CXXCoord.h"

#include "CXXCircle.h"

class CAtom;
using namespace std;

class CXXSphereNode {
private:
	CXXCoord theVertex;
	const CXXCircle *theIntersector;
//	std::vector<int> references;
	int shouldBeDrawn;
	CAtom *theAtom;
public:
		CXXSphereNode();
	CXXSphereNode(const CXXCoord &aCoord);
	const CXXCoord &vertex() const;
	int setDoDraw(const int yesNo);
	int doDraw() const;
	int setVertex(const CXXCoord &);
	void setIntersector (const CXXCircle *aCircle);
	const CXXCircle *getIntersector() const;
	int operator < (const CXXSphereNode &comparator) const {
		return (theVertex[2] < comparator.vertex()[2]);
	}
	void setAtom ( CAtom *anAtom) {
		theAtom = anAtom;
	};
	CAtom *getAtom() const {
		return theAtom;
	};
};

#endif

