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
#include "mmdb2/mmdb_manager.h"

using namespace std;

class CXXSphereNode {
private:
	CXXCoord<CXXCoord_ftype>theVertex;
	const CXXCircle *theIntersector;
//	std::vector<int> references;
	int shouldBeDrawn;
	mmdb::Atom* theAtom;
public:
		CXXSphereNode();
	CXXSphereNode(const CXXCoord<CXXCoord_ftype>&aCoord);
	const CXXCoord<CXXCoord_ftype>&vertex() const;
	int setDoDraw(const int yesNo);
	int doDraw() const;
	int setVertex(const CXXCoord<CXXCoord_ftype>&);
	void setIntersector (const CXXCircle *aCircle);
	const CXXCircle *getIntersector() const;
	int operator < (const CXXSphereNode &comparator) const {
		return (theVertex[2] < comparator.vertex()[2]);
	}
	void setAtom ( mmdb::Atom* anAtom) {
		theAtom = anAtom;
	};
	mmdb::Atom* getAtom() const {
		return theAtom;
	};
};

#endif

