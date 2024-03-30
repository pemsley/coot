/*
 * MoleculesToTriangles/CXXSurface/CXXSphereNode.h
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

