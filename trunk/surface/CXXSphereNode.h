/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
 *  CXXSphereNode.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  
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
	CXXSphereNode(CXXCoord &aCoord);
	CXXSphereNode(const CXXCoord &aCoord);
	~CXXSphereNode();
	const CXXCoord &vertex() const;
/*
	int nReferences() const;
	vector<int>::const_iterator beginReference() const;
	vector<int>::const_iterator endReference() const; 
	int reference(const int i) const;
	int eraseReferences();
	int addReference(int);
	int removeReference(int);
*/
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

