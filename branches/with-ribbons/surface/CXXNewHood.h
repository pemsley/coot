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
 *  CXXNewHood.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  
 *
 */

#ifndef CXXNewHood_included
#define CXXNewHood_included
#include <vector>
#include <math.h>
#ifndef  __MMDB_Manager__
#include "mmdb_manager.h"
#endif
#include "CXXCoord.h"
using namespace std;

class CXXCircle;
class CXXCircleNode;

class CXXNewHood{
private:
	PCAtom theAtomI;
	double theRadius;
	double theProbeRadius;
	CXXCoord theCentre;
	vector<CXXCircle>theCircles;
	vector<CXXCircleNode>theNodes;
	void init();
public:
		CXXNewHood();
	~CXXNewHood();
	CXXNewHood(PCAtom centralAtom, double radiusOfAtom1, double probeRadius);
	CXXNewHood(const CXXCircleNode &aNode, double probeRadius);
	
	int addAtom(PCAtom candidate, double radiusOfAtom2);
	int findSegments();
	
	const PCAtom getAtomI() const;
	double getRadius() const;
	double getProbeRadius() const;
	int getNCircles() const;

	const CXXCircle &getCircle(const int i) const;
	const CXXCoord &getCentre() const;
	
	int addNode(const CXXCircleNode &aNode);
	int addNodeAsAtom(const CXXCircleNode &aNode, int aNumber);
		
	int nCircles() const;
	
	int nNodes() const;
	const CXXCircleNode &getNode(const int iNode) const;
	
};
#endif

