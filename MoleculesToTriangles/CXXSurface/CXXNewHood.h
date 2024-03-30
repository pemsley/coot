/*
 * MoleculesToTriangles/CXXSurface/CXXNewHood.h
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

#ifndef CXXNewHood_included
#define CXXNewHood_included
#include <vector>
#include <list>
#include <map>
#include <math.h>
#ifndef  __MMDB_Manager__
#include "mmdb2/mmdb_manager.h"
#endif
#include "CXXCoord.h"
using namespace std;
#include "CXXBall.h"

class CXXCircle;
class CXXCircleNode;
class CXXSurface;
class CXXSphereElement;

class CXXNewHood{
private:
	mmdb::Atom* theAtomI;
	const CXXBall *theBall;
	double theRadius;
	double theProbeRadius;
	CXXCoord<CXXCoord_ftype>theCentre;
	std::list<CXXCircle  >theCircles;
	void init();
public:
	CXXNewHood();
	CXXNewHood(mmdb::Atom* centralAtom, double radiusOfAtom1, double probeRadius);
	CXXNewHood(const CXXCircleNode &aNode, double probeRadius);
	
	int addAtom(mmdb::Atom* candidate, double radiusOfAtom2);
	int findSegments();
	
	void initWith(mmdb::Atom* atomI, double radiusOfAtom1, double probeRadius) {
		theAtomI= atomI;
		theRadius = radiusOfAtom1 + probeRadius;
		theProbeRadius = probeRadius;
		theCentre = CXXCoord<CXXCoord_ftype>(atomI->x, atomI->y, atomI->z);
	};
    void initWith(const CXXCircleNode &aNode, double probeRadius);
	void initWith(const CXXBall *aBall);
    
	mmdb::Atom* getAtomI() const;
	double getRadius() const;
	double getProbeRadius() const;
	int getNCircles() const;
	
	const CXXCoord<CXXCoord_ftype>&getCentre() const;
	
	int addNode(const CXXCircleNode &aNode);
	int addBall(const CXXBall &aBall);
	
	size_t nCircles() const;
	
	int nNodes() const;
	const CXXCircleNode &getNode(const int iNode) const;
	
	static bool containsDrawable(const CXXNewHood &aHood);
    
	void identifyUniqueNodes(vector<CXXCircleNode  >&circleNodes, int selHnd) const;
	
    void triangulateAsRegularHoodInto(CXXSurface &aSurface, double delta, const CXXSphereElement *unitSphereAtOrigin) const;
	const std::list<CXXCircle  > &getCircles() const{
		return theCircles;
	};
	std::list<CXXCircle  > &getCircles() {
		return theCircles;
	};
	void triangulateAsBallHoodInto(CXXSurface &aSurface, double delta,
								   std::map<const CXXBall*, std::vector<CXXCoord<CXXCoord_ftype> > > &raggedEdges,
								   bool useEdges, int insideOrOutside, const CXXSphereElement &unitCellAtOriginForDelta) ;
	
};
#endif

