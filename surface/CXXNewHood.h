/*
 *  CXXNewHood.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
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

