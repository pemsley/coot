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
#include <list>
#include <map>
#include <math.h>
#ifndef  __MMDB_Manager__
#include "mmdb_manager.h"
#endif
#include "CXXCoord.h"
using namespace std;
#include "CXXAlloc.h"
#include "CXXBall.h"

class CXXCircle;
class CXXCircleNode;
class CXXSurface;
class CXXSphereElement;

class CXXNewHood{
private:
	PCAtom theAtomI;
	const CXXBall *theBall;
	double theRadius;
	double theProbeRadius;
	CXXCoord theCentre;
	std::list<CXXCircle, CXX::CXXAlloc<CXXCircle> >theCircles;
	void init();
public:
	CXXNewHood();
	CXXNewHood(PCAtom centralAtom, double radiusOfAtom1, double probeRadius);
	CXXNewHood(const CXXCircleNode &aNode, double probeRadius);
	
	int addAtom(PCAtom candidate, double radiusOfAtom2);
	int findSegments();
	
	void initWith(PCAtom atomI, double radiusOfAtom1, double probeRadius) {
		theAtomI= atomI;
		theRadius = radiusOfAtom1 + probeRadius;
		theProbeRadius = probeRadius;
		theCentre = CXXCoord(atomI->x, atomI->y, atomI->z);
	};
    void initWith(const CXXCircleNode &aNode, double probeRadius);
	void initWith(const CXXBall *aBall);
    
	const PCAtom getAtomI() const;
	double getRadius() const;
	double getProbeRadius() const;
	int getNCircles() const;
	
	const CXXCoord &getCentre() const;
	
	int addNode(const CXXCircleNode &aNode);
	int addBall(const CXXBall &aBall);
	
	int nCircles() const;
	
	int nNodes() const;
	const CXXCircleNode &getNode(const int iNode) const;
	
	static bool doesNotContainDrawable(const CXXNewHood &aHood);
    
	void identifyUniqueNodes(vector<CXXCircleNode, CXX::CXXAlloc<CXXCircleNode> >&circleNodes, int selHnd) const;
	
    void triangulateAsRegularHoodInto(CXXSurface *aSurface, double delta, const CXXSphereElement *unitSphereAtOrigin) const;
	const std::list<CXXCircle, CXX::CXXAlloc<CXXCircle> > &getCircles() const{
		return theCircles;
	};
	std::list<CXXCircle, CXX::CXXAlloc<CXXCircle> > &getCircles() {
		return theCircles;
	};
	void triangulateAsBallHoodInto(CXXSurface *aSurface, double delta,
								   std::map<const CXXBall*, std::vector<CXXCoord, CXX::CXXAlloc<CXXCoord> > > &raggedEdges, 
								   bool useEdges, int insideOrOutside) ;
	
};
#endif

