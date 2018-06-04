/*
 *  CXXNewHood.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXX_mot_CXXNewHood_included
#define CXX_mot_CXXNewHood_included
#include <vector>
#include <list>
#include <map>
#include <math.h>
#ifndef  __MMDB_Manager__
#include <mmdb2/mmdb_manager.h>
#endif
#include "CXXCoord.h"
using namespace std;
#include "CXXAlloc.h"
#include "CXXBall.h"

namespace CXX_mot {

class CXXCircle;
class CXXCircleNode;
class CXXSurface;
class CXXSphereElement;

class CXXNewHood{
private:
	mmdb::PAtom theAtomI;
	const CXXBall *theBall;
	double theRadius;
	double theProbeRadius;
	CXXCoord theCentre;
	std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >theCircles;
	void init();
public:
	CXXNewHood();
	CXXNewHood(mmdb::PAtom centralAtom, double radiusOfAtom1, double probeRadius);
	CXXNewHood(const CXXCircleNode &aNode, double probeRadius);
	
	int addAtom(mmdb::PAtom candidate, double radiusOfAtom2);
	int findSegments();
	
	void initWith(mmdb::PAtom atomI, double radiusOfAtom1, double probeRadius) {
		theAtomI= atomI;
		theRadius = radiusOfAtom1 + probeRadius;
		theProbeRadius = probeRadius;
		theCentre = CXXCoord(atomI->x, atomI->y, atomI->z);
	};
    void initWith(const CXXCircleNode &aNode, double probeRadius);
	void initWith(const CXXBall *aBall);
    
	const mmdb::PAtom getAtomI() const;
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
    
	void identifyUniqueNodes(vector<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >&circleNodes, int selHnd) const;
	
    void triangulateAsRegularHoodInto(CXXSurface *aSurface, double delta, const CXXSphereElement *unitSphereAtOrigin) const;
	const std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> > &getCircles() const{
		return theCircles;
	};
	std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> > &getCircles() {
		return theCircles;
	};
	void triangulateAsBallHoodInto(CXXSurface *aSurface, double delta,
								   std::map<const CXXBall*, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > > &raggedEdges, 
								   bool useEdges, int insideOrOutside) ;
	
};
}
#endif

