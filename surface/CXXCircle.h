/*
 *  CXXCircle.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXX_mot_CXXCircle_included
#define CXX_mot_CXXCircle_included

#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include "CXXCircleNode.h"
#include "CXXCoord.h"
#ifndef  __MMDB_Manager__
#include <mmdb2/mmdb_manager.h>
#endif

//#include "CXXAlloc.h"

using namespace std;

namespace CXX_mot {

class CXXNewHood;
class CXXBall;

class CXXCircle {
private:
	//Note that order now reflects position in 
	//initialization lists
	mmdb::PAtom theAtomJ;
	const CXXBall* theBallJ;
	CXXNewHood *theParent;
	CXXCoord centreOfSecondSphere;
	CXXCoord theNormal;
	double radiusOfSecondSphere;
	double radiusOfSphere;

	CXXCoord centreOfCircle;
	CXXCoord centreToCircle;
	CXXCoord referenceUnitRadius;
	
	double radiusOfCircle;
	
	
	list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >theNodes;
    int nIntersectingCircles;
	vector<CXXCircleNode *, CXX_old::CXXAlloc<CXXCircleNode *> >theStarts;
	vector<CXXCircleNode *, CXX_old::CXXAlloc<CXXCircleNode *> >theStops;
	
	int completelyEaten;
	
	int nodeNumber;
	
public:
		CXXCircle();
	CXXCircle (CXXNewHood *aHood, mmdb::PAtom atom2, double radiusOfAtom2, double probeRadius);
//	CXXCircle (const CXXCircle &oldOne);
	CXXCircle (CXXNewHood *aHood, const CXXBall &aBall);
	
    void performPrecalculations();
    
	int meetsCircle(const CXXCircle &otherCircle, vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > &nodeList) const;
	int isSomewhereInsideSphere(const CXXCoord &centre, const double radius) const;
	
	int sortNodes();	
	int newIdentifyArcs();
 
	const CXXCoord &getCentreOfSecondSphere() const;
	const CXXCoord &getCentreOfSphere() const;
	const CXXCoord &getCentreToCircle() const;
	const mmdb::PAtom    getAtomJ() const;
	const list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> > &getNodes() const { return theNodes;};
	const CXXBall *getBallJ() const {
		return theBallJ;
	};
	
	double getRadiusOfSphere() const;
	double getRadiusOfSecondSphere() const;
	double getRadiusOfCircle() const{
		return radiusOfCircle;
	};
	
	void setRadiusOfCircle(const double &value){
		radiusOfCircle = value;
	};
	const CXXCoord &getCentreOfCircle() const{
		return centreOfCircle;
	};
	void setCentreOfCircle (const CXXCoord &aCentre){
		centreOfCircle = aCentre;
	};
	
	const CXXCoord &getNormal() const;
	void setNormal(const CXXCoord &aNormal){
		theNormal = aNormal;
	};
	void setCentreToCircle(const CXXCoord &aNormal){
		centreToCircle = aNormal;
	};
	
	CXXNewHood *getParent() const;
	
	unsigned getNNodes() const;
	const CXXCircleNode &getNode(const int i) const;
	CXXCircleNode &getNode(const int i);
	int getEaten() const;	
	void setEaten(int flag);	
	 
	unsigned nSegments() const;
	CXXCircleNode* start(const int i) const;
	CXXCircleNode* stop(const int i) const;
	int addNode(const CXXCircleNode &aNode);
	double getRadiusOfVdWCircle() const;
	const CXXCoord getCentreOfVdWCircle() const;
	void dumpVdw() const;
	int vdwIsBehind(const CXXCoord &crd) const;
	const CXXCoord vdwPlaneIntersect(const CXXCoord &c1, const CXXCoord &c2) const;
	int accIsBehind(const CXXCoord &crd) const{
        CXXCoord diff(crd - getCentreOfCircle());
        return (diff*getNormal() < 0.);
    };    
	const CXXCoord accPlaneIntersect(const CXXCoord &c1, const CXXCoord &c2) const;
	const CXXCoord &getReferenceUnitRadius() const;
	void setArbitraryReference();
    int trimNodesBy(const CXXCircle &otherCircle);
    void trimOwnNodes();
    int countDrawnNodes() const;
    int getNIntersectingCircles() const{
        return nIntersectingCircles;
    };
    bool abBracketsC(const CXXCircleNode &nodea, 
                     const CXXCircleNode &nodeb, 
                     const CXXCircleNode &nodec) const;
    bool smallabBracketsC(const CXXCircleNode &nodea, 
                     const CXXCircleNode &nodeb, 
                     const CXXCircleNode &nodec) const;
    
};

}
#endif

