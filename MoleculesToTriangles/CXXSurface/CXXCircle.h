/*
 *  CXXCircle.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXXCircle_included
#define CXXCircle_included

#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include "CXXCircleNode.h"
#include "CXXCoord.h"
#ifndef  __MMDB_Manager__
#include "mmdb2/mmdb_manager.h"
#endif


using namespace std;
 
class CXXNewHood;
class CXXBall;

class CXXCircle {
private:
	//Note that order now reflects position in 
	//initialization lists
	mmdb::Atom* theAtomJ;
	const CXXBall* theBallJ;
	CXXNewHood *theParent;
	CXXCoord<CXXCoord_ftype>centreOfSecondSphere;
	CXXCoord<CXXCoord_ftype>theNormal;
	double radiusOfSecondSphere;
	double radiusOfSphere;

	CXXCoord<CXXCoord_ftype>centreOfCircle;
	CXXCoord<CXXCoord_ftype>centreToCircle;
	CXXCoord<CXXCoord_ftype>referenceUnitRadius;
	
	double radiusOfCircle;
	
	
	list<CXXCircleNode  >theNodes;
    int nIntersectingCircles;
	vector<CXXCircleNode *  >theStarts;
	vector<CXXCircleNode *  >theStops;
	
	int completelyEaten;
	
	int nodeNumber;
    
    int containsEatenNodes;
	
public:
		CXXCircle();
	CXXCircle (CXXNewHood *aHood, mmdb::Atom* atom2, double radiusOfAtom2, double probeRadius);
//	CXXCircle (const CXXCircle &oldOne);
	CXXCircle (CXXNewHood *aHood, const CXXBall &aBall);
	
    void performPrecalculations();
    
	int meetsCircle(const CXXCircle &otherCircle, vector<CXXCoord<CXXCoord_ftype> > &nodeList) const;
	int isSomewhereInsideSphere(const CXXCoord<CXXCoord_ftype>&centre, const double radius) const;
	
	int sortNodes();	
	int newIdentifyArcs();
 
	const CXXCoord<CXXCoord_ftype>&getCentreOfSecondSphere() const;
	const CXXCoord<CXXCoord_ftype>&getCentreOfSphere() const;
	const CXXCoord<CXXCoord_ftype>&getCentreToCircle() const;
	mmdb::Atom*    getAtomJ() const;
	const list<CXXCircleNode  > &getNodes() const { return theNodes;};
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
	const CXXCoord<CXXCoord_ftype>&getCentreOfCircle() const{
		return centreOfCircle;
	};
	void setCentreOfCircle (const CXXCoord<CXXCoord_ftype>&aCentre){
		centreOfCircle = aCentre;
	};
	
	const CXXCoord<CXXCoord_ftype>&getNormal() const;
	void setNormal(const CXXCoord<CXXCoord_ftype>&aNormal){
		theNormal = aNormal;
	};
	void setCentreToCircle(const CXXCoord<CXXCoord_ftype>&aNormal){
		centreToCircle = aNormal;
	};
	
	CXXNewHood *getParent() const;
	
	size_t getNNodes() const;
	const CXXCircleNode &getNode(const int i) const;
	CXXCircleNode &getNode(const int i);
	int getEaten() const;	
	void setEaten(int flag);	
    void setContainsEatenNodes(const int flag) {
        containsEatenNodes = flag;
    };
    int getContainsEatenNodes() const {
        return containsEatenNodes;
    };
	size_t nSegments() const;
	CXXCircleNode* start(const int i) const;
	CXXCircleNode* stop(const int i) const;
	int addNode(const CXXCircleNode &aNode);
	double getRadiusOfVdWCircle() const;
	const CXXCoord<CXXCoord_ftype>getCentreOfVdWCircle() const;
	void dumpVdw() const;
	int vdwIsBehind(const CXXCoord<CXXCoord_ftype>&crd) const;
	const CXXCoord<CXXCoord_ftype>vdwPlaneIntersect(const CXXCoord<CXXCoord_ftype>&c1, const CXXCoord<CXXCoord_ftype>&c2) const;
	int accIsBehind(const CXXCoord<CXXCoord_ftype>&crd) const{
        CXXCoord<CXXCoord_ftype>diff(crd - getCentreOfCircle());
        return (diff*getNormal() < 0.);
    };    
	const CXXCoord<CXXCoord_ftype>accPlaneIntersect(const CXXCoord<CXXCoord_ftype>&c1, const CXXCoord<CXXCoord_ftype>&c2) const;
	const CXXCoord<CXXCoord_ftype>&getReferenceUnitRadius() const;
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
#endif

