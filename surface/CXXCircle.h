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
 *  CXXCircle.h
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  
 *
 */
#ifndef CXXCircle_included
#define CXXCircle_included

#include <vector>
#include <algorithm>
#include <iostream>
#include "CXXCircleNode.h"
#include "CXXCoord.h"
#ifndef  __MMDB_Manager__
#include "mmdb_manager.h"
#endif

using namespace std;
 
class CXXNewHood;

class CXXCircle {
private:
	//Note that order now reflects position in 
	//initialization lists
	PCAtom theAtomJ;
	CXXNewHood *theParent;
	CXXCoord centreOfSphere;
	CXXCoord centreOfSecondSphere;
	CXXCoord theNormal;
	double radiusOfSecondSphere;
	double radiusOfSphere;

	CXXCoord centreOfCircle;
	CXXCoord centreToCircle;
	CXXCoord unitCentreToCircle;
	CXXCoord referenceUnitRadius;
	
	double radiusOfCircle;
	
	
	vector<CXXCircleNode>theNodes;
	vector<int>theStarts;
	vector<int>theStops;

	
	int completelyEaten;
	
	int nodeNumber;
	
public:
		CXXCircle();
	CXXCircle (CXXNewHood *aHood, PCAtom atom2, double radiusOfAtom2, double probeRadius);
	CXXCircle (CXXNewHood *aHood, const CXXCircleNode &aNode, double radiusOfAtom2, int nodeNumber);
	CXXCircle (const CXXCircle &oldOne);
	
	int meetsCircle(const CXXCircle &otherCircle, vector<CXXCoord> &nodeList) const;
	int isSomewhereInsideSphere(const CXXCoord &centre, const double radius) const;
	
	int sortNodes();
	int identifyArcs();
 
	const CXXCoord &getCentreOfSecondSphere() const;
	const CXXCoord &getCentreOfSphere() const;
	const CXXCoord &getCentreToCircle() const;
	const CXXCoord &getUnitCentreToCircle() const;
	const PCAtom    getAtomJ() const;
	
	double getRadiusOfSphere() const;
	double getRadiusOfSecondSphere() const;
	double getRadiusOfCircle() const;
	
	void setRadiusOfCircle(const double &value){
		radiusOfCircle = value;
	};
	const CXXCoord &getCentreOfCircle() const;
	void setCentreOfCircle (const CXXCoord &aCentre){
		centreOfCircle = aCentre;
	};
	
	const CXXCoord &getNormal() const;
	void setNormal(const CXXCoord &aNormal){
		theNormal = aNormal;
	};
	
	CXXNewHood *getParent() const;
	
	unsigned getNNodes() const;
	const CXXCircleNode &getNode(const int i) const;
	int getEaten() const;	
	void setEaten(int flag);	
	 
	unsigned nSegments() const;
	int start(const int i) const;
	int stop(const int i) const;
	int addNode(const CXXCircleNode &aNode);
	double getRadiusOfVdWCircle() const;
	const CXXCoord getCentreOfVdWCircle() const;
	void dumpVdw() const;
	int vdwIsBehind(const CXXCoord &crd) const;
	const CXXCoord vdwPlaneIntersect(const CXXCoord &c1, const CXXCoord &c2) const;
	int accIsBehind(const CXXCoord &crd) const;
	const CXXCoord accPlaneIntersect(const CXXCoord &c1, const CXXCoord &c2) const;
	const CXXCoord &getReferenceUnitRadius() const;
	void setArbitraryReference();
	int getNodeNumber() const {
		return nodeNumber;
	};
};
#endif

