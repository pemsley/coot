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
 *  CXXNewHood.cpp
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  
 *
 */

#if defined _MSC_VER
#include <string>
#endif

#include "CXXNewHood.h"
#include "mmdb_manager.h"
#include "mmdb_tables.h"
#include <CXXCircle.h>
#include <CXXCircleNode.h>

void CXXNewHood::init(){
	theAtomI = 0;
	theRadius = 0;
	theCircles.resize(0);
	theNodes.resize(0);
	theCentre = CXXCoord(0.,0.,0.);
}

CXXNewHood::CXXNewHood(){
	init();
}
CXXNewHood::~CXXNewHood(){}
CXXNewHood::CXXNewHood(PCAtom centralAtom, double radiusOfAtom1, double probeRadius) :
theAtomI(centralAtom), 
theRadius(radiusOfAtom1+probeRadius), 
theProbeRadius(probeRadius){
	theCentre = CXXCoord(theAtomI->x, theAtomI->y, theAtomI->z);
}

CXXNewHood::CXXNewHood(const CXXCircleNode &aNode, double probeRadius) :
theAtomI(0), 
theRadius(probeRadius), 
theProbeRadius(0),
theCentre(aNode.getCoord()) {
}

int CXXNewHood::addAtom(PCAtom anAtomJ, double radiusOfAtom2){
	if (anAtomJ->serNum == theAtomI->serNum) {
		//		std::cout << "Rejecting self " << anAtomJ->serNum << " " <<  theAtomI->serNum << endl;
		return 0; //Worried this might not be unique
	}
	double radiusOfAtomJ = radiusOfAtom2 + theProbeRadius;
	CXXCoord centreOfAtomJ = CXXCoord(anAtomJ->x, anAtomJ->y, anAtomJ->z);
	CXXCoord centreOfAtomI = theCentre;
	//This test is to deal with unlikely but possible event of being fed asphere with identical coordinates
	if (centreOfAtomI != centreOfAtomJ){ 
		CXXCoord vectorIJ = centreOfAtomJ - centreOfAtomI;
		double separ;
		if ((separ = vectorIJ.get3DLengthSq()) < (theRadius+radiusOfAtomJ)*(theRadius+radiusOfAtomJ)){
			if (separ > 0.00001){
				theCircles.push_back(CXXCircle(this, anAtomJ, radiusOfAtom2, theProbeRadius));
				return 1;
			}
		}
	}
	
	return 0;
}

int CXXNewHood::addNodeAsAtom(const CXXCircleNode &aNode, int aNumber){
	if (aNode.getCoord() != theCentre){
		CXXCoord vectorIJ = aNode.getCoord() - theCentre;
		double separ;
		if ((separ = vectorIJ.get3DLengthSq()) < (2.*theRadius)*(2.*theRadius)){
			if (separ > 0.00001){
				theCircles.push_back(CXXCircle(this, aNode, theRadius, aNumber));
				theCircles.back().setArbitraryReference();
				return 1;
			}
		}
	}
	return 0;
}

int CXXNewHood::findSegments(){
	//First cause circles to find their nodes
	
	for (unsigned int iCircle1=0;
		 iCircle1<theCircles.size();
		 iCircle1++){
		CXXCircle &circle1(theCircles[iCircle1]);
		if (!circle1.getEaten()){
			for (unsigned int iCircle2=0; 
				 iCircle2<theCircles.size() && !circle1.getEaten(); 
				 iCircle2++){
				CXXCircle &circle2(theCircles[iCircle2]);
				if (iCircle1!=iCircle2){
					vector<CXXCoord> theIntersections(2);
					int eatenBy=circle1.meetsCircle(circle2, theIntersections);
					if (eatenBy == 2){
						circle1.setEaten(1);
					}
					else if (eatenBy == 1){
						circle2.setEaten(1);
					}
					else if (eatenBy==0){
						circle1.addNode(CXXCircleNode(&circle1,
													  circle2.getAtomJ(),
													  theIntersections[0],
													  1));
						circle1.addNode(CXXCircleNode(&circle1,
													  circle2.getAtomJ(),
													  theIntersections[1],
													  2));
					}
				}
			}
		}
		// Here identify case where a circle has simply not intersected with another circle
		//and so is an intact 2*PI circle
		if (!(circle1.getEaten()) && circle1.getNNodes()==0){
			circle1.setArbitraryReference();
		}
	}
	
	for (unsigned int iCircle1=0; iCircle1<theCircles.size(); iCircle1++){
		//Don't bother with this bit if  the circle has been eaten
		if (!theCircles[iCircle1].getEaten()){
			//Then make the circles sort these
			theCircles[iCircle1].sortNodes();
			//Then make the circles upload their unique arcs to us
			theCircles[iCircle1].identifyArcs();
		}
	}
	return 0;
}

const CXXCoord &CXXNewHood::getCentre() const{
	return theCentre;
}

double CXXNewHood::getRadius() const{
	return theRadius;
}
int CXXNewHood::nCircles() const{
	return theCircles.size();
}
const PCAtom CXXNewHood::getAtomI() const {
	return theAtomI;
}
int CXXNewHood::addNode(const CXXCircleNode &node){
	int oldNNodes = theNodes.size();
	theNodes.push_back(node);
	return oldNNodes;
}
int CXXNewHood::nNodes() const{
	return theNodes.size();
}

const CXXCircleNode &CXXNewHood::getNode(const int iNode) const{
	return theNodes[iNode];
}
double CXXNewHood::getProbeRadius() const{
	return theProbeRadius;
}

const CXXCircle &CXXNewHood::getCircle(const int i) const{
	return theCircles[i];
}

