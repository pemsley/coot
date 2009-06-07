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
 *  CXXCircle.cpp
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  
 *
 */

#include <math.h>
#include "CXXCircle.h"
#include <CXXNewHood.h>
#include "mmdb_manager.h"
#include "mmdb_tables.h"

CXXCircle::CXXCircle() :
theAtomJ(0),
theParent(0),
completelyEaten(0)
{
}

CXXCircle::CXXCircle (CXXNewHood *aHood, PCAtom atom2, double radiusOfAtom2, double probeRadius) :
theAtomJ(atom2),
theParent(aHood),
centreOfSphere(aHood->getCentre()),
radiusOfSphere(aHood->getRadius()),
completelyEaten(0)
{
	centreOfSecondSphere = CXXCoord(theAtomJ->x, theAtomJ->y, theAtomJ->z); 
	theNormal = centreOfSecondSphere - centreOfSphere;
	theNormal.normalise();
	radiusOfSecondSphere = radiusOfAtom2 + probeRadius;
	
	// --------------------------------------- This calculates the centre ot cirlce as well as the radius of circle of an arc -----------//
	
	
	
	// first calculate the centre of the intersection circle - this starts at the centre of the central sphere and is parallel
	// to the connecting vector of the firstSphere and the scondSphere
	
	// calculate connecting vector centre sphere to centre circle is parallel to:
	
	centreToCircle = centreOfSecondSphere - centreOfSphere;
	
	// connecting vector has a radius of : (Rc*Rc + Rcs*Rcs - Rs*Rs)/2*Rcs 
	// where	Rc is the radius of the central atom
	//			Rs is the radius of the next atom in Sphere
	//			Rcs is the lenght of their distance vector <- dummy variables 
	
	double Rs,Rc,Rcs;					
	Rs = radiusOfSecondSphere;
	Rc = radiusOfSphere;
	Rcs = centreToCircle.get3DLength();
	
	// fisrt normalise
	centreToCircle /= Rcs;
	// then scale
	centreToCircle *= (Rc*Rc + Rcs*Rcs - Rs*Rs)/(2*Rcs);
	
//	cout << "[" << aHood->getAtomI()->element << "] " << (Rc*Rc + Rcs*Rcs - Rs*Rs)/(2*Rcs) << endl;
	//Now using the radius of the sphere and the intersection plane calculate radius of intersection circle
	radiusOfCircle = sqrt(radiusOfSphere*radiusOfSphere - centreToCircle.get3DLengthSq() );

	//Now precalculate centre of circle
	centreOfCircle = centreOfSphere + centreToCircle;

	unitCentreToCircle = centreToCircle;
	unitCentreToCircle.normalise();

}

CXXCircle::CXXCircle (CXXNewHood *aHood, const CXXCircleNode &aNode, double radiusOfAtom2, int aNumber) :
theAtomJ(0),
theParent(aHood),
centreOfSphere(theParent->getCentre()),
centreOfSecondSphere(aNode.getCoord()),
theNormal(centreOfSecondSphere-centreOfSphere),
radiusOfSecondSphere(radiusOfAtom2),
radiusOfSphere(aHood->getRadius()),
nodeNumber(aNumber)
{
	theNormal.normalise();
	
	// --------------------------------------- This calculates the centre ot cirlce as well as the radius of circle of an arc -----------//
	
	
	
	// first calculate the centre of the intersection circle - this starts at the centre of the central sphere and is parallel
	// to the connecting vector of the firstSphere and the scondSphere
	
	// calculate connecting vector centre sphere to centre circle is parallel to:
	
	centreToCircle = centreOfSecondSphere - centreOfSphere;
	
	// connecting vector has a radius of : (Rc*Rc + Rcs*Rcs - Rs*Rs)/2*Rcs 
	// where	Rc is the radius of the central atom
	//			Rs is the radius of the next atom in Sphere
	//			Rcs is the lenght of their distance vector <- dummy variables 
	
	double Rs,Rc,Rcs;					
	Rs = radiusOfSecondSphere;
	Rc = radiusOfSphere;
	Rcs = centreToCircle.get3DLength();
	
	// fisrt normalise
	centreToCircle.scale(1/Rcs);
	// then scale
	centreToCircle.scale((Rc*Rc + Rcs*Rcs - Rs*Rs)/(2*Rcs));
	
	//Now using the radius of the sphere and the intersection plane calculate radius of intersection circle
	radiusOfCircle = sqrt(radiusOfSphere*radiusOfSphere - centreToCircle.get3DLength()*centreToCircle.get3DLength() );
	
	//Now precalculate centre of circle
	centreOfCircle = centreOfSphere + centreToCircle;

	unitCentreToCircle = centreToCircle;
	unitCentreToCircle.normalise();

}

CXXCircle::CXXCircle (const CXXCircle &oldOne) :
theAtomJ(oldOne.getAtomJ()),
theParent(oldOne.getParent()),
centreOfSphere(oldOne.getCentreOfSphere()),
centreOfSecondSphere(oldOne.getCentreOfSecondSphere()),
theNormal(oldOne.getNormal()),
radiusOfSecondSphere(oldOne.getRadiusOfSecondSphere()),
radiusOfSphere(oldOne.getRadiusOfSphere()),
centreOfCircle(oldOne.centreOfCircle),
centreToCircle(oldOne.getCentreToCircle()),
unitCentreToCircle(oldOne.getUnitCentreToCircle()),
referenceUnitRadius(oldOne.getReferenceUnitRadius()),
radiusOfCircle(oldOne.getRadiusOfCircle()),
completelyEaten(oldOne.getEaten()),
nodeNumber(oldOne.getNodeNumber())
{
	theNodes.resize(oldOne.getNNodes());
	for (unsigned i=0; i<oldOne.getNNodes(); i++){
		theNodes[i] = oldOne.getNode(i);
		theNodes[i].setParent(this);
	}
	theStarts.resize(oldOne.nSegments());
	theStops.resize(oldOne.nSegments());
	for (unsigned i=0; i<oldOne.nSegments(); i++){
		theStarts[i] = oldOne.start(i);
		theStops[i] = oldOne.stop(i);
	}
}

int CXXCircle::meetsCircle(const CXXCircle &otherCircle, vector<CXXCoord> &nodeList) const{
		
	// check if there is an intersection between this circle and another circle
	// return 1 if there is intersection
	// 2 if the this circle is entirely swallowedd by the atom that generates the otherCircle, 
	// and 0 otherwise
	
	CXXCoord intersectA, intersectB;
	
	// some simplifications:	rj - centreToCircle of current (this) arc circle
	//							rk - centreToCircle of newArc circle
	
	CXXCoord rj(getCentreToCircle());
	CXXCoord rk(otherCircle.getCentreToCircle());
	CXXCoord rjxrk;
	
	// following Totrovs vector match - plane intersection of the two circles.
	// PlaneIntersectioon points at intersection of two this- and new- circlePlanes
	// with plane spanned by the centre of the two VDW spheres and the centre of thisArc
	
	// some dummy variables
	
	double rjrk, rjSquared, rkSquared;
	double a, b, distSq, radiusOfSphereSq, x;
	
	rjSquared = rj.dot(rj);
	rkSquared = rk.dot(rk);
	rjrk = rj.dot(rk);
	
	a = ((rjSquared - rjrk)*rkSquared)/(rjSquared*rkSquared - rjrk*rjrk);
	b = ((rkSquared - rjrk)*rjSquared)/(rjSquared*rkSquared - rjrk*rjrk);
	
	// with this now:
	CXXCoord planeIntersect((rj*a) + (rk*b));
		
	// if the distance between the intersection point of the three planes is closer then
	//the radius of the VDW + probe sphere, then the two circle of this arc and new arc intersect....
	
	// distance of intersectino point
	distSq = planeIntersect.dot(planeIntersect);			
	radiusOfSphereSq = radiusOfSphere*radiusOfSphere;
	
	
	// distance of (VDW + probe) sphere of thisAtom = radiusOfSphere therefore:
	if (distSq < radiusOfSphereSq) {
		
		// the CIRCLES corresponding to the new and the current (old) arc cross each other
		// the crossingpoints of the two CIRCLES are the points where the line defined by the 
		// intersection of the two arcplanes cuts throught the new circle 
		
		// the intersection line is perpendicular the plane spanned by rj and rk - that is parallel to:
		rjxrk = rj ^ rk;
		rjxrk.scale(1/rjxrk.get3DLength());
		
		
		
		//the intersection line is defined by this vector and one of its point - the planeIntersection vector:
		// => line planeIntersect + x*(rjxrk), x some real number
		// intersection points have x=+-0.5*length of the secante of circle therefore:
		x = sqrt(radiusOfSphereSq - distSq);
		
		rjxrk.scale(x);
		

		//We need to see whether the centreToCircle vectors of these circles indicate that the
		//corresponding torus is "behind" the centre of the first circle, or in front.  This determines
		//Whether the first or the second calculated intersection point represents a point to start drawing.
		//It boils down to a sort of parity thing.  Think about it !
		
	    int forwardOne = (theNormal.dot(unitCentreToCircle) > 0.?1:-1);
	    int forwardTwo = (otherCircle.getNormal().dot(otherCircle.getUnitCentreToCircle()) > 0.?1:-1);
		if (forwardOne*forwardTwo > 0){
			intersectA = planeIntersect+rjxrk;
			intersectB = planeIntersect-rjxrk;
		}
		else {
			intersectB = planeIntersect+rjxrk;
			intersectA = planeIntersect-rjxrk;
		}

	}
	else {
		//We get here if the two circles don't intersect
		//Now check whether this circle falls somewhere inside the
		//sphere indicated by the other circle.  If it does, then 
		//(given it doesn't intersect), it must *always* be inside
		// otherwise it must *always* be outside
		
		int oneEatenByTwo = 0;
		if (!getEaten()) {
			oneEatenByTwo = isSomewhereInsideSphere (otherCircle.getCentreOfSecondSphere(),
													 otherCircle.getRadiusOfSecondSphere());
		}
		// completelyEaten means circle1 is completely eaten by circle 2
		if (oneEatenByTwo)
			return 2;
		// otherwise circle2 may be completely eaten by circle 1
		int twoEatenByOne = 0;
		if (!otherCircle.getEaten()) {
			twoEatenByOne = otherCircle.isSomewhereInsideSphere(getCentreOfSecondSphere(),
																getRadiusOfSecondSphere());
		}
		if (twoEatenByOne) return 1;
		return 3;
	}
	
	nodeList[0] = intersectA + centreOfSphere;
	nodeList[1] = intersectB + centreOfSphere;
	
	return 0;
}

int CXXCircle::isSomewhereInsideSphere(const CXXCoord &centre, const double radius) const{
	//This can be checked by looking at any point on the circle rj and evaluating
	//if it is within the sphere
	CXXCoord xAxis(1.0, 0.0, 0.0);
	CXXCoord yAxis(0.0, 1.0, 0.0);
	
	//Need a vector orthogonal to centreToCircle.  Get this by crossing 
	CXXCoord radiusVector(0.,0.,0.);
	if (fabs(theNormal.dot(xAxis))<0.9999999) radiusVector = xAxis ^ theNormal;
	else radiusVector = yAxis ^ theNormal;
	radiusVector.normalise();
	radiusVector *= radiusOfCircle;
	CXXCoord pointOnCircle = getCentreOfCircle() + radiusVector;
	CXXCoord vectorToSphere = pointOnCircle - centre;
	double distanceToSphereSq = vectorToSphere.get3DLengthSq();
	if (distanceToSphereSq < radius*radius){
		return 1;
	}
	else {
		return 0;
	}
}

const CXXCoord &CXXCircle::getCentreOfSecondSphere() const{
	return centreOfSecondSphere;
}
const CXXCoord &CXXCircle::getCentreOfSphere() const{
	return centreOfSphere;
} 
const CXXCoord &CXXCircle::getCentreToCircle() const{
	return centreToCircle;
}
const CXXCoord &CXXCircle::getUnitCentreToCircle() const{
	return unitCentreToCircle;
}
const CXXCoord &CXXCircle::getCentreOfCircle() const{
	return centreOfCircle;
}
double CXXCircle::getRadiusOfSphere() const{
	return radiusOfSphere;
}
double CXXCircle::getRadiusOfSecondSphere() const{
	return radiusOfSecondSphere;
}
double CXXCircle::getRadiusOfCircle() const{
	return radiusOfCircle;
}
const PCAtom CXXCircle::getAtomJ() const{
	return theAtomJ;
}
const CXXCoord &CXXCircle::getNormal() const {
	return theNormal;
}

int CXXCircle::sortNodes(){

	if (theNodes.empty()) return 1;
	
	//Deal with the special case that this circle has two nodes introduced that correspond to
	//arbitrarily chosen start and end points on a complete circuit
	if (theNodes.size() == 2 &&
		theNodes[0].getFlag() == -1 &&
		theNodes[1].getFlag() == -1) {
		referenceUnitRadius = theNodes[0].getUnitRadius();
		theNodes[0].setAngle(0);
		theNodes[0].setFlag(2);
		theNodes[1].setAngle(2.*M_PI);
		theNodes[1].setFlag(1);
	}
	else {
		int iStartNode = -1;
		//March round till we find an "Entry" node (flag == 2)
		for (unsigned int iNode=0; iNode < theNodes.size() && iStartNode<0; iNode++){
			if (theNodes[iNode].getFlag() == 2) { 
				iStartNode = iNode;
			}
		}
		//For subsequent nodes, cause their corresponding angle around the circle to be
		//calculated with respect to this node
		referenceUnitRadius = theNodes[iStartNode].getUnitRadius();
		if (iStartNode==-1) cout << "Shouldn't happen in sortNodes" << endl;
		for (unsigned int iNode = 0; iNode< theNodes.size(); iNode++){
			if (iNode == (unsigned int) iStartNode) theNodes[iNode].setAngle(0.);
			else theNodes[iNode].setReference(referenceUnitRadius);
			//		cout << iNode << " " << *theNodes[iNode].getAngle() << " " << theNodes[iNode].getFlag() << endl;
		}
		
		
	}
	
	sort(theNodes.begin(), theNodes.end());
	
	//for (unsigned int iNode = 0; iNode< theNodes.size(); iNode++){
	//		cout << iNode << " " << *theNodes[iNode].getAngle() << " " << theNodes[iNode].getFlag() << endl;
	//}
	return 0;
}


int CXXCircle::identifyArcs() {
	
	int nNodes = theNodes.size();

	//No arcs to upload if completely eaten
	if (completelyEaten || nNodes==0) {
		theStarts.resize(0);
		theStops.resize(0);
		return 1;
	}
	

	//If this is an intact orbit, then atomK will be zero, and the first two nodes are trivially the start and
	//Stop points
	PCAtom atomK(theNodes[0].getAtomK());
	if (!atomK){
		theStarts.push_back(0);
		theStops.push_back(1);
		return 0;
	}

	//Prepare arrays of starts and corresponding ends modulated so that end omega is greater than start omega
	vector<double>nodesOne(nNodes/2);
	vector<double>nodesTwo(nNodes/2);


	unsigned int iCutter = 0;
	for (int i=0; i<nNodes; i++){
		const CXXCircleNode &theNode(theNodes[i]);
		if (theNode.getFlag() ==1){
			int iAtomK = theNode.getAtomK()->serNum;
			nodesOne[iCutter] = theNode.getAngle();
			int unMatched = 1;
			for ( int j=0; j<nNodes && unMatched; j++){
				const CXXCircleNode &theOtherNode(theNodes[j]);
				if (theOtherNode.getFlag() ==2 && theOtherNode.getAtomK()->serNum == iAtomK){
					double angle = theOtherNode.getAngle();
					while (angle < nodesOne[iCutter]) angle += M_PI*2.;
					nodesTwo[iCutter] = angle;
					unMatched = 0;
				}
			}
			iCutter++;
		}
	}
	
	int iStartNode, iEndNode;
	 int iNode = 0;
	
	while (iNode<nNodes){
		//March round until we are standing on a "2" and looking at a "1"
		while (iNode < nNodes && theNodes[iNode].getFlag() == 2 ) {
			iStartNode = iNode;
			iNode++;
		}
		iEndNode = iNode;
		
		if (iNode<nNodes) {
		//To decide how we proceed, we need to do the following:
		//Let w2_0 be the omega value associated with our node 2
		//let w1_0 be the omega value associated with our node1
		//increment w1_0 by 2_PI until w1_0> w2_0
		//This should only be drawn if for all other circles k
		// let w1k be omega associated with node 1
		// let w2k be omega associated with node 2
		// increment w2k by 2*pi until w2k > w1k
		// Now decrement w1k and w2k by 2*PI until w1k < w2_0
		// With this done, this segment is excluded if w2k> w1_0
			double w2_0 = theNodes[iStartNode].getAngle();
			double w1_0 = theNodes[iEndNode].getAngle();
			while (w1_0 < w2_0) w1_0 += 2.*M_PI;
			int doDraw = 1;
			for ( int k=0; k<nNodes/2 && doDraw; k++){
				double w1k = nodesOne[k];
				double w2k = nodesTwo[k];
				while (w1k > w2_0){
					w1k -= 2.*M_PI;
					w2k -= 2.*M_PI;
				}
				doDraw = (w1_0 > w2k);
			}
			
			if (doDraw){
				theStarts.push_back(iStartNode);
				theStops.push_back(iEndNode);
			//Flag that these nodes are suitable for including 
			//in the depiction of the reentrant surface
				theNodes[iStartNode].setDeleted(0);
				theNodes[iEndNode].setDeleted(0);
			}
			
			//Now advance round to the next node 2 and repeat
			while (iNode < nNodes && theNodes[iNode].getFlag() == 1) {
				iNode++;
			}
			iStartNode = iNode;
			
		}
		
	}
	return 0;
}

CXXNewHood *CXXCircle::getParent() const{
	return theParent;
}

unsigned CXXCircle::getNNodes() const{
	return theNodes.size();
}

const CXXCircleNode &CXXCircle::getNode(const int i) const{
	return theNodes[i];
}

int CXXCircle::addNode(const CXXCircleNode &aNode){
	theNodes.push_back(aNode);
	return 0;
}

int CXXCircle::getEaten() const{
	return completelyEaten;
}

void CXXCircle::setEaten(int flag){
	completelyEaten = flag;
}

unsigned CXXCircle::nSegments () const{
	return theStarts.size();
}

int CXXCircle::start(const int i) const {
	return theStarts[i];
}

int CXXCircle::stop(const int i) const{
	return theStops[i];
}

double CXXCircle::getRadiusOfVdWCircle() const{
	return getRadiusOfCircle() * 
	(getRadiusOfSphere()-theParent->getProbeRadius()) / getRadiusOfSphere();
}

const CXXCoord CXXCircle::getCentreOfVdWCircle() const{
	CXXCoord diff = getCentreToCircle();
	double frac = (getRadiusOfSphere()-theParent->getProbeRadius()) / getRadiusOfSphere();
	diff *= frac;
	return getCentreOfSphere()+diff;
}

void CXXCircle::dumpVdw() const{
	CXXCoord aRadius; 
	if (theNormal.dot(CXXCoord(1.,0.,0.)) < 0.9999999)
		aRadius = theNormal^CXXCoord(1.,0.,0.);
	else aRadius = theNormal^CXXCoord(0.,1.,0.);
	aRadius.normalise();
	CXXCoord bRadius = theNormal^aRadius;
	CXXCoord offset = theNormal;
	offset *= 0.001;
	CXXCoord cCent = getCentreOfVdWCircle()+offset;
	aRadius *= getRadiusOfVdWCircle();
	bRadius *= getRadiusOfVdWCircle();
	for (int i=0; i<12; i++){
		cout << "add triangle ";

		for (int j=0; j<3; j++){
			cout << cCent[j] << " ";
		}
		
		for (int j=0; j<3; j++){
			cout << cCent[j]
			+ cos(30.*M_PI/180.*double(i+1)) *aRadius[j]
			+ sin(30.*M_PI/180.*double(i+1)) *bRadius[j]
			<< " ";
		}	
		
		for (int j=0; j<3; j++){
			cout << cCent[j]
			+ cos(30.*M_PI/180.*double(i)) *aRadius[j]
			+ sin(30.*M_PI/180.*double(i)) *bRadius[j]
			<< " ";
		}		
		cout << endl;
	}
}

int CXXCircle::vdwIsBehind(const CXXCoord &crd) const{
	CXXCoord diff = crd - getCentreOfVdWCircle();
	return (diff.dot(getNormal()) < 0.);
}

const CXXCoord CXXCircle::vdwPlaneIntersect(const CXXCoord &A, const CXXCoord &B) const{
	CXXCoord AO = getCentreOfVdWCircle() - A;
	CXXCoord AB = B - A;
	double frac = AO.dot(getNormal()) / AB.dot(getNormal());
	AB *= frac;
	return A + AB;
}

int CXXCircle::accIsBehind(const CXXCoord &crd) const{
	CXXCoord diff(crd - getCentreOfCircle());
	return (diff.dot(getNormal()) < 0.);
}

const CXXCoord CXXCircle::accPlaneIntersect(const CXXCoord &A, const CXXCoord &B) const{
	CXXCoord AO(getCentreOfCircle() - A);
	CXXCoord AB(B - A);
	double frac = AO.dot(getNormal()) / AB.dot(getNormal());
	AB *= frac;
	CXXCoord result(A + AB);
	CXXCoord OC = result - getCentreOfCircle();
	double radialLength = OC.get3DLength();
	OC *= getRadiusOfCircle()/radialLength;
	return getCentreOfCircle() + OC;
}

const CXXCoord & CXXCircle::getReferenceUnitRadius() const{
	return referenceUnitRadius;
}

void CXXCircle::setArbitraryReference(){
	CXXCoord v1unit;
	CXXCoord xAxis(1.,0.,0.,0.);
	CXXCoord zAxis(0.,0.,1.,0.);
	if (fabs(getNormal().dot(zAxis)) < 0.9999){
		v1unit = (getNormal() ^ zAxis);
		v1unit.normalise();
	}
	else {
		v1unit = (getNormal() ^ xAxis);
		v1unit.normalise();
	}
	referenceUnitRadius=v1unit;
	v1unit *= getRadiusOfCircle();
	CXXCoord arbitraryEdgePoint(getCentreOfCircle() + v1unit);
	if (theParent){
		addNode(CXXCircleNode(this, 0, arbitraryEdgePoint, -1));
		addNode(CXXCircleNode(this, 0, arbitraryEdgePoint, -1));
	}
}
		
		
