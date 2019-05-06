/*
 *  CXXCircle.cpp
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include "CXXCircle.h"
#include "CXXNewHood.h"
#include <mmdb2/mmdb_manager.h>
#include <mmdb2/mmdb_tables.h>

CXX_mot::CXXCircle::CXXCircle() :
theAtomJ(0),
theBallJ(0),
theParent(0),
centreOfSecondSphere(CXXCoord(0.,0.,0.)),
theNormal(CXXCoord(0.,0.,0.)),
radiusOfSecondSphere(0.),
radiusOfSphere(0.),
centreOfCircle(CXXCoord(0.,0.,0.)),
centreToCircle(CXXCoord(0.,0.,0.)),
referenceUnitRadius(CXXCoord(0.,0.,0.)),
radiusOfCircle(0.),
nIntersectingCircles(0),
completelyEaten(0),
nodeNumber(0)
{
}

CXX_mot::CXXCircle::CXXCircle (CXXNewHood *aHood, mmdb::PAtom atom2, double radiusOfAtom2, double probeRadius) :
theAtomJ(atom2),
theBallJ(0),
theParent(aHood),
centreOfSecondSphere(CXXCoord(0.,0.,0.)),
theNormal(CXXCoord(0.,0.,0.)),
radiusOfSecondSphere(0.),
radiusOfSphere(aHood->getRadius()),
centreOfCircle(CXXCoord(0.,0.,0.)),
centreToCircle(CXXCoord(0.,0.,0.)),
referenceUnitRadius(CXXCoord(0.,0.,0.)),
radiusOfCircle(0.),
nIntersectingCircles(0),
completelyEaten(0),
nodeNumber(0)
{
	centreOfSecondSphere = CXXCoord(theAtomJ->x, theAtomJ->y, theAtomJ->z); 
	theNormal = centreOfSecondSphere - getCentreOfSphere();
	radiusOfSecondSphere = radiusOfAtom2 + probeRadius;
	
	performPrecalculations();
}

CXX_mot::CXXCircle::CXXCircle (CXXNewHood *aHood, const CXXBall &aBall) :
theAtomJ(aBall.getAtomI()),
theBallJ(&aBall),
theParent(aHood),
centreOfSecondSphere(aBall.getCoord()),
theNormal(CXXCoord(0.,0.,0.)),
radiusOfSecondSphere(aBall.getRadius()),
radiusOfSphere(aHood->getRadius()),
centreOfCircle(CXXCoord(0.,0.,0.)),
centreToCircle(CXXCoord(0.,0.,0.)),
referenceUnitRadius(CXXCoord(0.,0.,0.)),
radiusOfCircle(0.),
nIntersectingCircles(0),
completelyEaten(0),
nodeNumber(0)
{
    theNormal=centreOfSecondSphere-getCentreOfSphere();
	
    performPrecalculations();	
}

void CXX_mot::CXXCircle::performPrecalculations(){
    theNormal.normalise();
	
	// first calculate the centre of the intersection circle - this starts at the centre of the central sphere and is parallel
	// to the connecting vector of the firstSphere and the scondSphere
	
	// calculate connecting vector centre sphere to centre circle is parallel to:
	
	centreToCircle = centreOfSecondSphere - getCentreOfSphere();
	
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
	centreOfCircle = getCentreOfSphere() + centreToCircle;
	
	CXXCoord unitCentreToCircle = centreToCircle;
	unitCentreToCircle.normalise();
}
/*
CXX_mot::CXXCircle::CXXCircle (const CXXCircle &oldOne) :
theAtomJ(oldOne.getAtomJ()),
theBallJ(oldOne.getBallJ()),
theParent(oldOne.getParent()),
centreOfSecondSphere(oldOne.getCentreOfSecondSphere()),
theNormal(oldOne.getNormal()),
radiusOfSecondSphere(oldOne.getRadiusOfSecondSphere()),
radiusOfSphere(oldOne.getRadiusOfSphere()),
centreOfCircle(oldOne.centreOfCircle),
centreToCircle(oldOne.getCentreToCircle()),
referenceUnitRadius(oldOne.getReferenceUnitRadius()),
radiusOfCircle(oldOne.getRadiusOfCircle()),
theNodes(oldOne.getNodes()),
nIntersectingCircles(oldOne.getNIntersectingCircles()),
completelyEaten(oldOne.getEaten())
{
	for (unsigned i=0; i<oldOne.getNNodes(); i++){
		theNodes[i].setParent(this);
	}
	theStarts.resize(oldOne.nSegments());
	theStops.resize(oldOne.nSegments());
	for (unsigned i=0; i<oldOne.nSegments(); i++){
		theStarts[i] = oldOne.start(i);
		theStops[i] = oldOne.stop(i);
	}
}
*/
int CXX_mot::CXXCircle::meetsCircle(const CXXCircle &otherCircle, vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > &nodeList) const{
	
	// check if there is an intersection between this circle and another circle
	// return 0 if there is intersection
	// 2 if the this circle is entirely swallowedd by the atom that generates the otherCircle, 
	// 1 if the converse of the above applies
	// and 3 otherwise
	
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
	
	rjSquared = rj*rj;
	rkSquared = rk*rk;
	rjrk = rj*rk;
	
	a = ((rjSquared - rjrk)*rkSquared)/(rjSquared*rkSquared - rjrk*rjrk);
	b = ((rkSquared - rjrk)*rjSquared)/(rjSquared*rkSquared - rjrk*rjrk);
	
	// with this now:
	CXXCoord planeIntersect((rj*a) + (rk*b));
	
	// if the distance between the intersection point of the three planes is closer then
	//the radius of the VDW + probe sphere, then the two circle of this arc and new arc intersect....
	
	// distance of intersectino point
	distSq = planeIntersect*planeIntersect;			
	radiusOfSphereSq = radiusOfSphere*radiusOfSphere;
	
	// distance of (VDW + probe) sphere of thisAtom = radiusOfSphere therefore:
	if (distSq < radiusOfSphereSq) {
		
		// the CIRCLES corresponding to the new and the current (old) arc cross each other
		// the crossingpoints of the two CIRCLES are the points where the line defined by the 
		// intersection of the two arcplanes cuts throught the new circle 
		
		// the intersection line is perpendicular the plane spanned by rj and rk - that is parallel to:
		rjxrk = rj ^ rk;
		rjxrk.normalise();
		
		//the intersection line is defined by this vector and one of its point - the planeIntersection vector:
		// => line planeIntersect + x*(rjxrk), x some real number
		// intersection points have x=+-0.5*length of the secante of circle therefore:
		x = sqrt(radiusOfSphereSq - distSq);
		
		rjxrk.scale(x);
		
		//We need to see whether the centreToCircle vectors of these circles indicate that the
		//corresponding torus is "behind" the centre of the first circle, or in front.  This determines
		//Whether the first or the second calculated intersection point represents a point to start drawing.
		//It boils down to a sort of parity thing.  Think about it !
		
		CXXCoord unitCentreToCircle(centreToCircle);
		unitCentreToCircle.normalise();
	    int forwardOne = (theNormal*unitCentreToCircle > 0.?1:-1);
		CXXCoord otherCircleUnitCentreToCircle(otherCircle.getCentreToCircle());
		otherCircleUnitCentreToCircle.normalise();
	    int forwardTwo = (otherCircle.getNormal()*otherCircleUnitCentreToCircle > 0.?1:-1);
		if (forwardOne*forwardTwo > 0){
			intersectA = planeIntersect+rjxrk;
			intersectB = planeIntersect-rjxrk;
		}
		else {
			intersectB = planeIntersect+rjxrk;
			intersectA = planeIntersect-rjxrk;
		}
		nodeList[0] = intersectA + getCentreOfSphere();
		nodeList[1] = intersectB + getCentreOfSphere();
		
		return 0;
		
	}
	else {
		//We get here if the two circles don't intersect
		//Now check whether this circle falls somewhere inside the
		//sphere indicated by the other circle.  If it does, then 
		//(given it doesn't intersect), it must *always* be inside
		// otherwise it must *always* be outside
		
		if (otherCircle.accIsBehind(getCentreOfCircle()) == 0) return 2;
		if(accIsBehind(otherCircle.getCentreOfCircle()) == 0) return 1;
		return 3;
	}
	
}

int CXX_mot::CXXCircle::isSomewhereInsideSphere(const CXXCoord &centre, const double radius) const{
	//This can be checked by looking at any point on the circle rj and evaluating
	//if it is within the sphere
	CXXCoord xAxis(1.0, 0.0, 0.0);
	CXXCoord yAxis(0.0, 1.0, 0.0);
	
	//Need a vector orthogonal to centreToCircle.  Get this by crossing 
	CXXCoord radiusVector(0.,0.,0.);
	if (fabs(theNormal*xAxis)<0.9999999) radiusVector = xAxis ^ theNormal;
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

const CXX_mot::CXXCoord &CXX_mot::CXXCircle::getCentreOfSecondSphere() const{
	return centreOfSecondSphere;
}
const CXX_mot::CXXCoord &CXX_mot::CXXCircle::getCentreOfSphere() const{
    return theParent->getCentre();
};
const CXX_mot::CXXCoord &CXX_mot::CXXCircle::getCentreToCircle() const{
	return centreToCircle;
}
double CXX_mot::CXXCircle::getRadiusOfSphere() const{
	return radiusOfSphere;
}
double CXX_mot::CXXCircle::getRadiusOfSecondSphere() const{
	return radiusOfSecondSphere;
}
const mmdb::PAtom CXX_mot::CXXCircle::getAtomJ() const{
	return theAtomJ;
}

const CXX_mot::CXXCoord &CXX_mot::CXXCircle::getNormal() const {
	return theNormal;
}

int CXX_mot::CXXCircle::trimNodesBy(const CXXCircle &otherCircle)  {
    int deletedCount = 0;
    list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator lastNode = theNodes.end();
    for (list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator nodeIter = theNodes.begin();
         nodeIter != lastNode;
         ++nodeIter){
        CXXCircleNode &node(*nodeIter);
        if (node.getOtherCircle()!=this){
            if (!node.isDeleted() && node.getFlag()!=-1){
                if (!otherCircle.accIsBehind(node.getCoord())){
                    node.setDeleted(1);
                    deletedCount ++;
                }
            }
        }
    }
    return deletedCount;
}

bool CXX_mot::CXXCircle::abBracketsC(const CXXCircleNode &nodea, 
                            const CXXCircleNode &nodeb, 
                            const CXXCircleNode &nodec) const
{
    double aCrossBdotN = (nodea.getUnitRadius()^nodeb.getUnitRadius())*getNormal();
    double aCrossCdotN = (nodea.getUnitRadius()^nodec.getUnitRadius())*getNormal();
    double bCrossCdotN = (nodeb.getUnitRadius()^nodec.getUnitRadius())*getNormal();
    if (nodea.getUnitRadius().isNearly(nodec.getUnitRadius(), 0.000000000001) ||
        nodeb.getUnitRadius().isNearly(nodec.getUnitRadius(), 0.000000000001)){
        return false;
    }
    //    std::cout << "aC " << aCrossBdotN << " " << aCrossCdotN << " " << bCrossCdotN << std::endl;
    bool liesBetween = false;
    if (aCrossBdotN>0.){
        liesBetween = aCrossCdotN>0. && bCrossCdotN<0.;
    }
    else {
        liesBetween = aCrossCdotN>0. || bCrossCdotN<0.;
    }
    return liesBetween;
}

bool CXX_mot::CXXCircle::smallabBracketsC(const CXXCircleNode &nodea, 
                            const CXXCircleNode &nodeb, 
                            const CXXCircleNode &nodec) const
{
    double aCrossBdotN = (nodea.getUnitRadius()^nodeb.getUnitRadius())*getNormal();
    double aCrossCdotN = (nodea.getUnitRadius()^nodec.getUnitRadius())*getNormal();
    double bCrossCdotN = (nodeb.getUnitRadius()^nodec.getUnitRadius())*getNormal();
    bool liesBetween = false;
    if (aCrossBdotN>0.){
        liesBetween = aCrossCdotN>0. && bCrossCdotN<0.;
    }
    else {
        liesBetween = aCrossCdotN<0. && bCrossCdotN>0.;
    }
    return liesBetween;
}

void CXX_mot::CXXCircle::trimOwnNodes(){    
	
	//Delete all nodes that lie within any of the segments
    list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> > &otherCircles = theParent->getCircles();
	int nodesToDraw = countDrawnNodes();
	list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator otherCirclesEnd = otherCircles.end();
	for (list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator otherCircleIter = otherCircles.begin();
		 otherCircleIter != otherCirclesEnd && nodesToDraw > 0;
		 ++otherCircleIter){
		CXXCircle &otherCircle (*otherCircleIter);
		if (!otherCircle.getEaten()  &&
			&otherCircle != this){
			list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator nodesEnd = theNodes.end();
			for (list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator nodeIter = theNodes.begin();
				 nodeIter != nodesEnd && nodesToDraw>0;
				 ++nodeIter){
				CXXCircleNode &nodec(*nodeIter);
				if (nodec.getOtherCircle() != &otherCircle){
					if (nodec.isDeleted() == 0){
						if (!otherCircle.accIsBehind(nodec.getCoord())){
							nodec.setDeleted(1);
							nodesToDraw--;
						}
					}
				}
			}
		}
	}
	return;
}

int CXX_mot::CXXCircle::sortNodes(){
	if (theNodes.empty()) return 1;
	
	int initialNNodes = theNodes.size();
	//Deal with the special case that this circle has two nodes introduced that correspond to
	//arbitrarily chosen start and end points on a complete circuit
	if (initialNNodes == 2 &&
		theNodes.front().getFlag() == -1 &&
		theNodes.back().getFlag() == -1) {
		referenceUnitRadius = theNodes.begin()->getUnitRadius();
		theNodes.front().setAngle(0);
		theNodes.back().setFlag(2);
		theNodes.front().setAngle(2.*M_PI);
		theNodes.back().setFlag(1);
		return 0;
	}
    
    theNodes.remove_if(CXXCircleNode::shouldDelete);
	
    if (theNodes.size()%2) std::cout << "Seem to have non-even number of nodes (after trimming)\n";
	
	initialNNodes = theNodes.size();
    
    CXXCircleNode *startNode = 0;    
	//March round till we find an undeleted class2 node
    list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator nodesEnd = theNodes.end();
    for (list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator nodeIter = theNodes.begin();
         nodeIter != nodesEnd;
         ++nodeIter){
        CXXCircleNode &theNode(*nodeIter);
		if (theNode.isDeleted()==0) {
			if (theNode.getFlag()==2){
				startNode = &theNode;
			}
		}
	}
	if (startNode == 0) {
	   if (theNodes.size()>0) std::cout << theNodes.size() << "nodes but no startpoint\n";
	   theNodes.clear();
	   return 0;
	}	
	referenceUnitRadius = startNode->getUnitRadius();
    
    nodesEnd = theNodes.end();
    for (list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator nodeIter = theNodes.begin();
         nodeIter != nodesEnd;
         ++nodeIter){
        CXXCircleNode &theNode(*nodeIter);
		if (&theNode == startNode) {
			theNode.setAngle(0.);
		}
		else  theNode.setReference(referenceUnitRadius);
	}			
	
    theNodes.sort(CXXCircleNode::angleLessThan);
	
	return 0;
}

int CXX_mot::CXXCircle::newIdentifyArcs(){
	int nNodes = theNodes.size();
	
	//No arcs to upload if completely eaten
	if (getEaten() || nNodes==0) {
	   theStarts.clear();
	   theStops.clear();
	   return 1;
	}	
	
	//If this is an intact orbit, then atomK will be zero, and the first two nodes are trivially the start and
	//Stop points
	mmdb::PAtom atomK(theNodes.begin()->getAtomK());
	if (!atomK && theNodes.size()>1){
		theStarts.push_back(&(*theNodes.begin()));
		theStops.push_back(&theNodes.back());
		return 0;
	}
	
    list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator nodesEnd = theNodes.end();
    for (list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::iterator nodeIter = theNodes.begin();
         nodeIter != nodesEnd;
         ++nodeIter){
        CXXCircleNode &theNode(*nodeIter);

		//Tiptoe through the undeleted nodes
		if (theNode.isDeleted() == 0){
			if (theNode.getFlag() == 2){
				theStarts.push_back(&theNode);
			}
			else {
				theStops.push_back(&theNode);
			}
		}
	}
    if (theStarts.size() != theStops.size()) {
        std::cout << "uneven count of starts and stops\n";
        theStarts.clear();
    }
	return 0;
}

CXX_mot::CXXNewHood *CXX_mot::CXXCircle::getParent() const{
	return theParent;
}

unsigned CXX_mot::CXXCircle::getNNodes() const{
	return theNodes.size();
}
/*
CXXCircleNode &CXX_mot::CXXCircle::getNode(const int i) {
	return theNodes[i];
}

const CXXCircleNode &CXX_mot::CXXCircle::getNode(const int i) const{
	return theNodes[i];
}
*/
int CXX_mot::CXXCircle::addNode(const CXXCircleNode &aNode){
	theNodes.push_back(aNode);
	return 0;
}

int CXX_mot::CXXCircle::getEaten() const{
	return completelyEaten;
}

void CXX_mot::CXXCircle::setEaten(int flag){
	completelyEaten = flag;
    if (flag) theNodes.clear();
}

unsigned CXX_mot::CXXCircle::nSegments () const{
	return theStarts.size();
}

CXX_mot::CXXCircleNode* CXX_mot::CXXCircle::start(const int i) const {
	return theStarts[i];
}

CXX_mot::CXXCircleNode* CXX_mot::CXXCircle::stop(const int i) const{
	return theStops[i];
}

double CXX_mot::CXXCircle::getRadiusOfVdWCircle() const{
	return getRadiusOfCircle() * 
	(getRadiusOfSphere()-theParent->getProbeRadius()) / getRadiusOfSphere();
}

const CXX_mot::CXXCoord CXX_mot::CXXCircle::getCentreOfVdWCircle() const{
	CXXCoord diff = getCentreToCircle();
	double frac = (getRadiusOfSphere()-theParent->getProbeRadius()) / getRadiusOfSphere();
	diff *= frac;
	return getCentreOfSphere()+diff;
}

void CXX_mot::CXXCircle::dumpVdw() const{
	CXXCoord aRadius; 
	if (theNormal*CXXCoord(1.,0.,0.) < 0.9999999)
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

int CXX_mot::CXXCircle::vdwIsBehind(const CXXCoord &crd) const{
	CXXCoord diff = crd - getCentreOfVdWCircle();
	return (diff*getNormal() < 0.);
}

const CXX_mot::CXXCoord CXX_mot::CXXCircle::vdwPlaneIntersect(const CXXCoord &A, const CXXCoord &B) const{
	CXXCoord AO = getCentreOfVdWCircle() - A;
	CXXCoord AB = B - A;
	double frac = (AO*getNormal()) / (AB*getNormal());
	AB *= frac;
	return A + AB;
}

const CXX_mot::CXXCoord CXX_mot::CXXCircle::accPlaneIntersect(const CXXCoord &A, const CXXCoord &B) const{
	CXXCoord AO(getCentreOfCircle() - A);
	CXXCoord AB(B - A);
	double frac = (AO*getNormal()) / (AB*getNormal());
	AB *= frac;
	CXXCoord result(A + AB);
	CXXCoord OC = result - getCentreOfCircle();
	double radialLength = OC.get3DLength();
	OC *= getRadiusOfCircle()/radialLength;
	return getCentreOfCircle() + OC;
}

const CXX_mot::CXXCoord & CXX_mot::CXXCircle::getReferenceUnitRadius() const{
	return referenceUnitRadius;
}

void CXX_mot::CXXCircle::setArbitraryReference(){
	CXXCoord v1unit;
	CXXCoord xAxis(1.,0.,0.,0.);
	CXXCoord zAxis(0.,0.,1.,0.);
	if (fabs(getNormal()*zAxis) < 0.9999){
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
    addNode(CXXCircleNode(this, 0, arbitraryEdgePoint, -1));
    addNode(CXXCircleNode(this, 0, arbitraryEdgePoint, -1));
}

int CXX_mot::CXXCircle::countDrawnNodes() const
{
    int answer = 0;
    list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::const_iterator nodesEnd = theNodes.end();
    for (list<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >::const_iterator nodeIter = theNodes.begin();
         nodeIter != nodesEnd;
         ++nodeIter){
        const CXXCircleNode &theNode(*nodeIter);
        
        if (theNode.isDeleted()==0) answer++;
    }
    return answer;
}


