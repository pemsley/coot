/*
 *  CXXNewHood.cpp
 *  CXXSurface
 *
 *  Created by martin on Sat Feb 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <sstream>
#include <cstring>

#include "CXXSurface.h"
#include "CXXNewHood.h"
#include <mmdb2/mmdb_manager.h>
#include <mmdb2/mmdb_tables.h>
#include "CXXCircle.h"
#include "CXXCircleNode.h"
#include "CXXTorusElement.h"

void CXX_mot::CXXNewHood::init(){
	theAtomI = 0;
	theRadius = 0;
    theProbeRadius = 0;
	theCentre = CXXCoord(0.,0.,0.);
	theCircles.clear();
}

CXX_mot::CXXNewHood::CXXNewHood(){
	init();
}

CXX_mot::CXXNewHood::CXXNewHood(mmdb::PAtom centralAtom, double radiusOfAtom1, double probeRadius) :
theAtomI(centralAtom), 
theRadius(radiusOfAtom1+probeRadius), 
theProbeRadius(probeRadius){
	theCentre = CXXCoord(theAtomI->x, theAtomI->y, theAtomI->z);
}

void CXX_mot::CXXNewHood::initWith(const CXXCircleNode &aNode, double probeRadius) {
    theAtomI= 0;
    theRadius = probeRadius;
    theProbeRadius = 0;
    theCentre = aNode.getCoord();
}

void CXX_mot::CXXNewHood::initWith(const CXXBall *aBall){
	theBall = aBall;
	theRadius = aBall->getRadius();
	theProbeRadius = 0.;
	theCentre = aBall->getCoord();
};

int CXX_mot::CXXNewHood::addAtom(mmdb::PAtom anAtomJ, double radiusOfAtom2){
	if (anAtomJ->serNum == theAtomI->serNum) {
		//		std::cout << "Rejecting self " << anAtomJ->serNum << " " <<  theAtomI->serNum << endl;
		return 0; //Worried this might not be unique
	}
	double radiusOfAtomJ = radiusOfAtom2 + theProbeRadius;
	CXXCoord centreOfAtomJ = CXXCoord(anAtomJ->x, anAtomJ->y, anAtomJ->z);
	CXXCoord centreOfAtomI = theCentre;
	
	//This test is to deal with unlikely but possible event of being fed asphere with identical coordinates
	int result = 0;
	if (!centreOfAtomI.isNearly(centreOfAtomJ, 0.0001)){
        CXXCoord separ = centreOfAtomJ - centreOfAtomI;
        double sumrad = radiusOfAtomJ + theRadius;
        if (separ.get3DLengthSq()<(sumrad*sumrad)){
			CXXCircle aCircle(this, anAtomJ, radiusOfAtom2, theProbeRadius);
            theCircles.push_back(aCircle);
            result= 1;
        }
	}
	
	return result;
}

int CXX_mot::CXXNewHood::addBall(const CXXBall &aBall){
    if (!theCentre.isNearly(aBall.getCoord(), 0.0001)){
        CXXCoord separ = aBall.getCoord() - theCentre;
        double sumrad = theRadius + aBall.getRadius();
        if (separ.get3DLengthSq()<(sumrad*sumrad)){
            theCircles.push_back(CXXCircle(this, aBall));
            return 1;
        }
    }
    return 0;
}

int CXX_mot::CXXNewHood::findSegments(){
	vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > theIntersections(2);
    CXXCircleNode circleNode0;
    CXXCircleNode circleNode1;
    CXXCircleNode circleNode2;
    CXXCircleNode circleNode3;
	//First cause circles to find their nodes
	
	std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circleIter = theCircles.begin();
		 circleIter != circlesEnd;
		 ++circleIter){
		CXXCircle &circle1(*circleIter);
		if (!circle1.getEaten()){
			for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circle2Iter = circleIter;
				 circle2Iter != circlesEnd;
				 ++circle2Iter){
				CXXCircle &circle2(*circle2Iter);
				if (&circle1 != &circle2){
					int eatenBy=circle1.meetsCircle(circle2, theIntersections);
					if (eatenBy == 2){
						circle1.setEaten(1);
					}
					else if (eatenBy == 1){
						circle2.setEaten(1);
					}
					else if (eatenBy==0){
                        circleNode0.setParent(&circle1);
                        circleNode0.setOtherCircle(&circle2);
                        circleNode0.setCoord(theIntersections[0]);
                        circleNode0.setFlag(1);
						circle1.addNode(circleNode0);
                        /*CXXCircleNode(&circle1,
						 &circle2,
						 theIntersections[0],
						 1));*/
                        circleNode1.setParent(&circle1);
                        circleNode1.setOtherCircle(&circle2);
                        circleNode1.setCoord(theIntersections[1]);
                        circleNode1.setFlag(2);
						circle1.addNode(circleNode1);
                        /*CXXCircleNode(&circle1,
						 &circle2,
						 theIntersections[1],
						 2)); */
                        
                        circleNode2.setParent(&circle2);
                        circleNode2.setOtherCircle(&circle1);
                        circleNode2.setCoord(theIntersections[0]);
                        circleNode2.setFlag(2);
						circle2.addNode(circleNode2);
                        /*CXXCircleNode(&circle2,
						 &circle1,
						 theIntersections[0],
						 2));*/
                        
                        circleNode3.setParent(&circle2);
                        circleNode3.setOtherCircle(&circle1);
                        circleNode3.setCoord(theIntersections[1]);
                        circleNode3.setFlag(1);
						circle2.addNode(circleNode3);
                        /*CXXCircleNode(&circle2,
						 &circle1,
						 theIntersections[1],
						 1)); */
					}
				}
			}
		}
    }
	circlesEnd = theCircles.end();
	for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circleIter = theCircles.begin();
		 circleIter != circlesEnd;
		 ++circleIter){
		CXXCircle &circle1(*circleIter);
        bool completeOrbit = circle1.getNNodes() == 0;
        //Don't bother with this bit if  the circle has been eaten
		if (!circle1.getEaten()){
            // Here handle case where a circle has simply not intersected with another circle
            // and so is an intact 2*PI circle
            if (completeOrbit){
                circle1.setArbitraryReference();
            }		
            circle1.trimOwnNodes();
 			// Make the circle sort its nodes
			circle1.sortNodes();
			// Then make the circles upload their unique arcs to us
			circle1.newIdentifyArcs();
		}
	}
	return 0;
}

const CXX_mot::CXXCoord &CXX_mot::CXXNewHood::getCentre() const{
	return theCentre;
}

double CXX_mot::CXXNewHood::getRadius() const{
	return theRadius;
}
int CXX_mot::CXXNewHood::nCircles() const{
	return theCircles.size();
}
const mmdb::PAtom CXX_mot::CXXNewHood::getAtomI() const {
	return theAtomI;
}

double CXX_mot::CXXNewHood::getProbeRadius() const{
	return theProbeRadius;
}


bool CXX_mot::CXXNewHood::doesNotContainDrawable(const CXXNewHood &aHood) {
	//Identify how many of the circles of this hood actually have active segments
	//and hence surface elements to draw...no nodes means either no intersecting
	//spheres, or no accessible surfce patches 

	int nActiveCircles = 0;
	const std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >& theCircles(aHood.getCircles());
	std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circlesEnd = theCircles.end();
	bool result = (theCircles.size() == 0);
	for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circleIter = theCircles.begin();
		 circleIter != circlesEnd;
		 ++circleIter){
		const CXXCircle &theCircle(*circleIter);
		
		if (theCircle.nSegments()!=0){
			result = true;
			nActiveCircles++;
		}
	}
	return !result;
}

void CXX_mot::CXXNewHood::triangulateAsRegularHoodInto(CXXSurface *aSurface, double delta, const CXXSphereElement *unitSphereAtOrigin) const
{    
    const CXXNewHood &theNewHood(*this);
    
    //For now, copy the unit sphere to atom position, and translate and scale to become a vdw sphere
    CXXSphereElement vdwSphere(*unitSphereAtOrigin);
    vdwSphere.scaleBy(theNewHood.getRadius());		
    vdwSphere.translateBy (theNewHood.getCentre());
    vdwSphere.setAtom(theNewHood.getAtomI());
    
    //Use the patch edges to trim the triangles of the VDW sphere
	
	std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circleIter = theCircles.begin();
		 circleIter != circlesEnd && vdwSphere.getNDrawnTriangles();
		 ++circleIter){
		const CXXCircle &theCircle(*circleIter);        
        //All circles are needed to eat into the vdw sphere (except those that
        //were wholely subsumed by another circle
        if (!theCircle.getEaten()) {
            vdwSphere.trimBy(theCircle, 1);				
        }
    }                
    
    //We collect the sphere vertices that have been generated in the trimming process
    //so that we can introduce them into the corresponding torus later
    map<const CXXCircle *, vector<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >  > rawEdges;
    if (vdwSphere.nFlatTriangles()) {
        for (unsigned iVertex = 0; iVertex < vdwSphere.nVertices(); iVertex++){
            if (vdwSphere.vertex(iVertex).doDraw()){
                const CXXCircle *correspondingCircle = vdwSphere.vertex(iVertex).getIntersector();
                if (correspondingCircle != 0){
                    CXXCircleNode extraNode(correspondingCircle,  0, vdwSphere.vertex(iVertex).vertex(),iVertex);
                    extraNode.setReference(correspondingCircle->getReferenceUnitRadius());					
                    rawEdges[correspondingCircle].push_back(extraNode);					
                }
            }
        }
    }
    
    //Create torus elements for each of the segments, and elaborate them with the additional vertices
    //that correspond to ragged edges from the sphere
	circlesEnd = theCircles.end();
	for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circleIter = theCircles.begin();
		 circleIter != circlesEnd;
		 ++circleIter){
		const CXXCircle &theCircle(*circleIter);        
        if ((!theCircle.getEaten()) && theCircle.getNNodes() != 0){
            //Here cause the sphere to identify triangles that lie on each circle, and assign to them
            //"absolute" omega values:  this will be used later when addo=ing torus vertices
            vdwSphere.flagCutTriangles(theCircle);
            
            //Loop over segments of the circle, generating torus elements, and adding sphere nodes where
            //appropriate
            
            for (unsigned iTorus = 0; iTorus < theCircle.nSegments(); iTorus++){
                CXXTorusElement theTorus (theCircle, iTorus, delta, getProbeRadius());
                for (unsigned iRawEdge=0; iRawEdge < rawEdges[&theCircle].size(); iRawEdge++){
                    theTorus.addEdgeVertex(rawEdges[&theCircle][iRawEdge]);
                }
                theTorus.upload(aSurface);  
                //Here cause the sphere to subdivide the triangles that are interrupted by
                //nodes around this torus
                
                vdwSphere.addTorusVertices(theTorus);
            }	
        } 
    }        
    //Copy surface triangles from this sphere into the surface
    aSurface->upLoadSphere(vdwSphere, getProbeRadius(), CXXSphereElement::Contact);        
}

void CXX_mot::CXXNewHood::identifyUniqueNodes(vector<CXXCircleNode, CXX_old::CXXAlloc<CXXCircleNode> >&circleNodes, int selHnd) const {
    //Generate the equivalent hoods
	std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circleIter = theCircles.begin();
		 circleIter != circlesEnd;
		 ++circleIter){
		const CXXCircle &theCircle(*circleIter);
        if (!theCircle.getEaten()){
            //Unique Nodes will be segment starts, and should have a "otherCircle",
            //to distinguish them from the nodes added to a dummy orbit
            for (unsigned iSegment = 0; iSegment<theCircle.nSegments(); iSegment++){
                const CXXCircleNode *ends[2];
                ends[0] = theCircle.start(iSegment);
                ends[1] = theCircle.stop(iSegment);
                for (int iTerminus = 0; iTerminus<2; iTerminus++){
                    const CXXCircleNode &aNode(*ends[iTerminus]);
                    if (aNode.getOtherCircle()){
						bool includeNode = false;
                        //First easy option...all three nodes in selection...we should only include one of the
                        //three occurrences of this 
                        if (aNode.getAtomJ()->isInSelection(selHnd) &&
                            aNode.getAtomK()->isInSelection(selHnd)) {
                            if ( aNode.getAtomI() < aNode.getAtomJ() &&
                                aNode.getAtomJ() < aNode.getAtomK() ){
								includeNode = true;
                            }
                        }
                        //Next easy case: both j and k are not in selection...we 
                        //will only encounter this node once
                        else if (!aNode.getAtomJ()->isInSelection(selHnd) &&
                                 !aNode.getAtomK()->isInSelection(selHnd)) {
							if (aNode.getAtomJ() < aNode.getAtomK()) includeNode = true;
                        }
                        //Rolling between two atoms in the selection, towards a third atom 
                        //*not* in the selection
                        //This can happen in two possible ways, of which we accept only 1
                        else if (aNode.getAtomJ()->isInSelection(selHnd) &&
                                 !aNode.getAtomK()->isInSelection(selHnd)) {
                            if ( aNode.getAtomI() < aNode.getAtomJ() ){
								includeNode = true;
                            }
                        } 
                        //Rolling between one atom in selection and one atom not in, towards a third atom
                        //that is in the selection.  Every case where
                        //this can happen has an equivalent case above from which we keep the node
                        else if (aNode.getAtomK()->isInSelection(selHnd)) {
                        }
						if (includeNode) {
#pragma omp critical (circleNodes)
							circleNodes.push_back(aNode);	
						}
                    }
				}
            }
        }
    }
}

void CXX_mot::CXXNewHood::triangulateAsBallHoodInto(CXXSurface *aSurface, double delta, 
										   std::map<const CXXBall*, std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > > &raggedEdges, 
										   bool useEdges, int insideOrOutside) 
{    
	CXXSphereElement ballSphere;
	theBall->initSphereElement(ballSphere, delta);
	//Now use the circles identified in the above hood to trim the ball sphere element
	int wasTrimmed = 0;
	std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::const_iterator circleIter = theCircles.begin();
		 circleIter != circlesEnd && ballSphere.getNDrawnTriangles ();
		 ++circleIter){
		const CXXCircle &theCircle(*circleIter);		//our ball hood will contain circles that correspond to the initial edges of the
		//triangular patches: since they have no ballJ, we can exclude them from he trimming process
		if (theCircle.getBallJ()!=0) {
			if (ballSphere.trimBy(theCircle, 1)==1) {
				wasTrimmed = 1;
			}
		}
	}
	if (!useEdges && !wasTrimmed){
#pragma omp critical (mainTriangles)
		aSurface->upLoadSphere(ballSphere, theBall->getRadius(), insideOrOutside);
	}
	else if (!useEdges && wasTrimmed && ballSphere.getNDrawnTriangles()>0 ) {
		//We get here if this sphere was trimmed by intersection with other spheres...
		//this is the point at which to distribute the corresponding loose ends amongst
		//the probes that generated them
		std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circlesEnd = theCircles.end();
		for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circleIter = theCircles.begin();
			 circleIter != circlesEnd;
			 ++circleIter){
			CXXCircle &theCircle(*circleIter);
			if (theCircle.getBallJ()!=0) {
				ballSphere.identifyRaggedEdges(theCircle, raggedEdges);
			}
		}	
	}
	else if (useEdges && wasTrimmed && ballSphere.getNDrawnTriangles()>0 ) {
		std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circlesEnd = theCircles.end();
		for (std::list<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >::iterator circleIter = theCircles.begin();
			 circleIter != circlesEnd && ballSphere.getNDrawnTriangles()>0;
			 ++circleIter){
			CXXCircle &theCircle(*circleIter);
			const CXXBall *ballJ(theCircle.getBallJ());
			ballSphere.flagCutTriangles(theCircle);
			std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > &ballsEdges(raggedEdges[ballJ]);
			std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> >::iterator ballsEdgesEnd(ballsEdges.end());
			for (std::vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> >::iterator ballsEdge = ballsEdges.begin();
				 ballsEdge!= ballsEdges.end();
				 ballsEdge++){
				ballSphere.addVertex(CXXCircleNode(&theCircle, 0, *ballsEdge, 0));
			}
		}
#pragma omp critical (mainTriangles)
		aSurface->upLoadSphere(ballSphere, theBall->getRadius(), insideOrOutside);
	}
}

