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
#include "CXXSurfaceVertex.h"
#include "CXXSurface.h"
#include "CXXNewHood.h"
#include "mmdb2/mmdb_manager.h"
#include "mmdb2/mmdb_tables.h"
#include "CXXCircle.h"
#include "CXXCircleNode.h"
#include "CXXTorusElement.h"

void CXXNewHood::init(){
	theAtomI = 0;
	theRadius = 0;
    theProbeRadius = 0;
	theCentre = CXXCoord<CXXCoord_ftype>(0.,0.,0.);
	theCircles.resize(0);
}

CXXNewHood::CXXNewHood(){
	init();
}

CXXNewHood::CXXNewHood(mmdb::Atom* centralAtom, double radiusOfAtom1, double probeRadius) :
theAtomI(centralAtom), 
theRadius(radiusOfAtom1+probeRadius), 
theProbeRadius(probeRadius){
	theCentre = CXXCoord<CXXCoord_ftype>(theAtomI->x, theAtomI->y, theAtomI->z);
}

void CXXNewHood::initWith(const CXXCircleNode &aNode, double probeRadius) {
    theAtomI= 0;
    theRadius = probeRadius;
    theProbeRadius = 0;
    theCentre = aNode.getCoord();
}

void CXXNewHood::initWith(const CXXBall *aBall){
	theBall = aBall;
	theRadius = aBall->getRadius();
	theProbeRadius = 0.;
	theCentre = aBall->getCoord();
};

int CXXNewHood::addAtom(mmdb::Atom* anAtomJ, double radiusOfAtom2){
	if (anAtomJ->serNum == theAtomI->serNum) {
		//		std::cout << "Rejecting self " << anAtomJ->serNum << " " <<  theAtomI->serNum << endl;
		return 0; //Worried this might not be unique
	}
	double radiusOfAtomJ = radiusOfAtom2 + theProbeRadius;
	CXXCoord<CXXCoord_ftype>centreOfAtomJ = CXXCoord<CXXCoord_ftype>(anAtomJ->x, anAtomJ->y, anAtomJ->z);
	CXXCoord<CXXCoord_ftype>centreOfAtomI = theCentre;
	
	//This test is to deal with unlikely but possible event of being fed asphere with identical coordinates
	int result = 0;
	if (!centreOfAtomI.isNearly(centreOfAtomJ, 0.0001)){
        CXXCoord<CXXCoord_ftype>separ = centreOfAtomJ - centreOfAtomI;
        double sumrad = radiusOfAtomJ + theRadius;
        if (separ.get3DLengthSq()<(sumrad*sumrad)){
			CXXCircle aCircle(this, anAtomJ, radiusOfAtom2, theProbeRadius);
            theCircles.push_back(aCircle);
            result= 1;
        }
	}
	
	return result;
}

int CXXNewHood::addBall(const CXXBall &aBall){
    if (!theCentre.isNearly(aBall.getCoord(), 0.0001)){
        CXXCoord<CXXCoord_ftype>separ = aBall.getCoord() - theCentre;
        double sumrad = theRadius + aBall.getRadius();
        if (separ.get3DLengthSq()<(sumrad*sumrad)){
            theCircles.push_back(CXXCircle(this, aBall));
            return 1;
        }
    }
    return 0;
}

int CXXNewHood::findSegments(){
	vector<CXXCoord<CXXCoord_ftype> > theIntersections(2);
    CXXCircleNode circleNode;
	//First cause circles to find their nodes
	
	std::list<CXXCircle  >::iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle  >::iterator circle1Iter = theCircles.begin();
		 circle1Iter != circlesEnd;
		 ++circle1Iter){
		CXXCircle &circle1(*circle1Iter);
		if (!circle1.getEaten()){
			for (std::list<CXXCircle  >::iterator circle2Iter = circle1Iter;
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
                        
                        //Here I am going to add the corresponding nodes *only* if they are not nixed by another circle
						
                        for (int i=0; i<2; i++){
                            bool addNodes = true;
                            for (std::list<CXXCircle>::iterator circle3Iter = theCircles.begin();
                                 circle3Iter != circlesEnd && addNodes;
                                 ++circle3Iter){
                                if (circle3Iter != circle1Iter && circle3Iter != circle2Iter){
                                    CXXCircle &circle3(*circle3Iter);
                                    if (!circle3.getEaten()){
                                        if (!circle3.accIsBehind(theIntersections[i])) addNodes = false;
                                    }
                                }
                            }
                            if (addNodes) {
                                circleNode.setParent(&circle1);
                                circleNode.setOtherCircle(&circle2);
                                circleNode.setCoord(theIntersections[i]);
                                circleNode.setFlag((i==0?1:2));
                                circle1.addNode(circleNode);
                                circleNode.setParent(&circle2);
                                circleNode.setOtherCircle(&circle1);
                                circleNode.setCoord(theIntersections[i]);
                                circleNode.setFlag((i==0?2:1));
                                circle2.addNode(circleNode);
                                
                            }
                            else {
                                circle1.setContainsEatenNodes(1);   
                                circle2.setContainsEatenNodes(1);
                            }
                        }
                    }
                }
            }
        }
    }
    circlesEnd = theCircles.end();
    for (std::list<CXXCircle  >::iterator circleIter = theCircles.begin();
         circleIter != circlesEnd;
         ++circleIter){
        CXXCircle &circle1(*circleIter);
        bool completeOrbit = (circle1.getNNodes() == 0 && circle1.getContainsEatenNodes()==0);
        //Don't bother with this bit if  the circle has been eaten
        if (!circle1.getEaten()){
            // Here handle case where a circle has simply not intersected with another circle
            // and so is an intact 2*PI circle
            if (completeOrbit){
                circle1.setArbitraryReference();
            }		
            // Make the circle sort its nodes
            circle1.sortNodes();
            // Then make the circles upload their unique arcs to us
            circle1.newIdentifyArcs();
        }
    }
    return 0;
}

const CXXCoord<CXXCoord_ftype>&CXXNewHood::getCentre() const{
	return theCentre;
}

double CXXNewHood::getRadius() const{
	return theRadius;
}
size_t CXXNewHood::nCircles() const{
	return theCircles.size();
}
mmdb::Atom* CXXNewHood::getAtomI() const {
	return theAtomI;
}

double CXXNewHood::getProbeRadius() const{
	return theProbeRadius;
}


bool CXXNewHood::containsDrawable(const CXXNewHood &aHood) {
	//Identify how many of the circles of this hood actually have active segments
	//and hence surface elements to draw...no nodes means either no intersecting
	//spheres, or no accessible surfce patches
	const std::list<CXXCircle  >& theCircles(aHood.getCircles());
    bool result = (theCircles.size() == 0);
	std::list<CXXCircle  >::const_iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle  >::const_iterator circleIter = theCircles.begin();
		 circleIter != circlesEnd && !result;
		 ++circleIter){
		const CXXCircle &theCircle(*circleIter);
		
		if (theCircle.nSegments()!=0){
			result = true;
		}
	}
	return result;
}

void CXXNewHood::triangulateAsRegularHoodInto(CXXSurface &aSurface, double delta, const CXXSphereElement *unitSphereAtOrigin) const
{    
    const CXXNewHood &theNewHood(*this);
    
    //For now, copy the unit sphere to atom position, and translate and scale to become a vdw sphere
    CXXSphereElement vdwSphere(*unitSphereAtOrigin);
    vdwSphere.scaleBy(theNewHood.getRadius());		
    vdwSphere.translateBy (theNewHood.getCentre());
    vdwSphere.setAtom(theNewHood.getAtomI());
    
    //Use the patch edges to trim the triangles of the VDW sphere
	
	std::list<CXXCircle  >::const_iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle  >::const_iterator circleIter = theCircles.begin();
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
    map<const CXXCircle *, vector<CXXCircleNode  >  > rawEdges;
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
	for (std::list<CXXCircle  >::const_iterator circleIter = theCircles.begin();
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
                aSurface.uploadTorus(theTorus);
                //                theTorus.upload(aSurface);  
                //Here cause the sphere to subdivide the triangles that are interrupted by
                //nodes around this torus
                
                vdwSphere.addTorusVertices(theTorus);
            }	
        } 
    }        
    //Copy surface triangles from this sphere into the surface
    aSurface.upLoadSphere(vdwSphere, getProbeRadius(), CXXSphereElement::Contact);
}

void CXXNewHood::identifyUniqueNodes(vector<CXXCircleNode  >&circleNodes, int selHnd) const {
    //Generate the equivalent hoods
	std::list<CXXCircle  >::const_iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle  >::const_iterator circleIter = theCircles.begin();
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

void CXXNewHood::triangulateAsBallHoodInto(CXXSurface &aSurface, double delta,
										   std::map<const CXXBall*, std::vector<CXXCoord<CXXCoord_ftype> > > &raggedEdges,
										   bool useEdges, int insideOrOutside,
                                           const CXXSphereElement &unitCellAtOriginForDelta)
{    
    CXXSphereElement ballSphere;
    theBall->initSphereElement(ballSphere, delta, unitCellAtOriginForDelta);
	//Now use the circles identified in the above hood to trim the ball sphere element
	int wasTrimmed = 0;
	std::list<CXXCircle  >::const_iterator circlesEnd = theCircles.end();
	for (std::list<CXXCircle  >::const_iterator circleIter = theCircles.begin();
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
		aSurface.upLoadSphere(ballSphere, theBall->getRadius(), insideOrOutside);
	}
	else if (!useEdges && wasTrimmed && ballSphere.getNDrawnTriangles()>0 ) {
		//We get here if this sphere was trimmed by intersection with other spheres...
		//this is the point at which to distribute the corresponding loose ends amongst
		//the probes that generated them
		std::list<CXXCircle  >::iterator circlesEnd = theCircles.end();
		for (std::list<CXXCircle  >::iterator circleIter = theCircles.begin();
			 circleIter != circlesEnd;
			 ++circleIter){
			CXXCircle &theCircle(*circleIter);
			if (theCircle.getBallJ()!=0) {
				ballSphere.identifyRaggedEdges(theCircle, raggedEdges);
			}
		}	
	}
	else if (useEdges && wasTrimmed && ballSphere.getNDrawnTriangles()>0 ) {
		std::list<CXXCircle  >::iterator circlesEnd = theCircles.end();
		for (std::list<CXXCircle  >::iterator circleIter = theCircles.begin();
			 circleIter != circlesEnd && ballSphere.getNDrawnTriangles()>0;
			 ++circleIter){
			CXXCircle &theCircle(*circleIter);
			const CXXBall *ballJ(theCircle.getBallJ());
			ballSphere.flagCutTriangles(theCircle);
			std::vector<CXXCoord<CXXCoord_ftype> > &ballsEdges(raggedEdges[ballJ]);
			std::vector<CXXCoord<CXXCoord_ftype> >::iterator ballsEdgesEnd(ballsEdges.end());
			for (std::vector<CXXCoord<CXXCoord_ftype> >::iterator ballsEdge = ballsEdges.begin();
				 ballsEdge!= ballsEdgesEnd;
				 ballsEdge++){
				ballSphere.addVertex(CXXCircleNode(&theCircle, 0, *ballsEdge, 0));
			}
		}
#pragma omp critical (mainTriangles)
		aSurface.upLoadSphere(ballSphere, theBall->getRadius(), insideOrOutside);
	}
}

