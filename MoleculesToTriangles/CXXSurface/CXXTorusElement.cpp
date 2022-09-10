/*
 *  CXXTorusElement.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
#include "CXXSurfaceVertex.h"
#include "CXXTorusElement.h"
#include "CXXTorusNode.h"
#include "CXXCoord.h"
#include "CXXSurface.h"
#include "CXXTriangle.h"
#include "CXXCircle.h"
#include "CXXCircleNode.h"
#include "CXXNewHood.h"

CXXCircle CXXTorusElement::nullCircle = CXXCircle();

CXXTorusElement::CXXTorusElement() : 
theCircle(nullCircle),
debug(0)
{
	init();
}

CXXTorusElement::~CXXTorusElement()
{
}

void CXXTorusElement::init()
{
	flatTriangles.resize(0);
}

CXXTorusElement::CXXTorusElement(const CXXCircle &aCircle, int iEdge, double delta, double probeRadius) :
theCircle(aCircle),
rProbe(probeRadius),
debug(0)
{
	init();
	//	cout << "Hood Atom " << aCircle.getParent()->getAtomI()->element << " Circle Atom " << aCircle.getAtomJ()->element << endl;
	const CXXCoord<CXXCoord_ftype>xAxis(1.,0.,0.,0.);
	const CXXCoord<CXXCoord_ftype>yAxis(0.,1.,0.,0.);
	const CXXCoord<CXXCoord_ftype>zAxis(0.,0.,1.,0.);
	
	//Get centre, axis, and relevant radii of the torus; 
	torusCentre = theCircle.getCentreOfCircle(); 
	torusAxis = theCircle.getNormal();
	rTraj = theCircle.getRadiusOfCircle();
	
	const CXXCircleNode &startNode(*theCircle.start(iEdge));
	const CXXCircleNode &endNode(*theCircle.stop(iEdge));
	v1unit = startNode.getUnitRadius();
	v2unit = endNode.getUnitRadius();
	CXXCoord<CXXCoord_ftype>temp1 = startNode.getCoord();
	
	n1unit = v1unit ^ torusAxis;
	
	// Convert from position vectors to start and end range in theta
	CXXCoord<CXXCoord_ftype>diff = theCircle.getCentreOfSecondSphere()-temp1;
	diff.normalise();
	theta1 = acos(torusAxis*diff);
	diff = theCircle.getCentreOfSphere()-temp1;
	diff.normalise();
	theta2 = acos(torusAxis*diff);
	double halfWay1= (theta2+theta1)/2.;
	double halfWay2 = halfWay1;
	
	// Now in omega:  by definition, the omega value of vertex1 is zero, 
	absoluteStartOmega = startNode.getAngle();
	
	omega1 = 0;
	//startNode having an assigned otherCircle is the coded flag to indicate that this 
	//node is a real one coming from circle intersection, rather than one assigned
	//as an arbitrary referenceUnitRadius
	if (startNode.getOtherCircle()){
		omega2 = endNode.getAngle() - startNode.getAngle();
		while (omega2 < 0.) omega2 += 2.*M_PI;
	} 
	else {
		omega2 = 2.*M_PI;
	}
	
	deltaOmega=omega2-omega1;
	deltaTheta=theta2-halfWay2;
	
	nOmegaSteps = 1;
	while (fabs(deltaOmega)>delta){
		nOmegaSteps*=2;
		deltaOmega /= 2.;
	}
	
	nThetaSteps = 1;
	while (fabs(deltaTheta)>delta){
		nThetaSteps*=2;
		deltaTheta/=2.;
	}
	
	// Deal with case that the radius of the probe trajectory is less than probe Radius
	//Note I keep the same deltaTheta step so that it stitches more naturally with 
	//the corresponding nodes that flank it.  They will need to include the point at
	//the tip of the pointy bit into their edges, so I make space for and store this point
	if (rTraj <= probeRadius){
		halfWay1 = asin(rTraj/probeRadius);
		halfWay2 = M_PI - halfWay1;
		
		nThetaSteps = int((fabs(theta2-halfWay2)) / deltaTheta) + 1;
	}
	
	nodes.resize((nThetaSteps+1)*(nOmegaSteps+1));
	
	double omega = omega1;
	for (int iOmega=0; iOmega<=nOmegaSteps; iOmega++){
		double theta = theta2;
		for (int iTheta=0; iTheta<=nThetaSteps; iTheta++){
			CXXTorusNode aNode(theta,omega);
			CXXCoord<CXXCoord_ftype>xyz(coordFromThetaOmega(theta, omega));
			aNode.setCoord(xyz);
			nodes[iTheta+(iOmega*(nThetaSteps+1))] = aNode;
			theta-= deltaTheta;
			if (deltaTheta > 0.) theta = (theta > halfWay2 ? theta :halfWay2);
			else theta = (theta < halfWay2 ? theta :halfWay2);
		}
		omega+=deltaOmega;
	}
	
	for (int iOmega=0; iOmega<nOmegaSteps; iOmega++){
		for (int iTheta=0; iTheta<nThetaSteps; iTheta++){
			size_t iTriangle = flatTriangles.size();
			flatTriangles.push_back(CXXTriangle(((iOmega+1)*(nThetaSteps+1))+(iTheta+0), 
												((iOmega+0)*(nThetaSteps+1))+(iTheta+0),
												((iOmega+0)*(nThetaSteps+1))+(iTheta+1),
												iTriangle));
			if (iTheta == 0) edgeTriangles.push_back(&flatTriangles.back());
			iTriangle = flatTriangles.size();
			flatTriangles.push_back(CXXTriangle(((iOmega+0)*(nThetaSteps+1))+(iTheta+1),
												((iOmega+1)*(nThetaSteps+1))+(iTheta+1), 
												((iOmega+1)*(nThetaSteps+1))+(iTheta+0),
												iTriangle));
		}
	}
	//Add pointer to appropriate atom
	for (unsigned int iNode = 0; iNode<nodes.size(); iNode++){
		nodes[iNode].setAtom(theCircle.getParent()->getAtomI());
	}
}

size_t CXXTorusElement::addNode(CXXTorusNode &aNode){
	CXXTorusNode newNode(aNode);
	int debug = 0;
	
	CXXCoord<CXXCoord_ftype>xyz = coordFromThetaOmega(newNode.getTheta(), newNode.getOmega());
	newNode.setCoord(xyz);
	nodes.push_back(newNode);
	if (debug) cout << "Node number " << nodes.size() << endl;
	return (nodes.size() - 1);
}

const size_t CXXTorusElement::numberOfTorusNodes(void){
	return nodes.size();
}

const CXXCoord<CXXCoord_ftype>CXXTorusElement::probeAtOmega(double omega) const{
	CXXCoord<CXXCoord_ftype>temp1 = v1unit;
	temp1.scale(cos(omega));
	CXXCoord<CXXCoord_ftype>temp2 = n1unit;
	temp2.scale(sin(omega));
	CXXCoord<CXXCoord_ftype>temp3 = temp1 + temp2;
	temp3.scale(rTraj);
	return torusCentre + temp3;
}

const CXXCoord<CXXCoord_ftype>CXXTorusElement::normalToProbeAtTheta(CXXCoord<CXXCoord_ftype>&p, double theta) const {
	CXXCoord<CXXCoord_ftype>temp1 = torusAxis;
	temp1.scale(cos(theta));
	CXXCoord<CXXCoord_ftype>temp2 = torusCentre - p;
	temp2.normalise();
	temp2.scale(sin(theta));
	return temp1 + temp2;
}

const CXXCoord<CXXCoord_ftype>CXXTorusElement::coordFromThetaOmega(double theta, double omega) const{
	CXXCoord<CXXCoord_ftype>p(probeAtOmega(omega));
	CXXCoord<CXXCoord_ftype>normal(normalToProbeAtTheta(p, theta));
	normal.scale(rProbe);
	return  p + normal;
}

int CXXTorusElement::upload(CXXSurface *aSurface){
	//	Add vertices to surface
	size_t oldVertexCount;
	{
		double verticesBuffer[nodes.size()*3];// = new double[nodes.size()*3];
		double accessiblesBuffer[nodes.size()*3];// = new double[nodes.size()*3];
		double normalsBuffer[nodes.size()*3];// = new double[nodes.size()*3];
		for (unsigned int i=0; i< nodes.size(); i++){
			for (int j=0; j<3; j++) verticesBuffer[3*i+j] = nodes[i].coord().element(j);
			CXXCoord<CXXCoord_ftype>accessible = probeAtOmega(nodes[i].getOmega());
			for (int j=0; j<3; j++) accessiblesBuffer[3*i+j] = accessible.element(j);
			CXXCoord<CXXCoord_ftype>normal = normalToProbeAtTheta(accessible, nodes[i].getTheta());
			normal.scale(-1.);
			for (int j=0; j<3; j++) normalsBuffer[3*i+j] = normal.element(j);
		}
		oldVertexCount = aSurface->numberOfVertices();
		aSurface->updateWithVectorData(nodes.size(), "vertices", oldVertexCount, verticesBuffer);
		aSurface->updateWithVectorData(nodes.size(), "accessibles", oldVertexCount, accessiblesBuffer);
		aSurface->updateWithVectorData(nodes.size(), "normals",  oldVertexCount, normalsBuffer);
		//delete [] verticesBuffer;
		//delete [] accessiblesBuffer;
		//delete [] normalsBuffer;
	}	
	//Add atom pointers to the surface
	{
		void *atomBuffer[nodes.size()];// = new void*[nodes.size()];
		for (unsigned int i=0; i< nodes.size(); i++){
			atomBuffer[i] = (void *)nodes[i].getAtom();
		}
		aSurface->updateWithPointerData(nodes.size(), "atom", oldVertexCount, atomBuffer);
		//delete [] atomBuffer;
	}
	// Add triangles to surface
	{
		int triangleBuffer[flatTriangles.size()*3];// = new int[flatTriangles.size()*3];
		int nToDraw = 0;
        list<CXXTriangle  >::iterator trianglesEnd(flatTriangles.end());
        for (list<CXXTriangle  >::iterator triangle=flatTriangles.begin();
             triangle != trianglesEnd;
             ++triangle){
            CXXTriangle &flatTriangle(*triangle);
			if (flatTriangle.doDraw()){
				for (unsigned int j=0; j<3; j++){
					//Note the 2-j, this changes the sense of the triangle to reflect the fact that we
					//actually visualise the inside of the torus
					triangleBuffer[(3*nToDraw)+j] = int(flatTriangle[j] + oldVertexCount);
				}
				nToDraw++;
			}
		}
		aSurface->extendTriangles(triangleBuffer, nToDraw);
		//delete [] triangleBuffer;
	}
	return 0;
}	

const CXXTorusNode &CXXTorusElement::getNode(const int i) const{
	return nodes[i];
} 

void CXXTorusElement::addEdgeVertex(CXXCircleNode &aNode){
	double omega = aNode.getAngle() - absoluteStartOmega;
	while (omega < 0.) omega += 2.*M_PI;
	if (omega < omega2){
		//This vertex falls somewhere in the range of the segment:  find which of the edgeTriagles it falls within
		int triangleFound = 0;
		list<CXXTriangle *  >::iterator matchingTriangle;
		list<CXXTriangle *  >::iterator edgeTrianglesEnd = edgeTriangles.end();
		for (list<CXXTriangle *  >::iterator triangle = edgeTriangles.begin();
             triangle != edgeTrianglesEnd && !triangleFound;
             ++triangle){

			CXXTriangle &theFlatTriangle(**triangle);
			double omegaStart = nodes[theFlatTriangle[1]].getOmega();
			double omegaEnd   = nodes[theFlatTriangle[0]].getOmega();
			if (omega >= omegaStart && omega <= omegaEnd){
				triangleFound = 1;
				matchingTriangle = triangle;
			}
		}
		
		if (triangleFound){
			CXXTorusNode aNode(theta2,omega);
			CXXCoord<CXXCoord_ftype>xyz(coordFromThetaOmega(theta2, omega));
			aNode.setCoord(xyz);
			aNode.setAtom(theCircle.getParent()->getAtomI());
			nodes.push_back(aNode);
			
			CXXTriangle &theFT(**matchingTriangle);
			theFT.setDoDraw(0);
			edgeTriangles.erase(matchingTriangle);

            flatTriangles.push_back(CXXTriangle(theFT[0], nodes.size()-1, theFT[2]));
            edgeTriangles.push_back(&flatTriangles.back());
            flatTriangles.push_back(CXXTriangle(nodes.size()-1, theFT[1], theFT[2]));
            edgeTriangles.push_back(&flatTriangles.back());
		}
	}
}


