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
 *  CXXSphereElement.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  
 *
 */

#include "CXXSphereElement.h"
#include <CXXCoord.h>
#include <CXXSphereTriangle.h>
#include <CXXSphereTriangleEdge.h>
#include <CXXSphereNode.h>
#include <CXXSphereFlatTriangle.h>
#include <CXXSurface.h>
#include <CXXCircleNode.h>
#include <CXXCircle.h>
#include <CXXNewHood.h>
#include <CXXTorusElement.h>
#include <math.h>
#include "mmdb_manager.h"
#include "mmdb_tables.h"


void CXXSphereElement::init(){
	theVertices.reserve(720);
	theTriangles.reserve(720);
	flatTriangles.reserve(2000);
	theEdges.reserve(720);
	nDrawnTriangles=0;
}
CXXSphereElement::CXXSphereElement() : 
theAtom(0)
{
	init();
}

CXXSphereElement::~CXXSphereElement(){
}

CXXSphereElement::CXXSphereElement (PCAtom anAtom, double del) : 
theAtom(anAtom), 
deltaRadians(del)
{
	init();
	theCentre = CXXCoord(anAtom->x, anAtom->y, anAtom->z);
	calculate();
}

CXXSphereElement::CXXSphereElement(const CXXCoord &position, double radius, double del) :
theCentre(position), 
theRadius(radius),
deltaRadians(del) {
	init ();	
	calculate();
	for (unsigned i=0; i<flatTriangles.size(); i++){
		for (int j=0; j<3; j++){
			flatTriangles[i].setEdgeCircle(j,0);
		}
	}
}

int CXXSphereElement::calculate(){
	theVertices.resize(0);
	theTriangles.resize(0);
	flatTriangles.resize(0);
	
	CXXCoord xAxis(1.,0.,0.,0.);
	CXXCoord yAxis(0.,1.,0.,0.);
	CXXCoord zAxis(0.,0.,1.,0.);
	CXXCoord normal1(zAxis);
	CXXCoord normal2(xAxis);
	CXXCoord normal3(yAxis);
	xAxis.scale(theRadius);
	yAxis.scale(theRadius);
	zAxis.scale(theRadius);
	CXXSphereNode v1(theCentre + xAxis);
	int iV1 = addVertex(v1);
	CXXSphereNode v2(theCentre + yAxis);
	int iV2 = addVertex(v2);
	CXXSphereNode v3(theCentre + zAxis);
	int iV3 = addVertex(v3);
	
	//Make two surface patches that sum to one :  First inside of an octant		
	CXXSphereTriangleEdge edge1(normal1, iV1, iV2, 
								theCentre, theCentre, 
								theRadius, this);
	int iEdge1 = addEdge(edge1);
	CXXSphereTriangleEdge edge2(normal2, iV2, iV3, 
								theCentre, theCentre, 
								theRadius, this);
	int iEdge2 = addEdge(edge2);
	CXXSphereTriangleEdge edge3(normal3, iV3, iV1, 
								theCentre, theCentre,
								theRadius, this);
	int iEdge3 = addEdge(edge3);
	int vertices[] = {iV1, iV2, iV3};
	int edges[] = {iEdge1, iEdge2, iEdge3};
	CXXSphereTriangle aTriangle(this, vertices, edges, theRadius, theCentre);
	addTriangle(aTriangle);
	
	while (theTriangles.size()){
		if (!theTriangles.back().bisect(delta())){
			flattenLastTriangle();
		}
	}
	//Then outside of an octant
	normal1.scale(-1.);
	normal2.scale(-1.);
	normal3.scale(-1.);
	
	edge1 = CXXSphereTriangleEdge(normal1, iV2, iV1, 
								  theCentre, theCentre, 
								  theRadius, this);
	iEdge1 = addEdge(edge1);
	edge2 = CXXSphereTriangleEdge(normal2, iV3, iV2, 
								  theCentre, theCentre, 
								  theRadius, this);
	iEdge2 = addEdge(edge2);
	edge3 = CXXSphereTriangleEdge(normal3, iV1, iV3, 
								  theCentre, theCentre, 
								  theRadius, this);
	iEdge3 = addEdge(edge3);
	int vertices1[] = {iV1, iV3, iV2};
	int edges1[] = {iEdge3, iEdge2, iEdge1};
	aTriangle = CXXSphereTriangle(this, vertices1, edges1, theRadius, theCentre);
	addTriangle(aTriangle);
	
	while (theTriangles.size()){
		if (!theTriangles.back().bisect(delta())){
			flattenLastTriangle();
		}
	}
	nDrawnTriangles = countDrawnTriangles();
	return 0;
}

int CXXSphereElement::addTriangularPatch(const CXXCoord &u1_in, const CXXCoord &u2_in, 
										 const CXXCoord &u3_in, CAtom *anAtom,
										 vector <PointyBit> &pointyBits)
{	
	
	int oldNFlat = nFlatTriangles();
	for (int i=0; i<oldNFlat; i++){
		CXXSphereFlatTriangle &flat(flatTriangles[i]);
		flat.setEdgeCircle(0,0);
		flat.setEdgeCircle(1,0);
		flat.setEdgeCircle(2,0);
	}
	
	theAtom = anAtom;
	CXXCoord u1(u1_in);
	CXXCoord u2(u2_in);
	CXXCoord u3(u3_in);
	
	CXXCoord normal1(u1^u2);
	normal1.normalise();
	
	//This method works only for triangle for which all angles are < 180.  Here, we check this
	//by seeing whether the normal to u1->u2 points towards or away from u3, if not, then we reverse
	//the order of the vertices
	if (normal1.dot(u3)<0.) {
		CXXCoord temp(u1);
		u1 = u2;
		u2 = temp;
		normal1 *= -1.;
	}
	CXXCoord normal2(u2^u3);
	normal2.normalise();
	CXXCoord normal3(u3^u1);
	normal3.normalise();
	
	CXXCoord v1Coord(u1);
	CXXCoord v2Coord(u2);
	CXXCoord v3Coord(u3);
	v1Coord *=theRadius;
	v2Coord *= theRadius;
	v3Coord *= theRadius;
	CXXSphereNode v1(theCentre + v1Coord);
	int iV1 = addVertex(v1);
	CXXSphereNode v2(theCentre + v2Coord);
	int iV2 = addVertex(v2);
	CXXSphereNode v3(theCentre + v3Coord);
	int iV3 = addVertex(v3);
	
	CXXCircle edgeCircle1, edgeCircle2, edgeCircle3;
	edgeCircle1.setNormal(normal1);
	edgeCircle1.setCentreOfCircle(theCentre);
	edgeCircle1.setRadiusOfCircle(theRadius);
	
	edgeCircle2.setNormal(normal2);
	edgeCircle2.setCentreOfCircle(theCentre);
	edgeCircle2.setRadiusOfCircle(theRadius);
	
	edgeCircle3.setNormal(normal3);
	edgeCircle3.setCentreOfCircle(theCentre);
	edgeCircle3.setRadiusOfCircle(theRadius);
	
	CXXSphereTriangleEdge edge1(normal1, iV1, iV2, 
								theCentre, theCentre, 
								theRadius, this);
	if(!pointyBits[0].isNull) edge1.setCircle(&edgeCircle1);
	else edge1.setCircle(0);
	int iEdge1 = addEdge(edge1);
	
	
	CXXSphereTriangleEdge edge2(normal2, iV2, iV3, 
								theCentre, theCentre, 
								theRadius, this);
	if(!pointyBits[1].isNull) edge2.setCircle(&edgeCircle2);
	else edge2.setCircle(0);
	int iEdge2 = addEdge(edge2);
	CXXSphereTriangleEdge edge3(normal3, iV3, iV1, 
								theCentre, theCentre, 
								theRadius, this);
	if(!pointyBits[2].isNull) edge3.setCircle(&edgeCircle3);
	else edge3.setCircle(0);
	int iEdge3 = addEdge(edge3);
	
	int vertices[] = {iV1, iV2, iV3};
	int edges[] = {iEdge1, iEdge2, iEdge3};
	
	CXXSphereTriangle aTriangle(this, vertices, edges, theRadius, theCentre);
	aTriangle.setAtom(anAtom);
	addTriangle(aTriangle);
	
	while (theTriangles.size()){
		if (!theTriangles.back().bisect(delta())){
			flattenLastTriangle();
		}
	}
	int newNFlat = nFlatTriangles();
	for (int i=oldNFlat; i<newNFlat; i++){
		CXXSphereFlatTriangle &flat(flatTriangles[i]);
		theVertices[flat[0]].setAtom(anAtom);
		theVertices[flat[1]].setAtom(anAtom);
		theVertices[flat[2]].setAtom(anAtom);
	}		
	
	if (!pointyBits[0].isNull){
		edgeCircle1.setArbitraryReference();
		flagCutTriangles(edgeCircle1);
		CXXCircleNode cn1(&edgeCircle1, 0, pointyBits[0].coord, 0);
		cn1.setReference(edgeCircle1.getReferenceUnitRadius());
		addVertex(cn1);
	}
	
	if (!pointyBits[1].isNull){
		edgeCircle2.setArbitraryReference();
		flagCutTriangles(edgeCircle2);
		CXXCircleNode cn1(&edgeCircle2, 0, pointyBits[1].coord, 0);
		cn1.setReference(edgeCircle2.getReferenceUnitRadius());
//		addVertex(cn1);
	}
	
	if (!pointyBits[2].isNull){
		edgeCircle3.setArbitraryReference();
		flagCutTriangles(edgeCircle3);
		CXXCircleNode cn1(&edgeCircle3, 0, pointyBits[2].coord, 0);
		cn1.setReference(edgeCircle3.getReferenceUnitRadius());
		addVertex(cn1);
	}
	
	
	for (unsigned i=0; i<flatTriangles.size(); i++){
		CXXSphereFlatTriangle &flat(flatTriangles[i]);
		flat.setEdgeCircle(0,0);
		flat.setEdgeCircle(1,0);
		flat.setEdgeCircle(2,0);
	}
	clearCutFlags();
	nDrawnTriangles = countDrawnTriangles();
	
	
	return 0;
}

CXXSphereElement::CXXSphereElement (const CXXSphereElement &oldOne) :
theAtom ( oldOne.getAtom()),
theCentre ( oldOne.centre()),
theVertices(oldOne.getVertices()),
theTriangles(oldOne.getTriangles()),
flatTriangles(oldOne.getFlatTriangles()),
theRadius ( oldOne.radius()),
deltaRadians ( oldOne.delta()),
nDrawnTriangles (oldOne.getNDrawnTriangles())
{
	//Makes a complete independent copy of the oldOne;
	
	for (int i=0; i<oldOne.nTriangles(); i++){
		theTriangles[i].setSphereElement(this);
	}
	theEdges.resize(oldOne.nEdges());
	for (int i=0; i<oldOne.nEdges(); i++){
		theEdges[i]=CXXSphereTriangleEdge(oldOne.edge(i));
		theEdges[i].setSphereElement(this);
	}
}

CXXSphereElement::CXXSphereElement(const CXXCircleNode &aNode, double delta, double radius_in, int selHnd, vector<PointyBit> &pointyBits) :
theCentre(aNode.getCoord()), theRadius(radius_in), deltaRadians(delta){
	init();
	
	CAtom *atomK=aNode.getAtomK();
	CAtom *atomJ=aNode.getAtomJ();
	CAtom *atomI=aNode.getAtomI();
	
	CXXCoord u1(atomK->x, atomK->y, atomK->z);
	CXXCoord u2(atomJ->x, atomJ->y, atomJ->z);
	CXXCoord u3(atomI->x, atomI->y, atomI->z);
	
	u1 = u1 - theCentre;
	u1.normalise();
	u2 = u2-theCentre;
	u2.normalise();
	u3 = u3-theCentre;
	u3.normalise();
	
	CXXCoord u12;
	u12 = u1 + u2;
	u12.normalise();
	CXXCoord u23;
	u23 = u2 + u3;
	u23.normalise();
	CXXCoord u31;
	u31 = u3 + u1;
	u31.normalise();
	CXXCoord u123(u1 + u2 + u3);
	u123.normalise();
	
	
	
	vector<PointyBit> patchPointyBits(3);

	//Corner of triangle closest to atom I
	if (atomI->isInSelection(selHnd)) {
		for (unsigned i=0; i<3; i++) patchPointyBits[i].isNull = 1;
		for (unsigned i=0; i<pointyBits.size(); i++){
			if (pointyBits[i].atomI == atomI){
				if (pointyBits[i].atomJ == atomJ){
					patchPointyBits[1] = pointyBits[i];
				}
			}
		}
		addTriangularPatch(u123, u3,  u23, atomI, patchPointyBits);
		
		for (unsigned i=0; i<3; i++) patchPointyBits[i].isNull = 1;
		for (unsigned i=0; i<pointyBits.size(); i++){
			if (pointyBits[i].atomI == atomI){
				if (pointyBits[i].atomJ == atomK){
					patchPointyBits[1] = pointyBits[i];
				}
			}
		}
		addTriangularPatch(u123, u31, u3,  atomI, patchPointyBits);
	}
	
	//Corner of triangle closest to atom J
	if (atomJ->isInSelection(selHnd)){
		for (unsigned i=0; i<3; i++) patchPointyBits[i].isNull = 1;
		for (unsigned i=0; i<pointyBits.size(); i++){
			if (pointyBits[i].atomI == atomJ){
				if (pointyBits[i].atomJ == atomK){
					patchPointyBits[1] = pointyBits[i];
				}
			}
		}
		addTriangularPatch(u123, u2,  u12, atomJ, patchPointyBits);
		for (unsigned i=0; i<3; i++) patchPointyBits[i].isNull = 1;
		for (unsigned i=0; i<pointyBits.size(); i++){
			if (pointyBits[i].atomI == atomJ){
				if (pointyBits[i].atomJ == atomI){
					patchPointyBits[1] = pointyBits[i];
				}
			}
		}
		addTriangularPatch(u123, u23, u2,  atomJ, patchPointyBits);
	}
	
	//Corner of triangle closest to atom K
	if (atomK->isInSelection(selHnd)) {
		for (unsigned i=0; i<3; i++) patchPointyBits[i].isNull = 1;
		for (unsigned i=0; i<pointyBits.size(); i++){
			if (pointyBits[i].atomI == atomK){
				if (pointyBits[i].atomJ == atomI){
					patchPointyBits[1] = pointyBits[i];
				}
			}
		}
		addTriangularPatch(u123, u1,  u31, atomK, patchPointyBits);
		for (unsigned i=0; i<3; i++) patchPointyBits[i].isNull = 1;
		for (unsigned i=0; i<pointyBits.size(); i++){
			if (pointyBits[i].atomI == atomK){
				if (pointyBits[i].atomJ == atomJ){
					patchPointyBits[1] = pointyBits[i];
				}
			}
		}
		addTriangularPatch(u123, u12, u1,  atomK, patchPointyBits);
	}
	
}

int CXXSphereElement::addVertex(const CXXSphereNode &vert){
	int debug = 0;
	int iVertex = theVertices.size();
	theVertices.push_back(vert);
	if (debug) cout << "Added vertex number " << theVertices.size() << endl;
	return iVertex;
}

int CXXSphereElement::addTriangle(const CXXSphereTriangle &newTri){
	int debug = 0;
	int iTriangle = theTriangles.size();
	theTriangles.push_back(newTri);
	if (debug) cout << "Added triangle number " << theTriangles.size() << endl;
	return iTriangle;
}

int CXXSphereElement::addEdge(const CXXSphereTriangleEdge &newEdge){
	int iEdge = theEdges.size();
	theEdges.push_back(newEdge);
	return iEdge;
}

void CXXSphereElement::flattenLastTriangle(void){
	int debug = 0;
	if (debug)  cout << "Finished with Triangle " << flatTriangles.size() <<endl; 
	int i = theTriangles.back().vertex(0);
	int j = theTriangles.back().vertex(1);
	int k = theTriangles.back().vertex(2);
	CXXSphereFlatTriangle newTriangle(i, k, j, flatTriangles.size());
	newTriangle.setAtom (theTriangles.back().getAtom());
	newTriangle.setEdgeCircle(0, theEdges[theTriangles.back().edge(2)].getCircle());
	newTriangle.setEdgeCircle(1, theEdges[theTriangles.back().edge(1)].getCircle());
	newTriangle.setEdgeCircle(2, theEdges[theTriangles.back().edge(0)].getCircle());
	
	flatTriangles.push_back(newTriangle);
	theTriangles.pop_back();
	return;
} 

const CXXCoord &CXXSphereElement::centre() const{
	return theCentre;
}
const double CXXSphereElement::radius() const{
	return theRadius;
}
const unsigned CXXSphereElement::nVertices() const{
	return theVertices.size();
}
const int CXXSphereElement::nEdges() const {
	return theEdges.size();
}
const int CXXSphereElement::nTriangles() const{
	return theTriangles.size();
}
const int CXXSphereElement::nFlatTriangles() const{
	return flatTriangles.size();
}
const CXXSphereNode &CXXSphereElement::vertex(const int iVertex) const{
	return theVertices[iVertex];
}
void CXXSphereElement::moveVertex( const int iVertex, const CXXCoord &position){
	theVertices[iVertex].setVertex(position);
	return ;
}
const CXXSphereTriangle &CXXSphereElement::triangle(const int iTriangle) const{
	return theTriangles[iTriangle];
}
const CXXSphereTriangleEdge &CXXSphereElement::edge(const int iEdge) const{
	return theEdges[iEdge];
}
const CXXSphereFlatTriangle &CXXSphereElement::flatTriangle(const int iFlatTriangle) const{
	return flatTriangles[iFlatTriangle];
}
void CXXSphereElement::hideFlatTriangle(const int iFlatTriangle){
	flatTriangles[iFlatTriangle].setDoDraw(0);
	nDrawnTriangles--;
}

int CXXSphereElement::trimBy(const CXXCircle &aCircle){
	unsigned oldNTriangles = flatTriangles.size();
	for (unsigned int iTriangle = 0; iTriangle<oldNTriangles; iTriangle++){
		CXXSphereFlatTriangle &flatTriangle(flatTriangles[iTriangle]);
		if (flatTriangle.doDraw()){
			int iNodes[] = {flatTriangle[0], flatTriangle[1], flatTriangle[2]};
			
			int isInFront[] = {
				aCircle.accIsBehind(theVertices[iNodes[0]].vertex()),
				aCircle.accIsBehind(theVertices[iNodes[1]].vertex()),
				aCircle.accIsBehind(theVertices[iNodes[2]].vertex())};
			int nInFront = isInFront[0] + isInFront[1] + isInFront[2];
			
			//First deal with all vertices in front
			if (nInFront == 3){
			}
			//Next case is all vertices behind
			else if (nInFront == 0 ){
				hideFlatTriangle(iTriangle);
				for (int j=0; j<3; j++) theVertices[iNodes[j]].setDoDraw(0);
			}
			//Next case is one vertex in front, other two vertices behind
			else if (nInFront == 1) {
				int spun=0;
				for (int i=0; i<3 && !spun; i++){
					if (isInFront[i] && (!isInFront[(i+1)%3]) && (!isInFront[(i+2)%3])){
						spun = 1;
						hideFlatTriangle(iTriangle);
						
						theVertices[iNodes[(i+1)%3]].setDoDraw(0);
						theVertices[iNodes[(i+2)%3]].setDoDraw(0);
						
						CXXCoord newVertex1 = aCircle.accPlaneIntersect
							(theVertices[iNodes[i]].vertex(), 
							 theVertices[iNodes[(i+1)%3]].vertex());
						
						//If the triangle edge we are examining was already generated by trimming by another
						//circle, then we have to calculate the intersection vertex more carefully, bearing in mind
						//that these are circle intersections
						
						const CXXCircle *oldCircle = flatTriangle.getEdgeCircle(i);
						if (oldCircle){
							vector<CXXCoord> theIntersections(2);							 
							oldCircle->meetsCircle(aCircle, theIntersections);
							CXXCoord diff1(theIntersections[0]- newVertex1);
							CXXCoord diff2(theIntersections[1]- newVertex1);
							newVertex1 = (diff1.get3DLengthSq() < diff2.get3DLengthSq()?theIntersections[0]:theIntersections[1]);
						}
						
						CXXSphereNode newNode1(newVertex1);
						iNodes[(i+1)%3] = addVertex(newNode1);
						theVertices[iNodes[(i+1)%3]].setIntersector(&aCircle);
						
						CXXCoord newVertex2 = aCircle.accPlaneIntersect
							(theVertices[iNodes[i]].vertex(), 
							 theVertices[iNodes[(i+2)%3]].vertex());
						
						const CXXCircle *oldCircle1 = flatTriangle.getEdgeCircle((i+2)%3);
						if (oldCircle1){
							vector<CXXCoord> theIntersections(2);							 
							oldCircle1->meetsCircle(aCircle, theIntersections);
							CXXCoord diff1(theIntersections[0]- newVertex2);
							CXXCoord diff2(theIntersections[1]- newVertex2);
							newVertex2 = (diff1.get3DLengthSq() < diff2.get3DLengthSq()?theIntersections[0]:theIntersections[1]);
						}
						
						CXXSphereNode newNode2(newVertex2);
						iNodes[(i+2)%3] = addVertex(newNode2);
						theVertices[iNodes[(i+2)%3]].setIntersector(&aCircle);
						
						flatTriangles.push_back(CXXSphereFlatTriangle(iNodes[i], iNodes[(i+1)%3], iNodes[(i+2)%3]));
						flatTriangles.back().setAtom(flatTriangles[iTriangle].getAtom());
						flatTriangles.back().setEdgeCircle(0, oldCircle);
						flatTriangles.back().setEdgeCircle(1, &aCircle);
						flatTriangles.back().setEdgeCircle(2, oldCircle1);
						nDrawnTriangles++;
						
					}
				}
			}
			//Now deal with case that one coordinate is chopped off, and two are infront of plane
			else {
				int spun=0;
				for (int i=0; i<3 && !spun; i++){
					if (!isInFront[i] && isInFront[(i+1)%3] && isInFront[(i+2)%3]){
						
						spun = 1;
						hideFlatTriangle(iTriangle);
						theVertices[iNodes[i]].setDoDraw(0);
						
						CXXCoord newVertex1 = aCircle.accPlaneIntersect
							(theVertices[iNodes[(i+1)%3]].vertex(), 
							 theVertices[iNodes[i]].vertex());
						
						const CXXCircle *oldCircle = flatTriangles[iTriangle].getEdgeCircle(i);
						if (oldCircle){
							vector<CXXCoord> theIntersections(2);							 
							oldCircle->meetsCircle(aCircle, theIntersections);
							CXXCoord diff1(theIntersections[0]- newVertex1);
							CXXCoord diff2(theIntersections[1]- newVertex1);
							newVertex1 = (diff1.get3DLengthSq() < diff2.get3DLengthSq()?theIntersections[0]:theIntersections[1]);
						}
						
						CXXSphereNode newNode1(newVertex1);
						int iNewNode1 = addVertex(newNode1);
						theVertices[iNewNode1].setIntersector(&aCircle);
						
						flatTriangles.push_back(CXXSphereFlatTriangle(iNewNode1, iNodes[(i+1)%3], iNodes[(i+2)%3]));
						flatTriangles.back().setAtom(flatTriangles[iTriangle].getAtom());
						flatTriangles.back().setEdgeCircle(0, flatTriangles[iTriangle].getEdgeCircle(i));
						flatTriangles.back().setEdgeCircle(1, flatTriangles[iTriangle].getEdgeCircle((i+1)%3));
						flatTriangles.back().setEdgeCircle(2, 0);
						nDrawnTriangles++;
						
						CXXCoord newVertex2 = aCircle.accPlaneIntersect
							(theVertices[iNodes[(i+2)%3]].vertex(), 
							 theVertices[iNodes[i]].vertex());
						
						const CXXCircle *oldCircle1 = flatTriangles[iTriangle].getEdgeCircle((i+2)%3);
						if (oldCircle1){
							vector<CXXCoord> theIntersections(2);							 
							oldCircle1->meetsCircle(aCircle, theIntersections);
							CXXCoord diff1(theIntersections[0]- newVertex2);
							CXXCoord diff2(theIntersections[1]- newVertex2);
							newVertex2 = (diff1.get3DLengthSq() < diff2.get3DLengthSq()?theIntersections[0]:theIntersections[1]);
						} 
						
						CXXSphereNode newNode2(newVertex2);
						int iNewNode2 = addVertex(newNode2);
						theVertices[iNewNode2].setIntersector(&aCircle);
						
						flatTriangles.push_back(CXXSphereFlatTriangle(iNewNode2, iNewNode1, iNodes[(i+2)%3]));
						flatTriangles.back().setAtom(flatTriangles[iTriangle].getAtom());
						flatTriangles.back().setEdgeCircle(0, &aCircle);
						flatTriangles.back().setEdgeCircle(1, 0);
						flatTriangles.back().setEdgeCircle(2, flatTriangles[iTriangle].getEdgeCircle((i+2)%3));
						nDrawnTriangles++;
						
						
					}
				}
			}
			//No more cases I think
		}
	}
	
	return 0;
}

const double CXXSphereElement::delta() const{
	return deltaRadians;
}

int CXXSphereElement::translateBy (const CXXCoord &crd){
	theCentre = theCentre + crd;
	for (unsigned i=0; i<nVertices(); i++){
		const CXXCoord result = theVertices[i].vertex() + crd;
		theVertices[i].setVertex(result);
	}
	for (int i=0; i<nTriangles(); i++){
		CXXCoord result = theTriangles[i].centre() + crd;
		theTriangles[i].setCentre(result);		
	}
	for (int i=0; i<nEdges(); i++){
		CXXCoord result = theEdges[i].edgeCentre() + crd;
		theEdges[i].setEdgeCentre(result);
		result = theEdges[i].sphereCentre() + crd;
		theEdges[i].setSphereCentre(result);
	}
	return 0;
}

int CXXSphereElement::scaleBy (const double factor){
	double oldRadius = theRadius;
	theRadius = theRadius * factor;
	for (unsigned i=0; i<nVertices(); i++){
		CXXCoord result = theVertices[i].vertex() - theCentre;
		result.scale(theRadius / oldRadius);
		result = result + theCentre;
		theVertices[i].setVertex(result);
	}
	return 0;
}

int CXXSphereElement::setAtom(const CAtom *anAtom){
	theAtom = anAtom;
	return 0;
}

const CAtom *CXXSphereElement::getAtom() const{
	return theAtom;
}

const std::vector<CXXSphereNode> &CXXSphereElement::getVertices() const {
	return theVertices;
}
const std::vector<CXXSphereTriangle> &CXXSphereElement::getTriangles() const {
	return theTriangles;
}
const std::vector<CXXSphereTriangleEdge> &CXXSphereElement::getEdges() const {
	return theEdges;
}
const std::vector<CXXSphereFlatTriangle> &CXXSphereElement::getFlatTriangles() const {
	return flatTriangles;
}

int CXXSphereElement::addFlatTriangle(const CXXSphereFlatTriangle &theTriangle){
	flatTriangles.push_back(theTriangle);
	return 0;
}

CXXCoord CXXSphereElement::voronoiPoint(const CXXCoord &a, const CXXCoord &b, const CXXCoord &c) const{
	CXXCoord edge1(b-a);
	CXXCoord edge2(c-a);
	CXXCoord planeNormal(edge1^edge2);
	CXXCoord trivialResult(a+b+c);
	trivialResult /= 3.;
	if (planeNormal.get3DLengthSq()<1e-8){
		return trivialResult;
	}
	planeNormal.normalise();
	CXXCoord edgeNormal1(planeNormal^edge1);
	CXXCoord edgeNormal2(planeNormal^edge2);
	
	CXXCoord midPoint1(edge1);
	midPoint1 /= 2.;
	midPoint1 = midPoint1 + a;
	
	CXXCoord midPoint2(edge2);
	midPoint2 /= 2.;
	midPoint2 = midPoint2 + a;
	
	double k2 = (edgeNormal1[1]*(midPoint2[0]-midPoint1[0]) - edgeNormal1[0]*(midPoint2[1]-midPoint1[1])) /
		(edgeNormal1[0]*edgeNormal2[1]-edgeNormal2[0]*edgeNormal1[1]);
	
	CXXCoord result(edgeNormal2);
	result *= k2;
	result = result + midPoint2;
	//Confirm that the voronoi point lies on the triangle, otherwise, shift it to as close as
	//possible
	CXXCoord diff(result - a);
	double c1 = diff.dot(edge1);
	double c2 = diff.dot(edge2);
	if (c1> 0. && c1 < 1. && c2 > 0. && c2 < 1. && c1 + c2 < 1.) return result;
	else return trivialResult;
	
}

void CXXSphereElement::flagCutTriangles(const CXXCircle &theCircle){
	//First collect the list of triangles that have this circle as an edge, and evaluate their omega values
	//	edgeTriangles[&theCircle] = set <CXXSphereFlatTriangle *>;
	clearCutFlags();
	int nTriangles = flatTriangles.size();
	for (unsigned i=0; i<nTriangles; i++){
		CXXSphereFlatTriangle &flatTriangle(flatTriangles[i]);
		if (flatTriangle.doDraw()){
			int spun = 0;
			for (int j=0; j<3 && !spun; j++){
				if (flatTriangle.getEdgeCircle(j) == &theCircle){
					spun = 1;
					TriangleEdgePair edge(i,j);
					edgeTriangles.push_back(edge);
					
					const CXXCoord &v1(theVertices[flatTriangle[ j     ]].vertex());
					CXXCircleNode cn1(&theCircle, 0, v1, 0);
					cn1.setReference(theCircle.getReferenceUnitRadius());
					flatTriangle.setCircleNode( j     , cn1);
					
					const CXXCoord &v2(theVertices[flatTriangle[(j+1)%3]].vertex());
					CXXCircleNode cn2(&theCircle, 0, v2, 0);
					cn2.setReference(theCircle.getReferenceUnitRadius());
					flatTriangle.setCircleNode((j+1)%3, cn2);
				}
			}
		}
	}
}

void CXXSphereElement::addTorusVertices(const CXXTorusElement &aTorus){
	//Walk around the edge of this segment, evaluating location of extra vertices
	for (int iStep=0; iStep <= aTorus.getNOmegaSteps(); iStep++){
		double omega = aTorus.getAbsoluteStartOmega() + double(iStep) * aTorus.getDeltaOmega();
		while (omega < 0.) omega += 2.*M_PI;
		while (omega > 2.*M_PI) omega -= 2.*M_PI;
		const CXXCoord newNode(aTorus.probeAtOmega(omega-aTorus.getAbsoluteStartOmega()));
		CXXCircleNode extraNode(&aTorus.getCircle(), 0, newNode, 0);
		extraNode.setAngle(omega);
		addVertex(extraNode);
	}
}

void CXXSphereElement::addVertex(const CXXCircleNode &extraNode){
	double omega = extraNode.getAngle();
	const CXXCoord &newNode = extraNode.getCoord();
	const CXXCircle *aCircle = extraNode.getParent();
	
	list<TriangleEdgePair>::iterator cutTriangle;
	list<TriangleEdgePair>::iterator removeTriangle;
	int cutFound = 0;
	for (cutTriangle=edgeTriangles.begin(); 
		 cutTriangle!=edgeTriangles.end() && !cutFound; 
		 cutTriangle++){
		int iTriangle = (*cutTriangle).triangle;
		int iEdge = (*cutTriangle).edge;
		if (flatTriangles[iTriangle].doDraw()){
			if (flatTriangles[iTriangle].getEdgeCircle(iEdge) == aCircle){
				const CXXCircleNode &cn1 = flatTriangles[iTriangle].getCircleNode( iEdge     );
				double o1 = cn1.getAngle();
				const CXXCircleNode &cn2 = flatTriangles[iTriangle].getCircleNode((iEdge+1)%3);
				double o2 = cn2.getAngle();
				double omin, omax;
				if (o2>o1) {
					if (o2-o1 < M_PI){
						omin = o1;
						omax = o2;
					}
					else {
						omin = o2;
						omax = o1 + 2.*M_PI;
					}
				}
				else {
					if (o1-o2 < M_PI){
						omin = o2;
						omax = o1;
					}
					else {
						omin = o1;
						omax = o2 + 2.*M_PI;
					}
				}
				if (omega>omin && omega < omax){
					/*
					 cout << "select object banana add ball " << newNode[0] << " " << newNode[1] << " " << newNode[2] <<" 0.05\n";
					 cout << "#" <<omin <<" " <<omega<<" "<<omax<<endl;
					 CXXCoord projection(aTorus.coordFromThetaOmega(aTorus.getTheta2(),omega-aTorus.getAbsoluteStartOmega()));
					 cout << "apple " << omega <<" " << omin << " " << omax << endl; 
					 
					 cout << "select object banana add line " << newNode[0] << " " << newNode[1] << " " << newNode[2] << " "<<
					 projection[0] << " " << projection[1] << " " << projection[2] <<  "\n";
					 */
					
					cutFound = 1;
					removeTriangle = cutTriangle;
					int iv3 = addVertex(CXXSphereNode(newNode));
					
					flatTriangles.push_back(CXXSphereFlatTriangle(flatTriangles[iTriangle]));
					flatTriangles.back().setElement   ((iEdge+1)%3, iv3);
					flatTriangles.back().setCircleNode((iEdge+1)%3, extraNode);
					nDrawnTriangles++;
					edgeTriangles.push_back(TriangleEdgePair(flatTriangles.size()-1, iEdge));
					
					flatTriangles.push_back(CXXSphereFlatTriangle(flatTriangles[iTriangle]));
					flatTriangles.back().setElement   ( iEdge,      iv3);
					flatTriangles.back().setCircleNode( iEdge,      extraNode);
					nDrawnTriangles++;
					edgeTriangles.push_back(TriangleEdgePair(flatTriangles.size()-1, iEdge));
					
					hideFlatTriangle((*cutTriangle).triangle);
				}				
			}
		}
	}
	if (cutFound) edgeTriangles.erase(removeTriangle);
}

int CXXSphereElement::countDrawnTriangles() const{
	int nTri = nFlatTriangles();
	int nToDraw = 0;
	for (int i=0; i<nTri; i++){
		if (flatTriangles[i].doDraw()) nToDraw++;
	}
	return nToDraw;
}
