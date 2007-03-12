/*
 *  CXXSphereElement.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXXSphereElement_included
#define CXXSphereElement_included
#include <vector>
#include <list> 
#include <iostream>
#include "CXXCoord.h"

using namespace std;
class CXXSphereTriangleEdge;
class CXXSphereTriangle;
class CXXSphereNode;
class CAtom;
class CXXSphereFlatTriangle;
class CXXTriangle;
class CXXSurface;
class CXXCircleNode;
class CXXCircle;
class CXXTorusElement;
#include "CXXPointyBit.h"

class TriangleEdgePair{
public:
	int triangle;
	int edge;
	TriangleEdgePair(int i, int j) : triangle(i), edge(j){
	};
};

class CXXSphereElement{
private:
	const CAtom *theAtom;
	CXXCoord theCentre;
	vector<CXXSphereNode>theVertices;
	vector<CXXSphereTriangle>theTriangles;
	vector<CXXSphereTriangleEdge>theEdges;
	vector<CXXSphereFlatTriangle>flatTriangles;
	list<TriangleEdgePair> edgeTriangles;
	double theRadius;
	double deltaRadians;
	void init();
public: 
		CXXSphereElement();
	~CXXSphereElement();
	
	//Constructor to fully triangulate a sphere element at a given coordinate and radius
	CXXSphereElement(const CXXCoord &position, double radius, double del);
	//Create a fully triangulated Sphere elment covering an atom using the above
	CXXSphereElement (CAtom *anAtom, double delta);

	//Create a triangulated sphere element from a centr, a radius, 
	// and three unit vectors pointing to the surface. Assume all angles
	// Less than 180 degrees
	//CXXSphereElement (CXXCoord centre, double radius, CXXCoord u1, CXXCoord u2, CXXCoord u3, double del);

	//Constructur which uses the one above to create the surface patch that corresponds to a particular 
	//Node (i.e. point where probe is in contact with three atoms)
	CXXSphereElement (const CXXCircleNode &aNode, double del, double radius_in, int selHnd, vector<PointyBit> &pointyBits);
	
	//Copy constructor
	CXXSphereElement (const CXXSphereElement &oldOne);
	int addVertex(const CXXSphereNode &vert);
	int addTriangle(const CXXSphereTriangle &);
	int addEdge(const CXXSphereTriangleEdge &);
	int addFlatTriangle(const CXXSphereFlatTriangle &);
	
	void flattenLastTriangle(void);
	
	//Accessors to allow copy constructor
	const std::vector<CXXSphereNode> &getVertices() const;
	const std::vector<CXXSphereTriangle> &getTriangles() const;
	const std::vector<CXXSphereTriangleEdge> &getEdges() const;
	const std::vector<CXXSphereFlatTriangle> &getFlatTriangles() const;
	
	const CXXCoord &centre() const;
	const unsigned nVertices() const;
	const CXXSphereNode &vertex(const int iVertex) const;
	void moveVertex(const int iVertex, const CXXCoord &position);
	const int nTriangles() const ;
	const CXXSphereTriangle &triangle(const int iTriangle) const;
	const int nFlatTriangles() const;
	const CXXSphereFlatTriangle &flatTriangle(const int iFlatTriangle) const;	
	void hideFlatTriangle(const int iFlatTriangle);	
	const double radius() const;
	const int nEdges() const;
	const CXXSphereTriangleEdge &edge(const int iEdge) const;
	const CAtom *getAtom() const;
	const double delta() const;

	//This will calculate a triangulated surface from current position, radius, and delta
	int calculate();
		
	//Method to delete Flat Triangles that are on the wrong side ofa circle.  The circle is brought closer
	//to the centre of the sphere by a factor (radius-probeRadius) / radius
	int trimBy(const CXXCircle &aCircle);

	//Method to add Flat triangles to the surface as appropriate.  Sense decides wheither we 
	// define normal for inside or outside of the sphere.  Accessible surface positions are
	//displaced from surface positions by a vector length probeRadius in the direction of the normal
	int upLoad(CXXSurface *aSurface, double probeRadius, const int sense) const;
	
	//Methods to translate and scale a sphere element
	int translateBy (const CXXCoord &crd);
	int scaleBy (const double factor);
	
	// Might wish to assign an atom to the surface
	int setAtom(const CAtom *anAtom);
	
	int addTriangularPatch(const CXXCoord &u1, 
						   const CXXCoord &u2, 
						   const CXXCoord &u3, CAtom *, 
						   vector <PointyBit> &pointyBits);
	
	static const int Inside = 1;
	static const int Outside = 2;

	CXXCoord voronoiPoint(const CXXCoord &a, const CXXCoord &b, const CXXCoord &c) const;

	//flagCutTriangles: a stitching routine that inserts points around the edge of an intersection with a family of
	//segments of a torus
	
	void flagCutTriangles(const CXXCircle &aCircle);
	
	//addTorusNodes: a stitching routine to add vertices that arise from an intersecting torus into the sphere surface
	
	void addTorusVertices(const CXXTorusElement &aTorus);

	void clearCutFlags(){
		edgeTriangles.empty();
	};
	
	
	int nDrawnTriangles;
	int getNDrawnTriangles() const {
		return nDrawnTriangles;
	};
	int countDrawnTriangles() const;
	
	void addVertex(const CXXCircleNode &aCircle);
};
#endif

