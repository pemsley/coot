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
#include "mmdb2/mmdb_manager.h"

using namespace std;
class CXXSphereTriangleEdge;
#include "CXXSphereTriangle.h"
class CXXSphereNode;
#include "CXXSphereFlatTriangle.h"
class CXXTriangle;
class CXXSurface;
class CXXCircleNode;
#include "CXXCircle.h"
class CXXTorusElement;

class TriangleEdgePair{
public:
	int triangle;
	int edge;
	TriangleEdgePair(int i, int j) : triangle(i), edge(j){
	};
};

class CXXSphereElement{
private:
	mmdb::Atom* theAtom;
	CXXCoord<CXXCoord_ftype>theCentre;
	vector<CXXSphereNode>theVertices;
	vector<CXXSphereTriangle  >theTriangles;
	vector<CXXSphereTriangleEdge>theEdges;
	list<CXXSphereFlatTriangle  >flatTriangles;
    vector<vector<CXXCircle> >theCircles;
	map<CXXSphereFlatTriangle *, int> edgeTriangles;
	double theRadius;
	double deltaRadians;
	void init();
public: 
		CXXSphereElement();
//	~CXXSphereElement();
	
	//Constructor to fully triangulate a sphere element at a given coordinate and radius
	CXXSphereElement(const CXXCoord<CXXCoord_ftype>&position, double radius, double del);
	//Create a fully triangulated Sphere elment covering an atom using the above
	CXXSphereElement (mmdb::Atom* anAtom, double delta);

	//Create a triangulated sphere element from a centr, a radius, 
	// and three unit vectors pointing to the surface. Assume all angles
	// Less than 180 degrees
	//CXXSphereElement (CXXCoord<CXXCoord_ftype>centre, double radius, CXXCoord<CXXCoord_ftype>u1, CXXCoord<CXXCoord_ftype>u2, CXXCoord<CXXCoord_ftype>u3, double del);

	//Initialiser which uses the one above to create the surface patch that corresponds to a particular 
	//Node (i.e. point where probe is in contact with three atoms)
	void initWith(const CXXCircleNode &aNode, double del, double radius_in, bool *includeAtoms, int UseOrGenerate);
	//Copy constructor
	CXXSphereElement (const CXXSphereElement &oldOne);
	size_t addVertex(const CXXSphereNode &vert);
	size_t addTriangle(const CXXSphereTriangle &);
	size_t addEdge(const CXXSphereTriangleEdge &);
	int addFlatTriangle(const CXXSphereFlatTriangle &);
	
	void flattenLastTriangle(void);
	
	//Accessors to allow copy constructor
	const std::vector<CXXSphereNode> &getVertices() const {
		return theVertices;
	};
	const std::vector<CXXSphereTriangle  > &getTriangles() const {
		return theTriangles;
	};
	const std::vector<CXXSphereTriangleEdge> &getEdges() const{
		return theEdges;
	};
	const std::list<CXXSphereFlatTriangle  > &getFlatTriangles() const{
		return flatTriangles;
	};
	const vector<vector<CXXCircle > > &getCircles() const{
		return theCircles;
	};
	const CXXCoord<CXXCoord_ftype>&centre() const;
	const size_t nVertices() const;
	const CXXSphereNode &vertex(const size_t iVertex) const;
	void moveVertex(const int iVertex, const CXXCoord<CXXCoord_ftype>&position);
//	const int nTriangles() const ;
	const CXXSphereTriangle &triangle(const size_t iTriangle) const;
	const size_t nFlatTriangles() const;
//	const CXXSphereFlatTriangle &flatTriangle(const int iFlatTriangle) const;	
//	void hideFlatTriangle(const int iFlatTriangle);	
	const double radius() const;
	const size_t nEdges() const;
	const CXXSphereTriangleEdge &edge(const size_t iEdge) const;
	mmdb::Atom* getAtom() const;
	const double delta() const;

	//This will calculate a triangulated surface from current position, radius, and delta
	int calculate();
		
	//Method to delete Flat Triangles that are on the wrong side ofa circle.  The circle is brought closer
	//to the centre of the sphere by a factor (radius-probeRadius) / radius
	int trimBy(const CXXCircle &aCircle, int carefully);

	//Method to add Flat triangles to the surface as appropriate.  Sense decides wheither we 
	// define normal for inside or outside of the sphere.  Accessible surface positions are
	//displaced from surface positions by a vector length probeRadius in the direction of the normal
	int upLoad(CXXSurface *aSurface, double probeRadius, const int sense) const;
	
	//Methods to translate and scale a sphere element
	int translateBy (const CXXCoord<CXXCoord_ftype>&crd);
	int scaleBy (const double factor);
	
	// Might wish to assign an atom to the surface
	int setAtom(mmdb::Atom* anAtom);
	
	int addTriangularPatch(const CXXCoord<CXXCoord_ftype>&u1, 
						   const CXXCoord<CXXCoord_ftype>&u2, 
						   const CXXCoord<CXXCoord_ftype>&u3, mmdb::Atom* , 
						   vector<CXXCircle  >&circles, 
						   int UseOrGenerate);
	
	static const int Reentrant = 1;
	static const int Contact = 2;
	static const int VDW = 3;
	static const int Accessible = 4;
	static const int GenerateCircles = 0;
	static const int UseCircles = 1;
	
	CXXCoord<CXXCoord_ftype>voronoiPoint(const CXXCoord<CXXCoord_ftype>&a, const CXXCoord<CXXCoord_ftype>&b, const CXXCoord<CXXCoord_ftype>&c) const;

	
	void identifyRaggedEdges(CXXCircle &theCircle, std::map<const CXXCircleNode*, vector<CXXCoord<CXXCoord_ftype> > >&raggedEdges);
	
	void identifyRaggedEdges(CXXCircle &theCircle, std::map<const CXXBall*, vector<CXXCoord<CXXCoord_ftype> > >&raggedEdges);

	//flagCutTriangles: a stitching routine that inserts points around the edge of an intersection with a family of
	//segments of a torus
	
	size_t flagCutTriangles(const CXXCircle &aCircle);
	
	//addTorusNodes: a stitching routine to add vertices that arise from an intersecting torus into the sphere surface
	
	void addTorusVertices(const CXXTorusElement &aTorus);
    void addCircleVertices(const CXXCircle &theCircle, int iEdge, double delta);

	void clearCutFlags(){
		edgeTriangles.clear();
	};
	
	
	int nDrawnTriangles;
	int getNDrawnTriangles() const {
		return nDrawnTriangles;
	};
	int countDrawnTriangles() const;
	
	int addVertex(const CXXCircleNode &aCircle);

	void initWith(const CXXCoord<CXXCoord_ftype>&aCentre, mmdb::Atom* atomI, mmdb::Atom* atomJ, mmdb::Atom* atomK, 
									double delta, double radius_in, const bool *includeAtoms);
		
};
#endif

