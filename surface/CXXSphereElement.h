/*
 *  CXXSphereElement.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXX_mot_CXXSphereElement_included
#define CXX_mot_CXXSphereElement_included
#include <vector>
#include <list> 
#include <iostream>
#include "CXXCoord.h"

using namespace std;

// try class CXXSphereTriangleEdge;
#include "CXXSphereTriangle.h"

// try class CXXSphereNode;
#include "CXXSphereFlatTriangle.h"
// try class CXXTriangle;

// try class CXXSurface;
#include "CXXSurface.h"

// try class CXXCircleNode;
#include "CXXCircle.h"

namespace CXX_mot {

  class CXXTorusElement;
  
  class TriangleEdgePair{
  public:
     int triangle;
     int edge;
     TriangleEdgePair(int i, int j) : triangle(i), edge(j) {}
  };

class CXXSphereElement{
private:
	const mmdb::Atom *theAtom;
	CXXCoord theCentre;
	vector<CXXSphereNode, CXX_old::CXXAlloc<CXXSphereNode> >theVertices;
	vector<CXXSphereTriangle, CXX_old::CXXAlloc<CXXSphereTriangle> >theTriangles;
	vector<CXXSphereTriangleEdge, CXX_old::CXXAlloc<CXXSphereTriangleEdge> >theEdges;
	list<CXXSphereFlatTriangle, CXX_old::CXXAlloc<CXXSphereFlatTriangle> >flatTriangles;
    vector<vector<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >, CXX_old::CXXAlloc<vector<CXXCircle, CXX_old::CXXAlloc<CXXCircle> > > >theCircles;
	map<CXXSphereFlatTriangle *, int> edgeTriangles;
	double theRadius;
	double deltaRadians;
	void init();
public: 
		CXXSphereElement();
//	~CXXSphereElement();
	
	//Constructor to fully triangulate a sphere element at a given coordinate and radius
	CXXSphereElement(const CXXCoord &position, double radius, double del);
	//Create a fully triangulated Sphere elment covering an atom using the above
	CXXSphereElement (mmdb::Atom *anAtom, double delta);

	//Create a triangulated sphere element from a centr, a radius, 
	// and three unit vectors pointing to the surface. Assume all angles
	// Less than 180 degrees
	//CXXSphereElement (CXXCoord centre, double radius, CXXCoord u1, CXXCoord u2, CXXCoord u3, double del);

	//Initialiser which uses the one above to create the surface patch that corresponds to a particular 
	//Node (i.e. point where probe is in contact with three atoms)
	void initWith(const CXXCircleNode &aNode, double del, double radius_in, bool *includeAtoms, int UseOrGenerate);
	//Copy constructor
	CXXSphereElement (const CXXSphereElement &oldOne);
	int addVertex(const CXXSphereNode &vert);
	int addTriangle(const CXXSphereTriangle &);
	int addEdge(const CXXSphereTriangleEdge &);
	int addFlatTriangle(const CXXSphereFlatTriangle &);
	
	void flattenLastTriangle(void);
	
	//Accessors to allow copy constructor
	const std::vector<CXXSphereNode, CXX_old::CXXAlloc<CXXSphereNode> > &getVertices() const {
		return theVertices;
	};
	const std::vector<CXXSphereTriangle, CXX_old::CXXAlloc<CXXSphereTriangle> > &getTriangles() const {
		return theTriangles;
	};
	const std::vector<CXXSphereTriangleEdge, CXX_old::CXXAlloc<CXXSphereTriangleEdge> > &getEdges() const{
		return theEdges;
	};
	const std::list<CXXSphereFlatTriangle, CXX_old::CXXAlloc<CXXSphereFlatTriangle> > &getFlatTriangles() const{
		return flatTriangles;
	};
	const vector<vector<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >, CXX_old::CXXAlloc<vector<CXXCircle, CXX_old::CXXAlloc<CXXCircle> > > > &getCircles() const{
		return theCircles;
	};
	const CXXCoord &centre() const;
	const unsigned nVertices() const;
	const CXXSphereNode &vertex(const int iVertex) const;
	void moveVertex(const int iVertex, const CXXCoord &position);
//	const int nTriangles() const ;
	const CXXSphereTriangle &triangle(const int iTriangle) const;
	const int nFlatTriangles() const;
//	const CXXSphereFlatTriangle &flatTriangle(const int iFlatTriangle) const;	
//	void hideFlatTriangle(const int iFlatTriangle);	
	const double radius() const;
	const int nEdges() const;
	const CXXSphereTriangleEdge &edge(const int iEdge) const;
	const mmdb::Atom *getAtom() const;
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
	int translateBy (const CXXCoord &crd);
	int scaleBy (const double factor);
	
	// Might wish to assign an atom to the surface
	int setAtom(const mmdb::Atom *anAtom);
	
	int addTriangularPatch(const CXXCoord &u1, 
						   const CXXCoord &u2, 
						   const CXXCoord &u3, mmdb::Atom *, 
						   vector<CXXCircle, CXX_old::CXXAlloc<CXXCircle> >&circles, 
						   int UseOrGenerate);
	
	static const int Reentrant = 1;
	static const int Contact = 2;
	static const int VDW = 3;
	static const int Accessible = 4;
	static const int GenerateCircles = 0;
	static const int UseCircles = 1;
	
	CXXCoord voronoiPoint(const CXXCoord &a, const CXXCoord &b, const CXXCoord &c) const;

	
	void identifyRaggedEdges(CXXCircle &theCircle, std::map<const CXXCircleNode*, vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > >&raggedEdges);
	
	void identifyRaggedEdges(CXXCircle &theCircle, std::map<const CXXBall*, vector<CXXCoord, CXX_old::CXXAlloc<CXXCoord> > >&raggedEdges);

	//flagCutTriangles: a stitching routine that inserts points around the edge of an intersection with a family of
	//segments of a torus
	
	int flagCutTriangles(const CXXCircle &aCircle);
	
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

	void initWith(const CXXCoord &aCentre, mmdb::PAtom atomI, mmdb::PAtom atomJ, mmdb::PAtom atomK, 
									double delta, double radius_in, const bool *includeAtoms);
		
};
}
#endif

