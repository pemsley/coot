/*
 *  CXXTorusElement.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXXTorusElement_included
#define CXXTorusElement_included
//#include <CXXTorusTriangle.h>
#include "CXXTorusNode.h"
#include "CXXTriangle.h"
#include "CXXCircleNode.h"

#include <vector>
#include <list>

class CXXSurface;
class CXXTriangle;
class CXXCircle;
#include "CXXCoord.h"

using std::vector;
 
class CXXTorusElement {
private:
	static CXXCircle nullCircle;
	const CXXCircle &theCircle;
	vector <CXXTorusNode> nodes;
	list <CXXTriangle  > flatTriangles;
	//omega1 and omega2 are the limits of the arc on the surface of a torus
	//	that is in the plane of the axis of the torus
	double omega1;
	double omega2;
	//When triangulated, the torus is divided up into 2**N (N>=0) in omega and
	//theta.  The corresponding angular step sizes are stored here
	int nOmegaSteps;
	double deltaOmega;
	int nThetaSteps;
	double deltaTheta;
	//absoluteStartOmega is the value of the omega angle of the start
	//of this arc, defined with respect to the referenceUnitVector of the
	//parent circle
	double absoluteStartOmega;
	//theta1 and theta2 are the limits of the arc on the surface of a torus
	//	that is in the plane of the axis of the torus
	double theta1;
	double theta2;
	//Here a set of triangles that point into the flatTriangles array to identify
	//triangles that constitute the edge strip (i.e. last step in theta
	list<CXXTriangle *> edgeTriangles;
	// v1unit and v2unit are unit vectors from the centre of the circle that defines
	//the trajectory of the probe around a torus, to the start and end point around that
	//orbit
	CXXCoord<CXXCoord_ftype>v1unit;
	CXXCoord<CXXCoord_ftype>v2unit;
	CXXCoord<CXXCoord_ftype>n1unit;
	//torusUnit and torusCentre are the axis and the centre of the torus respectively
	CXXCoord<CXXCoord_ftype>torusAxis;
	CXXCoord<CXXCoord_ftype>torusCentre;
	//rProbe and rTraj are the two radiuses that define the shape of the torus
	double rProbe, rTraj;
	void init();
	int debug;
	
public:
		CXXTorusElement();
	~CXXTorusElement();
	CXXTorusElement(const CXXCircle &aCircle, int iEdge, double delta, double probeRadius);
	void deleteLastTriangle(void);
	size_t addNode(CXXTorusNode &aNode);
	const size_t numberOfTorusNodes(void);
	const CXXCoord<CXXCoord_ftype>coordFromThetaOmega(double theta, double omega) const;
	const CXXCoord<CXXCoord_ftype>normalToProbeAtTheta(CXXCoord<CXXCoord_ftype>&p, double theta) const;
	const CXXCoord<CXXCoord_ftype>probeAtOmega(double omega) const;
	const CXXTorusNode &getNode(const int i) const;
	int upload(CXXSurface *aSurface);
	void addEdgeVertex(CXXCircleNode &aNode);
	int getNOmegaSteps() const {
		return nOmegaSteps;
	};
	double getDeltaOmega() const {
		return deltaOmega; 
	};
	double getAbsoluteStartOmega() const {
		return absoluteStartOmega;
	};
	const CXXCircle &getCircle() const {
		return theCircle;
	};
	double getTheta2() const {
		return theta2;
	};
    const CXXTorusNode &node(size_t i) const {
        return nodes[i];
    };
    size_t nTorusNodes() const {
        return nodes.size();
    };
    size_t nFlatTriangles() const {
        return flatTriangles.size();
    };
    list <CXXTriangle  >::const_iterator firstTriangle() const {
        return flatTriangles.begin();
    };
    list <CXXTriangle  >::const_iterator endOfTriangles() const {
        return flatTriangles.end();
    };
};

#endif
