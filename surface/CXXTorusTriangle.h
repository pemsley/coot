/*
 *  CXXTorusTriangle.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXXTorusTriangle_included
#define CXXTorusTriangle_included
#include "CXXTorusNode.h"

class CXXTorusElement;

class CXXTorusTriangle {
private:
	CXXTorusElement *theTorusElement;
	int nodes[3];
	double deltaThetas[3];
	double deltaOmegas[3];
	void init();
public:
		CXXTorusTriangle();
	CXXTorusTriangle(const CXXTorusTriangle &target);
	CXXTorusTriangle(CXXTorusElement *aParent, 
					 int, int, int,
					 double, double, double, double, double, double);
	void setNode(const int i, int aNode);
	int getNode(const int i) const;
	CXXTorusElement *getTorusElement() const; 
	void setTorusElement(CXXTorusElement *anElement); 
	int bisect(double angleInDegrees);
	int tooLongEdge (double angle);
	double getDeltaTheta(const int iEdge) const;
	double getDeltaOmega(const int iEdge) const;
	void setDeltaTheta(const int iEdge, const double newTheta);
	void setDeltaOmega(const int iEdge, const double newOmega);
	void selfReport() const;
};

#endif
