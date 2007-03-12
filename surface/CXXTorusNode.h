/*
 *  CXXTorusNode.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXXTorusNode_included
#define CXXTorusNode_included
#include "CXXCoord.h"

class CAtom;

class CXXTorusNode{
private:
	CAtom *theAtom;
	double theta;
	double omega;
	CXXCoord crd;
	int ind;
	void init();
public:
		CXXTorusNode();
	~CXXTorusNode();
	CXXTorusNode(double the, double om);
	int setTheta(const double inTheta);
	int setOmega(const double inOmega);
	const CXXCoord &coord() const;
	const double getTheta() const;
	const double getOmega() const;
	int setIndex(int i);
	const int index(void) const;
	int setCoord(const CXXCoord &);
	int setAtom(CAtom *anAtom);
	CAtom *getAtom() const;
	CXXTorusNode(const CXXTorusNode &oldOne){
	  theAtom = oldOne.getAtom();
	  crd = oldOne.coord();
	  omega = oldOne.getOmega();
	  theta = oldOne.getTheta();
	  ind = oldOne.index();
	}
};

#endif
