/*
 *  CXXTorusNode.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <cstring>

#ifndef CXX_mot_CXXTorusNode_included
#define CXX_mot_CXXTorusNode_included
#include <mmdb2/mmdb_manager.h>

#include "CXXCoord.h"

// class mmdb::Atom;

namespace CXX_mot {

class CXXTorusNode{
private:
	mmdb::Atom *theAtom;
	double theta;
	double omega;
	CXXCoord crd;
	void init();
public:
	CXXTorusNode();
	CXXTorusNode(double the, double om);
	int setTheta(const double inTheta);
	int setOmega(const double inOmega);
	const CXXCoord &coord() const;
	const double getTheta() const;
	const double getOmega() const;
	int setCoord(const CXXCoord &);
	int setAtom(mmdb::Atom *anAtom);
	mmdb::Atom *getAtom() const;
};

}
#endif
