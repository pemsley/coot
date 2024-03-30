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
#include "mmdb2/mmdb_manager.h"

class CXXTorusNode{
private:
	mmdb::Atom* theAtom;
	double theta;
	double omega;
	CXXCoord<CXXCoord_ftype>crd;
	void init();
public:
		CXXTorusNode();
	CXXTorusNode(double the, double om);
	int setTheta(const double inTheta);
	int setOmega(const double inOmega);
	const CXXCoord<CXXCoord_ftype>&coord() const;
	const double getTheta() const;
	const double getOmega() const;
	int setCoord(const CXXCoord<CXXCoord_ftype>&);
	int setAtom(mmdb::Atom* anAtom);
	mmdb::Atom* getAtom() const;
};

#endif
