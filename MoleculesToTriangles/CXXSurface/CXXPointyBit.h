/*
 *  CXXPointyBit.h
 *  CXXSurface
 *
 *  Created by Martin Noble on 30/05/2005.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXXPointyBit_included
#define CXXPointyBit_included
#include "CXXCoord.h"

class const mmdb::Atom;

class PointyBit {
public:
	mmdb::Atom* atomI;
	mmdb::Atom* atomJ;
	CXXCoord<CXXCoord_ftype>coord;
	int   isNull;
	PointyBit() : atomI(0), atomJ(0), isNull(1) {};
        
	PointyBit(mmdb::Atom* _atomI, mmdb::Atom* _atomJ, const CXXCoord<CXXCoord_ftype>&_coord, int _isNull) :
	atomI(_atomI), atomJ(_atomJ), coord(_coord), isNull(_isNull) {};
};

#endif
