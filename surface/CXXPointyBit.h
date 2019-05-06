/*
 *  CXXPointyBit.h
 *  CXXSurface
 *
 *  Created by Martin Noble on 30/05/2005.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXX_mot_CXXPointyBit_included
#define CXX_mot_CXXPointyBit_included
#include "CXXCoord.h"

class mmdb::Atom;

class PointyBit {
public:
	mmdb::Atom *atomI;
	mmdb::Atom *atomJ;
	CXXCoord coord;
	int   isNull;
	PointyBit() : atomI(0), atomJ(0), isNull(1) {};
        
	PointyBit(mmdb::Atom *_atomI, mmdb::Atom *_atomJ, const CXXCoord &_coord, int _isNull) :
	atomI(_atomI), atomJ(_atomJ), coord(_coord), isNull(_isNull) {};
};

#endif
