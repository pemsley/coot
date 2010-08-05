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

class CAtom;

class PointyBit {
public:
	CAtom *atomI;
	CAtom *atomJ;
	CXXCoord coord;
	int   isNull;
	PointyBit() : atomI(0), atomJ(0), isNull(1) {};
        
	PointyBit(CAtom *_atomI, CAtom *_atomJ, const CXXCoord &_coord, int _isNull) :
	atomI(_atomI), atomJ(_atomJ), coord(_coord), isNull(_isNull) {};
};

#endif
