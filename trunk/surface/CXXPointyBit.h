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

typedef struct PointyBit_{
	CAtom* atomI;
	CAtom* atomJ;
	CXXCoord coord;
	int   isNull;
} PointyBit;

#endif
