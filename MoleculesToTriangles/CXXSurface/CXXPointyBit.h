/*
 * MoleculesToTriangles/CXXSurface/CXXPointyBit.h
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
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
