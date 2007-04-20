/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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
/*
 *  CXXSurfaceVertex.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  
 *
 */
#ifndef CXXSurfaceVertex_included
#define CXXSurfaceVertex_included
#include "CXXCoord.h"
#include <vector>
#include <deque>

using namespace std;

class CXXSurfaceVertex {
private:
	vector<void *> pointers;
	vector<CXXCoord> vectors;
	vector<double> scalars;
	void init();
public:
	CXXSurfaceVertex();
	~CXXSurfaceVertex();
	
	int setXyz(unsigned int coordType, double *xyz);
	int setXyzr(unsigned int coordType, double *xyzr);
	int setCoord(unsigned int coordType, const CXXCoord &crd);
	int setScalar(unsigned int scalarType, double value);
	int setPointer(unsigned int pointerType, void *value);

	double x(unsigned int coordType);
	double y(unsigned int coordType);
	double z(unsigned int coordType);
	double r(unsigned int coordType);
	double *xyzPntr(unsigned int coordType);
	int getXyz(unsigned int coordType, double *x);
	int getXyzr(unsigned int coordType, double *x);
	double scalar(unsigned int scalarType) const;
	void *pointer(unsigned int pointerType) const;
	const CXXCoord *coordPntr(unsigned int coordType) const;
	const CXXCoord &coordRef(unsigned int coordType) const;	
	int nPointers() const;
	int nVectors() const;
	int nScalars() const;
};
#endif




