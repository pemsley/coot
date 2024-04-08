/*
 * MoleculesToTriangles/CXXSurface/CXXSurfaceVertex.h
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
#ifndef CXXSurfaceVertex_included
#define CXXSurfaceVertex_included
#include "CXXCoord.h"
#include <vector>
#include <deque>

using namespace std;

class CXXSurfaceVertex {
private:
	vector<void *> pointers;
	vector<CXXCoord<CXXCoord_ftype>  > vectors;
	vector<double> scalars;
	void init();
public:
	CXXSurfaceVertex() {init();};
	
	int setXyz(size_t coordType, double *xyz);
	int setXyzr(size_t coordType, double *xyzr);
	int setCoord(size_t coordType, const CXXCoord<CXXCoord_ftype> &crd);
	int setScalar(size_t scalarType, double value);
	int setPointer(size_t pointerType, void *value);

	double x(size_t coordType);
	double y(size_t coordType);
	double z(size_t coordType);
	double r(size_t coordType);
	CXXCoord_ftype *xyzPntr(size_t coordType);
	int getXyz(size_t coordType, double *x);
	int getXyzr(size_t coordType, double *x);
	double scalar(size_t scalarType) const;
	void *pointer(size_t pointerType) const;
	const CXXCoord<CXXCoord_ftype> *coordPntr(size_t coordType) const;
	const CXXCoord<CXXCoord_ftype> &coordRef(size_t coordType) const;
	size_t nPointers() const;
	size_t nVectors() const;
	size_t nScalars() const;
};
#endif




