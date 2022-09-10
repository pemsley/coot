/*
 *  CXXSurfaceVertex.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
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




