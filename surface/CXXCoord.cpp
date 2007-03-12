/*
 *  CXXCoord.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Jan 24 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXCoord.h"
#include <math.h>

int CXXCoord::setXyz(const double *x){
	for (int i=0; i<3; i++) xyzr[i]= x[i];
	xyzr[3] = 0.;
	return 0;
}
int CXXCoord::getXyzr(double *x) const{
	for (int i=0; i<4; i++) x[i]= xyzr[i];
	return 0;
}
int CXXCoord::setXyzr(const double *x){
	for (int i=0; i<4; i++) xyzr[i]= x[i];
	return 0;
}

CXXCoord CXXCoord::scalarMultiply(const double factor){
	CXXCoord result(xyzr[0]*factor, xyzr[1]*factor, xyzr[2]*factor, xyzr[3]*factor);
	return result;
}

double CXXCoord::get3DLength() {
	
	return sqrt(xyzr[0]*xyzr[0]+xyzr[1]*xyzr[1]+xyzr[2]*xyzr[2]);
	
}

double CXXCoord::get3DLengthSq() {
	
	return xyzr[0]*xyzr[0]+xyzr[1]*xyzr[1]+xyzr[2]*xyzr[2];
	
}


int CXXCoord::isZero() const{
	
	if (xyzr[0] == xyzr[1] == xyzr[2] == 0)
		return 1;
	else
		return 0;
	
}

void CXXCoord::normalise(){
	double length = get3DLength();
	scale (1. / length);
}

double CXXCoord::angleBetween(const CXXCoord &v1, const CXXCoord &v2) const{
	double cosTheta = v1.dot(v2);
	CXXCoord v3 = v1^v2;
	double sinTheta = v3.get3DLength();
	double result = atan2(sinTheta, cosTheta);
	while (result < 0.) result += 2.*M_PI;
	if (dot(v3)< 0.) result = 2.*M_PI - result;
	return result;
}
