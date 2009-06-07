/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
 *  CXXCoord.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Jan 24 2004.
 *  
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
	
   if ((xyzr[0] == 0) && (xyzr[1] == 0) && (xyzr[2] == 0))
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
