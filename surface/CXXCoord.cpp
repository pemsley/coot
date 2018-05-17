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

int CXX_mot::CXXCoord::setXyz(const CXXCoord_ftype *x){
	for (int i=0; i<3; i++) xyzr[i]= x[i];
	xyzr[3] = 0.;
	return 0;
}
int CXX_mot::CXXCoord::getXyzr(CXXCoord_ftype *x) const{
	for (int i=0; i<4; i++) x[i]= xyzr[i];
	return 0;
}
int CXX_mot::CXXCoord::setXyzr(const CXXCoord_ftype *x){
	for (int i=0; i<4; i++) xyzr[i]= x[i];
	return 0;
}

CXX_mot::CXXCoord CXX_mot::CXXCoord::scalarMultiply(const CXXCoord_ftype factor){
	CXXCoord result(xyzr[0]*factor, xyzr[1]*factor, xyzr[2]*factor, xyzr[3]*factor);
	return result;
}

CXX_mot::CXXCoord_ftype CXX_mot::CXXCoord::get3DLength() {
	
	return sqrt(xyzr[0]*xyzr[0]+xyzr[1]*xyzr[1]+xyzr[2]*xyzr[2]);
	
}

CXX_mot::CXXCoord_ftype CXX_mot::CXXCoord::get3DLengthSq() {
		return xyzr[0]*xyzr[0]+xyzr[1]*xyzr[1]+xyzr[2]*xyzr[2];
}


int CXX_mot::CXXCoord::isZero() const{
	
   if ((xyzr[0] == 0) && (xyzr[1] == 0) && (xyzr[2] == 0))
      return 1;
   else
      return 0;
	
}

void CXX_mot::CXXCoord::normalise(){
	CXXCoord_ftype length = get3DLength();
	scale (1. / length);
}

CXX_mot::CXXCoord_ftype CXX_mot::CXXCoord::angleBetween(const CXXCoord &v1, const CXXCoord &v2) const{
	CXXCoord_ftype cosTheta = v1*v2;
	CXXCoord v3 = v1^v2;
	CXXCoord_ftype sinTheta = v3.get3DLength();
	CXXCoord_ftype result = atan2(sinTheta, cosTheta);
	if ( (*this)*v3< 0.) result *=-1.;
	while (result < 0.) result += 2.*M_PI;
	return result;
}

std::ostream &operator << ( std::ostream &out, const CXX_mot::CXXCoord &c )
{
    out << "[ ";
    for (int k=0; k< 4; k++) {
        out << c[k] << ' ';
    }
    out << ']';
	return out;
}
//CXXCoord_ftype CXXCoord::DEGTORAD = M_PI/180.;
//CXXCoord_ftype CXXCoord::RADTODEG = 180./M_PI;



