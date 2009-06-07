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
 *  CXXCoord.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Jan 24 2004.
 *  
 *
 */
#ifndef CXXCoord_included
#define CXXCoord_included
#include <iostream>


class CXXCoord {
private:
	double xyzr[4];
public:
	CXXCoord() {
	  for (int i=0; i<4; i++) xyzr[i] = 0.;
	}
	CXXCoord(const double *x){
          xyzr[0] = x[0]; xyzr[1]=x[1]; xyzr[2]=x[2]; xyzr[3]=x[3];
	}
	CXXCoord(const double x,const double y, const double z) {
	  xyzr[0]=x; xyzr[1]=y; xyzr[2]=z; xyzr[3]=0.;
	}	
	CXXCoord(const double x, const double y, const double z, const double r) {
	  xyzr[0]=x; xyzr[1]=y; xyzr[2]=z; xyzr[3]=r;
	}
	~CXXCoord(){}
	CXXCoord operator + (const CXXCoord &param) const{
		return CXXCoord(xyzr[0]+param[0], xyzr[1]+param[1], xyzr[2]+param[2], xyzr[3]+param[3]);
	}
	CXXCoord operator - (const CXXCoord &param) const{
		return CXXCoord(xyzr[0]-param[0], xyzr[1]-param[1], xyzr[2]-param[2], xyzr[3]-param[3]);
	}
	
	CXXCoord& operator *= (double factor){
		for (int i=0; i<4; i++){xyzr[i] = xyzr[i]*factor;}
	    return *this;
	}
	
	CXXCoord& operator /= (double invFactor){
		double factor = 1. / invFactor;
		for (int i=0; i<4; i++){xyzr[i] = xyzr[i]*factor;}
	    return *this;
	}
 	void dump() const{
		std::cout << xyzr[0] << " " << " " << xyzr[1]<< " " << xyzr[2]<< " " << xyzr[3] << std::endl;
	}
	double operator [] (int element) const{
		return xyzr[element];
	}
	int operator == (const CXXCoord &comparator) const{
		if (xyzr[0] != comparator[0]) return 0;
		if (xyzr[1] != comparator[1]) return 0;
		if (xyzr[2] != comparator[2]) return 0;
		return 1;
	}
	int operator == (CXXCoord &comparator) const{
		if (xyzr[0] != comparator[0]) return 0;
		if (xyzr[1] != comparator[1]) return 0;
		if (xyzr[2] != comparator[2]) return 0;
		return 1;
	}
	int operator != (const CXXCoord &comparator) const{
		if (xyzr[0] != comparator[0]) return 1;
		if (xyzr[1] != comparator[1]) return 1;
		if (xyzr[2] != comparator[2]) return 1;
		return 0;
	}
	double dot(const CXXCoord param) const{
		double result = 0.;
		for (int i=0; i<3; i++) result += xyzr[i]*param[i];
		return result;
	}
	CXXCoord operator ^ (const CXXCoord &param) const{
		CXXCoord temp( 
					   (y()*param.z() - param.y()*z()),
					   -(x()*param.z() - param.x()*z()) ,
					   (x()*param.y() - param.x()*y()) ,
					   0.
					   );
		return temp;
	}
	CXXCoord operator * (const double factor) const {
		CXXCoord temp(xyzr[0]*factor, xyzr[1]*factor, xyzr[2]*factor, xyzr[3]*factor);
		return temp;
	};
	CXXCoord scalarMultiply(const double factor);
	double x() const{
		return (xyzr[0]);
	}
	double y() const{
		return xyzr[1];
	}
	double z() const{
		return xyzr[2];
	}
	double r() const{
		return xyzr[3];
	}
	double element(int i) const{
		return xyzr[i];
	}
	double *xyzPntr(){
		return xyzr;
	}
	int isZero() const;
	int setX(const double x){
		xyzr[0]= x;
		return 0;
	}
	int setY(const double y){
		xyzr[1]= y;
		return 0;
	}
	int setZ(const double z){
		xyzr[2]= z;
		return 0;
	}
	int setR(const double r){
		xyzr[3]= r;
		return 0;
	}
	int getXyz(double *x) const{
		for (int i=0; i<3; i++) x[i]= xyzr[i];
		return 0;
	}
	int setXyz(const double *);
	int getXyzr(double *) const;
	int setXyzr(const double *);
	
	double angleBetween(const CXXCoord &v1, const CXXCoord &v2) const;
//Given two unit vectors v1 and v2, this evaluates the angle from v1 to v2 around the 
	//axis in a positive sense, returning it in the range 0 to 2PI

	double get3DLength();
	int scale(double factor){
		for (int i=0; i<3; i++) xyzr[i] *= factor;
		return 0;
	}
	double get3DLengthSq();
	void normalise();
};

#endif
