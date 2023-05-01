/*
 *  CXXCoord.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Jan 24 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXX_mot_CXXCoord_included
#define CXX_mot_CXXCoord_included

#include <math.h>

	
#include <iostream>

namespace CXX_mot {

typedef double CXXCoord_ftype;

class CXXCoord {
public:
/* 	static CXXCoord_ftype DEGTORAD; */
/* 	static CXXCoord_ftype RADTODEG; */
/* 	static double DEGTORAD; */
/* 	static double RADTODEG; */
	
	CXXCoord_ftype xyzr[4];
	CXXCoord() {
		for (int i=0; i<4; i++){xyzr[i] = 0.;}
    };
	CXXCoord(const CXXCoord_ftype *x){
          xyzr[0] = x[0]; xyzr[1]=x[1]; xyzr[2]=x[2]; xyzr[3]=x[3];
	};
	CXXCoord(const CXXCoord_ftype x,const CXXCoord_ftype y, const CXXCoord_ftype z) {
	  xyzr[0]=x; xyzr[1]=y; xyzr[2]=z; xyzr[3]=0.;
	};	
	CXXCoord(const CXXCoord_ftype x, const CXXCoord_ftype y, const CXXCoord_ftype z, const CXXCoord_ftype r) {
	  xyzr[0]=x; xyzr[1]=y; xyzr[2]=z; xyzr[3]=r;
	};
	CXXCoord operator + (const CXXCoord &param) const{
		return CXXCoord(xyzr[0]+param[0], xyzr[1]+param[1], xyzr[2]+param[2], xyzr[3]+param[3]);
	};
	CXXCoord operator - (const CXXCoord &param) const{
		CXXCoord_ftype result[4];
		for (int i=0; i<4; i++) result[i] = xyzr[i]-param[i];
		return CXXCoord(result);
	};
	
	CXXCoord& operator *= (CXXCoord_ftype factor){
	for (int i=0; i<4; i++) xyzr[i] *= factor;
	    return *this;
	};
	
	CXXCoord& operator /= (CXXCoord_ftype invFactor){
		CXXCoord_ftype factor = 1. / invFactor;
		for (int i=0; i<4; i++){xyzr[i] = xyzr[i]*factor;}
	    return *this;
	};
    
    CXXCoord operator / (CXXCoord_ftype invFactor){
        CXXCoord result(*this);
        result /= invFactor;
        return result;
    };
        
    CXXCoord operator * (CXXCoord_ftype factor){
        CXXCoord result(*this);
        result *= factor;
        return result;
    };
    
 	void dump() const{
		std::cout << xyzr[0] << " " << " " << xyzr[1]<< " " << xyzr[2]<< " " << xyzr[3] << std::endl;
	};
    const CXXCoord_ftype& operator [] (unsigned element) const {
		return xyzr[element];
	};
    CXXCoord_ftype& operator [] (unsigned element) {
		return xyzr[element];
	};
	int operator == (const CXXCoord &comparator) const{
		if (xyzr[0] != comparator[0]) return 0;
		if (xyzr[1] != comparator[1]) return 0;
		if (xyzr[2] != comparator[2]) return 0;
		return 1;
	};
	int operator == (CXXCoord &comparator) const{
		if (xyzr[0] != comparator[0]) return 0;
		if (xyzr[1] != comparator[1]) return 0;
		if (xyzr[2] != comparator[2]) return 0;
		return 1;
	};
	int operator != (const CXXCoord &comparator) const{
		if (xyzr[0] != comparator[0]) return 1;
		if (xyzr[1] != comparator[1]) return 1;
		if (xyzr[2] != comparator[2]) return 1;
		return 0;
	};
	CXXCoord operator ^ (const CXXCoord &param) const{
		CXXCoord temp( 
					   (y()*param.z() - param.y()*z()),
					   -(x()*param.z() - param.x()*z()) ,
					   (x()*param.y() - param.x()*y()) ,
					   0.
					   );
		return temp;
	};
	CXXCoord_ftype operator * (const CXXCoord &other) const {
		CXXCoord_ftype temp(xyzr[0]*other[0] + xyzr[1]*other[1] + xyzr[2]*other[2] + xyzr[3]*other[3]);
		return temp;
	};
	CXXCoord operator * (const CXXCoord_ftype factor) const {
		CXXCoord temp(xyzr[0]*factor, xyzr[1]*factor, xyzr[2]*factor, xyzr[3]*factor);
		return temp;
	};
	CXXCoord scalarMultiply(const CXXCoord_ftype factor);
	CXXCoord_ftype x() const{
		return (xyzr[0]);
	};
	CXXCoord_ftype y() const{
		return xyzr[1];
	};
	CXXCoord_ftype z() const{
		return xyzr[2];
	};
	CXXCoord_ftype r() const{
		return xyzr[3];
	};
	CXXCoord_ftype element(int i) const{
		return xyzr[i];
	};
	CXXCoord_ftype *xyzPntr(){
		return xyzr;
	};
	int isZero() const;
	int setX(const CXXCoord_ftype x){
		xyzr[0]= x;
		return 0;
	};
	int setY(const CXXCoord_ftype y){
		xyzr[1]= y;
		return 0;
	};
	int setZ(const CXXCoord_ftype z){
		xyzr[2]= z;
		return 0;
	};
	int setR(const CXXCoord_ftype r){
		xyzr[3]= r;
		return 0;
	};
	int getXyz(CXXCoord_ftype *x) const{
		for (int i=0; i<3; i++) x[i]= xyzr[i];
		return 0;
	};
	int setXyz(const CXXCoord_ftype *);
	int getXyzr(CXXCoord_ftype *) const;
	int setXyzr(const CXXCoord_ftype *);
	
	CXXCoord_ftype angleBetween(const CXXCoord &v1, const CXXCoord &v2) const;
//Given two unit vectors v1 and v2, this evaluates the angle from v1 to v2 around the 
	//axis in a positive sense, returning it in the range 0 to 2PI

	CXXCoord_ftype get3DLength();
	int scale(CXXCoord_ftype factor){
		for (int i=0; i<3; i++) xyzr[i] *= factor;
		return 0;
	};
	CXXCoord_ftype get3DLengthSq();
	void normalise();
    
    friend std::ostream &operator << ( std::ostream &out, const CXXCoord &c );
    
    bool isNearly(const CXXCoord &compare, double tolerance) const {
		CXXCoord diff(*this-compare);
        if (fabs(float(diff[0]))>tolerance) return false;
        if (fabs(float(diff[1]))>tolerance) return false;
        if (fabs(float(diff[2]))>tolerance) return false;
        if (fabs(float(diff[3]))>tolerance) return false;
        return true;
    };
};

}
#endif
