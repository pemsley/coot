/*
 *  CXXCoord.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Sat Jan 24 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CXXCoord_included
#define CXXCoord_included

#include <math.h>
#include <iostream>

typedef double CXXCoord_ftype;

template <typename T>
class CXXCoord{
public:
    static T CXX_DEGTORAD;
    static T CXX_RADTODEG;
    
    T xyzr[4];
    CXXCoord<T>() {
#pragma omp simd
        for (int i=0; i<4; i++){xyzr[i] = 0.;}
    };
    CXXCoord<T>(const T *x){
#pragma omp simd
        for (int i=0; i<4; i++){
            xyzr[i] = x[i];
        }
    };
    CXXCoord<T>(const T x,const T y, const T z) {
        xyzr[0]=x; xyzr[1]=y; xyzr[2]=z; xyzr[3]=0.;
    };
    CXXCoord<T>(const T x, const T y, const T z, const T r) {
        xyzr[0]=x; xyzr[1]=y; xyzr[2]=z; xyzr[3]=r;
    };
    CXXCoord<T>operator + (const CXXCoord<T>&param) const{
        T sum[4];
#pragma omp simd
        for (int i=0; i<4; i++){
            sum[i] = xyzr[i]+param[i];
        }
        return CXXCoord<T>(sum);
    };
    CXXCoord<T>operator - (const CXXCoord<T>&param) const{
        T diff[4];
#pragma omp simd
        for (int i=0; i<4; i++){
            diff[i] = xyzr[i]-param[i];
        }
        return CXXCoord<T>(diff);
    };
    
    CXXCoord& operator *= (T factor){
#pragma omp simd
        for (int i=0; i<4; i++) xyzr[i] *= factor;
        return *this;
    };
    
    CXXCoord& operator /= (T invFactor){
        T factor = 1. / invFactor;
        for (int i=0; i<4; i++){xyzr[i] = xyzr[i]*factor;}
        return *this;
    };
    
    CXXCoord<T>operator / (T invFactor){
        CXXCoord<T>result(*this);
        result /= invFactor;
        return result;
    };
    
    CXXCoord<T>operator * (T factor){
        CXXCoord<T>result(*this);
        result *= factor;
        return result;
    };
    
    void dump() const{
        std::cout << xyzr[0] << " " << " " << xyzr[1]<< " " << xyzr[2]<< " " << xyzr[3] << std::endl;
    };
    const T& operator [] (unsigned element) const {
        return xyzr[element];
    };
    T& operator [] (unsigned element) {
        return xyzr[element];
    };
    int operator == (const CXXCoord<T>&comparator) const{
        if (xyzr[0] != comparator[0]) return 0;
        if (xyzr[1] != comparator[1]) return 0;
        if (xyzr[2] != comparator[2]) return 0;
        return 1;
    };
    int operator == (CXXCoord<T>&comparator) const{
        if (xyzr[0] != comparator[0]) return 0;
        if (xyzr[1] != comparator[1]) return 0;
        if (xyzr[2] != comparator[2]) return 0;
        return 1;
    };
    int operator != (const CXXCoord<T>&comparator) const{
        if (xyzr[0] != comparator[0]) return 1;
        if (xyzr[1] != comparator[1]) return 1;
        if (xyzr[2] != comparator[2]) return 1;
        return 0;
    };
    CXXCoord<T>operator ^ (const CXXCoord<T>&param) const{
        CXXCoord<T>temp(
                            (y()*param.z() - param.y()*z()),
                            -(x()*param.z() - param.x()*z()) ,
                            (x()*param.y() - param.x()*y()) ,
                            0.
                            );
        return temp;
    };
    T operator * (const CXXCoord<T>&other) const {
        T result = 0.;
        for (int i=0; i<4; i++){
            result += xyzr[i]*other[i];
        }
        return result;
    };
    CXXCoord<T> & operator+=(const CXXCoord<T> &other){
        for (int i=0; i<4; i++){
            (*this)[i] += other[i];
        }
        return *this;
    };
    CXXCoord<T> & operator-=(const CXXCoord<T> &other){
        for (int i=0; i<4; i++){
            (*this)[i] -= other[i];
        }
        return *this;
    };
    CXXCoord<T>operator * (const T factor) const {
        CXXCoord<T>temp(xyzr[0]*factor, xyzr[1]*factor, xyzr[2]*factor, xyzr[3]*factor);
        return temp;
    };
    T x() const{
        return (xyzr[0]);
    };
    T y() const{
        return xyzr[1];
    };
    T z() const{
        return xyzr[2];
    };
    T r() const{
        return xyzr[3];
    };
    T element(int i) const{
        return xyzr[i];
    };
    T *xyzPntr(){
        return xyzr;
    };
    int setX(const T x){
        xyzr[0]= x;
        return 0;
    };
    int setY(const T y){
        xyzr[1]= y;
        return 0;
    };
    int setZ(const T z){
        xyzr[2]= z;
        return 0;
    };
    int setR(const T r){
        xyzr[3]= r;
        return 0;
    };
    int getXyz(T *x) const{
        for (int i=0; i<3; i++) x[i]= xyzr[i];
        return 0;
    };
    
    
    T get3DLength() const {return sqrt(get3DLengthSq());};
    int scale(T factor){
        for (int i=0; i<3; i++) xyzr[i] *= factor;
        return 0;
    };
    T get3DLengthSq() const {return *this * *this;};
    
    bool isNearly(const CXXCoord<T>&compare, double tolerance) const {
        CXXCoord<T>diff(*this-compare);
        if (fabs(float(diff[0]))>tolerance) return false;
        if (fabs(float(diff[1]))>tolerance) return false;
        if (fabs(float(diff[2]))>tolerance) return false;
        if (fabs(float(diff[3]))>tolerance) return false;
        return true;
    };
    
    int setXyz(const T *x){
        for (int i=0; i<3; i++) xyzr[i]= x[i];
        xyzr[3] = 0.;
        return 0;
    }
    int getXyzr(T *x) const{
        for (int i=0; i<4; i++) x[i]= xyzr[i];
        return 0;
    }
    int setXyzr(const T *x){
        for (int i=0; i<4; i++) xyzr[i]= x[i];
        return 0;
    }
    
    CXXCoord<T> scalarMultiply(const T factor){
        CXXCoord<T> result(xyzr[0]*factor, xyzr[1]*factor, xyzr[2]*factor, xyzr[3]*factor);
        return result;
    }
    
    int isZero() const{
        
        if (xyzr[0] == xyzr[1] == xyzr[2] == 0)
        return 1;
        else
        return 0;
        
    }
    
    void normalise(){
        T length = get3DLength();
        scale (1. / length);
    }
    
    T angleBetween(const CXXCoord<T>&v1, const CXXCoord<T>&v2) const{
        //Given two unit vectors v1 and v2, this evaluates the angle from v1 to v2 around the
        //axis in a positive sense, returning it in the range 0 to 2PI
        T cosTheta = v1*v2;
        CXXCoord<T>v3 = v1^v2;
        T sinTheta = v3.get3DLength();
        T result = atan2(sinTheta, cosTheta);
        if ( (*this)*v3< 0.) result *=-1.;
        while (result < 0.) result += 2.*M_PI;
        return result;
    }
    
    friend std::ostream &operator << ( std::ostream &out, const CXXCoord<T>&c )
    {
        out << "[ ";
        for (int k=0; k< 4; k++) {
            out << c[k] << ' ';
        }
        out << ']';
        return out;
    }
    
    
    
    
};

typedef CXXCoord<double> DCXXCoord;
typedef CXXCoord<float> FCXXCoord;

#endif
