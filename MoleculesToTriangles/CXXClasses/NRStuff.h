/*
 *  NRStuff.h
 *  MMDBRibbons
 *
 *  Created by Martin Noble on 17/07/2008.
 *  Copyright 2008 LMB, Oxford University. All rights reserved.
 *
 */

#ifndef NRStuff_h
#define NRStuff_h

#include <vector>
#include <iostream>
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

class NRSpline {
private:
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> y2;
    bool yDoublePrimeCalculated;
public:
    NRSpline() : yDoublePrimeCalculated(false) {};
	void clearSpline(){
        x.clear();
        y.clear();
        y2.clear();
	};
	~NRSpline(){
		clearSpline();
    }
    void addPair(float xVal, float yVal){
        x.push_back(xVal);
        y.push_back(yVal);
    }
    void calculateYDoublePrime(float ypLow, float ypHigh);
    float yForXEquals(const float xVal) ;
};

class CoordSpline {
public:
    NRSpline xSpline;
    NRSpline ySpline;
    NRSpline zSpline;
    NRSpline rSpline;
public:
    CoordSpline(){};
    void addPair(float xVal, const FCXXCoord  &coord) {
        xSpline.addPair(xVal, (float)coord[0]);
        ySpline.addPair(xVal, (float)coord[1]);
        zSpline.addPair(xVal, (float)coord[2]);
        rSpline.addPair(xVal, (float)coord[3]);
    };
    void calculateYDoublePrimes(float ypLow, float ypHigh) {
        xSpline.calculateYDoublePrime(ypLow, ypHigh);
        ySpline.calculateYDoublePrime(ypLow, ypHigh);
        zSpline.calculateYDoublePrime(ypLow, ypHigh);
        rSpline.calculateYDoublePrime(ypLow, ypHigh);
    };
    FCXXCoord  coordForXEquals(const float xVal) {
        FCXXCoord  result;
        result[0] = xSpline.yForXEquals(xVal);
        result[1] = ySpline.yForXEquals(xVal);
        result[2] = zSpline.yForXEquals(xVal);
        result[3] = rSpline.yForXEquals(xVal);
        return result;
    };
	void clearSpline(){
		xSpline.clearSpline();
		ySpline.clearSpline();
		zSpline.clearSpline();
		rSpline.clearSpline();
	};
};

#endif
