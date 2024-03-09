/*
 * MoleculesToTriangles/CXXClasses/NRStuff.h
 *
 * Copyright 2009 by Martin Noble, University of Oxford
 * Author: Martin Noble
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
