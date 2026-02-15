/*
 * MoleculesToTriangles/CXXClasses/NRStuff.h
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

#ifndef NRStuff_h
#define NRStuff_h

#include <vector>
#include <iostream>

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

class CoordSpline {
    std::vector<FCXXCoord> ctlPts;
    std::vector<FCXXCoord> spline;
    void DialASpline(float t, const std::vector<float> &a,  const std::vector<FCXXCoord> &p, int Cn, int interp, std::vector<FCXXCoord> &output, const int idx, std::vector<FCXXCoord> &work);
  public:
    CoordSpline(){
    };
    FCXXCoord  coordForXEquals(const float xVal) {
        int idx = int((xVal+1.0)/ctlPts.size() * spline.size()-1);
        int spline_size = spline.size();
        if (idx < 0)
            idx = 0;
        if (idx >= spline_size)
            idx = spline_size -1;
        return spline[idx];
    }
    void clearSpline(){
        ctlPts.clear();
        spline.clear();
    }
    void calculateYDoublePrimes(float ypLow, float ypHigh, int samplesPerSegment = 6) {
        spline = SplineCurve(ctlPts,(ctlPts.size()-1)*samplesPerSegment,3,1);
    }
    std::vector <FCXXCoord> SplineCurve(const std::vector<FCXXCoord> &ctlPts, int nsteps, int Cn, int iinterp);
    void addPair(float xVal, const FCXXCoord  &coord) {
        ctlPts.push_back(coord);
    }
};

#endif
