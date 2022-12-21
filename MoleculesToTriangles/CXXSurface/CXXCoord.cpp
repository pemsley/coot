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

// use clipper functions?
template<> double CXXCoord<double>::CXX_DEGTORAD = M_PI/180.;
template<> double CXXCoord<double>::CXX_RADTODEG = 180./M_PI;
template<> float CXXCoord<float>::CXX_DEGTORAD = M_PI/180.;
template<> float CXXCoord<float>::CXX_RADTODEG = 180./M_PI;


