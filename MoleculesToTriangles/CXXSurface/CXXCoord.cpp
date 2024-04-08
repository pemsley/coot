/*
 * MoleculesToTriangles/CXXSurface/CXXCoord.cpp
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
#define _USE_MATH_DEFINES
#include "CXXCoord.h"
#include <math.h>

// use clipper functions?
template<> double CXXCoord<double>::CXX_DEGTORAD = M_PI/180.;
template<> double CXXCoord<double>::CXX_RADTODEG = 180./M_PI;
template<> float CXXCoord<float>::CXX_DEGTORAD = M_PI/180.;
template<> float CXXCoord<float>::CXX_RADTODEG = 180./M_PI;


