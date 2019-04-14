/* density-contour/CIsoSurface.cpp
 * 
 * Copyright 2000 Paul Bourke 
 * Copyright 2000 Cory Gene Bloyd
 * Copyright 2000 Raghavendra Chandrashekara
 * Copyright 2005 The University of York
 * 
 * Author: Raghavendra Chandrashekara, Paul Bourke and Cory Gene Bloyd
 *         Paul Emsley and Kevin Cowtan
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
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

//
// File Name: CIsoSurface.cpp
// Last Modified: 5/8/2000
// Author: Raghavendra Chandrashekara (based on source code provided
// by Paul Bourke and Cory Gene Bloyd)
// Email: rc99@doc.ic.ac.uk, rchandrashekara@hotmail.com
//
// Description: This is the implementation file for the CIsoSurface class.
//
// Extra crystallograhic and export code added by Paul Emsley and
// Kevin Cowtan.
//

#include <stdlib.h>
#include <iomanip>
#include <fstream>

// #include "stdafx.h" commented by PE
#include <math.h>

#ifdef HAVE_CXX_THREAD
#include <thread>
#endif

// #include "clipper/mtz/mtz_io.h" // no need for this
//#include "map_io.h"


#include "CIsoSurface.h"

template <class T> const unsigned int CIsoSurface<T>::m_edgeTable[256] = {
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

// template <class T> const unsigned int CIsoSurface<T>::m_triTable[256][16] = {
template <class T> const int CIsoSurface<T>::m_triTable[256][16] = {
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
	{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
	{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
	{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
	{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
	{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
	{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
	{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
	{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
	{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
	{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
	{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
	{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
	{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
	{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
	{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
	{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
	{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
	{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
	{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
	{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
	{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
	{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
	{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
	{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
	{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
	{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
	{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
	{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
	{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
	{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
	{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
	{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
	{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
	{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
	{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
	{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
	{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
	{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
	{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
	{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
	{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
	{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
	{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
	{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
	{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};

template <class T> CIsoSurface<T>::CIsoSurface()
{
	m_fCellLengthX = 0;
	m_fCellLengthY = 0;
	m_fCellLengthZ = 0;
	m_nCellsX = 0;
	m_nCellsY = 0;
	m_nCellsZ = 0;
	m_nTriangles = 0;
	m_nNormals = 0;
	m_nVertices = 0;
	m_ppt3dVertices = NULL;
	m_piTriangleIndices = NULL;
	m_pvec3dNormals = NULL;
	m_ptScalarField = NULL;
	m_tIsoLevel = 0;
	m_bValidSurface = false;
}

template <class T> CIsoSurface<T>::~CIsoSurface()
{
	DeleteSurface();
}

template <class T> void CIsoSurface<T>::GenerateSurface(const T* ptScalarField, T tIsoLevel, unsigned int nCellsX, unsigned int nCellsY, unsigned int nCellsZ, float fCellLengthX, float fCellLengthY, float fCellLengthZ)
{

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_0 = std::chrono::high_resolution_clock::now();
#endif

   if (m_bValidSurface)
      DeleteSurface();

   m_tIsoLevel = tIsoLevel;
   m_nCellsX = nCellsX;
   m_nCellsY = nCellsY;
   m_nCellsZ = nCellsZ;
   m_fCellLengthX = fCellLengthX;
   m_fCellLengthY = fCellLengthY;
   m_fCellLengthZ = fCellLengthZ;
   m_ptScalarField = ptScalarField;

   unsigned int nPointsInXDirection = (m_nCellsX + 1);
   unsigned int nPointsInSlice = nPointsInXDirection*(m_nCellsY + 1);

   // Generate isosurface.
   for (unsigned int z = 0; z < m_nCellsZ; z++)
      for (unsigned int y = 0; y < m_nCellsY; y++)
	 for (unsigned int x = 0; x < m_nCellsX; x++) {
	    // Calculate table lookup index from those
	    // vertices which are below the isolevel.
	    unsigned int tableIndex = 0;
	    if (m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] < m_tIsoLevel)
	       tableIndex |= 1;
	    if (m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x] < m_tIsoLevel)
	       tableIndex |= 2;
	    if (m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + (x+1)] < m_tIsoLevel)
	       tableIndex |= 4;
	    if (m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)] < m_tIsoLevel)
	       tableIndex |= 8;
	    if (m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x] < m_tIsoLevel)
	       tableIndex |= 16;
	    if (m_ptScalarField[(z+1)*nPointsInSlice + (y+1)*nPointsInXDirection + x] < m_tIsoLevel)
	       tableIndex |= 32;
	    if (m_ptScalarField[(z+1)*nPointsInSlice + (y+1)*nPointsInXDirection + (x+1)] < m_tIsoLevel)
	       tableIndex |= 64;
	    if (m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + (x+1)] < m_tIsoLevel)
	       tableIndex |= 128;

	    // Now create a triangulation of the isosurface in this
	    // cell.
	    if (m_edgeTable[tableIndex] != 0) {
	       if (m_edgeTable[tableIndex] & 8) {
		  POINT3DID pt = CalculateIntersection(x, y, z, 3);
		  unsigned int id = GetEdgeID(x, y, z, 3);
		  m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
	       }
	       if (m_edgeTable[tableIndex] & 1) {
		  POINT3DID pt = CalculateIntersection(x, y, z, 0);
		  unsigned int id = GetEdgeID(x, y, z, 0);
		  m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
	       }
	       if (m_edgeTable[tableIndex] & 256) {
		  POINT3DID pt = CalculateIntersection(x, y, z, 8);
		  unsigned int id = GetEdgeID(x, y, z, 8);
		  m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
	       }
					
	       if (x == m_nCellsX - 1) {
		  if (m_edgeTable[tableIndex] & 4) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 2);
		     unsigned int id = GetEdgeID(x, y, z, 2);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
		  if (m_edgeTable[tableIndex] & 2048) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 11);
		     unsigned int id = GetEdgeID(x, y, z, 11);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
	       }
	       if (y == m_nCellsY - 1) {
		  if (m_edgeTable[tableIndex] & 2) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 1);
		     unsigned int id = GetEdgeID(x, y, z, 1);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
		  if (m_edgeTable[tableIndex] & 512) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 9);
		     unsigned int id = GetEdgeID(x, y, z, 9);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
	       }
	       if (z == m_nCellsZ - 1) {
		  if (m_edgeTable[tableIndex] & 16) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 4);
		     unsigned int id = GetEdgeID(x, y, z, 4);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
		  if (m_edgeTable[tableIndex] & 128) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 7);
		     unsigned int id = GetEdgeID(x, y, z, 7);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
	       }
	       if ((x==m_nCellsX - 1) && (y==m_nCellsY - 1))
		  if (m_edgeTable[tableIndex] & 1024) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 10);
		     unsigned int id = GetEdgeID(x, y, z, 10);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
	       if ((x==m_nCellsX - 1) && (z==m_nCellsZ - 1))
		  if (m_edgeTable[tableIndex] & 64) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 6);
		     unsigned int id = GetEdgeID(x, y, z, 6);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
	       if ((y==m_nCellsY - 1) && (z==m_nCellsZ - 1))
		  if (m_edgeTable[tableIndex] & 32) {
		     POINT3DID pt = CalculateIntersection(x, y, z, 5);
		     unsigned int id = GetEdgeID(x, y, z, 5);
		     m_i2pt3idVertices.insert(ID2POINT3DID::value_type(id, pt));
		  }
					
	       for (unsigned int i = 0; m_triTable[tableIndex][i] != -1; i += 3) {
		  TRIANGLE triangle;
		  unsigned int pointID0, pointID1, pointID2;
		  pointID0 = GetEdgeID(x, y, z, m_triTable[tableIndex][i]);
		  pointID1 = GetEdgeID(x, y, z, m_triTable[tableIndex][i+1]);
		  pointID2 = GetEdgeID(x, y, z, m_triTable[tableIndex][i+2]);
		  triangle.pointID[0] = pointID0;
		  triangle.pointID[1] = pointID1;
		  triangle.pointID[2] = pointID2;
		  m_trivecTriangles.push_back(triangle);
	       }
	    }
	 }
	
#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_1 = std::chrono::high_resolution_clock::now();
#endif
   RenameVerticesAndTriangles();
#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_2 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   std::cout << "   GenerateSurface() d10 " << d10 << "  d21 " << d21 << " microseconds\n";
#endif

   //CalculateNormals();
   m_bValidSurface = true;
}


// This is not used
//
// You might consider this to be inefficient.
//
// Heyho, I suppose strictly speaking it is.  However it takes 14ms to copy over
// a 1Mpt map, so speed doesn't matter, simplicity of writing and reading is
// more important (as is almost always the case of course).
//
// Having said that, it would be fine to have a version which didn't use
// this ptScalarField thing.  But that is substantially more work (I tried
// it first, actually).  Xmap and the grid needs to go directly to (a version of)
// GenerateSurface and instead of m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)
// etc, use map.get_data(c(x,y,z)) [note x,y,z are ints] which is nicer.
//
// Also note that Interpolate() needs to be fixed similarly.
//
//

template <class T> // vector<CartesianPair>
std::pair<int, int>
CIsoSurface<T>::rangeify(const clipper::Grid_map &grid, int isample_step,
			 int isection_start,
			 int isection_end, int n_sections) const {

   // we need to include the last section

   int gmin = grid.min().w();
   int gmax = grid.max().w();

   if (isample_step != 1) {

      // haven't worked this out yet
      return std::pair<int,int>(gmin, gmax);

   } else {

      int grange = gmax - gmin;

      float f1 = static_cast<float>(isection_start)/static_cast<float>(n_sections);
      float f2 = static_cast<float>(isection_end)/static_cast<float>(n_sections);

      // fg2 uses an additional +1 because to get:

      // rangeify input: 0 1 3 gmin 14 gmax 66   output 14 32
      // rangeify input: 1 2 3 gmin 14 gmax 66   output 31 49
      // rangeify input: 2 3 3 gmin 14 gmax 66   output 48 67

      // This covers the gap between sections - we need both edges

      int fg1 = grange * f1 + gmin;
      int fg2 = grange * f2 + gmin + 1;

      if (false)
	 std::cout << ".....rangeify input: " << isection_start << " " << isection_end
		   << " " << n_sections << " gmin " << gmin << " gmax " << gmax
		   << "   output " << fg1 << " " << fg2 << std::endl;

      return std::pair<int, int> (fg1, fg2);
   }
}


// The stardard usage of GenerateSurface_from_Xmap, generated usually from re-centring
// the graphics.
//
template <class T> // vector<CartesianPair>
coot::CartesianPairInfo
CIsoSurface<T>::GenerateSurface_from_Xmap(const clipper::Xmap<T>& crystal_map,
					  T tIsoLevel,
					  float box_radius, // half length
					  coot::Cartesian centre_point,
					  int isample_step,
					  int iream_start, int iream_end, int n_reams,
					  bool is_em_map) {

   // std::cout << "------ start GenerateSurface_from_Xmap() " << n_reams << std::endl;

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_0 = std::chrono::high_resolution_clock::now();
#endif

   // We need to convert the Cartesian centre_point to a grid_coord
   // using coord_orth we can use the xmap cell to give coord_frac 
   // and then use xmap::grid_sampling() to get a coord_grid.
   //
   // Similarly, to generate the extents of the map in which we are 
   // interested, we add and subtract box_radius from centre_point 
   // in all directions to give us u_start, u_end, v_start, v_end and
   // w_start, w_end.
   // 
   // When it comes to writing out the lines/triangles, we will need to add
   // the offset of the bottom left hand corner.
   //

   clipper::Coord_orth centre(centre_point.get_x(),
			      centre_point.get_y(),
			      centre_point.get_z());

   clipper::Coord_frac centref = centre.coord_frac(crystal_map.cell());

   clipper::Coord_frac box0(
			    centref.u() - box_radius/crystal_map.cell().descr().a(),
			    centref.v() - box_radius/crystal_map.cell().descr().b(),
			    centref.w() - box_radius/crystal_map.cell().descr().c() );
   clipper::Coord_frac box1(
			    centref.u() + box_radius/crystal_map.cell().descr().a(),
			    centref.v() + box_radius/crystal_map.cell().descr().b(),
			    centref.w() + box_radius/crystal_map.cell().descr().c() );

   //Note that this introduces a rounding step - is this what you want?
   clipper::Grid_map grid(box0.coord_grid(crystal_map.grid_sampling()),
			  box1.coord_grid(crystal_map.grid_sampling()));

   if (false) { // debug
      std::cout << "    tIsoLevel: " << tIsoLevel << std::endl;
      std::cout << "    box_radius " << box_radius << std::endl;
      std::cout << "    centre_point: " << centre_point << std::endl;
      std::cout << "    isample_step " << isample_step << std::endl;
      std::cout << "    iream_start " << iream_start << std::endl;
      std::cout << "    iream_end " << iream_end << std::endl;
      std::cout << "    n_reams " << n_reams << std::endl;
      std::cout << "    box0: " << box0.format() << std::endl;
      std::cout << "    box1: " << box1.format() << std::endl;
      std::cout << "    grid: " << grid.format() << std::endl;
   }

   std::pair<int, int> rt = rangeify(grid, isample_step, iream_start, iream_end, n_reams);
   clipper::Coord_grid base_grid = grid.min();
   base_grid.w() = rt.first;

   // std::cout << "debug:: allocating ptScalarField grid.size() " << grid.size() << std::endl;
   int grid_size = (rt.second-rt.first+1) * (grid.max().u() - grid.min().u() + 1) * (grid.max().v() - grid.min().v() + 1);
   T* ptScalarField = new T[grid_size];

   clipper::Xmap_base::Map_reference_coord ix(crystal_map);

   int nu = grid.max().u() - grid.min().u() + 1;
   int nv = grid.max().u() - grid.min().u() + 1;
   int icount = 0; // needs offset rt.first * nu * nv (base_grid)
   int w, v, u, ii;

   // note: old test is <=
   // for (w = grid.min().w(); w <= grid.max().w(); w+=isample_step ) {
   // and so rangeify gives a +1 on the second element.

   for (w = rt.first; w <= rt.second; w+=isample_step) {
      for (v = grid.min().v(); v <= grid.max().v(); v+=isample_step) {
	 ix.set_coord(clipper::Coord_grid( grid.min().u(), v, w ));
	 for (u = grid.min().u(); u <= grid.max().u(); u+= isample_step) {
	    // std::cout << "ix " << ix.coord().format() << " icount " << icount << std::endl;
	    if (icount < grid_size) {
	       ptScalarField[icount] = crystal_map[ ix ];
	    } else {
	       std::cout << "ERROR:: out of grid " << icount << " " << grid_size << " "
			 << ix.coord().format() << " min,max "
			 << grid.min().format() << " " << grid.max().format() << std::endl;
	    }
	    icount++;
	    for(ii=0; ii<isample_step; ii++)
	       ix.next_u();
	 }
      }
   }

#ifdef ANALYSE_CONTOURING_TIMING
  auto tp_3 = std::chrono::high_resolution_clock::now();
#endif

  /*
  GenerateSurface(ptScalarField, tIsoLevel,
		  (grid.nu()-1)/isample_step,
		  (grid.nv()-1)/isample_step,
		  (grid.nw()-1)/isample_step,
		  isample_step * 1.0, isample_step * 1.0, isample_step * 1.0);
  */
  GenerateSurface(ptScalarField, tIsoLevel,
		  (grid.nu()-1)/isample_step,
		  (grid.nv()-1)/isample_step,
		  (rt.second-rt.first-1)/isample_step,
		  isample_step * 1.0, isample_step * 1.0, isample_step * 1.0);

#ifdef ANALYSE_CONTOURING_TIMING
  auto tp_4 = std::chrono::high_resolution_clock::now();
#endif
  delete [] ptScalarField; // 7 ms
  
#ifdef ANALYSE_CONTOURING_TIMING
  auto tp_5 = std::chrono::high_resolution_clock::now();
#endif


  clipper::Coord_frac base_frc = base_grid.coord_frac(crystal_map.grid_sampling());
  coot::CartesianPairInfo cpi =
     returnTriangles(crystal_map,
		     base_frc,
		     box_radius, centre_point, is_em_map);
#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_6 = std::chrono::high_resolution_clock::now();

   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
   auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
   auto d65 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_6 - tp_5).count();

   std::cout << " GenerateSurface_from_Xmap"
	     << "  d10 " << d10 << "  d21 " << d21
	     << "  d32 " << d32 << "  d43 " << d43
	     << "  d54 " << d54 << "  d65 " << d65 
	     << " milliseconds\n";
#endif

   return cpi;
}

template <class T> // vector<CartesianPair>
coot::CartesianPairInfo
CIsoSurface<T>::GenerateSurface_from_NXmap(const clipper::NXmap<T>& nx_map,
					   T tIsoLevel,
					   float box_radius, // half length
					   coot::Cartesian centre_point,
					   int isample_step) {

   // When it comes to writing out the lines/triangles, we will need to add
   // the offset of the bottom left hand corner.
   //

   clipper::Coord_orth centre(centre_point.get_x(),
			      centre_point.get_y(),
			      centre_point.get_z());

  
   clipper::Coord_frac centre_fr(0.5, 0.5, 0.5); // some function of centre.

   clipper::Coord_frac box0(centre_fr.u() - 0.2,
			    centre_fr.v() - 0.2,
			    centre_fr.w() - 0.2);
   clipper::Coord_frac box1(centre_fr.u() + 0.2,
			    centre_fr.v() + 0.2,
			    centre_fr.w() + 0.2);

   clipper::Coord_frac grid_min_cf;
   clipper::Coord_grid grid_min;
   clipper::Coord_grid grid_max;
   clipper::Grid_range grid(grid_min, grid_max);

  
   T* ptScalarField = new T[grid.size()];

   std::cout << "box0: " << box0.format() << std::endl
	     << "box1: " << box1.format() << std::endl;

   clipper::NXmap_base::Map_reference_coord ix(nx_map); 
   int icount = 0; 
   for ( int w = grid.min().w(); w <= grid.max().w(); w+=isample_step ) { 
      for ( int v = grid.min().v(); v <= grid.max().v(); v+=isample_step ) { 
	 ix.set_coord( clipper::Coord_grid( grid.min().u(), v, w ) ); 
	 for ( int u = grid.min().u(); u <= grid.max().u(); u+= isample_step ) { 
	    ptScalarField[icount] = nx_map[ ix ]; 
	    icount++;
	    for(int ii=0; ii<isample_step; ii++) 
	       ix.next_u(); 
	 } 
      } 
   } 
   
   GenerateSurface(ptScalarField, tIsoLevel,
		   (grid.nu()-1)/isample_step,
		   (grid.nv()-1)/isample_step,
		   (grid.nw()-1)/isample_step,
		   isample_step * 1.0, isample_step * 1.0, isample_step * 1.0);
  
   delete [] ptScalarField;
   return returnTriangles(nx_map,
			  grid_min_cf,
			  box_radius,
			  centre_point);
} 




// ======================================================================
//            Triangles - for solid surface rendering
// ======================================================================
//
template <class T>
coot::density_contour_triangles_container_t
CIsoSurface<T>::GenerateTriangles_from_Xmap(const clipper::Xmap<T>& crystal_map,
					    T tIsoLevel,
					    float box_radius, // half length
					    coot::Cartesian centre_point,
					    int isample_step) {

   coot::density_contour_triangles_container_t tri_con;

   // We need to convert the Cartesian centre_point to a grid_coord
   // using coord_orth we can use the xmap cell to give coord_frac 
   // and then use xmap::grid_sampling() to get a coord_grid.
   //
   // Similarly, to generate the extents of the map in which we are 
   // interested, we add and subtract box_radius from centre_point 
   // in all directions to give us u_start, u_end, v_start, v_end and
   // w_start, w_end.
   // 
   // When it comes to writing out the lines/triangles, we will need to add
   // the offset of the bottom left hand corner.
   //

   clipper::Coord_orth centre( centre_point.get_x(), centre_point.get_y(),
			       centre_point.get_z() );
   double radius_sqd = box_radius * box_radius; // convert to double

   // clipper::Coord_frac centref = crystal_map.cell().to_frac( centre );
  
   clipper::Coord_frac centref = centre.coord_frac(crystal_map.cell() ); 

   clipper::Coord_frac box0(
			    centref.u() - box_radius/crystal_map.cell().descr().a(),
			    centref.v() - box_radius/crystal_map.cell().descr().b(),
			    centref.w() - box_radius/crystal_map.cell().descr().c() );
   clipper::Coord_frac box1(
			    centref.u() + box_radius/crystal_map.cell().descr().a(),
			    centref.v() + box_radius/crystal_map.cell().descr().b(),
			    centref.w() + box_radius/crystal_map.cell().descr().c() );

   // old style (early 2002) convertion operator:
   // 
   //clipper::Grid_map grid( crystal_map.grid_sampling().to_grid( box0 ),
   // crystal_map.grid_sampling().to_grid( box1 ) );

   // using this constructor (coords.h):
   //! constructor: takes grid limits
   // Grid_map( const Coord_grid& min, const Coord_grid& max );

   // question: we have a Coord_frac and want to convert it to a 
   // Coord_grid.  How?

   // note: 
   // Coord_frac::coord_grid(const Grid& g) returns a Coord_grid
   // but how do we get a Grid from an Xmap?
  
   //Note that this introduces a rounding step - is this what you want?
   clipper::Grid_map grid( box0.coord_grid(crystal_map.grid_sampling()),
			   box1.coord_grid(crystal_map.grid_sampling()));
			  
   //   cout << "INFO: centre_point is :" << centre_point << endl;
   //   cout << "INFO: box0         is :" << box0.format() << endl;
   //   cout << "INFO: box1         is :" << box1.format() << endl;



   T* ptScalarField = new T[grid.size()];

   //cout << "box0: " << box0.format() << endl
   //    << "box1: " << box1.format() << endl;

   clipper::Xmap_base::Map_reference_coord ix( crystal_map ); 
   int icount = 0; 
   for ( int w = grid.min().w(); w <= grid.max().w(); w+=isample_step ) { 
      for ( int v = grid.min().v(); v <= grid.max().v(); v+=isample_step ) { 
	 ix.set_coord( clipper::Coord_grid( grid.min().u(), v, w ) ); 
	 for ( int u = grid.min().u(); u <= grid.max().u(); u+= isample_step ) { 
	    ptScalarField[icount] = crystal_map[ ix ]; 
	    icount++;
	    for(int ii=0; ii<isample_step; ii++) 
	       ix.next_u(); 
	 } 
      } 
   } 
   
   GenerateSurface(ptScalarField, tIsoLevel,
		   (grid.nu()-1)/isample_step,
		   (grid.nv()-1)/isample_step,
		   (grid.nw()-1)/isample_step,
		   isample_step * 1.0, isample_step * 1.0, isample_step * 1.0);

   delete [] ptScalarField;

   // now fill tri_con
   clipper::Coord_frac base = grid.min().coord_frac(crystal_map.grid_sampling());
   T nu = crystal_map.grid_sampling().nu();
   T nv = crystal_map.grid_sampling().nv();
   T nw = crystal_map.grid_sampling().nw();

   // what is the maximum index in m_piTriangleIndices ?  (we
   // shouldn't need to do this - it should be clear(?) from other
   // code what this number is, c.f. check_max_min_vertices()
   //
   unsigned int max_index = 0;
   for (unsigned int i=0; i < m_nTriangles*3; i++) {
      if (m_piTriangleIndices[i] > max_index)
	 max_index = m_piTriangleIndices[i];
   }
   tri_con.points.resize(max_index+1);
   tri_con.normals.resize(max_index+1);
   
   //
   unsigned nt_for_index = 0;
   for (unsigned int nt=0; nt < m_nTriangles; nt++) {

      int i = nt*3;
      unsigned int j   = m_piTriangleIndices[i]; 
      unsigned int jp  = m_piTriangleIndices[i+1]; 
      unsigned int jp2 = m_piTriangleIndices[i+2];

      clipper::Coord_frac cf_1 = clipper::Coord_frac(m_ppt3dVertices[j][0]/nu,
						     m_ppt3dVertices[j][1]/nv,
						     m_ppt3dVertices[j][2]/nw) + base;
      clipper::Coord_frac cf_2 = clipper::Coord_frac(m_ppt3dVertices[jp][0]/nu,
						     m_ppt3dVertices[jp][1]/nv,
						     m_ppt3dVertices[jp][2]/nw) + base;
      clipper::Coord_frac cf_3 = clipper::Coord_frac(m_ppt3dVertices[jp2][0]/nu,
						     m_ppt3dVertices[jp2][1]/nv,
						     m_ppt3dVertices[jp2][2]/nw) + base;

      clipper::Coord_orth co_1 = cf_1.coord_orth(crystal_map.cell());
      clipper::Coord_orth co_2 = cf_2.coord_orth(crystal_map.cell());
      clipper::Coord_orth co_3 = cf_3.coord_orth(crystal_map.cell());
      tri_con.points[j  ] = co_1;
      tri_con.points[jp ] = co_2;
      tri_con.points[jp2] = co_3;

      clipper::Coord_orth sum_pt = co_1;
      sum_pt += co_2;
      sum_pt += co_3;

      TRIANGLE tri;
      tri.pointID[0] = j;
      tri.pointID[1] = jp;
      tri.pointID[2] = jp2;
      tri.mid_point = clipper::Coord_orth(0.333333333 * sum_pt.x(),
                                          0.333333333 * sum_pt.y(),
                                          0.333333333 * sum_pt.z());
      tri.back_front_projection_distance = 0; // Will be reset

      bool valid_co = true;

      // Don't add this triangle if it's outside the sphere
      //
      if ((tri.mid_point-centre).lengthsq() > radius_sqd) valid_co = false;
      // If you want to bring back "all" triangle, test on true
      if (valid_co) {

         // Note we apply a negation to get the normal pointing out of
         // the surface, so the shiny surface is on the outside.
         //
         tri.normal_for_flat_shading = clipper::Coord_orth(-clipper::Coord_orth::cross((co_2-co_1), (co_3-co_1)).unit());

         // If the contour level is negative then the normals need to
         // point in the other direction (c.f. a positive contour).  If
         // we don't do this, the bright shiny surfaces of the negative
         // level are on the inside.
         //
         if (tIsoLevel < 0.0)
	    tri.normal_for_flat_shading = - tri.normal_for_flat_shading;
      
         tri_con.point_indices.push_back(tri);
      }
   }

   tri_con.calculate_normals();
   return tri_con;
}



template <class T> bool CIsoSurface<T>::IsSurfaceValid()
{
	return m_bValidSurface;
}

template <class T> void CIsoSurface<T>::DeleteSurface()
{
	m_fCellLengthX = 0;
	m_fCellLengthY = 0;
	m_fCellLengthZ = 0;
	m_nCellsX = 0;
	m_nCellsY = 0;
	m_nCellsZ = 0;
	m_nTriangles = 0;
	m_nNormals = 0;
	m_nVertices = 0;
	if (m_ppt3dVertices != NULL) {
		delete[] m_ppt3dVertices;
		m_ppt3dVertices = NULL;
	}
	if (m_piTriangleIndices != NULL) {
		delete[] m_piTriangleIndices;
		m_piTriangleIndices = NULL;
	}
	if (m_pvec3dNormals != NULL) {
		delete[] m_pvec3dNormals;
		m_pvec3dNormals = NULL;
	}
	m_ptScalarField = NULL;
	m_tIsoLevel = 0;
	m_bValidSurface = false;
}

template <class T> int CIsoSurface<T>::GetVolumeLengths(float& fVolLengthX, float& fVolLengthY, float& fVolLengthZ)
{
	if (IsSurfaceValid()) {
		fVolLengthX = m_fCellLengthX*m_nCellsX;
		fVolLengthY = m_fCellLengthY*m_nCellsY;
		fVolLengthZ = m_fCellLengthZ*m_nCellsZ;
		return 1;
	}
	else
		return -1;
}

template <class T> unsigned int CIsoSurface<T>::GetEdgeID(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo)
{
	switch (nEdgeNo) {
	case 0:
		return GetVertexID(nX, nY, nZ) + 1;
	case 1:
		return GetVertexID(nX, nY + 1, nZ);
	case 2:
		return GetVertexID(nX + 1, nY, nZ) + 1;
	case 3:
		return GetVertexID(nX, nY, nZ);
	case 4:
		return GetVertexID(nX, nY, nZ + 1) + 1;
	case 5:
		return GetVertexID(nX, nY + 1, nZ + 1);
	case 6:
		return GetVertexID(nX + 1, nY, nZ + 1) + 1;
	case 7:
		return GetVertexID(nX, nY, nZ + 1);
	case 8:
		return GetVertexID(nX, nY, nZ) + 2;
	case 9:
		return GetVertexID(nX, nY + 1, nZ) + 2;
	case 10:
		return GetVertexID(nX + 1, nY + 1, nZ) + 2;
	case 11:
		return GetVertexID(nX + 1, nY, nZ) + 2;
	default:
		// Invalid edge no.
		return -1;
	}
}

template <class T> unsigned int CIsoSurface<T>::GetVertexID(unsigned int nX, unsigned int nY, unsigned int nZ)
{
	return 3*(nZ*(m_nCellsY + 1)*(m_nCellsX + 1) + nY*(m_nCellsX + 1) + nX);
}

template <class T> POINT3DID CIsoSurface<T>::CalculateIntersection(unsigned int nX, unsigned int nY, unsigned int nZ, unsigned int nEdgeNo)
{
	float x1, y1, z1, x2, y2, z2;
	unsigned int v1x = nX, v1y = nY, v1z = nZ;
	unsigned int v2x = nX, v2y = nY, v2z = nZ;
	
	switch (nEdgeNo)
	{
	case 0:
		v2y += 1;
		break;
	case 1:
		v1y += 1;
		v2x += 1;
		v2y += 1;
		break;
	case 2:
		v1x += 1;
		v1y += 1;
		v2x += 1;
		break;
	case 3:
		v1x += 1;
		break;
	case 4:
		v1z += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 5:
		v1y += 1;
		v1z += 1;
		v2x += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 6:
		v1x += 1;
		v1y += 1;
		v1z += 1;
		v2x += 1;
		v2z += 1;
		break;
	case 7:
		v1x += 1;
		v1z += 1;
		v2z += 1;
		break;
	case 8:
		v2z += 1;
		break;
	case 9:
		v1y += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 10:
		v1x += 1;
		v1y += 1;
		v2x += 1;
		v2y += 1;
		v2z += 1;
		break;
	case 11:
		v1x += 1;
		v2x += 1;
		v2z += 1;
		break;
	}

	x1 = v1x*m_fCellLengthX;
	y1 = v1y*m_fCellLengthY;
	z1 = v1z*m_fCellLengthZ;
	x2 = v2x*m_fCellLengthX;
	y2 = v2y*m_fCellLengthY;
	z2 = v2z*m_fCellLengthZ;

	unsigned int nPointsInXDirection = (m_nCellsX + 1);
	unsigned int nPointsInSlice = nPointsInXDirection*(m_nCellsY + 1);
	T val1 = m_ptScalarField[v1z*nPointsInSlice + v1y*nPointsInXDirection + v1x];
	T val2 = m_ptScalarField[v2z*nPointsInSlice + v2y*nPointsInXDirection + v2x];
	POINT3DID intersection = Interpolate(x1, y1, z1, x2, y2, z2, val1, val2);
	
	return intersection;
}

template <class T> POINT3DID CIsoSurface<T>::Interpolate(float fX1, float fY1, float fZ1, float fX2, float fY2, float fZ2, T tVal1, T tVal2)
{
	POINT3DID interpolation;
	interpolation.newID = 0;
	float mu;

	mu = float((m_tIsoLevel - tVal1))/(tVal2 - tVal1);
	interpolation.x = fX1 + mu*(fX2 - fX1);
	interpolation.y = fY1 + mu*(fY2 - fY1);
	interpolation.z = fZ1 + mu*(fZ2 - fZ1);

	return interpolation;
}

#include "utils/coot-utils.hh" // for get_max_number_of_threads()

template <class T> void CIsoSurface<T>::RenameVerticesAndTriangles() {

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_0 = std::chrono::high_resolution_clock::now();
#endif

   unsigned int nextID = 0;
   ID2POINT3DID::iterator mapIterator = m_i2pt3idVertices.begin();
   TRIANGLEVECTOR::iterator vecIterator = m_trivecTriangles.begin();

   // Rename vertices.
   while (mapIterator != m_i2pt3idVertices.end()) {
      mapIterator->second.newID = nextID;
      nextID++;
      mapIterator++;
   }

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_1 = std::chrono::high_resolution_clock::now();
#endif

   // Now rename triangles (don't do this with (now inner) threads)
   while (vecIterator != m_trivecTriangles.end()) {
      for (unsigned int i=0; i<3; i++) {
	 unsigned int newID = m_i2pt3idVertices.at(vecIterator->pointID[i]).newID;
	 vecIterator->pointID[i] = newID;
      }
      vecIterator++;
   }
   
#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_2 = std::chrono::high_resolution_clock::now();
#endif

   // Copy all the vertices and triangles into two arrays so that they
   // can be efficiently accessed.
   // Copy vertices.
   mapIterator = m_i2pt3idVertices.begin();
   m_nVertices = m_i2pt3idVertices.size();
   m_ppt3dVertices = new POINT3D[m_nVertices];

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_3 = std::chrono::high_resolution_clock::now();
#endif

   for (unsigned int i = 0; i < m_nVertices; i++, mapIterator++) {
      m_ppt3dVertices[i][0] = (*mapIterator).second.x;
      m_ppt3dVertices[i][1] = (*mapIterator).second.y;
      m_ppt3dVertices[i][2] = (*mapIterator).second.z;
   }

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_4 = std::chrono::high_resolution_clock::now();
#endif

   // Copy vertex indices which make triangles.
   vecIterator = m_trivecTriangles.begin();
   m_nTriangles = m_trivecTriangles.size();
   m_piTriangleIndices = new unsigned int[m_nTriangles*3];
   for (unsigned int i = 0; i < m_nTriangles; i++, vecIterator++) {
      m_piTriangleIndices[i*3  ] = (*vecIterator).pointID[0];
      m_piTriangleIndices[i*3+1] = (*vecIterator).pointID[1];
      m_piTriangleIndices[i*3+2] = (*vecIterator).pointID[2];
   }

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_5 = std::chrono::high_resolution_clock::now();
#endif

   m_i2pt3idVertices.clear();
   m_trivecTriangles.clear();

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_6 = std::chrono::high_resolution_clock::now();

   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
   auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
   auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
   auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
   auto d65 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_6 - tp_5).count();

   std::cout << "   RenameVerticesAndTriangles d10 " << d10 << "  d21 " << d21
	     << "  d32 " << d32 << "  d43 " << d43
	     << "  d54 " << d54 << "  d65 " << d65
	     << " milliseconds\n";
#endif
}

// static
template <class T> void CIsoSurface<T>::rename_tris_in_thread(const std::pair<unsigned int, unsigned int> &idx_range, TRIANGLEVECTOR &tv, const ID2POINT3DID &point_map) {

   for (std::size_t idx=idx_range.first; idx<idx_range.second; idx++) {
      for (unsigned int i=0; i<3; i++) {

	 unsigned int new_id = point_map.at(tv[idx].pointID[i]).newID;
	 tv[idx].pointID[i] = new_id;
      }
   }
}


// debugging function
template <class T> void CIsoSurface<T>::check_max_min_vertex_index_from_triangles() {

   unsigned int max_x = -999, max_y = -999, max_z = -999; 
   unsigned int min_x =  999, min_y =  999, min_z =  999;

   unsigned int v_max = 0;
   unsigned int v1; 

   std::cout << "checking m_nTriangles=" << m_nTriangles << " triangles\n"; 
   std::cout << "         m_nVertices =" << m_nVertices << " vertices\n"; 

   for (unsigned int i = 0; i < m_nTriangles; i++) {

      // recall that m_piTriangleIndices is a list of unsigned ints
      // (indexing the vertices).

      v1 = m_piTriangleIndices[i];

      if (v1 > v_max) {
	 v_max = v1;
      }
      
   }

   std::cout << "max vertex from triangle usage is: " << v_max << std::endl;

}


// debugging function
template <class T> void CIsoSurface<T>::check_max_min_vertices() {

   T max_x = -999, max_y = -999, max_z = -999; 
   T min_x =  999, min_y =  999, min_z =  999;

   std::cout << "checking m_nVertices=" << m_nVertices << " vertices\n"; 
   for (unsigned int i = 0; i < m_nVertices; i++) {

      if (m_ppt3dVertices[i][0] > max_x)
	 max_x = m_ppt3dVertices[i][0]; 
      if (m_ppt3dVertices[i][1] > max_y)
	 max_y = m_ppt3dVertices[i][1]; 
      if (m_ppt3dVertices[i][2] > max_z)
	 max_z = m_ppt3dVertices[i][2]; 

      if (m_ppt3dVertices[i][0] < min_x)
	 min_x = m_ppt3dVertices[i][0]; 
      if (m_ppt3dVertices[i][1] < min_y)
	 min_y = m_ppt3dVertices[i][1]; 
      if (m_ppt3dVertices[i][2] < min_z)
	 min_z = m_ppt3dVertices[i][2]; 

   }

   std::cout << "Debug: check_max_min_vertices (min and max x, y and z): \n"
	<< min_x << " " << max_x << "\n"
	<< min_y << " " << max_y << "\n"
	<< min_z << " " << max_z << "\n";
   
}
   

template <class T> void CIsoSurface<T>::CalculateNormals()
{
	m_nNormals = m_nVertices;
	m_pvec3dNormals = new VECTOR3D[m_nNormals];
	
	// Set all normals to 0.
	for (unsigned int i = 0; i < m_nNormals; i++) {
		m_pvec3dNormals[i][0] = 0;
		m_pvec3dNormals[i][1] = 0;
		m_pvec3dNormals[i][2] = 0;
	}

	// Calculate normals.
	for (unsigned int i = 0; i < m_nTriangles; i++) {
		VECTOR3D vec1, vec2, normal;
		unsigned int id0, id1, id2;
		id0 = m_piTriangleIndices[i*3];
		id1 = m_piTriangleIndices[i*3+1];
		id2 = m_piTriangleIndices[i*3+2];
		vec1[0] = m_ppt3dVertices[id1][0] - m_ppt3dVertices[id0][0];
		vec1[1] = m_ppt3dVertices[id1][1] - m_ppt3dVertices[id0][1];
		vec1[2] = m_ppt3dVertices[id1][2] - m_ppt3dVertices[id0][2];
		vec2[0] = m_ppt3dVertices[id2][0] - m_ppt3dVertices[id0][0];
		vec2[1] = m_ppt3dVertices[id2][1] - m_ppt3dVertices[id0][1];
		vec2[2] = m_ppt3dVertices[id2][2] - m_ppt3dVertices[id0][2];
		normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
		normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
		normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
		m_pvec3dNormals[id0][0] += normal[0];
		m_pvec3dNormals[id0][1] += normal[1];
		m_pvec3dNormals[id0][2] += normal[2];
		m_pvec3dNormals[id1][0] += normal[0];
		m_pvec3dNormals[id1][1] += normal[1];
		m_pvec3dNormals[id1][2] += normal[2];
		m_pvec3dNormals[id2][0] += normal[0];
		m_pvec3dNormals[id2][1] += normal[1];
		m_pvec3dNormals[id2][2] += normal[2];
	}

	// Normalize normals.
	for (unsigned int i = 0; i < m_nNormals; i++) {
		float length = sqrt(m_pvec3dNormals[i][0]*m_pvec3dNormals[i][0] + m_pvec3dNormals[i][1]*m_pvec3dNormals[i][1] + m_pvec3dNormals[i][2]*m_pvec3dNormals[i][2]);
		m_pvec3dNormals[i][0] /= length;
		m_pvec3dNormals[i][1] /= length;
		m_pvec3dNormals[i][2] /= length;
	}
}

#define TRUE  1
#define FALSE 0

// PE
template <class T> unsigned int CIsoSurface<T>::nTriangles(void) {

   return m_nTriangles;
}

// PE
template <class T> void CIsoSurface<T>::morphVertices(void) {

   // We look first for cluster of points that are close to each other.
   //
   // Where are the points?
   //
   // Well, looking at CalculateNormals(), we see that the list
   // m_piTriangleIndices() is used and that is used by pulling off 3
   // indexes at a time, these point correspond to a triangle.  It seems
   // that this does not use the TRIANGLE structure. 
   //
   // We can pop off m_nTriangles triangles.
   //
   // For TriangleIndex (actually point index) i, the vertex is at:
   // (m_ppt3dVertices[i][0], m_ppt3dVertices[i][1], m_ppt3dVertices[i][2])
   // 

   // small triangles vector with stl usage.
   //
   std::vector<int> small_triangles; 
   
   for (unsigned int i=0; i< m_nTriangles*3; i+=3) {
      if (isSmallTriangle(i)) {
	 //
	 adjustVertices(i);
	 small_triangles.push_back(i);

	 // we need to find who uses points i, i+1, i+2 and
	 // change their references to this new points.
	 //
	 // Also, we should remove this triangle from the list.
	 // Sounds like a job for std vector.
	 //
	 // The question is, should we delete this small triangle
	 // from the Vertices array.  We can decide that later,
	 // lets make a vector of them
      }
   }

   // now lets use that small_triangles vector
   //
   std::cout << "We found " << small_triangles.size() << " small triangles\n";
   //for (unsigned int i=0; i< small_triangles.size(); i++) {

      // cout << small_triangles[i] << " was the position of a small triangle\n";
   //}

   // copy over the established triangles into the new triangle list
   // stepping over the triangles that do not exist.
   // 
   // We shall keep a (+ve) offset (so that it can be declared unsigned)
   // 

   // I really want to do this with recursion (this construction is a 2-liner
   // in scheme).
   //
   // I should take a moment to work out how to do recursion in c++.
   // 
   unsigned int offset=0;
   int n_deleted = 0;
   int debug_counter = 0;
   
   // unsigned int *newTriangles = new unsigned int [m_nTriangles*3 ];
   // 3*small_triangles.size()];
   

//    for (unsigned int i=0; i< m_nTriangles*3; i+=3) {

//       // cout << "comparing " << i << " and "
//       //      << small_triangles[offset] << endl;
      
//       if (i == small_triangles[offset]) {
// 	 // we skip the copy of this triangle
// 	 //
// 	 // and update to the next small_triangle
// 	 offset++;
// 	 n_deleted++; 
//       } else {
// 	 // OK to copy (the usual case, of course).
// 	 newTriangles[i-3*offset  ] = m_piTriangleIndices[i  ];
// 	 newTriangles[i-3*offset+1] = m_piTriangleIndices[i+1];
// 	 newTriangles[i-3*offset+2] = m_piTriangleIndices[i+2];

// 	 if (debug_counter < i-3*offset+2) {
// 	    debug_counter = i-3*offset+2;
// 	    // cout << "updating debug_counter to " << debug_counter << endl;
// 	 }
//       }
//    }

   // cout << "we found " << n_deleted << " triangles\n";
   //cout << "maximum m_piTriangleIndices can be indexed to is "
   // << debug_counter << endl;

   // now adjust the number of triangles, now that we have deleted
   // some of them.
   //
   //cout << "Reducing m_nTriangles (" << m_nTriangles
   //	<< ") by 3*" << small_triangles.size() << endl;
   
   // m_nTriangles -= 3*small_triangles.size();
   // and finally delete the old m_piTriangleIndices and
   // reassign it to the new array.

   //delete [] m_piTriangleIndices;
   //m_piTriangleIndices = newTriangles; 
   
}

// This is the easy/cheap way of doing this because we are not
// creating/deleting new vertices (and hence not re-ording the 
// indicies of the triangles that use them).  
// 
template <class T> void CIsoSurface<T>::adjustVertices(unsigned int i) {

   T t1_x, t1_y, t1_z;
   T t2_x, t2_y, t2_z;
   T t3_x, t3_y, t3_z;

   unsigned int j0, j1, j2; 

   j0 = m_piTriangleIndices[i]; 
   j1 = m_piTriangleIndices[i+1]; 
   j2 = m_piTriangleIndices[i+2]; 
   
   t1_x = m_ppt3dVertices[j0][0]; 
   t1_y = m_ppt3dVertices[j0][1]; 
   t1_z = m_ppt3dVertices[j0][2]; 

   t2_x = m_ppt3dVertices[j1][0]; 
   t2_y = m_ppt3dVertices[j1][1]; 
   t2_z = m_ppt3dVertices[j1][2]; 

   t3_x = m_ppt3dVertices[j2][0]; 
   t3_y = m_ppt3dVertices[j2][1]; 
   t3_z = m_ppt3dVertices[j2][2];

   // POINT3DID point;

   float x = (t1_x + t2_x + t3_x)/3.0; 
   float y = (t1_y + t2_y + t3_y)/3.0; 
   float z = (t1_z + t2_z + t3_z)/3.0;

   // first point
    m_ppt3dVertices[j0][0] = x;
    m_ppt3dVertices[j0][1] = y;
    m_ppt3dVertices[j0][2] = z;

    // let the other points be this point too.
    // 
    m_piTriangleIndices[i+1] = j0;
    m_piTriangleIndices[i+2] = j0;


}
      

   

//
template <class T> bool CIsoSurface<T>::isSmallTriangle(unsigned int i) {

   T t1_x, t1_y, t1_z;
   T t2_x, t2_y, t2_z;
   T t3_x, t3_y, t3_z;

   unsigned int j; 

   j = m_piTriangleIndices[i]; 
   
   t1_x = m_ppt3dVertices[j][0]; 
   t1_y = m_ppt3dVertices[j][1]; 
   t1_z = m_ppt3dVertices[j][2]; 

   j = m_piTriangleIndices[i+1]; 

   t2_x = m_ppt3dVertices[j][0]; 
   t2_y = m_ppt3dVertices[j][1]; 
   t2_z = m_ppt3dVertices[j][2]; 

   j = m_piTriangleIndices[i+2]; 

   t3_x = m_ppt3dVertices[j][0]; 
   t3_y = m_ppt3dVertices[j][1]; 
   t3_z = m_ppt3dVertices[j][2];

   T small_dist = 0.1;

   if ( fabsf((float) (t1_x-t2_x)) < small_dist) { 
      if ( fabsf((float) (t1_y-t2_y)) < small_dist) { 
	 if ( fabsf((float) (t1_z-t2_z)) < small_dist) {
	    if ( fabsf((float) (t1_x-t3_x)) < small_dist) { 
	       if ( fabsf((float) (t1_y-t3_y)) < small_dist) { 
		  if ( fabsf((float) (t1_z-t3_z)) < small_dist) {

		     //cout << "we found a small triangle:\n";
		     //cout << "(" << t1_x << ", " << t1_y << ", " << t1_z << ")\n";
		     //cout << "(" << t2_x << ", " << t2_y << ", " << t2_z << ")\n";
		     //cout << "(" << t3_x << ", " << t3_y << ", " << t3_z << ")\n";

		     return TRUE; 
		  }
	       }
	    }
	 }
      }
   }

   return FALSE; 

}

// This function is not used.
// 
template <class T> void CIsoSurface<T>::writeTriangles(std::string filename) {

   T t1_x, t1_y, t1_z;
   T t2_x, t2_y, t2_z;
   T t3_x, t3_y, t3_z;

   unsigned int j;

   // question: what is the maximum m_piTriangleIndices can be indexed to? 

   std::cout << "In writeTriangles, m_nVertices is " << m_nVertices
	<< " and m_nTriangles is " << m_nTriangles << std::endl;

   check_max_min_vertices(); // debugging
   
   int i_tri_out = 0;

   std::ofstream outfile(filename.c_str());; // c_str() is needed.

   if (! outfile) {
      std::cout << "Could not open " << filename.c_str() << " for some reason\n";
   }

   done_line_list_t done_line_list;
   to_vertex_list_t to_vertex_list; 

   for (unsigned int i=0; i < m_nTriangles*3; i+=3) {
   
      j = m_piTriangleIndices[i]; 
      
      t1_x = m_ppt3dVertices[j][0]; 
      t1_y = m_ppt3dVertices[j][1]; 
      t1_z = m_ppt3dVertices[j][2]; 
      
      j = m_piTriangleIndices[i+1];

      t2_x = m_ppt3dVertices[j][0]; 
      t2_y = m_ppt3dVertices[j][1]; 
      t2_z = m_ppt3dVertices[j][2]; 
      
      j = m_piTriangleIndices[i+2]; 
      
      t3_x = m_ppt3dVertices[j][0]; 
      t3_y = m_ppt3dVertices[j][1]; 
      t3_z = m_ppt3dVertices[j][2];

      // to_vertex_list = done_line_list.to_vertices[m_piTriangleIndices[i]];
      //
      //to_vertex_list = done_line_list.to_vertices.at(m_piTriangleIndices[i]);
      
      // if (to_vertex_list.vertex_list.at(m_piTriangleIndices[i+1]) == 1) {

      //cout << "done bond: " << m_piTriangleIndices[i]
      //      << " to " << m_piTriangleIndices[i+1] << endl;

	 //} else { 
	 outfile << i << "\n";

	 outfile.setf(std::ios::scientific); //, ios::floatfield);
	 
	 outfile << t1_x << " " << t1_y << " " << t1_z << "\n";
	 outfile << t2_x << " " << t2_y << " " << t2_z << "\n";
	 outfile << t3_x << " " << t3_y << " " << t3_z << "\n";
	 
	 i_tri_out++;

	 // now mark it as done
	 //ndone_line_list[m_piTriangleIndices[i]].assign(m_piTriangleIndices[i+1],1);
	 
	 //      }
   }

   outfile.close();

   std::cout << "we wrote " << i_tri_out << " triangles. 3*" << i_tri_out
	<< "=" << 3*i_tri_out << " to " << filename.c_str() << std::endl;

}

bool
do_line(done_line_list_t &done_line_list, int j, int jp) {

   if (done_line_list.done_before(j,jp) == 1) {
      //cout << "top: " << j << "," << jp << " done before!" << endl;
      return 0;
   } else {
      // _was_ new (now marked as done). 
      // cout << "top: " << j << "," << jp << " was new" << endl;
      return 1;
   }
}


template <class T>
coot::CartesianPairInfo
CIsoSurface<T>::returnTriangles(const clipper::Xmap<T>& xmap,
				const clipper::Coord_frac& base,
				float radius,
				coot::Cartesian centre,
				bool is_em_map) const {

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_0 = std::chrono::high_resolution_clock::now();
#endif
   coot::CartesianPairInfo result_wrapper;

   // result_wrapper.data = result;
   // result_wrapper.size = line_index; 

   T nu = xmap.grid_sampling().nu();
   T nv = xmap.grid_sampling().nv();
   T nw = xmap.grid_sampling().nw();

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_1 = std::chrono::high_resolution_clock::now();
#endif
   result_wrapper.data = new coot::CartesianPair[m_nTriangles*3]; // at most it can be this
   result_wrapper.size = 0; // indexer of the result array.
#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_2 = std::chrono::high_resolution_clock::now();
#endif

   clipper::Coord_frac cf;
   clipper::Coord_orth co1, co2, co3;
   float radius_sqd = radius * radius;
   clipper::Coord_orth centre_clipper(centre.x(), centre.y(), centre.z());

   done_line_list_t done_line_list;

   int face_count_1 = 0, face_count_2 = 0, face_count_3 = 0;
   int not_passed_back_count = 0;
   short int face;

   unsigned int done_count = 0, d1_2, d2_3, d1_3;

   coot::Cartesian co1_c;
   coot::Cartesian co2_c;
   coot::Cartesian co3_c;

   bool valid_co_1 = true;
   bool valid_co_2 = true;
   bool valid_co_3 = true;

   double max_x=0, max_y=0, max_z=0;
   if (is_em_map) {
      max_x = xmap.cell().descr().a();
      max_y = xmap.cell().descr().b();
      max_z = xmap.cell().descr().c();
   }

#ifdef ANALYSE_CONTOURING_TIMING
   auto tp_3 = std::chrono::high_resolution_clock::now();
#endif

   unsigned int m_nTriangles_x_3 = m_nTriangles*3; // do the multiply once

   for (unsigned int i=0; i < m_nTriangles_x_3; i+=3) {

      unsigned int j   = m_piTriangleIndices[i]; 
      unsigned int jp  = m_piTriangleIndices[i+1]; 
      unsigned int jp2 = m_piTriangleIndices[i+2];

      cf = clipper::Coord_frac(m_ppt3dVertices[j][0]/nu,
			       m_ppt3dVertices[j][1]/nv,
			       m_ppt3dVertices[j][2]/nw) + base;
      co1 = cf.coord_orth(xmap.cell());
      co1_c = coot::Cartesian(co1.x(), co1.y(), co1.z());

      cf = clipper::Coord_frac(m_ppt3dVertices[jp][0]/nu,
			       m_ppt3dVertices[jp][1]/nv,
			       m_ppt3dVertices[jp][2]/nw) + base;
      co2 = cf.coord_orth(xmap.cell());
      co2_c = coot::Cartesian(co2.x(), co2.y(), co2.z());

      cf = clipper::Coord_frac(m_ppt3dVertices[jp2][0]/nu,
			       m_ppt3dVertices[jp2][1]/nv,
			       m_ppt3dVertices[jp2][2]/nw) + base;
      co3 = cf.coord_orth(xmap.cell());
      co3_c =  coot::Cartesian( co3.x(), co3.y(), co3.z());

      valid_co_1 = true;
      valid_co_2 = true;
      valid_co_3 = true;

      if ((co1_c-centre).amplitude_squared() > radius_sqd)
	 valid_co_1 = false;
      if ((co2_c-centre).amplitude_squared() > radius_sqd)
	 valid_co_2 = false;
      if ((co3_c-centre).amplitude_squared() > radius_sqd)
	 valid_co_3 = false;

      if (valid_co_1 && valid_co_2)
	 result_wrapper.data[result_wrapper.size++] = coot::CartesianPair(co1_c, co2_c);
      if (valid_co_1 && valid_co_3)
	 result_wrapper.data[result_wrapper.size++] = coot::CartesianPair(co1_c, co3_c); 
      if (valid_co_2 && valid_co_3)
	 result_wrapper.data[result_wrapper.size++] = coot::CartesianPair(co3_c, co2_c);


   }

   return result_wrapper;
}

template <class T> // vector<CartesianPair>
coot::CartesianPairInfo
CIsoSurface<T>::returnTriangles( const clipper::NXmap<T>& nx_map,
				 const clipper::Coord_frac& base,
				 float radius,
				 coot::Cartesian centre) const {

   coot::CartesianPairInfo result_wrapper;
   
   coot::CartesianPair *result = new coot::CartesianPair[m_nTriangles*3]; // at most it can be this
   int line_index = 0; // indexer of the result array.
   
   result_wrapper.data = result;
   result_wrapper.size = line_index; 
      
   return result_wrapper;
}


// -----------------------------------------------------------------
// testing stuff
// 

// i is the from index, j is the to index.
//
// We'll check once and mark both ways
bool
done_line_list_t::done_before(int i, int j) {

   //cout << "i=" << i << " j=" << j
   //	<< " max_from_vertex=" << max_from_vertex << endl;

   
   int itmp = i > j ? i : j;

   // Because, if this had been done before, the array would have
   // been of the right size.
   // 
   if (itmp >= from_vertices_size) { 
      resize_and_copy(itmp);
      mark_as_done(i,j);
      return 0;
   }

   // Same as above, but this time we do not need to extend to array.
   //
   if (itmp > max_from_vertex) {
      mark_as_done(i,j);
      return 0;
   }

   if  ( (from_vertices[i]).contains(j) == 1) {
      return 1;
   } else {
      mark_as_done(i,j);
      return 0;
   } 

}

done_line_list_t::done_line_list_t() {

   int start_size = 40000; // ahem... not 10.

   from_vertices = new to_vertex_list_t[start_size];
   from_vertices_size = start_size;  // the size of the array
   max_from_vertex = -1;             // the maximum vertex encountered so far.

}

done_line_list_t::~done_line_list_t() {

   //cout << "destroying a done_line_list_t" << endl;
   //cout << "from_vertices_size is " << from_vertices_size << endl;
   //cout << "max_from_vertex is " << max_from_vertex << endl;
   delete [] from_vertices;
   //cout << "done deleting from_vertices" << endl;

}

void
done_line_list_t::resize_and_copy(int i) {

   //cout << "resize and copy to fix start vertex " << i << endl;
   //cout << "resize: from_vertices_size was " << from_vertices_size << endl;
   
   int new_size = int (rint(from_vertices_size +
			    (i - from_vertices_size + 500 )*1.5));

   //cout << "resize: new_size is " << new_size << endl;

   to_vertex_list_t *new_list = new to_vertex_list_t[new_size];

   for (int ii=0; i<max_from_vertex; i++)
      new_list[ii] = from_vertices[ii];

   max_from_vertex = i;
   delete [] from_vertices;      // out with the old
   from_vertices_size = new_size;// in with the new.
   from_vertices = new_list; 

   // cout << "resize: max_from_vertex is " << max_from_vertex << endl;
}

// We come here only when there is not already a mark.
// 
void
done_line_list_t::mark_as_done(int i, int j) {

   // need to mark both i and j indices first.
   //
   // First mark i first:
   //
   //
   //cout << "Marking as done: first way: " << i << "," << j << endl;
   //
   to_vertex_list_t *v = &from_vertices[i]; // don't copy!
   v->add(j);

   //cout << "Marking as done: secon way: " << j << "," << i << endl;
   //
   v = &from_vertices[j];  // don't copy
   v->add(i);

   max_from_vertex = (max_from_vertex > i) ? max_from_vertex : i;
   max_from_vertex = (max_from_vertex > j) ? max_from_vertex : j;

}

//
to_vertex_list_t
done_line_list_t::getVertex(unsigned int i) const {

   return from_vertices[i];

} 

//
void
to_vertex_list_t::add(int i) {

   //cout << "add: n_vertices is currently: " << n_vertices << endl;
   //cout << "add: adding vertex: " << i << endl;
   
   if ( n_vertices < vertex_list_size ) {
      
      //cout << "add: no need for a resize as " << n_vertices
      //	   << " < " << vertex_list_size << endl;
      vertex_list[n_vertices] = i;
      n_vertices++;      
      //cout << "add: now n_vertices is " << n_vertices << endl;
      
   } else {
      //cout << "add: vertex_list resizing" << endl;
     int new_size = vertex_list_size ? vertex_list_size + 2 : 4;
      int *new_list = new int[new_size ];
      // copy across the old data
      //
      for (int ii=0; ii<n_vertices; ii++)
	 new_list[ii] = vertex_list[ii];
      
      vertex_list_size = new_size;
      delete [] vertex_list;
      vertex_list = new_list;
      vertex_list[n_vertices] = i;
      n_vertices++;

      //cout << "add: vertex_list_size expanded to "
      //	   << vertex_list_size << endl;
      //cout << "add: now n_vertices is " << n_vertices << endl;
   }



}

// copy constuctor
//
to_vertex_list_t::to_vertex_list_t(const to_vertex_list_t &a) {

   //cout << "making a default to_vertex_list_t" << endl;
   Copy(a);
}

//      
void
to_vertex_list_t::Copy(const to_vertex_list_t &a) {
   
   // cout << "to_vertex_list_t Copy" << endl;
   //
   // int *new_vertex_list = new int[a.vertex_list_size];
   // for (int ii=0; ii<a.vertex_list_size; ii++)
   //    new_vertex_list[ii] = a.vertex_list[ii];

   n_vertices       = a.n_vertices;
   vertex_list_size = a.vertex_list_size;
   std::cout << "post Copy(): vertex_list_size = " << vertex_list_size << std::endl;
   std::cout << "post Copy(): n_vertices = " << n_vertices << std::endl;
      
}

// 
const to_vertex_list_t&
to_vertex_list_t::operator=(const to_vertex_list_t &a) {

   Copy(a);

   return *this;
}

   

// This is a question asked of the class
//
bool
to_vertex_list_t::contains(int i_test_vertex) {

   //cout << "looking for vertex " << i_test_vertex << " in "
   //	<< n_vertices << " vertices" << endl;
   for (int i=0; i< n_vertices; i++) {
      //cout << "for index i=" << i << " comparing " << i_test_vertex
      //   << " and " << vertex_list[i] << endl;
      if (vertex_list[i] == i_test_vertex) {
	 return 1;
	 
      }
   }

   return 0;
}

to_vertex_list_t::to_vertex_list_t() {

   vertex_list = 0; /* Allocate memory as needed for this list */

   vertex_list_size = 0;
   n_vertices = 0;
}

to_vertex_list_t::~to_vertex_list_t() {

   //cout << "~to_v_l_t: vertex_list_size is " << vertex_list_size << endl;
   //cout << "~to_v_l_t: n_vertices is " << n_vertices << endl;
   
   if (vertex_list_size > 0) { 
      delete [] vertex_list;
   }
   //cout << "~to_v_l_t: done deleting" << endl;

}



// Instantiate template(s)
// 
// template class CIsoSurface<short>;
// template class CIsoSurface<unsigned short>;
template class CIsoSurface<float>;


