/* density-contour/Vectors.h
 * 
 * Copyright 2000 Raghavendra Chandrashekara
 * Copyright 2005 The University of York
 * 
 * Author: Raghavendra Chandrashekara, Paul Bourke and Cory Gene Bloyd
 *         Paul Emsley and Kevin Cowtan
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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
// File Name: Vectors.h
// Last Modified: 9/8/2000
// Author: Raghavendra Chandrashekara
// Email: rc99@doc.ic.ac.uk, rchandrashekara@hotmail.com
//
// Description: This is the implementation file for POINT3DXYZ class.

#ifndef VECTORS_H
#define VECTORS_H
// File Name: Vectors.h
// Last Modified: 7/8/2000
// Author: Raghavendra Chandrashekara
// Email: rc99@doc.ic.ac.uk, rchandrashekara@hotmail.com
//
// Description: This file contains some useful structures.

typedef float POINT3D[3];
typedef float VECTOR3D[3];

struct POINT3DXYZ {
	float x, y, z;
	friend POINT3DXYZ operator+(const POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2);
	friend POINT3DXYZ operator-(const POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2);
	friend POINT3DXYZ operator*(const POINT3DXYZ& pt3dPoint, float fScale);
	friend POINT3DXYZ operator*(float fScale, const POINT3DXYZ& pt3dPoint);
	friend POINT3DXYZ operator/(const POINT3DXYZ& pt3dPoint, float fScale);
	friend POINT3DXYZ& operator*=(POINT3DXYZ& pt3dPoint, float fScale);
	friend POINT3DXYZ& operator/=(POINT3DXYZ& pt3dPoint, float fScale);
	friend POINT3DXYZ& operator+=(POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2);
	friend POINT3DXYZ& operator-=(POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2);

};

typedef POINT3DXYZ VECTOR3DXYZ;
#endif // VECTORS_H
