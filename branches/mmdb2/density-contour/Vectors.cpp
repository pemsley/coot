/* density-contour/Vectors.cpp
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

// #include "stdafx.h"  commented by Paul Emsley 26-12-2001 
#include "Vectors.h"

POINT3DXYZ operator+(const POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2)
{
	POINT3DXYZ result;

	result.x = pt3dPoint1.x + pt3dPoint2.x;
	result.y = pt3dPoint1.y + pt3dPoint2.y;
	result.z = pt3dPoint1.z + pt3dPoint2.z;

	return result;
}

POINT3DXYZ operator-(const POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2)
{
	POINT3DXYZ result;

	result.x = pt3dPoint1.x - pt3dPoint2.x;
	result.y = pt3dPoint1.y - pt3dPoint2.y;
	result.z = pt3dPoint1.z - pt3dPoint2.z;

	return result;
}

POINT3DXYZ operator*(const POINT3DXYZ& pt3dPoint, float fScale)
{
	POINT3DXYZ result;

	result.x = pt3dPoint.x*fScale;
	result.y = pt3dPoint.y*fScale;
	result.z = pt3dPoint.z*fScale;

	return result;
}

POINT3DXYZ operator*(float fScale, const POINT3DXYZ& pt3dPoint)
{
	POINT3DXYZ result;

	result.x = pt3dPoint.x*fScale;
	result.y = pt3dPoint.y*fScale;
	result.z = pt3dPoint.z*fScale;

	return result;
}

POINT3DXYZ operator/(const POINT3DXYZ& pt3dPoint, float fScale)
{
	POINT3DXYZ result;

	result.x = pt3dPoint.x/fScale;
	result.y = pt3dPoint.y/fScale;
	result.z = pt3dPoint.z/fScale;
	
	return result;
}

POINT3DXYZ& operator*=(POINT3DXYZ& pt3dPoint, float fScale)
{
	pt3dPoint.x *= fScale;
	pt3dPoint.y *= fScale;
	pt3dPoint.z *= fScale;

	return pt3dPoint;
}

POINT3DXYZ& operator/=(POINT3DXYZ& pt3dPoint, float fScale)
{
	pt3dPoint.x /= fScale;
	pt3dPoint.y /= fScale;
	pt3dPoint.z /= fScale;

	return pt3dPoint;
}

POINT3DXYZ& operator+=(POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2)
{
	pt3dPoint1.x += pt3dPoint2.x;
	pt3dPoint1.y += pt3dPoint2.y;
	pt3dPoint1.z += pt3dPoint2.z;

	return pt3dPoint1;
}

POINT3DXYZ& operator-=(POINT3DXYZ& pt3dPoint1, const POINT3DXYZ& pt3dPoint2)
{
	pt3dPoint1.x -= pt3dPoint2.x;
	pt3dPoint1.y -= pt3dPoint2.y;
	pt3dPoint1.z -= pt3dPoint2.z;
	
	return pt3dPoint1;
}
