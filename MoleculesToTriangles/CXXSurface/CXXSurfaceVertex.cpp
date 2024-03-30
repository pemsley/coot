/*
 * MoleculesToTriangles/CXXSurface/CXXSurfaceVertex.cpp
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

#include "CXXSurfaceVertex.h"

void CXXSurfaceVertex::init(){
  pointers.reserve(1);
  vectors.reserve(4);
  scalars.reserve(2);
}

int CXXSurfaceVertex::setXyz(size_t coordType, double *xyz){
	double xyzr[4];
	
	for (int i=0; i<3; i++) xyzr[i] = xyz[i];
	xyzr[3] = 0.;
	setXyzr(coordType, xyzr);
	return 0;
}


int CXXSurfaceVertex::setXyzr(size_t coordType, double *xyzr){
	if (coordType < vectors.size()){
		vectors[coordType] = CXXCoord<CXXCoord_ftype>(xyzr[0],xyzr[1],xyzr[2],xyzr[3]);
	}
	else {
		vectors.resize(coordType+1);
		vectors[coordType] =  CXXCoord<CXXCoord_ftype>(xyzr[0],xyzr[1],xyzr[2],xyzr[3]);
	}
	return 0;
}

int CXXSurfaceVertex::setCoord(size_t coordType, const CXXCoord<CXXCoord_ftype>&crd){
	if (coordType < vectors.size()){
		vectors[coordType] = crd;
	}
	else {
		vectors.resize(coordType+1);
		vectors[coordType] = crd;
	}
	return 0;
}

double CXXSurfaceVertex::x(size_t coordType){
	if (coordType < vectors.size()){
		return vectors[coordType].element(0);
	}
	else return 0.;
}

double CXXSurfaceVertex::y(size_t coordType){
	if (coordType < vectors.size()){
		return vectors[coordType].element(1);
	}
	else return 0.;
}

double CXXSurfaceVertex::z(size_t coordType){
	if (coordType < vectors.size()){
		return vectors[coordType].element(2);
	}
	else return 0.;
}

double CXXSurfaceVertex::r(size_t coordType){
	if (coordType < vectors.size()){
		return vectors[coordType].element(3);
	}
	else return 0.;
}

CXXCoord_ftype *CXXSurfaceVertex::xyzPntr(size_t coordType){
	if (coordType < vectors.size()){
		return vectors[coordType].xyzPntr();
	}
	else return 0;
}

int CXXSurfaceVertex::setScalar(size_t scalarType, double value){
	if (scalarType < scalars.size()){
		scalars[scalarType] = value;
	}
	else {
		scalars.resize(scalarType+1);
		scalars[scalarType] = value;
	}
	return 0;
}

int CXXSurfaceVertex::setPointer(size_t pointerType, void *value){
	if (pointerType < pointers.size()){
		pointers[pointerType] = value;
	}
	else {
		pointers.resize(pointerType+1);
		pointers[pointerType] = value;
	}
	return 0;
}

void *CXXSurfaceVertex::pointer(size_t pointerType) const{
	if (pointerType>=pointers.size()){
		return 0;
	}
	else {
		return pointers[pointerType];
	}
}

const CXXCoord<CXXCoord_ftype>&CXXSurfaceVertex::coordRef(size_t coordType) const{
	return vectors[coordType];
}

double CXXSurfaceVertex::scalar(size_t scalarType) const{
	if (scalarType>=scalars.size()){
		return 0;
	}
	else {
		return scalars[scalarType];
	}
}

const CXXCoord<CXXCoord_ftype>*CXXSurfaceVertex::coordPntr(size_t vectorType) const{
	if (vectorType>=vectors.size()){
		return 0;
	}
	else {
		return &vectors[vectorType];
	}
}

size_t CXXSurfaceVertex::nPointers() const{
	return pointers.size();
}

size_t CXXSurfaceVertex::nVectors() const{
	return vectors.size();
}

size_t CXXSurfaceVertex::nScalars() const{
	return scalars.size();
}
