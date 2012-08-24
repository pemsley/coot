/*
 *  CXXSurfaceVertex.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXSurfaceVertex.h"

CXX::CXXAlloc<CXXSurfaceVertex> CXXSurfaceVertex::allocator = CXX::CXXAlloc<CXXSurfaceVertex>();

void CXXSurfaceVertex::init(){
  pointers.reserve(1);
  vectors.reserve(4);
  scalars.reserve(2);
}

int CXXSurfaceVertex::setXyz(unsigned int coordType, double *xyz){
	double xyzr[4];
	
	for (int i=0; i<3; i++) xyzr[i] = xyz[i];
	xyzr[3] = 0.;
	setXyzr(coordType, xyzr);
	return 0;
}


int CXXSurfaceVertex::setXyzr(unsigned int coordType, double *xyzr){
	if (coordType <= vectors.size()){
		vectors[coordType-1] = CXXCoord(xyzr[0],xyzr[1],xyzr[2],xyzr[3]);
	}
	else {
		vectors.resize(coordType);
		vectors[coordType-1] =  CXXCoord(xyzr[0],xyzr[1],xyzr[2],xyzr[3]);
	}
	return 0;
}

int CXXSurfaceVertex::setCoord(unsigned int coordType, const CXXCoord &crd){
	if (coordType <= vectors.size()){
		vectors[coordType-1] = crd;
	}
	else {
		vectors.resize(coordType);
		vectors[coordType-1] = crd;
	}
	return 0;
}

double CXXSurfaceVertex::x(unsigned int coordType){
	if (coordType>0 && coordType <= vectors.size()){
		return vectors[coordType-1].element(0);
	}
	else return 0.;
}

double CXXSurfaceVertex::y(unsigned int coordType){
	if (coordType>0 && coordType <= vectors.size()){
		return vectors[coordType-1].element(1);
	}
	else return 0.;
}

double CXXSurfaceVertex::z(unsigned int coordType){
	if (coordType>0 && coordType <= vectors.size()){
		return vectors[coordType-1].element(2);
	}
	else return 0.;
}

double CXXSurfaceVertex::r(unsigned int coordType){
	if (coordType>0 && coordType <= vectors.size()){
		return vectors[coordType-1].element(3);
	}
	else return 0.;
}

CXXCoord_ftype *CXXSurfaceVertex::xyzPntr(unsigned int coordType){
	if (coordType>0 && coordType <= vectors.size()){
		return vectors[coordType-1].xyzPntr();
	}
	else return 0;
}

int CXXSurfaceVertex::setScalar(unsigned int scalarType, double value){
	if (scalarType <= scalars.size()){
		scalars[scalarType-1] = value;
	}
	else {
		scalars.resize(scalarType);
		scalars[scalarType-1] = value;
	}
	return 0;
}

int CXXSurfaceVertex::setPointer(unsigned int pointerType, void *value){
	if (pointerType <= pointers.size()){
		pointers[pointerType-1] = value;
	}
	else {
		pointers.resize(pointerType);
		pointers[pointerType-1] = value;
	}
	return 0;
}

void *CXXSurfaceVertex::pointer(unsigned int pointerType) const{
	if (pointerType<1||pointerType>pointers.size()){
		return 0;
	}
	else {
		return pointers[pointerType-1];
	}
}

const CXXCoord &CXXSurfaceVertex::coordRef(unsigned int coordType) const{
	return vectors[coordType-1];
}

double CXXSurfaceVertex::scalar(unsigned int scalarType) const{
	if (scalarType<1||scalarType>scalars.size()){
		return 0;
	}
	else {
		return scalars[scalarType-1];
	}
}

const CXXCoord *CXXSurfaceVertex::coordPntr(unsigned int vectorType) const{
	if (vectorType<1||vectorType>vectors.size()){
		return 0;
	}
	else {
		return &vectors[vectorType-1];
	}
}

int CXXSurfaceVertex::nPointers() const{
	return pointers.size();
}

int CXXSurfaceVertex::nVectors() const{
	return vectors.size();
}

int CXXSurfaceVertex::nScalars() const{
	return scalars.size();
}
