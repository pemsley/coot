/*
 *  CXXTriangle.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "CXXTriangle.h"

int CXXTriangle::setIjk(const int *lijk){
	for (int i=0; i<3; i++) ijk[i] = lijk[i];
	return 0;
}

int CXXTriangle::setIjk(const int i, const int j, const int k){
	ijk[0] = i;
	ijk[1] = j;
	ijk[2] = k;
	return 0;
}

int CXXTriangle::getIjk(int *lijk) const{
	for (int i=0; i<3; i++) lijk[i] = ijk[i];
	return 0;
}

int CXXTriangle::setDoDraw(const int yesNo) {
	shouldBeDrawn = yesNo;
	return shouldBeDrawn;
}

CAtom *CXXTriangle::getAtom() const {
	return theAtom;
}

void CXXTriangle::setAtom(CAtom *anAtom) {
	theAtom = anAtom;
}

