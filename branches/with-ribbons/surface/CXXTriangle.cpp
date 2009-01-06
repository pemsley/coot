/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
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
/*
 *  CXXTriangle.cpp
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  
 *
 */

#include "CXXTriangle.h"

CXXTriangle::~CXXTriangle()
{
}

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

int CXXTriangle::i() const{
	return ijk[0];
}

int CXXTriangle::j() const{
	return ijk[1];
}

int CXXTriangle::k() const{
	return ijk[2];
}

int CXXTriangle::element(const int iElement) const{
	return ijk[iElement];
}

int CXXTriangle::setDoDraw(const int yesNo) {
	shouldBeDrawn = yesNo;
	return shouldBeDrawn;
}

int CXXTriangle::doDraw() const {
	return shouldBeDrawn;
}

CAtom *CXXTriangle::getAtom() const {
	return theAtom;
}

void CXXTriangle::setAtom(CAtom *anAtom) {
	theAtom = anAtom;
}

