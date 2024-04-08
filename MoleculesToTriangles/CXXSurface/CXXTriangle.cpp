/*
 * MoleculesToTriangles/CXXSurface/CXXTriangle.cpp
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

#include "CXXTriangle.h"

int CXXTriangle::setIjk(const size_t *lijk){
	for (int i=0; i<3; i++) ijk[i] = lijk[i];
	return 0;
}

int CXXTriangle::setIjk(const size_t i, const size_t j, const size_t k){
	ijk[0] = i;
	ijk[1] = j;
	ijk[2] = k;
	return 0;
}

int CXXTriangle::getIjk(size_t *lijk) const{
	for (int i=0; i<3; i++) lijk[i] = ijk[i];
	return 0;
}

int CXXTriangle::setDoDraw(const int yesNo) {
	shouldBeDrawn = yesNo;
	return shouldBeDrawn;
}

mmdb::Atom* CXXTriangle::getAtom() const {
	return theAtom;
}

void CXXTriangle::setAtom(mmdb::Atom* anAtom) {
	theAtom = anAtom;
}

