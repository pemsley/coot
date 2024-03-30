/*
 * MoleculesToTriangles/CXXSurface/CXXTriangle.h
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
#ifndef CXXTriangle_included
#define CXXTriangle_included

#include "CXXSurfaceVertex.h"
#include "mmdb2/mmdb_manager.h"

class CXXSurfaceVertex;

class CXXTriangle {
private:
	friend class CXXFlatTriangle;
	size_t ijk[4];
	mmdb::Atom* theAtom;
	int shouldBeDrawn;
public:
	CXXTriangle() : theAtom(0), shouldBeDrawn(1) {
	  ijk[0]=ijk[1]=ijk[2]=ijk[3]=0;
	};
	CXXTriangle(const size_t &i, const size_t &j, const size_t &k) : theAtom(0), shouldBeDrawn(1){
	  ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=0;
	};
	CXXTriangle(const size_t &i, const size_t &j, const size_t &k, const size_t &l) : theAtom(0), shouldBeDrawn(1){
	  ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=l;
	};
	CXXTriangle(const size_t &i, const size_t &j, const size_t &k, mmdb::Atom* anAtom) : theAtom(anAtom), shouldBeDrawn(1){
	  ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=0;
	};
	CXXTriangle(const int *ijk_in) : theAtom(0), shouldBeDrawn(1){ 
	  for (int i=0; i<3; i++) ijk[i] = ijk_in[i];
	};
	int setIjk(const size_t i, const size_t j, const size_t k);
	int getIjk(size_t *lijk) const;
	int setIjk(const size_t *lijk);
	int getIjk(size_t *i, size_t *j, size_t *k) const;
	size_t *ijkPntr();
	int setDoDraw(const int);
	int doDraw() const {
		return shouldBeDrawn;
	};
	mmdb::Atom* getAtom() const;
	void setAtom(mmdb::Atom* anAtom);
	void setElement(const size_t target, const size_t value) {
		ijk[target] = value;
	};
    const size_t& operator [] (size_t element) const {
		return ijk[element];
	};
    size_t& operator [] (size_t element) {
		return ijk[element];
	};
    static bool doNotDraw(CXXTriangle &aNode){
        return aNode.doDraw() == 0;
    };
};

#endif
