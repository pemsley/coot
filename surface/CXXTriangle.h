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
 *  CXXTriangle.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  
 *
 */
#ifndef CXXTriangle_included
#define CXXTriangle_included

#include "CXXSurfaceVertex.h"

class CXXSurfaceVertex;

class CXXTriangle {
private:
	int ijk[4];
	class CAtom *theAtom;
	int shouldBeDrawn;
public:
	CXXTriangle() : theAtom(0), shouldBeDrawn(1) {
	  ijk[0]=ijk[1]=ijk[2]=ijk[3]=0;
	};
	CXXTriangle(const int &i, const int &j, const int &k) : theAtom(0), shouldBeDrawn(1){
	  ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=0;
	};
	CXXTriangle(const int &i, const int &j, const int &k, const int &l) : theAtom(0), shouldBeDrawn(1){
	  ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=l;
	};
	CXXTriangle(const int &i, const int &j, const int &k, CAtom *anAtom) : theAtom(anAtom), shouldBeDrawn(1){
	  ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=0;
	};
	CXXTriangle(const int *ijk_in) : theAtom(0), shouldBeDrawn(1){ 
	  for (int i=0; i<3; i++) ijk[i] = ijk_in[i];
	};
	~CXXTriangle();
	int setIjk(const int i, const int j, const int k);
	int getIjk(int *lijk) const;
	int setIjk(const int *lijk);
	int getIjk(int *i, int *j, int *k) const;
	int *ijkPntr();
	int i() const;
	int j() const;
	int k() const;
	int element(const int iElement) const;
	int setDoDraw(const int);
	int doDraw() const;
	CAtom *getAtom() const;
	void setAtom(CAtom *anAtom);
	int operator [] (int element) const{
		return ijk[element];
	};
	void setElement(const int target, const int value) {
		ijk[target] = value;
	};
	
};

#endif
