/*
 *  CXXTriangle.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Fri Jan 23 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
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
