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

class CXXSurfaceVertex;

class CXXTriangle {
private:
	friend class CXXFlatTriangle;
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
	int setIjk(const int i, const int j, const int k);
	int getIjk(int *lijk) const;
	int setIjk(const int *lijk);
	int getIjk(int *i, int *j, int *k) const;
	int *ijkPntr();
	int setDoDraw(const int);
	int doDraw() const {
		return shouldBeDrawn;
	};
	CAtom *getAtom() const;
	void setAtom(CAtom *anAtom);
	void setElement(const int target, const int value) {
		ijk[target] = value;
	};
    const int& operator [] (unsigned element) const {
		return ijk[element];
	};
    int& operator [] (unsigned element) {
		return ijk[element];
	};
    static bool doNotDraw(CXXTriangle &aNode){
        return aNode.doDraw() == 0;
    };
};

#endif
