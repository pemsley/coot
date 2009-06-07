/* 
 * 
 * Copyright 2004 by The University of Oxford
 * Author: Martin Noble, Jan Gruber
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
 *  CXXTorusNode.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  
 *
 */

#ifndef CXXTorusNode_included
#define CXXTorusNode_included
#include "CXXCoord.h"

class CAtom;

class CXXTorusNode{
private:
	CAtom *theAtom;
	double theta;
	double omega;
	CXXCoord crd;
	int ind;
	void init();
public:
		CXXTorusNode();
	~CXXTorusNode();
	CXXTorusNode(double the, double om);
	int setTheta(const double inTheta);
	int setOmega(const double inOmega);
	const CXXCoord &coord() const;
	const double getTheta() const;
	const double getOmega() const;
	int setIndex(int i);
	const int index(void) const;
	int setCoord(const CXXCoord &);
	int setAtom(CAtom *anAtom);
	CAtom *getAtom() const;
	CXXTorusNode(const CXXTorusNode &oldOne){
	  theAtom = oldOne.getAtom();
	  crd = oldOne.coord();
	  omega = oldOne.getOmega();
	  theta = oldOne.getTheta();
	  ind = oldOne.index();
	}
};

#endif
