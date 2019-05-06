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
 *  CXXTorusTriangle.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  
 *
 */
#ifndef CXX_mot_CXXTorusTriangle_included
#define CXX_mot_CXXTorusTriangle_included
#include "CXXTorusNode.h"

class CXXTorusElement;

class CXXTorusTriangle {
private:
	CXXTorusElement *theTorusElement;
	int nodes[3];
	double deltaThetas[3];
	double deltaOmegas[3];
	void init();
public:
		CXXTorusTriangle();
	CXXTorusTriangle(const CXXTorusTriangle &target);
	CXXTorusTriangle(CXXTorusElement *aParent, 
					 int, int, int,
					 double, double, double, double, double, double);
	void setNode(const int i, int aNode);
	int getNode(const int i) const;
	CXXTorusElement *getTorusElement() const; 
	void setTorusElement(CXXTorusElement *anElement); 
	int bisect(double angleInDegrees);
	int tooLongEdge (double angle);
	double getDeltaTheta(const int iEdge) const;
	double getDeltaOmega(const int iEdge) const;
	void setDeltaTheta(const int iEdge, const double newTheta);
	void setDeltaOmega(const int iEdge, const double newOmega);
	void selfReport() const;
};

#endif
