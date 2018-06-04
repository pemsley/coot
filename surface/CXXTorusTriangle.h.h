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
 *  CXXTorusTriangle.h.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  
 *
 */

#ifndef CXX_mot_CXXTorusTriangle_included
#define CXX_mot_CXXTorusTriangle_included
#include <CXXTorusNode.h>
#include <CXXCoord.h>
#include <vector>

class CXXTorusTriangle {
private:
	static CXXTorusNode nodes[3];
public:
	CXXTorusTriangle();
	~CXXTorusTriangle();
	CXXTorusTriangle(const CXXTorusTriangle &another);
	CXXTorusTriangle(CXXTorusNode &, CXXTorusNode &, CXXTorusNode &);
	CXXTorusNode &element(int i) const;
};

#endif
