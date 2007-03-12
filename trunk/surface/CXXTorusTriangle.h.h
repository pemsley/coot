/*
 *  CXXTorusTriangle.h.h
 *  CXXSurface
 *
 *  Created by Martin Noble on Mon Feb 09 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CXXTorusTriangle_included
#define CXXTorusTriangle_included
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
