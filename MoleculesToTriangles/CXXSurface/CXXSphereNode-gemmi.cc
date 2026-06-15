/*
 * MoleculesToTriangles/CXXSurface/CXXSphereNode-gemmi.cc
 *
 * gemmi-native twin of CXXSphereNode.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "CXXSphereNode-gemmi.hh"
#include "CXXCoord.h"
#include "CXXTriangle-gemmi.hh"
#include <list>

coot::m2t::CXXSphereNode::CXXSphereNode() :
   theVertex(CXXCoord<CXXCoord_ftype>(0., 0., 0.)),
   theIntersector(0), shouldBeDrawn(1), theAtom(0) {}

coot::m2t::CXXSphereNode::CXXSphereNode(const CXXCoord<CXXCoord_ftype>&aCoord) :
   theVertex(aCoord), theIntersector(0), shouldBeDrawn(1), theAtom(0) {}

const CXXCoord<CXXCoord_ftype>& coot::m2t::CXXSphereNode::vertex() const { return theVertex; }
int coot::m2t::CXXSphereNode::doDraw() const { return shouldBeDrawn; }
int coot::m2t::CXXSphereNode::setDoDraw(const int yesNo) { shouldBeDrawn = yesNo; return shouldBeDrawn; }
int coot::m2t::CXXSphereNode::setVertex(const CXXCoord<CXXCoord_ftype>&crd) { theVertex = crd; return 0; }
void coot::m2t::CXXSphereNode::setIntersector(const CXXCircle *aCircle) { theIntersector = aCircle; }
const coot::m2t::CXXCircle *coot::m2t::CXXSphereNode::getIntersector() const { return theIntersector; }
