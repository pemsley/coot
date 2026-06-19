/*
 * MoleculesToTriangles/CXXSurface/CXXTorusNode-gemmi.cc
 *
 * gemmi-native twin of CXXTorusNode.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "CXXTorusNode-gemmi.hh"

void coot::m2t::CXXTorusNode::init() {
   theAtom = 0;
   theta = omega = 0.;
   crd = CXXCoord<CXXCoord_ftype>(0., 0., 0.);
}
coot::m2t::CXXTorusNode::CXXTorusNode() { init(); }
coot::m2t::CXXTorusNode::CXXTorusNode(double inTheta, double inOmega) {
   init();
   theta = inTheta;
   omega = inOmega;
}
int coot::m2t::CXXTorusNode::setTheta(const double inTheta) { theta = inTheta; return 0; }
int coot::m2t::CXXTorusNode::setOmega(const double inOmega) { omega = inOmega; return 0; }
int coot::m2t::CXXTorusNode::setCoord(const CXXCoord<CXXCoord_ftype>&aCoord) { crd = aCoord; return 0; }
const CXXCoord<CXXCoord_ftype>& coot::m2t::CXXTorusNode::coord() const { return crd; }
const double coot::m2t::CXXTorusNode::getTheta() const { return theta; }
const double coot::m2t::CXXTorusNode::getOmega() const { return omega; }
int coot::m2t::CXXTorusNode::setAtom(const gemmi::Atom* anAtom) { theAtom = anAtom; return 0; }
const gemmi::Atom* coot::m2t::CXXTorusNode::getAtom() const { return theAtom; }
