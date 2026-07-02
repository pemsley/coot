/*
 * MoleculesToTriangles/CXXSurface/CXXTriangle-gemmi.cc
 *
 * gemmi-native twin of CXXTriangle.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include "CXXTriangle-gemmi.hh"

int coot::m2t::CXXTriangle::setIjk(const size_t *lijk) {
   for (int i=0; i<3; i++) ijk[i] = lijk[i];
   return 0;
}
int coot::m2t::CXXTriangle::setIjk(const size_t i, const size_t j, const size_t k) {
   ijk[0] = i; ijk[1] = j; ijk[2] = k;
   return 0;
}
int coot::m2t::CXXTriangle::getIjk(size_t *lijk) const {
   for (int i=0; i<3; i++) lijk[i] = ijk[i];
   return 0;
}
int coot::m2t::CXXTriangle::setDoDraw(const int yesNo) {
   shouldBeDrawn = yesNo;
   return shouldBeDrawn;
}
const gemmi::Atom* coot::m2t::CXXTriangle::getAtom() const {
   return theAtom;
}
void coot::m2t::CXXTriangle::setAtom(const gemmi::Atom* anAtom) {
   theAtom = anAtom;
}
