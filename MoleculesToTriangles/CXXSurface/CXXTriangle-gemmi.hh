/*
 * MoleculesToTriangles/CXXSurface/CXXTriangle-gemmi.hh
 *
 * gemmi-native twin of CXXTriangle.h. Atom-identity token is const gemmi::Atom*.
 * Type is coot::m2t::CXXTriangle. CXXSurfaceVertex is mmdb-free (reused).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXTriangle_gemmi_included
#define CXXTriangle_gemmi_included

#include "CXXSurfaceVertex.h"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXTriangle {
      private:
         // friend class CXXFlatTriangle; // (re-add when the sphere-flat-triangle twin lands)
         size_t ijk[4];
         const gemmi::Atom* theAtom;
         int shouldBeDrawn;
      public:
         CXXTriangle() : theAtom(0), shouldBeDrawn(1) {
            ijk[0]=ijk[1]=ijk[2]=ijk[3]=0;
         }
         CXXTriangle(const size_t &i, const size_t &j, const size_t &k) : theAtom(0), shouldBeDrawn(1) {
            ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=0;
         }
         CXXTriangle(const size_t &i, const size_t &j, const size_t &k, const size_t &l) : theAtom(0), shouldBeDrawn(1) {
            ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=l;
         }
         CXXTriangle(const size_t &i, const size_t &j, const size_t &k, const gemmi::Atom* anAtom) : theAtom(anAtom), shouldBeDrawn(1) {
            ijk[0]=i; ijk[1]=j; ijk[2]=k; ijk[3]=0;
         }
         CXXTriangle(const int *ijk_in) : theAtom(0), shouldBeDrawn(1) {
            for (int i=0; i<3; i++) ijk[i] = ijk_in[i];
         }
         int setIjk(const size_t i, const size_t j, const size_t k);
         int getIjk(size_t *lijk) const;
         int setIjk(const size_t *lijk);
         int setDoDraw(const int);
         int doDraw() const { return shouldBeDrawn; }
         const gemmi::Atom* getAtom() const;
         void setAtom(const gemmi::Atom* anAtom);
         void setElement(const size_t target, const size_t value) { ijk[target] = value; }
         const size_t& operator [] (size_t element) const { return ijk[element]; }
         size_t& operator [] (size_t element) { return ijk[element]; }
         static bool doNotDraw(CXXTriangle &aNode) { return aNode.doDraw() == 0; }
      };
   }
}

#endif
