/*
 * MoleculesToTriangles/CXXSurface/CXXTorusNode-gemmi.hh
 *
 * gemmi-native twin of CXXTorusNode.h. Atom token const gemmi::Atom*.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXTorusNode_gemmi_included
#define CXXTorusNode_gemmi_included
#include "CXXCoord.h"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXTorusNode {
      private:
         const gemmi::Atom* theAtom;
         double theta;
         double omega;
         CXXCoord<CXXCoord_ftype> crd;
         void init();
      public:
         CXXTorusNode();
         CXXTorusNode(double the, double om);
         int setTheta(const double inTheta);
         int setOmega(const double inOmega);
         const CXXCoord<CXXCoord_ftype>& coord() const;
         const double getTheta() const;
         const double getOmega() const;
         int setCoord(const CXXCoord<CXXCoord_ftype>&);
         int setAtom(const gemmi::Atom* anAtom);
         const gemmi::Atom* getAtom() const;
      };
   }
}

#endif
