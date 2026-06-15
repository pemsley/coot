/*
 * MoleculesToTriangles/CXXClasses/SurfacePrimitive-gemmi.hh
 *
 * gemmi-native twin of SurfacePrimitive.h. Consumes the coot::m2t CXXSurfaceMaker
 * (driven by gemmi selection sets) and bakes per-vertex colour from the coot::m2t
 * ColorScheme via each vertex's atom -> CRA -> colorForAtom. Renders as the reused
 * VertexColorNormalPrimitive (atomArray left null - the colour is baked in).
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */
#ifndef SurfacePrimitive_gemmi_hh
#define SurfacePrimitive_gemmi_hh
#include <memory>
#include <set>
#include "VertexColorNormalPrimitive.h"
#include "ColorScheme-gemmi.hh"
#include "MoleculesToTriangles/CXXSurface/CXXSurfaceMaker-gemmi.hh"
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class SurfacePrimitive : public VertexColorNormalPrimitive {
      private:
         CXXSurfaceMaker *cxxSurfaceMaker;
         std::shared_ptr<ColorScheme> colorScheme;
      public:
         enum SurfaceType { AccessibleSurface, VdWSurface, MolecularSurface };
         SurfacePrimitive();
         // selSet = central (chunk) atoms, contextSet = full selection (neighbours)
         SurfacePrimitive(gemmi::Structure *structure, gemmi::Model &model,
                          const std::set<const gemmi::Atom*> &selSet,
                          const std::set<const gemmi::Atom*> &contextSet,
                          std::shared_ptr<ColorScheme> _colorScheme, enum SurfaceType type,
                          float probeRadius, float radiusMultiplier);
         virtual ~SurfacePrimitive() {
            if (cxxSurfaceMaker) delete cxxSurfaceMaker;
         }
         virtual void generateArrays();
         CXXSurfaceMaker *getCXXSurfaceMaker() { return cxxSurfaceMaker; }
      };
   }
}

#endif
