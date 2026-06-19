/*
 * MoleculesToTriangles/CXXSurface/CXXSurfaceMaker-gemmi.hh
 *
 * gemmi-native twin of CXXSurfaceMaker.h. Driven by a gemmi::Structure/Model +
 * selection sets (std::set<const gemmi::Atom*>) instead of an mmdb::Manager + selHnd.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef CXXSurfaceMaker_gemmi_included
#define CXXSurfaceMaker_gemmi_included

#include "CXXSurface-gemmi.hh"
#include <vector>
#include <set>
#include <map>
#include <gemmi/model.hpp>

namespace coot {
   namespace m2t {

      class CXXVCN;
      class CXXBall;
      class CXXCircleNode;
      class CXXSphereElement;

      class CXXSurfaceMaker {
      private:
         gemmi::Structure *allAtomsManager;
         std::vector<CXXSurface> childSurfaces;
      public:
         CXXSurfaceMaker() { allAtomsManager = 0; }
         ~CXXSurfaceMaker() {}
         // Primary gemmi entry point (molecular / accessible surface of a selection).
         CXXSurfaceMaker(gemmi::Structure *structure, gemmi::Model &model,
                         const std::set<const gemmi::Atom*> &selSet,
                         const std::set<const gemmi::Atom*> &contextSet,
                         double delta, double probeRadius, bool blend_edges);

         std::vector<CXXSurface> &getChildSurfaces() { return childSurfaces; }
         gemmi::Structure *getStructure() const { return allAtomsManager; }

         int calculateFromAtoms(gemmi::Structure *structure, gemmi::Model &model,
                                const std::set<const gemmi::Atom*> &selSet,
                                const std::set<const gemmi::Atom*> &contextSet,
                                double delta, double probeRadius, bool blend_edges);

         // VdW / solvent-accessible (ball-only, no reentrant probes)
         int calculateVDWFromAtoms(gemmi::Structure *structure, gemmi::Model &model,
                                   const std::set<const gemmi::Atom*> &selSet,
                                   const std::set<const gemmi::Atom*> &contextSet,
                                   double delta, double probeRadius, double radiusMultiplier);
         int calculateAccessibleFromAtoms(gemmi::Structure *structure, gemmi::Model &model,
                                          const std::set<const gemmi::Atom*> &selSet,
                                          const std::set<const gemmi::Atom*> &contextSet,
                                          double delta, double probeRadius, double radiusMultiplier);

         SurfaceParameters measuredProperties();
         double getAtomRadius(const gemmi::Atom *atom, const gemmi::Residue &residue);
         int assignAtom(gemmi::Structure *structure, const std::set<const gemmi::Atom*> &selSet);

         void generateArrays(std::vector<CXXVCN> &vcns, std::vector<CXXCoord<float> > &accessibles,
                             std::vector<const gemmi::Atom*> &atoms);

         void memberHandleCentralAtoms(
            const int atomNr,
            const std::vector<const CXXBall *> *vdwBallPntrs,
            CXXSurface *elementSurfacesArray,
            const float radiusMultiplier,
            const float probeRadius,
            const float delta,
            const CXXSphereElement *unitSphereAtOrigin,
            const std::map<const CXXBall *, std::vector<const CXXBall *> > *contactMap,
            std::vector<CXXCircleNode> *splitReentrantProbesArray,
            const std::set<const gemmi::Atom*> &selSet);
      };
   }
}

#endif
