/*
 * MoleculesToTriangles/CXXClasses/DiscreteSegment-gemmi.hh
 *
 * gemmi-native twin of DiscreteSegment.h (mmdb->gemmi migration). Stores a
 * gemmi::CRA per C-alpha (carries chain/residue context for colour and SSE).
 * Lives alongside the original; type is coot::m2t::DiscreteSegment.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef DISCRETE_SEGMENT_GEMMI_HH
#define DISCRETE_SEGMENT_GEMMI_HH

#include <gemmi/model.hpp>   // gemmi::CRA, gemmi::Atom
#include "NRStuff.h"

#include <vector>
#include <tuple>
#include <map>
#include <memory>

namespace coot {
   namespace m2t {

      class ColorScheme;  // gemmi colour scheme (Phase 3)
      class ColorRule;

      class DiscreteSegment {
      private:
         std::vector<gemmi::CRA> calphas;
         std::vector<FCXXCoord> calphaCoords;
         std::vector<FCXXCoord> normalOnes;
         std::vector<FCXXCoord> normalTwos;
         std::vector<std::tuple<float, float, float>> anisoValues;
         CoordSpline coordSpline;
         CoordSpline normalOnesSpline;
         CoordSpline normalTwosSpline;
         CoordSpline anisoSpline;
         CoordSpline colorSpline;
      public:
         ~DiscreteSegment() {
            calphas.clear();
            normalOnes.clear();
            normalTwos.clear();
            anisoValues.clear();
         }
         void addCalpha(const gemmi::CRA &calpha) {
            calphas.push_back(calpha);
            const gemmi::Position &p = calpha.atom->pos;
            calphaCoords.push_back(FCXXCoord(p.x, p.y, p.z, 1.0f));
            anisoValues.push_back(std::make_tuple(1.0f, 1.0f, 1.0f));
         }
         void addCalpha(const gemmi::CRA &calpha, float radius) {
            calphas.push_back(calpha);
            const gemmi::Position &p = calpha.atom->pos;
            calphaCoords.push_back(FCXXCoord(p.x, p.y, p.z, radius));
            anisoValues.push_back(std::make_tuple(1.0f, 1.0f, 1.0f));
         }
         void addCalpha(const gemmi::CRA &calpha, float ax, float ay, float az) {
            calphas.push_back(calpha);
            const gemmi::Position &p = calpha.atom->pos;
            calphaCoords.push_back(FCXXCoord(p.x, p.y, p.z, 1.0f));
            anisoValues.push_back(std::make_tuple(ax, ay, az));
         }
         std::tuple<float, float, float> anisoFor(float xVal) {
            FCXXCoord aniso = anisoSpline.coordForXEquals(xVal);
            return std::make_tuple(aniso.x(), aniso.y(), aniso.z());
         }
         FCXXCoord operator [] (int i) {
            const gemmi::Position &p = calphas[i].atom->pos;
            return FCXXCoord(p.x, p.y, p.z);
         }
         int nCalphas() { return int(calphas.size()); }
         gemmi::CRA calpha(int iCalpha) { return calphas[iCalpha]; }
         void evaluateNormals();
         void smoothBetas();
         void evaluateSplines(int samplesPerSegment = 6);
         // gemmi: no selection-handle map - the colour scheme tests matches(CRA) itself.
         void evaluateColors(std::shared_ptr<ColorScheme> colorScheme);
         FCXXCoord coordFor(float xVal);
         FCXXCoord normalOneFor(float xVal);
         FCXXCoord normalTwoFor(float xVal);
         FCXXCoord colorFor(float xVal);
      };
   }
}

#endif // DISCRETE_SEGMENT_GEMMI_HH
