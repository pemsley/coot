/*
 * MoleculesToTriangles/CXXClasses/DiscreteSegment-gemmi.cc
 *
 * gemmi-native twin of DiscreteSegment.cpp.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include <string>
#include "DiscreteSegment-gemmi.hh"
#include "ColorScheme-gemmi.hh"

extern "C" {
#include "lfit.h"
}

namespace {
   // Secondary structure is stamped on gemmi::Residue::flag ('H'/'E'/'L') by
   // coot::m2t::assign_secondary_structure (DSSP). 'E' == strand (smoothed).
   bool segment_residue_is_strand(const gemmi::CRA &cra) {
      return cra.residue && cra.residue->flag == 'E';
   }
}

void coot::m2t::DiscreteSegment::evaluateNormals()
{
   normalOnes.resize(calphas.size());
   normalTwos.resize(calphas.size());
   FCXXCoord xAxis(1., 0., 0.);
   FCXXCoord yAxis(0., 1., 0.);
   normalOnes[0] = xAxis;
   normalOnes[normalOnes.size()-1] = xAxis;
   normalTwos[0] = yAxis;
   normalTwos[normalOnes.size()-1] = yAxis;
   if (normalOnes.size() < 3) return;

   for (int i=1; i<(int)(calphas.size()-1); i++) {
      FCXXCoord a(calphaCoords[i-1]);
      FCXXCoord b(calphaCoords[i]);
      FCXXCoord c(calphaCoords[i+1]);
      FCXXCoord ab = b-a;
      FCXXCoord cb = b-c;
      FCXXCoord ac = c-a;
      FCXXCoord normal = ab^cb;
      if (normal*normalOnes[i-1] < 0.) normal *= -1.;
      normal.normalise();
      normalOnes[i] = normal;
      FCXXCoord normal2 = normal^ac;
      normal2.normalise();
      normalTwos[i] = normal2;
   }
   normalOnes[0] = normalOnes[1];
   normalOnes[normalOnes.size()-1] = normalOnes[normalOnes.size()-2];
   normalTwos[0] = normalTwos[1];
   normalTwos[normalTwos.size()-1] = normalTwos[normalTwos.size()-2];
}

void coot::m2t::DiscreteSegment::smoothBetas()
{
   float *xes = new float[calphas.size()];
   float *yes[3];
   float *n1s[3];
   float *n2s[3];
   for (int i=0; i<3; i++) {
      yes[i] = new float[calphas.size()];
      n1s[i] = new float[calphas.size()];
      n2s[i] = new float[calphas.size()];
   }
   int *fixeds = new int[calphas.size()];

   for (int iCalpha=0; iCalpha<(int)calphas.size(); iCalpha++) {
      int nInStrand = 0;
      while (segment_residue_is_strand(calphas[iCalpha+nInStrand])
             && iCalpha+nInStrand < (int)calphas.size()) {
         fixeds[nInStrand] = 0;
         xes[nInStrand] = nInStrand;
         for (int i=0; i<3; i++) {
            yes[i][nInStrand] = calphaCoords[nInStrand + iCalpha][i];
            n1s[i][nInStrand] = normalOnes[nInStrand + iCalpha][i];
            n2s[i][nInStrand] = normalTwos[nInStrand + iCalpha][i];
         }
         nInStrand++;
      }
      if (nInStrand >= 2) {
         fixeds[0] = fixeds[nInStrand-1] = 1;
         for (int i=0; i<3; i++) {
            int order = nInStrand - 1;
            if (order > 3) order = 3;
            ForcePolynomial(order, xes, yes[i], fixeds, nInStrand, 0, nInStrand-1);
            ForcePolynomial(order, xes, n1s[i], fixeds, nInStrand, 0, nInStrand-1);
            ForcePolynomial(order, xes, n2s[i], fixeds, nInStrand, 0, nInStrand-1);
         }
         for (int i=0; i<nInStrand; i++) {
            for (int j=0; j<3; j++) {
               calphaCoords[iCalpha + i][j] = yes[j][i];
               normalOnes[iCalpha + i][j] = n1s[j][i];
               normalTwos[iCalpha + i][j] = n2s[j][i];
            }
         }
      }
      iCalpha += nInStrand;
   }

   for (int i=0; i<3; i++) { delete [] yes[i]; delete [] n1s[i]; delete [] n2s[i]; }
   delete [] xes;
   delete [] fixeds;
}

void coot::m2t::DiscreteSegment::evaluateSplines(int samplesPerSegment)
{
   for (int i=0; i<(int)calphas.size(); i++) {
      coordSpline.addPair((float)i, calphaCoords[i]);
      normalOnesSpline.addPair((float)i, normalOnes[i]);
      normalTwosSpline.addPair((float)i, normalTwos[i]);
      auto [ax, ay, az] = anisoValues[i];
      anisoSpline.addPair((float)i, FCXXCoord(ax, ay, az, 1.0f));
   }
   coordSpline.calculateYDoublePrimes(1e30f, 1e30f, samplesPerSegment);
   normalOnesSpline.calculateYDoublePrimes(1e30f, 1e30f, samplesPerSegment);
   normalTwosSpline.calculateYDoublePrimes(1e30f, 1e30f, samplesPerSegment);
   anisoSpline.calculateYDoublePrimes(1e30f, 1e30f);
}

FCXXCoord coot::m2t::DiscreteSegment::coordFor(float xVal)     { return coordSpline.coordForXEquals(xVal); }
FCXXCoord coot::m2t::DiscreteSegment::normalOneFor(float xVal) { return normalOnesSpline.coordForXEquals(xVal); }
FCXXCoord coot::m2t::DiscreteSegment::normalTwoFor(float xVal) { return normalTwosSpline.coordForXEquals(xVal); }

void coot::m2t::DiscreteSegment::evaluateColors(std::shared_ptr<ColorScheme> colorScheme)
{
   for (int i=0; i<(int)calphas.size(); i++) {
      FCXXCoord color = colorScheme->colorForAtom(calphas[i]);
      colorSpline.addPair((float)i, color);
   }
   colorSpline.calculateYDoublePrimes(1e30f, 1e30f);
}

FCXXCoord coot::m2t::DiscreteSegment::colorFor(float xVal) { return colorSpline.coordForXEquals(xVal); }
