/*
 * MoleculesToTriangles/CXXClasses/AtomPropertyRampColorRule-gemmi.hh
 *
 * gemmi-native twin of AtomPropertyRampColorRule.h. Same HSV/spline ramp maths;
 * colorForAtom reads the residue seq number / atom B-factor from a gemmi::CRA.
 * Type is coot::m2t::AtomPropertyRampColorRule.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef AtomPropertyRampColorRule_gemmi_hh
#define AtomPropertyRampColorRule_gemmi_hh

#include "ColorRule-gemmi.hh"
#include "NRStuff.h"
#include <cmath>
#include <memory>

namespace coot {
   namespace m2t {

      class AtomPropertyRampColorRule : public ColorRule {
      private:
         FCXXCoord startHSV, middleHSV, endHSV;
         float startValue, endValue;
         int rampType;
         CoordSpline spline;
         int ramp_points_size;
      public:
         enum RampType { BFactor, ResidueNumber };

         AtomPropertyRampColorRule() : startValue(1.), endValue(1000.0) {
            rank = 1.0;
            ramp_points_size = 6;
            type = AtomPropertyRamp;
            rampType = ResidueNumber;
            compoundSelection = std::make_shared<compound_selection_t>("/*/*/*.*/*");
            startHSV  = FCXXCoord(240., 0.6, 0.6, 1.);
            middleHSV = FCXXCoord(120., 0.6, 0.6, 1.);
            endHSV    = FCXXCoord(0.,   0.6, 0.6, 1.);
            updateSpline();
         }
         float getStartValue() const { return startValue; }
         void setStartValue(const float value) { startValue = value; updateSpline(); }
         float getEndValue() const { return endValue; }
         void setEndValue(float value) { endValue = value; updateSpline(); }
         FCXXCoord getStartRGB() { return hsvToRGB(startHSV); }
         void setStartRGB(FCXXCoord color) {
            startHSV = rgbToHSV(color);
            while (middleHSV[0] > startHSV[0]) middleHSV[0] -= 360.;
            while (endHSV[0] > middleHSV[0]) endHSV[0] -= 360.;
            updateSpline();
         }
         FCXXCoord getMiddleRGB() { return hsvToRGB(middleHSV); }
         void setMiddleRGB(FCXXCoord color) {
            middleHSV = rgbToHSV(color);
            while (middleHSV[0] > startHSV[0]) middleHSV[0] -= 360.;
            while (endHSV[0] > middleHSV[0]) endHSV[0] -= 360.;
            updateSpline();
         }
         FCXXCoord getEndRGB() { return hsvToRGB(endHSV); }
         void setEndRGB(FCXXCoord color) {
            endHSV = rgbToHSV(color);
            while (middleHSV[0] > startHSV[0]) middleHSV[0] -= 360.;
            while (endHSV[0] > middleHSV[0]) endHSV[0] -= 360.;
            updateSpline();
         }
         int getRampType() { return rampType; }
         void setRampType(int type_in) { rampType = type_in; }
         void setNumberOfRampPoints(int _accu) { ramp_points_size = _accu; }

         void updateSpline() {
            spline.clearSpline();
            for (int i=0; i<ramp_points_size; i++) {
               float f1 = (1.0 - float(i)/ramp_points_size);
               float f2 = (float(i)/ramp_points_size);
               FCXXCoord hsv = middleHSV;
               if (float(i)/ramp_points_size < 0.5) {
                  hsv[0] = f1*startHSV[0] + f2*middleHSV[0];
               } else {
                  hsv[0] = f1*middleHSV[0] + f2*endHSV[0];
               }
               spline.addPair(f1*startValue+f2*endValue, hsv);
            }
            spline.calculateYDoublePrimes(1e30, 1e30);
         }

         FCXXCoord hsvToRGB(FCXXCoord hsv) {
            FCXXCoord rgb;
            rgb[3] = hsv[3];
            int i;
            float f, p, q, t;
            if (hsv[1] == 0.) { rgb[0] = rgb[1] = rgb[2] = hsv[2]; return rgb; }
            hsv[0] /= 60.f;
            i = fmod(floor(hsv[0]), 6.f);
            f = hsv[0] - i;
            p = hsv[2] * (1 - hsv[1]);
            q = hsv[2] * (1 - hsv[1] * f);
            t = hsv[2] * (1 - (hsv[1] * (1 - f)));
            switch (i) {
               case 0: rgb[0] = hsv[2]; rgb[1] = t; rgb[2] = p; break;
               case 1: rgb[0] = q; rgb[1] = hsv[2]; rgb[2] = p; break;
               case 2: rgb[0] = p; rgb[1] = hsv[2]; rgb[2] = t; break;
               case 3: rgb[0] = p; rgb[1] = q; rgb[2] = hsv[2]; break;
               case 4: rgb[0] = t; rgb[1] = p; rgb[2] = hsv[2]; break;
               default: rgb[0] = hsv[2]; rgb[1] = p; rgb[2] = q; break;
            }
            return rgb;
         }
         FCXXCoord rgbToHSV(FCXXCoord rgb) {
            FCXXCoord hsv;
            hsv[3] = rgb[3];
            float r = rgb[0], g = rgb[1], b = rgb[2];
            float mn, mx, delta;
            mn = ((r<g?r:g)<b?(r<g?r:g):b);
            mx = ((r>g?r:g)>b?(r>g?r:g):b);
            hsv[2] = mx;
            delta = mx - mn;
            if (mx != 0) hsv[1] = delta / mx;
            else { hsv[1] = 0; hsv[0] = -1; return hsv; }
            if (fabs(delta) < 1e-6) hsv[0] = 0.;
            else if (r == mx) hsv[0] = (g - b) / delta;
            else if (g == mx) hsv[0] = 2 + (b - r) / delta;
            else hsv[0] = 4 + (r - g) / delta;
            hsv[0] *= 60.;
            if (hsv[0] < 0) hsv[0] += 360.;
            return hsv;
         }

         static AtomPropertyRampColorRule *defaultRampRule() { return new AtomPropertyRampColorRule(); }

         virtual FCXXCoord colorForAtom(const gemmi::CRA &cra) override {
            float value;
            if (rampType == ResidueNumber) {
               value = (cra.residue && cra.residue->seqid.num) ? float(*cra.residue->seqid.num) : startValue;
               if (value < startValue) value = startValue;
               if (value > endValue) value = endValue;
            } else if (rampType == BFactor) {
               value = cra.atom ? cra.atom->b_iso : startValue;
               if (value < startValue) value = startValue;
               if (value > endValue) value = endValue;
            } else {
               value = 0.5;
            }
            return colorForValue(value);
         }

         FCXXCoord colorForValue(float value) {
            if (value < startValue) value = startValue;
            if (value > endValue) value = endValue;
            // CoordSpline::addPair() discards the x-coordinate and
            // coordForXEquals() interprets its argument in control-point-index
            // space (~0..nCtlPts), not data space. Passing a raw residue number
            // (range up to endValue) saturates to the last control point, so
            // the ramp collapses to its end colour. Map value (in
            // [startValue,endValue]) onto that space, inverting the index
            // formula in coordForXEquals so frac 0..1 spans spline[0..n-1].
            float frac = (endValue > startValue) ? (value - startValue) / (endValue - startValue) : 0.f;
            int nC = spline.nCtlPts();
            int nS = spline.nSplinePts();
            float x = value;
            if (nC > 0 && nS > 0)
               x = (frac * float(nS - 1) + 1.0f) * float(nC) / float(nS) - 1.0f;
            FCXXCoord hsv = spline.coordForXEquals(x);
            while (hsv[0] <   0.) hsv[0] += 360.;
            while (hsv[0] > 360.) hsv[0] -= 360.;
            hsv[1] = (hsv[1]<0.?0.:hsv[1]); hsv[1] = (hsv[1]>1.?1.:hsv[1]);
            hsv[2] = (hsv[2]<0.?0.:hsv[2]); hsv[2] = (hsv[2]>1.?1.:hsv[2]);
            hsv[3] = (hsv[3]<0.?0.:hsv[3]); hsv[3] = (hsv[3]>1.?1.:hsv[3]);
            return hsvToRGB(hsv);
         }
      };
   }
}

#endif // AtomPropertyRampColorRule_gemmi_hh
