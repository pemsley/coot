/* src/c-interface-sharpen.cc
 * 
 * Copyright 2013 by Bernhard Lohkamp
 * Copyright 2016 by Medical Research Council
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifdef USE_PYTHON
#ifndef PYTHONH
#define PYTHONH
#include <Python.h>
#endif
#endif

#include "graphics-info.h"
#include "c-interface.h"

/*  ----------------------------------------------------------------------- */
/*           Map kurtosis B factor optimization                             */
/*  ----------------------------------------------------------------------- */
float optimal_B_kurtosis(int imol) {
// 
// CALCULATES THE OPTIMAL BFACTOR:
// PERFORMS A GOLDEN SECTION SEARCH ON
// THE KURTOSIS OF THE ENTIRE SUPPLIED MAP
// RICHARDTJÖRNHAMMAR 2016-05
   // TOL is for finding maxium with golden search
   // TOLB is for search window for local maximum
   // substract a linear background to find local maxima
   //
   float sharpening_limit = graphics_info_t::map_sharpening_scale_limit;
   float golden_ratio = (sqrt(5.0)-1.0)*0.5;
   float kurtosis = 0.0, B_optimal = 0.0;
   float a = -1.0*sharpening_limit, b = 1.0*sharpening_limit, TOL = 1E-2, TOLB = 40E0;
   float fc = 0.0, fd = 0.0, k = 0.0, m = 0.0;
   float c = b-golden_ratio*(b-a);
   float d = a+golden_ratio*(b-a);
   float a0 = a;
   bool what;

   if (is_valid_map_molecule(imol)) {
      if (graphics_info_t::molecules[imol].sharpen_b_factor_kurtosis_optimised() < -999.0) {

         graphics_info_t::molecules[imol].sharpen(a, false, 0);
         map_statistics_t ms01 = graphics_info_t::molecules[imol].map_statistics();
         fc = ms01.kurtosis;
         graphics_info_t::molecules[imol].sharpen(b, false, 0);
         map_statistics_t ms02 = graphics_info_t::molecules[imol].map_statistics();
         fd = ms02.kurtosis;

         k = (fd-fc)/(b-a);
         m = fc;

         while( d-c > TOL )
         {
            graphics_info_t::molecules[imol].sharpen(c, false, 0);
            map_statistics_t ms1 = graphics_info_t::molecules[imol].map_statistics();
            what = d-c>TOLB;
            if (what) {
               fc = ms1.kurtosis/(k*(c-a0)+m);
            } else {
               fc = ms1.kurtosis;
            }
            graphics_info_t::molecules[imol].sharpen(d, false, 0);
            map_statistics_t ms2 = graphics_info_t::molecules[imol].map_statistics();
            if (what) {
               fd = ms2.kurtosis/(k*(d-a0)+m);
            } else {
               fd = ms2.kurtosis;
            }

            if( fc > fd ) { // FIND MAXIMUM
               b = d; d = c;
               c = b - golden_ratio*( b - a );
            } else {
               a = c; c = d;
               d = a + golden_ratio*( b - a );
            }
         }
         B_optimal    = (c+d)*0.5;
         graphics_info_t::molecules[imol].set_sharpen_b_factor_kurtosis_optimised(B_optimal);
      } else {
         // have already a calculated one, so use that one
         B_optimal = graphics_info_t::molecules[imol].sharpen_b_factor_kurtosis_optimised();
      }
   }
   return B_optimal;
}
