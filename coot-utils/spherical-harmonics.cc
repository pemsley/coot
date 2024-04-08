/*
 * coot-utils/spherical-harmonics.cc
 *
 * Copyright 2018 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <iostream>
// #include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_legendre.h>
#include "spherical-harmonics.hh"


void
coot::spherical_harmonics(clipper::NXmap<float> &nxmap) {

   /*  -- GSL 2 or so?
   const gsl_sf_legendre_t norm = GSL_SF_LEGENDRE_SCHMIDT;
   const size_t lmax = 3;

   double x = 0.5;

   size_t dim = gsl_sf_legendre_array_n(lmax);
   double *p = malloc(sizeof(double) * dim);

   gsl_sf_legendre_array(norm, lmax, x, p);
   */


   gsl_sf_result r;
   int i = gsl_sf_legendre_P1_e(-0.5, &r);
   std::cout << "gsl_sf " << i << " " << r.val << " " << r.err << std::endl;

}

