/* src/molecule-class-info.h
 * 
 * Copyright 2008 The University of Oxford
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


// Having to set up the include files like this so that
// molecule-class-info.h can be parsed, is silly.

#include <string>
#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "Cartesian.h"
#include "mmdb-crystal.h"
#include "molecule-class-info.h"

void
molecule_class_info_t::sharpen(float b_factor) {

   int n_data = 0;
   int n_tweaked = 0;

//    clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis = original_fphis;
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis;
   clipper::HKL_info::HKL_reference_index hri;   
   for (hri = fphis.first(); !hri.last(); hri.next()) {
      n_data++;

      std::cout << " " << hri.invresolsq() << std::endl;

      fphis[hri].f() = 0.0;
      n_tweaked++;
   }
}

