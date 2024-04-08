/*
 * coot-utils/coot_shiftfield.h
 *
 * Copyright 2018 by Kevin Cowtan & The University of York
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


#ifndef COOT_SHIFTFIELD
#define COOT_SHIFTFIELD


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


namespace coot {

   // update the temperature-factors of the atoms in mol
   void
   shift_field_b_factor_refinement(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &fo0,
                                   const clipper::HKL_data<clipper::data32::Flag> &free,
                                   mmdb::Manager *mol, int ncycles);

   // update the coordinates of the atoms in mol
   void shift_field_xyz_refinement(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &fo0,
                                   const clipper::HKL_data<clipper::data32::Flag> &free,
                                   mmdb::Manager *mol,
                                   float resolution);

}

#endif
