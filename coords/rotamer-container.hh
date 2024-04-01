/*
 * coords/rotamer-container.hh
 * 
 * Copyright 2017 by Medical Research Council
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
 */


#ifndef ROTAMER_CONTAINER_HH
#define ROTAMER_CONTAINER_HH

#include <clipper/core/coords.h>
#include "utils/colour-holder.hh"

#include "ligand/rotamer.hh"

class rotamer_markup_container_t {

public:
   coot::residue_spec_t spec;
   clipper::Coord_orth pos;
   coot::colour_holder col;
   coot::rotamer_probability_info_t rpi;

   rotamer_markup_container_t(const coot::residue_spec_t &spec_in,
                              const clipper::Coord_orth &pos_in,
			      const coot::colour_holder &col_in,
                              const coot::rotamer_probability_info_t &prob_in) :
      spec(spec_in), pos(pos_in), col(col_in), rpi(prob_in) {}
   rotamer_markup_container_t() {}
};

#endif // ROTAMER_CONTAINER_HH

