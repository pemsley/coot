/*
 * coot-utils/c-beta-deviations.hh
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

#include <map>
#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>
#include "mini-mol/atom-quads.hh"

namespace coot {

   class c_beta_deviation_t {

      // there can be several c_beta deviations per residue (alt confs)

   public:
      c_beta_deviation_t(mmdb::Atom *at_in, const clipper::Coord_orth &pos_ideal_in,  double dist_in) :
	 at(at_in), pos_ideal(pos_ideal_in), dist(dist_in) { }
      c_beta_deviation_t() : at(nullptr), pos_ideal(clipper::Coord_orth(-1,-1,-1)), dist(0.0) {} // needed for map
      mmdb::Atom *at;
      clipper::Coord_orth pos_ideal;
      double dist;
   };


   std::map<mmdb::Residue *, std::map<std::string, c_beta_deviation_t> >
   get_c_beta_deviations(mmdb::Manager *mol);
   // maybe other, more fine-grained arguments needed?

   std::map<std::string, c_beta_deviation_t> get_c_beta_deviations(mmdb::Residue *residue_p);

   clipper::Coord_orth make_CB_ideal_pos(const atom_quad &q, const std::string &res_name);
   
}
