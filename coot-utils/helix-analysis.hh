/*
 * coot-utils/helix-analysis.hh
 *
 * Copyright 2012 by University of York
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

#include <stdexcept>
#include <vector>

#include <clipper/core/clipper_types.h>
#include <clipper/core/coords.h>

#include "mini-mol/atom-quads.hh"

namespace coot {

   // All angles internal in radians here.
   // 
   class helix_params_t {
      clipper::Mat33<double> A;
      clipper::Coord_orth B;
      void calc_A();
      void calc_B();
   public:
      helix_params_t(int resno_start_in, atom_quad quad_in, double t) {
	 resno_start = resno_start_in;
	 quad = quad_in;
	 torsion = clipper::Util::d2rad(t);
	 calc_A();
	 calc_B();
      }
      int resno_start; // start of this quad
      atom_quad quad;
      double torsion;
   }; 

   class helix_params_container_t {
      // get the set of atoms starting from the given residue serial number.
      atom_quad get_quad(const std::string &atom_name, mmdb::Chain *chain_p, int res_serial_no);
      mmdb::Manager *mol;
	 
   public:
      helix_params_container_t() {}
      std::vector<helix_params_t> params;
      void make(mmdb::Manager *mol, const std::string chain_id, int resno_start, int resno_end);
   };

}
