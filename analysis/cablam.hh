/*
 * analysis/cablam.hh
 *
 * Copyright 2009 by University of Oxford
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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


#include "mini-mol/atom-quads.hh"
#include <clipper/core/coords.h>


namespace coot {

   // Given a residue selection: 
   // 
   // Find pairs (like Ramachandran)
   // Return diagnostics:
   //   a vector of cablam torsion angles (and scores when we
   //   have a table to look them up).

   class cablam_pseudo_torsion_info {
   public:
      torsion_atom_quad quad;
   };


   class cablam {
      clipper::Coord_orth get_closest_CA_CA_approach(const coot::torsion_atom_quad &quad) const;
      clipper::Coord_orth get_closest_CA_CA_approach(const clipper::Coord_orth &CA_pos_p,
						     const clipper::Coord_orth &CA_pos_t,
						     const clipper::Coord_orth &O_pos_p) const;
   public:
      cablam(mmdb::PResidue *residues, int n_sel_residues);
   };

}
