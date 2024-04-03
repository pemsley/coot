/*
 * ligand/ca-torsion-info.hh
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

#ifndef CA_TORSION_INFO_HH
#define CA_TORSION_INFO_HH

#include "angles/AngleInfo.h"
#include "mini-mol/mini-mol.hh"

namespace coot {

   class CA_torsion_info_t {

      // this class works in radians

   public:
      CA_torsion_info_t(const AngleInfo &ai_in) : ai(ai_in) {
	 status = false;
      }
      const AngleInfo &ai;
      bool status;
      double angle;
      double torsion;
      // seqnum is the residue number of the first residue to be added (N)
      // if offset = 1, building forwards on the C-terminus:
      // Let's score residues N-2, N-1, N, N+1.
      // We could also score N-3, N-2, N-1, N - which is perhaps more interesting?
      std::pair<bool, double> score_fragment(const minimol::fragment &frag,
					     mmdb::Residue *res_p,
					     mmdb::Residue *res_prev_p,
					     int seqnum, int offset);
   };

}



#endif // CA_TORSION_INFO_HH
