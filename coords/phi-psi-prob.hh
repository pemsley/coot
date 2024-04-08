/*
 * coords/phi-psi-prob.hh
 * 
 * Copyright 2022 by Medical Research Council
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

#ifndef PHI_PSI_PROB_HH
#define PHI_PSI_PROB_HH

#include "coot-utils/coot-rama.hh"
#include "coords/Cartesian.h"
#include "coords/ramachandran-container.hh"

namespace coot {
   class phi_psi_prob_t {
      public:
         phi_psi_prob_t() {}
      phi_psi_prob_t(const util::phi_psi_t &pp, const Cartesian &pos, const ramachandrans_container_t &rama_container);
      util::phi_psi_t phi_psi;
      Cartesian position;
      double probability;
      bool is_allowed_flag;
      bool is_allowed() const { return is_allowed_flag; }
      std::string residue_name() const { return phi_psi.residue_name(); }
   };
}


#endif // PHI_PSI_PROB_HH

