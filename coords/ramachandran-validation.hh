/*
 * coords/ramachandran-container.hh
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

#ifndef RAMACHANDRAN_VALIDATION_HH
#define RAMACHANDRAN_VALIDATION_HH

#include <vector>

#include "compat/coot-sysdep.h"
#include "geometry/residue-and-atom-specs.hh"
#include "rama-plot-phi-psi.hh"
#include "phi-psi-prob.hh"

namespace coot {

   std::vector<phi_psi_prob_t> ramachandran_validation(mmdb::Manager *mol, const ramachandrans_container_t &rc);

}


#endif // RAMACHANDRAN_VALIDATION_HH
