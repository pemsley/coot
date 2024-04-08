/*
 * coords/phi-psi-prob.cc
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


#include "phi-psi-prob.hh"

coot::phi_psi_prob_t::phi_psi_prob_t(const coot::util::phi_psi_t &pp, const coot::Cartesian &pos, const ramachandrans_container_t &rc) {

   is_allowed_flag = true;
   phi_psi = pp;
   position = pos;
   const clipper::Ramachandran *rama = &rc.rama;

   if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
   if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

   // Use the right version of clipper!

   // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
   // if (phi_psi.is_pre_pro())
   // if (phi_psi.residue_name() != "GLY")
   // rama = &rc.rama_pre_pro;

   probability = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                   clipper::Util::d2rad(phi_psi.psi()));

   if (probability < 0.002) // from src/rama-plot.hh
      is_allowed_flag = false;

}
