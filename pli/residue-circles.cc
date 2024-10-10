/* lbg/residue-circles.cc
 * 
 * Author: Paul Emsley
 * Copyright 2010 by The University of Oxford
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

#ifdef USE_PYTHON
#include <Python.h>
#endif

#include "flev.hh"
#include "optimise-residue-circles.hh"
#include "flev-scale-factor.hh"


std::vector<int>
flev_t::get_primary_indices() const {

   std::vector<int> primary_indices;  // this primary_indices needs to
   // get passed to the
   // primary_indices used in
   // residue cirlce optimization.

   for(unsigned int ic=0; ic<residue_circles.size(); ic++) {
      if (residue_circles[ic].is_a_primary_residue()) {
         primary_indices.push_back(ic);
      }
   }
   return primary_indices;
}

void
flev_t::initial_residues_circles_layout() {

   std::cout << "HHHHHHHHHHHHHere A with residue_circles.size() " << residue_circles.size() << std::endl;

   if (true) {
      std::cout << "------------------- flev_t::initial_residues_circles_layout() residue circles ------------" << std::endl;
      for (unsigned int i=0; i<residue_circles.size(); i++) {
         const auto &rc = residue_circles[i];
         std::cout << "   " << std::setw(2) << i << " : "
                 << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.x << " "
                 << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.y << std::endl;
      }
   }


   // when we move a primary, we want to know it's index in
   // residue_circles, because that's what we really want to move.
   //
   // std::vector<std::pair<int, lbg_info_t::residue_circle_t> > primaries;

   // now a class data member because it is used in the layout penalty
   // function (we want to have nice bond lengths for residues bonded
   // to the ligand).
   //
   std::vector<int> primary_indices;  // this primary_indices needs to
   // get passed to the
   // primary_indices used in
   // residue circle optimization.

   for(unsigned int ic=0; ic<residue_circles.size(); ic++) {
      if (residue_circles[ic].is_a_primary_residue()) {
         primary_indices.push_back(ic);
      }
   }

   // primaries get placed first and are checked for non-crossing
   // ligand interaction bonds (they have penalty scored added if they
   // cross).
   //
   try {
      std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair =
         mol.ligand_extents();
      ligand_grid grid(l_e_pair.first, l_e_pair.second);
      grid.fill(mol);

      std::cout << "HHHHHHHHHHHHHere B with residue_circles.size() " << residue_circles.size() << std::endl;
      if (true) {
         std::cout << "------------------- flev_t::initial_residues_circles_layout() post B residue circles ------------"
                   << std::endl;
         for (unsigned int i=0; i<residue_circles.size(); i++) {
            const auto &rc = residue_circles[i];
            std::cout << "   " << std::setw(2) << i << " : "
                    << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.x << " "
                    << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.y << std::endl;
         }
      }

      for (unsigned int iprimary=0; iprimary<primary_indices.size(); iprimary++) {
         int idx = primary_indices[iprimary];
         std::vector<std::pair<lig_build::pos_t, double> > attachment_points =
            residue_circles[idx].get_attachment_points(mol);
         initial_primary_residue_circles_layout(grid, idx, attachment_points);
      }
      std::cout << "HHHHHHHHHHHHHere C with residue_circles.size() " << residue_circles.size() << std::endl;
      if (true) {
         std::cout << "------------------- flev_t::initial_residues_circles_layout() post C residue circles ------------"
                   << std::endl;
         for (unsigned int i=0; i<residue_circles.size(); i++) {
            const auto &rc = residue_circles[i];
            std::cout << "   " << std::setw(2) << i << " : "
                    << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.x << " "
                    << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.y << std::endl;
         }
      }
      // show_mol_ring_centres();

      position_non_primaries(grid, primary_indices); // untrap residues as needed.

      std::cout << "HHHHHHHHHHHHHere D with residue_circles.size() " << residue_circles.size() << std::endl;
      if (true) {
         std::cout << "------------------- flev_t::initial_residues_circles_layout() post D residue circles ------------"
                   << std::endl;
         for (unsigned int i=0; i<residue_circles.size(); i++) {
            const auto &rc = residue_circles[i];
            std::cout << "   " << std::setw(2) << i << " : "
                    << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.x << " "
                    << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.y << std::endl;
         }
      }

   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
   }
}


void
flev_t::refine_residue_circle_positions() { // changes the positions of residue_circles

   std::cout << "--------------- refine_residue_circle_positions() --- start --- "
             << residue_circles.size() << std::endl;

   std::vector<int> primary_indices = get_primary_indices();
   initial_residues_circles_layout(); // twiddle residue_circles
   std::vector<residue_circle_t> current_circles = residue_circles;

   for (int iround=0; iround<30; iround++) {
      std::cout << "flev_t::refine_residue_circle_positions(): iround      " << iround << std::endl;
      std::pair<int, std::vector<residue_circle_t> > new_c =
	 optimise_residue_circle_positions(residue_circles, current_circles, primary_indices);
      current_circles = new_c.second;
      if (new_c.first == GSL_ENOPROG)
	 break;
      if (new_c.first == GSL_SUCCESS) {
	 break;
      }
   }
   residue_circles = current_circles;
}

std::pair<bool,lig_build::pos_t>
flev_t::get_residue_circles_top_left() const {

   bool status = false;
   lig_build::pos_t p(999999,999999);

   if (residue_circles.size()) {
      status = true;
      for (unsigned int i=0; i<residue_circles.size(); i++) {
	 if (residue_circles[i].pos.x < p.x) p.x = residue_circles[i].pos.x;
	 if (residue_circles[i].pos.y < p.y) p.y = residue_circles[i].pos.y;
      }
   }
   return std::pair<bool, lig_build::pos_t> (status, p);
}




// Return the status and updated positions.
//
std::pair<int, std::vector<residue_circle_t> >
pli::optimise_residue_circles::solution() const {
   return std::pair<int, std::vector<residue_circle_t> > (status, current_circles);
}


void
pli::optimise_residue_circles::setup_angles() {

   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      std::vector<int> residue_indexes;
      for (unsigned int ires=0; ires<current_circles.size(); ires++) {
	 for (unsigned int ibond=0;
	      ibond<current_circles[ires].bonds_to_ligand.size();
	      ibond++) {
	    if (current_circles[ires].bonds_to_ligand[ibond].ligand_atom_name == mol.atoms[iat].get_atom_name()) {
	       residue_indexes.push_back(ires);
	    }
	 }
      }
      if (residue_indexes.size() > 1) {
	 // a bit complicated.
	 angle a1(iat, residue_indexes[0], residue_indexes[1]);
	 angles.push_back(a1);
	 if (residue_indexes.size() > 2) {
	    angle a2(iat, residue_indexes[1], residue_indexes[2]);
	    angles.push_back(a2);
	    angle a3(iat, residue_indexes[0], residue_indexes[2]);
	    angles.push_back(a3);
	 }
      }
   }

}


