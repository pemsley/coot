/* src/rotamer.hh
 * 
 * Copyright 2001, 2002, 2003, 2004, 2006 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#ifndef ROTAMER_HH
#define ROTAMER_HH

#include "chi-angles.hh"

namespace coot {


   class rotamer : public chi_angles {

      std::pair<short int, double> probability_of_this_rotamer(const std::vector<double> &chi_angles,
							       const std::vector<coot::simple_rotamer> &rots) const;
      std::vector<std::vector<std::string> >
      coot::rotamer::rotamer_atoms(const std::string &residue_name) const;
      std::vector<std::vector<int> > rotamer_atom_names_to_indices(const std::vector<std::vector<std::string> > &residue_rotamer_atoms, PCAtom *residue_atoms, int n_residue_atoms) const;
      double chi_torsion(const std::vector<int> &chi_angle_atom_indices,
			 PCAtom *residue_atoms);
      std::vector<coot::simple_rotamer>
      get_all_rotamers(const std::string &res_type) const;

      short int similar_rotamer_chi(double target, double model) const {
	 short int is = 0;
	 double diff = target - model;
	 while (diff > 180.0)
	    diff -= 360.0;
	 while (diff < -180.0)
	    diff += 360.0;

	 if (fabs(diff) < 40.0)
	    is = 1;

	 return is;
      }
      
   public:

      rotamer(CResidue *res, short int add_extra_PHE_and_TYR_rotamers_flag ) :
	 chi_angles(res, add_extra_PHE_and_TYR_rotamers_flag) {
      }

      // For use with Z-score (which is analysis only: we don't move anything)
      // 
      rotamer(CResidue *residue) : chi_angles(residue, 0) {}

      // LEU, VAL, THR have "nomenclature" (or real) chiral centres -
      // they are not dealt with here.
      // 
      // We deal with bifurcated symmetric non-chiral side chains (PHE, ASP,
      // GLU, THR)
      // 
      int optimize_rotamer_by_atom_names();

      std::pair<short int, double> probability_of_this_rotamer(); // can't const - mmdb
                                                                  // CResidue issues...


   };
}

#endif // ROTAMER_HH
