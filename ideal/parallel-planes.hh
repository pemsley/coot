/* ideal/parallel-planes.hh
 * 
 * Copyright 2013, 2015 by Medical Research Council
 * Author: Paul Emsley
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

#ifndef PARALLEL_PLANES_HH
#define PARALLEL_PLANES_HH

#include "coot-utils/coot-coord-utils.hh"

namespace coot {

   // Actually, the planes don't have to be parallel, there can be any
   // target angle - but typically the target angle is 0.
   // 
   
   class parallel_plane_atoms_t {
   public:
      residue_spec_t res_spec;
      std::vector<std::string> atom_names;
      std::string alt_conf;
      parallel_plane_atoms_t() {}
      parallel_plane_atoms_t(const residue_spec_t &res_spec_in,
			     const std::vector<std::string> &atom_names_in,
			     const std::string &alt_conf_in) :
         res_spec(res_spec_in), atom_names(atom_names_in), alt_conf(alt_conf_in) {}
      unsigned int size() const { return atom_names.size(); }
   };
   
   class parallel_planes_t {
      int parse_2nd_plane(const std::vector<std::string> &words, int offset);
      int parse_dist_and_type(const std::vector<std::string> &words, int offset);
   public:
      parallel_plane_atoms_t plane_1_atoms;
      parallel_plane_atoms_t plane_2_atoms;
      double target_angle;
      double sigma_angle;
      double sigma_distance;
      double sigma_combined_planes;  // because the target function tries to
                                     // flatten the atoms of translated plane pairs
                                     // needs synthetic sigma.
      std::pair<bool, double> distance;
      bool matches;
      parallel_planes_t(const residue_spec_t &res_spec_1,
			const residue_spec_t &res_spec_2,
			const std::vector<std::string> &ap1_names,
			const std::vector<std::string> &ap2_names,
			const std::string &alt_conf_1,
			const std::string &alt_conf_2) {
	 plane_1_atoms = parallel_plane_atoms_t(res_spec_1, ap1_names, alt_conf_1);
	 plane_2_atoms = parallel_plane_atoms_t(res_spec_2, ap2_names, alt_conf_2);
	 target_angle  = 0;
	 sigma_angle = 5;
	 sigma_distance = 0.2;
	 sigma_combined_planes = 0.06;
	 matches = false; // not parsed from a line (shouldn't matter)
	 distance = std::pair<bool, double> (false,0);
      }
      explicit parallel_planes_t(const std::string &line);
      void set_target_angle(double t) { target_angle = t; }
      void set_sigma_angle(double s) { sigma_angle = s; }
      void set_distance(double d) { distance = std::pair<bool,double> (true,d); }
      friend std::ostream &operator<<(std::ostream &s, parallel_planes_t pp);
   };
   std::ostream &operator<<(std::ostream &s, parallel_planes_t pp);
}

#endif // PARALLEL_PLANES_HH
