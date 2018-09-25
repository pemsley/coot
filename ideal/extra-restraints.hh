/* ideal/extra-restraints.hh
 * 
 * Copyright 2013 by Medical Research Council
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

#include "coot-utils/coot-coord-utils.hh"

namespace coot {
   
   class extra_restraints_t {

      bool matches_bond_template_p(const std::vector<std::string> &words) const;
      bool matches_angle_template_p(const std::vector<std::string> &words) const;

   public:

      class extra_bond_restraint_t {

      public:
	 atom_spec_t atom_1;
	 atom_spec_t atom_2;
	 double bond_dist;
	 double esd;
	 extra_bond_restraint_t() {}
	 extra_bond_restraint_t(const atom_spec_t &a1, const atom_spec_t &a2, double d, double e) {
	    atom_1 = a1;
	    atom_2 = a2;
	    bond_dist = d;
	    esd = e;
	 }
	 bool operator==(const residue_spec_t &rs) const {
	    if (residue_spec_t(atom_1) == rs)
	       return true;
	    if (residue_spec_t(atom_2) == rs)
	       return true;
	    return false;
	 }
	 bool is_deviant(const double &real_dist, const double &n_sigma) const {
	    return (fabs(real_dist-bond_dist)/esd > n_sigma);
	 }
      };

      class bond_eraser {
      public:
	 double n_sigma_lim;
	 
	 std::map<std::pair<atom_spec_t, atom_spec_t>, double,
		  bool(*)(const std::pair<atom_spec_t, atom_spec_t> &,
			  const std::pair<atom_spec_t, atom_spec_t> &) > dist_map;
	 
	 
 	 bond_eraser(const std::map<std::pair<atom_spec_t, atom_spec_t>, double,
 		     bool(*)(const std::pair<atom_spec_t, atom_spec_t> &,
 		             const std::pair<atom_spec_t, atom_spec_t> &)
		     > &dist_map_in, double nsi) {
 	    n_sigma_lim = nsi;
 	    dist_map = dist_map_in;
 	 } 
	 bool operator()(const extra_bond_restraint_t &br) {
	    std::pair<atom_spec_t, atom_spec_t> p(br.atom_1, br.atom_2);
	    std::map<std::pair<atom_spec_t, atom_spec_t>, double>::const_iterator it =
	       dist_map.find(p);
	    if (it != dist_map.end()) {
	       bool v =  (fabs(br.bond_dist - it->second)/br.esd >= n_sigma_lim);
	       std::cout << "comparing (" << br.bond_dist << " - " << it->second << ")/" << br.esd
			 << " = " << fabs(br.bond_dist - it->second)/br.esd << " vs " << n_sigma_lim
			 << " -> " << v << std::endl;
	       return (fabs(br.bond_dist - it->second)/br.esd >= n_sigma_lim);
	    } else {
	       std::cout << "not found bond restraint between " << br.atom_1 <<  " " << br.atom_2 << std::endl;
	       return false;
	    } 
	 } 
      };
        
      class extra_angle_restraint_t {
      public:
	 atom_spec_t atom_1;
	 atom_spec_t atom_2;
	 atom_spec_t atom_3;
	 double angle;
	 double esd;
	 extra_angle_restraint_t(const atom_spec_t &a1, const atom_spec_t &a2,
				   const atom_spec_t &a3, 
				   double angle_in, double esd_in) {
	    atom_1 = a1;
	    atom_2 = a2;
	    atom_3 = a3;
	    angle = angle_in;
	    esd = esd_in;
	 }
      };

      class extra_torsion_restraint_t {
      public:
	 atom_spec_t atom_1;
	 atom_spec_t atom_2;
	 atom_spec_t atom_3;
	 atom_spec_t atom_4;
	 double torsion_angle;
	 double esd;
	 int period;
	 extra_torsion_restraint_t(const atom_spec_t &a1, const atom_spec_t &a2,
				   const atom_spec_t &a3, const atom_spec_t &a4,
				   double torsion_in, double esd_in, int period_in) {
	    atom_1 = a1;
	    atom_2 = a2;
	    atom_3 = a3;
	    atom_4 = a4;
	    torsion_angle = torsion_in;
	    esd = esd_in;
	    period = period_in;
	 }
      };
      
      class extra_start_pos_restraint_t {
      public:
	 atom_spec_t atom_1;
	 double esd;
	 extra_start_pos_restraint_t(const atom_spec_t &a1, double e) {
	    atom_1 = a1;
	    esd = e;
	 }
      };

      std::vector<extra_bond_restraint_t> bond_restraints;
      std::vector<extra_angle_restraint_t> angle_restraints;
      std::vector<extra_torsion_restraint_t> torsion_restraints;
      std::vector<extra_start_pos_restraint_t> start_pos_restraints;
      std::vector<parallel_planes_t> parallel_plane_restraints;

      void read_refmac_extra_restraints(const std::string &file_name);

      bool has_restraints() const {
	 if (bond_restraints.size() > 0)
	    return true;
	 if (angle_restraints.size() > 0)
	    return true;
	 else if (torsion_restraints.size() > 0)
	    return true;
	 else if (start_pos_restraints.size() > 0)
	    return true;
	 else if (parallel_plane_restraints.size() > 0)
	    return true;
	 else 
	    return false;
      }
      
      void clear() {
	 bond_restraints.clear();
	 angle_restraints.clear();
	 torsion_restraints.clear();
	 start_pos_restraints.clear();
         parallel_plane_restraints.clear();
      }
      void add_restraints(const extra_restraints_t &r) {
	 for (unsigned int i=0; i<r.bond_restraints.size(); i++)
	    bond_restraints.push_back(r.bond_restraints[i]);
	 for (unsigned int i=0; i<r.angle_restraints.size(); i++)
	    angle_restraints.push_back(r.angle_restraints[i]);
	 for (unsigned int i=0; i<r.torsion_restraints.size(); i++)
	    torsion_restraints.push_back(r.torsion_restraints[i]);
	 for (unsigned int i=0; i<r.start_pos_restraints.size(); i++)
	    start_pos_restraints.push_back(r.start_pos_restraints[i]);
	 for (unsigned int i=0; i<r.parallel_plane_restraints.size(); i++)
	    parallel_plane_restraints.push_back(r.parallel_plane_restraints[i]);
      }
      void delete_restraints_for_residue(const residue_spec_t &rs);
      // updates restraint on atom if it can, else adds
      void add_start_pos_restraint(const atom_spec_t &atom_1_in, double esd_in) {
	 bool already_exists = false;
	 for (unsigned int i=0; i<start_pos_restraints.size(); i++) {
	    if (start_pos_restraints[i].atom_1 == atom_1_in) {
	       start_pos_restraints[i].esd = esd_in;
	       already_exists = true;
	       break;
	    } 
	 }
	 if (! already_exists) {
	    extra_start_pos_restraint_t e(atom_1_in, esd_in);
	    start_pos_restraints.push_back(e);
	 }
      }

      // We want to interpolate proSMART restraints from start to final model.
      // We have proSMART restraints for both models.
      // 
      // Let's return a list of extra bond restraint indices that are
      // between the bond restraints of this (presumably start) and
      // final.
      // 
      std::vector<std::pair<unsigned int, unsigned int> >
      find_pair_indices(const extra_restraints_t &final) const;

      // n_path_points is the number of points along the trajectory.
      // It should be at least 2.
      // 
      void
      write_interpolated_restraints(const extra_restraints_t &final,
				    unsigned int n_path_points,
				    std::string file_name_stub) const;

      std::map<mmdb::Atom *, clipper::Coord_orth>
      position_point_map(mmdb::Manager *mol_running, mmdb::Manager *mol_ref) const;
      

      void
      write_interpolated_models_and_restraints(const extra_restraints_t &final,
					       mmdb::Manager *mol_1, // corresponds to this
					       mmdb::Manager *mol_2, // corresponds to final
					       unsigned int n_path_points,
					       std::string file_name_stub) const;

      void write_interpolated_restraints(std::ofstream &f,
					 const std::vector<extra_bond_restraint_t> &final_bonds_restraints,
					 double frac,
					 unsigned int idx_1,
					 unsigned int idx2) const;

      void write_interpolated_models(mmdb::Manager *mol_running,
				     const std::map<mmdb::Atom *, clipper::Coord_orth> &matching_atoms_1,
				     const std::map<mmdb::Atom *, clipper::Coord_orth> &matching_atoms_2,
				     unsigned int n_path_points,
				     std::string file_name_stub) const;
      
   };

}
