/* coot-utils/atom-overlaps.hh
 * 
 * Copyright 2015 by Medical Research Council
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


#ifndef ATOM_OVERLAPS_HH
#define ATOM_OVERLAPS_HH

// #include <mmdb2/mmdb_manager.h>
// #include <vector>

#include "compat/coot-sysdep.h"
#include "geometry/protein-geometry.hh"
#include <unordered_map>

namespace coot {

   class atom_overlaps_dots_container_t {
   public:
      class spikes_t {
      public:
	 std::string type;
	 std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > positions;
	 const std::pair<clipper::Coord_orth, clipper::Coord_orth> &operator[](unsigned int idx) const {
	    return positions[idx];
	 }
	 unsigned int size() const { return positions.size(); }
      };
      class dot_t {
      public:
	 double overlap;
	 clipper::Coord_orth pos;
	 std::string col;
	 dot_t(double o, const std::string &col_in, const clipper::Coord_orth &pos_in) : pos(pos_in), col(col_in) {
	    overlap = o;
	 }
      };
      atom_overlaps_dots_container_t() {
	 // I think this speeds things up a bit.
	 dots["big-overlap"  ].reserve(2500);
	 dots["small-overlap"].reserve(2500);
	 dots["close-contact"].reserve(2500);
	 dots["wide-contact" ].reserve(2500);
	 dots["H-bond"       ].reserve(2500);
	 dots["vdw-surface"  ].reserve(2500);
      }

      // 6,000 a/t -> size 145,000
      // 1,000 a/t -> size  27,000
      explicit atom_overlaps_dots_container_t(unsigned int n_atoms_per_thread) {

	 dots["big-overlap"  ].reserve(25 * n_atoms_per_thread);
	 dots["small-overlap"].reserve(25 * n_atoms_per_thread);
	 dots["close-contact"].reserve(25 * n_atoms_per_thread);
	 dots["wide-contact" ].reserve(25 * n_atoms_per_thread);
	 dots["H-bond"       ].reserve(25 * n_atoms_per_thread);
	 dots["vdw-surface"  ].reserve(25 * n_atoms_per_thread);
      }

      std::unordered_map<std::string, std::vector<dot_t> > dots;
      spikes_t clashes;
      void add(const atom_overlaps_dots_container_t &other) {
	 std::unordered_map<std::string, std::vector<dot_t> >::const_iterator it;
	 for (it=other.dots.begin(); it!=other.dots.end(); ++it)
	    if (it->second.size())
	       dots[it->first].insert(dots[it->first].end(),it->second.begin(), it->second.end());
	 if (other.clashes.positions.size())
	    clashes.positions.insert(clashes.positions.end(),
				     other.clashes.positions.begin(),
				     other.clashes.positions.end());
      }
      double score() const {
	 std::unordered_map<std::string, std::vector<dot_t> >::const_iterator it;
	 // do these match the types in overlap_delta_to_contact_type()?
	 double r = 0;
	 it = dots.find("H-bond");
	 if (it != dots.end()) r += it->second.size();
	 it = dots.find("wide-contact");
	 if (it != dots.end()) r += 0.1 * it->second.size();
	 it = dots.find("close-contact");
	 if (it != dots.end()) r -= 0.0 * it->second.size();
	 it = dots.find("small-overlap");
	 if (it != dots.end()) r -= 0.1 * it->second.size();
	 it = dots.find("big-overlap");
	 if (it != dots.end()) r -= 0.6 * it->second.size();
	 r -= clashes.size();
	 return r;
      }
      void debug() const {
	 std::unordered_map<std::string, std::vector<atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
	 for (it=dots.begin(); it!=dots.end(); ++it)
	    std::cout << " contact dot map " << it->first << " size " << it->second.size() << std::endl;
      }

   };

   class atom_overlap_t {
   public:
      atom_overlap_t(mmdb::Atom *a1, mmdb::Atom *a2) {
	 atom_1 = a1;
	 atom_2 = a2;
	 overlap_volume = -1;
	 r_1 = -1;
	 r_2 = -1;
         is_h_bond = false;
         ligand_atom_index = -1;
      }
      atom_overlap_t(int ligand_atom_index_in,
		     mmdb::Atom *a1, mmdb::Atom *a2, const double &r_1_in, const double &r_2_in,
		     const double &o) {
	 ligand_atom_index = ligand_atom_index_in;
	 atom_1 = a1;
	 atom_2 = a2;
	 r_1 = r_1_in;
	 r_2 = r_2_in;
         is_h_bond = false;
	 overlap_volume = o;
      }
      int ligand_atom_index;
      double r_1, r_2;
      mmdb::Atom *atom_1; 
      mmdb::Atom *atom_2; 
      double overlap_volume;
      bool is_h_bond;
   };

   class atom_overlaps_container_t {

      void init();
      void init_for_all_atom();
      enum overlap_mode_t { CENTRAL_RESIDUE, ALL_ATOM };
      overlap_mode_t overlap_mode;
      mmdb::Manager *mol;
      bool have_dictionary; // for central residue (or should it be all residues?)
      bool molecule_has_hydrogens;
      mmdb::Residue *res_central;
      std::vector<mmdb::Residue *> neighbours;
      int udd_h_bond_type_handle;
      int udd_residue_index_handle;
      double probe_radius;
      bool ignore_water_contacts_flag;

      // for energy types -> vdw radius and h-bond type
      std::map<std::string, double> type_to_vdw_radius_map;
      std::map<mmdb::Atom *, double> central_residue_atoms_vdw_radius_map; // ligand atoms
      std::map<mmdb::Atom *, double> neighbour_atoms_vdw_radius_map; // neighbouring atoms
      dictionary_residue_restraints_t central_residue_dictionary;
      // for ligand and environment neighbours
      std::vector<dictionary_residue_restraints_t> neighb_dictionaries;
      // for all atom
      std::map<std::string, dictionary_residue_restraints_t> dictionary_map;

      double get_vdw_radius_ligand_atom(mmdb::Atom *at);
      double get_vdw_radius_neighb_atom(mmdb::Atom *at, unsigned int idx_neighb_res);
      double get_vdw_radius_neighb_atom(int idx_neighb_atom) const;
      double get_overlap_volume(const double &dist, const double &r_1, const double &r_2) const; // in A^3
      const protein_geometry *geom_p;

      bool clashable_alt_confs(mmdb::Atom *at_1, mmdb::Atom *at_2) const;
      std::vector<double> env_residue_radii;
      void setup_env_residue_atoms_radii(int i_sel_hnd_env_atoms); // fill above
      // which calls:
      double type_energy_to_radius(const std::string &te) const;

      void add_residue_neighbour_index_to_neighbour_atoms();
      std::vector<double> neighb_atom_radius;

      // first is yes/no, second is if the H is on the ligand.
      // also allow water to be return true values.
      static
      std::pair<bool, bool> is_h_bond_H_and_acceptor(mmdb::Atom *ligand_atom,
						     mmdb::Atom *env_atom,
						     int udd_h_bond_type_handle);

      // more general/useful version of the above
      class h_bond_info_t {
      public:
	 bool is_h_bond_H_and_acceptor;
	 bool is_h_bond_donor_and_acceptor;
	 bool H_is_first_atom_flag;
	 bool H_is_second_atom_flag;
	 bool donor_is_second_atom_flag;
	 h_bond_info_t(mmdb::Atom *ligand_atom,
		       mmdb::Atom *env_atom,
		       int udd_h_bond_type_handle);
         std::string format() const;
      };

      hb_t get_h_bond_type(mmdb::Atom *at);
      // store the results of a contact search.
      //
      // when we draw the surface of the ligand, for each atom, we don't want to draw
      // surface points that are nearer to another atom than the "central" atom.  So keep
      // a quick store of what's close to what and the radius.
      //
      std::map<int, std::vector<std::pair<mmdb::Atom *, double> > > ligand_atom_neighbour_map;
      std::map<int, std::vector<int> > ligand_to_env_atom_neighbour_map;
      // adds radii too.
      void fill_ligand_atom_neighbour_map();
      // include inner cusps (ugly/simple)
      bool is_inside_another_ligand_atom(int idx,
					 const clipper::Coord_orth &pt_idx_at) const;
      // exclude inner cusps (modern/pretty)
      bool is_inside_another_ligand_atom(int idx,
					 const clipper::Coord_orth &probe_pos,
					 const clipper::Coord_orth &pt_idx_at) const;

      // for all-atom contacts
      static
      bool is_inside_another_atom_to_which_its_bonded(int atom_idx,
						      mmdb::Atom *at,
						      const clipper::Coord_orth &pt_on_surface,
						      const std::vector<int> &bonded_neighb_indices,
						      mmdb::Atom **atom_selection,
						      const std::vector<double> &neighb_atom_radius);
      bool is_inside_an_env_atom_to_which_its_bonded(int idx,
						     const std::vector<int> &bonded_neighb_indices,
						     mmdb::Atom **env_residue_atoms,
						     const clipper::Coord_orth &pt_at_surface);
      double clash_spike_length;
      void mark_donors_and_acceptors();
      void mark_donors_and_acceptors_central_residue(int udd_h_bond_type_handle);
      void mark_donors_and_acceptors_for_neighbours(int udd_h_bond_type_handle);
      // return a contact-type and a colour
      static std::pair<std::string, std::string> overlap_delta_to_contact_type(double delta, bool is_h_bond);
      static std::pair<std::string, std::string> overlap_delta_to_contact_type(double delta, const h_bond_info_t &hbi,
                                                                               bool molecule_has_hydrogens_flag);
      static void test_get_type(double delta, bool is_h_bond, std::string *c_type_p, std::string *col);
      // can throw std::exception
      const dictionary_residue_restraints_t &get_dictionary(mmdb::Residue *r, unsigned int idx) const;
      // where BONDED here means bonded/1-3-angle/ring related
      enum atom_interaction_type { CLASHABLE, BONDED, IGNORED };
      atom_interaction_type
      bonded_angle_or_ring_related(mmdb::Manager *mol,
				   mmdb::Atom *at_1,
				   mmdb::Atom *at_2,
				   bool exclude_mainchain_also,
				   std::map<std::string, std::vector<std::pair<std::string, std::string> > > *bonded_neighbours,
				   std::map<std::string, std::vector<std::vector<std::string> > > *ring_list_map);
      bool are_bonded_residues(mmdb::Residue *res_1, mmdb::Residue *res_2) const;
      bool in_same_ring(mmdb::Atom *at_1, mmdb::Atom *at_2,
			std::map<std::string, std::vector<std::vector<std::string> > > &ring_list_map) const;
      // check LINK records
      bool is_linked(mmdb::Atom *at_1, mmdb::Atom *at_2) const;
      bool is_angle_related_via_link(mmdb::Atom *at_1, mmdb::Atom *at_2,
                                     const std::vector<std::pair<std::string, std::string> > &bonds_for_at_1,
                                     const std::vector<std::pair<std::string, std::string> > &bonds_for_at_2) const;
      bool is_ss_bonded_or_CYS_CYS_SGs(mmdb::Atom *at_1, mmdb::Atom *at_2) const;
      bool is_ss_bonded(mmdb::Residue *residue_p) const;
//       bool in_same_ring(const std::string &atom_name_1,
// 			const std::string &atom_name_2,
// 			const std::vector<std::vector<std::string> > &ring_list) const;
      std::vector<std::vector<std::string> > phe_ring_list() const;
      std::vector<std::vector<std::string> > his_ring_list() const;
      std::vector<std::vector<std::string> > trp_ring_list() const;
      std::vector<std::vector<std::string> > pro_ring_list() const;

      atom_overlaps_dots_container_t all_atom_contact_dots_internal_multi_thread(double dot_density_in,
										 mmdb::Manager *mol,
										 int i_sel_hnd_1,
										 int i_sel_hnd_2,
										 mmdb::realtype min_dist,
										 mmdb::realtype max_dist,
										 bool make_vdw_surface);

      atom_overlaps_dots_container_t
      all_atom_contact_dots_internal_single_thread(double dot_density_in,
						   mmdb::Manager *mol,
						   int i_sel_hnd_1,
						   int i_sel_hnd_2,
						   mmdb::realtype min_dist,
						   mmdb::realtype max_dist,
						   bool make_vdw_surface);

      static bool overlap_sorter(const atom_overlap_t &ao1, const atom_overlap_t &ao2);
      void sort_overlaps();
      bool kludge_filter(mmdb::Atom *at_1, mmdb::Atom *at_2) const;
      
   public:
      // we need mol to use UDDs to mark the HB donors and acceptors (using coot-h-bonds.hh)
      atom_overlaps_container_t(mmdb::Residue *res_central_in,
				const std::vector<mmdb::Residue *> &neighbours_in,
				mmdb::Manager *mol,
				const protein_geometry *geom_p_in);
      atom_overlaps_container_t(mmdb::Residue *res_central_in,
				mmdb::Residue *neighbour,
				mmdb::Manager *mol,
				const protein_geometry *geom_p_in);
      // this one for contact dots (around central ligand)
      atom_overlaps_container_t(mmdb::Residue *res_central_in,
				const std::vector<mmdb::Residue *> &neighbours_in,
				mmdb::Manager *mol,
				const protein_geometry *geom_p_in,
				double clash_spike_length_in,
				double probe_radius_in = 0.25);
      // all atom contact dots and atom overlaps
      atom_overlaps_container_t(mmdb::Manager *mol_in,
				const protein_geometry *geom_p_in,
				bool ignore_water_contacts_flag,
				double clash_spike_length_in = 0.5,
				double probe_radius_in = 0.25);

      // If there are no overlaps, is it because there was no dictionary for one or more residues?
      bool get_have_dictionary() const { return have_dictionary; }
      std::vector<atom_overlap_t> overlaps;
      void make_overlaps();
      void make_all_atom_overlaps();
      void contact_dots_for_overlaps() const; // old
      atom_overlaps_dots_container_t contact_dots_for_ligand(double dot_density_in = 1.02);
      // this should be a vector or derived symmetry_atom class really.
      std::vector<atom_overlap_t> symmetry_contacts(float d);
      atom_overlaps_dots_container_t all_atom_contact_dots(double dot_density = 0.5,
							   bool make_vdw_surface = false);

      // public because thread -> static
      //
      static
      atom_overlaps_dots_container_t
      contacts_for_atom(int iat,
			mmdb::Atom **atom_selection,
			const std::map<int, std::vector<int> > &contact_map,
			const std::map<int, std::vector<int> > &bonded_map,
			const std::vector<double> &neighb_atom_radius,
			int udd_h_bond_type_handle,
			bool molecule_has_hydrogens,
			double probe_radius,
			double dot_density_in,
			double clash_spike_length,
			bool make_vdw_surface);
      static
      void
      contacts_for_atoms(int iat_start, int iat_end,
			 mmdb::Atom **atom_selection,
			 const std::map<int, std::vector<int> > &contact_map,
			 const std::map<int, std::vector<int> > &bonded_map,
			 const std::vector<double> &neighb_atom_radius,
			 int udd_h_bond_type_handle,
			 bool molecule_has_hydrogens,
			 double probe_radius,
			 double dot_density_in,
			 double clash_spike_length,
			 bool make_vdw_surface,
			 atom_overlaps_dots_container_t *ao_results); // fill this

      float score(); // not const because calls all_atom_contact_dots()

      static void contacts_for_atom_test();
      static void contacts_for_atom_test_1(int iat);
      static void contacts_for_atom_test_2(int iat, mmdb::Atom **atom_selection);
      static void contacts_for_atom_test_3(int iat,
					   mmdb::Atom **atom_selection,
					   const std::map<int, std::vector<int> > &contact_map);
      static void contacts_for_atom_test_4(int iat,
					   mmdb::Atom **atom_selection,
					   const std::map<int, std::vector<int> > &contact_map,
					   const std::map<int, std::vector<int> > &bonded_map);
      static void contacts_for_atom_test_6(int iat,
					   mmdb::Atom **atom_selection,
					   const std::map<int, std::vector<int> > &contact_map,
					   const std::map<int, std::vector<int> > &bonded_map,
					   double probe_radius, double dot_density_in);
      static void contacts_for_atom_test_7(int iat,
					   mmdb::Atom **atom_selection,
					   const std::map<int, std::vector<int> > &contact_map,
					   const std::map<int, std::vector<int> > &bonded_map,
					   const std::vector<double> &neighb_atom_radius,
					   double probe_radius, double dot_density_in);
      static void contacts_for_atom_test_all(int iat,
					     mmdb::Atom **atom_selection,
					     const std::map<int, std::vector<int> > &contact_map,
					     const std::map<int, std::vector<int> > &bonded_map,
					     const std::vector<double> &neighb_atom_radius,
					     int udd_h_bond_type_handle,
					     double probe_radius,
					     double dot_density_in,
					     double clash_spike_length,
					     bool make_vdw_surface);
   };

}



#endif // ATOM_OVERLAPS_HH
