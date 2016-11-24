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

namespace coot {

   class atom_overlaps_dots_container_t {
   public:
      class spikes_t {
      public:
	 std::string type;
	 std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > positions;
	 const std::pair<clipper::Coord_orth, clipper::Coord_orth> &operator[](unsigned int idx) {
	    return positions[idx];}
	 unsigned int size() const { return positions.size(); }
      };
      atom_overlaps_dots_container_t() {}
      std::map<std::string, std::vector<clipper::Coord_orth> > dots;
      spikes_t spikes;
   };

   class atom_overlap_t {
   public:
      atom_overlap_t(mmdb::Atom *a1, mmdb::Atom *a2) {
	 atom_1 = a1;
	 atom_2 = a2;
	 overlap_volume = -1;
	 r_1 = -1;
	 r_2 = -1;
      }
      atom_overlap_t(int ligand_atom_index_in,
		     mmdb::Atom *a1, mmdb::Atom *a2, const double &r_1_in, const double &r_2_in,
		     const double &o) {
	 ligand_atom_index = ligand_atom_index_in;
	 atom_1 = a1;
	 atom_2 = a2;
	 r_1 = r_1_in;
	 r_2 = r_2_in;
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
      mmdb::Manager *mol;
      bool have_dictionary; // for central residue (or should it be all residues?)
      mmdb::Residue *res_central;
      std::vector<mmdb::Residue *> neighbours;
      int udd_h_bond_type_handle;
      int udd_residue_index_handle;
      double probe_radius;
      
      // for energy types -> vdw radius and h-bond type
      std::map<std::string, double> type_to_vdw_radius_map;
      std::map<mmdb::Atom *, double> central_residue_atoms_vdw_radius_map; // ligand atoms
      std::map<mmdb::Atom *, double> neighbour_atoms_vdw_radius_map; // neighbouring atoms
      dictionary_residue_restraints_t central_residue_dictionary;
      std::vector<dictionary_residue_restraints_t> neighb_dictionaries;
      double get_vdw_radius_ligand_atom(mmdb::Atom *at);
      double get_vdw_radius_neighb_atom(mmdb::Atom *at, unsigned int idx_neighb_res);
      double get_vdw_radius_neighb_atom(int idx_neighb_atom) const;
      double get_overlap_volume(const double &dist, const double &r_1, const double &r_2) const; // in A^3
      const protein_geometry *geom_p;

      std::vector<double> env_residue_radii;
      void setup_env_residue_atoms_radii(int i_sel_hnd_env_atoms); // fill above
      void add_residue_neighbour_index_to_neighbour_atoms();
      std::vector<double> neighb_atom_radius;

      // first is yes/no, second is if the H is on the ligand
      // 
      std::pair<bool, bool> is_h_bond_H_and_acceptor(mmdb::Atom *ligand_atom,
						     mmdb::Atom *env_atom) const;

      hb_t get_h_bond_type(mmdb::Atom *at);
      // store the results of a contact search.
      //
      // when we draw the surface of the ligand, for each atom, we don't want to draw
      // surface points that are nearer to another atom than the "central" atom.  So keep
      // a quick store of what's close to what and the radius.
      //
      std::map<int, std::vector<std::pair<mmdb::Atom *, double> > > ligand_atom_neighbour_map;
      // std::map<int, std::vector<mmdb::Atom *> > ligand_to_env_atom_neighbour_map;
      std::map<int, std::vector<int> > ligand_to_env_atom_neighbour_map;
      void fill_ligand_atom_neighbour_map();
      // include inner cusps (ugly/simple)
      bool is_inside_another_ligand_atom(int idx,
					 const clipper::Coord_orth &pt_idx_at) const;
      // exclude inner cusps (modern/pretty)
      bool is_inside_another_ligand_atom(int idx,
					 const clipper::Coord_orth &probe_pos,
					 const clipper::Coord_orth &pt_idx_at) const;
      double clash_spike_length;
      void mark_donors_and_acceptors();
      std::string overlap_delta_to_contact_type(double delta, bool is_h_bond) const;
      

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
      // this one for contact dots
      atom_overlaps_container_t(mmdb::Residue *res_central_in,
				const std::vector<mmdb::Residue *> &neighbours_in,
				mmdb::Manager *mol,
				const protein_geometry *geom_p_in,
				double clash_spike_length_in,
				double probe_radius_in = 0.25);

      std::vector<atom_overlap_t> overlaps;
      void make_overlaps();
      void contact_dots_for_overlaps() const; // old
      atom_overlaps_dots_container_t contact_dots();
   };

}



#endif // ATOM_OVERLAPS_HH
