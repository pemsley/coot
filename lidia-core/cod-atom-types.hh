/* lidia-core/cod-atom-types.hh
 * 
 * Copyright 2016 by Medical Research Council
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

#include <vector>
#include <string>

#include "rdkit-interface.hh"

#include "third-neighbour-info-t.hh"

#include "cod-atom-type-t.hh"
#include "primes.hh"

namespace cod {

   // put these in a class?
   static std::map<std::string, std::pair<unsigned int, unsigned int> > element_period_group_map;
   void fill_element_period_group_map();
   
   class atom_types_t {
      void handle_bigger_rings_from_fused_rings(RDKit::ROMol &rdkm,
						const std::vector<std::vector<int> > &fused_rings);

      bool is_ring_member(unsigned int iat_ui,
			  const std::vector<std::vector<int> > &fused_rings);

      // nb_level = 0 mean this atom neighbour info - the whole thing
      // 
      // can throw a std::runtime_error
      //
      // std::pair<std::string, std::list<third_neighbour_info_t> >
      atom_type_t
      get_cod_atom_type(const RDKit::Atom *atom_base_p,
			const RDKit::Atom *atom_nb_1_p,
			const RDKit::Atom *atom_nb_2_p,
			const RDKit::Atom *atom_nb_3_p,
			const RDKit::ROMol &rdkit_mol,
			int nb_level=0);

      std::vector<std::string> sort_neighbours(const std::vector<std::string> &neighbours_in,
					       int level);
      
      static bool neighbour_sorter(const std::string &a, const std::string &b);

      static bool fei_neighb_sorter(const std::string &a, const std::string &b);
      
   
      cod::third_neighbour_info_t
      get_cod_nb_3_type(const RDKit::Atom *atom_base_p, // the parent of atom_p
			const RDKit::Atom *atom_nb_1,
			const RDKit::Atom *atom_nb_2,
			const RDKit::Atom *atom_nb_3,
			const RDKit::ROMol &rdkit_mol);

      bool check_for_3rd_nb_info(const RDKit::Atom *atom_base_p,
				 const RDKit::Atom *atom_nb_1,
				 const RDKit::Atom *atom_nb_2,
				 const RDKit::Atom *atom_nb_3,
				 const RDKit::ROMol &rdkm);

      bool related_via_angle(const RDKit::Atom *atom_in_1_p,
			     const RDKit::Atom *atom_in_2_p,
			     const RDKit::ROMol &rdkm) const;

      // first is 3rd level (without 3rd neighbour info) and second is full (with neighbour info)
      // 
      std::pair<std::string, std::string>
      make_cod_level_3_and_4_atom_type(const RDKit::Atom *base_atom_p,
				       const std::string &atom_ele,
				       const std::vector<std::string> &neighbour_types,
				       const std::list<third_neighbour_info_t> &tniv,
				       int level);
      // which calls:
      std::string make_cod_3rd_neighb_info_type(const std::list<third_neighbour_info_t> &tniv);

      unsigned int get_smallest_ring_info(const RDKit::Atom *atom_p) const;

      // neighbour info atom type: e.g. "N(CCS)(SN)" -> "C-3:S-2:"
      // 
      std::string make_cod_level_2_atom_type(const RDKit::Atom *base_atom_p, const RDKit::ROMol &rdkm);

      //
      unsigned int make_hash_index(const RDKit::Atom *base_atom_p) const;

      unsigned int make_hash_index(const RDKit::Atom *base_atom_p, const primes &p) const;
      
      // which calls
      std::pair<unsigned int, unsigned int> get_period_group(const RDKit::Atom *at) const;
      

      std::vector<std::vector<int> > trace_path(unsigned int idx,
						const std::map<int, std::vector<int> > &bond_map,
						unsigned int n_max_bonds);
      std::vector<std::vector<int> > 
      trace_path(unsigned int idx,
		 std::vector<int> in_path_indices,
		 unsigned int target_idx,
		 const std::map<int, std::vector<int> > &bond_map,
		 unsigned int level);

      static bool atomRingSorter(const std::vector<int> &r1, const std::vector<int> &r2);


      // int hybridization_to_int(RDKit::Atom::HybridizationType) const;

   public:
      // can throw a std::runtime_error
      //
      // rdkit_mol is not const because there is no const beginAtoms() operator.
      std::vector<atom_type_t>
      get_cod_atom_types(RDKit::ROMol &rdkm, bool add_name_as_property = true);
   };

}
