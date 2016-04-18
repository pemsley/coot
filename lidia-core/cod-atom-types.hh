
#include <vector>
#include <string>

#include "rdkit-interface.hh"

#include "third-neighbour-info-t.hh"

namespace cod {

   class atom_types_t {
      void handle_bigger_rings_from_fused_rings(RDKit::ROMol &rdkm,
						const std::vector<std::vector<int> > &fused_rings);

      bool is_ring_member(unsigned int iat_ui,
			  const std::vector<std::vector<int> > &fused_rings);

      // nb_level = 0 mean this atom neighbour info - the whole thing
      // 
      // can throw a std::runtime_error
      //
      std::pair<std::string, std::list<third_neighbour_info_t> >
      get_cod_atom_type(RDKit::Atom *atom_base_p,
			RDKit::Atom *atom_nb_1_p,
			RDKit::Atom *atom_nb_2_p,
			RDKit::Atom *atom_nb_3_p,
			const RDKit::ROMol &rdkit_mol,
			int nb_level=0);

      std::vector<std::string> sort_neighbours(const std::vector<std::string> &neighbours_in,
					       int level);
      
      static bool neighbour_sorter(const std::string &a, const std::string &b);

      static bool fei_neighb_sorter(const std::string &a, const std::string &b);
      
   
      cod::third_neighbour_info_t
      get_cod_nb_3_type(RDKit::Atom *atom_base_p, // the parent of atom_p
			RDKit::Atom *atom_nb_1,
			RDKit::Atom *atom_nb_2,
			RDKit::Atom *atom_nb_3,
			const RDKit::ROMol &rdkit_mol);

      bool check_for_3rd_nb_info(RDKit::Atom *atom_base_p,
				 RDKit::Atom *atom_nb_1,
				 RDKit::Atom *atom_nb_2,
				 RDKit::Atom *atom_nb_3,
				 const RDKit::ROMol &rdkm);

      bool related_via_angle(RDKit::Atom *atom_in_1_p,
			     RDKit::Atom *atom_in_2_p,
			     const RDKit::ROMol &rdkm) const;

      std::string make_cod_type(RDKit::Atom *base_atom_p,
				const std::string &atom_ele,
				const std::vector<std::string> &neighbour_types,
				const std::list<third_neighbour_info_t> &tniv,
				int level);
      // which calls:
      std::string make_cod_3rd_neighb_info_type(const std::list<third_neighbour_info_t> &tniv);

      std::vector<std::vector<int> > trace_path(unsigned int idx,
						const std::map<int, std::vector<int> > &bond_map,
						unsigned int n_max_bonds);
      std::vector<std::vector<int> > 
      trace_path(unsigned int idx,
		 std::vector<int> in_path_indices,
		 unsigned int target_idx,
		 const std::map<int, std::vector<int> > &bond_map,
		 unsigned int level);

   public:
      // can throw a std::runtime_error
      //
      // rdkit_mol is not const because there is no const beginAtoms() operator.
      std::vector<std::string>
      get_cod_atom_types(RDKit::ROMol &rdkm, bool add_name_as_property = true);
   };

}
