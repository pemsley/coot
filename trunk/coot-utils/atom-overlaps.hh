
#ifndef ATOM_OVERLAPS_HH
#define ATOM_OVERLAPS_HH

// #include <mmdb2/mmdb_manager.h>
// #include <vector>

#include "compat/coot-sysdep.h"
#include "geometry/protein-geometry.hh"

namespace coot {

   class atom_overlap_t {
   public:
      atom_overlap_t(mmdb::Atom *a1, mmdb::Atom *a2) {
	 atom_1 = a1;
	 atom_2 = a2;
      }
      atom_overlap_t(mmdb::Atom *a1, mmdb::Atom *a2, const double &o) {
	 atom_1 = a1;
	 atom_2 = a2;
	 overlap_volume = o;
      }
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
      
      // for energy types -> vdw radius and h-bond type
      std::map<std::string, double> type_to_vdw_radius_map;
      std::map<mmdb::Atom *, double> central_residue_atoms_vdw_radius_map; // ligand atoms
      std::map<mmdb::Atom *, double> neighbour_atoms_vdw_radius_map; // neighbouring atoms
      dictionary_residue_restraints_t central_residue_dictionary;
      std::vector<dictionary_residue_restraints_t> neighb_dictionaries;
      double get_vdw_radius_ligand_atom(mmdb::Atom *at);
      double get_vdw_radius_neighb_atom(mmdb::Atom *at, unsigned int idx_neighb);
      double get_overlap_volume(const double &dist, const double &r_1, const double &r_2) const; // in A^3
      const protein_geometry *geom_p;

      // first is yes/no, second is if the H is on the ligand
      // 
      std::pair<bool, bool> is_h_bond_H_and_acceptor(mmdb::Atom *ligand_atom,
				    mmdb::Atom *env_atom,
				    const double &d) const;
      hb_t get_h_bond_type(mmdb::Atom *at);
      
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
      std::vector<atom_overlap_t> overlaps;
      void make_overlaps();
   };

}


#endif // ATOM_OVERLAPS_HH
