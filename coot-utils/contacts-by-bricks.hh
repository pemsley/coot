
#ifndef CONTACTS_BY_BRICK_HH
#define CONTACTS_BY_BRICK_HH

#include <vector>
#include <set>
#include <string>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   // The returned vector is partially asymmetric: indices are added if the base/index atom is not fixed.
   // so contains M-M and M-F, - and *not* F-F or F-M.
   // For moving atoms, there is symmetry: M1-M2 exists so does M2-M1. Use the index values
   // to add them uniquely.

   // For good speed, this needs to exist for the lifetime of a restraints_container_t

   class contacts_by_bricks {
      float dist_nbc_max; // 8.0 say
      mmdb::PAtom *atoms;
      int n_atoms;
      bool only_between_different_residues_flag;
      std::vector<bool> fixed_flags; // a flag for every atom
      float lower_left[3]; // lowest coordinates (x,y,z) of bricking system
      float brick_size; // 6.0
      int range[3]; // how many bricks along x, y z, yes use int.
      unsigned int get_brick_index(mmdb::Atom *at) const;
      unsigned int idx_3d_to_idx_1d(int idx_3d[3]) const;
      void set_lower_left_and_range(mmdb::PAtom *atoms_in, int n_atoms);
      std::vector<std::set<unsigned int> > atoms_in_bricks;
      std::vector<std::vector<unsigned int> > thread_index_sets;
      void find_the_contacts_in_bricks(std::vector<std::set<unsigned int> > *vec,
				       bool only_between_different_residues_flag) const; // fill vec
      void find_the_contacts_between_bricks(std::vector<std::set<unsigned int> > *vec,
					    bool only_between_different_residues_flag) const; // fill vec
      void find_the_contacts_between_bricks_simple(std::vector<std::set<unsigned int> > *vec,
						   bool only_between_different_residues_flag) const; // fill vec
      void find_the_contacts_between_bricks_multi_thread(std::vector<std::set<unsigned int> > *vec,
							 bool only_between_different_residues_flag) const; // fill vec
      // the function for the thread:
      static void find_the_contacts_between_bricks_multi_thread_workpackage(std::vector<std::set<unsigned int> > *vec,
									    const std::vector<unsigned int> &index_set,
									    const std::vector<std::set<unsigned int> > &atoms_in_bricks,
									    const std::vector<bool> &fixed_flags,
									    const int brick_range[3],
									    mmdb::PAtom *atoms,
									    int brick_index_max,
									    float dist_max,
									    bool only_between_different_residues_flag);
      void fill_the_bricks();

   public:
      // Fill (or edit) contacts vec - which may or may not be filled already.
      // The fixed atom indices should also be passed, so that
      // fixed to fixed atom indices don't get added to the contact vector
      contacts_by_bricks(mmdb::PAtom *, int n_atoms, const std::set<unsigned int> &fixed_atom_indices);
      void set_dist_max(float dist_max_in);
      void find_the_contacts(std::vector<std::set<unsigned int> > *vec,
			     bool only_between_different_residues_flag=false);
   };

}

#endif // CONTACTS_BY_BRICK_HH
