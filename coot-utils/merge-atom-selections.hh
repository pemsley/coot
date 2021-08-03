
#ifndef MERGE_ATOM_SELECTIONS_HH
#define MERGE_ATOM_SELECTIONS_HH

#include <utility>
#include <vector>
#include <string>
#include <mmdb2/mmdb_manager.h>


namespace coot {


   class delete_a_chain_t {
   public:
      enum delete_a_chain_type_t { NONE, DELETE_FIRST_CHAIN, DELETE_SECOND_CHAIN};
      bool chains_were_mergeable;
      bool short_fragment_is_in_first_selection;
      bool short_fragment_is_upstream_fragment;
      delete_a_chain_type_t delete_type;
   public:
      delete_a_chain_t(bool a, bool b, bool c) :
         chains_were_mergeable(a), short_fragment_is_in_first_selection(b), short_fragment_is_upstream_fragment(c) {
         delete_type = NONE;
      }
   };

   // make this a member function of delete_a_chain_t?
   void delete_the_short_overlapping_chain(delete_a_chain_t dac,
                                           mmdb::Manager *mol,
                                           const std::string &chain_id_i_chain,
                                           const std::string &chain_id_j_chain);

   class match_container_for_residues_t {
      void meld_residues(std::vector<mmdb::Residue *> res_vec, mmdb::Residue *residue_2,
			                int res_no_delta, mmdb::Chain *to_chain_p, mmdb::Manager *mol);
   public:
      mmdb::Residue *residue_1;
      mmdb::Residue *residue_2;

      std::vector<mmdb::Residue *> fragment_1_res_vec;
      std::vector<mmdb::Residue *> fragment_2_res_vec;
      std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > atom_pairs;
      match_container_for_residues_t(mmdb::Residue *r1, mmdb::Residue *r2);
      match_container_for_residues_t() : residue_1(NULL), residue_2(NULL) {}
      void add(mmdb::Atom *at_1, mmdb::Atom *at_2);
      // atom_selection_1(true) vs atom_selection_2(false) and upstream(true) vs downstream (false)
      delete_a_chain_t
      find_short_fragment_around_overlap(mmdb::Manager *mol,
                                         int selection_handle_1,
                                         int selection_handle_2) const;
      void delete_upstream(mmdb::Manager *mol, bool from_first, int selection_handle_1, int selection_handle_2);
      void delete_downstream(mmdb::Manager *mol, bool from_first, int selection_handle_1, int selection_handle_2);
      // merge_flags used as in find_short_fragment_around_overlap()
      void meld(mmdb::Manager *mol, std::pair<bool, bool> merge_flags);
      std::vector<mmdb::Residue *> residue_vector_from_residue(mmdb::Manager *mol, mmdb::Residue *residue_p) const;
      void debug() const;
   };


   class match_container_t {
      public:
      std::vector<match_container_for_residues_t> matches;
      void add(mmdb::Atom *at_1, mmdb::Atom *at_2);
      // return a null for residue_1 on failure
      match_container_for_residues_t find_best_match() const;
   };

   // tests if the 2 selections have overlapping atoms
   // how about std::pair<bool, std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > >?
   // At one particular residue, the atoms overlap, all residues downstream in one selection
   // and upstream in the other should be deleted - and the atoms of the merging residue
   // should be averaged.
   std::pair<bool, match_container_for_residues_t>
   mergeable_atom_selections(mmdb::Manager *mol, int selection_handle_1, int selection_handle_2);

   // merge selection 2 into 1 and renumber if necessary - delete overlapping atoms.
   delete_a_chain_t
   merge_atom_selections(mmdb::Manager *mol, int selection_handle_1, int selection_handle_2);

   void merge_atom_selections(mmdb::Manager *mol);

   std::vector<mmdb::Residue *> atom_selection_to_residue_vector(mmdb::Manager *mol, int selection_handle);

   void delete_the_matched_residues_matched_residue(mmdb::Manager *mol, match_container_for_residues_t m,
                                                    bool short_fragment_is_in_first_selection);

   // maybe be a regular coot-util function?
   void renumber_chains_start_at_least_at_1(mmdb::Manager *mol);

}

#endif // MERGE_ATOM_SELECTIONS_HH

