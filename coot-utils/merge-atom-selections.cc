
#include <string>
#include <set>
#include <queue>
#include <algorithm>

#include <clipper/core/coords.h>

#include "utils/coot-utils.hh"

#include "merge-atom-selections.hh"
#include "coot-coord-utils.hh"

coot::match_container_for_residues_t::match_container_for_residues_t(mmdb::Residue *r1, mmdb::Residue *r2) {

   residue_1 = r1;
   residue_2 = r2;

}

// tests if the 2 selections have overlapping atoms
// how about std::pair<bool, std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > >?
std::pair<bool, coot::match_container_for_residues_t>
coot::mergeable_atom_selections(mmdb::Manager *mol, int selection_handle_1, int selection_handle_2) {

   bool status = false;
   match_container_for_residues_t m;

   mmdb::realtype max_dist = 2.2;

   if (mol) {
      mmdb::Contact *pscontact = NULL;
      int n_contacts;
      float min_dist = 0.01;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
         for (int j=0; j<4; j++)
            my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mmdb::Atom **atom_selection_1 = 0;
      mmdb::Atom **atom_selection_2 = 0;
      int n_atoms_1;
      int n_atoms_2;
      mol->GetSelIndex(selection_handle_1, atom_selection_1, n_atoms_1);
      mol->GetSelIndex(selection_handle_2, atom_selection_2, n_atoms_2);

      mol->SeekContacts(atom_selection_1, n_atoms_1,
                        atom_selection_2, n_atoms_2,
                        0, max_dist,
                        0, // in same residue
                        pscontact, n_contacts,
                        0, &my_matt, i_contact_group);
      if (n_contacts > 0) {
         // std::cout << "in mergeable_atom_selections() n_contacts is " << n_contacts << std::endl;
         if (pscontact) {
            match_container_t match_set;
            for (int i=0; i<n_contacts; i++) {
               mmdb::Atom *at_1 = atom_selection_1[pscontact[i].id1];
               mmdb::Atom *at_2 = atom_selection_2[pscontact[i].id2];
               std::string atom_name_1 = at_1->GetAtomName();
               std::string atom_name_2 = at_2->GetAtomName();
               if (atom_name_1 == atom_name_2) {
                  std::pair<mmdb::Atom *, mmdb::Atom *> p(at_1, at_2);
                  // std::cout << "Adding pair " << atom_spec_t(at_1) << " and " << atom_spec_t(at_2) << std::endl;
                  match_set.add(at_1, at_2);
               }
            }
            // match_set is filled.

            std::cout << "in mergeable_atom_selections() match_set size " << match_set.matches.size() << std::endl;

            // which is the closest matching residue with more than 2 matches?
            // find_best_match() returns null for res_1 on failure
            match_container_for_residues_t best_match = match_set.find_best_match();
            // std::cout << "in mergeable_atom_selections() best_match " << best_match.residue_1 << std::endl;
            // best_match.debug();
            if (best_match.residue_1) {
               m = best_match;
               status = true;
               if (true)
                  std::cout << "in mergeable_atom_selections() best_match "
                            << residue_spec_t(best_match.residue_1) << " " << residue_spec_t(best_match.residue_2) << std::endl;
            }
         }
      }
   }
   return std::pair<bool, match_container_for_residues_t> (status, m);
}

void
coot::match_container_t::add(mmdb::Atom *at_1, mmdb::Atom *at_2) {

   mmdb::Residue *res_1 = at_1->residue;
   mmdb::Residue *res_2 = at_2->residue;
   if (res_1) {
      if (res_2) {
         bool added = false;
         for (unsigned int i=0; i<matches.size(); i++) {
            if (matches[i].residue_1 == res_1) {
               if (matches[i].residue_2 == res_2) {
                   matches[i].add(at_1, at_2);
                   added = true;
                   break;
               }
            }
         }
         if (! added) {
            // make a new one
            match_container_for_residues_t m(res_1, res_2);
            m.add(at_1, at_2);
            matches.push_back(m);
         }
      }
   }
}

void
coot::match_container_for_residues_t::add(mmdb::Atom *at_1, mmdb::Atom *at_2) {

   std::pair<mmdb::Atom *, mmdb::Atom *> p(at_1, at_2);
   atom_pairs.push_back(p);

}

coot::match_container_for_residues_t
coot::match_container_t::find_best_match() const {

   match_container_for_residues_t rm; // with null

   double best_av_devi = 999999999.9;
   for (unsigned int i=0; i<matches.size(); i++) {
      const match_container_for_residues_t &m = matches[i];
      if (m.atom_pairs.size() > 2) {
         double sum_devi = 0;
         for (unsigned int iat=0; iat<m.atom_pairs.size(); iat++) {
            mmdb::Atom *at_1 = m.atom_pairs[iat].first;
            mmdb::Atom *at_2 = m.atom_pairs[iat].second;
            clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
            clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
            double dd = (pt_1-pt_2).lengthsq();
            sum_devi += sqrt(dd);
         }
         double av_devi = sum_devi/static_cast<double>(m.atom_pairs.size());
         if (av_devi < best_av_devi) {
            best_av_devi = av_devi;
            rm = m;
         }
      }
   }
   return rm;
}


// atom_selection_1(true) vs atom_selection_2(false)  and  upstream(true) vs downstream (false)
// @return short-fragment-is-in-first-selection, short_fragment_is_upstream_fragment
//
std::pair<bool, bool>
coot::match_container_for_residues_t::find_short_fragment_around_overlap(mmdb::Manager *mol, int selection_handle_1, int selection_handle_2) const {

   bool is_first_selection = true;
   bool is_upstream = true;

   mmdb::Atom **atom_selection_1 = 0;
   mmdb::Atom **atom_selection_2 = 0;
   int n_atoms_1;
   int n_atoms_2;
   mol->GetSelIndex(selection_handle_1, atom_selection_1, n_atoms_1);
   mol->GetSelIndex(selection_handle_2, atom_selection_2, n_atoms_2);

   // how many atoms below the matched atoms?
   int atom_sel_1_n_above = 0;
   int atom_sel_1_n_below = 0;
   int atom_sel_2_n_above = 0;
   int atom_sel_2_n_below = 0;

   bool hit_matches = false;
   for (int iat=0; iat<n_atoms_1; iat++) {
      mmdb::Atom *at = atom_selection_1[iat];
      // bleugh
      if (! hit_matches) {
         for (unsigned int ip=0; ip<atom_pairs.size(); ip++) {
            if (atom_pairs[ip].first == at) {
               hit_matches = true;
               break;
            }
         }
      }
      if (! hit_matches)
         atom_sel_1_n_above++;
      else
         atom_sel_1_n_below++;
   }
   hit_matches = false;
   for (int iat=0; iat<n_atoms_2; iat++) {
      mmdb::Atom *at = atom_selection_2[iat];
      if (! hit_matches) {
         for (unsigned int ip=0; ip<atom_pairs.size(); ip++) {
            if (atom_pairs[ip].second == at) {
               hit_matches = true;
               break;
            }
         }
      }
      if (! hit_matches)
         atom_sel_2_n_above++;
      else
         atom_sel_2_n_below++;
   }

   // idx == 1: is_in_first = true;  is_upstream = true
   // idx == 2: is_in_first = true;  is_upstream = false
   // idx == 3: is_in_first = false; is_upstream = true
   // idx == 4: is_in_first = false; is_upstream = false

   int idx = 1;
   if (atom_sel_1_n_below <= atom_sel_1_n_above)
      if (atom_sel_1_n_below <= atom_sel_2_n_above)
         if (atom_sel_1_n_below <= atom_sel_2_n_below)
            idx = 2;
   if (atom_sel_2_n_above <= atom_sel_1_n_above)
      if (atom_sel_2_n_above <= atom_sel_1_n_below)
         if (atom_sel_2_n_above <= atom_sel_2_n_below)
            idx = 3;
   if (atom_sel_2_n_below <= atom_sel_1_n_above)
      if (atom_sel_2_n_below <= atom_sel_1_n_below)
         if (atom_sel_2_n_below <= atom_sel_2_n_above)
            idx = 4;

   if (idx==2) { is_upstream = false; }
   if (idx==3) { is_first_selection = false; }
   if (idx==4) { is_first_selection = false; is_upstream = false;}

   if (false)
      std::cout << "debug in find_short_fragment_around_overlap() "
                << "atom_sel_1_n_below " << atom_sel_1_n_below << " "
                << "atom_sel_1_n_above " << atom_sel_1_n_above << " "
                << "atom_sel_2_n_below " << atom_sel_2_n_below << " "
                << "atom_sel_2_n_above " << atom_sel_2_n_above
                << " with idx " << idx << "\n";

   return std::pair<bool, bool>(is_first_selection, is_upstream);

}

// merge selection 2 into 1 and renumber if necessary - delete overlapping atoms.
bool
coot::merge_atom_selections(mmdb::Manager *mol, int selection_handle_1, int selection_handle_2) {

   bool done_merge = false;
   std::pair<bool, match_container_for_residues_t> m = mergeable_atom_selections(mol, selection_handle_1, selection_handle_2);

   // std::cout << "debug in merge_atom_selections(): for " << selection_handle_1 << " " << selection_handle_2 << " m first is " << m.first << std::endl;

   if (m.first) {
      // either upstream or downstream of selection_handle_1 or selection_handle_2 is a short fragment.
      // delete that, which means that we keep the other fragment for that selection
      // which means that we know we keep the same stream of the other selection (and delete the other)
      std::pair<bool,bool> r = m.second.find_short_fragment_around_overlap(mol, selection_handle_1, selection_handle_2);

      if (r.first) {
         if (r.second) {
            m.second.delete_upstream(mol,   true,  selection_handle_1);
            m.second.delete_downstream(mol, false, selection_handle_2);
            // m.second.delete_downstream(mol, true,  selection_handle_1);
            // m.second.delete_upstream(mol,   false, selection_handle_2);
            m.second.meld(mol, r);
         } else {
            m.second.delete_downstream(mol, true, selection_handle_1);
            m.second.delete_upstream(mol,  false, selection_handle_2);
            m.second.meld(mol, r);
         }
      } else {

         if (r.second) {
            m.second.delete_upstream(mol,   false, selection_handle_2);
            m.second.delete_downstream(mol,  true, selection_handle_1);
            m.second.meld(mol, r);
         } else {
            m.second.delete_downstream(mol, false, selection_handle_2);
            m.second.delete_upstream(mol,    true, selection_handle_1);
            m.second.meld(mol, r);
         }
      }
   }

   done_merge = m.first;
   return done_merge;
}

void
coot::match_container_for_residues_t::delete_downstream(mmdb::Manager *mol, bool from_first, int selection_handle_1) {

   bool debug = false;

   // delete downstream in selection
   bool found_matchers = false;
   std::set<mmdb::Residue *> delete_these_residues;
   mmdb::Atom **atom_selection_1 = 0;
   int n_atoms_1;
   mol->GetSelIndex(selection_handle_1, atom_selection_1, n_atoms_1);
   mmdb::Chain *chain_p = atom_selection_1[0]->GetChain();
   mmdb::Residue *matchers_residue = 0;
   for (int i=0; i<n_atoms_1; i++) {
      mmdb::Atom *at = atom_selection_1[i];
      if (from_first) {
	 for (unsigned int ip=0; ip<atom_pairs.size(); ip++) {
	    if (atom_pairs[ip].first == at) {
	       found_matchers = true;
	       matchers_residue = at->residue;
	       break;
	    }
	 }
      } else {
	 for (unsigned int ip=0; ip<atom_pairs.size(); ip++) {
	    if (atom_pairs[ip].second == at) {
	       found_matchers = true;
	       matchers_residue = at->residue;
	       break;
	    }
	 }
      }

      // if we are *past* the matching atoms (not in them)
      if (found_matchers)
	 if (at->residue != matchers_residue)
	    delete_these_residues.insert(at->residue);
   }
   if (delete_these_residues.size() > 0) {
      std::set<mmdb::Residue *>::iterator it;
      for (it=delete_these_residues.begin(); it!=delete_these_residues.end(); it++) {
	 mmdb::Residue *r = *it;
         if (debug)
	    std::cout << "in delete_downstream() delete " << residue_spec_t(r) << std::endl;
      }
      for (it=delete_these_residues.begin(); it!=delete_these_residues.end(); it++) {
	 mmdb::Residue *r = *it;
	 // chain_p->DeleteResidue(r->GetSeqNum(), r->GetInsCode());
         delete r;
      }
      chain_p->TrimResidueTable();
   }
}

void
coot::match_container_for_residues_t::debug() const {

   std::cout << "debug this match_container_for_residues_t ";
   if (residue_1)
      std::cout << residue_spec_t(residue_1) << " ";
   else
      std::cout << "residue-1 null ";
   if (residue_2)
      std::cout << residue_spec_t(residue_2) << " ";
   else
      std::cout << "residue-2 null ";
   std::cout << "with " << atom_pairs.size() << " atom pairs" << std::endl;
   for (unsigned int i=0; i<atom_pairs.size(); i++)
      std::cout << "    " << atom_spec_t(atom_pairs[i].first) << " "
                << atom_spec_t(atom_pairs[i].second) << "\n";

}


// #include "coot-coord-extras.hh"


// fails with insertion codes - but who cares about that?
std::vector<mmdb::Residue *>
coot::match_container_for_residues_t::residue_vector_from_residue(mmdb::Manager *mol, mmdb::Residue *residue_in_p) const {

   // simple_residue_tree() will select residues of the other fragment - that is not what we want.
   // float close_dist_max = 2.0;
   // std::vector<mmdb::Residue *> v = simple_residue_tree(residue_p, mol, close_dist_max);

   std::vector<mmdb::Residue *> v;
   mmdb::Chain *chain_p = residue_in_p->GetChain();

   std::queue<mmdb::Residue *> q;
   std::map<residue_spec_t, mmdb::Residue *> m;
   q.push(residue_in_p);
   while (! q.empty()) {

      // process the first element and remove it from the queue
      mmdb::Residue *residue_p = q.front();
      int res_no = residue_p->GetSeqNum();
      q.pop();
      m[residue_spec_t(residue_p)] = residue_p;

      mmdb::Residue *rp = chain_p->GetResidue(res_no-1, "");
      mmdb::Residue *rn = chain_p->GetResidue(res_no+1, "");
      if (rp) {
         residue_spec_t spec(rp);
         if (m.find(spec) == m.end())
            q.push(rp);
      }
      if (rn) {
         residue_spec_t spec(rn);
         if (m.find(spec) == m.end())
            q.push(rn);
      }
   }
   std::map<residue_spec_t, mmdb::Residue *>::const_iterator it;
   for (it=m.begin(); it!=m.end(); it++) {
      v.push_back(it->second);
   }
   return v;
}


void
coot::match_container_for_residues_t::meld(mmdb::Manager *mol, std::pair<bool, bool> merge_flags) {

   // the small fragment info is in merge_flags

   if (true) {
      // merge 2 into 1
      // do we want to add residue from the second selection onto the
      // beginning or end of this?
      if (merge_flags.first) {

         int res_no_delta = residue_1->GetSeqNum() - residue_2->GetSeqNum();
         // std::cout << "debug in meld() A res_no_delta " << res_no_delta << std::endl;

         // upstream of selection 1 was deleted, so
         // merge "upstream" of the second selection, that is,
         // merge all of it, except the overlaped residue
         // "upstream" means that the residue numbers should be
         // based on the merging residue of selection 1.

         // selection_handle_2  and 1 is out of date after a deletion of atoms.
         //
         // std::vector<mmdb::Residue *> res_vec = atom_selection_to_residue_vector(mol, selection_handle_2);

         mmdb::Chain *to_chain_p = residue_1->GetChain();

         // select the residues of the fragment - looking for residues that are contiguous with the passed
         // residue.
         std::vector<mmdb::Residue *> res_vec = residue_vector_from_residue(mol, residue_2);

         meld_residues(res_vec, residue_2, res_no_delta, to_chain_p);

         // Now delete the matching residue from residue_2;
         delete residue_2; // but use the pointer for testing in the for loop


      } else {

         int res_no_delta = residue_2->GetSeqNum() - residue_1->GetSeqNum();
         // std::cout << "debug in meld() B res_no_delta " << res_no_delta << std::endl;
         std::vector<mmdb::Residue *> res_vec = residue_vector_from_residue(mol, residue_1);
         for (unsigned int i=0; i<res_vec.size(); i++) {
            mmdb::Residue *residue_p = res_vec[i];
            residue_p->seqNum += res_no_delta;
         }

         mmdb::Chain *to_chain_p = residue_2->GetChain();
         // select the residues of the fragment - looking for residues that are contiguous with the passed
         // residue.

         // Now delete the matching residue residue_1 (high value comment)
         delete residue_1;

         meld_residues(res_vec, residue_1, 0, to_chain_p);

         int n_residues = to_chain_p->GetNumberOfResidues();
         if (n_residues > 0) {
            int res_no_first = to_chain_p->GetResidue(0)->GetSeqNum();
            if (res_no_first < 1) {
               int res_no_delta = 1 - res_no_first;
               for (int i=0; i<n_residues; i++) {
                  mmdb::Residue *r = to_chain_p->GetResidue(i);
                  r->seqNum += res_no_delta;
               }
            }
         }
      }
   }
}


void
coot::match_container_for_residues_t::meld_residues(std::vector<mmdb::Residue *> res_vec, mmdb::Residue *residue_2,
						    int res_no_delta, mmdb::Chain *to_chain_p) {

   for (unsigned int i=0; i<res_vec.size(); i++) {
      mmdb::Residue *residue_p = res_vec[i];
      if (! residue_p) continue;
      if (residue_p != residue_2) {
         residue_spec_t spec_pre(residue_p);
         residue_p->seqNum += res_no_delta;
         residue_spec_t spec_post(residue_p);
         if (true) // debug
            std::cout << "in meld() res_no_delta " << res_no_delta << " " << " residue " << spec_pre << " becomes "
                      << spec_post << std::endl;
         int this_res_seq_num = residue_p->GetSeqNum();

         // don't make a new chain for this residue
         mmdb::Residue *residue_copy = util::deep_copy_this_residue_add_chain(residue_p, "", true, false);

         if (! residue_copy) {
            std::cout << "WARNING:: deep_copy_this_residue_add_chain() returned NULL for "
                      << residue_spec_t(residue_p) << std::endl;
         } else {

            int n_chain_residues;
            mmdb::PResidue *chain_residues;
            to_chain_p->GetResidueTable(chain_residues, n_chain_residues);
            int best_diff = 99999;
            int target_res_serial_number = -1;
            for (int iserial=0; iserial<n_chain_residues; iserial++) {
               mmdb::Residue *r = chain_residues[iserial];
               int chain_residue_seq_num = r->GetSeqNum();
               int this_diff = chain_residue_seq_num - this_res_seq_num;
               // std::cout << "   Here with iserial " << iserial << " and this_diff " << this_diff << std::endl;
               if (this_diff > 0) {
                  if (this_diff < best_diff) {
                     best_diff = this_diff;
                     target_res_serial_number = iserial;
                  }
               }
            }

            if (false)
               std::cout << "   debug in meld() for " << residue_spec_t(residue_p) << " target_res_serial_number "
                         << target_res_serial_number << " best_diff " << best_diff << std::endl;

            if (target_res_serial_number >= 0) {
               if (false)
                  std::cout << "   InsResidue() on " << residue_spec_t(residue_copy) << " to chain \"" << to_chain_p->GetChainID()
                            << "\"" << std::endl;
               to_chain_p->InsResidue(residue_copy, target_res_serial_number);
            } else {
               if (false)
                  std::cout << "   AddResidue() on " << residue_spec_t(residue_copy) << " to chain \"" << to_chain_p->GetChainID()
                            << "\"" << std::endl;
               to_chain_p->AddResidue(residue_copy);
            }

            delete residue_p;
         }
      }
   }
}

void
coot::match_container_for_residues_t::delete_upstream(mmdb::Manager *mol, bool from_first, int selection_handle_1) {

   if (false) { // debugging
      std::cout << "in delete_upstream here are the atom pairs " << std::endl;
      for (unsigned int ip=0; ip<atom_pairs.size(); ip++) {
         std::cout << "   " << atom_spec_t(atom_pairs[ip].first) << " " << atom_spec_t(atom_pairs[ip].second) << std::endl;
      }
   }

   // delete upstream in selection_1
   bool found_matchers = false;
   std::set<mmdb::Residue *> delete_these_residues;
   mmdb::Atom **atom_selection_1 = 0;
   int n_atoms_1;
   mol->GetSelIndex(selection_handle_1, atom_selection_1, n_atoms_1);
   mmdb::Chain *chain_p = atom_selection_1[0]->GetChain();
   for (int i=0; i<n_atoms_1; i++) {
      mmdb::Atom *at = atom_selection_1[i];
      for (unsigned int ip=0; ip<atom_pairs.size(); ip++) {
	 if (from_first) {
	    if (atom_pairs[ip].first == at) {
	       found_matchers = true;
	       break;
	    }
	 } else {
	    if (atom_pairs[ip].second == at) {
	       found_matchers = true;
	       break;
	    }
	 }
      }
      if (found_matchers)
	 break;

      if (false) { // debug
         mmdb::Residue *residue_p = at->residue;
         std::set<mmdb::Residue *>::const_iterator it = delete_these_residues.find(residue_p);
         if (it == delete_these_residues.end())
         std::cout << "debug in delete_upstream() i = " << i << " inserting residue for atom "
                   << atom_spec_t(at) << std::endl;
      }
      delete_these_residues.insert(at->residue);
   }
   if (delete_these_residues.size() > 0) {
      std::set<mmdb::Residue *>::iterator it;

      for (it=delete_these_residues.begin(); it!=delete_these_residues.end(); it++) {
	 mmdb::Residue *r = *it;
         if (false)
	    std::cout << "in delete_upstream() delete residue " << residue_spec_t(r) << std::endl;
      }

      for (it=delete_these_residues.begin(); it!=delete_these_residues.end(); it++) {
	 mmdb::Residue *r = *it;
	 // chain_p->DeleteResidue(r->GetSeqNum(), r->GetInsCode());
	 delete r;
      }
      chain_p->TrimResidueTable();
   }
}

void
coot::merge_atom_selections(mmdb::Manager *mol) {

   // make atom selections from chains/fragments and call merge_atom_selections(mol, handle_1, handle_2)

   // when a merge has been made, the selections need to be reset - this will probably need a
   // while (continue_merging) {..} construction.

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {

      int  merge_count = 0; // for debugging filename construction
      bool continue_merging = true;
      while (continue_merging) {
	 continue_merging = false;
	 int n_chains = model_p->GetNumberOfChains();
	 std::vector<int> atom_selection_handles;
	 std::vector<std::string> chain_ids; // for debugging
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int sel_hnd = mol->NewSelection();

	    mol->Select(sel_hnd, mmdb::STYPE_ATOM, 0,
			chain_p->GetChainID(),
			mmdb::ANY_RES, "*",  // starting res
			mmdb::ANY_RES, "*",  // ending res
			"*",  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			mmdb::SKEY_NEW // selection key
			);
	    atom_selection_handles.push_back(sel_hnd);
            chain_ids.push_back(chain_p->GetChainID());
	 }

	 bool r = false;
	 std::vector<int> atom_selections_changed;
	 for (unsigned int i=0; i<atom_selection_handles.size(); i++) {
	    for (unsigned int j=i; j<atom_selection_handles.size(); j++) {
	       if (i != j) {
		  // tests for mergeability first, did we do a merge?
		  r = merge_atom_selections(mol, atom_selection_handles[i], atom_selection_handles[j]);

                  if (false) { // debugging
                     if (r) {
                        std::string debug_file_name = "debug-merge-" + util::int_to_string(merge_count) + ".pdb";
                        mol->WritePDBASCII(debug_file_name.c_str());
                        merge_count++;
                     }
                  }
		  if (r) {
		     continue_merging = true;
		     atom_selections_changed.push_back(atom_selection_handles[i]);
		     atom_selections_changed.push_back(atom_selection_handles[j]);
		     break;
		  }
	       }
	    }
	    if (r)
	       break;
	 }

	 if (true) // crashes because molecule was modified? (20190112-PE doesn't seem to crash)
	    for (unsigned int i=0; i<atom_selection_handles.size(); i++)
	       if (std::find(atom_selections_changed.begin(), atom_selections_changed.end(), atom_selection_handles[i]) ==
		   atom_selections_changed.end())
		  mol->DeleteSelection(atom_selection_handles[i]);

	 if (r) {
	    mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
	    util::pdbcleanup_serial_residue_numbers(mol);
	    mol->FinishStructEdit();
	 }
      }
   }

}
