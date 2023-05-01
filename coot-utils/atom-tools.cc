
#include "geometry/residue-and-atom-specs.hh"
#include "atom-tools.hh"

// So that we don't draw the atoms in the "static" molecules that have intermediate atoms.
// - should make the view more clear.
//
// maybe we should return a vector of bools, a vector that corresponds to the atoms in mol_atom_sel
//
std::set<int>
coot::atom_indices_in_other_molecule(atom_selection_container_t mol_atom_sel,
				     atom_selection_container_t moving_atom_sel) {

   std::set<int> s;
   bool debug = false;

   if (false) {
      std::cout << "Here with moving_atom_sel SelectionHandle    "    << moving_atom_sel.SelectionHandle << std::endl;
      std::cout << "Here with moving_atom_sel UDDAtomIndexHandle "    << moving_atom_sel.UDDAtomIndexHandle << std::endl;
      std::cout << "Here with moving_atom_sel UDDOldAtomIndexHandle " << moving_atom_sel.UDDOldAtomIndexHandle << std::endl;
      
      std::cout << "Here with mol_atom_sel SelectionHandle    "    << mol_atom_sel.SelectionHandle << std::endl;
      std::cout << "Here with mol_atom_sel UDDAtomIndexHandle "    << mol_atom_sel.UDDAtomIndexHandle << std::endl;
      std::cout << "Here with mol_atom_sel UDDOldAtomIndexHandle " << mol_atom_sel.UDDOldAtomIndexHandle << std::endl;
   }

   int mol_udd_atom_index_handle = mol_atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
   
   for (int iat=0; iat<moving_atom_sel.n_selected_atoms; iat++) {
      mmdb::Atom *at = moving_atom_sel.atom_selection[iat];
      if (debug)
	 std::cout << "debug:: in atom_indices_in_other_molecule() with moving atom " << atom_spec_t(at) << std::endl;
      int idx = -1;
      if (at->GetUDData(moving_atom_sel.UDDOldAtomIndexHandle, idx) == mmdb::UDDATA_Ok) {
	 if ((idx < mol_atom_sel.n_selected_atoms) && (idx != -1)) {
	    mmdb::Atom *at_mol = mol_atom_sel.atom_selection[idx];
	    atom_spec_t moving_atom_spec(at);
	    atom_spec_t mol_atom_spec(at_mol);
	    if (moving_atom_spec.atom_name == mol_atom_spec.atom_name) {
	       if (moving_atom_spec.res_no == mol_atom_spec.res_no) {
		  if (moving_atom_spec.chain_id == mol_atom_spec.chain_id) {
		     int udd_idx_atom;
		     if (debug)
			std::cout << "debug:: here with mol_udd_atom_index_handle " << mol_udd_atom_index_handle << std::endl;
		     if (at_mol->GetUDData(mol_udd_atom_index_handle, udd_idx_atom) == mmdb::UDDATA_Ok) {
			if (udd_idx_atom == idx) {
			   s.insert(idx);
			   if (debug)
			      std::cout << "INFO:: atom_indices_in_other_molecule " << iat << " "
					<< moving_atom_spec << " " << idx << std::endl;
			} else {
			   std::cout << "WARNING:: atom_indices_in_other_molecule() rejecting atom from set because "
				     << udd_idx_atom << " is not " << idx << std::endl;
			}
		     } else {
			std::cout << "WARNING:: atom_indices_in_other_molecule() GetUDData failure " << mol_udd_atom_index_handle
				  << std::endl;
		     }
		  } else {
		     std::cout << "WARNING:: atom_indices_in_other_molecule not same chain id " << moving_atom_spec << std::endl;
		  }
	       } else {
		  std::cout << "WARNING:: atom_indices_in_other_molecule not same res_no " << moving_atom_spec << std::endl;
	       }
	    } else {
	       std::cout << "WARNING:: atom_indices_in_other_molecule not same atom_name "
			 << moving_atom_spec << " " << mol_atom_spec << std::endl;
	    }
	 } else {
	    std::cout << "WARNING:: atom_indices_in_other_molecule - bad atom index " << idx << " "
		      << mol_atom_sel.n_selected_atoms << std::endl;
	 }
      } else {
	 if (debug)
	    std::cout << "WARNING:: atom_indices_in_other_molecule " << iat << " " << atom_spec_t(at) << " idx fail" << std::endl;
      }
   }

   return s;
}


#include "coot-coord-utils.hh" // for co

void
coot::find_out_of_register_errors(mmdb::Manager *mutations_mol, mmdb::Manager *ref_mol) {

   // we have done mutations post alignment. There might be out-of-register errors
   // (probably are). I want to find them.
   // The ref_mol is the reference/deposited structure.
   // The mutated molecule will not have the correct residue numbering
   // and will have mis-built loops (too many or too few residues).
   //
   // So find 5-residue fragments in the reference, and try to find that
   // sequence in the mutations mol. Then find the difference of the positions
   // based around the central residue (say).
   //
   // Which chains? The first chains in both molecules

   if (! mutations_mol) return;
   if (! ref_mol) return;

   class sequence_info_t {
   public:
      sequence_info_t(mmdb::Residue *r, int rn, const std::string &resname) :
         residue_p(r), res_no(rn), res_name(resname) {}
      mmdb::Residue *residue_p;
      int res_no;
      std::string res_name;
   };

   auto fill_sequence = [] (mmdb::Manager *mol) {

                           std::vector<sequence_info_t> sequence;
                           int imodel = 1;
                           mmdb::Model *model_p = mol->GetModel(imodel);
                           if (model_p) {
                              int n_chains = model_p->GetNumberOfChains();
                              for (int ichain=0; ichain<n_chains; ichain++) {
                                 mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                 int nres = chain_p->GetNumberOfResidues();
                                 for (int ires=0; ires<nres; ires++) {
                                    mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                    int res_no = residue_p->GetSeqNum();
                                    std::string rn = residue_p->GetResName();
                                    sequence_info_t si(residue_p, res_no, rn);
                                    sequence.push_back(si);
                                 }
                                 if (! sequence.empty())
                                    break;
                              }
                           }
                           return sequence;
                        };

   // return empty vector on failure
   auto make_a_5_residue_fragment = [] (int i, const std::vector<sequence_info_t> &si) {
                                       std::vector<sequence_info_t> s;
                                       int si_size = si.size();
                                       int i_start = i-2;
                                       if (i_start >= 0) {
                                          int i_last = i + 2;   // last index in the fragment
                                          int i_end = i_last + 1;  // first invalid index
                                          if (i_end < si_size) {
                                             std::vector<sequence_info_t>::const_iterator it_begin = si.begin() + i_start;
                                             std::vector<sequence_info_t>::const_iterator it_end   = si.begin() + i_end;
                                             s = std::vector<sequence_info_t>(it_begin, it_end);
                                          }
                                       }
                                       return s;
                                    };
   // return empty vector on failure
   auto find_fragment_in_sequence = [] (const std::vector<sequence_info_t> &search_frag,
                                        const std::vector<sequence_info_t> &sequence) {
                                       std::vector<sequence_info_t> s;
                                       for (unsigned int i=2; i<sequence.size(); i++) {
                                          unsigned int n_matches = 0;
                                          if (i < (sequence.size()-2)) {
                                             for (unsigned int j=0; j<search_frag.size(); j++) {
                                                if (search_frag[j].res_name == sequence[i+j-2].res_name)
                                                   n_matches++;
                                             }
                                          }
                                          if (n_matches == 5) {
                                             if (s.empty()) {
                                                std::vector<sequence_info_t>::const_iterator it_begin = sequence.begin() + i-2;
                                                std::vector<sequence_info_t>::const_iterator it_end   = sequence.begin() + i+3;
                                                s = std::vector<sequence_info_t>(it_begin, it_end);
                                             } else {
                                                std::cout << "Double hit for fragment sequence returning first hit " << std::endl;
                                             }
                                          }
                                       }
                                       return s;
                                    };

   auto get_positions_difference = [] (const std::vector<sequence_info_t> &frag_1,
                                       const std::vector<sequence_info_t> &frag_2) {
                                      float d = 0.0f;
                                      unsigned int n_hits = 0;
                                      for (unsigned int i=0; i<frag_1.size(); i++) {
                                         mmdb::Residue *res = frag_1[i].residue_p;
                                         mmdb::Atom *at_1 = res->GetAtom(" CA ");
                                         if (at_1) {
                                            res = frag_2[i].residue_p;
                                            mmdb::Atom *at_2 = res->GetAtom(" CA ");
                                            if (at_2) {
                                               clipper::Coord_orth pos_1 = co(at_1);
                                               clipper::Coord_orth pos_2 = co(at_2);
                                               double dd = (pos_2-pos_1).lengthsq();
                                               d += sqrtf(dd);
                                               std::cout << "adding d " << sqrtf(dd) << std::endl;
                                               n_hits++;
                                            }
                                         }
                                      }
                                      if (n_hits != 5)
                                         return static_cast<float>(-1.0);
                                      return d;
                                    };

   std::vector<sequence_info_t> ref_sequence = fill_sequence(ref_mol);
   std::vector<sequence_info_t> mut_sequence = fill_sequence(mutations_mol);

   if (ref_sequence.empty()) {
      std::cout << "empty ref sequence " << std::endl;
      return;
   }
   if (mut_sequence.empty()) {
      std::cout << "empty mut sequence " << std::endl;
      return;
   }

   std::vector<std::pair<mmdb::Residue *, float> > position_difference_results;
   for (unsigned int i=2; i<ref_sequence.size(); i++) {
      std::vector<sequence_info_t> five_res_frag = make_a_5_residue_fragment(i, ref_sequence);
      if (five_res_frag.size() == 5) {
         std::vector<sequence_info_t> related_frag = find_fragment_in_sequence(five_res_frag,
                                                                               mut_sequence);
         if (related_frag.size() == 5) {
            mmdb::Residue *r = five_res_frag[2].residue_p;
            std::cout << "Checking differences for residue " << residue_spec_t(r) << std::endl;
            float positions_difference = get_positions_difference(five_res_frag, related_frag);
            float d = positions_difference;
            position_difference_results.push_back(std::pair<mmdb::Residue *, float>(r, d));
         }
      }
   }

   for (unsigned int i=0; i<position_difference_results.size(); i++) {
      std::cout << i << "   " << residue_spec_t(position_difference_results[i].first)
                << " " << position_difference_results[i].second << std::endl;
   }
   
}
