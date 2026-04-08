/*
 * src/molecule-class-info-alignment.cc
 *
 * Copyright 2018 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include "molecule-class-info.h"

// 20180302 add PIR alignment parsing
void
molecule_class_info_t::associate_pir_alignment(const std::string &chain_id, const std::string &alignment) {

   if (alignment.size()) {
      pir_alignments[chain_id] = coot::pir_alignment_t(alignment);
   }
}

std::pair<bool,int>
molecule_class_info_t::max_res_no_in_chain(mmdb::Chain *chain_p) const {

   bool status = false;
   int res_no_max = -9999;
   if (chain_p) {
      int n = chain_p->GetNumberOfResidues();
      for (int i=0; i<n; i++) { // do this backwards?
	 mmdb::Residue *residue_p = chain_p->GetResidue(i);
	 int resno_this = residue_p->GetSeqNum();
	 if (resno_this > res_no_max) {
	    res_no_max = resno_this;
	    status = true;
	 }
      }
   }
   return std::pair<bool, int> (status, res_no_max);
}

std::pair<bool,int>
molecule_class_info_t::max_res_no_in_chain(const std::string &chain_id) const {

   bool status = false;
   int res_no_max = -9999;
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::string chain_id_local = chain_p->GetChainID();
         if (chain_id_local == chain_id) {
            std::pair<bool, int> p = max_res_no_in_chain(chain_p);
            status = p.first;
            res_no_max = p.second;
         }
      }
   }
   return std::pair<bool, int> (status, res_no_max);
}

std::pair<bool,int>
molecule_class_info_t::min_res_no_in_chain(const std::string &chain_id) const {

   bool status = false;
   int res_no_min = 99999;
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         std::string chain_id_local = chain_p->GetChainID();
         if (chain_id_local == chain_id) {
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int res_no = residue_p->GetSeqNum();
                  if (res_no < res_no_min) {
                     res_no_min = res_no;
                     status = true;
                  }
               }
            }
         }
      }
   }
   return std::pair<bool, int> (status, res_no_min);
}



// apply the alignment
void
molecule_class_info_t::apply_pir_alignment(const std::string &chain_id) {

   bool debug = true;
   std::map<std::string, coot::pir_alignment_t>::const_iterator it;
   // std::cout << "INFO:: apply_pir_alignment " << chain_id << std::endl;

   it = pir_alignments.find(chain_id);
   if (it == pir_alignments.end()) {
      std::cout << "WARNING:: apply_pir_alignment() No chain \"" << chain_id << "\" found in "
                << pir_alignments.size() << " alignments" << std::endl;
   } else {

      // Happy path
      const coot::pir_alignment_t &a = it->second;

      std::cout << "DEBUG:: in apply_pir_alignment() with matches.size() " << a.matches.size() << std::endl;
 
      if (a.matches.size() > 0) {

	 mmdb::Chain *chain_p = 0;
	 int imod = 1;
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
 	 if (model_p) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       mmdb::Chain *chain_p_this = model_p->GetChain(ichain);
	       std::string this_chain_id = chain_p_this->GetChainID();
	       if (this_chain_id == chain_id) {
		  chain_p = chain_p_this;
		  break;
	       }
	    }
	 }

         std::cout << "DEBUG:: in apply_pir_alignment() with chain_p " << chain_p << std::endl;

	 if (chain_p) {

	    // so now we have the alignment and the chain to which it should be applied

	    make_backup(__FUNCTION__);
	    int bs = backups_state();
	    turn_off_backup();

	    mmdb::PResidue *residues = 0;
	    int n_residues;
	    chain_p->GetResidueTable(residues, n_residues);
	    int i_res = 0; // to start with
	    std::vector<mmdb::Residue *> deletables;

            std::cout << "DEBUG:: in apply_pir_alignment() with a.size() " << a.size() << std::endl;

	    if (a.size() > 0) {
	       if (a.size(0) > 0) {
		  std::vector<coot::pir_alignment_t::matched_residue_t> matches = a.get_matches(0);

                  std::cout << "INFO:: in apply_pir_alignment() need to apply " << matches.size()
                            << " alignment matches" << std::endl;

		  for (std::size_t i_align_pair=0; i_align_pair<matches.size(); i_align_pair++) {

		     const coot::pir_alignment_t::matched_residue_t &mr = matches[i_align_pair];

		     if (mr.aligned == '-') {

			// That means there is no residue in the model for this pair.
			// We will need to build this loop. We can't mutate it. Move along
			// the sequence

		     } else {

			std::string pir_res_type = coot::util::single_letter_to_3_letter_code(mr.aligned);

                        if (debug)
                           std::cout << "i_align_pair: " << i_align_pair << "  " << mr << std::endl;

			if (i_res < n_residues) {
			   mmdb::Residue *residue_p = residues[i_res];

			   bool found_matching_residue = false;
			   while (!found_matching_residue) {
			      if (residue_p) {
				 std::string mol_res_type = residue_p->GetResName();
				 if (debug)
				    std::cout << "looking for \"" << mr.aligned << "\" \""
					      << pir_res_type << "\" found "
					      << mol_res_type << std::endl;
                                 if (residue_p->GetSeqNum() < a.resno_start) {
				    // we haven't found the sarting residue yet
				    if (debug)
				       std::cout << "We haven't found the starting residue yet "
						 << residue_p->GetSeqNum() << " " << a.resno_start
						 << std::endl;
                                 } else {
                                     // happy Path
				    if (pir_res_type == mol_res_type) {
				       found_matching_residue = true;

                                       if (mr.target == '-') {
                                          // there was a residue in the model that we don't want in the final sequence
                                          // delete it!
					  std::string aligned_res_type = coot::util::single_letter_to_3_letter_code(mr.aligned);
					  std::string current_res_type = residue_p->GetResName();
					  if (aligned_res_type == current_res_type) {
					     deletables.push_back(residue_p);
					  } else {
					     std::cout << "Something strange on Delete "
						       << coot::residue_spec_t(residue_p) << std::endl;
					  }
                                       } else {
				          // now actually mutate (if needed)
				          std::string to = "to";
				          if (mr.aligned == mr.target) to = "..";
					  if (debug)
					     std::cout << "INFO:: mutate " << coot::residue_spec_t(residue_p)
						       << " " << residue_p->GetResName() << " from "
						       << mr.aligned << " " << to << " " << mr.target << std::endl;
				          if (mr.aligned != mr.target) {
				             std::string new_residue_type =
						coot::util::single_letter_to_3_letter_code(mr.target);
					     bool verbose_mutate = false;
				             mutate(residue_p, new_residue_type, verbose_mutate);
                                         }
                                       }
				    }
				 }
			      }

			      // try something new

			      i_res++;
			      // std::cout << "trying next residue " << i_res << std::endl;
			      if (i_res < n_residues) {
				 residue_p = residues[i_res];
			      } else {
				 break; // the while loop
			      }
			   }

			} else {
			   break; // we have gone past the last residue of the chain
			}
		     }
		  }
	       }
	    }

	    if (deletables.size()) {
	       std::vector<coot::residue_spec_t> specs;
	       for (unsigned int jj=0; jj<deletables.size(); jj++)
		  specs.push_back(coot::residue_spec_t(deletables[jj]));
	       delete_residues(specs);
	    }

	    apply_pir_renumber(a, chain_p);

	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	    if (bs)
	       turn_on_backup();
	 }
      }
   }
}

void
molecule_class_info_t::apply_pir_renumber(const coot::pir_alignment_t &a_in, mmdb::Chain *chain_p) {

   coot::pir_alignment_t a = a_in;

   std::cout << "----------------- now apply resno_start_structure from " << a.resno_start << " "
	     << a.resno_end << " to (new-start) " << a.resno_start_structure << std::endl;

   if (a.resno_start != -1) {
      if (a.resno_start_structure != -1) {

	 int resno_end = a.resno_end;

	 if (a.resno_end == -1) {
	    // set to be the C-terminus
	    std::pair<bool,int> mr = max_res_no_in_chain(chain_p);
	    if (mr.first) {
	       resno_end = mr.second;
	    }
	 }

	 int change_by_delta = a.resno_start_structure - a.resno_start;
	 std::string chain_id(chain_p->GetChainID());
	 std::cout << "apply_pir_renumber " << a.resno_start_structure << " "
		   << resno_end << " " << change_by_delta<< std::endl;
	 renumber_residue_range(chain_id, a.resno_start, resno_end, change_by_delta);
      }
   }
}
