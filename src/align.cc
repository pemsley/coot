/* src/align.cc
 * 
 * Copyright 2006 The University of York
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

#include <iostream>
#include <string>
#include <vector>

#include "coot-utils/coot-coord-utils.hh"

#include <mmdb2/mmdb_math_align.h>
#include <mmdb2/mmdb_tables.h>

#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "utils/coot-utils.hh"
#include "coot-align.hh"

void do_alignment(atom_selection_container_t asc);
coot::chain_mutation_info_container_t align_on_chain(mmdb::Chain *chain_p, mmdb::PResidue *SelResidues, int nSelResidues, const std::string &target);
std::string make_model_string(mmdb::PResidue *SelResidues, int nSelResidues);

int
main(int argc, char **argv) {

    int i = int (30 * 0.41);
    std::cout << "DEBUG:: i " << i << std::endl;
    i = int(-0.1); 
    std::cout << "DEBUG:: i " << i << std::endl;
    i = int(-0.8); 
    std::cout << "DEBUG:: i " << i << std::endl;
    i = int(0.8); 
    std::cout << "DEBUG:: i " << i << std::endl;
    i = int(rint(0.6)); 
    std::cout << "DEBUG:: i " << i << std::endl;
    i = int(rint(-0.6)); 
    std::cout << "DEBUG:: i " << i << std::endl;
    

    if (argc > 1) { 
       atom_selection_container_t asc = get_atom_selection(argv[1], true, false, false);
       if (asc.n_selected_atoms > 0) { 
	  do_alignment(asc);
       }
    }
    return 0;
}


void
do_alignment(atom_selection_container_t asc) {

   std::string target("DTSGTVCLSLPPEAIASDGPFPYSQDGTERFDSCVNCAWQRTGVVFQNRESVLPTQSYGYYHEYTVITP");

   mmdb::Manager *mol = asc.mol;
   int n_models = mol->GetNumberOfModels();
   // int max_chain_length = coot::util::max_number_of_residues_in_chain(mol);
   for (int imodel = 1; imodel <= n_models; imodel++) { 
      mmdb::Model *model_p = mol->GetModel(imodel);
      mmdb::Chain *chain_p;
      const char *chain_id;
      int n_chains = model_p->GetNumberOfChains();

      for (int ich=0; ich<n_chains; ich++) {
	 chain_p = model_p->GetChain(ich);
	 chain_id = chain_p->GetChainID();
	 std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
	 if (m.first) {
	    // int offset = m.second;
	    int selHnd = mol->NewSelection();
	    mmdb::PResidue *SelResidues = NULL;
	    int nSelResidues;

	    mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
			chain_id,
			mmdb::ANY_RES, "*",
			mmdb::ANY_RES, "*",
			"*",  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			mmdb::SKEY_NEW // selection key
			);
	    mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

	    if (nSelResidues > 0) {

	       align_on_chain(chain_p, SelResidues, nSelResidues, target);
	    }
	    mol->DeleteSelection(selHnd);
	 }
      }
   }
}

coot::chain_mutation_info_container_t
align_on_chain(mmdb::Chain *chain_p, mmdb::PResidue *SelResidues, int nSelResidues,
	       const std::string &target)  {

   coot::chain_mutation_info_container_t ch_info(chain_p->GetChainID());

   std::cout << "\n\n----------------------------------------------------" << std::endl;
   std::cout << "           chain " << chain_p->GetChainID() << std::endl;
   std::cout << "----------------------------------------------------\n\n" << std::endl;

   // std::string model ("DVSGTVCLSALPPEATDTLNLIASDGPFPYSQDFQNRESVLPTQSYGYYHEYTVITP");
   std::string model = make_model_string(SelResidues, nSelResidues);


   mmdb::math::Alignment align;

   // default values (it seems)
   mmdb::realtype wgap = 0.0;
   mmdb::realtype wspace = -1.0;

   std::string stripped_target = coot::util::remove_whitespace(target);
   
   align.Align(model.c_str(), stripped_target.c_str());

   std::cout << "model : " << align.GetAlignedS() << std::endl;
   std::cout << "target: " << align.GetAlignedT() << std::endl;
   std::cout << "score " << align.GetScore() << std::endl;

   // get something like:
   // DVSGTVCLSALPPEATDTLNLIASDGPFPYSQD----F--------Q-------NRESVLPTQSYGYYHEYTVITP
   // DTSGTVCLS-LPPEA------IASDGPFPYSQDGTERFDSCVNCAWQRTGVVFQNRESVLPTQSYGYYHEYTVITP

   // So GetAlignedS() and GetAlignedT() return strings of the same length
   //
   // We run across the length of the returned string, looking for differences:

   // model  has residue             :  Simple mutate.
   // target has different residue   :  

   // model  has residue             :  There is an insertion in the model.
   // target has "-"                 :  Delete that residue, offset residues numbers 
   //                                :  to the right by -1.
   // 
   // model  has "-"                 :  The is a deletion in the model.
   // target has residue             :  Insert a gap for the residue by offsetting 
   //                                :  residue numbers to the right by +1.

   std::string s=align.GetAlignedS();
   std::string t=align.GetAlignedT();


   if (s.length() == t.length()) {

      std::vector<int> resno_offsets(s.length(), 0);
      for (unsigned int ires=0; ires<s.length(); ires++) {
	 if (s[ires] != t[ires]) {

	    // Case 1: (simple mutate)
	    if ((s[ires] != '-') && t[ires] != '-') {
	       std::cout << "mutate ires " << ires << " "
			 << s[ires] << " to " << t[ires] << std::endl;
	       std::string target_type =
		  coot::util::single_letter_to_3_letter_code(t[ires]);
	       coot::residue_spec_t res_spec(ires);
	       ch_info.add_mutation(res_spec, target_type);
	    }

	    // Case 2: model had insertion
	    if ((s[ires] != '-') && t[ires] == '-') {
	       for (unsigned int i=ires+1; i<s.length(); i++)
		  resno_offsets[i] -= 1;
	       std::cout << "Delete ires " << ires << " " << s[ires] << std::endl;
	       coot::residue_spec_t res_spec(ires);
	       ch_info.add_deletion(res_spec);
	    }

	    // Case 3: model has a deletion
	    if ((s[ires] == '-') && t[ires] != '-') {
	       for (unsigned int i=ires+1; i<s.length(); i++)
		  resno_offsets[i] += 1;
	       coot::residue_spec_t res_spec(ires);
	       std::string target_type =
		  coot::util::single_letter_to_3_letter_code(t[ires]);
	       ch_info.add_insertion(res_spec, target_type);
	       std::cout << "Insert ires " << ires << " " << t[ires] << std::endl;
	    }
	 }
      }
   }
   ch_info.rationalize_insertions();
   return ch_info;
}

std::string make_model_string(mmdb::PResidue *SelResidues, int nSelResidues) {

   std::string s;
   std::string this_residue;
   char r[10];
   
   for (int i=0; i<nSelResidues; i++) {
      
      this_residue = "X";
      mmdb::pstr rn = SelResidues[i]->GetResName();
      std::string residue_name(rn);
      mmdb::Get1LetterCode(rn, r);
      this_residue = r[0];
      if (residue_name != "HOH") 
	 s += this_residue;
   }

   std::cout << "DEBUG:: make_model_string returns " << s << std::endl;
   return s;
}
