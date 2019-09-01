
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm> // std::find()
#include <stdlib.h> // abs()
#include "secondary-structure-headers.hh"
#include "geometry/residue-and-atom-specs.hh"

coot::secondary_structure_header_records::secondary_structure_header_records(mmdb::Manager *mol, bool needsCalcSecStructure) {

   if (mol) {
      int ss_status = 0; // OK

      if (false) {
	 int sse_as_int = mmdb::SSE_Strand;
	 std::cout << "debug: SSE_Strand " << sse_as_int << std::endl;
	 sse_as_int = mmdb::SSE_Bulge;
	 std::cout << "debug: SSE_Bulge " << sse_as_int << std::endl;
	 sse_as_int = mmdb::SSE_3Turn;
	 std::cout << "debug: SSE_3Turn " << sse_as_int << std::endl;
	 sse_as_int = mmdb::SSE_4Turn;
	 std::cout << "debug: SSE_4Turn " << sse_as_int << std::endl;
	 sse_as_int = mmdb::SSE_5Turn;
	 std::cout << "debug: SSE_5Turn " << sse_as_int << std::endl;
	 sse_as_int = mmdb::SSE_Helix;
	 std::cout << "debug: SSE_Helix " << sse_as_int << std::endl;
	 sse_as_int = mmdb::SSE_None;
	 std::cout << "debug: SSE_None " << sse_as_int << std::endl;
      }

      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
	 mmdb::Model *model_p = mol->GetModel(imod);
	 if (model_p) {
	    //
	    if (needsCalcSecStructure) {
	       ss_status = model_p->CalcSecStructure(1);
	       // c.f. mmdb::SSERC_Ok
	       // std::cout << "debug:: ss_status: " << ss_status << std::endl;
	    }
	    std::vector<helix_info_t> helices;
	    std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > strands;
	    std::vector<std::vector<mmdb::Residue *> > helices_with_residues;
	    std::vector<std::vector<mmdb::Residue *> > strands_with_residues;
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(ichain);
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::Residue *running_res_p = NULL;
	       std::vector<mmdb::Residue *> running_residues;
	       mmdb::Residue *start_res_p = NULL;
	       // weirdness - the SSE types in atom.h are mmdb::SSE_FLAG enums
	       // but the SSE attribute of a residue is a unsigned char.
	       // mmdb::SSE_FLAG current_type = mmdb::SSE_None;
	       unsigned char current_type = mmdb::SSE_None;
	       int current_type_as_int = 0; // None
	       for (int ires=0; ires<nres; ires++) {

		  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		  int sse_as_int = residue_p->SSE;

		  if (sse_as_int == 2) // bulges are strands
		     sse_as_int = 1;

		  if (ires == 0) start_res_p = residue_p;

		  if (sse_as_int != current_type_as_int) {
		     if (running_residues.size()) {

			// output previous, if start_res_p was not set, then
			// we encounted the 
			//
			if ((current_type_as_int == 6) || (current_type_as_int == 1))
			   if (false)
			      std::cout << "finished a thing " << residue_spec_t(start_res_p) << " "
					<< residue_spec_t(running_res_p) << " " << current_type_as_int
					<< std::endl;
			if (current_type_as_int == 6) {
			   helix_info_t h(start_res_p, running_res_p, running_residues.size());
			   helices.push_back(h);
			   helices_with_residues.push_back(running_residues);
			   running_residues.clear();
			}
			if (current_type_as_int == 1) {
			   if (running_residues.size() > 1) {
			      std::pair<mmdb::Residue *, mmdb::Residue *> p(start_res_p, running_res_p);
			      strands.push_back(p);
			      // std::cout << "saving a strand of size " << running_residues.size() << std::endl;
			      strands_with_residues.push_back(running_residues);
			   }
			   running_residues.clear();
			}

		     }
		     start_res_p = residue_p;
		  }

		  // now handle the last residue
		  //
		  if (ires == (nres - 1)) {
		     if (current_type_as_int == 6) {
			helix_info_t h(start_res_p, running_res_p, running_residues.size());
			helices.push_back(h);
			helices_with_residues.push_back(running_residues);
		     }
		     if (current_type_as_int == 1) {
			if (running_residues.size() > 1) {
			   std::pair<mmdb::Residue *, mmdb::Residue *> p(start_res_p, running_res_p);
			   strands.push_back(p);
			   strands_with_residues.push_back(running_residues);
			}
		     }
		  }

		  // setup for next
		  //
		  if (sse_as_int == 1 || sse_as_int == 6)
		     running_residues.push_back(residue_p);
		  running_res_p = residue_p;
		  current_type_as_int = sse_as_int;
	       }
	    }

	    // std::cout << "for model number " << imod << std::endl;
	    for (std::size_t ihel=0; ihel<helices_with_residues.size(); ihel++) {
               if (false)
	          std::cout << "  helix number  " << ihel << " with "
			    << helices_with_residues[ihel].size() << " residues" << std::endl;
	       if (false)
		  for (std::size_t ii=0; ii<helices_with_residues[ihel].size(); ii++)
		     std::cout << "      " << residue_spec_t(helices_with_residues[ihel][ii]) << std::endl;
	    }
	    for (std::size_t istrand=0; istrand<strands_with_residues.size(); istrand++) {
	       if (false)
	          std::cout << "   strand number " << istrand << " with "
			    << strands_with_residues[istrand].size() << " residues "
			    << " from " << residue_spec_t(strands_with_residues[istrand].front()) << " to "
			    << residue_spec_t(strands_with_residues[istrand].back())
			    << std::endl;
	       if (false)
		  for (std::size_t ii=0; ii<strands_with_residues[istrand].size(); ii++)
		     std::cout << "      " << residue_spec_t(strands_with_residues[istrand][ii]) << std::endl;
	    }
	    make_sheets(mol, model_p, strands_with_residues);
	    make_helices(mol, model_p, helices);
	 }
	 //
      }
   }

}

void
coot::secondary_structure_header_records::make_helices(mmdb::Manager *mol,
						       mmdb::Model *model_p,
						       const std::vector<helix_info_t> &helices_in) {

   unsigned int n_helices = helices_in.size();

   // static not dynamic
   access_model *am = static_cast<access_model *> (model_p);

   for (std::size_t i=0; i<n_helices; i++) {
      mmdb::Helix *h = new mmdb::Helix;
      mmdb::Residue *first_res = helices_in[i].start_res;
      mmdb::Residue *last_res  = helices_in[i].end_res;
      std::string helix_id = helix_index_to_helix_id(i);
      int helix_class = 1; // FIXME

      std::size_t size_of_ICode   = 10;
      std::size_t size_of_ChainID = 10;
      std::size_t size_of_ResName = 20;

      h->serNum = i+1;
      strcpy(h->helixID, helix_id.c_str());
      int ll = helix_id.size();
      if (ll < 20)
	 h->helixID[ll] = 0;

      for (std::size_t ii=0; ii<size_of_ResName; ii++)
	 h->initResName[ii] = 0;
      for (std::size_t ii=0; ii<size_of_ChainID; ii++)
	 h->initChainID[ii] = 0;
      for (std::size_t ii=0; ii<size_of_ICode; ii++)
	 h->initICode[ii] = 0;
      for (std::size_t ii=0; ii<size_of_ResName; ii++)
	 h->endResName[ii] = 0;
      for (std::size_t ii=0; ii<size_of_ChainID; ii++)
	 h->endChainID[ii] = 0;
      for (std::size_t ii=0; ii<size_of_ICode; ii++)
	 h->endICode[ii] = 0;

      strcpy(h->initResName, first_res->GetResName());
      strcpy(h->initChainID, first_res->GetChainID());
      h->initSeqNum  = first_res->GetSeqNum();
      strcpy(h->initICode, first_res->GetInsCode());

      strcpy(h->endResName, last_res->GetResName());
      strcpy(h->endChainID, last_res->GetChainID());
      h->endSeqNum  = last_res->GetSeqNum();
      strcpy(h->endICode, last_res->GetInsCode());
      h->helixClass = helix_class;
      h->length = helices_in[i].length;

      am->add_helix(h);
   }
}


void
coot::secondary_structure_header_records::make_sheets(mmdb::Manager *mol,
						      mmdb::Model *model_p,
						      const std::vector<std::vector<mmdb::Residue *> > &strands_with_residues) {

   std::vector<std::vector<strand_relation_t> > order = get_sheet_order(mol, model_p, strands_with_residues);

   mmdb::Sheets *sheets = new mmdb::Sheets;
   sheets->nSheets = order.size();
   sheets->sheet = new mmdb::PSheet[order.size()];
   unsigned int n_added_sheets = 0;

   for (std::size_t isheet=0; isheet<order.size(); isheet++) {
      const std::vector<strand_relation_t> &sheet_order = order[isheet];
      std::string id = sheet_index_to_sheet_id(n_added_sheets);
      // std::cout << "made an empty new sheet..." << " for isheet " << n_added_sheets << " with id " << id << std::endl;
      if (true) {
         mmdb::Sheet *sheet = new mmdb::Sheet;
         std::vector<mmdb::Strand *> vec_of_strands;
         // std::cout << "debug:: sheet_order has size " << sheet_order.size() << std::endl;
         for (std::size_t istrand=0; istrand<sheet_order.size(); istrand++) {
	    mmdb::Strand *strand = new mmdb::Strand;
	    strcpy(strand->sheetID, id.c_str());
	    strand->strandNo = istrand+1;
	    mmdb::Residue *first_res = NULL;
	    mmdb::Residue *last_res = NULL;
            if (false)
	       std::cout << "   debug:: strands_with_residues[sheet_order[" << istrand
		         << "].strand_idx] has size "
		         << strands_with_residues[sheet_order[istrand].strand_idx].size() << std::endl;
	    if (strands_with_residues[sheet_order[istrand].strand_idx].size() > 0) {
	       first_res = strands_with_residues[sheet_order[istrand].strand_idx][0];
	       last_res  = strands_with_residues[sheet_order[istrand].strand_idx].back();

	       strcpy(strand->initResName, first_res->GetResName());
	       strcpy(strand->initChainID, first_res->GetChainID());
	       strand->initSeqNum  = first_res->GetSeqNum();
	       strcpy(strand->initICode, first_res->GetInsCode());

	       strcpy(strand->endResName, last_res->GetResName());
	       strcpy(strand->endChainID, last_res->GetChainID());
	       strand->endSeqNum  = last_res->GetSeqNum();
	       strcpy(strand->endICode, last_res->GetInsCode());

	       strand->sense = strand_relation_t::sense_to_pdb_sense(sheet_order[istrand].sense);

	       // std::cout << "added strand to vec of strands " << std::endl;
	       vec_of_strands.push_back(strand);
	    }
	 }
         // add strands to sheet
         mmdb::PStrand *sheet_strands = new mmdb::PStrand[vec_of_strands.size()];
         for (std::size_t istrand=0; istrand<vec_of_strands.size(); istrand++)
	    sheet_strands[istrand] = vec_of_strands[istrand];
         sheet->strand = sheet_strands;
         strcpy(sheet->sheetID, id.c_str());
         sheet->nStrands = vec_of_strands.size();

         // sheet needs to be added to sheets
         sheets->sheet[isheet] = sheet;
         n_added_sheets++;
      }
   }

   // add sheet to model_p
   //
   // can't do this (sheets is protected): model_p->sheets = sheets;

   // access_model *am = access_model(model_p);

   access_model *am = static_cast<access_model *> (model_p);
   am->add_sheets(sheets);
}

std::vector<std::vector<coot::secondary_structure_header_records::strand_relation_t> >
coot::secondary_structure_header_records::get_sheet_order(mmdb::Manager *mol,
							  mmdb::Model *model_p,
							  const std::vector<std::vector<mmdb::Residue *> > &strands_with_residues) {

   std::vector<std::vector<strand_relation_t> > ordered_strands; // fill and return this

   int model_number = model_p->GetSerNum();
   std::vector<mmdb::Atom *> N_atoms_vec;
   std::vector<mmdb::Atom *> O_atoms_vec;

   for (std::size_t i=0; i<strands_with_residues.size(); i++) {
      for (std::size_t ii=0; ii<strands_with_residues[i].size(); ii++) {
	 mmdb::Residue *residue_p = strands_with_residues[i][ii];
	 mmdb::Atom **residue_atoms = 0;
	 int n_residue_atoms;
	 residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    mmdb::Atom *at = residue_atoms[iat];
	    std::string atom_name = at->GetAtomName();
	    if (atom_name == " O  ")
	       O_atoms_vec.push_back(at);
	    if (atom_name == " N  ")
	       N_atoms_vec.push_back(at);
	 }
      }
   }

   mmdb::realtype max_dist = 3.5;

   mmdb::Atom **N_atoms = new mmdb::PAtom [N_atoms_vec.size()];
   mmdb::Atom **O_atoms = new mmdb::PAtom [O_atoms_vec.size()];
   int n_N_atoms = N_atoms_vec.size();
   int n_O_atoms = O_atoms_vec.size();
   for (std::size_t i=0; i<N_atoms_vec.size(); i++)
      N_atoms[i] = N_atoms_vec[i];
   for (std::size_t i=0; i<O_atoms_vec.size(); i++)
      O_atoms[i] = O_atoms_vec[i];

   // every sheet has a strand relation set
   std::vector<std::set<strand_relation_t> > strand_ordering(strands_with_residues.size());

   {
      mmdb::Contact *pscontact = NULL;
      int n_contacts;
      float min_dist = 0.1;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
	 for (int j=0; j<4; j++)
	    my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
      mol->SeekContacts(N_atoms, n_N_atoms,
			O_atoms, n_O_atoms,
			0, max_dist,
			0, // in same residue
			pscontact, n_contacts,
			0, &my_matt, i_contact_group); // makes reverses also
      if (n_contacts > 0) {
	 if (pscontact) {
	    for (int ic=0; ic<n_contacts; ic++) {
	       mmdb::Atom *at_1 = N_atoms[pscontact[ic].id1];
	       mmdb::Atom *at_2 = O_atoms[pscontact[ic].id2];
	       mmdb::Residue *r_1 = at_1->residue;
	       mmdb::Residue *r_2 = at_2->residue;
	       if (r_1 != r_2) {
		  int rn_1 = r_1->GetSeqNum();
		  int rn_2 = r_2->GetSeqNum();
		  if (abs(rn_1 - rn_2) > 2) {
		     // std::cout << "    " << atom_spec_t(at_1) << " " << atom_spec_t(at_2) << std::endl;

		     // in which strands are those residues?
		     for (std::size_t i=0; i<strands_with_residues.size(); i++) {
			if (std::find(strands_with_residues[i].begin(),
				      strands_with_residues[i].end(),
				      r_1) != strands_with_residues[i].end()) {
			   for (std::size_t j=0; j<strands_with_residues.size(); j++) {
			      if (j != i) {
				 if (std::find(strands_with_residues[j].begin(),
					       strands_with_residues[j].end(),
					       r_2) != strands_with_residues[j].end()) {
				    if (false)
				       std::cout << "strand i " << i << " strand j " << j
						 << " connected by " << residue_spec_t(r_1)
						 << " and " << residue_spec_t(r_2)
						 << std::endl;
				    bool found = false;
				    // std::cout << "calling get_strand_sense() for " << i << " " << j << std::endl;
				    strand_relation_t::sense_t sense = strand_relation_t::get_strand_sense(strands_with_residues[i], strands_with_residues[j]);
				    strand_relation_t sr_1(j, sense);
				    strand_relation_t sr_2(i, sense);
				    
				    // either both or neither are found
				    
				    if (! found) {
				       strand_ordering[i].insert(sr_1);
				       strand_ordering[j].insert(sr_2);
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   delete [] N_atoms;
   delete [] O_atoms;

   // for a sane sheet, there should be 2 strands with 1 neighbour and all the others should have 2.

   if (false) { // debug
      for (std::size_t i=0; i<strand_ordering.size(); i++) {
	 std::cout << "strand_ordering i: " << i << " ";
	 std::set<strand_relation_t>::const_iterator it;
	 for (it=strand_ordering[i].begin(); it!=strand_ordering[i].end(); it++)
	    std::cout << "  (" << it->strand_idx << " " << it->sense << ")";
	 std::cout << std::endl;
      }
   }

   std::vector<unsigned int> done_strands;

   // search the strands for a strand with 1 neighbour.  When you find it - start a sheet.
   //
   for (std::size_t i=0; i<strand_ordering.size(); i++) {

      // thought needed for barrels - strands of which all have 2 neighbours.
      // Use a function get_first_strand() which tries to find
      // a strand with 1 neighbour, but if it fails will try any with 2.
      // What about multiple barrels? Hmm. Quite tricky.
      std::vector<strand_relation_t> sheet_order_local;

      // note: strand_ordering[i] is a strand-relation set. Every strand i has a strand_ordering[i]
      //
      if (strand_ordering[i].size() == 1) { // start of a sheet
	 if (std::find(done_strands.begin(),
		       done_strands.end(), i) == done_strands.end()) {
	    strand_relation_t sr_first(i, strand_relation_t::FIRST);
	    sheet_order_local.push_back(sr_first);
	    done_strands.push_back(i);
	    unsigned int this_strand = i;
	    bool cont = true;
	    while (cont) {
	       const std::set<strand_relation_t> &sr = strand_ordering[this_strand];
	       // find the strand in the set that is not already in ordered_strands -
	       // that's where we will go next
	       std::set<strand_relation_t>::const_iterator it;
	       bool found_one = false;
	       for (it=sr.begin(); it!=sr.end(); it++) {

		  bool found = false;
		  if (std::find(sheet_order_local.begin(), sheet_order_local.end(), *it) ==
		  sheet_order_local.end()) {

		     sheet_order_local.push_back(*it);
		     done_strands.push_back(it->strand_idx);
		     this_strand = it->strand_idx;
		     found_one = true;
		     break;
		  }
	       }
	       if (! found_one) {
		  cont = false;
	       }
	    }
	 }
      }
      if (sheet_order_local.size() > 0) {
         // std::cout << "pushing back a sheet of size " << sheet_order_local.size() << std::endl;
         ordered_strands.push_back(sheet_order_local);
      }
   }
   return ordered_strands;
}

#include <clipper/core/coords.h>
#include "coot-coord-utils.hh"

// static
coot::secondary_structure_header_records::strand_relation_t::sense_t coot::secondary_structure_header_records::strand_relation_t::get_strand_sense(const std::vector<mmdb::Residue *> &strand_1,
																		   const std::vector<mmdb::Residue *> &strand_2) {

   strand_relation_t::sense_t s(coot::secondary_structure_header_records::strand_relation_t::NO_RESULT);

   if (strand_1.size() > 1) {
      if (strand_2.size() > 1) {
	 mmdb::Residue *sr10 = strand_1[0];
	 mmdb::Residue *sr1e = strand_1.back();
	 mmdb::Residue *sr20 = strand_2[0];
	 mmdb::Residue *sr2e = strand_2.back();

	 clipper::Coord_orth sr10_pt = util::average_position(sr10);
	 clipper::Coord_orth sr1e_pt = util::average_position(sr1e);
	 clipper::Coord_orth sr20_pt = util::average_position(sr20);
	 clipper::Coord_orth sr2e_pt = util::average_position(sr2e);

	 clipper::Coord_orth v1 = sr1e_pt - sr10_pt;
	 clipper::Coord_orth v2 = sr2e_pt - sr20_pt;
	 clipper::Coord_orth v1u(v1.unit());
	 clipper::Coord_orth v2u(v2.unit());

	 double cos_theta = clipper::Coord_orth::dot(v1u, v2u);
	 // std::cout << "cos_theta for strands " << cos_theta << std::endl;
	 if (cos_theta > 0)
	    s = strand_relation_t::PARALLEL;
	 else 
	    s = strand_relation_t::ANTI_PARALLEL;
      }
   }
   return s;

}
   
