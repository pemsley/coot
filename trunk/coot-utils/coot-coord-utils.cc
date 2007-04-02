
#include <algorithm>
#include "coot-utils.hh"
#include "coot-coord-utils.hh"
#include "mmdb_tables.h"  // for Get1LetterCode()
#include "mmdb_graph.h" // for graph matching

std::vector<std::string>
coot::util::residue_types_in_molecule(CMMDBManager *mol) { 

   std::vector<std::string> v;

   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in trim molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in residues_in_molecule: "
			    << std::endl;
	       } else { 
		  int nres = chain->GetNumberOfResidues();
		  CResidue *residue_p;
		  for (int ires=0; ires<nres; ires++) { 
		     residue_p = chain->GetResidue(ires);

//		     int n_atoms = residue_p->GetNumberOfAtoms();
// 		     for (int iat=0; iat<n_atoms; iat++) {
// 			CAtom *atom_p = residue_p->GetAtom(iat);
// 		     }

		     std::string resname = residue_p->name;

		     if (! is_member_p(v, resname)) { 
			v.push_back(resname);
		     } 
		  }
	       }
	    }
	 }
      }
   }
   return v;
}

clipper::RTop_orth
coot::util::matrix_convert(mat44 mat) {
   
   clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
				      mat[1][0], mat[1][1], mat[1][2],
				      mat[2][0], mat[2][1], mat[2][2]);
   clipper::Coord_orth cco(mat[0][3], mat[1][3], mat[2][3]);
   return clipper::RTop_orth(clipper_mat, cco);
} 


short int
coot::is_main_chain_p(CAtom *at) { 

   std::string mol_atom_name(at->name);
   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " O  ") {
      return 1;
   } else {
      return 0;
   } 
}

short int
coot::is_main_chain_or_cb_p(CAtom *at) { 

   std::string mol_atom_name(at->name);
   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " CB " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " O  ") {
      return 1;
   } else {
      return 0;
   } 
}

// return 0 or 1
bool coot::is_hydrogen_p(CAtom *at) {

   std::string mol_atom_ele(at->element);
   if (mol_atom_ele == " H" ||
       mol_atom_ele == " D") {
      return 1;
   } else {
      return 0;
   }
}



std::vector<std::string>
coot::util::residue_types_in_chain(CChain *chain_p) {

   std::vector<std::string> v;
   
   PCResidue residue_p;
   int nres = chain_p->GetNumberOfResidues();
   for (int ires=0; ires<nres; ires++) { 
      residue_p = chain_p->GetResidue(ires);
      if (residue_p) {
	 std::string n(residue_p->name);
	 if (! is_member_p(v, n))
	    v.push_back(n);
      }
   }
   return v;
}

// Return -1 on badness 
int
coot::util::max_number_of_residues_in_chain(CMMDBManager *mol) {

   int max_number_of_residues_in_chain = -1;
   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in max_number_of_residues_in_chain "
		      << nchains << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in residues_in_molecule: "
			    << std::endl;
	       } else { 
		  int nres = chain->GetNumberOfResidues();
		  if (nres > max_number_of_residues_in_chain) {
		     max_number_of_residues_in_chain = nres;
		  }
	       }
	    }
	 }
      }
   }
   return max_number_of_residues_in_chain;
}

// Return -1 on badness.
// 
// So that we can calculate the lenght of the graph x axis - there may
// be gaps, which is why max_number_of_residues_in_chain would fail.
// 
int
coot::util::max_min_max_residue_range(CMMDBManager *mol) {

   int max_min_max = -1;
   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in max_min_max_residue_range "
		      << nchains << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in residues_in_molecule: "
			    << std::endl;
	       } else { 
		  int nres = chain->GetNumberOfResidues();
		  if (nres > 0) { 
		     int min_resno = 99999;
		     int max_resno = -99999;
		     int this_resno;
		     for (int i=0; i<nres; i++) {
			this_resno = chain->GetResidue(i)->GetSeqNum();
			if (this_resno > max_resno)
			   max_resno = this_resno;
			if (this_resno < min_resno)
			   min_resno = this_resno;
		     }
		     int range = max_resno - min_resno + 1; // residues 1, 2 -> 2
		     if (range > max_min_max)
			max_min_max = range;
		  }
	       }
	    }
	 }
      }
   }
   return max_min_max;
}

// Return a vector of residue that are in this fragment.
// Fragments are marked by consecutively numbered residues.  A
// gap in the sequence numbers marks the end/beginning of a
// fragment.
std::vector<PCResidue>
coot::util::get_residues_in_fragment(CChain *chain_p,
				     coot::residue_spec_t clicked_residue) {

   std::vector<PCResidue> r;
   int nSelResidues;
   PCResidue *SelResidues;

   chain_p->GetResidueTable(SelResidues, nSelResidues);
   for (int i=0; i<nSelResidues; i++)
      r.push_back(SelResidues[i]);

   return r; 
} 


// Return -1 on badness (actually, number of chains in the last model)
int
coot::util::number_of_chains(CMMDBManager *mol) {

   int nchains = -1;

   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 // run over chains of the existing mol
	 nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in number_of_chains "
		      << nchains << std::endl;
	 }
      }
   }
   return nchains;
}

// Return -1 on badness:
int
coot::get_selection_handle(CMMDBManager *mol, const atom_spec_t &at) {

   int SelHnd = -1;
   if (mol) { 
      SelHnd = mol->NewSelection();
      mol->SelectAtoms (SelHnd, 0, (char *) at.chain.c_str(),
			at.resno, // starting resno, an int
			(char *) at.insertion_code.c_str(), // any insertion code
			at.resno, // ending resno
			(char *) at.insertion_code.c_str(), // ending insertion code
			"*", // any residue name
			(char *) at.atom_name.c_str(), // atom name
			"*", // elements
			(char *) at.alt_conf.c_str()  // alt loc.
			);
   }
   return SelHnd;
}



std::vector<std::string>
coot::util::chains_in_molecule(CMMDBManager *mol) { 

   std::vector<std::string> v;

   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in trim molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in residues_in_molecule: "
			    << std::endl;
	       } else {
		  v.push_back(chain->GetChainID());
	       }
	    }
	 }
      }
   }
   return v;
}

std::pair<short int, int>
coot::util::min_resno_in_chain(CChain *chain_p) {

   short int found_residues = 0;
   int min_resno = 99999999;
   
   if (chain_p == NULL) {  
      // This should not be necessary. It seem to be a
      // result of mmdb corruption elsewhere - possibly
      // DeleteChain in update_molecule_to().
      std::cout << "NULL chain in residues_in_molecule: "
		<< std::endl;
   } else { 
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      int resno;
      for (int ires=0; ires<nres; ires++) {
	 residue_p = chain_p->GetResidue(ires);
	 resno = residue_p->seqNum;
	 if (resno < min_resno) {
	    min_resno = resno;
	    found_residues = 1;
	 }
      }
   }
   return std::pair<short int, int>(found_residues, min_resno);
}


std::pair<short int, int>
coot::util::max_resno_in_chain(CChain *chain_p) {

   short int found_residues = 0;
   int max_resno = -31999;
   
   if (chain_p == NULL) {  
      // This should not be necessary. It seem to be a
      // result of mmdb corruption elsewhere - possibly
      // DeleteChain in update_molecule_to().
      std::cout << "NULL chain in residues_in_molecule: "
		<< std::endl;
   } else { 
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      int resno;
      if (nres > 0) { 
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    resno = residue_p->seqNum;
	    if (resno > max_resno) {
	       max_resno = resno;
	       found_residues = 1;
	    }
	 }
      } else {
	 // there was a chain, but no residues in it. which is what
	 // happens in pointer_atom_molecule() when we create a
	 // Pointer Atom molecule.
	 //
	 // In that case, we want the added atom to have residue
	 // number 1.
	 max_resno = 0;
      } 
   }
//    std::cout << "DEBUG:: max_resno_in_chain returning " << found_residues
// 	     << " " << max_resno << std::endl;
   return std::pair<short int, int>(found_residues, max_resno);
}

int
coot::util::number_of_residues_in_molecule(CMMDBManager *mol) {

   int nres = 0;
   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain_p;
	 CResidue *residue_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in trim molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in number residues_in_molecule: "
			    << std::endl;
	       } else { 
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) { 
		     residue_p = chain_p->GetResidue(ires);
		     if (residue_p != NULL) {
			nres++;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return nres;
}


// Enough kludgey code.  Utilities.
//
// On failure, return "";
std::string
coot::util::single_letter_to_3_letter_code(char code) {

   if (code == 'G') return std::string("GLY");
   if (code == 'A') return std::string("ALA");
   if (code == 'V') return std::string("VAL");
   if (code == 'S') return std::string("SER");
   if (code == 'N') return std::string("ASN");
   if (code == 'P') return std::string("PRO");
   if (code == 'D') return std::string("ASP");
   if (code == 'C') return std::string("CYS");
   if (code == 'Q') return std::string("GLN");
   if (code == 'E') return std::string("GLU");
   if (code == 'H') return std::string("HIS");
   if (code == 'I') return std::string("ILE");
   if (code == 'L') return std::string("LEU");
   if (code == 'K') return std::string("LYS");
   if (code == 'M') return std::string("MET");
   if (code == 'F') return std::string("PHE");
   if (code == 'T') return std::string("THR");
   if (code == 'W') return std::string("TRP");
   if (code == 'Y') return std::string("TYR");
   if (code == 'R') return std::string("ARG");

   return std::string("");
}


// Match on graph
// 
// Return the orientation matrix moving res_moving to res_reference
// and a flag letting us know that the match worked OK.
// 
std::pair<bool, clipper::RTop_orth>
coot::graph_match(CResidue *res_moving,
		  CResidue *res_reference) {

   clipper::RTop_orth rtop;
   bool success = 0;

   CGraph graph1;
   CGraph graph2;

   graph1.MakeGraph(res_moving);
   graph2.MakeGraph(res_reference);

   int build_status1 = graph1.Build(1);
   int build_status2 = graph2.Build(1);

   if (build_status1 != 0) {
      std::cout << "ERROR:: build_status1: " << build_status1 << std::endl;
   } else { 
      if (build_status2 != 0) {
	 std::cout << "ERROR:: build_status2: " << build_status2 << std::endl;
      } else { 
	 CGraphMatch match;

	 match.MatchGraphs(&graph1, &graph2, 4, 0);
	 int n_match = match.GetNofMatches();
	 match.PrintMatches();

	 double best_match_sum = 1e20;
	 int best_match = -1;
	 clipper::RTop_orth best_rtop;
	 for (int imatch=0; imatch<n_match; imatch++) {
	    int n;
	    realtype p1, p2;
	    ivector FV1, FV2;
	    match.GetMatch(imatch, FV1, FV2, n, p1, p2); // n p1 p2 set
// 	    For understanding only.  
// 	    std::cout << "Match number: " << imatch << "  " << p1*100 << "% "
// 		      << p2*100 << "% "<< std::endl;
	    std::vector<clipper::Coord_orth> coords_1_local;
	    std::vector<clipper::Coord_orth> coords_2_local;
	    for (int ipair=1; ipair<=n; ipair++) {
	       PCVertex V1 = graph1.GetVertex ( FV1[ipair] );
	       PCVertex V2 = graph2.GetVertex ( FV2[ipair] );
	       if ((!V1) || (!V2))  {
		  std::cout << "Can't get vertices for match "
			    << ipair << std::endl;
	       } else  {
// 		  printf(" %4i.  [%4s] <-> [%4s]\n",
// 			 ipair, V1->GetName(), V2->GetName());
		  CAtom *at1 = res_moving->atom[V1->GetUserID()];
		  CAtom *at2 = res_reference->atom[V2->GetUserID()];
		  coords_1_local.push_back(clipper::Coord_orth(at1->x, at1->y, at1->z));
		  coords_2_local.push_back(clipper::Coord_orth(at2->x, at2->y, at2->z));
	       }
	    }
	    clipper::RTop_orth rtop_local(coords_1_local, coords_2_local);
	    
	    double dist_sum = 0.0;
	    for (unsigned int i=0; i<coords_1_local.size(); i++) {
	       dist_sum += clipper::Coord_orth::length(coords_2_local[i],
						       coords_1_local[i].transform(rtop_local));
	    }
	    if (dist_sum < best_match_sum) {
	       // Debugging
	       // std::cout << "better dist_sum: " << dist_sum << std::endl;
	       best_rtop = rtop_local;
	       best_match_sum = dist_sum;
	       best_match = imatch;
	    }
	 } // imatch loop

	 if (best_match != -1) {
	    rtop = best_rtop;
	    success = 1;
	 }
      }
   }
   return std::pair<bool, clipper::RTop_orth> (success, rtop);
}


float
coot::util::median_temperature_factor(PPCAtom atom_selection,
				      int n_atoms,
				      float low_cutoff,
				      float high_cutoff,
				      short int apply_low_cutoff,
				      short int apply_high_cutoff) {

   float this_b;
   float median = 0;
   std::vector<float> b;
   for (int i=0; i<n_atoms; i++) {
      this_b = atom_selection[i]->tempFactor;
      if ((apply_low_cutoff && (this_b > low_cutoff)) ||
	  !apply_low_cutoff) {
	 if ((apply_high_cutoff && (this_b > high_cutoff)) ||
	     !apply_high_cutoff) {
	    b.push_back(this_b);
	 }
      }
   }

   if (b.size() > 0) { 
      std::sort(b.begin(), b.end());
      int mid_pos = b.size()/2;
      median = b[mid_pos];
   }
   return median;
} 

float
coot::util::average_temperature_factor(PPCAtom atom_selection,
				       int n_atoms,
				       float low_cutoff,
				       float high_cutoff,
				       short int apply_low_cutoff,
				       short int apply_high_cutoff) {

   float this_b = 0.0;
   float b_sum = 0.0;
   int n_sum = 0;

   for (int i=0; i<n_atoms; i++) {
      this_b = atom_selection[i]->tempFactor;
      if ((apply_low_cutoff && (this_b > low_cutoff)) ||
	  !apply_low_cutoff) {
	 if ((apply_high_cutoff && (this_b > high_cutoff)) ||
	     !apply_high_cutoff) {
	    b_sum += this_b;
	    n_sum++;
	 }
      }
   }

   float mean = 0.0;
   if (n_atoms > 0)
      mean = b_sum/float(n_atoms);
   return mean;
} 




// The flanking residues (if any) are in the residue selection (SelResidues).
// The flags are not needed now we have made adjustments in the calling
// function.
// 
// create_mmdbmanager_from_res_selection must make adjustments
// 
// Note: there is now a molecule-class-info version of this - perhaps
// we should call it?  Next bug fix here: move over to the function call.
// 
// We need to pass orig_mol because of atom index transfer
std::pair<CMMDBManager *, int>
coot::util::create_mmdbmanager_from_res_selection(CMMDBManager *orig_mol,
						  PCResidue *SelResidues, 
						  int nSelResidues, 
						  int have_flanking_residue_at_start,
						  int have_flanking_residue_at_end, 
						  const std::string &altconf,
						  const std::string &chain_id_1,
						  short int residue_from_alt_conf_split_flag) {

   int start_offset = 0;
   int end_offset = 0;
   
//    if (have_flanking_residue_at_start)
//       start_offset = -1;
//    if (have_flanking_residue_at_end)
//       end_offset = +1; 

   CMMDBManager *residues_mol = new CMMDBManager;
   CModel *model = new CModel;
   CChain *chain = new CChain;
   short int whole_res_flag = 0; // not all alt confs, only this one ("A") and "".

   // For the active residue range (i.e. not the flanking residues) we only want
   // to refine the atoms that have the alt conf the same as the picked atom
   // (and that is altconf, passed here).
   // 
   // However, for *flanking residues* it's different.  Say we are refining a
   // non-split residue with alt conf "": Say that residue has a flanking
   // residue that is completely split, into A and B.  In that case we want
   // either "" or "A" for the flanking atoms.
   // 
   // And say we want to refine the A alt conf of a completely split residue
   // that has a flanking neighbour that is completely unsplit (""), we want
   // atoms that are either "A" or "".
   // 
   // So let's try setting whole_res_flag to 1 for flanking residues.

   CResidue *r;
   int atom_index_handle = residues_mol->RegisterUDInteger(UDR_ATOM, "mol's atom index");
   int afix_handle_orig = orig_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
   int afix_handle_new_mol = -1;
   if (afix_handle_orig >= 0) { 
      afix_handle_new_mol = residues_mol->RegisterUDInteger(UDR_ATOM, "shelx afix");
      int udd_afix_handle_test = residues_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
//       std::cout << "DEBUG:: in create_mmdbmanager_from_res_selection, "
// 		<< "test : " << afix_handle_new_mol << " vs " << udd_afix_handle_test << std::endl;
   }

   for (int ires=start_offset; ires<(nSelResidues + end_offset); ires++) { 

      if ( (ires == 0) || (ires == nSelResidues -1) ) { 
	 if (! residue_from_alt_conf_split_flag)
	    whole_res_flag = 1;
      } else { 
	 whole_res_flag = 0;
      }

      r = coot::util::deep_copy_this_residue_with_atom_index_and_afix_transfer(orig_mol, SelResidues[ires], altconf, whole_res_flag, atom_index_handle, afix_handle_new_mol);
      
      chain->AddResidue(r);
      r->seqNum = SelResidues[ires]->GetSeqNum();
      r->SetResName(SelResidues[ires]->GetResName());
   }
   chain->SetChainID(chain_id_1.c_str());
   model->AddChain(chain);
   residues_mol->AddModel(model);
   residues_mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   residues_mol->FinishStructEdit();

   if (afix_handle_orig >= 0) { 
      afix_handle_new_mol = residues_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
      int imod = 1;
      
      CModel *model_p = residues_mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    
// 	    for (int iat=0; iat<n_atoms; iat++) {
// 	       at = residue_p->GetAtom(iat);
// 	       int check_afix_number;
// 	       if (at->GetUDData(afix_handle_new_mol, check_afix_number) == UDDATA_Ok) {
// 		  std::cout << "atom " << at << " has afix handle " << check_afix_number
// 			    << std::endl;
// 	       } else {
// 		  std::cout << "Failed to get afix number right after set! "
// 			    << at << std::endl;
//	       }
//          }
	 }
      }
   }

//    int udd_afix_handle_inter = residues_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
//    std::cout << "DEBUG:: about to return from create_mmdbmanager_from_res_selection, "
// 	     << "udd_afix_handle_inter : " << udd_afix_handle_inter << std::endl;

   return std::pair<CMMDBManager *, int>(residues_mol, atom_index_handle);
}


// ignore atom index transfer
CMMDBManager *
coot::util::create_mmdbmanager_from_atom_selection(CMMDBManager *orig_mol,
						   int SelectionHandle) {
   
   CMMDBManager *atoms_mol = new CMMDBManager;
   CModel *model = new CModel;
   atoms_mol->AddModel(model);

   PCAtom *atoms = 0;
   int n_selected_atoms;
   orig_mol->GetSelIndex(SelectionHandle, atoms, n_selected_atoms);

   for (int iatom=0; iatom<n_selected_atoms; iatom++) {

      int atom_residue_selection_handle = atoms_mol->NewSelection();
      CAtom *at = atoms[iatom];

      // Does this new atoms residue exist in new atoms_mol?
      // if yes,
      //   add a copy of this atom to that residue
      // else
      //   does this new atom's chain exist in new atoms_mol?
      //   if yes
      //      create a new residue in that chain
      //      add a copy of this atom to that new residue
      //   else
      //      create a new chain
      //      create a new residue in that chain
      //      add a copy of this atom to that new residue
      // 
      atoms_mol->Select (atom_residue_selection_handle,
			STYPE_RESIDUE, 0, // .. TYPE, iModel
			at->GetChainID(), // Chain(s)
			at->GetSeqNum(), at->GetInsCode(),  // starting res
			at->GetSeqNum(), at->GetInsCode(),  // ending res
			at->GetResName(),  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			SKEY_NEW // selection key
			);
      PCResidue *sel_residues;
      int n_sel_residues;
      atoms_mol->GetSelIndex(atom_residue_selection_handle,
			     sel_residues, n_sel_residues);

      if (n_sel_residues > 0) {
	 CAtom *new_atom = new CAtom;
	 new_atom->Copy(at);
	 // std::cout << "      simply adding atom: " << new_atom << std::endl;
	 sel_residues[0]->AddAtom(new_atom);
      } else {
	 // residue was not found!
	 // Now we select on the chain
	 int atom_chain_selection_handle = atoms_mol->NewSelection();
	 atoms_mol->Select (atom_chain_selection_handle,
			   STYPE_CHAIN, 0, // .. TYPE, iModel
			    at->GetChainID(), // Chain(s)
			    ANY_RES, "*",  // starting res
			    ANY_RES, "*",  // ending res
			    "*",  // residue name
			    "*",  // Residue must contain this atom name?
			    "*",  // Residue must contain this Element?
			    "*",  // altLocs
			    SKEY_NEW // selection key
			    );
	 PCChain *sel_chains;
	 int n_sel_chains;
	 atoms_mol->GetSelIndex(atom_chain_selection_handle, sel_chains, n_sel_chains);

	 if (n_sel_chains > 0) {
// 	    std::cout << "creating residue in chain " << sel_chains[0]->GetChainID()
// 		      << " with specs " <<   at->GetResName() << " "
// 		      << at->GetSeqNum() << " "
// 		      << at->GetInsCode() << std::endl;
	    CResidue *residue = new CResidue(sel_chains[0],
					     at->GetResName(),
					     at->GetSeqNum(),
					     at->GetInsCode());
	    CAtom *new_atom = new CAtom; 
	    new_atom->Copy(at);
	    // std::cout << "   adding atom: " << new_atom << std::endl;
	    residue->AddAtom(new_atom);
	    
	 } else {
	    // There was no chain even.
	    CChain *chain = new CChain(model, at->GetChainID());
	    
// 	    std::cout << "   adding residue..." << std::endl;
// 	    std::cout << "creating residue in chain " << chain->GetChainID()
// 		      << " with specs " <<   at->GetResName() << " "
// 		      << at->GetSeqNum() << " "
// 		      << at->GetInsCode() << std::endl;
	    CResidue *residue = new CResidue(chain,
					     at->GetResName(),
					     at->GetSeqNum(),
					     at->GetInsCode());
	    CAtom *new_atom = new CAtom; 
	    new_atom->Copy(at);
	    // std::cout << "   adding atom: " << new_atom << std::endl;
	    residue->AddAtom(new_atom);
	 } 
      }
      atoms_mol->DeleteSelection(atom_residue_selection_handle);
   }

   // now set the spacegroup and cell:
   realtype a[6];
   realtype vol;
   int orthcode;
   orig_mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
   char *sg = orig_mol->GetSpaceGroup();
   size_t l = strlen(sg+1);
   char *sgc = new char[l];
   strcpy(sgc, sg);
   atoms_mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
   atoms_mol->SetSpaceGroup(sgc);
   atoms_mol->FinishStructEdit();
   return atoms_mol;
} 



// Note, we also create a chain and add this residue to that chain.
// We do this so that we have a holder for the segid.
// 
// whole_residue_flag: only copy atoms that are either in this altLoc,
// or has an altLoc of "".
// 
CResidue *
coot::util::deep_copy_this_residue(const CResidue *residue,
				   const std::string &altconf,
				   short int whole_residue_flag) {

   // Horrible casting to CResidue because GetSeqNum and GetAtomTable
   // are not const functions.
   // 
   CResidue *rres = new CResidue;
   CChain   *chain_p = new CChain;
   chain_p->SetChainID(((CResidue *)residue)->GetChainID());
   rres->seqNum = ((CResidue *)residue)->GetSeqNum();
   strcpy(rres->name, residue->name);

   PPCAtom residue_atoms;
   int nResidueAtoms;
   ((CResidue *)residue)->GetAtomTable(residue_atoms, nResidueAtoms);
   CAtom *atom_p;
   
   for(int iat=0; iat<nResidueAtoms; iat++) {
      std::string this_atom_alt_loc(residue_atoms[iat]->altLoc);
      if (whole_residue_flag ||
	  this_atom_alt_loc  == altconf || this_atom_alt_loc == "") { 
	 atom_p = new CAtom;
	 atom_p->Copy(residue_atoms[iat]);
	 rres->AddAtom(atom_p);
      }
   }
   chain_p->AddResidue(rres);
   return rres;
}

// Note, we also create a chain and add this residue to that chain.
// We do this so that we have a holder for the segid.
// 
// whole_residue_flag: only copy atoms that are either in this altLoc,
// or has an altLoc of "".
//
//
// Don't transfer afix number if afix_udd_handle (for the new atoms) is -1;
CResidue *
coot::util::deep_copy_this_residue_with_atom_index_and_afix_transfer(CMMDBManager *std_mol, 
								     const CResidue *residue,
								     const std::string &altconf,
								     short int whole_residue_flag,
								     int atom_index_udd_handle,
								     int afix_udd_handle) {

   // Note we pass the atom_index_udd_handle from the new residues
   // molecule and get the old one by using GetUDDHandle on the
   // old/std mol.

   // Horrible casting to CResidue because GetSeqNum and GetAtomTable
   // are not const functions.
   // 
   CResidue *rres = new CResidue;
   CChain   *chain_p = new CChain;
   chain_p->SetChainID(((CResidue *)residue)->GetChainID());
   rres->seqNum = ((CResidue *)residue)->GetSeqNum();
   strcpy(rres->name, residue->name);

   PPCAtom residue_atoms;
   int nResidueAtoms;
   ((CResidue *)residue)->GetAtomTable(residue_atoms, nResidueAtoms);
   CAtom *atom_p;

   // We are not passed the handle, we have to look it up.
   int mol_atom_index_handle = std_mol->GetUDDHandle(UDR_ATOM, "atom index");
   int mol_afix_handle = -1;
   int afix_number;
   if (afix_udd_handle >= 0) {
      mol_afix_handle = std_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
   }
      
   for(int iat=0; iat<nResidueAtoms; iat++) {
      std::string this_atom_alt_loc(residue_atoms[iat]->altLoc);
      if (whole_residue_flag ||
	  this_atom_alt_loc  == altconf || this_atom_alt_loc == "") { 
	 atom_p = new CAtom;
	 atom_p->Copy(residue_atoms[iat]);
	 int mol_atom_index = -1;
	 residue_atoms[iat]->GetUDData(mol_atom_index_handle, mol_atom_index);
	 atom_p->PutUDData(atom_index_udd_handle, mol_atom_index);
	 // and shelx afix data:
	 if (mol_afix_handle >= 0) {
	    if (residue_atoms[iat]->GetUDData(mol_afix_handle, afix_number) == UDDATA_Ok) { 
// 	       std::cout << "DEBUG:: Transfering udd afix: " << afix_number
// 			 << " using std mol_udd_handle " << mol_afix_handle
// 			 << " to new atom with afix_udd_handle: " << afix_udd_handle
// 			 << std::endl;
	       atom_p->PutUDData(afix_udd_handle, afix_number);
// 	    } else {
// 	       std::cout << "Ooops! Can get afix number from original atom"
// 			 << std::endl;
	    } 
	 } 
	 rres->AddAtom(atom_p);
      }
   }
   chain_p->AddResidue(rres);
   return rres;
}


// transform atoms in residue
void
coot::util::transform_atoms(CResidue *res, const clipper::RTop_orth &rtop) {

   PPCAtom residue_atoms;
   int natoms;
   clipper::Coord_orth co;
   clipper::Coord_orth trans_pos; 
   res->GetAtomTable(residue_atoms, natoms);
   for (int iatom=0; iatom<natoms; iatom++) {
      co = clipper::Coord_orth(residue_atoms[iatom]->x, 
			       residue_atoms[iatom]->y, 
			       residue_atoms[iatom]->z);
      trans_pos = co.transform(rtop);
      residue_atoms[iatom]->x = trans_pos.x();
      residue_atoms[iatom]->y = trans_pos.y();
      residue_atoms[iatom]->z = trans_pos.z();
   }
} 



CAtom *
coot::util::intelligent_this_residue_mmdb_atom(CResidue *res_p) {

   PCAtom *residue_atoms;
   int nResidueAtoms;
   
   res_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int i=0; i<nResidueAtoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == " CA ") {
	 return residue_atoms[i];
      }
   }

   if (nResidueAtoms > 0) {
      return residue_atoms[0];
   }

   // failure
   return NULL;

}

float
coot::util::occupancy_sum(PCAtom *atoms, int n_atoms) {

   float os = 0.0;

   for (int i=0; i<n_atoms; i++) {
      os += atoms[i]->occupancy;
   }
   return os;
} 

short int
coot::util::is_nucleotide(CResidue *residue_p) {
   
   short int nuc = 0;

   if (residue_p) { 
      std::string type(residue_p->name); // all spaces cut

      if (type == "Ad") { 
	 nuc = 1;
      } else { 
	 if (type == "Cd") { 
	    nuc = 1;
	 } else { 
	    if (type == "Gd") { 
	       nuc = 1;
	    } else { 
	       if (type == "Td") { 
		  nuc = 1;
	       } else { 
		  if (type == "Ar") { 
		     nuc = 1;
		  } else { 
		     if (type == "Cr") { 
			nuc = 1;
		     } else { 
			if (type == "Gr") { 
			   nuc = 1;
			} else { 
			   if (type == "Ur") { 
			      nuc = 1;
			   } else {
			      if (type == "Ud") {  // happens sometimes, e.g. brna.pdb
				 nuc = 1;
			      } else {
				 if (type == "Tr") {  // not very likely
				    nuc = 1;
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
   return nuc;
} 


// Return NULL on no such chain:
// 
CChain *
coot::util::chain_only_of_type(CMMDBManager *mol, const std::string &residue_type) {

   CChain *chain_p = NULL;
   CChain *single_type_chain_p = NULL;
   CModel *model_p = mol->GetModel(1);
   CResidue *residue_p;
   
   if (model_p) { 
      int nchains = model_p->GetNumberOfChains();
      for (int ich=0; ich<nchains; ich++) {
	 chain_p = model_p->GetChain(ich);
	 int nres = chain_p->GetNumberOfResidues();
	 short int all_same_type_flag = 1; 
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    std::string resname(residue_p->name);
	    if (! (resname == residue_type)) {
	       all_same_type_flag = 0;
	       break;
	    }
	 }
	 if (all_same_type_flag) { 
	    single_type_chain_p = chain_p;
	    break;
	 }
      }
   } 
   return single_type_chain_p;
} 


std::string
coot::util::three_letter_to_one_letter(const std::string &resname) {

   std::string n;
   char r[10];
   //    std::cout << "DEBUG:: resname length: " << resname << " " << resname.length() << std::endl;

   short int done_locally = 0;
   if (resname.length() == 2) {
      if ((resname[1] == 'd') ||
	  (resname[1] == 'r')) {
	 n = resname.substr(0, 1);
	 done_locally = 1;
      }
   }
   

   if (! done_locally) { 
      // This casting should not be required with a modern mmdb.
      Get1LetterCode( (char *) resname.c_str(), r);
      n = r[0];
   }

   return n;
} 

std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::util::extents(CMMDBManager *mol) {


   int selHnd = mol->NewSelection();
   mol->SelectAtoms(selHnd, 0, "*", ANY_RES, "*", ANY_RES, "*",
		    "*", "*", "*", "*");
   std::pair<clipper::Coord_orth, clipper::Coord_orth> p =
      coot::util::extents(mol, selHnd);
   mol->DeleteSelection(selHnd);
   return p;
}

std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::util::extents(CMMDBManager *mol,
		     int SelectionHandle) {

   std::pair<clipper::Coord_orth, clipper::Coord_orth> p;

   PCAtom *atoms = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(SelectionHandle, atoms, n_selected_atoms);
   float most_x = -99999;
   float most_y = -99999;
   float most_z = -99999;
   float least_x = 99999;
   float least_y = 99999;
   float least_z = 99999;
   for (int i=0; i<n_selected_atoms; i++) {
      if (atoms[i]->x > most_x) most_x = atoms[i]->x;
      if (atoms[i]->y > most_y) most_y = atoms[i]->y;
      if (atoms[i]->z > most_z) most_z = atoms[i]->z;

      if (atoms[i]->x < least_x) least_x = atoms[i]->x;
      if (atoms[i]->y < least_y) least_y = atoms[i]->y;
      if (atoms[i]->z < least_z) least_z = atoms[i]->z;
   }

   clipper::Coord_orth p1( most_x,  most_y,  most_z);
   clipper::Coord_orth p2(least_x, least_y, least_z);

   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (p2, p1);
}


// pair.second = 0 for failure
// pair.first  = 1 for success
// 
std::pair<clipper::RTop_orth, short int>
coot::util::get_ori_to_this_res(CResidue *res) {

   std::pair<clipper::RTop_orth, short int>  pair;
   // clipper::RTop_orth op;

   PPCAtom residue_atoms;
   int nResidueAtoms;
   res->GetAtomTable(residue_atoms, nResidueAtoms);
   if (nResidueAtoms == 0) {
      std::cout << "ERROR:: something broken in atom residue selection before ";
      std::cout << "get_ori_to_this_res, 0 atoms in given residue" << std::endl;
   } else {

      clipper::Coord_orth ca, c, n;
      int n_found = 0;
      bool found_ca = 0;
      bool found_c  = 0;
      bool found_n  = 0;
      
      for(int iat=0; iat<nResidueAtoms; iat++) {
	 std::string atom_name = residue_atoms[iat]->name;
	 if (atom_name == " CA ") {
	    n_found++;
	    found_ca = 1;
	    ca = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
	 if (atom_name == " C  ") {
	    n_found++;
	    found_c = 1;
	    c  = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
	 if (atom_name == " N  ") {
	    n_found++;
	    found_n = 1;
	    n  = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
      }

      if (n_found != 3) {
	 std::cout << "DISASTER! Not all necessary atoms found in residue ";
	 std::cout << res->GetChainID() << " " << res->GetSeqNum() << std::endl;
	 if (found_ca == 0)
	    std::cout << "    CA is missing " << std::endl;
	 if (found_c == 0)
	    std::cout << "    C is missing " << std::endl;
	 if (found_n == 0)
	    std::cout << "    N is missing " << std::endl;
	 pair.second = 0; // failure
	 return pair;
      }

      // testing
//       ca = clipper::Coord_orth(0,0,0);
//       n  = clipper::Coord_orth(0.87, 0, 1.23);
//       c  = clipper::Coord_orth(0.83, 0, -1.18);
//       clipper::Coord_orth cb(-1.03, -1.11, 0);

      
      // now get the vectors of the orientation:
      //
      clipper::Coord_orth can_unit = clipper::Coord_orth((n - ca).unit());
      clipper::Coord_orth cac_unit = clipper::Coord_orth((c - ca).unit());

      clipper::Coord_orth bisector ((can_unit + cac_unit).unit());
      clipper::Coord_orth diff_unit((can_unit - cac_unit).unit());

      clipper::Coord_orth cross_prod(clipper::Coord_orth::cross(diff_unit,bisector));
      clipper::Coord_orth cpu = clipper::Coord_orth(cross_prod.unit());

//       std::cout << "bisector   " << bisector.format() << std::endl;
//       std::cout << "diff_unit  " << diff_unit.format() << std::endl;
//       std::cout << "cross prod " << cross_prod.format() << std::endl;

      // bisector   -> new x axis
      // diff_unit  -> new z axis
      // cross_prod -> new y axis

      clipper::Mat33<double> m(  bisector.x(),   bisector.y(),   bisector.z(),
			              cpu.x(),        cpu.y(),        cpu.z(),
			        diff_unit.x(),  diff_unit.y(),  diff_unit.z());

      pair.first = clipper::RTop_orth(m.transpose(), ca);
      pair.second = 1;

   }
   return pair;
} 


// For sequence/sequence alignment


// Take into account the insertion code too:
std::vector<std::pair<CResidue *, int> >
coot::util::sort_residues_by_seqno(PCResidue *residues,
				   int nResidues) {

   std::vector<std::pair<CResidue *, int> > v;

   // construct the vector
   for (int i=0; i<nResidues; i++)
      v.push_back(std::pair<CResidue *, int> (residues[i], i));

//    // test
//    std::vector<int> vi;
//    vi.push_back(10);
//    vi.push_back(1);
//    vi.push_back(0);

//    std::vector<int>::iterator start = vi.begin();
//    std::vector<int>::iterator end   = vi.end();

//    sort(start, end);

   // and now with our lovely data (not int)
   // 
   std::vector<std::pair<CResidue *, int> >::iterator start = v.begin();
   std::vector<std::pair<CResidue *, int> >::iterator end   = v.end();
   std::sort(start, end, coot::util::compare_residues);
   // sort(start, end) for things that have implicit comparison function.
      
   return v;
}

bool
coot::util::compare_residues(const std::pair<CResidue *, int> &a,
			     const std::pair<CResidue *, int> &b) {

   int r1 = a.first->GetSeqNum();
   int r2 = b.first->GetSeqNum();

   if (r1 < r2) {
      return 1;
   } else {
      if (r1 > r2) {
	 return 0;
      } else {
	 std::string ins1(a.first->GetInsCode());
	 std::string ins2(b.first->GetInsCode());
	 if (ins1 > ins2) {
	    return 0;
	 } else {
	    return 1; // check
	 }
      }
   }

   return 1;
} 


// Use the results of the above to give us a sequence string:
std::string
coot::util::model_sequence(const std::vector<std::pair<CResidue *, int> > &sa) {

   std::string s;
   char r[10];
   for (unsigned int i=0; i<sa.size(); i++) {
      std::string this_residue = "X";
      Get1LetterCode(sa[i].first->GetResName(), r);
      this_residue = r[0];
      s += this_residue;
   }
   return s;
}


// residues with insertion codes
std::vector<CResidue *>
coot::util::residues_with_insertion_codes(CMMDBManager *mol) {

   std::vector<CResidue *> v;
   
   int imod = 1;
      
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 std::string inscode(residue_p->GetInsCode());
	 if (inscode != "") {
	    v.push_back(residue_p);
	 }
      }
   }
   return v;
}

// return the status too, torsion in radians.
std::pair<short int, double>
coot::util::omega_torsion(CResidue *C_residue, CResidue *N_residue, const std::string &altconf) {

   double omega = 0;
   short int istatus = 0; // initially unset
   PPCAtom C_residue_atoms = NULL;
   int nCResidueAtoms;
   C_residue->GetAtomTable(C_residue_atoms, nCResidueAtoms);
   CAtom *C_residue_CA_atom_p = NULL;
   CAtom *C_residue_C_atom_p = NULL; 

   PPCAtom N_residue_atoms = NULL;
   int nNResidueAtoms;
   N_residue->GetAtomTable(N_residue_atoms, nNResidueAtoms);
   CAtom *N_residue_CA_atom_p = NULL;
   CAtom *N_residue_N_atom_p = NULL;

   for (int i=0; i<nCResidueAtoms; i++) {
      std::string atom_name = C_residue_atoms[i]->name;
      std::string altconf_atom = C_residue_atoms[i]->altLoc;
      if (atom_name == " CA ")
	 if (altconf_atom == altconf)
	    C_residue_CA_atom_p = C_residue_atoms[i];
      if (atom_name == " C  ")
	 if (altconf_atom == altconf)
	    C_residue_C_atom_p = C_residue_atoms[i];
   }

   for (int i=0; i<nNResidueAtoms; i++) {
      std::string atom_name = N_residue_atoms[i]->name;
      std::string altconf_atom = N_residue_atoms[i]->altLoc;
      if (atom_name == " CA ")
	 if (altconf_atom == altconf)
	    N_residue_CA_atom_p = N_residue_atoms[i];
      if (atom_name == " N  ")
	 if (altconf_atom == altconf)
	    N_residue_N_atom_p = N_residue_atoms[i];
   }

   if (C_residue_CA_atom_p && C_residue_C_atom_p && N_residue_N_atom_p && N_residue_CA_atom_p) {
      clipper::Coord_orth ca1(C_residue_CA_atom_p->x,
			      C_residue_CA_atom_p->y,
			      C_residue_CA_atom_p->z);
      clipper::Coord_orth c1(C_residue_C_atom_p->x,
			     C_residue_C_atom_p->y,
			     C_residue_C_atom_p->z);
      clipper::Coord_orth ca2(N_residue_CA_atom_p->x,
			      N_residue_CA_atom_p->y,
			      N_residue_CA_atom_p->z);
      clipper::Coord_orth n2(N_residue_N_atom_p->x,
			     N_residue_N_atom_p->y,
			     N_residue_N_atom_p->z);

      omega = clipper::Coord_orth::torsion(ca1, c1, n2, ca2);
      istatus = 1;
   }
   return std::pair<short int, double> (istatus, omega);
}


// This function is based on mmdb's CATom::MakePDBAtomName(), so we share
// the copyright (with Eugene/EBI).
//
// 
// 
std::string
coot::pad_atom_name(const std::string &atom_id, const std::string &element) { 

   std::string name = atom_id;
   std::string new_name; // = atom_id;

   // std::cout << "DEBUG:: :" << atom_id << ": :" << element << ":" << std::endl;

   if (name.length() == 4) {

      new_name = name;

   } else {
   
      if (element == "") {
	 if (name.length() == 1) {
	    new_name = " ";
	    new_name += name;
	    new_name += "  ";
	 }
      
      } else { // element was defined

	 // consider the following cases of 
	 // atom_id and element:
	 //   "CA"  "C"  -> " CA "
	 //   "AC8" "C"  -> "AC8 " : in this case we post apply
	 //                          the space because the second char of
	 //                          the atom_id is equal to the
	 //                          element
	 //   "NN1" "N"  -> "NN1 " :
	 //   "NP"  "P"  -> "NP  "
	 
	 if (element.length() == 1) {
	    int k = atom_id.length();
	    if (k == 3) {
	       // std::cout << "comparing " << atom_id.substr(1,1) << " and " << element << std::endl;
	       if (atom_id.substr(1,1) == element) {
		  new_name = atom_id + " ";
	       } else {
		  new_name = " " + atom_id;
	       }
	    } else { 
	       // promote the characters one space
	       if (k==2) {
		  // e.g "NP" "P", or "CA" "C"
		  if ((atom_id.substr(1,1) == element) && (element != "H")) {
		     new_name = atom_id + "  ";
		  } else {
		     new_name = " " + atom_id;
		     new_name += " ";
		  }
	       } else {
		  new_name = " ";
		  new_name += atom_id;
		  // fill the rest with spaces
		  new_name += "  ";
	       }
	    }

	 } else {
	    // element was 2 (or more) characters
	    if (element[0] == ' ') { // unusual from dict, but usual for SHELXL reading
	       if (element[1] != name[1]) {
		  // promote the characters one space
		  new_name = " ";
		  new_name += atom_id;
		  int k = atom_id.length();
		  // fill the rest with spaces
		  if (k == 1) 
		     new_name += "  ";
		  if (k == 2)
		     new_name += " ";
	       } else {
		  int k = atom_id.length();
		  new_name = " ";
		  new_name += atom_id;
		  if (k == 1) 
		     new_name += "  ";
		  if (k == 2)
		     new_name += " ";
	       } 
	    } else {
	       // This is the usual case, there is no leading space
	       //
	       // left justify the name and pad with spaces
	       new_name = atom_id;
	       if (atom_id.size() == 1) // never happens?
		  new_name += "   ";
	       if (atom_id.size() == 2)
		  new_name += "  ";
	       if (atom_id.size() == 3)
		  new_name += " ";
	       
	    }
	 }
      }
   }

   // debug
//    CAtom at;
//    at.SetAtomName(atom_id.c_str());
//    at.SetElementName(type_symbol.c_str());
//    at.MakePDBAtomName();
//    if (new_name != std::string(at.name)) {
//       std::cout << "name pad failure, mmdb coot :" << at.name << ": :"
// 		<< new_name << ": for in_atom :" << atom_id.c_str()
// 		<< ": element :" << type_symbol.c_str() << ":" << std::endl;
//    } else {
//       std::cout << "name pad match " << new_name << std::endl;
//    }

//    std::cout << "new atom name :" << new_name << ": from :"
// 	     << atom_id << ": :" << element << ":" << std::endl;
   
   return new_name;
}

double
coot::lsq_plane_deviation(const std::vector<clipper::Coord_orth> &v,
			  const clipper::Coord_orth &pt) {

   int n_atoms = v.size();
   clipper::Coord_orth sum(0,0,0);
   for (int i=0; i<n_atoms; i++)
      sum += v[i];
   double factor = 1/double(n_atoms);
   clipper::Coord_orth midpoint(sum.x()*factor, sum.y()*factor, sum.z()*factor);
   
   clipper::Matrix<double> mat(3,3);
   for (int i=0; i<n_atoms; i++) {
      mat(0,0) += (v[i].x() - midpoint.x()) * (v[i].x() - midpoint.x());
      mat(1,1) += (v[i].y() - midpoint.y()) * (v[i].y() - midpoint.y());
      mat(2,2) += (v[i].z() - midpoint.z()) * (v[i].z() - midpoint.z());
      mat(0,1) += (v[i].x() - midpoint.x()) * (v[i].y() - midpoint.y());
      mat(0,2) += (v[i].x() - midpoint.x()) * (v[i].z() - midpoint.z());
      mat(1,2) += (v[i].y() - midpoint.y()) * (v[i].z() - midpoint.z());
   }
   mat(1,0) = mat(0,1);
   mat(2,0) = mat(0,2);
   mat(2,1) = mat(1,2);
   std::vector<double> eigens = mat.eigen(true);
   // Let's now extract the values of a,b,c normalize them
   std::vector<double> abcd;
   abcd.push_back(mat(0,0));
   abcd.push_back(mat(1,0));
   abcd.push_back(mat(2,0));
   double sqsum = 1e-20;
   
   for (int i=0; i<3; i++)
      sqsum += abcd[i] * abcd[i];
   for (int i=0; i<3; i++)
      abcd[i] /= sqsum;
   
   //set d, recall di = Axi+Byi+Czi-D, so xi = x_cen, yi = y_cen, zi = z_cen:
   abcd.push_back( abcd[0]*midpoint.x() + abcd[1]*midpoint.y() + abcd[2]*midpoint.z() );
   
   double val = abcd[0]*pt.x() + abcd[1]*pt.y() + abcd[2]*pt.z() - abcd[3];
   return val;
}

bool
coot::compare_atom_specs_user_float(const coot::atom_spec_t &a1, const coot::atom_spec_t &a2) {

   return a1.float_user_data < a2.float_user_data ? 1 : 0;

} 


      // For use with interesting-things-gui, make the list argument
      // from a vector of atom specs.
      // 
      // use the user data in the atom spec to give us the molecule number
      // and the button label
      // 
std::string
coot::util::interesting_things_list(const std::vector<atom_spec_t> &v) {

   // e.g. (list) for empty v
   // (list (list "button label" imol-no chain-id resno atom-name)
   //       (list "button label" imol-no chain-id resno atom-name)
   // )

   std::string r = " (list ";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].chain;
      atom_str += "\" ";
      atom_str += int_to_string(v[i].resno);
      atom_str += " \"";
      atom_str += v[i].insertion_code;
      atom_str += "\" \"";
      atom_str += v[i].atom_name;
      atom_str += "\" \"";
      atom_str += v[i].alt_conf;
      atom_str += " \"";

      std::string button_label("Clash gap: ");
      button_label += float_to_string(v[i].float_user_data);
      button_label += " : ";
      button_label += v[i].chain;
      button_label += " ";
      button_label += int_to_string(v[i].resno);
      button_label += " ";
      if (v[i].insertion_code != "") {
	 button_label += v[i].insertion_code;
 	 button_label += " ";
      }
      button_label += v[i].atom_name;
      if (v[i].alt_conf != "") {
	 button_label += ",";
	 button_label += v[i].alt_conf;
	 button_label += " ";
      }

      std::string s = "(list ";
      s += single_quote(button_label);
      s += " ";
      s += int_to_string(v[i].int_user_data);
      s += " ";
      s += atom_str;
      s += ")\n";
      
      r += s;
   }

   r += ")";
   return r;
}

// error_type is e.g. "Z score", "Clash gap"
// 
std::string
coot::util::interesting_things_list_with_fix(const std::vector<coot::util::atom_spec_and_button_info_t> &v,
					     const std::string error_type) {

   // e.g. (list) for empty v
   // (list (list "button label" imol-no chain-id resno atom-name)
   //       (list "button label" imol-no chain-id resno atom-name)
   // )
   //
   // if we have a fix, the callback function is not "" and then the
   // returned thing becomes:
   // 
   // (list (list "button label" imol-no chain-id resno atom-name callback-func)
   //       (list "button label" imol-no chain-id resno atom-name callback-func)
   // )
   //
   // where callback-func is e.g. (lambda() (do-180-degree-side-chain-flip 0 "A" 45 "" ""))

   std::string r = " (list ";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].as.chain;
      atom_str += "\" ";
      atom_str += int_to_string(v[i].as.resno);
      atom_str += " \"";
      atom_str += v[i].as.insertion_code;
      atom_str += "\" \"";
      atom_str += v[i].as.atom_name;
      atom_str += "\" \"";
      atom_str += v[i].as.alt_conf;
      atom_str += " \"";

      std::string button_label = v[i].button_label;

      std::string s = "(list ";
      s += single_quote(button_label);
      s += " ";
      s += int_to_string(v[i].as.int_user_data);
      s += " ";
      s += atom_str;

      if (v[i].callback_func != "") {
	 s += " ";
	 s +=  v[i].callback_func;
      }
      
      s += ")\n";
      
      r += s;
   }

   r += ")";
   return r;
}


// Here std_base is at some arbitary position when passed.
// 
void
coot::util::mutate_base(CResidue *residue, CResidue *std_base) {

   std::vector<std::string> adenine;  // Pirimidine
   adenine.push_back(" N9 ");
   adenine.push_back(" C8 ");
   adenine.push_back(" N7 ");
   adenine.push_back(" C5 ");
   adenine.push_back(" C4 ");
   // 
   adenine.push_back(" N1 ");
   adenine.push_back(" C2 ");
   adenine.push_back(" N3 ");
   adenine.push_back(" C6 ");
   adenine.push_back(" N6 ");

   std::vector<std::string> guanine; // Pirimidine
   guanine.push_back(" N9 ");
   guanine.push_back(" C8 ");
   guanine.push_back(" N7 ");
   guanine.push_back(" C5 ");
   guanine.push_back(" C4 ");
   //
   guanine.push_back(" N1 ");
   guanine.push_back(" C2 ");
   guanine.push_back(" N3 ");
   guanine.push_back(" C6 ");
   guanine.push_back(" O6 ");
   guanine.push_back(" N2 "); // No matcher for this in adenine

   std::vector<std::string> thymine;  // Purine
   thymine.push_back(" N1 ");
   thymine.push_back(" C2 ");
   thymine.push_back(" N3 ");
   thymine.push_back(" C4 ");
   thymine.push_back(" C5 ");
   thymine.push_back(" C6 ");
   // 
   thymine.push_back(" O2 ");
   thymine.push_back(" O4 ");
   thymine.push_back(" C5M");
   
   std::vector<std::string> cytosine;  // Purine
   cytosine.push_back(" N1 ");
   cytosine.push_back(" C2 ");
   cytosine.push_back(" N3 ");
   cytosine.push_back(" C4 ");
   cytosine.push_back(" C5 ");
   cytosine.push_back(" C6 ");
   // 
   cytosine.push_back(" O2 ");
   cytosine.push_back(" N4 ");
   
   std::vector<std::string> uracil;  // Purine
   uracil.push_back(" N1 ");
   uracil.push_back(" C2 ");
   uracil.push_back(" N3 ");
   uracil.push_back(" C4 ");
   uracil.push_back(" C5 ");
   uracil.push_back(" C6 ");
   // 
   uracil.push_back(" O2 ");
   uracil.push_back(" O4 ");
   

   // These next 2 are in match order, don't change it.
   std::vector<std::string> purine; // A and G
   purine.push_back(" N9 ");
   purine.push_back(" C4 ");
   purine.push_back(" C5 ");
   purine.push_back(" N7 ");
   purine.push_back(" C8 ");

   std::vector<std::string> pyrimidine; // T, U and C
   pyrimidine.push_back(" N1 ");
   pyrimidine.push_back(" C2 ");
   pyrimidine.push_back(" N3 ");
   pyrimidine.push_back(" C5 ");
   pyrimidine.push_back(" C6 ");
   pyrimidine.push_back(" C4 ");

   // We need to know whether we have purine or pyrimidine for both
   // the molecule base and the std_base.
   //
   // We need to get all (5 for pyrimidine, 6 for purine) the
   // coordinates for both bases.
   // 
   // If they are either or both are pyrimidine we match 5 atoms,
   // If they are both purine we match 6 atoms.


   // So what are the input base types?
   //
   // These for flags should be set to something after our test
   short int mol_base_is_pyrimidine = -1;
   short int mol_base_is_purine     = -1;
   short int std_base_is_pyrimidine = -1;
   short int std_base_is_purine     = -1;

   std::string mol_base_name = residue->GetResName();
   std::string std_base_name = std_base->GetResName();

   if (mol_base_name == "Ar" || mol_base_name == "Ad" ||
       mol_base_name == "Gr" || mol_base_name == "Gd") {
      mol_base_is_purine = 1;
      mol_base_is_pyrimidine = 0;
   }

   if (mol_base_name == "Cr" || mol_base_name == "Cd" ||
       mol_base_name == "Ur" || mol_base_name == "Ud" ||
       mol_base_name == "Tr" || mol_base_name == "Td") {
      mol_base_is_pyrimidine = 1;
      mol_base_is_purine = 0;
   }

   if (std_base_name == "Ar" || std_base_name == "Ad" ||
       std_base_name == "Gr" || std_base_name == "Gd") {
      std_base_is_purine = 1;
      std_base_is_pyrimidine = 0;
   }

   if (std_base_name == "Cr" || std_base_name == "Cd" ||
       std_base_name == "Tr" || std_base_name == "Td" ||
       std_base_name == "Ur" || std_base_name == "Ud") {
      std_base_is_pyrimidine = 1;
      std_base_is_purine = 0;
   }

   if ((mol_base_is_pyrimidine == -1) || (mol_base_is_purine == -1) || 
       (std_base_is_pyrimidine == -1) || (std_base_is_purine == -1) ) {

      std::cout << "ERROR:: unassigned type "
		<< "mol_base_is_pyrimidine:" << " "
		<< mol_base_is_pyrimidine << " "
		<< "mol_base_is_purine: " << " "
		<< mol_base_is_purine << " "
		<< "std_base_is_pyrimidine: " << " "
		<< std_base_is_pyrimidine << " "
		<< "std_base_is_purine: " << " "
		<< std_base_is_purine << " "
		<< residue->GetResName() << " " << std_base->GetResName()
		<< std::endl;

      
   } else {
   
//       std::cout << "DEBUG:: assigned types "
// 		<< "mol_base_is_pyrimidine:" << " "
// 		<< mol_base_is_pyrimidine << " "
// 		<< "mol_base_is_purine: " << " "
// 		<< mol_base_is_purine << " "
// 		<< "std_base_is_pyrimidine: " << " "
// 		<< std_base_is_pyrimidine << " "
// 		<< "std_base_is_purine: " << " "
// 		<< std_base_is_purine << " "
// 		<< residue->GetResName() << " " << std_base->GetResName()
// 		<< std::endl;

      int n_match_atoms = 5;
      if (mol_base_is_pyrimidine && std_base_is_pyrimidine)
	 n_match_atoms = 6;

      std::vector<std::string> moving_name_vector;
      std::vector<std::string> refrce_name_vector;

      if (std_base_is_purine)
	 moving_name_vector = purine;
      else
	 moving_name_vector = pyrimidine;

      if (mol_base_is_purine)
	 refrce_name_vector = purine;
      else
	 refrce_name_vector = pyrimidine;
      
      PCAtom *std_base_atoms;
      int n_std_base_atoms;

      PCAtom *mol_base_atoms;
      int n_mol_base_atoms;
      
      residue->GetAtomTable( mol_base_atoms, n_mol_base_atoms);
      std_base->GetAtomTable(std_base_atoms, n_std_base_atoms);

      std::vector<clipper::Coord_orth> refrce_atom_positions;
      std::vector<clipper::Coord_orth> moving_atom_positions;

//       for (unsigned int i=0; i<refrce_name_vector.size(); i++)
// 	 std::cout << "ref base search atom :" << refrce_name_vector[i]
// 		   << ":" << std::endl;
//       for (unsigned int i=0; i<moving_name_vector.size(); i++)
// 	 std::cout << "mov base search atom :" << moving_name_vector[i]
// 		   << ":" << std::endl;
      
      for (int j=0; j<n_match_atoms; j++) {
	 for (int i=0; i<n_mol_base_atoms; i++) {
	    std::string atom_name = mol_base_atoms[i]->name;
	    if (refrce_name_vector[j] == atom_name) {
	       refrce_atom_positions.push_back(clipper::Coord_orth(mol_base_atoms[i]->x,
								   mol_base_atoms[i]->y,
								   mol_base_atoms[i]->z));
	       // std::cout << "Found " << atom_name << " in reference " << std::endl;
	    }
	 }
      }

      for (int j=0; j<n_match_atoms; j++) {
	 for (int i=0; i<n_std_base_atoms; i++) {
	 std::string atom_name = std_base_atoms[i]->name;
	    if (moving_name_vector[j] == atom_name) {
	       moving_atom_positions.push_back(clipper::Coord_orth(std_base_atoms[i]->x,
								   std_base_atoms[i]->y,
								   std_base_atoms[i]->z));
	       // std::cout << "Found " << atom_name << " in moving (std) base " << std::endl;
	    }
	 }
      }

      if (int(refrce_atom_positions.size()) != n_match_atoms) {
	 std::cout << "ERROR:: wrong number of reference atoms found! "
		   << refrce_atom_positions.size() << std::endl;
      } else {

	 if (int(moving_atom_positions.size()) != n_match_atoms) {
	    std::cout << "ERROR:: wrong number of moving atoms found! "
		   << moving_atom_positions.size() << std::endl;

	 } else {

	    clipper::RTop_orth rtop (moving_atom_positions, refrce_atom_positions);
      
	    double sum_dist = 0.0;
	    double sum_dist2 = 0.0;
	    double mind =  999999999.9;
	    double maxd = -999999999.9;
	    double d;
	    for (unsigned int i=0; i<refrce_atom_positions.size(); i++) {
	       d = clipper::Coord_orth::length(refrce_atom_positions[i],
					       clipper::Coord_orth(moving_atom_positions[i].transform(rtop)));
	       sum_dist  += d;
	       sum_dist2 += d*d;
	       if (d>maxd)
		  maxd = d;
	       if (d<mind)
		  mind = d;
	    }
	    double mean = sum_dist/double(moving_atom_positions.size());
	    double var  = sum_dist2/double(moving_atom_positions.size()); // no mean*mean
	    std::cout << "INFO:: " << moving_atom_positions.size() << " matched atoms had: \n"
		      << "   mean devi: " << mean << "\n"
		      << "    rms devi: " << sqrt(var) << "\n"
		      << "    max devi: " << maxd << "\n"
		      << "    min devi: " << mind << std::endl;

	    CAtom *at;
	    // We are going to delete the current atoms of the residue
	    // and add the std_base ones.  First what *are* the atom
	    // names we what to add or delete?
	    // 
	    std::vector<std::string> mol_base_atom_names;
	    if (mol_base_name == "Ar" || mol_base_name == "Ad")
	       mol_base_atom_names = adenine;
	    if (mol_base_name == "Gr" || mol_base_name == "Gd")
	       mol_base_atom_names = guanine;
	    if (mol_base_name == "Cr" || mol_base_name == "Cd")
	       mol_base_atom_names = cytosine;
	    if (mol_base_name == "Tr" || mol_base_name == "Td")
	       mol_base_atom_names = thymine;
	    if (mol_base_name == "Ur" || mol_base_name == "Ud")
	       mol_base_atom_names = uracil;

	    if (mol_base_atom_names.size() == 0) {
	       std::cout << "ERROR:: failed to find mol_base_name for mol_base_atom_names\n";
	    } else {
	       
	       std::vector<std::string> std_base_atom_names;
	       if (std_base_name == "Ar" || std_base_name == "Ad")
		  std_base_atom_names = adenine;
	       if (std_base_name == "Gr" || std_base_name == "Gd")
		  std_base_atom_names = guanine;
	       if (std_base_name == "Cr" || std_base_name == "Cd")
		  std_base_atom_names = cytosine;
	       if (std_base_name == "Tr" || std_base_name == "Td")
		  std_base_atom_names = thymine;
	       if (std_base_name == "Ur" || std_base_name == "Ud")
		  std_base_atom_names = uracil;
	    
	       if (std_base_atom_names.size() == 0) {
		  std::cout << "ERROR:: failed to find std_base_name for std_base_atom_names\n";
	       } else {

		  // now find the atoms of the given residue and apply
		  // the transformation and add them to the residue;
		  short int have_deleted = 0;
		  for (unsigned int iat=0; iat<mol_base_atom_names.size(); iat++) {
		     for (int i=0; i<n_mol_base_atoms; i++) {
			if (mol_base_atoms[i]) { 
			   if (mol_base_atom_names[iat] == mol_base_atoms[i]->name) {
// 			      std::cout << ".... Deleting Atom " << mol_base_atoms[i]->name
// 					<< std::endl;
			      residue->DeleteAtom(i);
			      mol_base_atoms[i] = NULL;
			      have_deleted = 1;
			      break;
			   }
			}
		     }
		  }
 		  if (have_deleted)
 		     residue->TrimAtomTable();
		  
		  for (unsigned int iat=0; iat<std_base_atom_names.size(); iat++) {
		     for (int i=0; i<n_std_base_atoms; i++) {
			if (std_base_atom_names[iat] == std_base_atoms[i]->name) {
			   clipper::Coord_orth p(std_base_atoms[i]->x,
						 std_base_atoms[i]->y,
						 std_base_atoms[i]->z);
			   clipper::Coord_orth pt = p.transform(rtop);
			   std::string ele = std_base_atom_names[iat].substr(0,2);
			   at = new CAtom;
// 			   std::cout << ".... Adding Atom " << std_base_atoms[i]->name
// 				     << std::endl;
			   at->SetCoordinates(pt.x(), pt.y(), pt.z(), 1.0, 20.0);
			   at->SetAtomName(std_base_atoms[i]->name);
			   at->SetElementName(ele.c_str());
			   residue->AddAtom(at);
			}
		     }
		  }
	    
	    
		  //       
		  std::string new_base_name = std_base_name;

		  if (mol_base_name.length() != 2) {
		     new_base_name = mol_base_name;
		  } else {
		     // the normal case
		     new_base_name = std_base_name.substr(0,1) +
			mol_base_name.substr(1,2);
		  }
		  
		  strcpy(residue->name, new_base_name.c_str());
		  residue->TrimAtomTable();

	       }
	    }
	 }
      }
   }
} 

clipper::Spacegroup
coot::util::get_spacegroup_from_symops(CMMDBManager *mol) {

   mat44 my_matt;
   clipper::Spacegroup sgp; // initially null.
   clipper::Spgr_descr sgd;
   int ierr = mol->GetTMatrix(my_matt, 0, 0, 0, 0);
   if (ierr == 0) {
      // Good, we have symm info
      std::string symmetry_ops("");
      std::string s = "";
      for ( int i = 0; i < mol->GetNumberOfSymOps(); i++ ) {
	 char* symop_str = mol->GetSymOp(i); // shallow copy, I think
	 s += std::string( symop_str );
	 s += ";";
      }
      short int spacegroup_ok = 1; // set OK initially, then try:
      try {
	 sgd = clipper::Spgr_descr(s, clipper::Spgr_descr::Symops);
	 sgp.init(sgd);
      } catch ( clipper::Message_base exc ) {
	 spacegroup_ok = 0;
      }
   }
   return sgp;
} 

std::vector<std::pair<coot::atom_spec_t, std::string> >
coot::util::gln_asn_b_factor_outliers(CMMDBManager *mol) {

   std::vector<std::pair<coot::atom_spec_t, std::string> > v;
   int imod = 1;
      
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;

   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      CAtom *at;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 int n_atoms = residue_p->GetNumberOfAtoms();
	 std::string residue_name(residue_p->GetResName());
	 if ((residue_name == "ASN") ||
	     (residue_name == "GLN")) {
	    float b_sum = 0.0;
	    float b_sum_sq = 0.0;
	    CAtom *oatom = NULL;
	    CAtom *natom = NULL;
	    CAtom *go_to_atom = NULL;
	    int n_residue_atoms = 0;
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       std::string altloc(at->altLoc);
	       if (altloc == "") {
		  std::string atom_name(at->GetAtomName());
		  
		  if (((residue_name == "GLN") &&
		       ((atom_name == " OE1") ||
			(atom_name == " NE2"))) ||
		      ((residue_name == "ASN") &&
		       ((atom_name == "OD1") ||
			(atom_name == "ND2")))) {
		     
		     if (residue_name == "GLN") { 
			if (atom_name == " OE1") {
			   oatom = at;
			}
			if (atom_name == " NE2") {
			   natom = at;
			}
		     }
		     if (residue_name == "ASN") { 
			if (atom_name == " OD1") {
			   oatom = at;
			}
			if (atom_name == " ND2") {
			   natom = at;
			}
		     }
		  } else {
		     // is a normal atom of the residue:
		     b_sum += at->tempFactor;
		     b_sum_sq += at->tempFactor * at->tempFactor;
		     n_residue_atoms++;
		  }
		  // find the atom to centre on when the button is
		  // clicked (that's why we need an atom spec vector,
		  // not a residue spec.
		  if (residue_name == "GLN") {
		     if (atom_name == " CD ")
			go_to_atom = at;
		  }
		  if (residue_name == "ASN") {
		     if (atom_name == " CG ")
			go_to_atom = at;
		  }
	       }
	    }
	    // end of residue atoms loop
	    if (oatom) {
	       if (natom) {
		  if (n_residue_atoms > 2) {
		     float mean = b_sum/float(n_residue_atoms);
		     float var = b_sum_sq/float(n_residue_atoms) - mean*mean;
		     float std_dev = sqrt(var);
		     float diff = (oatom->tempFactor - natom->tempFactor)/2.0;
		     // we are only interested in cases that have the
		     // O atom B-factor greater than the N atom
		     // B-factor because only they can be fixed by
		     // flipping.
		     if (diff > 0.0) { 
			float z = fabs(diff/std_dev);
			// std::cout << z << std::endl;
			// 		     std::cout << diff << "/sqrt(" << var << ") = " << z
			// 			       << "           " << b_sum << " " << b_sum_sq << " "
			// 			       << n_residue_atoms << std::endl;
			
			if (go_to_atom) { 
			   if (z > 2.25) {
			      coot::atom_spec_t as(go_to_atom);
			      std::string button_label = "Z score = ";
			      button_label += coot::util::float_to_string(z);
			      button_label += "   ";
			      button_label += go_to_atom->GetChainID();
			      button_label += " ";
			      button_label += int_to_string(go_to_atom->GetSeqNum());
			      button_label += " ";
			      button_label += go_to_atom->GetResName();
			      as.float_user_data = z;
			      std::pair<coot::atom_spec_t, std::string> p(as, button_label);
			      v.push_back(p);
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   if (v.size() > 1) {
      // sort v;
      std::sort(v.begin(), v.end(), coot::compare_atom_specs_user_float_in_pair);
   }
   return v;
}

bool
coot::compare_atom_specs_user_float_in_pair(const std::pair<coot::atom_spec_t, std::string> &a,
					    const std::pair<coot::atom_spec_t, std::string> &b) {

   return b.first.float_user_data < a.first.float_user_data ? 1 : 0;
}

bool
coot::util::residue_has_hydrogens_p(CResidue *res) {

   bool result = 0;

   if (res) { 
      PPCAtom residue_atoms;
      int natoms;
      res->GetAtomTable(residue_atoms, natoms);
      for (int iat=0; iat<natoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 std::string ele(at->element);
	 if ((ele == " H") || (ele == " D")) {
	    result = 1;
	    break;
	 }
      }
   }
   return result;
}
