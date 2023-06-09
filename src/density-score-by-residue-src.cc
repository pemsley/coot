/* src/density-score-residue.cc
 * 
 * Copyright 2004 by The University of York
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */



// Portability (getopt) gubbins
#ifndef _MSC_VER
#include <unistd.h> // for getopt(3)
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>


#include <mmdb2/mmdb_math_align.h>
#include <mmdb2/mmdb_tables.h>

#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.hh"

class scored_atom_t {
public:
   std::string atom_name;
   float d;
   scored_atom_t(const std::string &a, float d_in) : atom_name(a) {
      d = d_in;
   }
};

class scored_residue_t {
public:
   std::vector<scored_atom_t> scored_atoms;
   std::string ins_code;
   int resno;
   float d;
   scored_residue_t(int resno_in, const std::string &ins_code_in, float d_in,
		    const std::vector<scored_atom_t> &scored_atoms_in) : scored_atoms(scored_atoms_in), ins_code(ins_code_in) {
      resno = resno_in;
      d = d_in;
   }
   // no scored atoms:
   scored_residue_t(int resno_in, const std::string &ins_code_in, float d_in) : ins_code(ins_code_in) { 
      resno = resno_in;
      d = d_in;
   }
};
   
class scored_chain_t {
public:
   std::string chain_id;
   std::vector<scored_residue_t> scores;
   scored_chain_t(const std::string &chain_id_in, const std::vector<scored_residue_t> &scores_in) :
      chain_id(chain_id_in), scores(scores_in) {
   }
};

// return -111.0 on failure to find level
// 
float estimate_cut_level(const std::vector<scored_chain_t> &density_score_results) {

   float cl = -111.0;
   float sum = 0;
   float sum_sq = 0;
   int n = 0;
   for(unsigned int i=0; i<density_score_results.size(); i++) {
      std::cout << "Density Scores for chain " << density_score_results[i].chain_id << std::endl;
      for(unsigned int ires=0; ires<density_score_results[i].scores.size(); ires++) {
	 sum += density_score_results[i].scores[ires].d;
	 sum_sq += density_score_results[i].scores[ires].d * density_score_results[i].scores[ires].d;
	 n++;
      }
   }
   if (n > 0) {
      float mean = sum/float(n);
      float var = sum_sq/float(n) - mean * mean;
      cl = mean - 2.8 * sqrt(var);
      if (cl > 0) { 
	 std::cout << "INFO:: Auto-Cut level estimated at: " << cl << std::endl;
      } else {
	 std::cout << "INFO:: Auto-Cut level estimated at: 0.0 (reset from " << cl << ")" << std::endl;
	 cl = 0.0;
      } 
   }
   
   return cl;
}



void cut_residues(atom_selection_container_t asc, const std::vector<scored_chain_t> &density_score_results,
		  const std::string &pdb_out_filename, float cut_level) {


   short int made_a_deletion = 0;
   for(unsigned int i=0; i<density_score_results.size(); i++) {
      for(unsigned int ires=0; ires<density_score_results[i].scores.size(); ires++) {
	 if (density_score_results[i].scores[ires].d < cut_level) {
	    std::cout << "Cutting residue " << density_score_results[i].chain_id.c_str()
		      << density_score_results[i].scores[ires].resno << " with density score "
		      << density_score_results[i].scores[ires].d << std::endl;
	    // the -1 here is a real hack!  I don't know why it is needed.
	    asc.mol->DeleteResidue(1, density_score_results[i].chain_id.c_str(),
				   density_score_results[i].scores[ires].resno-1);
	    made_a_deletion = 1;
	 } 
      }
   }
   
   if (made_a_deletion) {
	 asc.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
         asc.mol->FinishStructEdit();
   }
   asc.mol->WritePDBASCII(pdb_out_filename.c_str());
} 

void density_score_molecule(std::string pdb_filename,
			    std::string mtz_filename,
			    std::string pdb_out_filename,
			    std::string f_col,
			    std::string phi_col,
			    std::string w_col,
			    short int use_weights,
			    short int score_atoms_flag,
			    short int do_cut_flag,
			    float cut_level) {


   atom_selection_container_t asc = get_atom_selection(pdb_filename, true, false, false);
   if (asc.n_selected_atoms > 0) {

      std::vector<scored_chain_t> density_score_results;
      clipper::Xmap<float> xmap;
      coot::util::map_fill_from_mtz(&xmap, mtz_filename, f_col, phi_col, w_col,
				    use_weights, 0);
      
      int n_models = asc.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 
      
	 mmdb::Model *model_p = asc.mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in ... " << std::endl;
	       } else { 
		  int nres = chain_p->GetNumberOfResidues();
		  mmdb::PResidue residue_p;
		  std::vector<scored_residue_t> scored_res;
		  
		  // For each residue in chain
		  for (int ires=0; ires<nres; ires++) { 
		     residue_p = chain_p->GetResidue(ires);
		     int nResidueAtoms;
		     mmdb::PPAtom residue_atoms;

		     residue_p->GetAtomTable(residue_atoms, nResidueAtoms);

		     float f = coot::util::map_score(residue_atoms, nResidueAtoms, xmap, 1);
		     f /= coot::util::occupancy_sum(residue_atoms, nResidueAtoms);
		     std::vector<scored_atom_t> sc_at;
		     if (score_atoms_flag) {
			for (int iat=0; iat<nResidueAtoms; iat++) { 
			   float d = coot::util::map_score_atom(residue_atoms[iat], xmap);
			   sc_at.push_back(scored_atom_t(std::string(residue_atoms[iat]->name), d));
			}
		     }
		     scored_residue_t scored_residue(residue_p->GetSeqNum(), residue_p->GetInsCode(), f, sc_at);
		     scored_res.push_back(scored_residue);
		  }
		  density_score_results.push_back(scored_chain_t(chain_p->GetChainID(),
								 scored_res));
	       }
	    }
	 }
      }

      for(unsigned int i=0; i<density_score_results.size(); i++) {
	 std::cout << "Density Scores for chain " <<
	    density_score_results[i].chain_id << std::endl;
	 for(unsigned int ires=0; ires<density_score_results[i].scores.size(); ires++) {
	    std::cout << "Residue " << density_score_results[i].chain_id << " " 
		      << density_score_results[i].scores[ires].resno << " "
		      << density_score_results[i].scores[ires].d << std::endl;
	    if (score_atoms_flag) {
	       for (unsigned int iat=0; iat<density_score_results[i].scores[ires].scored_atoms.size(); iat++) {
		  std::cout << "    atom: "
			    << density_score_results[i].scores[ires].scored_atoms[iat].atom_name
			    << "  "
			    << density_score_results[i].scores[ires].scored_atoms[iat].d
			    << std::endl;
	       }
	    }
	 }
      }

      if (do_cut_flag) {
	 float cl = estimate_cut_level(density_score_results);
	 if (cl > -100) { // didn't fail
	    cut_residues(asc, density_score_results, pdb_out_filename, cl);
	 }
      }
   }
}

int
main(int argc, char **argv) {

   if (argc < 9) { 
      std::cout << "Usage: " << argv[0] 
		<< " --pdbin pdb-in-filename" << " --hklin mtz-filename"
		<< " --f f_col_label"
		<< " --phi phi_col_label"
		<< " --weight weight_col_label"
		<< " --score-atoms"
		<< " --cut"
		<< " --cut-level (absolute) level"
		<< " --pdbout pdb-out-filename"
		<< "\n";
      std::cout << "     where pdbin is the protein (typically)\n";

   } else {

      std::string pdb_file_name;
      std::string pdb_out_filename; 
      std::string  mtz_filename;
      std::string         f_col;
      std::string       phi_col;
      std::string         w_col;
      int n_used_args = 0;
      short int use_weights = 0;
      short int score_atoms_flag = 0;
      short int do_cut_flag = 0; // cut out badly scoring residues
      float cut_level = 0.1; // absolute level below which we (can) cut out residues.

      const char *optstr = "i:h:f:p:w";
      struct option long_options[] = {
	 {"pdbin", 1, 0, 0},
	 {"hklin", 1, 0, 0},
	 {"f", 1, 0, 0},
	 {"phi", 1, 0, 0},
	 {"weight", 1, 0, 0},
	 {"score-atoms", 0, 0, 0},
	 {"cut", 0, 0, 0},
	 {"cut-level", 1, 0, 0},
	 {"pdbout", 1, 0, 0},
	 {0, 0, 0, 0}
      };

      int ch;
      int option_index = 0;
      while ( -1 != 
	      (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

	 // std::cout << "DEBUG:: " << option_index << " " << coot_optarg << std::endl;

	 switch(ch) { 
	    
	 case 0:
	    if (coot_optarg) { 
	       std::string arg_str = long_options[option_index].name;

	       if (arg_str == "pdbin") { 
		  pdb_file_name = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "pdbout") { 
		  pdb_out_filename = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "hklin") { 
		  mtz_filename = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "f") { 
		  f_col = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "phi") {
		  phi_col = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "weight") {
		  w_col = coot_optarg;
		  n_used_args += 2;
		  use_weights = 1; // yes, do use them
	       }
	       if (arg_str == "cut-level") {
		  n_used_args += 2;
		  cut_level = atof(coot_optarg);
	       } 
	       
	       
	    } else { 
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "score-atoms") {
		  std::cout << "scoring individual atoms\n";
		  score_atoms_flag = 1;
		  n_used_args++;
	       } else {
		  if (arg_str == "cut") {
		     std::cout << "ENABLED:: cutting bad residues\n";
		     do_cut_flag = 1;
		  } else { 
		     std::cout << "Malformed option: "
			       << long_options[option_index].name << std::endl;
		  }
	       }
	    }
	    break;

	 case 'i':
	    pdb_file_name = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 case 'h':
	    mtz_filename = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 case 'f':
	    f_col = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 case 'p':
	    phi_col = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 default:
	    std::cout << "default coot_optarg: " << coot_optarg << std::endl;
	    break;
	 }
      }

      // short int do_it = 0;
      if (pdb_file_name.length() == 0) { 
	 std::cout << "Missing input PDB file\n";
	 exit(1);
      } else { 
	 if (mtz_filename.length() == 0) { 
	    std::cout << "Missing MTZ file\n";
	    exit(1);
	 } else { 
	    if (f_col.length() == 0) { 
	       std::cout << "Missing F column name\n";
	       exit(1);
	    } else {
	       if (phi_col.length() == 0) { 
		  std::cout << "Missing PHI column name\n";
		  exit(1);
	       } else { 
// 		  std::cout << "argc: " << argc << " n_used_args: " 
// 			    << n_used_args << std::endl;
		  density_score_molecule(pdb_file_name,
					 mtz_filename,
					 pdb_out_filename, // often dummy
					 f_col,
					 phi_col,
					 w_col,
					 use_weights,
					 score_atoms_flag,
					 do_cut_flag, cut_level);
	       }
	    }
	 }
      }
   }
   return 0;
}
