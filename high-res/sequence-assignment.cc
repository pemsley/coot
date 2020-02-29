/* high-res/sequence-assignment.cc
 * 
 * Copyright 2003, 2004  by Paul Emsley, The University of York
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
#include <algorithm>
#include <stdio.h>
#include <math.h>

#include "coords/Cartesian.h"
#include "coords/mmdb-extras.h"
#include "coords/mmdb-crystal.h" // for atom_selection_container_t usage

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "ligand/dunbrack.hh"
#include "sequence-assignment.hh"


void
coot::sequence_assignment::scored_chain_info_t::debug() const {

   // Down the page we have the residue number, across the page we
   // have the scores for the different residue types
   //

   std::cout << "Chain: " << chain_name_ << "\n";
   std::cout << " Resno GLY   ALA   SER   VAL   THR   PRO   ASN   ASP   "
	     << "CYS   GLN   GLU   "
	     << "HIS   ILE   LEU   LYS   MET   PHE   TYR   TRP   ARG\n";
   for (unsigned int i=0; i<residue_side_chain_score.size(); i++) {
      printf(" %3d ", i);
      for (unsigned int j=0; j<residue_side_chain_score[i].size(); j++) {
	 printf(" %1.3f", residue_side_chain_score[i][j]);
      }
      std::cout << "\n";
   }
}

void
coot::sequence_assignment::scored_chain_info_t::debug_with_cout() const {

   // Down the page we have the residue number, across the page we
   // have the scores for the different residue types
   //

   std::cout << "Chain: " << chain_name_ << "\n";
   std::cout << " Resno GLY    ALA    SER    VAL     THR    PRO    ASN    ASP    "
	     << "CYS    GLN    GLU "
	     << "HIS    ILE    LEU    LYS    MET    PHE    TYR    TRP    ARG\n";
   for (unsigned int i=0; i<residue_side_chain_score.size(); i++) {
      std::cout << " ";
      std::cout.width(3);
      std::cout << i << " ";
      for (unsigned int j=0; j<residue_side_chain_score[i].size(); j++) {
	 std::cout.width(3);
	 std::cout.precision(3);
	 std::cout << residue_side_chain_score[i][j] << " ";
      }
      std::cout << "\n";
   }

}

std::string
coot::sequence_assignment::side_chain_name_index_to_name(const side_chain_name_index &idx) {

   if (idx == coot::sequence_assignment::GLY) return "GLY";
   if (idx == coot::sequence_assignment::ALA) return "ALA";
   if (idx == coot::sequence_assignment::SER) return "SER";
   if (idx == coot::sequence_assignment::VAL) return "VAL";
   if (idx == coot::sequence_assignment::THR) return "THR";
   if (idx == coot::sequence_assignment::PRO) return "PRO";
   if (idx == coot::sequence_assignment::ASN) return "ASN";
   if (idx == coot::sequence_assignment::ASP) return "ASP";
   if (idx == coot::sequence_assignment::CYS) return "CYS";
   if (idx == coot::sequence_assignment::GLN) return "GLN";
   if (idx == coot::sequence_assignment::GLU) return "GLU";
   if (idx == coot::sequence_assignment::HIS) return "HIS";
   if (idx == coot::sequence_assignment::ILE) return "ILE";
   if (idx == coot::sequence_assignment::LEU) return "LEU";
   if (idx == coot::sequence_assignment::LYS) return "LYS";
   if (idx == coot::sequence_assignment::MET) return "MET";
   if (idx == coot::sequence_assignment::PHE) return "PHE";
   if (idx == coot::sequence_assignment::TYR) return "TYR";
   if (idx == coot::sequence_assignment::TRP) return "TRP";
   if (idx == coot::sequence_assignment::ARG) return "ARG";

   
//    GLY, ALA, SER, VAL, THR, PRO, ASN, ASP, CYS,
//       GLN, GLU, HIS, ILE, LEU, LYS, MET, PHE, TYR,
//       TRP, ARG

   return "";
}



void
coot::sequence_assignment::scored_chain_info_t::add_score(int resno,
							  int sc_name_idx,
							  double score) {

   if (resno > max_resno()) {
      std::cout << "unable to assign to residue " << resno
		<< " in chain " << chain_name_ << std::endl;
   } else {
      if (resno <= int(residue_side_chain_score.size()) ) {
	 if (residue_side_chain_score[resno].size() == 0) {
	    residue_side_chain_score[resno].resize(20);
	 }
	 residue_side_chain_score[resno][sc_name_idx] = score;
      } else { 
	 std::cout << "residue number out of range: " 
		   << resno << " " << residue_side_chain_score.size()
		   << std::endl;
      } 
   } 
}

void
coot::sequence_assignment::side_chain_score_t::add_score(int chain_number,
							 const std::string &chain_id,
							 int residue_number,
							 int max_resno_in_chain,
							 int residue_idx,
							 double score) { 

   if (int(side_chain_score.size()) <= chain_number) { 
      // side_chain_score element needs to be initialized with
      // chain_id and number of residues in the chain.
      side_chain_score.resize(chain_number+1, 
			      scored_chain_info_t(chain_id, 
						  max_resno_in_chain));
   }

   side_chain_score[chain_number].add_score(residue_number, residue_idx, score);
}


void
coot::sequence_assignment::side_chain_score_t::debug() const {

   for (unsigned int i=0; i<side_chain_score.size(); i++) {
      std::cout << "\n";
      side_chain_score[i].debug();
   }
}

// http://ngfnblast.gbf.de/docs/fasta.html
// 
// For those programs that use amino acid query sequences (BLASTP and
// TBLASTN), the accepted amino acid codes are:
// 
//     A  alanine                         P  proline
//     B  aspartate or asparagine         Q  glutamine
//     C  cystine                         R  arginine
//     D  aspartate                       S  serine
//     E  glutamate                       T  threonine
//     F  phenylalanine                   U  selenocysteine
//     G  glycine                         V  valine
//     H  histidine                       W  tryptophane
//     I  isoleucine                      Y  tyrosine
//     K  lysine                          Z  glutamate or glutamine
//     L  leucine                         X  any
//     M  methionine                      *  translation stop
//     N  asparagine                      -  gap of indeterminate length

// sequence
void
coot::sequence_assignment::side_chain_score_t::add_fasta_sequence(const std::string &sequence_chain_id_in, const std::string &seq_in) {

   // format "> name\n <sequence>", we ignore everything that is not a
   // letter after the newline.

   // sequence is member data.  Let's fill it.

   std::string seq;

   int nchars = seq_in.length();
   short int found_greater = 0;
   short int found_newline = 0;
   std::string t;

   for (int i=0; i<nchars; i++) {

      //       std::cout << "checking character: " << seq_in[i] << std::endl;

      if (found_newline && found_greater) {
	 t = toupper(seq_in[i]);
	 if (is_fasta_aa(t)) {
	    // std::cout << "adding character: " << seq_in[i] << std::endl;
	    seq += t;
	 }
      }
      if (seq_in[i] == '>') {
	 // std::cout << "DEBUG:: " << seq_in[i] << " is > (greater than)\n";
	 found_greater = 1;
      }
      if (seq_in[i] == '\n') { 
	 if (found_greater) {
	    // std::cout << "DEBUG:: " << seq_in[i] << " is carriage return\n";
	    found_newline = 1;
	 }
      }
   }
   
   if (seq.length() > 0) { 
      std::cout << "storing sequence: " << seq << " for chain id: "
		<< sequence_chain_id_in << std::endl;
      input_sequence.push_back(std::pair<std::string, std::string> (sequence_chain_id_in,seq));   
      // make sequence_as_indices from input_sequence:
      // sequence_as_indices.push_back(std::pair<std::string, std::vector<int> > (sequence_chain_id_in, convert_seq_to_indices(seq)));
      sequence_infos.push_back(coot::sequence_assignment::sequence_info_t(sequence_chain_id_in, convert_seq_to_indices(seq)));
				    
   } else {
      std::cout << "WARNING:: no sequence found or improper fasta sequence format\n";
   }
}

 
// Note vector can return 1000 on unfound residue type.
// 
std::vector<coot::sequence_assignment::side_chain_name_index>
coot::sequence_assignment::side_chain_score_t::convert_seq_to_indices(const std::string &seq) const { 

   std::vector<coot::sequence_assignment::side_chain_name_index> r;
   std::string w = "WARNING:: The following codes were not comprehensible:\n";
   int incomprehensible = 0;
   
   for (unsigned int i=0; i<seq.length(); i++) {

      std::string s = seq.substr(i,1);
      coot::sequence_assignment::side_chain_name_index t = convert_slc_to_index(s);
      if (is_valid_residue_type(t))
	 r.push_back(t);
      else {
	 incomprehensible++;
	 w += s;
	 w += " at position ";
	 w += coot::util::int_to_string(i);
	 w += "\n";
      }
   }
   if (incomprehensible) {
      // and do a GUI warning?
      std::cout << w << std::endl;
   }
   return r;
}

// Return -1 on unfound residue.
coot::sequence_assignment::side_chain_name_index
coot::sequence_assignment::side_chain_score_t::convert_slc_to_index(const std::string &code) const {

   coot::sequence_assignment::side_chain_name_index i =
      coot::sequence_assignment::side_chain_name_index(1000); // invalid return
   
   if (code == "A") i = coot::sequence_assignment::ALA;
   else 
      if (code == "G") i = coot::sequence_assignment::GLY;
      else
	 if (code == "A") i = coot::sequence_assignment::ALA;
	 else 
	    if (code == "V") i = coot::sequence_assignment::VAL;
	    else
	       if (code == "S") i = coot::sequence_assignment::SER;
	       else 
		  if (code == "N") i = coot::sequence_assignment::ASN;
		  else 
		     if (code == "P") i = coot::sequence_assignment::PRO;
		     else 
			if (code == "D") i = coot::sequence_assignment::ASP;
			else 
			   if (code == "C") i = coot::sequence_assignment::CYS;
			   else 
			      if (code == "Q") i = coot::sequence_assignment::GLN;
			      else 
				 if (code == "E") i = coot::sequence_assignment::GLU;
				 else 
				    if (code == "H") i = coot::sequence_assignment::HIS; 
				    else 
				       if (code == "I") i = coot::sequence_assignment::ILE; 
				       else 
					  if (code == "L") i = coot::sequence_assignment::LEU; 
					  else 
					     if (code == "K") i = coot::sequence_assignment::LYS; 
					     else 
						if (code == "M") i = coot::sequence_assignment::MET; 
						else 
						   if (code == "F") i = coot::sequence_assignment::PHE; 
						   else 
						      if (code == "T") i = coot::sequence_assignment::THR; 
						      else 
							 if (code == "W") i = coot::sequence_assignment::TRP;
							 else 
							    if (code == "Y") i = coot::sequence_assignment::TYR; 
							    else 
							       if (code == "R") i = coot::sequence_assignment::ARG;
   return i;
}


short int
coot::sequence_assignment::side_chain_score_t::is_fasta_aa(const std::string &a) const { 

   short int r = 0;
   
   if (a == "A" || a == "G" ) { 
      r = 1;
   } else { 
      if (a == "B" 
	  || a == "C" || a == "D" || a == "E" || a == "F" || a == "H" || a == "I"
	  || a == "K" || a == "L" || a == "M" || a == "N" || a == "P" || a == "Q" 
	  || a == "R" || a == "S" || a == "T" || a == "U" || a == "V" || a == "W" 
	  || a == "Y" || a == "Z" || a == "X" || a == "*" || a == "-") {
	 r = 1;
      }
   }
   return r;
}


void
coot::sequence_assignment::side_chain_score_t::slider() const { 

   std::cout << "Sliding!\n";
   int s;

   for (unsigned int i=0; i<side_chain_score.size(); i++) { 
      for (unsigned int j=0; j<input_sequence.size(); j++) {
	 // s = side_chain_score[i].slider_hit(sequence_as_indices[j].second);
	 if (j < sequence_infos.size()) { 
	    s = side_chain_score[i].slider_hit(sequence_infos[j].residue_info);
	    if (s != -1) { 
	       std::cout << "We found a hit\n";
	    }
	 } else {
	    std::cout << "ERROR:: Trapped indexing error (slider) " << j << std::endl;
	 }
      }
   }
}

// return -1 on no outstanding slider position
int
coot::sequence_assignment::scored_chain_info_t::outstanding_slider_position(const std::vector<float> &s) const {

   int ir = -1;

   // identify the highest value in s
   // create a new vector r, without this value
   // 
   if (s.size() > 0) { 
      float best_val = -99999999.9;
      int best_val_idx = -1;
      for (unsigned int i=0; i<s.size(); i++) {
	 if (s[i] > best_val) { 
	    best_val = s[i];
	    best_val_idx = i;
	 }
      }
      
      std::vector<float> r;
      for (unsigned int i=0; i<s.size(); i++) {
	 if (int(i) != best_val_idx) { 
	    r.push_back(s[i]);
	 }
      }
      
      if (r.size() > 0) { 
	 float sum = 0;
	 float sum_sq = 0;
	 for (unsigned int i=0; i<r.size(); i++) { 
	    sum += r[i];
	    sum_sq += r[i]*r[i];
	 }
	 // float mean = sum/float(r.size());
	 float std_dev = sqrt(sum_sq/(float(r.size())));

	 float smallest_diff = 9e50;
	 for (unsigned int i=0; i<r.size(); i++) {
	    if ( (best_val - r[i]) < smallest_diff) { 
	       smallest_diff = best_val - r[i];
	    }
	 }

	 if (smallest_diff/std_dev > 3.0) { 
	    std::cout << "Found an outstanding matcher! at index"
		      << best_val_idx << std::endl;
	    ir = best_val_idx;
	 } 
	 
	 
      } else { 
	 ir = best_val_idx;
      }
   }

   return ir;
}

// Generate a set of scores for each offset
// 
// We want to pick the best score in this vector if it is outstanding
// 
int 
coot::sequence_assignment::scored_chain_info_t::slider_hit(const std::vector<std::pair<side_chain_name_index, float> > &seq) const { 

   // slide seq along residue_side_chain_score.  

   // What is the length of the returned object:
   // 
   // This table for fullies:
   //    <- SEQ ->                {length m}
   // RRRRRRRRRRRRRRRRRRRRRR      {length n}
   // 
   // residue_side_chain_score length  seq length    returned obj length
   // 
   //               2                      1              2
   //               4                      2              3
   //               n                      m              n-m+1
   // 
   // This table for partial overhangs (i.e. seq overhangs the end of
   // residue_side_chain_score) (either end)
   // 
   //    <- SSSSSSSSSSSEEEEQQQQQQUEEEENNNNCCCCEEEE ->
   // RRRRRRRRRRRRRRRRRRRRRR
   // 
   // residue_side_chain_score length  seq length    returned obj length
   // 
   //               2                      1              2
   //               4                      2              3
   //               n                      m              n+m-1
   // 
   // offset range:  for SEQ on RRRR
   // 
   // SEQ
   //   RRRR   -2     offset: 1-m
   // 
   //    SEQ
   // RRRR      3     offset: n-1
   
   std::vector<float> r;

   int m = seq.size();
   int n = residue_side_chain_score.size();
   float sum_score;
   int residues_table_index; // index of an index/residue_type score
   int sc_name_idx;

   std::cout << "DEBUG:: residue_side_chain_score has size "
	     << residue_side_chain_score.size() << std::endl;
   for (int offset=1-m; offset<=(n-1); offset++) {
      
      sum_score = 0.0;
      for (int i_in_seq=0; i_in_seq<m; i_in_seq++) { // e.g. m = 1315

	 // We want to pulll out the score of this residue type from
	 // the table residue_side_chain_score. 
	 // 
	 sc_name_idx = seq[i_in_seq].first;
	 residues_table_index = offset + i_in_seq;
	 if (residues_table_index >=0 && residues_table_index < n) {
	    if (int(residue_side_chain_score[residues_table_index].size()) > sc_name_idx) { 
	       sum_score += residue_side_chain_score[residues_table_index][sc_name_idx];
	    } else {
	       std::cout << "ERROR:: Trapped indexing problem (slider_hit): table index "
			 << residues_table_index << " sc idx: " << sc_name_idx
			 << " but size: "
			 << residue_side_chain_score[residues_table_index].size() << std::endl;
	    }
	 } 
      }
      
      r.push_back(sum_score);
   } 

   // The position of this call should be sorted out later.
   // Consider what this function should return to it container class.
   return outstanding_slider_position(r);

} 

void
coot::sequence_assignment::side_chain_score_t::generate_scores(mmdb::Manager *mol_in,
							       const clipper::Xmap<float> &xmap) {

   // How do I approach this, then?
   //
   //
   // We have to find an unassigned portion of the sequence.  How do
   // we do that?  We need to attach to each sequence residue an "is
   // assigned" flag.
   //
   // We find a bit of structure that is unassigned.  How do we do
   // that?  Attach a flag of unassigned-ness as an mmdb UDD.  So, we
   // have a question: how do we mark unassignedness?  Well, we'll
   // pass a argument POLY_ALA and ALL which tells the marking part
   // how to mark the structure it has.
   // 
   // So now we have structure that can be heterogeneous... e.g. there
   // is a stretch of poly ALA (or a separate chain) which will be
   // needing assignment.
   //
   // We need to find such sequence.  This is not easy.
   // find_unassigned_regions() which should return a list of
   // [chain-id + residue range]s.
   //
   // Let's imagine that it does, then what?  We get on to the
   // side_chain_score_t functions...
   //
   // 1) Select an unassigned_region from the list (vector)
   // 
   // 2) Find in the sequence contiguous stretches that are long
   //   enough to cover the length of the unassigned region.  Consider
   //   caching (or pre-analysis of) these results because will will do
   //   the same for the next unassigned region.
   //
   // 3) Now we have an unassigned sequence and an unassigned residue range
   //    so let's start scoring side chains...
   //
   // 4) Simply, for each side chain position, mutate and autofit
   // rotamer and return a score.
   // 

   mol = mol_in;
   mark_unassigned_residues(); // uses mol

   float pr_cut = 0.1; 
   std::vector<residue_range_t> urv = find_unassigned_regions(pr_cut);
   std::cout << "There were " << urv.size() << " unassigned regions\n";
   
      
}


float
coot::sequence_assignment::side_chain_score_t::auto_fit_score(const std::string &chain_id,
							      int resno,
							      const coot::sequence_assignment::side_chain_name_index &idx) {

   float f = 0; // unset
   
   // First find the standard residue matched to the position of resno.
   // Then do an auto-fit rotamer on it - which returns a score.

   // We want 2 major functions:
   //
   // mutate_to_target_type_and_rotamer_score which is given a model
   // mmdb::Residue (ie. what it currently is, the main-chain coordinates
   // of which are used to do the mutation).  We allow this function
   // to have control to run over the rotamers of dunbrack.
   //
   // auto-fit-residue, which we give a mmdb::Residue that has been mutated
   // to the correct type and runs over the dunbrack rotamers and
   // returns the best score (c.f. auto-fit-rotamer in
   // molecule_class_info_t).
   //
   // Between them these functions will fill the table.
   //
   
   return f;
}


float
coot::sequence_assignment::side_chain_score_t::auto_fit_score(mmdb::Residue *current_ala_res,
							      const coot::sequence_assignment::side_chain_name_index &idx,
							      const coot::dictionary_residue_restraints_t &rest,
							      const clipper::Xmap<float> &xmap) {

   // Find reference points in res (CA, N, C) used to make the fit.
   // (this is done in get_ori_to_this_res())
   //
   // Mutate from current residue type to residue type idx
   //
   // Run over rotamers and return the best score
   // (c.f. auto_fit_rotamer in molecule_class_info_t)
   //

   std::map<std::string, clipper::RTop_orth> rtops = 
      coot::util::get_ori_to_this_res(current_ala_res);

   mmdb::Residue *rot_res = get_standard_residue(idx);

   std::map<std::string, clipper::RTop_orth>::const_iterator it = rtops.find("");
   if (it != rtops.end()) { 
      move_std_res_to_this_res_pos(it->second, rot_res); // move rot_res atoms
   }

   // OK, so now we have a residue (rot_res) (which is a
   // rotated/translated standard residue) in the place of
   // current_ala_res.
   //
   // What we want to do now is something like auto-fit-rotamer, but
   // don't move the atoms (except for generating test postions), just
   // return the score of the best rotamer for this residue type
   // (idx).

   return best_rotamer_score(xmap, rest, rot_res);
   
}

// move std_res atoms
void
coot::sequence_assignment::side_chain_score_t::move_std_res_to_this_res_pos(const clipper::RTop_orth &rtop,
									    mmdb::Residue *std_residue) { 
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;

   std_residue->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int iat=0; iat<nResidueAtoms; iat++) {
      clipper::Coord_orth co(residue_atoms[iat]->x,
			     residue_atoms[iat]->y,
			     residue_atoms[iat]->z);
      clipper::Coord_orth rotted = co.transform(rtop);
      residue_atoms[iat]->x = rotted.x();
      residue_atoms[iat]->y = rotted.y();
      residue_atoms[iat]->z = rotted.z();
   }
}


mmdb::Residue *
coot::sequence_assignment::side_chain_score_t::get_standard_residue(const coot::sequence_assignment::side_chain_name_index &idx) const {

   return standard_residues[idx];
}

short int
coot::sequence_assignment::side_chain_score_t::cache_standard_residues() {

   atom_selection_container_t st_res_asc; // = read_standard_residues();
   mmdb::Residue *res;
   int nSelResidues;
   mmdb::PResidue *SelResidue = NULL;
   
   if (st_res_asc.read_success) {
      standard_residues.resize(20);
      for (int i=0; i<20; i++) {
	 int selHnd = st_res_asc.mol->NewSelection();
	 coot::sequence_assignment::side_chain_name_index idx =
	    coot::sequence_assignment::side_chain_name_index(i);
	 std::string residue_type = side_chain_name_index_to_name(idx);
	 st_res_asc.mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
				 "*", // Chain(s) it's "A" in this case.
				 mmdb::ANY_RES,"*",  // starting res
				 mmdb::ANY_RES,"*",  // ending res
				 (char *) residue_type.c_str(),  // residue name
				 "*",  // Residue must contain this atom name?
				 "*",  // Residue must contain this Element?
				 "*",  // altLocs
				 mmdb::SKEY_NEW // selection key
				 );
	 st_res_asc.mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
	 if (nSelResidues < 0) {
	    std::cout << "ERROR:: failed to find standard residue " << residue_type
		      << " in cache_standard_residues" << std::endl;
	 } else {
	    mmdb::Residue *r = coot::util::deep_copy_this_residue(SelResidue[0]);
	    standard_residues[i] = r;
	 }
      }
   }
   return 1; 
}


// Mark up the input structure as unassigned.  Use Residue UDD to do
// the assignment.
// 
std::vector<coot::residue_range_t>
coot::sequence_assignment::side_chain_score_t::find_unassigned_regions(float pr_cut) { 

   std::vector<coot::residue_range_t> v;
   int istate; // for UDD data reading

   int n_models = mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) { 
      
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      if (nchains <= 0) { 
	 std::cout << "bad nchains in molecule " << nchains
		   << std::endl;
      } else {
	 int in_ala_range_flag = 0;
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::PResidue residue_p;
	    mmdb::PResidue stop_residue;
	    mmdb::PResidue previous_residue = NULL;
	    int consecutive_ala_count = 0;
	    std::string chain_id = chain_p->GetChainID();
	    int iassigned;
	    int start_resno = -1; // unset
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int this_resno = residue_p->seqNum;
	       //
	       std::string restype = residue_p->name;
	       //
	       istate = residue_p->GetUDData(udd_assigned_handle, iassigned);
	       if (istate == mmdb::UDDATA_Ok) {
		  //std::cout << "UDdata was OK for " << residue_p->GetChainID() << " "
		  // 			    << residue_p->GetSeqNum() << std::endl;
		  // std::cout << "data: " << iassigned << std::endl;
		  if (iassigned == 1) {
		     if (in_ala_range_flag == 0) {
			in_ala_range_flag = 1;
			start_resno = this_resno;
		     } else {
			// just another residue in the unassigned
			// range.  Do nothing...
		     }
		  } else { // iassigned = 0;
		     // interesting things can happen if we *were* in
		     // a range...
		     if (in_ala_range_flag == 1) {
			if (previous_residue) { 
			   v.push_back(coot::residue_range_t(chain_id,
							     start_resno,
							     previous_residue->seqNum));
			   in_ala_range_flag = 0;
			}
		     }
		     in_ala_range_flag = 0;
		  }
	       } else {
		  std::cout << "ERROR:: can't get UDData!\n";
	       }

	       previous_residue = residue_p; // setup for next run.
	    }
	    // We have got to the end of the chain: were we in a ALA
	    // range? If so, use previous_residue to push back a range info
	    if (in_ala_range_flag)
	       if (previous_residue)
		  v.push_back(coot::residue_range_t(chain_id,
						    start_resno,
						    previous_residue->seqNum));
	 }
      }
   }
   return v;

}


void
coot::sequence_assignment::side_chain_score_t::mark_unassigned_residues() { 
   // Be conservative and mark up ALA ALA ALA (+ ALA *)
   // only

   int istate;
   udd_assigned_handle = mol->RegisterUDInteger (mmdb::UDR_RESIDUE, "assigned residue info");
   if (!udd_assigned_handle) {
      std::cout << "ERROR getting udd_assigned_handle\n";
   } 

   int n_models = mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) { 
      
      mmdb::Model *model_p = mol->GetModel(imod);
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
	       mmdb::PResidue save_res1 = 0; // unset
	       mmdb::PResidue save_res2 = 0; // unset
	       int consecutive_ala_count = 0;
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);

		  //
		  std::string restype = residue_p->name;
		  if (restype == "ALA")
		     consecutive_ala_count++;
		  else
		     consecutive_ala_count = 0;

		  // save the previous 2 residues if are starting a
		  // range:
		  if (consecutive_ala_count > 0) {
		     if (consecutive_ala_count < 3) { 
			if (consecutive_ala_count < 2) { 
			   save_res1 = residue_p;
			} else {
			   save_res2 = residue_p;
			}
		     }
		  }
			     

		  if (consecutive_ala_count > 2) {

		     // mark up this residue and if this is the third, the previous 2
		     //
		     istate = residue_p->PutUDData(udd_assigned_handle, 1);
		     if (istate == mmdb::UDDATA_WrongUDRType)
			std::cout << "ERROR::  mmdb:UDDATA_WrongUDRType in find_unassigned_regions 1" << std::endl;

		     // if this is the 3rd mark the previous 2:
		     if (consecutive_ala_count == 3) {
			istate = save_res1->PutUDData(udd_assigned_handle, 1);
			if (istate == mmdb::UDDATA_WrongUDRType)
			   std::cout << "ERROR::  mmdb:UDDATA_WrongUDRType in find_unassigned_regions 2" << std::endl;
			istate = save_res2->PutUDData(udd_assigned_handle, 1);
			if (istate == mmdb::UDDATA_WrongUDRType)
			   std::cout << "ERROR::  mmdb:UDDATA_WrongUDRType in find_unassigned_regions 3" << std::endl;
		     }
		     
		  } else {

		     // set UDD as not unassigned.
		     istate = residue_p->PutUDData(udd_assigned_handle, 0);
		     if (istate == mmdb::UDDATA_WrongUDRType)
			std::cout << "ERROR::  mmdb:UDDATA_WrongUDRType in find_unassigned_regions 4" << std::endl;
		  }
	       }
	    }
	 }
      }
   }
}


// Often we come here with a whole chain (fragment) (which has a
// nominal chain_id) to be assigned: e.g. residues 1 to 13 of chain G.
// We want to return all unassigned sequence ranges that are long
// enough for 13-1+1 residues.
//
// We make sequence_range_ts from sequence_info_t by using
// residue_range_ts.
// 
// Old notes:
// Find unassigned bits of the sequences that are long enough
// to accomodate a_residue_range.  So if the sequence is 100
// aa and we have already assigned a range 40-50, there are 2
// returned sequences: 1 to 39 and 51 to 100.
// 
std::vector<coot::sequence_assignment::sequence_range_t>
coot::sequence_assignment::side_chain_score_t::find_unassigned_sequence(const coot::residue_range_t &a_residue_range) const {
   
   float v_crit = 0.1;
   unsigned int required_range_size = a_residue_range.length();
   std::vector<coot::sequence_assignment::sequence_range_t> v;

   // well, we often don't know what the real chain is the want to
   // find this residue range, so lets look over all chains looking
   // for approprite stretches.
   //
   for (unsigned int is=0; is<sequence_infos.size(); is++) {
      if (sequence_infos[is].residue_info.size() > required_range_size) {
	 int in_unass_range_flag = 0;
	 int start_resno = -999;
	 int end_resno;
	 int previous_resno = -999;
	 
	 for (unsigned int ir=0; ir<sequence_infos[is].residue_info.size(); ir++) {
	    if (sequence_infos[is].residue_info[ir].second < v_crit) {
	       if (in_unass_range_flag) {
	       } else {
		  start_resno = 1;
		  in_unass_range_flag = 1;
	       }
	    } else {

	       // if we *were* in a range, add that to the returned
	       // vector:
	       if (in_unass_range_flag) {
		  in_unass_range_flag = 0; 
		  end_resno = ir - 1; // is it?  Checkme.
		  if ( (end_resno - start_resno + 1) >= int(required_range_size)) { 
		     v.push_back(coot::sequence_assignment::sequence_range_t(is,
									     start_resno, end_resno));
		  }
	       }
	    }
	    previous_resno = ir;
	 }
	 if (in_unass_range_flag) {
	    end_resno = previous_resno;
	    if ((end_resno - start_resno + 1) >= int(required_range_size)) { 
	       v.push_back(coot::sequence_assignment::sequence_range_t(is,
								       start_resno, end_resno));
	    }
	 }
      }
   }
   return v;
}
	 


// ----
void
coot::sequence_assignment::side_chain_score_t::test_residue_range_marking() {

   std::cout << "in test_residue_range_marking" << std::endl;
   std::vector<coot::residue_range_t> v = find_unassigned_regions(0.1);

   std::cout << "There were " << v.size() << " unassigned regions:\n";
   for (unsigned int i=0; i<v.size(); i++) {
      std::cout << "Range # " << i << " chain " << v[i].chain_id
		<< " " << v[i].start_resno << " "
		<< v[i].end_resno << std::endl;
   }

   coot::residue_range_t rr("A", 20, 30);

   // assign some sequence: that is apply some high proproability to
   // some sequence_infos data.

   unsigned int start_ass = 20;
   unsigned int end_ass = 40;
   for (unsigned int i=start_ass; i<end_ass && i<sequence_infos[0].residue_info.size(); i++) {
      sequence_infos[0].residue_info[i].second = 1.0;
   }
   
   std::vector<coot::sequence_assignment::sequence_range_t> sv =
      find_unassigned_sequence(rr);
   std::cout << "There were " << sv.size() << " available sequence regions:\n";
   for (unsigned int i=0; i<sv.size(); i++) {
      std::cout << i << " chain_id_idx: " << sv[i].chain_id_index << " "
		<< input_sequence[i].first << " "
		<< sv[i].start_sequence_resno << " "
		<< sv[i].end_sequence_resno << std::endl;
   }
}

// best-fit-residue, which we give a mmdb::Residue that has been mutated
// to the correct type and runs over the dunbrack rotamers and
// returns the best score (c.f. auto-fit-rotamer in
// molecule_class_info_t).
float
coot::sequence_assignment::side_chain_score_t::best_rotamer_score(const clipper::Xmap<float> &xmap,
								  const coot::dictionary_residue_restraints_t &rest,
								  mmdb::Residue *res) const {

   float best_score = 0.0;
   float score_this_rotamer;
   
   // c.f. molecule_class_info_t molecule-class-info-other.cc: auto_fit_best_rotamer
   //
   const std::string alt_conf = ""; // needs fixing?
   coot::dunbrack d(res, alt_conf);
   mmdb::Residue *rotamer_res;
   std::vector<float> probabilities = d.probabilities();
   if (probabilities.size() > 0) {
      for (unsigned int i=0; i<probabilities.size(); i++) {
	 std::cout << "--- Rotamer number " << i << " ------"  << std::endl;
	 rotamer_res = d.GetResidue(rest, i);
	 int n_selected_atoms;
	 mmdb::PAtom *residue_atoms;
	 rotamer_res->GetAtomTable(residue_atoms, n_selected_atoms);
	 score_this_rotamer = coot::util::map_score(residue_atoms, n_selected_atoms,
						    xmap, 1);
	 if (score_this_rotamer > best_score)
	    best_score = score_this_rotamer;
      }
   }
   return best_score;
}
