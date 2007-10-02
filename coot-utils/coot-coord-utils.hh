/* coot-utils/coot-coord-utils.hh
 * 
 * Copyright 2006, by The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#ifndef HAVE_COOT_COORD_UTILS_HH
#define HAVE_COOT_COORD_UTILS_HH

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include "mmdb_manager.h"

#include "clipper/core/coords.h"
#include "coot-lsq-types.h"

namespace coot {

   // Perhaps this should be a class function of a class derived from CMMDBManager?
   int write_coords_pdb(CMMDBManager *mol, const std::string &file_name);

   std::string pad_atom_name(const std::string &atom_name_in,
			     const std::string &element);

   class atom_spec_t {
   public:
      std::string chain;
      int resno;
      std::string insertion_code;
      std::string atom_name;
      std::string alt_conf;
      int int_user_data;
      float float_user_data;
      std::string string_user_data;
      atom_spec_t() {}
      atom_spec_t(const std::string &chain_in,
		  int resno_in,
		  const std::string &insertion_code_in,
		  const std::string &atom_name_in,
		  const std::string &alt_conf_in) {
	 chain = chain_in;
	 resno = resno_in;
	 insertion_code = insertion_code_in;
	 atom_name = atom_name_in;
	 alt_conf = alt_conf_in;
      }
      // This presumes at is a member of a coordinate hierarchy.
      atom_spec_t(CAtom *at) {
	 chain = at->GetChainID();
	 resno = at->GetSeqNum();
	 insertion_code = at->GetInsCode();
	 atom_name = at->name;
	 alt_conf = at->altLoc;
      }
      // This presumes at is a member of a coordinate hierarchy.
      atom_spec_t(CAtom *at, const std::string &user_data_string) {
	 chain = at->GetChainID();
	 resno = at->GetSeqNum();
	 insertion_code = at->GetInsCode();
	 atom_name = at->name;
	 alt_conf = at->altLoc;
	 string_user_data = user_data_string;
      }
      void selectatoms(CMMDBManager *mol, int SelHnd) {
	 mol->SelectAtoms(SelHnd, 0, (char *) chain.c_str(),
			  resno, (char *) insertion_code.c_str(),
			  resno, (char *) insertion_code.c_str(),
			  "*",
			  (char *) atom_name.c_str(),
			  "*",
			  (char *) alt_conf.c_str());
      }

      // Presumes that atom can get to SeqNum() and InsCode()? Need
      // tested against a residue not in a hierarchy.
      bool matches_spec(CAtom *atom) const;

      friend std::ostream& operator<< (std::ostream& s, const atom_spec_t &spec);
   };
   
   bool compare_atom_specs_user_float(const coot::atom_spec_t &a1,
				      const coot::atom_spec_t &a2);
   bool compare_atom_specs_user_float_in_pair(const std::pair<coot::atom_spec_t, std::string> &a,
					      const std::pair<coot::atom_spec_t, std::string> &b);
 


   class residue_spec_t {
   public:
      std::string chain;
      int resno;
      std::string insertion_code;
      residue_spec_t(int r) {
	 resno = r;
	 chain = "";
	 insertion_code = "";
      }
      residue_spec_t(const std::string &chain_in, int r) {
	 resno = r;
	 chain = chain_in;
	 insertion_code = "";
      }
      residue_spec_t(const std::string &chain_in, int r,
		     const std::string &ins_code_in) {
	 resno = r;
	 chain = chain_in;
	 insertion_code = ins_code_in;
      }
      residue_spec_t(CResidue *res) {
	 chain = res->GetChainID();
	 resno = res->GetSeqNum();
	 insertion_code = res->GetInsCode();
      } 
      // This one for coot_wrap_guile
      residue_spec_t() {
	 resno = -9999;
	 chain = "";
	 insertion_code = "";
      }
      short int unset_p() const {
	 short int u = 1;
	 if (resno != -9999)
	    u = 0;
	 return u;
      }
   };

   class torsion {
   public:
      // imol torsion pairs
      std::pair<int, atom_spec_t> atom_1;
      std::pair<int, atom_spec_t> atom_2;
      std::pair<int, atom_spec_t> atom_3;
      std::pair<int, atom_spec_t> atom_4;
      torsion(const std::pair<int, atom_spec_t> &atom_1_in,
	      const std::pair<int, atom_spec_t> &atom_2_in,
	      const std::pair<int, atom_spec_t> &atom_3_in,
	      const std::pair<int, atom_spec_t> &atom_4_in) {
	 atom_1 = atom_1_in;
	 atom_2 = atom_2_in;
	 atom_3 = atom_3_in;
	 atom_4 = atom_4_in;
      }
      torsion(int imol,
	      const atom_spec_t &atom_1_in,
	      const atom_spec_t &atom_2_in,
	      const atom_spec_t &atom_3_in,
	      const atom_spec_t &atom_4_in) {
	 atom_1 = std::pair<int, atom_spec_t> (imol, atom_1_in);
	 atom_2 = std::pair<int, atom_spec_t> (imol, atom_2_in);
	 atom_3 = std::pair<int, atom_spec_t> (imol, atom_3_in);
	 atom_4 = std::pair<int, atom_spec_t> (imol, atom_4_in);
      }
      
      // Find 4 atoms in residue that match the torsion spec.  If the
      // returning vector is not of size 4, then this function has
      // failed.
      std::vector<CAtom *> matching_atoms(CResidue *residue);
   };

   class lsq_range_match_info_t {
   public:
      int to_reference_start_resno;
      int to_reference_end_resno;
      int from_matcher_start_resno;
      int from_matcher_end_resno; // is this used? // for validation?
      std::string reference_chain_id;
      std::string matcher_chain_id;
      int match_type_flag; // CA/Main/All
      lsq_range_match_info_t() {};
      lsq_range_match_info_t(int to_reference_start_resno_in,
			     int to_reference_end_resno_in,
			     std::string reference_chain_id_in,
			     int from_matcher_start_resno_in,
			     int from_matcher_end_resno_in,
			     std::string matcher_chain_id_in,
			     short int match_type_flag_in) {
	 match_type_flag = match_type_flag_in;
	 to_reference_start_resno = to_reference_start_resno_in;
	 to_reference_end_resno = to_reference_end_resno_in;
	 from_matcher_start_resno = from_matcher_start_resno_in;
	 from_matcher_end_resno = from_matcher_end_resno_in;
	 reference_chain_id = reference_chain_id_in;
	 matcher_chain_id = matcher_chain_id_in;
      }
      
   };
   
   // return -1 on badness
   int get_selection_handle(CMMDBManager *mol, const atom_spec_t &at);

   //
   double lsq_plane_deviation(const std::vector<clipper::Coord_orth> &v,
			      const clipper::Coord_orth &pt);

   // return 0 or 1
   short int is_main_chain_p(CAtom *at);

   // return 0 or 1
   bool is_hydrogen_p(CAtom *at);

   // return 0 or 1
   short int is_main_chain_or_cb_p(CAtom *at);

   class graph_match_info_t {
   public:
      std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > matching_atom_names;
      bool success;
      // atom_match_names: ((atom_name_ref alt_conf_ref) (atom_name_wrk alt_conf_work))
      std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > atom_match_names; 
      clipper::RTop_orth rtop;
      int n_match;
      double dist_score;
      graph_match_info_t() {
	 n_match = 0;
	 success = 0;
	 dist_score = 0;
      }
   };
   
   // Match on graph
   // 
   // Return the orientation matrix moving res_moving to res_reference
   // and a flag letting us know that the match worked OK.  Also
   // return the vector of matching atom names, which will be used in
   // the ligand analysis (where the rtop is not applied)
   //
   // Sometimes, we don't want to do an automatic rotation/translation
   // of a test ligand onto the reference structure, we want to
   // measure the distances from exactly where they are (e.g. in
   // testing ligand placement).
   //
   // Usually though apply_rtop_flag will be 1.
   // 
   coot::graph_match_info_t graph_match(CResidue *res_moving, CResidue *res_reference, bool apply_rtop_flag);

   namespace util {

      class peptide_torsion_angles_info_t {
      public:
	 short int status;
	 // angles in radians from clipper::Coord_orth::torsion
	 double omega;
	 double phi;
	 double psi;
	 std::string altconf;
      };

      class atom_spec_and_button_info_t {
      public:
	 atom_spec_t as;
	 std::string button_label;
	 std::string callback_func;
	 atom_spec_and_button_info_t(atom_spec_t as_in,
				     std::string button_label_in,
				     std::string callback_func_in) {
	    as = as_in;
	    button_label = button_label_in;
	    callback_func = callback_func_in;
	 }
      };

      class cis_peptide_info_t {
      public:

	 int serial_number;
	 std::string chain_id_1;
	 std::string residue_name_1;
	 int resno_1;
	 std::string ins_code_1;
	 std::string chain_id_2;
	 std::string residue_name_2;
	 int resno_2;
	 std::string ins_code_2;
	 int model_number;
	 float omega_torsion_angle;

	 // normal constructor used by count_cis_peptides():
	 cis_peptide_info_t(const std::string &chain_id,
			    residue_spec_t res1,
			    residue_spec_t res2,
			    int model_number_in) {
	    model_number = model_number_in;
	    serial_number = -1; // unset
	    chain_id_1 = chain_id;
	    chain_id_2 = chain_id;
	    resno_1 = res1.resno;
	    resno_2 = res2.resno;
	    ins_code_1 = res1.insertion_code;
	    ins_code_2 = res2.insertion_code;
	 } 

	 // Full constructor
	 cis_peptide_info_t(int serial_number_in,
			    std::string chain_id_1_in,
			    std::string residue_name_1_in,
			    int resno_1_in,
			    std::string ins_code_1_in,
			    std::string chain_id_2_in,
			    std::string residue_name_2_in,
			    int resno_2_in,
			    std::string ins_code_2_in,
			    int model_number_in,
			    float omega_torsion_angle_in) {

	    serial_number = serial_number_in;
	    chain_id_1 = chain_id_1_in;
	    residue_name_1 = residue_name_1_in;
	    resno_1 = resno_1_in;
	    ins_code_2 = ins_code_1_in;
	    chain_id_2 = chain_id_2_in;
	    residue_name_2 = residue_name_2_in;
	    resno_2 = resno_2_in;
	    ins_code_2 = ins_code_2_in;
	    model_number = model_number_in;
	    omega_torsion_angle = omega_torsion_angle_in;
	 }

	 // Full from mmdb structure
	 cis_peptide_info_t(CCisPep *cis) {
	    serial_number = cis->serNum;
	    chain_id_1 = cis->chainID1;
	    residue_name_1 = cis->pep1;
	    resno_1 = cis->seqNum1;
	    ins_code_1 = cis->icode1;
	    chain_id_2 = cis->chainID2;
	    residue_name_2 = cis->pep2;
	    resno_2 = cis->seqNum2;
	    ins_code_2 = cis->icode2;
	    model_number = cis->modNum;
	    omega_torsion_angle = cis->measure;
	 }

	 bool operator==(const cis_peptide_info_t &a) {
	    bool r = 0;

	    // The model number in pdb files usually bogus because of badly formed CISPEP cards
	    
//  	    std::cout << "comparing "  // << model_number 
//  		      << " :" << chain_id_1 << ": " << resno_1 << " :" << ins_code_1 << ": :"
//  		      <<         chain_id_2 << ": " << resno_2 << " :" << ins_code_2 << ": "
//  		      << "\nto\n"
//  		      << "          " // << a.model_number
//  		      << " :" << a.chain_id_1 << ": " << a.resno_1 << " :" << a.ins_code_1 << ": :"
//  		      <<         a.chain_id_2 << ": " << a.resno_2 << " :" << a.ins_code_2 << ": "
//  		      << std::endl;
	    // if (a.model_number == model_number) { 
	       if (a.chain_id_1 == chain_id_1) { 
		  if (a.chain_id_2 == chain_id_2) {
		     if (a.resno_1 == resno_1) { 
			if (a.resno_2 == resno_2) {
			   if (a.ins_code_1 == ins_code_1) { 
			      if (a.ins_code_2 == ins_code_2) { 
				 r = 1;
			      }
			   }
			}
		     }
		  }
	       }
            // }
	    return r;
	 }
	 
      };
	 
      std::string single_letter_to_3_letter_code(char code);

      short int is_nucleotide(CResidue *r); // by refmac restraint naming
                                            // e.g. Td Ur Gr Ad etc.

      bool nucleotide_is_DNA(CResidue *r);  // test for presence of O2'

      bool residue_has_hydrogens_p(CResidue *res);

      std::vector<std::string> residue_types_in_molecule(CMMDBManager *mol);
      std::vector<std::string> residue_types_in_chain(CChain *chain_p);
      std::vector<std::string> chains_in_molecule(CMMDBManager *mol);
      int number_of_residues_in_molecule(CMMDBManager *mol);
      // Return -1 on badness
      int max_number_of_residues_in_chain(CMMDBManager *mol);
      // Return NULL on no such chain:
      CChain *chain_only_of_type(CMMDBManager *mol, const std::string &residue_type);

      clipper::RTop_orth matrix_convert(mat44 mat);
      
      // Return -1 on badness.
      // 
      // So that we can calculate the length of the graph x axis - there may
      // be gaps, which is why max_number_of_residues_in_chain would fail.
      // 
      int max_min_max_residue_range(CMMDBManager *mol);

      // return first == 0 if residues not found in chain
      std::pair<short int, int> min_resno_in_chain(CChain *chain_p);
      std::pair<short int, int> max_resno_in_chain(CChain *chain_p);

      // Return -1 on badness (actually, number of chains in the last model)
      int number_of_chains(CMMDBManager *mol);

      // The sum of all occupancies:
      float occupancy_sum(PCAtom *atoms, int n_atoms); 

      float median_temperature_factor(PPCAtom atom_selection,
				      int n_atoms,
				      float low_cutoff,
				      float high_cutoff,
				      short int apply_low_cutoff,
				      short int apply_high_cuttoff);
      float average_temperature_factor(PPCAtom atom_selection,
				       int n_atoms,
				       float low_cutoff,
				       float high_cutoff,
				       short int apply_low_cutoff,
				       short int apply_high_cuttoff);

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
      // 
      std::pair<CMMDBManager *, int>
      create_mmdbmanager_from_res_selection(CMMDBManager *orig_mol, 
					    PCResidue *SelResidues, 
					    int nSelResidues, 
					    int have_flanking_residue_at_start,
					    int have_flanking_residue_at_end, 
					    const std::string &altconf,
					    const std::string &chain_id_1,
					    short int residue_from_alt_conf_split_flag);

      // ignore atom index transfer, return NULL on error.
      CMMDBManager *create_mmdbmanager_from_residue(CMMDBManager *orig_mol,
						    CResidue *res);

      CMMDBManager *create_mmdbmanager_from_atom_selection(CMMDBManager *orig_mol,
							   int SelectionHandle);
      
      // utility function for above:
      CResidue* deep_copy_this_residue(const CResidue *residue,
				       const std::string &altconf,
				       short int whole_residue_flag);

      CResidue *copy_and_delete_hydrogens(CResidue *residue_in);
      
      // utility function for above:
      CResidue* deep_copy_this_residue_with_atom_index_and_afix_transfer(CMMDBManager *std_mol, 
									 const CResidue *residue,
									 const std::string &altconf,
									 short int whole_residue_flag,
									 int new_atom_index_udd_handle,
									 int new_afix_handle);

      // return a vector of residue that are in this fragment.
      // Fragments are marked by consecutively numbered residues.  A
      // gap in the sequence numbers marks the end/beginning of a
      // fragment.
      std::vector<PCResidue> get_residues_in_fragment(CChain *clicked_residue_chain_p,
						      residue_spec_t clicked_residue);

      // deleted by calling process
      std::pair<CMMDBManager *, std::vector<coot::residue_spec_t> > 
      get_fragment_from_atom_spec(const coot::atom_spec_t &atom_spec,
						CMMDBManager *mol);


      // transform atoms in residue
      void transform_atoms(CResidue *res, const clipper::RTop_orth &rtop);

      // transform all the atom in mol
      void transform_mol(CMMDBManager *mol, const clipper::RTop_orth &rtop);


      // Rotate position round vector, return a position.
      // 
      clipper::Coord_orth rotate_round_vector(const clipper::Coord_orth &direction,
					      const clipper::Coord_orth &position,
					      const clipper::Coord_orth &origin_shift,
					      double angle);
      

      // A useful function that was (is) in molecule_class_info_t
      //
      // If this residue has a CA, return a pointer to that atoms, if
      // not, return the first atom, if no atoms, return NULL.
      // 
      CAtom* intelligent_this_residue_mmdb_atom(CResidue *res_p);

      std::string three_letter_to_one_letter(const std::string &resname);


      // for sequence/sequence alignment
      //
      // Take into account the insertion code too:
      std::vector<std::pair<CResidue *, int> > sort_residues_by_seqno(PCResidue *residues,
								      int nResidues);
      // Use the results of the above to give us a sequence string:
      std::string model_sequence(const std::vector<std::pair<CResidue *, int> > &sa);
      bool compare_residues(const std::pair<CResidue *, int> &a,
			    const std::pair<CResidue *, int> &b);

      // extents
      std::pair<clipper::Coord_orth, clipper::Coord_orth> extents(CMMDBManager *mol);
      std::pair<clipper::Coord_orth, clipper::Coord_orth> extents(CMMDBManager *mol,
								  int SelectionHandle);
      // pair.second = 0 for failure
      // pair.first  = 1 for success
      //
      // How to convert a standard orientation residue to res position
      // 
      std::pair<clipper::RTop_orth, short int> get_ori_to_this_res(CResidue *res);

      // residues with insertion codes
      std::vector<CResidue *> residues_with_insertion_codes(CMMDBManager *mol); 

      // LSQing
      //
      std::pair<short int, clipper::RTop_orth>
      get_lsq_matrix(CMMDBManager *mol1,
		     CMMDBManager *mol2,
		     const std::vector<coot::lsq_range_match_info_t> &matches);
      // used by above
      // On useful return, first.length == second.length and first.length > 0.
      // 
      std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> > 
      get_matching_indices(CMMDBManager *mol1,
			   CMMDBManager *mol2,
			   int SelHnd1,
			   int SelHnd2,
			   const coot::lsq_range_match_info_t &match);

      // Return the status in the first position
      // Return the angle in radians.
      std::pair<short int, double> omega_torsion(CResidue *C_residue, CResidue *N_residue,
						 const std::string &altconf);
      // and along the same lines:
      // Return the angles in radians.
      peptide_torsion_angles_info_t peptide_torsions(CResidue *C_residue, CResidue *N_residue,
						     const std::string &altconf);

      // Return the RTop that matches moving to reference.  Don't move
      // moving though.
      std::pair<bool,clipper::RTop_orth> base_to_base(CResidue *reference,
						      CResidue *moving);

      // Return the RTop that matches moving to reference.  Don't move
      // moving.  Lite above, but add phosphate and furanose atoms.
      std::pair<bool,clipper::RTop_orth> nucleotide_to_nucleotide(CResidue *reference,
								  CResidue *moving);

      int mutate(CResidue *res, CResidue *std_res_unoriented, short int shelx_flag);

      // given a std residue oriented over residue, make the mutation
      // to std_residue
      void mutate_internal(CResidue *residue, CResidue *std_residue_oriented,
			   short int is_from_shelx_ins_flag);
      
      void mutate_base(CResidue *residue, CResidue *std_base);
      
      // For use with interesting-things-gui, make the list argument
      // from a vector of atom specs.
      // 
      // use the user data in the atom spec to give us the molecule number
      // and the button label
      // 
      std::string interesting_things_list(const std::vector<atom_spec_t> &v);

      // error_type is e.g. "Z score", "Clash gap"
      std::string
      interesting_things_list_with_fix(const std::vector<atom_spec_and_button_info_t> &v,
				       const std::string error_type);

      // 
      clipper::Spacegroup get_spacegroup_from_symops(CMMDBManager *mol);

      // return a set of residue specifiers and a string to put on the
      // button.  Atsushi Nakagawa told me about the problem of these
      // flips and how they could be identified by B factor
      // outlierness.
      std::vector<std::pair<coot::atom_spec_t, std::string> >
      gln_asn_b_factor_outliers(CMMDBManager *mol);

      // return the number of cis peptides in mol:
      int count_cis_peptides(CMMDBManager *mol);

      // more info on the real cis peptides derived from atom positions:
      std::vector<cis_peptide_info_t>
      cis_peptides_info_from_coords(CMMDBManager *mol);

      // remove wrong cis_peptides
      void remove_wrong_cis_peptides(CMMDBManager *mol);
      

   }
   
}

#endif // HAVE_COOT_COORD_UTILS_HH
