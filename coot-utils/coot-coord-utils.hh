/* coot-utils/coot-coord-utils.hh
 * 
 * Copyright 2006, 2007, by The University of York
 * Copyright 2008, 2009, 2010 by The University of Oxford
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

#ifndef HAVE_COOT_COORD_UTILS_HH
#define HAVE_COOT_COORD_UTILS_HH

#include <stdexcept>

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

// How should I do this better?
#define CXX_UNSET_CHARGE -99.8

namespace coot {

   // a generally useful class to be used with std::map where the
   // return type is an index (int).  We use this to determine if the
   // map returns a useful result (if the map used a starndard int, it
   // would return a int (typically 0) when the map failed (no indexer
   // in map).
   // 
   class map_index_t {
      int index_;
   public:
      enum index_type { UNASSIGNED = -1 };
      map_index_t() { index_ = UNASSIGNED; }
      map_index_t(int i) { index_ = i; }
      int index() const { return index_; }
      bool is_assigned() { return (index_ != UNASSIGNED); }
      bool operator==(const map_index_t &ti) const {
	 return (ti.index() == index_);
      }
   };


   // for canonical_base_name();
   enum base_t { RNA, DNA };

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
      int model_number;
      atom_spec_t() {
	 chain = "unset";
	 resno = MinInt4;
	 insertion_code = "";
	 model_number = -1;
	 int_user_data = -1;
      }
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
	 model_number = 1; 
      }
      // This presumes at is a member of a coordinate hierarchy.
      atom_spec_t(CAtom *at) {
	 chain = at->GetChainID();
	 resno = at->GetSeqNum();
	 insertion_code = at->GetInsCode();
	 model_number = at->GetModelNum();
	 atom_name = at->name;
	 alt_conf = at->altLoc;
	 int_user_data = -1; // mark as "unset" (better than not setting it)
      }
      // This presumes at is a member of a coordinate hierarchy.
      atom_spec_t(CAtom *at, const std::string &user_data_string) {
	 model_number = at->GetModelNum();
	 chain = at->GetChainID();
	 resno = at->GetSeqNum();
	 insertion_code = at->GetInsCode();
	 atom_name = at->name;
	 alt_conf = at->altLoc;
	 string_user_data = user_data_string;
      }
      void selectatoms(CMMDBManager *mol, int SelHnd) {
	 const char *chainid = chain.c_str();
	 const char *inscode = insertion_code.c_str();
	 const char *atname  = atom_name.c_str(); // atom name
	 const char *altconf = alt_conf.c_str();
	 
	 mol->SelectAtoms(SelHnd, 0, chainid, resno, inscode, resno, inscode,
			  "*", atname, "*", altconf);
      }

      // Presumes that atom can get to SeqNum() and InsCode()? Need
      // tested against a residue not in a hierarchy.
      bool matches_spec(CAtom *atom) const;

      bool operator==(const atom_spec_t &matcher) const {
	 if (matcher.model_number == model_number) { 
	    if (matcher.chain == chain) {
	       if (matcher.resno == resno) {
		  if (matcher.insertion_code == insertion_code) {
		     if (matcher.atom_name == atom_name) {
			if (matcher.alt_conf == alt_conf) {
			   return 1; 
			}
		     }
		  }
	       }
	    }
	 }
	 return 0;
      }

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
      residue_spec_t(const atom_spec_t &atom_spec) { 
         chain = atom_spec.chain;
         resno = atom_spec.resno;
         insertion_code = atom_spec.insertion_code;
      }
      // This one for coot_wrap_guile
      residue_spec_t() {
	 resno = MinInt4;
	 chain = "";
	 insertion_code = "";
      }
      bool unset_p() const {
	 short int u = 1;
	 if (resno != MinInt4)
	    u = 0;
	 return u;
      }
      bool operator==(const residue_spec_t &matcher) const {
	 if (matcher.chain == chain) {
	    if (matcher.resno == resno) {
	       if (matcher.insertion_code == insertion_code) {
		  return 1; 
	       }
	    }
	 }
	 return 0;
      }
      bool operator<(const residue_spec_t &matcher) const{
	 if (matcher.chain == chain) {
	    if (matcher.resno == resno) {
	       if (matcher.insertion_code == insertion_code) {
		  return 0; 
	       } else {
		  if (matcher.insertion_code < insertion_code)
		     return 0;
		  else
		     return 1;
	       }
	    } else {
	       if (matcher.resno < resno)
		  return 0;
	       else
		  return 1;
	    } 
	 } else {
	    if (matcher.chain < chain)
	       return 0;
	    else
	       return 1;
	 } 
	 return 0;
      } 

      friend std::ostream& operator<< (std::ostream& s, const residue_spec_t &spec);
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


   // Return dist in Angstroms, can throw an exception if any of the
   // atoms is null.
   // 
   double distance(CAtom *at_1, CAtom *at_2);
   
   // Return angle in degrees, can throw an exception if any of the
   // atoms is null.
   // 
   double angle(CAtom *at_1, CAtom *at_2, CAtom *at_3);

   class lsq_range_match_info_t {
   public:
      int to_reference_start_resno;
      int to_reference_end_resno;
      int from_matcher_start_resno;
      int from_matcher_end_resno; // is this used? // for validation?
      std::string reference_chain_id;
      std::string matcher_chain_id;
      int match_type_flag; // CA/Main/All
      bool is_single_atom_match;
      std::string reference_atom_name;
      std::string reference_alt_conf;
      std::string matcher_atom_name;
      std::string matcher_alt_conf;
      lsq_range_match_info_t() { is_single_atom_match = 0;};
      lsq_range_match_info_t(int to_reference_start_resno_in,
			     int to_reference_end_resno_in,
			     std::string reference_chain_id_in,
			     int from_matcher_start_resno_in,
			     int from_matcher_end_resno_in,
			     std::string matcher_chain_id_in,
			     short int match_type_flag_in) {
	 is_single_atom_match = 0;
	 match_type_flag = match_type_flag_in;
	 to_reference_start_resno = to_reference_start_resno_in;
	 to_reference_end_resno = to_reference_end_resno_in;
	 from_matcher_start_resno = from_matcher_start_resno_in;
	 from_matcher_end_resno = from_matcher_end_resno_in;
	 reference_chain_id = reference_chain_id_in;
	 matcher_chain_id = matcher_chain_id_in;
      }
      
      lsq_range_match_info_t(std::string reference_chain_id_in,
			     int reference_resno_in,
			     std::string reference_insertion_code_in,
			     std::string reference_atom_name_in,
			     std::string reference_alt_conf_in,
			     std::string matcher_chain_id_in,
			     int matcher_resno_in,
			     std::string matcher_insertion_code_in,
			     std::string matcher_atom_name_in,
			     std::string matcher_alt_conf_in) {
	 is_single_atom_match = 1;
	 match_type_flag = COOT_LSQ_ALL;
	 to_reference_start_resno = reference_resno_in;
	 to_reference_end_resno   = reference_resno_in;
	 from_matcher_start_resno = matcher_resno_in;
	 from_matcher_end_resno   = matcher_resno_in;
	 reference_chain_id       = reference_chain_id_in;
 	 matcher_chain_id         = matcher_chain_id_in;
	 reference_atom_name      = reference_atom_name_in;
	 reference_alt_conf       = reference_alt_conf_in;
	 matcher_atom_name        = matcher_atom_name_in;
	 matcher_alt_conf         = matcher_alt_conf_in;
      }
      
   };
   
   bool sort_chains_util(const std::pair<CChain *, std::string> &a,
			 const std::pair<CChain *, std::string> &b);

   // return -1 on badness
   int get_selection_handle(CMMDBManager *mol, const atom_spec_t &at);

   // return the lsq deviation of the pt atom and (in second) the rms
   // deviation of the atoms in the plane (not the including pt of
   // course)
   std::pair<double, double>
   lsq_plane_deviation(const std::vector<clipper::Coord_orth> &v,
		       const clipper::Coord_orth &pt);

   // return 0 or 1
   bool is_main_chain_p(CAtom *at);

   // return 0 or 1
   bool is_hydrogen_p(CAtom *at);

   // return 0 or 1
   bool is_main_chain_or_cb_p(CAtom *at);

   bool is_member_p(const std::vector<CResidue *> &v, CResidue *a);


   // Throw an exception if there is no consistent seg id for the
   // atoms in the given residue.
   std::string residue_atoms_segid(CResidue *residue_p);

   // Throw an exception if there is no consistent seg id for the
   // atoms in the given chain.
   std::string chain_atoms_segid(CChain *chain_p);

   // Use the above function the get the segid and insert it into all
   // the atoms of receiver.
   bool copy_segid(CResidue *provider, CResidue *receiver);

   // residues are in increasing number order (if they are not, this
   // causes the residue selection of mmdb to fail at the misplaced
   // residue).
   //
   bool residues_in_order_p(CChain *chain_p);

   // return success status as first element
   // 
   std::pair<bool, clipper::Coord_orth> centre_of_molecule(CMMDBManager *mol);

   // convert atoms in residue to HETATMs if residue is not of
   // standard PDB type
   int hetify_residue_atoms_as_needed(CResidue *res);

   // convert atoms in residue to HETATMs.  Return the number of HET
   // atoms.
   int hetify_residue_atoms(CResidue *res);

   // return residue specs for residues that have atoms that are
   // closer than radius Angstroems to any atom in the residue
   // specified by res_in.
   // 
   std::vector<residue_spec_t> residues_near_residue(const residue_spec_t &res_in,
						     CMMDBManager *mol,
						     float radius);

   std::vector<CResidue *> residues_near_residue(CResidue *res_ref, CMMDBManager *mol,
						 float radius);

   std::vector<CResidue *> residues_near_position(const clipper::Coord_orth &pt,
						  CMMDBManager *mol,
						  double radius);

   // a trivial class to hold the residue and the solvent exposure,
   // including and not including the ligand.
   class solvent_exposure_difference_helper_t {
   public:
      residue_spec_t res_spec;
      double exposure_fraction_holo;
      double exposure_fraction_apo;
      solvent_exposure_difference_helper_t(residue_spec_t res_spec_in, double h, double a) {
	 res_spec = res_spec_in;
	 exposure_fraction_holo = h;
	 exposure_fraction_apo  = a;
      }
   };


   // Don't include residues that are HOH residues that are not bonded
   // to the protein (if bonding to res_ref but not protein, then
   // reject.
   //
   std::vector<CResidue *> filter_residues_by_solvent_contact(CResidue *res_ref,
							      CMMDBManager *mol,
							      const std::vector<CResidue *> &residues,
							      const double &water_dist_max);

   // filter out residues that are HOH and are not in contact with
   // res_ref (i.e. have an atom closer than water_dist_max)
   // 
   // std::vector<CResidue *>
   // filter_residues_by_solvent_contact_to_residue(CResidue *res_ref,
   // const std::vector<CResidue *> &residues,
   // const double &water_dist_max);
   //
   // we don't do this here - we'll do it in the viewer - were we can
   // toggle on an off the waters (if they are not rejected at an
   // early stage).
   

   // Return a pair, the bool of which is set if the float is sensible.
   // 
   std::pair<bool,float> closest_approach(CMMDBManager *mol,
					  CResidue *r1, CResidue *r2);

   CResidue *nearest_residue_by_sequence(CMMDBManager *mol,
					 const residue_spec_t &spec);

  // create a new molecule.  rtop_frac comes from the symmetry
  // operator.  If pre_shift_abc is not of size 3 then don't apply it
  // (if it is, do so, of course).
  CMMDBManager *mol_by_symmetry(CMMDBManager *mol, 
				clipper::Cell cell,
				clipper::RTop_frac rtop_frac,
				std::vector<int> pre_shift_abc);

   class close_residues_from_different_molecules_t {
      // Interacting Residues: Return all residues of mol1, mol2 that
      // have atoms that are closer that dist to atoms of mol2/mol1.
      //
      CMMDBManager *combined_mol;
   public:
      close_residues_from_different_molecules_t() {};
      std::pair<std::vector<CResidue *>, std::vector<CResidue *> >
	 close_residues(CMMDBManager *mol1, CMMDBManager *mol2, float dist);
      void clean_up() { delete combined_mol; } 
   };


   // Fiddle with mol. 
   // 
   // sort chains in lexographical order
   void sort_chains(CMMDBManager *mol);


   // Pukka puckers?
   //
   // Throw an exception if it is not possible to generate pucker info
   // 
   class pucker_analysis_info_t {
      enum PUCKERED_ATOM_T { NONE, C1_PRIME, C2_PRIME, C3_PRIME, C4_PRIME, O4_PRIME};
      PUCKERED_ATOM_T puckered_atom_;
      std::string altconf;
      void assign_base_atom_coords(CResidue *residue_p);
      CAtom *N1_or_9;
      CAtom *C1_prime; 
   public:
      float plane_distortion;
      float out_of_plane_distance;
      std::vector<clipper::Coord_orth> ribose_atoms_coords;
      std::vector<clipper::Coord_orth> base_atoms_coords; // the perpendicular distance of the
                                                          // following phosphate is to the *BASE*
                                                          // atoms (stupid boy).

      // Throw an exception if it is not possible to generate pucker info
      // 
      pucker_analysis_info_t(CResidue *res, std::string altconf_in);
      
      // Use the 3' phosphate of the following residue to calculate
      // its out of plane distance.  Decide from that if this should
      // have been 3' or 2'.  Check vs the actual puckering.
      // 
      float phosphate_distance(CResidue *following_res);
      float phosphate_distance_to_base_plane(CResidue *following_res);
      std::string puckered_atom() const;
   };

   

   class graph_match_info_t {
   public:
      // atom_match_names: ((atom_name_wrk alt_conf_work) (atom_name_ref alt_conf_ref))
      std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > matching_atom_names;
      bool success;
      clipper::RTop_orth rtop;
      int n_match;
      double dist_score;
      graph_match_info_t() {
	 n_match = 0;
	 success = 0;
	 dist_score = 0;
      }
      // Change the names in res_moving_names to match those in
      // res_reference as much as possible.  When there is a name
      // collision (i.e. the name maped from the res_reference is
      // already in the res_moving_names (and that is not to be
      // replace by anything)), invent a new name for the atom.
      // Use the internal matching_atom_names.
      void match_names(CResidue *res_moving_names);
   private:
      // Find a new name for name_in that is not already in the residue
      // 
      std::string invent_new_name(const std::string &name_in,
				  const std::string &ele,
				  const std::vector<std::string> &residue_atom_name) const;
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
   graph_match_info_t graph_match(CResidue *res_moving, CResidue *res_reference, bool apply_rtop_flag);


   //
   bool mol_has_symmetry(CMMDBManager *mol);

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
			    int model_number_in,
			    float tors_in) {
	    model_number = model_number_in;
	    serial_number = -1; // unset
	    chain_id_1 = chain_id;
	    chain_id_2 = chain_id;
	    resno_1 = res1.resno;
	    resno_2 = res2.resno;
	    ins_code_1 = res1.insertion_code;
	    ins_code_2 = res2.insertion_code;
	    omega_torsion_angle = tors_in;
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

#ifdef HAVE_MMDB_WITH_CISPEP
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
#endif // HAVE_MMDB_WITH_CISPEP	 

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

      class quaternion {
      public:
	 float q0, q1, q2, q3;
	 quaternion(const float &q0in, const float &q1in,
		    const float &q2in, const float &q3in) {
	    q0 = q0in;
	    q1 = q1in;
	    q2 = q2in;
	    q3 = q3in;
	 } 
	 quaternion(const clipper::Mat33<double> &mat_in);
	 clipper::Mat33<double> matrix() const;
	 float convert_sign(const float &x, const float &y) const;
	 friend std::ostream&  operator<<(std::ostream&  s, const quaternion &q);
	 friend std::ofstream& operator<<(std::ofstream& s, const quaternion &q);

	 static void test_quaternion(); // test yourself
	 static bool close_float_p (const float &f1, const float &f2) { //testing func
	    float d = fabsf(f1-f2);
	    if (d < 0.001)
	       return 1;
	    else
	       return 0;
	 }
	 bool is_similar_p(const quaternion &q) {
	    bool r = 0;
	    if (close_float_p(q.q0, q0) &&
		close_float_p(q.q1, q1) && 
		close_float_p(q.q2, q2) && 
		close_float_p(q.q3, q3)) {
	       r = 1;
	    }
	    return r;
	 }
      };

      class chain_id_residue_vec_helper_t { 
      public: 
	 std::vector<CResidue *> residues;
	 std::string chain_id;
	 void sort_residues();
	 static bool residues_sort_func(CResidue *first, CResidue *second);
	 bool operator<(const chain_id_residue_vec_helper_t &c) const;
      };
	 
      std::string single_letter_to_3_letter_code(char code);

      short int is_nucleotide(CResidue *r); // by refmac restraint naming
                                            // e.g. Td Ur Gr Ad etc.

      bool nucleotide_is_DNA(CResidue *r);  // test for presence of O2'

      bool residue_has_hydrogens_p(CResidue *res);

      // Return NULL on residue not found in this molecule.
      // 
      CResidue *get_residue(const std::string &chain_id, int res_no,
			    const std::string &insertion_code,
			    CMMDBManager *mol);

      CResidue *get_first_residue(CMMDBManager *mol);

      // convenience interface to above
      CResidue *get_residue(const residue_spec_t &rs, CMMDBManager *mol);

      // Return NULL on residue not found in this molecule.
      // 
      CResidue *get_following_residue(const residue_spec_t &rs, 
				      CMMDBManager *mol);

      std::pair<bool, clipper::Coord_orth> get_residue_centre(CResidue *res);
      
      std::vector<std::string> get_residue_alt_confs(CResidue *res);


      std::vector<std::string> residue_types_in_molecule(CMMDBManager *mol);
      // non-standard means not one of the standard protein amino acid
      // residues.
      std::vector<std::string> non_standard_residue_types_in_molecule(CMMDBManager *mol);
      // a utility function 
      std::vector<std::string> standard_residue_types();

      std::vector<std::string> PDB_standard_residue_types();
      
      std::vector<std::string> residue_types_in_chain(CChain *chain_p);
      std::vector<std::string> residue_types_in_residue_vec(const std::vector<CResidue *> &residues);

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
      std::pair<bool, int> min_resno_in_chain(CChain *chain_p);
      std::pair<bool, int> max_resno_in_chain(CChain *chain_p);
      std::pair<bool, int> max_resno_in_molecule(CMMDBManager *mol);
      

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
      float standard_deviation_temperature_factor(PPCAtom atom_selection,
						  int n_atoms,
						  float low_cutoff,
						  float high_cutoff,
						  short int apply_low_cutoff,
						  short int apply_high_cuttoff);

      class contact_atoms_info_t {
      public:
	 enum ele_index_t { ELE_UNASSIGNED, ELE_NA, ELE_K, ELE_MG2, ELE_LI, ELE_CA2 };
	 class contact_atom_t {
	 public:
	    double dist;
	    CAtom *at;
	    contact_atom_t(CAtom *contactor, CAtom *central_atom);
	 };
      private:
	 std::vector<contact_atom_t> contact_atoms;
	 CAtom *at;
	 ele_index_t metal_type;
      public:
	 contact_atoms_info_t() {
	    at = NULL;
	    metal_type = ELE_UNASSIGNED;
	 }
	 contact_atoms_info_t(CAtom *at_central_in, const contact_atom_t &con_at) {
	    at = at_central_in;
	    contact_atoms.push_back(con_at);
	 }
	 unsigned int size() const {
	    return contact_atoms.size();
	 }
	 contact_atom_t &operator[](int i) {
	    return contact_atoms[i];
	 }
	 bool matches_atom(CAtom *at_in) {
	    return at_in == at;
	 }
	 void add(const contact_atom_t &con_at) {
	    contact_atoms.push_back(con_at);
	 }
	 bool has_contacts(int n_contacts, double dist_max) const {
	    bool r = 0; 
	    if (size() >= n_contacts) {
	       int n_local = 0;
	       for (unsigned int i=0; i<contact_atoms.size(); i++) {
		  if (contact_atoms[i].dist <= dist_max) {
		     n_local++;
		  }
	       }
	       if (n_local >= n_contacts)
		  r = 1;
	    }
	    return r;
	 }
	 CAtom *central_atom() const { return at; }
	 double smallest_contact_dist() const {
	    if (contact_atoms.size() == 0)
	       throw std::runtime_error("zero contacts");
	    double d = 999999999999.9;
	    for (unsigned int i=0; i<contact_atoms.size(); i++) {
	       if (contact_atoms[i].dist < d)
		  d = contact_atoms[i].dist;
	    }
	    return d;
	 }
	 bool test_for_na() const;
	 bool test_for_ele(ele_index_t ele_index) const;
	 static std::string ele_to_string(ele_index_t ele) {
	    std::string r = "Unknown";
	    if (ele == ELE_UNASSIGNED)
	       r = "Unassigned";
	    if (ele == ELE_NA)
	       r = "Sodium (Na+)";
	    if (ele == ELE_K)
	       r = "Potassium (K+)";
	    if (ele == ELE_LI)
	       r = "Lithium (Li+)";
	    if (ele == ELE_MG2)
	       r = "Magnesium (Mg+2)";
	    if (ele == ELE_CA2)
	       r = "Calcium (Ca+2)";
	    return r;
	 }
      };

      // This does nieve symmetry expansion, if the mol has symmetry
      // then this class should be used with (typically) a copy of the
      // molecule that has been moved as close as possible to the
      // origin.
      // 
      class water_coordination_t {
	 std::vector<contact_atoms_info_t> atom_contacts;
	 void add_contact(CAtom *atom_contactor, CAtom *atom_central);
	 
	 void add_contacts(CMMDBManager *mol,
			   PCAtom *water_selection, int n_water_atoms, 
			   PCAtom *atom_selection, int n_selected_atoms,
			   realtype min_dist, realtype max_dist,
			   const mat44 &my_mat);
	 void sort_contacts(std::vector<contact_atoms_info_t> *v) const;
	 static bool sort_contacts_func(const contact_atoms_info_t &first,
					const contact_atoms_info_t &second);
      public:
	 water_coordination_t(CMMDBManager *mol, realtype radius_limit);
	 water_coordination_t() {}; 

	 // I don't know what the "get" interface to
	 // water_coordination_t should be. So I'll make one up:
	 //
	 std::vector<contact_atoms_info_t>
	 get_highly_coordinated_waters(int n_contacts,  // at least n_contacts
				       double dist_max) const; // within dist_max

	 // This checks against build-in values from the literature
	 //
	 std::vector<std::pair<coot::util::contact_atoms_info_t, coot::util::contact_atoms_info_t::ele_index_t> > metals() const;
	 
      };

      // a simple full copy, caller deletes.
      CMMDBManager *copy_molecule(CMMDBManager *mol);

      // copy cell, symm, origin and scale cards from m1 to m2 (if possible)
      // return success status.
      bool copy_cell_and_symm_headers(CMMDBManager *m1, CMMDBManager *m2);
      

      // The flanking residues (if any) are in the residue selection (SelResidues).
      // The flags are not needed now we have made adjustments in the calling
      // function.
      // 
      // create_mmdbmanager_from_res_selection must make adjustments
      // 
      // Note: there is now a molecule-class-info version of this - perhaps
      // we should call it?  Next bug fix here: move over to the function call.
      // 
      // We pass the original molecule because here we do atom index transfer.
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

      // We don't mess with the chain ids (give as we get), but also
      // return the handle for the atom index transfer.
      std::pair<CMMDBManager *, int> create_mmdbmanager_from_mmdbmanager(CMMDBManager *);


     std::pair<bool, CMMDBManager *>
     create_mmdbmanager_from_residue_vector(const std::vector<CResidue *> &res_vec);

      // ignore atom index transfer, return NULL on error.
      CMMDBManager *create_mmdbmanager_from_residue(CMMDBManager *orig_mol,
						    CResidue *res);

      CMMDBManager *create_mmdbmanager_from_atom_selection(CMMDBManager *orig_mol,
							   int SelectionHandle,
							   bool invert_selection_flag=0);
      // uses the following:
      CMMDBManager *create_mmdbmanager_from_atom_selection_straight(CMMDBManager *orig_mol,
								    int SelectionHandle);
      // Beware: This destroys (inverts) the atom selection as passed.
      CMMDBManager *create_mmdbmanager_from_inverted_atom_selection(CMMDBManager *orig_mol,
								    int SelectionHandle);

      CMMDBManager *create_mmdbmanager_from_atom(CAtom *at);

      // a new residue for each point.  Caller deletes.
      // 
      CMMDBManager *create_mmdbmanager_from_points(const std::vector<clipper::Coord_orth> &pts);

      // calling function deletes
      // 
      CMMDBManager *create_mmdbmanager_from_residue_specs(std::vector<coot::residue_spec_t> &r1,
							  CMMDBManager *mol);

      void add_copy_of_atom(CMMDBManager *mol, CAtom *atom);

      // return success status, 1 is good, 0 is fail.  Use clipper::Coord_orth constructor
      // 
      bool add_atom(CResidue *res,
		    const std::string &atom_name_1,
		    const std::string &atom_name_2,
		    const std::string &atom_name_3,
		    const std::string &alt_conf, 
		    double length,
		    double angle, // degrees
		    double torsion,
		    const std::string &new_atom_name,
		    const std::string &new_atom_ele,
		    float new_atom_occ,
		    float new_atom_b_factor); // degrees
      
      // utility function for above:
      CResidue* deep_copy_this_residue_add_chain(CResidue *residue,
						 const std::string &altconf,
						 bool whole_residue_flag,
						 bool attach_to_new_chain_flag);
      CResidue *deep_copy_this_residue(CResidue *residue);
      

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
      std::pair<CMMDBManager *, std::vector<residue_spec_t> > 
      get_fragment_from_atom_spec(const atom_spec_t &atom_spec, CMMDBManager *mol);


      // transform atoms in residue
      void transform_atoms(CResidue *res, const clipper::RTop_orth &rtop);

      // transform all the atom in mol
      void transform_mol(CMMDBManager *mol, const clipper::RTop_orth &rtop);


      // Rotate position round vector, return a position.
      //
      // Note that angle is in radians.
      // 
      clipper::Coord_orth rotate_round_vector(const clipper::Coord_orth &direction,
					      const clipper::Coord_orth &position,
					      const clipper::Coord_orth &origin_shift,
					      double angle);

      // This presumes that a_residue_p and b_residue_p are valid.
      std::vector<std::pair<int, int> > pair_residue_atoms(CResidue *a_residue_p,
							   CResidue *b_residue_p);

      // A useful function that was (is) in molecule_class_info_t
      //
      // If this residue has a CA, return a pointer to that atoms, if
      // not, return the first atom, if no atoms, return NULL.
      // 
      CAtom* intelligent_this_residue_mmdb_atom(CResidue *res_p);

      // A function to find the previous (next) residue.  Find residue
      // by previous (next) serial number.
      // Return NULL if prevous (next) resiude not found.
      // this_residue must be part of a molecule hierarchy.
      //
      // (c.f. get_following_residue()).
      // 
      // (used in pepflip)
      CResidue *previous_residue(CResidue *this_residue);
      CResidue *next_residue(CResidue *this_residue);

      // normal sequence codes, X for non-protein
      std::string three_letter_to_one_letter(const std::string &resname);
      // as above, but allow specials, currently HOH -> "~"
      std::string three_letter_to_one_letter_with_specials(const std::string &resname);


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
		     const std::vector<lsq_range_match_info_t> &matches,
		     int every_nth);
      // used by above
      // On useful return, first.length == second.length and first.length > 0.
      // 
      std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> > 
      get_matching_indices(CMMDBManager *mol1,
			   CMMDBManager *mol2,
			   int SelHnd1,
			   int SelHnd2,
			   const lsq_range_match_info_t &match,
			   int every_nth);

      // Return the status in the first position
      // Return the angle in radians.
      std::pair<short int, double> omega_torsion(CResidue *C_residue, CResidue *N_residue,
						 const std::string &altconf);
      // and along the same lines:
      // Return the angles in radians.
      peptide_torsion_angles_info_t peptide_torsions(CResidue *C_residue, CResidue *N_residue,
						     const std::string &altconf);

      // return "" on no canonical name found
      std::string canonical_base_name(const std::string res_name_in, base_t rna_or_dna);

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
      std::string interesting_things_list_py(const std::vector<atom_spec_t> &v);

      // error_type is e.g. "Z score", "Clash gap"
      std::string
      interesting_things_list_with_fix(const std::vector<atom_spec_and_button_info_t> &v,
				       const std::string error_type);
      std::string
      interesting_things_list_with_fix_py(const std::vector<atom_spec_and_button_info_t> &v,
				          const std::string error_type);


      // return a set of residue specifiers and a string to put on the
      // button.  Atsushi Nakagawa told me about the problem of these
      // flips and how they could be identified by B factor
      // outlierness.
      std::vector<std::pair<atom_spec_t, std::string> >
      gln_asn_b_factor_outliers(CMMDBManager *mol);

      // return the number of cis peptides in mol:
      int count_cis_peptides(CMMDBManager *mol);

      // more info on the real cis peptides derived from atom positions:
      std::vector<cis_peptide_info_t>
      cis_peptides_info_from_coords(CMMDBManager *mol);

      // remove wrong cis_peptides
      void remove_wrong_cis_peptides(CMMDBManager *mol);

      // move waters round protein, fiddle with mol.
      // return the number of moved waters.
      int move_waters_around_protein(CMMDBManager *mol);
      // above uses this more primitive interface
      // 
      std::vector<std::pair<CAtom *, clipper::Coord_orth> >
      symmetry_move_atoms(const std::vector<clipper::Coord_orth> &protein_coords,
			  const std::vector<std::pair<CAtom*, clipper::Coord_orth> > &water_atoms,
			  clipper::Cell cell,
			  clipper::Spacegroup spacegroup);
	 

      // Throw an std::runtime_error exception on
      // not-able-to-extract-cell/symm-info.  (In such a case, we convert a
      // clipper::Message_base to a std::runtime_error).
      // 
      std::pair<clipper::Cell, clipper::Spacegroup> get_cell_symm(CMMDBManager *mol);

      // shove a cell from a clipper cell into the passed mol.
      bool set_mol_cell(CMMDBManager *mol, clipper::Cell cell);


      // caller must check that others has some points in it.
      // 
      double min_dist_to_points(const clipper::Coord_orth &pt,
				const std::vector<clipper::Coord_orth> &others);

      // Return the fractional shift needed to translate the protein
      // as close as possible to the origin (do not apply the shift).
      //
      // Throw an exception (e.g. no cell or symmetry).
      // 
      clipper::Coord_frac shift_to_origin(CMMDBManager *mol);
      clipper::Coord_frac shift_to_origin(const std::vector<clipper::Coord_orth> &protein_coords,
					  clipper::Cell cell,
					  clipper::Spacegroup spacegroup);

      //
      clipper::Mat33<double> residue_orientation(CResidue *residue_p,
						 const clipper::Mat33<double> &orientation_in);

      clipper::Coord_orth average_position(std::vector<clipper::Coord_orth> &pts);

      // Return the median position.  Throw an exception on failure
      // (e.g no atoms).
      // 
      clipper::Coord_orth median_position(CMMDBManager *mol);
      // also throws an exception
      clipper::Coord_orth median_position(const std::vector<clipper::Coord_orth> &pts);
      //
      clipper::Coord_orth
      translate_close_to_origin(const clipper::Coord_orth water_pos_pre,
				const clipper::Cell &cell);

      void translate_close_to_origin(CMMDBManager *mol);

      // Print secondary structure info:
      void print_secondary_structure_info(CModel *model_p);
      
   } // namespace util
   std::ostream&  operator<<(std::ostream&  s, const util::quaternion &q);
   std::ofstream& operator<<(std::ofstream& s, const util::quaternion &q);
   std::ostream& operator<<(std::ostream& s, const atom_spec_t &spec);
   std::ostream& operator<<(std::ostream& s, const residue_spec_t &spec);

} // namespace coot

#endif // HAVE_COOT_COORD_UTILS_HH
