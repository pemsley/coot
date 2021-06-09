/* high-res/sequence-assignment.hh
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


#ifndef HAVE_SEQUENCE_ASSIGNMENT_HH
#define HAVE_SEQUENCE_ASSIGNMENT_HH

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include "clipper/core/xmap.h"
#include <mmdb2/mmdb_manager.h>
#include "geometry/protein-geometry.hh"

namespace coot {

   // We get a vector of these things, they are info where there is
   // structure (poly-ala) but no assigned sequence - they refer to
   // the structure.
   // 
   class residue_range_t {
   public:
      std::string chain_id;
      int chain_id_index;
      int start_resno;
      int end_resno; // inclusive
      residue_range_t(const std::string &chain_id_in,
		      int start_resno_in,
		      int end_resno_in) : chain_id(chain_id_in) {
	 start_resno = start_resno_in;
	 end_resno = end_resno_in;
         chain_id_index = 0;
      }
      int length() const { return end_resno - start_resno + 1; }
   };

   namespace sequence_assignment { 

      enum side_chain_name_index { GLY, ALA, SER, VAL, THR, PRO, ASN, ASP, CYS,
				   GLN, GLU, HIS, ILE, LEU, LYS, MET, PHE, TYR,
				   TRP, ARG};

      std::string side_chain_name_index_to_name(const side_chain_name_index &idx);

      // We get a vector of these things which are info about a range
      // of sequence that has not been associcated with structure
      // (applies to sequence info).  Used in the construction of the
      // slider table.
      // 
      class sequence_range_t {
      public:
	 int chain_id_index;
	 int start_sequence_resno;
	 int end_sequence_resno; // inclusive
	 sequence_range_t(int chain_id_index_in,
			  int start_sequence_resno_in,
			  int end_sequence_resno_in) {
	    chain_id_index = chain_id_index_in;
	    start_sequence_resno = start_sequence_resno_in;
	    end_sequence_resno = end_sequence_resno_in;
	 } 
      };

      // Sequence info: We have one of these for each chain.  It
      // contains sequence info, each element has a chain id and for
      // each residue in the sequence there is a residue_type_as_int
      // (for fast comparisons) and associated with that residue also
      // a score of sureness of the assignemnt (0 to 1). Used in the
      // construction of the slider table.
      // 
      class sequence_info_t {
      public:
	 std::string chain_id;
	 // for each sequence position we set a type
	 // (side_chain_name_index) and a score for it's certainty of
	 // assignment (1.0 is sure).
	 // 
	 std::vector<std::pair<side_chain_name_index, float> > residue_info;
	 sequence_info_t(const std::string &chain_id_in,
			 const std::vector<side_chain_name_index> &a) : chain_id(chain_id_in) {
	    residue_info.resize(a.size());
	    for (unsigned int i=0; i<a.size(); i++) {
	       residue_info[i].first = a[i];
	       residue_info[i].second = 0.0; // nothing in the structure
					     // assigned to this sequence
					     // initially
	    }
	 }
      };

      // This is an older class - a helper class for
      // side_chain_score_t and is used by the slider (I think).
      // 
      class scored_chain_info_t {
	 std::vector<std::vector<float> > residue_side_chain_score;
	 std::string chain_name_;
	 int max_resno() { return residue_side_chain_score.size()-1;}

	 // return -1 on no outstanding slider position.
	 int outstanding_slider_position(const std::vector<float> &s) const;
      public:
	 scored_chain_info_t(const std::string &chin, int n_residues) : chain_name_(chin) {
	    residue_side_chain_score.resize(n_residues+1);
	 }
	 // address the residue number (resno) directly.  Index 0 is
	 // there but not used.
	 // 
	 void add_score(int resno, 
			int residue_index,
			double score);
	 void debug() const; // print score table
	 void debug_with_cout() const;
	 std::string chain_name() const { return chain_name_;}
	 // return -1 on no outstanding slider position:
	 int slider_hit(const std::vector<std::pair<side_chain_name_index, float> > &seq) const;
      };
      
      class side_chain_score_t {

	 // create a copy of the generate scores input data
	 mmdb::Manager *mol;
	 clipper::Xmap<float> *xmap;
	 std::vector <mmdb::Residue *> standard_residues;
	 
	 // pairs of (chain_id and sequence)
	 std::vector<std::pair<std::string, std::string> > input_sequence;
	 // pairs of (chain_id as sequence as int)s
	 // std::vector<std::pair<std::string, std::vector<int> > > sequence_as_indices;
	 // replaces the above:
	 std::vector<sequence_info_t> sequence_infos;
	 short int is_fasta_aa(const std::string &a) const;
	 std::vector<side_chain_name_index> convert_seq_to_indices(const std::string &seq) const;
	 side_chain_name_index convert_slc_to_index(const std::string &code) const;
	 // and a test to see if the above was valid:
	 bool is_valid_residue_type(const side_chain_name_index &i) const {
	    // return (i<999 ? 1 : 0);  20140308 : what did I mean!?
            if (i < 0)
               return false;
	    return true;
	 }
	 // isn't auto_fit_score using a mmdb::Residue a better fit?
	 float auto_fit_score(const std::string &chain_id,
			      int resno,
			      const side_chain_name_index &idx);
	 
	 float auto_fit_score(mmdb::Residue *res, // current res used to fit new one.
			      const side_chain_name_index &idx,
			      const coot::dictionary_residue_restraints_t &rest,
			      const clipper::Xmap<float> &xmap);

	 short int cache_standard_residues();
	 mmdb::Residue *get_standard_residue(const side_chain_name_index &idx) const;

	 void move_std_res_to_this_res_pos(const clipper::RTop_orth &rtop,
					   mmdb::Residue *std_res); // move std_res atoms

	 float best_rotamer_score(const clipper::Xmap<float> &xmap,
				  const dictionary_residue_restraints_t &rest,
				  mmdb::Residue *res) const;

	 // not const because it adds UDD data to mol
	 std::vector<residue_range_t> find_unassigned_regions(float pr_cut);
	 int udd_assigned_handle; // the UDD handle

	 // Find unassigned bits of the sequences that are long enough
	 // to accomodate a_residue_range.  So if the sequence is 100
	 // aa and we have already assigned a range 40-50, there are 2
	 // returned sequences: 1 to 39 and 51 to 100.
	 // 
	 std::vector<sequence_range_t>
	 find_unassigned_sequence(const residue_range_t &a_residue_range) const;
	 // And here is how we mark them up initially when we get them:
	 // (notice that we can also pass a critical value (0.0->1.0) below
	 // which the residues are considered unassigned.  Typical value of
	 // pr_crit is 0.1.
	 void mark_unassigned_residues();

      public:
      // Let's make chains 0-indexed and residues 1-indexed.
	 std::vector<scored_chain_info_t> side_chain_score;
	 side_chain_score_t() { udd_assigned_handle = -1; xmap = 0; mol = 0; };
	 void add_new_chain(const std::string &chain_name, int n_residues) {
	    side_chain_score.push_back(scored_chain_info_t(chain_name, n_residues));
	 }
	 void add_score(int chain_number,
			const std::string &str,
			int residue_number,
			int max_resno_in_chain,
			int residue_idx, // enum side_chain_name_index
			double score);

	 void debug() const;
	 void add_fasta_sequence(const std::string &sequence_chain_id_in,
				 const std::string &seq);
	 void slider() const;
	 void generate_scores(mmdb::Manager *mol, const clipper::Xmap<float> &xmap);

	 void test_residue_range_marking();
      };
   }
}

#endif // HAVE_SEQUENCE_ASSIGNMENT_HH
