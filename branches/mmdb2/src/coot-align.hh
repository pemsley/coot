#ifndef COOT_ALIGN_HH
#define COOT_ALIGN_HH

#include "coot-utils/coot-coord-utils.hh"

namespace coot {

   class mutate_insertion_range_info_t {
   public:
      int start_resno;
      std::vector<std::string> types;
      int end_resno() const {
	 return start_resno + types.size();
      }
      mutate_insertion_range_info_t(int st, const std::vector<std::string> &types_in) {
	 start_resno = st;
	 types = types_in;
      }
      friend std::ostream& operator <<(std::ostream &s, mutate_insertion_range_info_t &r);
   };
   std::ostream& operator <<(std::ostream &s, mutate_insertion_range_info_t &r);

   class chain_mutation_info_container_t {
   public:
      std::string chain_id;
      std::string alignedS;
      std::string alignedT;
      std::string alignedS_label;
      std::string alignedT_label;
      std::pair<bool, float> alignment_score;
      std::vector<mutate_insertion_range_info_t> insertions;
      std::vector<std::pair<residue_spec_t, std::string> > single_insertions;
      std::vector<residue_spec_t> deletions;
      std::vector<std::pair<residue_spec_t, std::string> > mutations;
      chain_mutation_info_container_t(const std::string &chain_id_in) {
	 chain_id = chain_id_in;
      }
      chain_mutation_info_container_t() {
	 alignment_score = std::pair<bool, float> (0,0);
      }
      void add_deletion(const residue_spec_t &res_spec) {
	 residue_spec_t r = res_spec;
	 r.chain = chain_id; // in case it was not set by the function
		             //  that called add_insertion.
	 deletions.push_back(r);
      }
      void add_mutation(const residue_spec_t &res_spec, const std::string &target_type) {
	 residue_spec_t r = res_spec;
	 r.chain = chain_id; // in case it was not set by the function
		             //  that called add_insertion.
	 mutations.push_back(std::pair<residue_spec_t, std::string> (r, target_type));
      }
      void add_insertion(const residue_spec_t &res_spec, const std::string &target_type) {
	 residue_spec_t r = res_spec;
	 r.chain = chain_id; // in case it was not set by the function
		             //  that called add_insertion.
	 single_insertions.push_back(std::pair<residue_spec_t, std::string> (r, target_type));
      }
      void rationalize_insertions(); // move stuff from single
				     // insertions to insertions.

      // throw an execption if there is no residue type to return for
      // the given spec.
      std::string get_residue_type(const residue_spec_t &spec) const;
      
      void print() const;

      double dissimilarity_score() const;
   };

}

#endif // COOT_ALIGN_HH
