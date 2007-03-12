
#include "coot-coord-utils.hh"

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
   };

   class chain_mutation_info_container_t {
   public:
      std::string chain_id;
      std::string alignedS;
      std::string alignedT;
      std::vector<mutate_insertion_range_info_t> insertions;
      std::vector<std::pair<residue_spec_t, std::string> > single_insertions;
      std::vector<residue_spec_t> deletions;
      std::vector<std::pair<residue_spec_t, std::string> > mutations;
      chain_mutation_info_container_t(const std::string &chain_id_in) {
	 chain_id = chain_id_in;
      }
      void add_deletion(const residue_spec_t &res_spec) {
	 deletions.push_back(res_spec);
      }
      void add_mutation(const residue_spec_t &res_spec, const std::string &target_type) {
	 mutations.push_back(std::pair<residue_spec_t, std::string> (res_spec, target_type));
      }
      void add_insertion(const residue_spec_t &res_spec, const std::string &target_type) {
	 single_insertions.push_back(std::pair<residue_spec_t, std::string> (res_spec, target_type));
      }
      void rationalize_insertions(); // move stuff from single
				     // insertions to insertions.
   };

}
