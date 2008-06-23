
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
      friend std::ostream& operator <<(std::ostream &s, mutate_insertion_range_info_t &r);
   };
   std::ostream& operator <<(std::ostream &s, mutate_insertion_range_info_t &r);

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
      
      void print() const {
	 std::cout << "The alignment resulted in the following" << std::endl;
	 std::cout << "   Insertions (coalesced):" << std::endl;
	 for (unsigned int i_ins=0; i_ins<insertions.size(); i_ins++) {
	    std::cout << "       from " << insertions[i_ins].start_resno << " to "
		      << insertions[i_ins].end_resno() << " ";
	    for (unsigned int it=0; it<insertions[i_ins].types.size(); it++) {
	       std::cout << insertions[i_ins].types[it] << " ";
	    }
	    std::cout << std::endl;
	 }
	 std::cout << "   Insertions (singles):" << std::endl;
	 for (unsigned int i_ins=0; i_ins<single_insertions.size(); i_ins++) {
	    std::cout << "    " << single_insertions[i_ins].first << " -> "
		      << single_insertions[i_ins].second << std::endl;
	 }
	 
	 std::cout << "   Deletions:" << std::endl;
	 for (unsigned int i_del=0; i_del<deletions.size(); i_del++) {
	    std::cout << "      " << deletions[i_del] << std::endl;
	 }
	 std::cout << "   Mutations:" << std::endl;
	 for (unsigned int i_mut=0; i_mut<mutations.size(); i_mut++) {
	    std::cout << "      " << mutations[i_mut].first << " -> "
		      << mutations[i_mut].second << std::endl;
	 }
      }
   };

}
