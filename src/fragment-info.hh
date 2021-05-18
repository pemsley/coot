
namespace coot {

   class fragment_info_t {
   public:
      class fragment_range_t {
      public:
	 residue_spec_t start_res;
	 residue_spec_t end_res;
	 fragment_range_t(const residue_spec_t &r1, const residue_spec_t &r2) {
	    start_res = r1;
	    end_res = r2;
	 }
      };
      std::string chain_id;
      std::vector<fragment_range_t> ranges;
      fragment_info_t() {}
      fragment_info_t(const std::string chain_id_in) { chain_id = chain_id_in; } 
      void add_range(const fragment_range_t &r) {
	 ranges.push_back(r);
      } 
   };


}
