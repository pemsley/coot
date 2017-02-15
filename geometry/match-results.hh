

namespace coot {

   // a container for the results of the comparison vs CCP4SRS graph matching.
   //
   class match_results_t {
   public:
      bool success;
      std::string name;
      std::string comp_id;
      mmdb::Residue *res;
      std::vector<std::pair<int, int> > graph_match_atom_indices;
      // clipper::RTop_orth
      match_results_t(const std::string &comp_id_in, const std::string &name_in, mmdb::Residue *res_in) {
	 name = name_in;
	 comp_id = comp_id_in;
	 res = res_in;
	 if (res_in)
	    success = true;
	 else
	    success = false;
      }
   };

}
