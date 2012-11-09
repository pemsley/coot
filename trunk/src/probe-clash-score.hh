
#ifndef PROBE_CLASH_SCORE_HH
#define PROBE_CLASH_SCORE_HH

namespace coot {

   // " B  72 CYS  HG A: B  53 HIS  H   "
   // 
   class probe_atom_spec_t : public atom_spec_t {
   public: 
      probe_atom_spec_t(const std::string &s) : atom_spec_t() {
	 if (s.length() > 14) { 
	    std::string chain_local = s.substr(0,2);
	    std::string res_no_str = s.substr(2, 4);
	    std::string atom_name_local = s.substr(11, 4);
	    try {
	       int resno_local = coot::util::string_to_int(res_no_str);
	       if (chain_local[0] == ' ')
		  if (chain_local.length() > 1)
		     chain = std::string(chain_local.substr(1));
	       else 
		  chain = chain_local;
	       resno = resno_local;
	       atom_name = atom_name_local;
	    }
	    catch (const std::exception &e) {
	       std::cout << "WARNING:: " << e.what() << std::endl;
	    }
	 }
      }
   };

   class probe_clash_score_t {
   public:
      bool filled;
      int n_bad_overlaps;
      int n_hydrogen_bonds;
      int n_small_overlaps;
      int n_close_contacts;
      int n_wide_contacts;
      probe_clash_score_t() {
	 filled = false;
      }
      probe_clash_score_t(const std::string &dots_file_name);
   }; 

   
   // couldn't get this to work - so I did it long hand.
   class spec_eraser {
   public:
      std::map<std::pair<probe_atom_spec_t, probe_atom_spec_t>, bool> ref_specs;
      spec_eraser(const std::map<std::pair<probe_atom_spec_t, probe_atom_spec_t>, bool> &ref_specs_in) {
	 ref_specs = ref_specs_in;
      }
      bool operator() (const std::pair<probe_atom_spec_t, probe_atom_spec_t> &s) const {
	 return true;
      } 
   };
}


#endif // PROBE_CLASH_SCORE_HH
