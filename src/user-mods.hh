
#include <string>
#include <vector>

#ifdef USE_GUILE
#include <libguile.h>
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
// no fix up needed 
#else    
#endif // SCM version
#endif // USE_GUILE

#ifdef USE_PYTHON
#include "c-interface-python.hh"
#endif 


namespace coot { 

   class flips_container {

      class flip {
      public:
	 coot::atom_spec_t atom_spec;
	 std::string residue_type;
	 std::string info_string;
	 std::string set_string;
	 float score;
	 flip (const coot::atom_spec_t &spec,
	       const std::string &set_string_in,
	       const std::string &info_string_in,
	       const std::string &res_type_in,
	       float score_in) {
	    atom_spec = spec;
	    residue_type = res_type_in;
	    info_string = info_string_in;
	    set_string = set_string_in;
	    score = score_in;
	 }
	 void print() const;
      };

      class no_adjust {
	 std::string info_string_;
      public:
	 std::vector<coot::atom_spec_t> atom_specs;
	 no_adjust(const std::vector<coot::atom_spec_t> &specs_in,
		   const std::string &info_string_in) {
	    atom_specs = specs_in;
	    info_string_ = info_string_in;
	 }
	 std::string info_string() const { return info_string_;}
      };

      std::vector<flip> flips;
      std::vector<no_adjust> no_adjustments;
      void parse_set_or_single(const std::string &line);
      void parse_no_adj(const std::string &line);
   
   public:
      flips_container(const std::string &file_name);
      void handle_user_mods(const std::string &filename);
      std::vector<std::string>
      get_user_mods(const std::string &filename) const;
      void store(const std::vector<std::string> &lines);
      unsigned int n_flips() const { return flips.size();}
      flip get_flip(unsigned int i) { return flips[i]; }
#ifdef USE_GUILE
      SCM user_mods() const;
#endif // USE_GUILE      
#ifdef USE_PYTHON
      PyObject *user_mods_py() const;
#endif // USE_PYTHON
   };
}
