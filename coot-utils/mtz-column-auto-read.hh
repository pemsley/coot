
#include <string>

namespace coot { 
   class mtz_column_trials_info_t {
   public:
      std::string f_col;
      std::string phi_col;
      std::string w_col;
      bool use_weights;
      bool is_diff_map;
      mtz_column_trials_info_t(const std::string &f_col_in,
			       const std::string &phi_col_in,
			       const std::string &w_col_in,
			       bool w_in, bool d_in) : f_col(f_col_in), phi_col(phi_col_in), w_col(w_col_in) {
	 use_weights = w_in;
	 is_diff_map = d_in;
      }
      mtz_column_trials_info_t(const std::string &f_col_in,
			       const std::string &phi_col_in,
			       bool d_in) : f_col(f_col_in), phi_col(phi_col_in), w_col(std::string("")) {
	 use_weights = false;
	 is_diff_map = d_in;
      }
   };
}
