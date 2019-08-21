
#ifndef UPDATING_MAP_PARAMS_T_HH
#define UPDATING_MAP_PARAMS_T_HH

class updating_map_params_t {
public:
   updating_map_params_t() {
      imol = -1;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
   }
   updating_map_params_t(int imol_in,
			 const std::string &mtz_file_name_in,
			 const std::string &f_col_in,
			 const std::string &phi_col_in,
			 const std::string &w_col_in,
			 bool use_weights_in, bool is_difference_map_in) {
      imol = imol_in;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
      mtz_file_name = mtz_file_name_in;
      f_col = f_col_in;
      phi_col = phi_col_in;
      weight_col = w_col_in;
      use_weights = use_weights_in;
      is_difference_map = is_difference_map_in;
      ctime.tv_sec = 0;
      ctime.tv_nsec = 0;
   }
   int imol;
   std::string mtz_file_name;
   std::string f_col;
   std::string phi_col;
   std::string weight_col;
   bool use_weights;
   bool is_difference_map;
   timespec ctime;
};

#endif // UPDATING_MAP_PARAMS_T_HH
