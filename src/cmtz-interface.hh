
#ifndef CMTZ_INTERFACE_HH
#define CMTZ_INTERFACE_HH

#include <string>
#include <vector>

namespace coot {

   // 2019, the info needed for make_and_draw_map
	//       first example used in extraction of info from refmac updating maps json file
	class mtz_to_map_info_t {
		public:
		std::string mtz_file_name;
		std::string f_col;
		std::string phi_col;
		std::string w_col;
		std::string id; // used to find the molecule name
		bool use_weights;
		bool is_difference_map;
		mtz_to_map_info_t() { use_weights = false; is_difference_map = false; }
	};

   // --------------------- This is messy ----------------------------
   // what was I thinking?

   class mtz_type_label {
   public:
      char column_type; 
      std::string column_label;
      int column_position;
      mtz_type_label(const std::string &column_label_in, char column_type_in, int column_position_in) {
	 column_type = column_type_in;
	 column_label = column_label_in;
	 column_position = column_position_in;
      };
      mtz_type_label() {
	 column_type = 0;
	 column_position = -1;
      }
   };


   class mtz_column_types_info_t {
   public:
      std::string mtz_filename;
      short int read_success; 
      std::vector<mtz_type_label> f_cols;
      std::vector<mtz_type_label> sigf_cols; // contains sigIs as well
      std::vector<mtz_type_label> d_cols;
      std::vector<mtz_type_label> sigd_cols;
      std::vector<mtz_type_label> phi_cols;
      std::vector<mtz_type_label> weight_cols;
      std::vector<mtz_type_label> r_free_cols;
      std::vector<mtz_type_label> hl_cols;
      std::vector<mtz_type_label> fpm_cols;
      std::vector<mtz_type_label> sigfpm_cols;
      std::vector<mtz_type_label> i_cols;
      std::vector<mtz_type_label> ipm_cols;
      std::vector<mtz_type_label> sigipm_cols;
      int selected_f_col; 
      int selected_phi_col;
      int selected_weight_col;
      int selected_refmac_fobs_col;
      int selected_refmac_sigfobs_col;
      int selected_refmac_r_free_col;
      int selected_refmac_phi_col;
      int selected_refmac_fom_col;
      int selected_refmac_hla_col;
      int selected_refmac_hlb_col;
      int selected_refmac_hlc_col;
      int selected_refmac_hld_col;
      int selected_refmac_fp_col;
      int selected_refmac_sigfp_col;
      int selected_refmac_fm_col;
      int selected_refmac_sigfm_col;
      int selected_refmac_iobs_col;
      int selected_refmac_sigiobs_col;
      int selected_refmac_ip_col;
      int selected_refmac_sigip_col;
      int selected_refmac_im_col;
      int selected_refmac_sigim_col;
      int use_weights;

      int get_prefered_label_idx(const std::vector<mtz_type_label> &cols,
				 const std::string &prefered_suffix) const {
	 int prf_idx = 0;
	 for (std::size_t i=0; i<cols.size(); i++) {
	    const std::string &col_lab = cols[i].column_label;
	    if (col_lab.length() >= prefered_suffix.length()) {
	       std::size_t ll = col_lab.length();
	       std::size_t lp = prefered_suffix.length();
	       std::string tail = col_lab.substr(ll-lp, lp);
	       if (tail == prefered_suffix) {
		  prf_idx = i;
		  break;
	       }
	    }
	 }
	 return prf_idx;
      }

      int get_prefered_f_col_idx() const {
	 return get_prefered_label_idx(f_cols, "/FWT");
      }
      int get_prefered_phi_col_idx() const {
	 return get_prefered_label_idx(phi_cols, "/PHWT");
      }
   };

   mtz_column_types_info_t get_mtz_columns(const std::string &filename);

   std::vector<std::string> get_f_cols(const std::string &mtz_file_name); 
   std::vector<std::string> get_sigf_cols(const std::string &mtz_file_name); 
   std::vector<std::string> get_r_free_cols(const std::string &mtz_file_name); 
   std::vector<std::string> get_phi_cols(const std::string &mtz_file_name); 
   std::vector<std::string> get_weight_cols(const std::string &mtz_file_name); 
   std::vector<std::string> get_d_cols(const std::string &mtz_file_name);

}

// 20220723-PE perhaps I should elide this whole file?
#ifdef EMSCRIPTEN_THING
void
my_combo_box_text_add_items(GtkComboBox *combobox,
			    const std::vector<coot::mtz_type_label> &labels,
			    int active_label_index);
#endif

#endif // CMTZ_INTERFACE_HH
