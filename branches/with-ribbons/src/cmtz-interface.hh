
#ifndef CMTZ_INTERFACE_HH
#define CMTZ_INTERFACE_HH

#include <string>
#include <vector>

namespace coot {

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
   };

   mtz_column_types_info_t get_f_phi_columns(const std::string &filename);

   std::vector<std::string> get_f_cols(const std::string &mtz_file_name); 
   std::vector<std::string> get_phi_cols(const std::string &mtz_file_name); 
   std::vector<std::string> get_weight_cols(const std::string &mtz_file_name); 
   std::vector<std::string> get_d_cols(const std::string &mtz_file_name); 

   // return 0 on not a valid mtz/data file
   GtkWidget *column_selector_using_cmtz(const std::string &filename);
   void fill_f_optionmenu(GtkWidget *optionmenu_f, short int is_expert_mode_flag);
   void setup_refmac_parameters(GtkWidget *window, 
				const coot::mtz_column_types_info_t &col_labs);
   void setup_refmac_parameters_from_file(GtkWidget *window);
}

GtkWidget *make_menu_item( gchar         *name,
                           GtkSignalFunc  callback,
                           gpointer       data );


void f_button_select(GtkWidget *item, GtkPositionType pos);
void phase_button_select(GtkWidget *item, GtkPositionType pos);
void weight_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_f_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_sigf_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_r_free_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_dialog_f_button_select(GtkWidget *item, GtkPositionType pos);
//void refmac_dialog_sigf_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_dialog_fpm_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_dialog_i_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_dialog_ipm_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_dialog_r_free_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_dialog_phases_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_dialog_fom_button_select(GtkWidget *item, GtkPositionType pos);
void refmac_dialog_hl_button_select(GtkWidget *item, GtkPositionType pos);

   

   



#endif // CMTZ_INTERFACE_HH
