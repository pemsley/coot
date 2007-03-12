
#ifndef CMTZ_INTERFACE_HH
#define CMTZ_INTERFACE_HH

#include <string>
#include <vector>

namespace coot {

   class mtz_type_label {
   public:
      char column_type; 
      std::string column_label;
      mtz_type_label(const std::string &column_label_in, char column_type_in) {
	 column_type = column_type_in;
	 column_label = column_label_in;
      } 
   };


   class mtz_column_types_info_t {
   public:
      std::string mtz_filename;
      short int read_success; 
      std::vector<mtz_type_label> f_cols;
      std::vector<mtz_type_label> sigf_cols;
      std::vector<mtz_type_label> d_cols;
      std::vector<mtz_type_label> sigd_cols;
      std::vector<mtz_type_label> phi_cols;
      std::vector<mtz_type_label> weight_cols;
      std::vector<mtz_type_label> r_free_cols;
      int selected_f_col; 
      int selected_phi_col;
      int selected_weight_col;
      int selected_refmac_fobs_col;
      int selected_refmac_sigfobs_col;
      int selected_refmac_r_free_col;
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

   

   



#endif // CMTZ_INTERFACE_HH
