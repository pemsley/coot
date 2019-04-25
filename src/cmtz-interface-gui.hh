
namespace coot { 
   // return 0 on not a valid mtz/data file
   GtkWidget *column_selector_using_cmtz(const std::string &filename);
   void fill_f_optionmenu(GtkWidget *optionmenu_f, short int is_expert_mode_flag);
   void setup_refmac_parameters(GtkWidget *window, 
				const coot::mtz_column_types_info_t &col_labs);
   void setup_refmac_parameters_from_file(GtkWidget *window);

}

GtkWidget *make_menu_item( gchar         *name,
                           GCallback  callback,
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

   

