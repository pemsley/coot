

class lsq_dialog_values_t {

public:
   int reference_molecule_number;
   int moving_molecule_number;
   int ref_res_range_start;
   int ref_res_range_end;
   int mov_res_range_start;
   int mov_res_range_end;
   std::string chain_id_ref;
   std::string chain_id_mov;
   lsq_dialog_values_t() {
      reference_molecule_number = -1;
      moving_molecule_number = -1;
      ref_res_range_start = 1;
      ref_res_range_end = 999;
      mov_res_range_start = 1;
      mov_res_range_end = 999;
   }
   void update(GtkWidget *dialog) {
      if (dialog) {
	 GtkWidget *mov_option_menu = lookup_widget(dialog, "least_squares_moving_molecule_optionmenu");
	 GtkWidget *ref_option_menu = lookup_widget(dialog, "least_squares_reference_molecule_optionmenu");
	 GtkWidget *ref_res_range_1 = lookup_widget(dialog, "least_squares_reference_range_1_entry");
	 GtkWidget *ref_res_range_2 = lookup_widget(dialog, "least_squares_reference_range_2_entry");
	 GtkWidget *mov_res_range_1 = lookup_widget(dialog, "least_squares_moving_range_1_entry");
	 GtkWidget *mov_res_range_2 = lookup_widget(dialog, "least_squares_moving_range_2_entry");

	 const char *txt_r_1 = gtk_entry_get_text(GTK_ENTRY(ref_res_range_1));
	 const char *txt_r_2 = gtk_entry_get_text(GTK_ENTRY(ref_res_range_2));
	 const char *txt_m_1 = gtk_entry_get_text(GTK_ENTRY(mov_res_range_1));
	 const char *txt_m_2 = gtk_entry_get_text(GTK_ENTRY(mov_res_range_2));
	 int ir1 = clipper::String(txt_r_1).i();
	 int ir2 = clipper::String(txt_r_2).i();
	 int im1 = clipper::String(txt_m_1).i();
	 int im2 = clipper::String(txt_m_2).i();

	 ref_res_range_start = ir1;
	 ref_res_range_end   = ir2;
	 mov_res_range_start = im1;
	 mov_res_range_end   = im2;
      }
   }
   void update(GtkWidget *dialog, const std::string &chain_id_1, const std::string &chain_id_2) {
      update(dialog);
      chain_id_ref = chain_id_1;
      chain_id_mov = chain_id_2;
   }
};
