#ifndef GTK_MANUAL_HH
#define GTK_MANUAL_HH

#include <gtk/gtk.h>
#include <string>

/* ------------------------------------------------------------------------------ */
/* Additional Representation Handling                                             */
/* ------------------------------------------------------------------------------ */

GtkWidget *display_control_add_reps_container(GtkWidget *display_control_window_glade,
					      int imol_no);
void display_control_add_reps(GtkWidget *add_reps_vbox,
			      int imol_no, int add_rep_no, 
			      bool show_it,
			      int bonds_box_type, const std::string &name);
void add_rep_toggle_button_toggled(GtkToggleButton       *button,
				   gpointer         user_data);

void add_add_reps_frame_and_vbox(GtkWidget *display_control_window_glade, 
				 GtkWidget *hbox_for_single_molecule, int imol_no);

#endif // GTK_MANUAL_HH


