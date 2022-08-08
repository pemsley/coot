#ifndef GTK_MANUAL_HH
#define GTK_MANUAL_HH

#include <gtk/gtk.h>
#include <string>

/* ------------------------------------------------------------------------------ */
/* Additional Representation Handling                                             */
/* ------------------------------------------------------------------------------ */

/* create a frame/combo_box and add it to the
   display_control_molecule_combo_box (or the vbox thereof). */
void display_control_molecule_combo_box(const std::string &name, int n, bool show_add_reps_frame_flag);

GtkWidget *display_control_add_reps_container(GtkWidget *display_control_window_glade,
                                              int imol_no);

GtkWidget *
display_control_add_reps_all_on_check_button(GtkWidget *display_control_window_glade,
                                             int imol_no);


void display_control_add_reps(GtkWidget *add_reps_vbox,
                              int imol_no, int add_rep_no,
                              bool show_it,
                              int bonds_box_type, const std::string &name);
void add_rep_toggle_button_toggled(GtkToggleButton       *button,
                                   gpointer         user_data);

void add_add_reps_frame_and_vbox(GtkWidget *display_control_window_glade,
                                 GtkWidget *hbox_for_single_molecule, int imol_no,
                                 bool show_add_reps_frame_flag);

// bring the interface into C++
void display_control_map_combo_box(const std::string &name, int imol);

#endif // GTK_MANUAL_HH

