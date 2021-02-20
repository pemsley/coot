
#ifndef RAMA_MOUSEY_HH
#define RAMA_MOUSEY_HH


#include <gtk/gtk.h>

#include <goocanvas.h>
#include "rama_plot.hh"
 
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS extern
#define END_C_DECLS
#endif

BEGIN_C_DECLS

// Callback Functions from declared in setup_canvas
// (i.e. can be c++ functions)
// 
gint rama_button_press (GtkWidget *widget, GdkEventButton *event);
gint rama_motion_notify(GtkWidget *widget, GdkEventMotion *event);
gint rama_key_press_event(GtkWidget *widget, GdkEventKey *event); // not used.
gint rama_key_release_event(GtkWidget *widget, GdkEventKey *event); 

//menu
void
on_rama_open_menuitem_activate(GtkMenuItem *item, gpointer user_data);
void
on_rama_print_menuitem_activate(GtkMenuItem *item, gpointer user_data);
void
on_rama_save_as_png_menuitem_activate(GtkMenuItem *item, gpointer user_data);
void
on_rama_close_menuitem_activate(GtkMenuItem *item, gpointer user_data);
void
on_rama_radiomenuitem_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data);
void
on_kleywegt_radiomenuitem_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data);
void
on_zoom_in_activate(GtkMenuItem *item, gpointer user_data);
void
on_zoom_out_activate(GtkMenuItem *item, gpointer user_data);
void
on_zoom_resize_menuitem_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data);
void
on_rama_about_menuitem_activate(GtkMenuItem *item, gpointer user_data);

//buttons
void
on_dynarama_selection_checkbutton_toggled(GtkToggleButton *button, gpointer user_data);
void
on_dynarama_selection_entry_activate(GtkEntry *entry, gpointer  user_data);
void
on_dynarama_selection_apply_button_clicked(GtkButton *button, gpointer user_data);
void
on_dynarama2_outliers_only_togglebutton_toggled(GtkToggleButton *button, gpointer user_data);
void
on_dynarama2_zoom_resize_togglebutton_toggled(GtkToggleButton *button, gpointer user_data);

//psi axis
void
on_psi_axis_classic_radioitem_toggled(GtkToggleButton *button, gpointer user_data);
void
on_psi_axis_paule_radioitem_toggled(GtkToggleButton *button, gpointer user_data);

//about dialog
void
on_rama_aboutdialog1_close(GtkDialog *dialog, gpointer user_data);
void
on_rama_aboutdialog1_response(GtkDialog *dialog, gint response_id, gpointer user_data);

// file chooser
void
on_rama_export_as_pdf_filechooserdialog_close(GtkDialog *dialog, gpointer user_data);
void
on_rama_export_as_pdf_filechooserdialog_response(GtkDialog *dialog, gint response_id, gpointer user_data);
void
on_rama_export_as_png_filechooserdialog_close(GtkDialog *dialog, gpointer user_data);
void
on_rama_export_as_png_filechooserdialog_response(GtkDialog *dialog, gint response_id, gpointer user_data);
void
on_rama_open_filechooserdialog_close(GtkDialog *dialog, gpointer user_data);
void
on_rama_export_as_png_filechooserdialog_response(GtkDialog *dialog, gint response_id, gpointer user_data);

// new functions
gboolean rama_item_button_press (GooCanvasItem *item, GooCanvasItem *target,
                                 GdkEventButton *event, gpointer data);
gboolean rama_item_button_release (GooCanvasItem *item, GooCanvasItem *target,
                                 GdkEventButton *event, gpointer data);
gboolean rama_item_enter_event (GooCanvasItem *item, GooCanvasItem *target,
                                 GdkEventCrossing *event, gpointer data);
gboolean rama_item_motion_event (GooCanvasItem *item, GooCanvasItem *target,
                                 GdkEventMotion *event, gpointer data);

// Extern C Functions called from callbacks.c
//
// (i.e.) must be extern c functions
void rama_show_preferences();
void rama_zoom_in(GtkWidget *widget);
void rama_zoom_out(GtkWidget *widget); 
gboolean rama_configure_event(GtkWidget *widget, GdkEventConfigure *event, gpointer data);

END_C_DECLS


#endif // RAMA_MOUSEY_HH

