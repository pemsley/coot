#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gtk/gtk.h>

#include "callbacks.h"
#include "interface.h"
#include "support.h"

#include "restraints-editor-c.h"


void
on_restraint_editor_add_restraint_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  restraints_editor_add_restraint_by_widget(w);
}


void
on_restraints_editor_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  gtk_widget_destroy(w);

}


void
on_restraints_editor_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{

  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  apply_restraint_by_widget(w);
}



void
on_restraint_editor_delete_restraint_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data)
{
  GtkWidget *w = lookup_widget(GTK_WIDGET(button), "restraints_editor_dialog");
  restraints_editor_delete_restraint_by_widget(w);
}

