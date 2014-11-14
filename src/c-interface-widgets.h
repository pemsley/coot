
/* I think that this should this be a .hh file */

void add_ligand_builder_menu_item_maybe();
void start_ligand_builder_gui(GtkMenuItem     *menuitem,
			      gpointer         user_data);


/* ------------------------------------------------------------------------- */
/*                    Generic Objects                                        */
/* ------------------------------------------------------------------------- */

GtkWidget *wrapped_create_generic_objects_dialog();
/* and this uses callbacks: */
void on_generic_objects_dialog_object_toggle_button_toggled(GtkButton       *button,
							    gpointer         user_data);
/* and... */
void
generic_objects_dialog_table_add_object_internal(const coot::generic_display_object_t &gdo,
						 GtkWidget *dialog,
						 GtkWidget *table,
						 int io);

/* return a new object number (so that we can set it to be displayed). */
int add_generic_display_object(const coot::generic_display_object_t &gdo); 
