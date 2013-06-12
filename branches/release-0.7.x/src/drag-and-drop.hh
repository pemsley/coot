
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS extern
#define END_C_DECLS     
#endif


BEGIN_C_DECLS

enum { TARGET_STRING };

gboolean
on_gl_canvas_drag_drop(GtkWidget *widget,
		       GdkDragContext *context,
		       gint x, gint y,
		       guint time,
		       gpointer user_data);


void
on_drag_data_received (GtkWidget *widget, 
		       GdkDragContext *context, 
		       gint x, gint y,
		       GtkSelectionData *selection_data, 
		       guint target_type, 
		       guint time,
		       gpointer data);

int handle_drag_and_drop_single_item(const std::string &file_name);

END_C_DECLS

