#ifndef VALIDATION_GRAPH_WIDGET_HH
#define VALIDATION_GRAPH_WIDGET_HH
#include <gtk/gtk.h>



G_BEGIN_DECLS   

#define COOT_VALIDATION_GRAPH_TYPE (coot_validation_graph_get_type ())
G_DECLARE_FINAL_TYPE  (CootValidationGraph, coot_validation_graph, COOT, COOT_VALIDATION_GRAPH, GtkWidget)


CootValidationGraph *coot_validation_graph_new();

G_END_DECLS

#endif // VALIDATION_GRAPH_WIDGET_HH