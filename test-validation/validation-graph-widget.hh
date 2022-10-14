#ifndef VALIDATION_GRAPH_WIDGET_HH
#define VALIDATION_GRAPH_WIDGET_HH
#include <gtk/gtk.h>
#include "validation-information.hh"
#include <memory>

G_BEGIN_DECLS   

#define COOT_VALIDATION_GRAPH_TYPE (coot_validation_graph_get_type ())
G_DECLARE_FINAL_TYPE  (CootValidationGraph, coot_validation_graph, COOT, COOT_VALIDATION_GRAPH, GtkWidget)


CootValidationGraph *coot_validation_graph_new();

G_END_DECLS

void coot_validation_graph_set_validation_information(CootValidationGraph* self, std::unique_ptr<coot::validation_information_t> vi);

#endif // VALIDATION_GRAPH_WIDGET_HH