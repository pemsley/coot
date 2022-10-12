#include "validation-graph-widget.hh"

G_BEGIN_DECLS
struct _CootValidationGraph {
    GtkWidget parent;
};

G_DEFINE_TYPE(CootValidationGraph, coot_validation_graph, GTK_TYPE_WIDGET)

// not sure what this is for or whether it is going to be needed at all
// struct _CootValidationGraphClass {
//     GObjectClass parent_class;
// };

static void coot_validation_graph_init(CootValidationGraph* ) {
    // not sure what to put here yet
}

static void coot_validation_graph_class_init(CootValidationGraphClass* klass) {
    // not sure what to put here yet
}

CootValidationGraph* 
coot_validation_graph_new()
{
  return COOT_COOT_VALIDATION_GRAPH(g_object_new (COOT_VALIDATION_GRAPH_TYPE, NULL));
}



G_END_DECLS
