#ifndef VALIDATION_GRAPH_WIDGET_HH
#define VALIDATION_GRAPH_WIDGET_HH
#include <gtk/gtk.h>
#include "validation-information.hh"

#if __cplusplus > 201402L
    // good times

    #include <memory>
#else
    // cope

    #include <memory>
    // Copied from: https://gist.github.com/chinmaygarde/970fd5bbd124754b7d36
    // Thank you kind man
    namespace std {
        template <typename T, typename... Args>
        unique_ptr<T> make_unique(Args&&... args) {
            return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
        }
    }
#endif


G_BEGIN_DECLS   

#define COOT_VALIDATION_GRAPH_TYPE (coot_validation_graph_get_type ())
G_DECLARE_FINAL_TYPE  (CootValidationGraph, coot_validation_graph, COOT, COOT_VALIDATION_GRAPH, GtkWidget)


CootValidationGraph *coot_validation_graph_new();
void coot_validation_graph_set_horizontal_zoom_scale(CootValidationGraph* self, float scale);

/// Enable single-chain mode for the given chain ID.
/// Single-chain mode is turned off if the `chain_id` is a null pointer
void coot_validation_graph_set_single_chain_mode(CootValidationGraph* self, const char* chain_id);

G_END_DECLS

void coot_validation_graph_set_validation_information(CootValidationGraph* self, std::shared_ptr<coot::validation_information_t> vi);


#endif // VALIDATION_GRAPH_WIDGET_HH