
#ifndef VALIDATION_GRAPHS_HH
#define VALIDATION_GRAPHS_HH

#include <gtk/gtk.h>
#include <type_traits>
namespace coot {


    enum class validation_graph_type : unsigned char {
        rota,
        temp_factor,
        density_fit,
        rama,
        omega,
        geometry,
        ncs
    };
    /// For using as unsigned char like so:
    /// `static_cast<validation_graph_type_repr_t>(enum_instance)`
    typedef std::underlying_type<validation_graph_type>::type validation_graph_type_repr_t;


    //todo: remove this class later
    class validation_graphs_t {

    public:
        GtkWidget *dynarama_is_displayed;
        GtkWidget *sequence_view_is_displayed;
        GtkWidget *geometry_graph;
        GtkWidget *b_factor_variance_graph;
        ////B B GRAPH
        GtkWidget *b_factor_graph;
        ////E B GRAPH
        GtkWidget *residue_density_fit_graph;
        GtkWidget *omega_distortion_graph;
        GtkWidget *rotamer_graph;
        GtkWidget *ncs_diffs_graph;
        validation_graphs_t() { init(); }

        void init() {
        dynarama_is_displayed = NULL;
        sequence_view_is_displayed = NULL;
        geometry_graph = NULL;
        b_factor_variance_graph = NULL;
        b_factor_graph = NULL;
        residue_density_fit_graph = NULL;
        omega_distortion_graph = NULL;
        rotamer_graph = NULL;
        ncs_diffs_graph = NULL;
        }
    };
} // namespace coot

#endif // VALIDATION_GRAPHS_HH

