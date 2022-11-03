
#ifndef VALIDATION_GRAPHS_HH
#define VALIDATION_GRAPHS_HH

#include <gtk/gtk.h>
#include <type_traits>
#include <string>
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
    std::string validation_graph_type_to_human_name(validation_graph_type graph_type);

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
        validation_graphs_t();

        void init();
    };
} // namespace coot

#endif // VALIDATION_GRAPHS_HH

