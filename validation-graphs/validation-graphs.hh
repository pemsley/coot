
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
        density_correlation,
        rama,
        omega,
        geometry,
        ncs
    };
    /// For using as unsigned char like so:
    /// `static_cast<validation_graph_type_repr_t>(enum_instance)`
    typedef std::underlying_type<validation_graph_type>::type validation_graph_type_repr_t;
    std::string validation_graph_type_to_human_name(validation_graph_type graph_type);
} // namespace coot

#endif // VALIDATION_GRAPHS_HH

