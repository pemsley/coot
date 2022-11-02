#include "validation-graphs.hh"


void coot::validation_graphs_t::init() {
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

coot::validation_graphs_t::validation_graphs_t() {
    init(); 
}

std::string coot::validation_graph_type_to_human_name(coot::validation_graph_type graph_type) {
    switch (graph_type) {
        case validation_graph_type::density_fit: {
            return "Density fit";
            break;
        }
        case validation_graph_type::geometry: {
            return "Geometry";
            break;
        }
        case validation_graph_type::ncs: {
            return "NCS";
            break;
        }
        case validation_graph_type::omega: {
            return "Omega";
            break;
        }
        case validation_graph_type::rama: {
            return "Rama";
            break;
        }
        case validation_graph_type::rota: {
            return "Rota";
            break;
        }
        case validation_graph_type::temp_factor: {
            return "Temp factor";
            break;
        }
        default: {
            return "unknown";
            break;
        }
    }
}