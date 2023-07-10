#include "validation-graphs.hh"

std::string coot::validation_graph_type_to_human_name(coot::validation_graph_type graph_type) {
    switch (graph_type) {
        case validation_graph_type::density_fit: {
            return "Density fit";
            break;
        }
        case validation_graph_type::geometry: {
            return "Geometry analysis";
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
            return "Ramachandran";
            break;
        }
        case validation_graph_type::rota: {
            return "Rotamer analysis";
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
