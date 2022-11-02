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