#include <mmdb2/mmdb_manager.h>
#include "coot-utils/simple-mesh.hh"

//! This uses residue->SSE == mmdb::SSE_Helix to remove helices.
//! so be sure that you have residue SSEs, e.g. by calling
//! secondary_structure_header_to_residue_sse(mol).
//!
//! Cn is 3
//! accuracy = 12.
//!
coot::simple_mesh_t
make_tubes_representation(mmdb::Manager *mol,
                          const std::string &atom_selection_str,
                          const std::string &colour_scheme,
                          float radius_for_coil,
                          int Cn_for_coil, int accuracy_for_coil,
                          unsigned int n_slices_for_coil,
                          int secondaryStructureUsageFlag);

//! Typically we might call this function for every chain.
//! For a bendy-helix representation, we would we don't want a
//! coil/tube where the helices are - so, in that case,
//! remove_trace_for_helices = true.
//!
coot::simple_mesh_t
make_coil_for_tubes_representation(mmdb::Manager *mol,
                                   const std::string &atom_selection_str,
                                   float radius_for_coil,
                                   int Cn_for_coil, int accuracy_for_coil,
                                   unsigned int n_slices_for_coil,
                                   bool remove_trace_for_helices);
coot::simple_mesh_t
make_mesh_for_helical_representation(mmdb::Manager *mol,
                                     const std::string &atom_selection_str,
                                     float radius_for_helices,
                                     unsigned int n_slices_for_helices);
