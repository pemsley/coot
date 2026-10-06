#include <mmdb2/mmdb_manager.h>
#include "m2t-mesh.hh"

//! This uses residue->SSE == mmdb::SSE_Helix to remove helices.
//! so be sure that you have residue SSEs, e.g. by calling
//! secondary_structure_header_to_residue_sse(mol).
//!
//! Cn is 3
//! accuracy = 12.
//!
//! helix_template_pdb_file_name locates the theoretical Z-aligned poly-ALA helix
//! (data/pdb-templates/theor-helix-z-ori-v2.pdb) that real helix segments get
//! superposed onto. This library can't depend on coot-utils to resolve an
//! installed data-directory path (it's low in the link order), so by default
//! this is a bare filename that only resolves if the current working directory
//! happens to contain it - pass an absolute path (e.g. built from
//! coot::package_data_dir() in the caller) to make this work regardless of cwd.
//! Unused when straight_helices is true (see make_mesh_for_straight_helical_representation()).
//!
//! straight_helices: false (default) superposes the reference helix onto every
//! residue triplet independently (make_mesh_for_helical_representation()) - this
//! tracks real local backbone bending/irregularity faithfully but looks visibly
//! segmented ("wormy") rather than a smooth rod. true fits one single straight axis
//! through the whole helix's CA atoms instead and draws one plain capped cylinder
//! along it (make_mesh_for_straight_helical_representation()), smoothing away that
//! local wobble.
coot::m2t::simple_mesh_t
make_tubes_representation(mmdb::Manager *mol,
                          const std::string &atom_selection_str,
                          const std::string &colour_scheme,
                          float radius_for_coil,
                          int Cn_for_coil, int accuracy_for_coil,
                          unsigned int n_slices_for_coil,
                          int secondaryStructureUsageFlag,
                          const std::string &helix_template_pdb_file_name = "theor-helix-z-ori-v2.pdb",
                          bool straight_helices = false);

//! Typically we might call this function for every chain.
//! For a bendy-helix representation, we would we don't want a
//! coil/tube where the helices are - so, in that case,
//! remove_trace_for_helices = true.
//!
coot::m2t::simple_mesh_t
make_coil_for_tubes_representation(mmdb::Manager *mol,
                                   const std::string &atom_selection_str,
                                   float radius_for_coil,
                                   int Cn_for_coil, int accuracy_for_coil,
                                   unsigned int n_slices_for_coil,
                                   bool remove_trace_for_helices);
coot::m2t::simple_mesh_t
make_mesh_for_helical_representation(mmdb::Manager *mol,
                                     const std::string &atom_selection_str,
                                     float radius_for_helices,
                                     unsigned int n_slices_for_helices);

//! Just the straight-cylinder helix geometry (one PCA-fit axis per helix, see
//! make_mesh_for_straight_helical_representation() in tubes.cc) - no coil/strand.
//! Meant to be merged with a Ribbon representation that has its "hideHelixGeometry"
//! parameter set, so the two don't overlap - see molecule_class_info_t's "TubeHelices"
//! style in src/molecule-class-info-mol-tris.cc.
coot::m2t::simple_mesh_t
make_straight_cylinder_helices_mesh(mmdb::Manager *mol,
                                    const std::string &atom_selection_str,
                                    float radius_for_helices,
                                    unsigned int n_slices_for_helices,
                                    int secondaryStructureUsageFlag);
