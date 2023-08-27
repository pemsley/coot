#ifndef MAKE_INSTANCED_GRAPHICAL_BONDS_HH
#define MAKE_INSTANCED_GRAPHICAL_BONDS_HH

#include "coords/graphical-bonds-container.hh"
#include "instancing.hh"
#include "bond-colour.hh"

void
make_instanced_graphical_bonds_spherical_atoms(coot::instanced_mesh_t &m, // add to this
                                               const graphical_bonds_container &gbc,
                                               coot::api_bond_colour_t bonds_box_type,
                                               int udd_handle_bonded_type,
                                               float base_atom_radius,
                                               float base_bond_radius,
                                               unsigned int num_subdivisions,
                                               const std::vector<glm::vec4> &colour_table);


void
make_instanced_graphical_bonds_hemispherical_atoms(coot::instanced_mesh_t &m, // add to this
                                                   const graphical_bonds_container &gbc,
                                                   coot::api_bond_colour_t bonds_box_type,
                                                   int udd_handle_bonded_type,
                                                   float atom_radius,
                                                   float bond_radius,
                                                   unsigned int num_subdivisions,
                                                   const std::vector<glm::vec4> &colour_table);

void
make_instanced_graphical_bonds_bonds(coot::instanced_mesh_t &m,
                                     const graphical_bonds_container &gbc,
                                     float bond_radius,
                                     unsigned int n_slices,
                                     unsigned int n_stacks,
                                     const std::vector<glm::vec4> &colour_table);

void make_graphical_bonds_spherical_atoms_with_vdw_radii_instanced(coot::instanced_mesh_t &m, const graphical_bonds_container &gbc,
                                                                   unsigned int num_subdivisions,
                                                                   const std::vector<glm::vec4> &colour_table,
                                                                   const coot::protein_geometry &geom);

void
make_instanced_graphical_bonds_bonds(coot::instanced_mesh_t &m,
                                     const graphical_bonds_container &gbc,
                                     float bond_radius,
                                     unsigned int n_slices,
                                     unsigned int n_stacks,
                                     const std::vector<glm::vec4> &colour_table);

void
make_graphical_bonds_cis_peptides(coot::simple_mesh_t &m,
                                  const graphical_bonds_container &gbc);

#endif // MAKE_INSTANCED_GRAPHICAL_BONDS_HH
