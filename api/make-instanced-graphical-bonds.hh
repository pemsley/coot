#ifndef MAKE_INSTANCED_GRAPHICAL_BONDS_HH
#define MAKE_INSTANCED_GRAPHICAL_BONDS_HH

#include "coords/graphical-bonds-container.hh"
#include "instancing.hh"
#include "bond-colour.hh"
#include "coords/mmdb-crystal.hh" // for symm_trans_t

void
make_instanced_graphical_bonds_spherical_atoms(coot::instanced_mesh_t &m, // add to this
                                               const graphical_bonds_container &gbc,
                                               coot::api_bond_colour_t bonds_box_type,
                                               float base_atom_radius,
                                               float base_bond_radius,
                                               bool show_atoms_as_aniso_flag,
                                               float aniso_probability, // 0.0 to 1.0
                                               bool show_aniso_atoms_as_ortep,
                                               bool show_atoms_as_empty_flag,
                                               unsigned int num_subdivisions,
                                               const std::vector<glm::vec4> &colour_table);


void
make_instanced_graphical_bonds_hemispherical_atoms(coot::instanced_mesh_t &m, // add to this
                                                   const graphical_bonds_container &gbc,
                                                   coot::api_bond_colour_t bonds_box_type,
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
                                                                   const coot::protein_geometry &geom, int imol=coot::protein_geometry::IMOL_ENC_ANY);

void
make_instanced_graphical_bonds_bonds(coot::instanced_mesh_t &m,
                                     const graphical_bonds_container &gbc,
                                     float bond_radius,
                                     unsigned int n_slices,
                                     unsigned int n_stacks,
                                     const std::vector<glm::vec4> &colour_table);

void
make_instanced_graphical_symmetry_bonds_bonds(coot::instanced_mesh_t &m,
                                              std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > symmetry_bonds_box,
                                              float bond_radius,
                                              unsigned int n_slices,
                                              unsigned int n_stacks,
                                              const std::vector<glm::vec4> &colour_table);


void
make_graphical_bonds_cis_peptides(coot::simple_mesh_t &m,
                                  const graphical_bonds_container &gbc);

#endif // MAKE_INSTANCED_GRAPHICAL_BONDS_HH
