/*
 * src/molecular-mesh-generator.hh
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */


#ifndef GM_INFO_T_HH
#define GM_INFO_T_HH

// is this a molecular mesh generator?

#include <vector>
#include <chrono>
#include <string>
#include <memory>
#include <map>
#include <mmdb2/mmdb_manager.h>

#include "utils/colour-holder.hh"

#ifdef USE_MOLECULES_TO_TRIANGLES
//  #include <MoleculesToTriangles/CXXClasses/RendererGL.h>
#include <MoleculesToTriangles/CXXClasses/Light.h>
#include <MoleculesToTriangles/CXXClasses/Camera.h>
#include <MoleculesToTriangles/CXXClasses/SceneSetup.h>
#include <MoleculesToTriangles/CXXClasses/ColorScheme.h>
#include <MoleculesToTriangles/CXXClasses/MyMolecule.h>
#include <MoleculesToTriangles/CXXClasses/RepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/VertexColorNormalPrimitive.h>
#endif

#include "generic-vertex.hh"
#include "coot-utils/g_triangle.hh"
#ifdef THIS_IS_HMT
#else
#include "coot-utils/cis-peptide-info.hh"
#endif
#include "molecular-triangles-mesh.hh"

#include "geometry/protein-geometry.hh"

class molecular_mesh_generator_t {
   std::vector<glm::vec3> generate_spline_points(std::vector<glm::vec3> &control_points, unsigned int n=20);
   void smooth_vertices(std::vector<s_generic_vertex> *v_p, unsigned int idx_begin, unsigned int idx_end);
   float get_torsion_angle(const std::vector<glm::vec3> &cis_pep_quad) const;
   std::vector<std::vector<std::pair<mmdb::Residue *, glm::vec3> > > make_CA_fragments(mmdb::Manager *mol) const;
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
   dodec_at_position(float scale, const glm::vec3 &position) const;
   void add_to_mesh(std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > *vp, // update vp
                    const std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > &vp_new) const;
   void add_to_mesh(std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > *vp, // update vp
                    const std::vector<s_generic_vertex> &gv, const std::vector<g_triangle> &tris) const;

#ifdef USE_MOLECULES_TO_TRIANGLES
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
   molecular_representation_instance_to_mesh(std::shared_ptr<MolecularRepresentationInstance> mol_rep_inst,
                                                         const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                                         const std::vector<std::pair<std::string, int> >   &M2T_int_params);
#endif
   int get_max_resno_for_polymer(mmdb::Chain *chain_p) const;
   std::vector<std::pair<std::string, std::string> > selection_colours;

public:
   molecular_mesh_generator_t() {
      t_previous = std::chrono::high_resolution_clock::now();
      t_start = t_previous;
      fill_atom_positions();
      update_mats_and_colours();
   }

   // This exists so that we can send over the information about the use-picked range-pair (for highlighting)
   class range_t {
   public:
      range_t() : is_valid(false) {};
      range_t(const std::string &c, int r1, int r2) : is_valid(true), chain_id(c), resno_start(r1), resno_end(r2) {};
      bool is_valid;
      std::string chain_id;
      int resno_start;
      int resno_end;
   };
   std::chrono::time_point<std::chrono::high_resolution_clock> t_previous;
   std::chrono::time_point<std::chrono::high_resolution_clock> t_start;
   std::vector<glm::vec3> atom_positions;
   std::vector<glm::mat4> bond_instance_matrices;
   std::vector<glm::vec4> bond_instance_colours;
   std::vector<float> original_bond_lengths;
   std::vector<glm::mat4> atom_instance_matrices;
   std::vector<glm::vec4> atom_instance_colours;
   // the atoms that contribute to the ith GM restraint
   std::map<unsigned int, std::tuple<unsigned int, unsigned int, float, bool> > bond_atom_pair_indices;
   std::vector<std::tuple<unsigned int, unsigned, float> > atom_atom_pair_indices;
   void fill_atom_positions();
   void update_mats_and_colours();
   void move_the_atoms_and_update_the_instancing_matrices();
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > get_cis_peptides(const std::string &pdb_file_name);
   std::map<int, std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > > get_cis_peptides_mesh(mmdb::Manager *mol);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > get_test_cis_peptides();

   // key is the model.
   // type: 0 uset
   // type: 1 cis
   // type: 2 pre-PRO cis
   // type: 3 twisted-trans

   std::map<int, std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > > make_cis_peptide_quads_mesh(mmdb::Manager *mol);

#ifdef THIS_IS_HMT
#else
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > make_cis_peptide_geom(const std::vector<glm::vec3> &cis_pep_quad,
                                                                                            coot::util::cis_peptide_quad_info_t::type_t type);
#endif
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > get_test_twisted_trans_peptides();
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > make_twisted_trans_peptide_geom(const std::vector<glm::vec3> &cis_pep_quad);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > get_test_molecular_triangles_mesh(mmdb::Manager *mol, int secondary_structure_usage_flag);

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
   get_molecular_triangles_mesh(mmdb::Manager *mol,
                                mmdb::Chain *chain_p,
                                const std::string &colour_scheme,
                                const std::string &style,
                                int secondary_structure_usage_flag,
                                const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                const std::vector<std::pair<std::string, int> > &M2T_int_params);
   void add_selection_and_colour(const std::string &sel, const std::string &col);

   std::vector<molecular_triangles_mesh_t>
   get_molecular_triangles_mesh(mmdb::Manager *mol,
                                const std::string &selection_string, // mmdb-format
                                const std::string &colour_scheme,
                                const std::string &style,
                                int secondary_structure_usage_flag,
                                const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                const std::vector<std::pair<std::string, int> > &M2T_int_params);

   //! selection mode is one of 1: residue, 2: sphere, 3: big sphere, 4: chain
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
   get_molecular_triangles_mesh_for_active_residue(int imol, mmdb::Manager *mol, int model_number,
                                                   mmdb::Residue *residue_p,
                                                   const coot::protein_geometry *geom_in,
                                                   int bond_width,
                                                   range_t range,
                                                   int selection_mode);


   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
   get_molecular_triangles_mesh_for_ribbon_with_user_defined_residue_colours(mmdb::Manager *mol, mmdb::Chain *chain_p,
                                                                             const std::vector<std::pair<unsigned int, coot::colour_holder> > &user_defined_colours,
                                                                             int secondary_structure_usage_flag,
                                                                             const std::vector<std::pair<std::string, float> > &M2T_float_params,
                                                                             const std::vector<std::pair<std::string, int> > &M2T_int_params);


};

#endif // GM_INFO_T

