/*
 * src/molecule-class-info-mol-tris.cc
 *
 * Copyright 2017 by Medical Research Council
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

#ifdef USE_PYTHON
#include "Python.h"
#endif

#include <algorithm>

#include "graphics-info.h"
#include "molecule-class-info.h"

// make and add to the scene
int
molecule_class_info_t::make_molecularrepresentationinstance(const std::string &atom_selection,
                                                            const std::string &colour_scheme,
                                                            const std::string &style) {

   int status = 0;
   return status;

}

//   std::vector<std::pair<std::string, float> > M2T_float_params;
//   std::vector<std::pair<std::string, float> > M2T_int_params;

//! Update float parameter for MoleculesToTriangles molecular mesh
void
molecule_class_info_t::M2T_updateFloatParameter(const std::string &param_name, float value) {

   M2T_float_params.push_back(std::make_pair(param_name, value));
}

//! Update int parameter for MoleculesToTriangles molecular mesh
void
molecule_class_info_t::M2T_updateIntParameter(const std::string &param_name, int value) {

   M2T_int_params.push_back(std::make_pair(param_name, value));
}



void
molecule_class_info_t::set_mol_triangles_is_displayed(int state) {

#ifdef USE_MOLECULES_TO_TRIANGLES
   if (molrepinsts.size()) {
      if (state) {
         for (auto mri : molrepinsts)
            graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(mri);
      } else {
         for (auto mri : molrepinsts)
            graphics_info_t::mol_tri_scene_setup->removeRepresentationInstance(mri);
      }
   }
#endif // USE_MOLECULES_TO_TRIANGLES

}

#include "molecular-mesh-generator.hh"
#include "MoleculesToTriangles/CXXClasses/tubes.hh"
int
molecule_class_info_t::add_molecular_representation(const std::string &atom_selection,
                                                    const std::string &colour_scheme,
                                                    const std::string &style,
                                                    int secondary_structure_usage_flag) {
   int status = 0;

   if (false)
      std::cout << "DEBUG:: in mcit::add_molecular_representation() atom_selection: \"" << atom_selection << "\""
                << " colour_scheme: \"" << colour_scheme << "\" style: \"" << style << "\"" << std::endl;

   GLenum err = glGetError();
   if (err)
      std::cout << "GL ERROR:: add_molecular_representation() --- start --- " << err << std::endl;

   if (! atom_sel.mol)  return 0;
   if (atom_sel.n_selected_atoms == 0)  return 0;

   gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0])); // needed?
   gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
   molecular_mesh_generator_t mmg;

   // user-facing label shown against this representation in the Display Manager -
   // for Worms/Bendix/Tube (helices), match the Draw > Molecule menu item name
   // exactly, rather than the generic "style: colour_scheme_label" used below.
   std::string name;
   if (style == "Tubes" && secondary_structure_usage_flag == 1) {
      name = "Worms";
   } else if (style == "Tubes") {
      name = "Bendix";
   } else if (style == "TubeHelices") {
      name = "Tube";
   } else {
      std::string colour_scheme_label = colour_scheme;
      if (colour_scheme == "Chain" || colour_scheme == "colorChainsScheme") colour_scheme_label = "By Chain";
      if (colour_scheme == "colorRampChainsScheme")                        colour_scheme_label = "Rainbow";
      if (colour_scheme == "colorBySecondaryScheme" || colour_scheme == "Secondary") colour_scheme_label = "Sec. Struct.";
      if (colour_scheme == "colorByElementScheme" || colour_scheme == "Element")     colour_scheme_label = "By Element";
      name = style + ": " + colour_scheme_label;
   }

   // identifies the "slot" that this representation occupies (independent of the
   // display label above), so that re-requesting the exact same representation replaces
   // itself rather than being added on top of its old self - without this, two
   // near-identical overlapping meshes z-fight and the new colouring never becomes
   // visible on screen. Worms and Bendix are both style "Tubes" (they differ only by
   // secondary_structure_usage_flag) but are meant to coexist as separate, independently
   // toggleable representations rather than replacing each other, so the flag is folded
   // into the key for that one style to keep their slots distinct.
   std::string style_key = style;
   if (style == "Tubes") style_key += std::to_string(secondary_structure_usage_flag);
   std::string representation_key = atom_selection + "\x1f" + style_key;
   meshes.erase(std::remove_if(meshes.begin(), meshes.end(),
      [&representation_key] (const Mesh &m) { return m.representation_key == representation_key; }),
      meshes.end());

   Material material;

   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: add_molecular_representation() pos-B " << err << std::endl;

   material.do_specularity = true;        // 20210905-PE make these user settable. Perhaps they are? I should check.
   material.shininess = 256.0;
   material.specular_strength = 0.56;

   if (style == "Tubes") { // bendy-helix "worm" representation

      // "Worms" (DONT_USE, flag 1 - no SSE computed, so the whole backbone is one
      // uncomputed run) want a thick uniform tube; "Bendix" (flag 0 or 2 - helices
      // drawn as cylinders) wants a thin coil to match ribbonStyleCoilThickness's
      // default in Ribbon mode, so the coil between helices doesn't dominate.
      float radius_for_coil = (secondary_structure_usage_flag == 1) ? 0.8f : 0.3f;
      int Cn_for_coil = 2;
      int accuracy_for_coil = 12;
      unsigned int n_slices_for_coil = 12;
      // MoleculesToTriangles can't resolve the installed data directory itself (it can't
      // depend on coot-utils), so resolve the helix reference template's absolute path
      // here instead - otherwise it only works when coot happens to be run from a
      // directory that happens to contain theor-helix-z-ori-v2.pdb.
      std::string helix_template_pdb_file_name = coot::package_data_dir() + "/theor-helix-z-ori-v2.pdb";
      coot::m2t::simple_mesh_t tubes_mesh =
         make_tubes_representation(atom_sel.mol, atom_selection, colour_scheme, radius_for_coil, Cn_for_coil,
                                   accuracy_for_coil, n_slices_for_coil, secondary_structure_usage_flag,
                                   helix_template_pdb_file_name);

      // tubes_mesh is coot::m2t::simple_mesh_t (MoleculesToTriangles can't depend on
      // coot-utils, so it uses its own mesh_vertex_t/mesh_triangle_t) - convert to the
      // s_generic_vertex/g_triangle that Mesh is built from everywhere else in this function.
      std::vector<s_generic_vertex> vertices;
      vertices.reserve(tubes_mesh.vertices.size());
      for (const auto &v : tubes_mesh.vertices)
         vertices.push_back(s_generic_vertex(v.pos, v.normal, v.color));
      std::vector<g_triangle> triangles;
      triangles.reserve(tubes_mesh.triangles.size());
      for (const auto &t : tubes_mesh.triangles)
         triangles.push_back(g_triangle(t.point_id[0], t.point_id[1], t.point_id[2]));

      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > verts_and_tris(vertices, triangles);
      Mesh mesh(verts_and_tris);
      mesh.set_name(name);
      mesh.set_representation_key(representation_key);
      meshes.push_back(mesh);
      meshes.back().setup(material);

   } else if (style == "TubeHelices") {

      // A normal Ribbon representation (strand = arrow, coil = thin tube, via the
      // same drawRibbon() everything else uses) with its helix geometry suppressed
      // (hideHelixGeometry), merged with a straight-cylinder helix mesh (one PCA-fit
      // axis per helix - see make_straight_cylinder_helices_mesh()) to fill the gap -
      // so helices come out as smooth rods instead of drawRibbon()'s native per-residue
      // elliptical sweep, while strand/coil keep their normal cartoon look untouched.

      // local copy: must not mutate the molecule's persistent M2T_int_params, or every
      // future plain Ribbon representation on this molecule would also lose its helices.
      std::vector<std::pair<std::string, int> > local_int_params = M2T_int_params;
      local_int_params.push_back(std::make_pair(std::string("hideHelixGeometry"), 1));

      std::vector<molecular_triangles_mesh_t> mtm =
         mmg.get_molecular_triangles_mesh(atom_sel.mol, atom_selection, colour_scheme, "Ribbon",
                                          secondary_structure_usage_flag,
                                          M2T_float_params, local_int_params);
      molecular_triangles_mesh_t meshes_together;
      for (unsigned int i=0; i<mtm.size(); i++)
         meshes_together.add_to_mesh(mtm[i].vertices, mtm[i].triangles);

      std::vector<s_generic_vertex> vertices = meshes_together.vertices;
      std::vector<g_triangle> triangles = meshes_together.triangles;

      float radius_for_helices = 2.5;
      unsigned int n_slices_for_helices = 16;
      coot::m2t::simple_mesh_t helix_mesh =
         make_straight_cylinder_helices_mesh(atom_sel.mol, atom_selection, radius_for_helices,
                                             n_slices_for_helices, secondary_structure_usage_flag);

      // offset the helix mesh's (0-based) triangle indices so they index correctly
      // into the combined vertex buffer once appended after the ribbon's vertices, then
      // convert from m2t's mesh_vertex_t/mesh_triangle_t to s_generic_vertex/g_triangle.
      unsigned int idx_base = vertices.size();
      vertices.reserve(vertices.size() + helix_mesh.vertices.size());
      for (const auto &v : helix_mesh.vertices)
         vertices.push_back(s_generic_vertex(v.pos, v.normal, v.color));
      triangles.reserve(triangles.size() + helix_mesh.triangles.size());
      for (auto t : helix_mesh.triangles) {
         t.rebase(idx_base);
         triangles.push_back(g_triangle(t.point_id[0], t.point_id[1], t.point_id[2]));
      }

      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > verts_and_tris(vertices, triangles);
      Mesh mesh(verts_and_tris);
      mesh.set_name(name);
      mesh.set_representation_key(representation_key);
      meshes.push_back(mesh);
      meshes.back().setup(material);

   } else if (colour_scheme == "colorRampChainsScheme") {

      std::cout << "---------------------------------------  Rainbow ----------------------" << std::endl;
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         molecular_triangles_mesh_t meshes_together;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            if (n_res > 1) {
               std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > verts_and_tris =
                  mmg.get_molecular_triangles_mesh(atom_sel.mol, chain_p, colour_scheme, style,
                                                   secondary_structure_usage_flag,
                                                   M2T_float_params, M2T_int_params);
               Mesh mesh(verts_and_tris);
               mesh.set_name(name);
               mesh.set_representation_key(representation_key);
               meshes.push_back(mesh);
               meshes.back().setup(material); // do I need the shader to do this!?
            }
         }
      }

   } else {

      if (false)
         std::cout << "DEBUG:: in mcit::add_molecular_representation() atom_selection: \"" << atom_selection << "\""
                   << " colour_scheme: \"" << colour_scheme << "\" style: \"" << style << "\""
                   << " non-colour-ramp path " << std::endl;

      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: add_molecular_representation() non-colour-ramp-path " << err << std::endl;

      std::vector<molecular_triangles_mesh_t> mtm =
         mmg.get_molecular_triangles_mesh(atom_sel.mol, atom_selection, colour_scheme, style,
                                          secondary_structure_usage_flag,
                                          M2T_float_params, M2T_int_params);

      // Mesh mesh(mtm);
      // meshes.push_back(mesh);
      // meshes.back().setup(&molecular_triangles_shader, material);

      {
         // hacketty hack! This is to make the "old" mesh method work again in Coot, without
         // moving to the Model method
         molecular_triangles_mesh_t meshes_together;
         for (unsigned int i=0; i<mtm.size(); i++) {
            if (false)
               std::cout << "meshes_together " << i << " " << mtm[i].vertices.size() << " "
                         << mtm[i].triangles.size() << std::endl;
            meshes_together.add_to_mesh(mtm[i].vertices, mtm[i].triangles);
         }
         std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
            meshes_together_pair(meshes_together.vertices, meshes_together.triangles);
         Mesh mesh(meshes_together_pair);
         mesh.set_name(name);
         mesh.set_representation_key(representation_key);
         meshes.push_back(mesh);
         // meshes.back().setup(&molecular_triangles_shader, material); 20210910-PE
         meshes.back().setup(material);
         // meshes.back().debug_to_file();
      }
   }

   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: add_molecular_representation() --- end --- " << err << std::endl;

   return status;
}



void
molecule_class_info_t::remove_molecular_representation(int idx) {

   if (idx >= 0) {

      // this will shuffle the indices of the other molecule representations, hmm...
      // molrepinsts.erase();
      if (molrepinsts.size() > 0) {
         std::vector<std::shared_ptr<MolecularRepresentationInstance> >::iterator it = molrepinsts.end();
         it --;
         molrepinsts.erase(it);
         std::cout << "erased - now molrepinsts size " << molrepinsts.size() << std::endl;
      }
   }

}
