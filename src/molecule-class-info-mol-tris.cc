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
   std::string name = atom_selection + " " + colour_scheme + " " + style;
   Material material;

   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: add_molecular_representation() pos-B " << err << std::endl;

   material.do_specularity = true;        // 20210905-PE make these user settable. Perhaps they are? I should check.
   material.shininess = 256.0;
   material.specular_strength = 0.56;

   // if (colour_scheme == "Rainbow") {
   if (colour_scheme == "colorRampChainsScheme") {

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
               mesh.set_name(atom_selection + " " + colour_scheme + " Rainbow Ribbons");
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
