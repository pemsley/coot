/*
 * src/gaussian-surface.cc
 *
 * Copyright 2023 by Medical Research Council
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


#include "cc-interface.hh"
#include "graphics-info.h"

#include "density-contour/gaussian-surface.hh"
#include "c-interface-generic-objects.h"

//! \brief set the sigma for gaussian surface
void set_gaussian_surface_sigma(float s) {
   graphics_info_t::gaussian_surface_sigma = s;
}

//! \brief set the contour_level for gaussian surface
void set_gaussian_surface_contour_level(float s) {
   graphics_info_t::gaussian_surface_contour_level = s;
}

//! \brief set the box_radius for gaussian surface
void set_gaussian_surface_box_radius(float s) {
   graphics_info_t::gaussian_surface_box_radius = s;
}

//! \brief set the grid_scale for gaussian surface
void set_gaussian_surface_grid_scale(float s) {
   graphics_info_t::gaussian_surface_grid_scale = s;
}

//! \brief set the fft B-factor for gaussian surface. Use 0 for no B-factor (default 100)
void set_gaussian_surface_fft_b_factor(float f) {
   graphics_info_t::gaussian_surface_fft_b_factor = f;
}

void set_gaussian_surface_chain_colour_mode(short int mode) {
   graphics_info_t::gaussian_surface_chain_colour_mode = mode;
}

//! \brief set the opacity for a given molecule's gaussian_surface
//!
//! @param imol the molecule index
//! @param opacity between 0. and 1.0
void set_gaussian_surface_opacity(int imol, float opacity) {

   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].gaussian_surface_opacity = opacity;
      g.graphics_draw();
   }
}


#include "c-interface.h" // for first_coords_imol();

void show_gaussian_surface_overlay() {

   auto get_model_molecule_vector = [] () {
                                       graphics_info_t g;
                                       std::vector<int> vec;
                                       int n_mol = g.n_molecules();
                                       for (int i=0; i<n_mol; i++)
                                          if (g.is_valid_model_molecule(i))
                                             vec.push_back(i);
                                       return vec;
                                    };

   GtkWidget *w = widget_from_builder("gaussian_surface_frame");
   GtkWidget *mol_chooser_combobox = widget_from_builder("gaussian_surface_molecule_chooser_combobox");
   GtkWidget *e_sigma          = widget_from_builder("gaussian_surface_sigma_entry");
   GtkWidget *e_radius         = widget_from_builder("gaussian_surface_radius_entry");
   GtkWidget *e_contour_level  = widget_from_builder("gaussian_surface_contour_level_entry");
   GtkWidget *e_b_factor       = widget_from_builder("gaussian_surface_b_factor_entry");
   GtkWidget *e_chain_col_mode = widget_from_builder("gaussian_surface_chain_colour_entry");

   gtk_editable_set_text(GTK_EDITABLE(e_sigma),         coot::util::float_to_string_using_dec_pl(graphics_info_t::gaussian_surface_sigma,         1).c_str());
   gtk_editable_set_text(GTK_EDITABLE(e_radius),        coot::util::float_to_string_using_dec_pl(graphics_info_t::gaussian_surface_box_radius,    1).c_str());
   gtk_editable_set_text(GTK_EDITABLE(e_contour_level), coot::util::float_to_string_using_dec_pl(graphics_info_t::gaussian_surface_contour_level, 2).c_str());
   gtk_editable_set_text(GTK_EDITABLE(e_b_factor),      coot::util::float_to_string_using_dec_pl(graphics_info_t::gaussian_surface_fft_b_factor,  0).c_str());
   gtk_editable_set_text(GTK_EDITABLE(e_chain_col_mode), std::to_string(graphics_info_t::gaussian_surface_chain_colour_mode).c_str());

   graphics_info_t g;
   int imol_active = first_coords_imol();
   auto mv = get_model_molecule_vector();
   g.fill_combobox_with_molecule_options(mol_chooser_combobox, NULL, imol_active,mv);

   gtk_widget_set_visible(w, TRUE);

}


int gaussian_surface(int imol) {

   auto make_an_ncs_ghost_surface = [] (int imol, mmdb::Manager *mol,
                                        unsigned int i_ch, const std::string &chain_id,
                                        const std::vector<coot::ghost_molecule_display_t> &gi,
                                        bool colour_by_ncs_ghost,
                                        const std::map<std::string, int> &chain_id_map) {

      graphics_info_t g;
      coot::colour_holder ch(0.66, 0.44, 0.44);
      float opacity = graphics_info_t::molecules[imol].gaussian_surface_opacity;

      if (colour_by_ncs_ghost) {
         for (const auto &ghost : gi) {
            if (ghost.chain_id == chain_id) {
               const std::string &target_chain_id = ghost.target_chain_id;
               std::map<std::string, int>::const_iterator it = chain_id_map.find(target_chain_id);
               if (it != chain_id_map.end()) {
                  int i_ch_for_target = it->second;
                  ch.rotate_by(0.22 * i_ch_for_target);
               }
            }
         }
      } else {
         ch.rotate_by(0.22 * i_ch);
      }
      glm::vec4 col = colour_holder_to_glm(ch);
      if (opacity != 1.0f) col.a = opacity;

      coot::gaussian_surface_t gauss_surf(mol, chain_id);
      coot::simple_mesh_t smesh = gauss_surf.get_surface();
      std::vector<s_generic_vertex> vertices(smesh.vertices.size());
      for (unsigned int i = 0; i < smesh.vertices.size(); i++) {
         vertices[i] = s_generic_vertex(smesh.vertices[i].pos,
                                        smesh.vertices[i].normal,
                                        smesh.vertices[i].color);
         vertices[i].color = col;
         // std::cout << i << " " << glm::to_string(vertices[i].pos) << "\n";
      }

      g.attach_buffers();

      std::string object_name("Gaussian Surface #");
      object_name += std::to_string(imol);
      object_name += std::string(" Chain ");
      object_name += chain_id;
      int obj_mesh = new_generic_object_number(object_name);
      meshed_generic_display_object &obj = g.generic_display_objects[obj_mesh];
      obj.imol = imol;
      obj.mesh.name = object_name;
      obj.mesh.set_draw_mesh_state(true);
      obj.mesh.import(vertices, smesh.triangles);
      obj.mesh.set_material_specularity(0.7, 128);
      obj.mesh.setup_buffers();
      g.graphics_draw();
   };


   auto make_an_ncs_chain_surface = [] (int imol, mmdb::Manager *mol,
                                        mmdb::Chain  *chain_p,
                                        const std::vector<std::vector<mmdb::Chain *> > &ncs_chains,
                                        float sigma, float contour_level, float box_radius,
                                        float grid_scale, float b_factor) {

      coot::colour_holder ch(0.66, 0.44, 0.44);
      int chain_set_idx = -1;
      for (unsigned int i=0; i<ncs_chains.size(); i++) {
         const auto &vc = ncs_chains[i];
         for (const auto &c : vc) {
            if (c == chain_p) {
               chain_set_idx = i;
               break;
            }
         }
      }
      ch.rotate_by(0.22 * chain_set_idx);

      glm::vec4 col = colour_holder_to_glm(ch);
      std::string chain_id = chain_p->GetChainID();

      // gaussian surface optional args:
      // float sigma=4.4, float contour_level=4.0, float box_radius=5.0, float grid_scale=0.7);
      //
      { // this code block is common
         coot::gaussian_surface_t gauss_surf(mol, chain_id, sigma, contour_level, box_radius, grid_scale, b_factor);
         coot::simple_mesh_t smesh = gauss_surf.get_surface();
         std::vector<s_generic_vertex> vertices(smesh.vertices.size());
         for (unsigned int i = 0; i < smesh.vertices.size(); i++) {
            vertices[i] = s_generic_vertex(smesh.vertices[i].pos,
                                           smesh.vertices[i].normal,
                                           smesh.vertices[i].color);
            vertices[i].color = col;
         }
         graphics_info_t g;
         g.attach_buffers();
         std::string object_name("Gaussian Surface #");
         object_name += std::to_string(imol);
         object_name += std::string(" Chain ");
         object_name += chain_id;
         int obj_mesh = new_generic_object_number(object_name);
         meshed_generic_display_object &obj = g.generic_display_objects[obj_mesh];
         obj.imol = imol;
         obj.mesh.name = object_name;
         obj.mesh.set_draw_mesh_state(true);
         obj.mesh.import(vertices, smesh.triangles);
         obj.mesh.set_material_specularity(0.7, 128);
         obj.mesh.setup_buffers();
      }
   };

   auto make_a_chain_surface = [] (int imol, mmdb::Manager *mol,
                                   int i_ch, const std::string &chain_id,
                                   float sigma, float contour_level, float box_radius,
                                   float grid_scale, float b_factor) {

      coot::colour_holder ch(0.66, 0.44, 0.44);
      ch.rotate_by(0.22 * i_ch);
      glm::vec4 col = colour_holder_to_glm(ch);
      { // this code block is common
         coot::gaussian_surface_t gauss_surf(mol, chain_id, sigma, contour_level, box_radius, grid_scale, b_factor);
         coot::simple_mesh_t smesh = gauss_surf.get_surface();
         std::vector<s_generic_vertex> vertices(smesh.vertices.size());
         for (unsigned int i = 0; i < smesh.vertices.size(); i++) {
            vertices[i] = s_generic_vertex(smesh.vertices[i].pos,
                                           smesh.vertices[i].normal,
                                           smesh.vertices[i].color);
            vertices[i].color = col;
         }
         graphics_info_t g;
         g.attach_buffers();
         std::string object_name("Gaussian Surface #");
         object_name += std::to_string(imol);
         object_name += std::string(" Chain ");
         object_name += chain_id;
         int obj_mesh = new_generic_object_number(object_name);
         meshed_generic_display_object &obj = g.generic_display_objects[obj_mesh];
         obj.imol = imol;
         obj.mesh.name = object_name;
         obj.mesh.set_draw_mesh_state(true);
         obj.mesh.import(vertices, smesh.triangles);
         obj.mesh.set_material_specularity(0.7, 128);
         obj.mesh.setup_buffers();
      }
   };

   int status = 0;
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      int imodel = 1;
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      mmdb::Model *model_p = mol->GetModel(imodel);
      if (model_p) {
         std::vector<std::vector<mmdb::Chain *> > ncs_chains = coot::ncs_related_chains(mol, imodel);
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id(chain_p->GetChainID());
            std::cout << "INFO:: Calculating Gaussian surface for chain " << chain_p->GetChainID()
                      << " with chain-colour mode " << graphics_info_t::gaussian_surface_chain_colour_mode << std::endl;
            // float sigma=4.4, float contour_level=4.0, float box_radius=5.0, float grid_scale=0.7);
            float sigma         = graphics_info_t::gaussian_surface_sigma;
            float contour_level = graphics_info_t::gaussian_surface_contour_level;
            float box_radius    = graphics_info_t::gaussian_surface_box_radius;
            float grid_scale    = graphics_info_t::gaussian_surface_grid_scale;
            float b_factor      = graphics_info_t::gaussian_surface_fft_b_factor;
            if (graphics_info_t::gaussian_surface_chain_colour_mode == 1)
               make_a_chain_surface(imol, mol, ichain, chain_id, sigma, contour_level, box_radius, grid_scale, b_factor);
            else
               make_an_ncs_chain_surface(imol, mol, chain_p, ncs_chains, sigma, contour_level, box_radius, grid_scale, b_factor);
         }
      }
      g.graphics_draw();
   }
   return status;
}
