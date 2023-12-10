

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



int gaussian_surface(int imol) {

   auto make_an_ncs_ghost_surface = [] (int imol, mmdb::Manager *mol,
                                        unsigned int i_ch, const std::string &chain_id,
                                        const std::vector<coot::ghost_molecule_display_t> &gi,
                                        bool colour_by_ncs_ghost,
                                        const std::map<std::string, int> &chain_id_map) {

      graphics_info_t g;
      coot::colour_holder ch(0.66, 0.44, 0.44);

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
                                        float grid_scale) {

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
      coot::gaussian_surface_t gauss_surf(mol, chain_id, sigma, contour_level, box_radius, grid_scale);
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
            std::cout << "INFOO:: Calculating Gaussian surface for chain " << chain_p->GetChainID() << std::endl;
            // float sigma=4.4, float contour_level=4.0, float box_radius=5.0, float grid_scale=0.7);
            float sigma = graphics_info_t::gaussian_surface_sigma;
            float contour_level = graphics_info_t::gaussian_surface_contour_level;
            float box_radius = graphics_info_t::gaussian_surface_box_radius;
            float grid_scale = graphics_info_t::gaussian_surface_grid_scale;
            make_an_ncs_chain_surface(imol, mol, chain_p, ncs_chains, sigma, contour_level, box_radius, grid_scale);
         }
      }
      g.graphics_draw();
   }
   return status;
}
