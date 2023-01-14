

#include "cc-interface.hh"
#include "graphics-info.h"

#include "density-contour/gaussian-surface.hh"
#include "c-interface-generic-objects.h"

int gaussian_surface(int imol) {

   auto make_a_surface = [] (int imol, mmdb::Manager *mol,
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
      object_name += " Chain ";
      object_name += chain_id;
      int obj_mesh = new_generic_object_number(object_name);
      meshed_generic_display_object &obj = g.generic_display_objects[obj_mesh];
      obj.mesh.name = object_name;
      obj.mesh.set_draw_mesh_state(true);
      obj.mesh.import(vertices, smesh.triangles);
      obj.mesh.set_material_specularity(0.7, 128);
      obj.mesh.setup_buffers();
      g.graphics_draw();
   };


   int status = 0;
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {

      std::vector<coot::ghost_molecule_display_t> gi = g.molecules[imol].NCS_ghosts();

      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      std::vector<std::string> chain_ids = g.molecules[imol].get_chain_ids();
      std::map<std::string, int> chain_id_map;
      for (unsigned int i_ch=0; i_ch<chain_ids.size(); i_ch++) {
         auto chain_id = chain_ids[i_ch];
         chain_id_map[chain_id] = i_ch;
      }

      if (false) {
         std::cout << ":::::::::::::::::::::::: ghosts " << gi.size() << std::endl;
         for (unsigned int ighost=0; ighost<gi.size(); ighost++) {
            const auto &ghost = gi[ighost];
            std::cout << "   " << ighost << " " << ghost.name << " " << ghost.chain_id << " " << ghost.target_chain_id
                      << "\n" << ghost.rtop.format() << std::endl;
         }
      }

      for (unsigned int i_ch=0; i_ch<chain_ids.size(); i_ch++) {
         auto chain_id = chain_ids[i_ch];
         // gi is used if colour_by_ncs_ghost is true
         bool colour_by_ncs_ghost = false;
         make_a_surface(imol, mol, i_ch, chain_id, gi, colour_by_ncs_ghost, chain_id_map);
         status = 1;
      }
   }
   return status;
}
