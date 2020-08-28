

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

#ifdef USE_MOLECULES_TO_TRIANGLES
int
molecule_class_info_t::add_molecular_representation(const std::string &atom_selection,
                                                    const std::string &colour_scheme,
                                                    const std::string &style) {

   gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0]));
   gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
   int status = 0;
   molecular_mesh_generator_t mmg;
   std::string name = "Ribbons";
   Material material;
   material.shininess = 1.5;
   material.specular_strength = 1.12;
   Shader molecular_triangles_shader;
   molecular_triangles_shader.init("moleculestotriangles.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);
   Mesh mesh(name);
   meshes.push_back(mesh);
   meshes.back().import(mmg.get_molecular_triangles_mesh(atom_sel.mol, atom_selection, colour_scheme, style));
   meshes.back().setup(&molecular_triangles_shader, material);

   std::cout << "........... now for molecule " << imol_no << " we have "
             << meshes.size() << " meshes " << std::endl;
   return status;
}

#endif // USE_MOLECULES_TO_TRIANGLES


#ifdef USE_MOLECULES_TO_TRIANGLES
void
molecule_class_info_t::remove_molecular_representation(int idx) {

   if (idx >= 0) {

      // this will shuffle the indices of the other molecule representations, hmm...
      // molrepinsts.erase();
      if (molrepinsts.size() > 0) {
    std::vector<std::shared_ptr<MolecularRepresentationInstance> >::const_iterator it = molrepinsts.end();
    it --;
    molrepinsts.erase(it);
    std::cout << "erased - now molrepinsts size " << molrepinsts.size() << std::endl;
      }
   }
}

#endif // USE_MOLECULES_TO_TRIANGLES
