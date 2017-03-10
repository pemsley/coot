
#include "Python.h"
#include "graphics-info.h"
#include "molecule-class-info.h"

// make and add to the scene
void
molecule_class_info_t::make_molecularrepresentationinstance() {

#ifdef USE_MOLECULES_TO_TRIANGLES
#ifdef HAVE_CXX11
   if (atom_sel.mol) {

      auto sp = std::make_shared<MyMolecule> (atom_sel.mol);

      std::cout << "making a molrepinst" << std::endl;

      auto cs = ColorScheme::colorBySecondaryScheme();
      molrepinst = MolecularRepresentationInstance::create(sp, cs, "ALL", "Ribbon");

      std::cout << "made a molrepinst " << molrepinst << std::endl;

      graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(molrepinst);

   }
#endif // HAVE_CXX11
#endif // USE_MOLECULES_TO_TRIANGLES
}
