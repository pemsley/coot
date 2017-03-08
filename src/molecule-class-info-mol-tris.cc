
#include "Python.h"
#include "graphics-info.h"
#include "molecule-class-info.h"

void
molecule_class_info_t::make_molecularrepresentationinstance() {

#ifdef HAVE_CXX11
   if (atom_sel.mol) {
      
      auto mm=new MyMolecule(atom_sel.mol);

      auto sp = std::shared_ptr<MyMolecule>(mm);

      std::cout << "making a molpreinst" << std::endl;
      // should be rep not pre
      auto mri = new MolecularRepresentationInstance;
      mri->create(sp, ColorScheme::colorBySecondaryScheme(), "ALL","Ribbon");
      std::cout << "making a molpreinst done" << std::endl;
      molrepinst=mri;

   }
#endif // AX_CXX_COMPILE_STDCXX_11   
}
