

// don't forget to git add this

#include "molecule-class-info.h"

void
molecule_class_info_t::add_molecular_symmetry(const clipper::Mat33<double> &mol_symm,
                                              const clipper::Coord_orth &mol_origin) {

   // std::cout << "#### pushing back mol_origin " << mol_origin.format() << std::endl;
   molecular_symmetry_matrices.push_back(std::make_pair(mol_symm, mol_origin));

}
