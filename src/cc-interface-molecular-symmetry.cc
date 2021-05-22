

// don't forget to git add this file

#ifdef USE_PYTHON  // fixes Python include file order problems
#include "Python.h"
#endif

#include "cc-interface.hh"
#include "c-interface.h"

#include "graphics-info.h"

void add_molecular_symmetry(int imol,
                            double r_00, double r_01, double r_02,
                            double r_10, double r_11, double r_12,
                            double r_20, double r_21, double r_22,
                            double about_origin_x,
                            double about_origin_y,
                            double about_origin_z) {

   if (is_valid_model_molecule(imol)) {
      clipper::Coord_orth molecule_origin(about_origin_x, about_origin_y, about_origin_z);
      clipper::Mat33<double> mol_symm_matrix(r_00, r_01, r_02,
                                             r_10, r_11, r_12,
                                             r_20, r_21, r_22);
      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      m.add_molecular_symmetry(mol_symm_matrix, molecule_origin);
   }

}
