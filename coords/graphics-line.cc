
#include "graphics-line.hh"

void
graphics_line_t::update(mmdb::Atom **atom_selection, int n_atoms) {

   if (atom_index_1 < n_atoms) {
      if (atom_index_1 >= 0) {
         if (atom_index_2 < n_atoms) {
            if (atom_index_2 >= 0) {
               mmdb::Atom *at_1 = atom_selection[atom_index_1];
               mmdb::Atom *at_2 = atom_selection[atom_index_2];
               coot::Cartesian p1(at_1->x, at_1->y, at_1->z);
               coot::Cartesian p2(at_2->x, at_2->y, at_2->z);
               positions = coot::CartesianPair(p1, p2);
            }
         }
      }
   }

}
