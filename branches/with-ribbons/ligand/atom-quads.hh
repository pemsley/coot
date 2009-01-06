
#ifndef LIGAND_ATOM_QUAD_HH
#define LIGAND_ATOM_QUAD_HH

#ifndef HAVE_STRING
#include <string>
#define HAVE_STRING
#endif

namespace coot { 

   // 4 atom names for a torsion:
   class atom_name_quad {
   public:
      std::string atom1;
      std::string atom2;
      std::string atom3;
      std::string atom4;
      atom_name_quad(const std::string &atom_name_1,
                     const std::string &atom_name_2,
                     const std::string &atom_name_3,
                     const std::string &atom_name_4) {
         atom1 = atom_name_1;
         atom2 = atom_name_2;
         atom3 = atom_name_3;
         atom4 = atom_name_4;
      }
   };

   class atom_index_quad {
   public:
      int index1;
      int index2;
      int index3;
      int index4;
      atom_index_quad() { // unassigned pair
         index1 = -1;
         index2 = -1;
         index3 = -1;
         index4 = -1;
      }
      atom_index_quad(int i1, int i2, int i3, int i4) {
         index1 = i1;
         index2 = i2;
         index3 = i3;
         index4 = i4;
      }
   };

}

#endif // LIGAND_ATOM_QUAD_HH

