
#ifndef ACEDRG_TYPES_FOR_RESIDUE_HH
#define ACEDRG_TYPES_FOR_RESIDUE_HH

#include <string>
#include <vector>

namespace coot {

   class acedrg_types_for_bond_t {
   public:
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_type_1;
      std::string atom_type_2;
      double bond_length;
      acedrg_types_for_bond_t() : bond_length(-1) {}
      acedrg_types_for_bond_t(const std::string &a11, const std::string &a12, const std::string &a21, const std::string &a22,
                              double bl) : atom_id_1(a11), atom_id_2(a12), atom_type_1(a21), atom_type_2(a22),
                                           bond_length(bl) {}
   };

   class acedrg_types_for_angle_t {
   public:
      std::string atom_id_1;
      std::string atom_id_2;
      std::string atom_id_3;
      std::string atom_type_1;
      std::string atom_type_2;
      std::string atom_type_3;
      double angle; // radians
      acedrg_types_for_angle_t(const std::string &a11, const std::string &a12, const std::string &a13,
                               const std::string &a21, const std::string &a22, const std::string &a23,
                               double a) : atom_id_1(a11), atom_id_2(a12), atom_id_3(a13),
                                           atom_type_1(a21), atom_type_2(a22), atom_type_3(a23),
                                           angle(a) {}
   };

   class acedrg_types_for_residue_t {
   public:
      std::vector<acedrg_types_for_bond_t> bond_types;
   };

}

#endif // ACEDRG_TYPES_FOR_RESIDUE_HH
