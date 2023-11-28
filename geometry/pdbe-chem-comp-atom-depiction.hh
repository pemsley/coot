#ifndef PDBE_CHEM_COMP_ATOM_DEPICTION_HH
#define PDBE_CHEM_COMP_ATOM_DEPICTION_HH

#include <vector>
#include <string>

namespace coot {

   class depiction_atom_t {
   public:
      depiction_atom_t(const std::string &atom_name, const std::string &element, double xx, double yy, int idx) :
         atom_name(atom_name), element(element), x(xx), y(yy), pdbe_ordinal(idx) {}
      std::string atom_name;
      std::string element;
      double x, y;
      int pdbe_ordinal;
   };

   class chem_comp_atom_depiction_t {
      std::string comp_id;
   public:
      chem_comp_atom_depiction_t() : comp_id("unset") {}
      chem_comp_atom_depiction_t(const std::string &comp_id_in, const std::vector<depiction_atom_t> &atoms_in) :
         comp_id(comp_id_in), atoms(atoms_in) {}
      std::vector<depiction_atom_t> atoms;
      bool empty() const { return atoms.empty(); }
   };

}

#endif // PDBE_CHEM_COMP_ATOM_DEPICTION_HH
