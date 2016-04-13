
#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#ifndef COD_TYPES_HH
#define COD_TYPES_HH

#include <string>
#include "use-rdkit.hh"

#include "bond-table-record-t.hh"
#include "bond-record-container-t.hh"

namespace cod {

   // we need info for neighbours (count) and aromatic status
   class ring_info_t {
      std::vector<int> atom_indices; // mirrors the atomRings() data (type) in RDKit
      bool aromaticity;
   public:
      ring_info_t(const std::vector<int> &ai) { atom_indices = ai; aromaticity = false; }
      ring_info_t(const std::vector<int> &ai, bool aromaticity_in)
      { atom_indices = ai; aromaticity = aromaticity_in; }
      void add_atom(unsigned int i) { atom_indices.push_back(i); }
      void set_aromaticity(bool arom) { aromaticity = arom; }
      unsigned int size() const { return atom_indices.size(); }
      bool get_aromaticity() const { return aromaticity; }
   };

   bond_record_container_t read_acedrg_table(const std::string &file_name);
   
} 

#endif // COD_TYPES_HH

#endif // MAKE_ENHANCED_LIGAND_TOOLS
