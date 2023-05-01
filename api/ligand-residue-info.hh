#ifndef LIGAND_RESIDUE_INFO_HH
#define LIGAND_RESIDUE_INFO_HH

#include "geometry/residue-and-atom-specs.hh"
#include "coords/Cartesian.h"

namespace coot {

   class ligand_residue_info_t {
   public:
      residue_spec_t residue_spec;
      std::string label;
      coot::Cartesian position;
      std::string residue_name;
   };

   class ligand_residues_container_t {
   public:
      std::vector<ligand_residue_info_t> ligand_residues;
   };

}

#endif // LIGAND_RESIDUE_INFO_HH
