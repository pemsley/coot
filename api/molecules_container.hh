#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <vector>
#include "coot_molecule.hh"

class molecules_container_t {

   std::vector<coot_molecule_t> molecules;

public:
   molecules_container_t() {}

   bool is_valid_model_molecule(int);
   bool is_valid_map_molecule(int);

   int flipPeptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf);
   int read_pdb(const std::string &file_name);
};

#endif // MOLECULES_CONTAINER_HH
