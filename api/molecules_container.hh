#ifndef MOLECULES_CONTAINER_HH
#define MOLECULES_CONTAINER_HH

#include <vector>
#include "coot_molecule.hh"
#include "validation-information.hh"

class molecules_container_t {

   std::vector<coot::molecule_t> molecules;

public:
   molecules_container_t() {}

   bool is_valid_model_molecule(int);
   bool is_valid_map_molecule(int);

   int flipPeptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf);
   int flipPeptide(int imol, const std::string &cid, const std::string &alt_conf);
   int read_pdb(const std::string &file_name);
   int read_mtz(const std::string &file_name, const std::string &f, const std::string &phi, const std::string &weight,
                bool use_weight, bool is_a_difference_map);
   coot::validation_information_t density_fit_analysis(int imol_model, int imol_map);

};

#endif // MOLECULES_CONTAINER_HH
