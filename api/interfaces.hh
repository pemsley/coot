#ifndef INTERFACES_HH
#define INTERFACES_HH

#include <string>
#include "geometry/residue-and-atom-specs.hh"

std::string flipPeptide(const std::string &pdb_file_name_in, const coot::residue_spec_t &rs,
                        const std::string &pdb_file_name_out);

int flipPeptide_mmdb(mmdb::Manager *mol, const coot::residue_spec_t &rs, const std::string &alt_conf);

#endif // INTERFACES_HH
