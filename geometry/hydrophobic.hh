#ifndef HYDROPHOBIC_HH
#define HYDROPHOBIC_HH

#include <string>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   bool is_hydrophobic_atom(const std::string &residue_name, const std::string &atom_name);

   bool is_hydrophobic_atom(mmdb::Atom *at);
}

#endif
