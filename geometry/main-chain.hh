
#ifndef MAIN_CHAIN_HH
#define MAIN_CHAIN_HH

#include <string>
#include <mmdb2/mmdb_manager.h>

namespace coot {
   
   // return 0 or 1
   bool is_main_chain_p(mmdb::Atom *at);

   // return 0 or 1
   bool is_hydrogen_p(mmdb::Atom *at);

   // return 0 or 1
   bool is_main_chain_or_cb_p(mmdb::Atom *at);

   // return 0 or 1
   bool is_main_chain_p(const std::string &atom_name);

   // return 0 or 1
   bool is_main_chain_or_cb_p(const std::string &atom_name);

}


#endif // MAIN_CHAIN_HH

