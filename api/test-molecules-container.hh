#ifndef TEST_MOLECULES_CONTAINER_HH
#define TEST_MOLECULES_CONTAINER_HH

#include <string>
#include <mmdb2/mmdb_manager.h>
#include "coords/Cartesian.h"

void starting_test(const char *func);
std::string reference_data(const std::string &file);
coot::Cartesian atom_to_cartesian(mmdb::Atom *at);

#endif // TEST_MOLECULES_CONTAINER_HH
