
#include "utils/coot-utils.hh"
#include "file-system-utils.hh"

int
make_directory_maybe(const std::string &dir) {
   return coot::util::create_directory(dir);
}

