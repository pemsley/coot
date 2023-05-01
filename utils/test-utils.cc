
#include "coot-utils.hh"

int main(int argc, char **argv) {

   if (true) {

      std::string file_name = "test.pdb";
      std::string c = "copied.pdb";
      coot::copy_file(file_name, c);

      // now read both files and check that they are the same.

   }

   return 0;
}
