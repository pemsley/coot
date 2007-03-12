
#include <iostream>
#include "dunbrack.hh"
#include "coot-utils.hh"

int main(int argc, char **argv) {

   CResidue *r = 0; 
   coot::dunbrack d(r);

//    std::string test_string = "45%";
//    std::vector<std::string> parts = coot::util::split_string(test_string, "%");
//    std::cout << test_string << " has " << parts.size() << " parts" << std::endl;
//    test_string = "45";
//    parts = coot::util::split_string(test_string, "%");
//    std::cout << test_string << " has " << parts.size() << " parts" << std::endl;

   if (argc > 1) {
      std::string n(argv[1]);
      std::cout << "Reading Library file: " << n << std::endl;
      d.read_penultimate_library(n);
   } else {
      std::cout << "Usage: " << argv[0] << " dictionary-file-name " << std::endl;
   }


   return 0;
}
