
#include "coot-coord-utils.hh"

int main(int argc, char **argv) {

   if (argc > 1) {
      std::string file_name = argv[1];
      CMMDBManager *t_mol = new CMMDBManager;
      int status = t_mol->ReadPDBASCII(file_name.c_str());
      if (status != Error_NoError) {
	 std::cout << "ERROR:: on reading " << file_name << std::endl;
      } else {
      }
   }
   return 0;
} 
