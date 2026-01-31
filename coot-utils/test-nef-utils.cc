
#include "nef-utils.hh"

int main(int argc, char **argv) {

   int status = 0;
   if (argc >1) {
      std::string filepath = argv[0];
      nef::NEFParser nef_parser = nef::NEFParser::from_file(filepath);
      status = 1;
   }
   return status;
}
