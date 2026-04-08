
#include "nef-utils.hh"
#include <iostream>

int main(int argc, char **argv) {

   int status = 0;
   if (argc > 1) {
      std::string filepath = argv[1];
      try {
         nef::NEFParser nef_parser = nef::NEFParser::from_file(filepath);
         nef_parser.parse();

         const auto& restraints = nef_parser.restraints();
         std::cout << "\nDistance restraints (first 10):\n";
         std::cout << std::string(70, '-') << "\n";
         size_t count = 0;
         for (const auto& r : restraints) {
            if (count++ >= 10) break;
            std::cout << r.to_string() << "\n";
         }
         std::cout << "\nTotal distance restraints: " << restraints.size() << "\n";
         status = 1;
      } catch (const std::exception& e) {
         std::cerr << "Error: " << e.what() << "\n";
         status = 0;
      }
   } else {
      std::cerr << "Usage: " << argv[0] << " <nef_file>\n";
   }
   return status;
}
