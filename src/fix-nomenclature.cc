
#include <stdlib.h>
#include <iostream> // fixes undefined strchr, strchrr problems

#include "geometry/protein-geometry.hh"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.hh"

#include "coot-nomenclature.hh"

int
main(int argc, char **argv) {

//    int a = 5;
//    std::cout << a << "\n";
//    a = a | 4; 
//    std::cout << a << "\n";
//    a = a | 8; 
//    std::cout << a << "\n";
//    a = a | 16; 
//    std::cout << a << "\n";
//    a = a | 16; 
//    std::cout << a << "\n";

//    a = 1;
//    a = a << 1;
//    std::cout << a << "\n";
//    a = a << 1;
//    std::cout << a << "\n";
//    a = a <<= 1;
//    std::cout << a << "\n";


   if (argc < 3) {
      std::cout << "Usage: " << argv[0] << " in-filename out-filename\n";
      exit(1);
   } else {
      coot::protein_geometry geom;
      std::string filename = argv[1];
      atom_selection_container_t asc = get_atom_selection(filename, true, false, false);
      coot::nomenclature n(asc.mol);
      asc.mol->WritePDBASCII(argv[2]);
      n.fix(&geom);
   }
   return 0;
}
