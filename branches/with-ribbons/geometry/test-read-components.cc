

#include <iostream>
#include "protein-geometry.hh"
#include "test-read-components.hh"

int main(int argc, char ** argv) {

   int r = 0;

   if (argc > 1) {
      read_components_file(argv[1]);
   } else {
      std::cout << "Usage: " << argv[0] << " cif-dict-file-name" << std::endl;
   }
   return r;
} 


void read_components_file(std::string cif_dictionary_filename) {

   CMMCIFFile ciffile;
   int ierr = ciffile.ReadMMCIFFile((char *)cif_dictionary_filename.c_str());
   
   if (ierr!=CIFRC_Ok) {
      std::cout << "dirty mmCIF file? " << cif_dictionary_filename.c_str()
		<< std::endl;
      std::cout << "    Bad CIFRC_Ok on ReadMMCIFFile" << std::endl;
   } else {
      std::cout << "There are " << ciffile.GetNofData()
		<< " data in " << cif_dictionary_filename << std::endl; 
      
      for(int idata=0; idata<ciffile.GetNofData(); idata++) {
	 PCMMCIFData data = ciffile.GetCIFData(idata);
	 std::string s = data->GetDataName();
	 // std::cout << idata << " " << s << std::endl;
      }
   }
}

//  1 fails
//  2 fails         500000
//  3 passes 295     50000
//  4 fails         200000
//  5 passes 610    100000
//  6 passes 917    150000
//  7 fails         179948
//  8 passed 1002   164870
//  9 passes 1032   169935
// 10 fails         174884
// 11 passes 1042   171864
// 12 passes 1048   172955
// 13 passes 1053   173860
// 14 passes 1055   174256
// 15 passed 1057   174573
// 16 fails         175060 (we knew that, 10 fails)
// 17 fails         174884 (is 10)
// 18 passes 1057   174573
// -> problematic monomer 773

// 20 fails        1368395

// -------------------
// now for testing 23.cif
//
