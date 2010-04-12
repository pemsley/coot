
#include "lbg-molfile.hh"

int
main(int argc, char **argv) {

   if (argc > 1) {
      lig_build::molfile_molecule_t mm;
      mm.read(argv[1]);
      lbg_info_t l;
      l.import(mm);
   }
   return 0;
}
