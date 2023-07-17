
#include <stdio.h>
#include "string.h"

#include <gtk/gtk.h>		/* necessary for c-interface.h */

#include "c-interface.h"
#include "read-cif.h"
#include "utils/coot-utils.hh"

#include "read-molecule.hh" // now with std::string args


/*  return 0 if does not have ".cif" or a ".fcf" extention and 1 if it
    does. */
/* That is read and display the map */
int
try_read_cif_file(const char *filename) {

   int imol_new = -1; // bad

   std::string ext = coot::util::file_name_extension(filename);
   if (ext == ".fcf") {
      imol_new = handle_shelx_fcf_file_internal(filename);
   }
   if (ext == ".cif") {
      imol_new = auto_read_cif_data_with_phases(filename);
   }
   return imol_new;

}

int
try_read_cif_file_and_calc_sfs(const char *filename, int imol) {

   std::string ext = coot::util::file_name_extension(filename);
   if (ext == ".cif") {
      read_cif_data(filename, imol);
   }
   return 0;
}



int
try_read_cns_data_file(const char *filename, int imol) {

   printf ("INFO trying to read %s as a CNS/X-PLOR data file\n", filename);
   return handle_cns_data_file(filename, imol);

}
