
#include <stdio.h>
#include "string.h"

#include <gtk/gtk.h>		/* necessary for c-interface.h */

#include "c-interface.h"
#include "read-cif.h"

/*  return 0 if does not have ".cif" or a ".fcf" extention and 1 if it
    does. */
/* That is read and display the map */
int 
try_read_cif_file(const char *filename) { 

  int imol_new = -1; // bad
  char *f_pt; 
  short int is_valid = 0;
  
   f_pt = strrchr(filename, '.');
   if (f_pt != NULL) { 
      if (! (strncmp(f_pt, ".fcf", 4))) { 
	 printf ("INFO trying to read %s as a SHELX fcf file\n", filename);
	 imol_new = handle_shelx_fcf_file_internal(filename);
      }
      if (! (strncmp(f_pt, ".cif", 4))) { 
	 printf ("%s is a mmCIF file\n", filename); 
	 imol_new = auto_read_cif_data_with_phases(filename);
      }
   }
   return imol_new;
}

int 
try_read_cif_file_and_calc_sfs(const char *filename, int imol) { 

   char *f_pt; 
   f_pt = strrchr(filename, '.');
   if (f_pt != NULL) { 
      if (! (strncmp(f_pt, ".cif",4))) { 
	 printf ("%s is a mmCIF file\n", filename); 

	 /* let's call that c-interface function then */
	 
	 read_cif_data(filename, imol);
	 
	 return 1; /* true */
      }
   }
   return 0; 
}



int 
try_read_cns_data_file(const char *filename, int imol) { 

   printf ("INFO trying to read %s as a CNS/X-PLOR data file\n", filename);
   return handle_cns_data_file(filename, imol);
  
}
