
#include <iostream>

#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS extern
#define END_C_DECLS     
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS

int mmcif_sfs_to_mtz(const char *cif_file_name, const char *mtz_out_file_name);

END_C_DECLS

int main(int argc, char **argv) {

   int status = 0;

   if (argc != 3) {
      std::cout << "Usage: " << argv[0] << "<input-cif-file-name> <output-mtz-file-name>"<< std::endl;
   } else {

      int success = mmcif_sfs_to_mtz(argv[1], argv[2]);
      if (success != 1)
	 status = 1; // fail.
   } 
   return status;
} 
