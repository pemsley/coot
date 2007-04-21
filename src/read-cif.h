
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


/*  return 0 if does not have ".cif" extention and 1 if it does. */
/* That is read and display the map */
int try_read_cif_file(const char *filename); 

int try_read_cif_file_and_calc_sfs(const char *filename, int imol);

int try_read_cns_data_file(const char *filename, int imol);


END_C_DECLS
