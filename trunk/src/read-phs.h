
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


int 
try_read_phs_file(const char *filename); 


int phs_pdb_cell_symm(); 

void test_read_coords(const gchar *filename); 

int phs_pdb_cell_symm();


void do_phs_cell_choice_window(); 

END_C_DECLS
