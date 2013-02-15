
#ifndef C_INTERFACE_MOGUL_HH
#define C_INTERFACE_MOGUL_HH
#ifdef USE_GUILE
#include <cstdio> // for std::FILE in gmp.h for libguile.h
#include <libguile.h>		/* for SCM type (returned by safe_scheme_command) */
#endif /*  USE_GUILE */


void mogul_markup(int imol, const char *chain_id, int res_no, const char *ins_code,
		  const char *mogul_out_file_name);

int update_restraints_using_mogul(int imol, const char *chain_id, int res_no, const char *ins_code,
				  const char *residue_type, 
				  const char *mogul_out_file_name);

void set_mogul_max_badness(float b);
float get_mogul_max_badness();

#ifdef __cplusplus
#ifdef USE_GUILE
SCM mogul_results_scm(const char *mogul_out_file_name);
#endif 
#endif 

#endif // C_INTERFACE_MOGUL_HH

