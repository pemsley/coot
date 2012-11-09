
#ifndef C_INTERFACE_MOGUL_HH
#define C_INTERFACE_MOGUL_HH
#ifdef USE_GUILE
#include <cstdio> // for std::FILE in gmp.h for libguile.h
#include <libguile.h>		/* for SCM type (returned by safe_scheme_command) */
#endif /*  USE_GUILE */

#ifdef __cplusplus
#ifdef USE_GUILE
SCM mogul_results_scm(const char *mogul_out_file_name);
#endif 
#endif 

#endif // C_INTERFACE_MOGUL_HH

