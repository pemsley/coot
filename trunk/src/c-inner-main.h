
#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS extern
#define END_C_DECLS     
#endif

#endif /* BEGIN_C_DECLS */

#if defined (USE_PYTHON)
#include "Python.h" 
#endif 

BEGIN_C_DECLS

// Bernhard turns these off for now.
#ifdef UNDERSTAND_GUILE_ON_WINDOWS
/* debugging with CYGWIN voodoo */
#if (defined(_WIN32) || defined(__CYGWIN__))
__declspec(dllimport) extern scm_t_option scm_debug_opts[];
__declspec(dllimport) extern scm_t_option scm_read_opts[];
__declspec(dllimport) extern scm_t_option scm_evaluator_trap_table[];
#endif
#endif

void c_inner_main(void *closure, int argc, char** argv); 
struct command_line_data;
void c_wrapper_scm_boot_guile(int argc, char** argv, struct command_line_data* pcld);
char* does_file_exist (const char     *directory,
		       const char     *filename); 
void start_command_line_python_maybe(char **argv);

END_C_DECLS
