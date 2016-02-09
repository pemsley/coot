
void my_wrap_scm_boot_guile(int argc, char** argv);
void start_command_line_python_maybe(char **argv);
void try_load_dot_coot_and_preferences();

#if defined (USE_GUILE) || defined (USE_PYTHON)
extern "C" {
   void SWIG_init();
}
#endif

