
/*  ----------------------------------------------------------------------- */
/*                  scripting                                               */
/*  ----------------------------------------------------------------------- */

// Not including Python.h here beacuse it needs to be first in the sources
// 20210820-PE (but maybe that won't affect it being added here... Hmm)
//
#include <vector>
#include <string>

#ifdef USE_GUILE
SCM safe_scheme_command_test(const char *cmd);
SCM safe_scheme_command(const std::string &scheme_command);
#else 
void safe_scheme_command(const std::string &scheme_command); /* do nothing */
#endif // USE_GUILE

#ifdef USE_PYTHON
void safe_python_command(const std::string &python_command); 
void safe_python_command_by_char_star(const char *python_command);
PyObject *py_clean_internal(PyObject *obj);
PyObject *safe_python_command_with_return(const std::string &python_cmd);
PyObject *safe_python_command_test(const char *cmd);
void safe_python_command_with_unsafe_thread(const char *cmd);
#endif // PYTHON
/*  Is this a repeat of something?  I don't know. */
void run_generic_script(const std::vector<std::string> &cmd_strings);

#ifdef USE_GUILE
SCM coot_has_python_p();
#endif


#ifdef USE_GUILE
SCM test_mol_triangles_scm(SCM i_scm, SCM j_scm);
#endif

#ifdef USE_PYTHON
// maybe a better name when things are more established.
// key can be a single-letter string or an int
void add_key_binding_gtk3_py(PyObject *key, int ctrl_key_flag, PyObject *func,
                             const std::string &description);
void set_light_position_py(int light_id, float x, float y, float z);
#endif

#ifdef USE_GUILE
void add_key_binding_gtk3_scm(int key, int ctrl_key, SCM thunk, const std::string &description);
#endif

void reload_shaders();

