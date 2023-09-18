#include <string>
#include "Python.h"

void layla_setup_python_basic(int argc, char **argv);
void layla_setup_python_coot_module();
void layla_setup_python_module(const std::string &module_name);
PyObject *layla_safe_python_command_with_return(const std::string &python_cmd);
std::string get_drug_via_wikipedia_and_drugbank_py(const std::string &drugname);