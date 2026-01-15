
#include <Python.h>
#include <string>

class execute_python_results_container_t {
   public:
   PyObject *result;
   std::string error_message;
   std::string stdout;
   execute_python_results_container_t() {
      result = nullptr;
   }
};

execute_python_results_container_t execute_python_code_with_result_internal(const std::string &code);
execute_python_results_container_t execute_python_multiline_code_with_result_internal(const std::string &code);
