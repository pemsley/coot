
#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "graphics-info.h"
#include "cc-interface.hh" // for pythonize_command_name()
#include "c-interface-scm.hh"

coot::command_arg_t
coot::scripting_function(const std::string &function_name,
			 const std::vector<coot::command_arg_t> &args) {

   coot::command_arg_t r;
   if (graphics_info_t::prefer_python) {
#ifdef USE_PYTHON      
      std::string c = pythonize_command_name(function_name);
      std::vector<std::string> command_strings;
      command_strings.push_back(c);
      for (unsigned int i=0; i<args.size(); i++) { 
	 command_strings.push_back(args[i].as_string());
      }
      std::string s = graphics_info_t::pythonize_command_strings(command_strings);
      PyObject *o = safe_python_command_with_return(s);
      if (o) {
	 if (PyBool_Check(o)) {
	    r.type = coot::command_arg_t::BOOL;
	    r.b = PyInt_AsLong(o);
	 } 
	 if (PyFloat_Check(o)) {
	    r.type = coot::command_arg_t::FLOAT;
	    r.f = PyFloat_AsDouble(o);
	 } 
	 if (PyInt_Check(o)) {
	    r.type = coot::command_arg_t::INT;
	    r.i = PyInt_AsLong(o);
	 } 
	 if (PyString_Check(o)) {
	    r.type = coot::command_arg_t::STRING;
	    r.s = PyString_AsString(o);
	 } 
      } 
#endif      
      
   } else {
#ifdef USE_GUILE
      std::string c = schemize_command_name(function_name);
      std::vector<std::string> command_strings;
      command_strings.push_back(c);
      for (unsigned int i=0; i<args.size(); i++) { 
	 command_strings.push_back(args[i].as_string());
      }
      std::string s = graphics_info_t::schemize_command_strings(command_strings);
      SCM ss = safe_scheme_command(s.c_str());

//       std::cout << "debug:: scripting_function() returns "
//  		<< scm_to_locale_string(display_scm(ss)) 
// 		<< std::endl;
       
      if (scm_is_true(scm_boolean_p(ss))) {
	 r.type = coot::command_arg_t::BOOL;
	 r.b = scm_to_bool(ss);
      }
      if (scm_is_true(scm_integer_p(ss))) {
	 r.type = coot::command_arg_t::INT;
	 r.i = scm_to_int(ss);
      } else { 
	 if (scm_is_true(scm_number_p(ss))) {
	    r.type = coot::command_arg_t::FLOAT;
	    r.f = scm_to_double(ss);
	 }
      }
      if (scm_is_true(scm_string_p(ss))) {
	 r.type = coot::command_arg_t::STRING;
	 r.s = scm_to_locale_string(ss);
      }
#endif 
   } 
   return r;
}
