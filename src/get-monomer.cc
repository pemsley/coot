

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "graphics-info.h"
#include "cc-interface-scripting.hh"

#include "get-monomer.hh"




// return a new molecule number
int get_monomer_molecule_by_network_and_dict_gen(const std::string &text) {

   int imol = -1;

   if (graphics_info_t::prefer_python) { 
#if defined USE_PYTHON
      std::string python_command("get_pdbe_cif_for_comp_id('");
      python_command += text;
      python_command += "')";
      safe_python_command(python_command);
#endif // USE_PYTHON
   } else { 
#if defined USE_GUILE
      std::string scheme_command("(get-pdbe-cif-for-comp-id \"");
      scheme_command += text;
      scheme_command += "\")";
      SCM r = safe_scheme_command(scheme_command);
      if (scm_is_true(scm_string_p(r))) {
	 scheme_command = "(generate-molecule-from-mmcif \"";
	 scheme_command += text;
	 scheme_command += "\"";
	 scheme_command += " ";
	 scheme_command += "\"";
	 scheme_command += scm_to_locale_string(r);
	 scheme_command += "\")";
	 r = safe_scheme_command(scheme_command);
	 if (scm_is_true(scm_integer_p(r))) {
	    int ir = scm_to_int(r);
	    imol = ir;
	 }
      }
#endif // USE_GUILE
   }

   return imol;
} 
