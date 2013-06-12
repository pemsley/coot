
#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "mogul-interface.hh"
#include "c-interface-mogul.hh"

#ifdef USE_GUILE
SCM
mogul_results_scm(const char *mogul_out_file_name) {

   SCM r = SCM_BOOL_F;
   coot::mogul m;
   m.parse(mogul_out_file_name);
   if (m.n_items() > 0) {
      r = SCM_EOL;
      for (unsigned int i=0; i<m.n_items(); i++) { 
	 const coot::mogul_item &item = m[i];
	 // SCM scm_item = scm_list_2(scm_double2num(item.z),
	 // scm_double2num(item.z));
	 SCM scm_item = scm_double2num(item.z);
	 r = scm_cons(scm_item, r);
      }
   }
   return r;
}
#endif 
