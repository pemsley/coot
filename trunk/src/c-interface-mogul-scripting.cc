
#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "coot-utils.hh"
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


#ifdef USE_GUILE
SCM mogul_results_process_many(const char *glob_dir, const char *glob_str) {

   SCM r = SCM_EOL;

   std::vector<std::string> files = coot::util::glob_files(glob_dir, glob_str);

   std::map<std::string, bool> done_comp_ids;
   for (unsigned int i=0; i<files.size(); i++) { 
      std::vector<std::string> bits = coot::util::split_string(files[i],  "-");
      if (bits.size() > 4) {
	 std::string comp_id        = bits[3];
	 std::string accession_code = bits[2];

	 // for each accession code, write out the worse z value given
	 // capped sigma values.
	 // 
	 if (1) {
	    coot::mogul m(files[i]);
	    std::pair<float, coot::mogul_item> mzb = m.get_max_z_badness(coot::mogul_item::BOND_OR_ANGLE);
	    if (mzb.first > 0) {
	       std::cout << "   max_z: " << files[i] << " " <<  accession_code << " " << mzb.first
			 << " " << mzb.second << std::endl;
	    }
	    SCM l = scm_list_2(scm_from_locale_string(accession_code.c_str()), scm_double2num(mzb.first));
	    r = scm_cons(l, r);
	 }


	 // Write out the counts and standard deviations for all bonds
	 // (or angles) so that we can see how the data are
	 // distributed:
	 // 
	 if (0) { 
	    std::map<std::string, bool>::const_iterator it = done_comp_ids.find(comp_id);
	    if (it == done_comp_ids.end()) {
	       std::cout << "show results for " << files[i] << std::endl;
	       coot::mogul m(files[i]);
	       for (unsigned int i_item=0; i_item<m.n_items(); i_item++) { 
		  // if (m[i_item].type == coot::mogul_item::BOND) {
		  if (m[i_item].type == coot::mogul_item::ANGLE) {
		     std::cout << "   mogul: " << m[i_item].counts << " " << m[i_item].std_dev
			       <<  "   "  << m[i_item].mean << std::endl;
		  }
	       }
	       done_comp_ids[comp_id] = true;
	    }
	 }
      }
   }

   return r;
}
#endif
