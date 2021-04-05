
#include <iostream>
// for stat()
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <Python.h>
#include <gtk/gtk.h>

#include "c-interface.h" // setup_symm_lib()
#include "graphics-info.h"

// put these functions in init_coot_as_python_module.cc

#include "coot-surface/rgbreps.h"

void check_reference_structures_dir() {

   char *coot_reference_structures = getenv("COOT_REF_STRUCTS");
   if (coot_reference_structures) {
      struct stat buf;
      int status = stat(coot_reference_structures, &buf);
      if (status != 0) { // file was not found in default location either
	 std::cout << "WARNING:: The reference structures directory "
		   << "(COOT_REF_STRUCTS): "
		   << coot_reference_structures << " was not found." << std::endl;
	 std::cout << "          Ca->Mainchain will not be possible." << std::endl;
      }
   } else {

      // check in the default place: pkgdatadir = $prefix/share/coot
      std::string pkgdatadir = PKGDATADIR;
      std::string ref_structs_dir = pkgdatadir;
      ref_structs_dir += "/";
      ref_structs_dir += "reference-structures";
      struct stat buf;
      int status = stat(ref_structs_dir.c_str(), &buf);
      if (status != 0) { // file was not found in default location either
	 std::cout << "WARNING:: No reference-structures found (in default location)."
		   << "          and COOT_REF_STRUCTS was not defined." << std::endl;
	 std::cout << "          Ca->Mainchain will not be possible." << std::endl;
      }
   }

}


void setup_symm_lib() {

   // How should the setting up of symmetry work?
   //
   // First we check the environment variable SYMINFO.
   //
   // If that is not set, then we look in the standard (hardwired at
   // compile time) place
   //
   // If both these fail then give an error message.

   char *syminfo = getenv("SYMINFO");
   if (!syminfo) {

      std::string symstring("SYMINFO=");

      // using PKGDATADIR will work for those who compiler, not the
      // binary users:
      std::string standard_file_name = PKGDATADIR; // xxx/share/coot
      standard_file_name += "/";
      standard_file_name += "syminfo.lib";

      struct stat buf;
      int status = stat(standard_file_name.c_str(), &buf);
      if (status != 0) { // standard-residues file was not found in default location

	 // This warning is only sensible for those who compile (or
	 // fink).  So let's test if SYMINFO was set before we write it
	 //
	 std::cout << "WARNING:: Symmetry library not found at "
		   << standard_file_name
		   << " and environment variable SYMINFO is not set." << std::endl;
	 std::cout << "WARNING:: Symmetry will not be possible\n";

      } else {

	 symstring += standard_file_name;

	 // Mind bogglingly enough, the string is part of the environment
	 // and malleable, so const char * of a local variable is not
	 // what we want at all.
	 //
	 //  We fire and forget, we don't want to change s.
	 //
	 char * s = new char[symstring.length() + 1];
	 strcpy(s, symstring.c_str());
	 putenv(s);
	 // std::cout << "DEBUG:: SYMINFO set/not set? s is " << s <<std::endl;
      }
   }
}


void
init_coot_as_python_module() {

   // graphics_info_t::coot_is_a_python_module is set true initially,
   // in main() it is set to false.
   // when we import coot from python, main() is not executed, so we come here
   // with graphics_info_t::coot_is_a_python_module as true

   if (!graphics_info_t::coot_is_a_python_module) return; // coot gui/embedded mode

#ifdef USE_LIBCURL
   curl_global_init(CURL_GLOBAL_NOTHING); // nothing extra (e.g. ssl or WIN32)
#endif

   setup_symm_lib();
   check_reference_structures_dir();
   graphics_info_t g;
   g.use_graphics_interface_flag = 0;
   g.init();

}
   
