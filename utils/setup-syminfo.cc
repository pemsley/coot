
#include <iostream>
#include <string.h>

#include "coot-utils.hh"
#include "setup-syminfo.hh"
// for stat()
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void setup_syminfo() {

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
      std::string standard_file_name = coot::package_data_dir();
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
	 // We fire and forget, we don't want to change s.
	 //
         int ll = symstring.length();
	 char * s = new char[ll + 1];
         s[ll] = 0;
	 strcpy(s, symstring.c_str());
	 putenv(s);
	 // std::cout << "DEBUG:: SYMINFO set/not set? s is " << s <<std::endl;
      }
   }
}
