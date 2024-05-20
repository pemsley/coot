/*
 * utils/setup-syminfo.cc
 *
 * Copyright 2009 by University of Oxford
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include <iostream>
#include <string.h>

#include "coot-utils.hh"
#include "setup-syminfo.hh"
#ifdef _MSC_VER
  // bypass on Windows - thank you Charles Ballard
#else
// for stat()
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

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
      standard_file_name += "/data/";
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
