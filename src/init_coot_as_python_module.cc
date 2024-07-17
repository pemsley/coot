/*
 * src/init_coot_as_python_module.cc
 *
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <iostream>
// for stat()
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <Python.h>
#include <gtk/gtk.h>

#include "c-interface.h" // setup_symm_lib()
#include "graphics-info.h"

#include "utils/setup-syminfo.hh"

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
      std::string pkgdatadir = coot::package_data_dir();
      std::string ref_structs_dir = pkgdatadir;
      ref_structs_dir += "/";
      ref_structs_dir += "reference-structures";
      struct stat buf;
      int status = stat(ref_structs_dir.c_str(), &buf);
      if (status != 0) { // file was not found in default location either
	 std::cout << "WARNING:: No reference-structures found (in default location)."
                   << "   " << ref_structs_dir
		   << " and COOT_REF_STRUCTS was not defined." << std::endl;
	 std::cout << "          Ca->Mainchain will not be possible." << std::endl;
      }
   }

}




void
init_coot_as_python_module() {


   std::cout << "###################################################### init_coot_as_python_module() " << std::endl;

   // 20230526-PE this is called by coot_wrap_python() SWIG_init().

   // graphics_info_t::coot_is_a_python_module is set true initially,
   // in main() it is set to false.
   // when we import coot from python, main() is not executed, so we come here
   // with graphics_info_t::coot_is_a_python_module as true

   if (!graphics_info_t::coot_is_a_python_module) return; // coot gui/embedded mode

#ifdef USE_LIBCURL
   curl_global_init(CURL_GLOBAL_NOTHING); // nothing extra (e.g. ssl or WIN32)
#endif

   // 20240717-PE why is this stuff here - in python-related code?
   //       It is so that "import coot" from python init's coot properly

   mmdb::InitMatType();
   setup_syminfo();
   check_reference_structures_dir();
   graphics_info_t g;
   g.use_graphics_interface_flag = 0;
   g.init();

}
