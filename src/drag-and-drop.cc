/* src/drag-and-drop.cc
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2013 by Medical Research Council
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
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <iostream>
#include <gtk/gtk.h>
#include "drag-and-drop.hh"

#ifdef USE_GUILE   
#include <cstdio> // for std::FILE in gmp.h for libguile.h
#include <libguile.h>		/* for SCM type (returned by safe_scheme_command) */
#endif

#include "cc-interface.hh"
#include "c-interface.h"
#include "read-molecule.hh" // now with std::string args
#include "curl-utils.hh"


// why doesn't this live in a header?
void network_get_accession_code_entity(const std::string &text, int mode);

//! \brief handle the string that get when an URL (or text from an url) is dropped.
//
// uri string can be a concatenation of string with \ns between them
// (and a \n to end)
// 
int handle_drag_and_drop_string(const std::string &uri_in) {

   int handled = FALSE;
   bool tried_already = false;
   std::string uri = uri_in;
   std::string url = uri_in;

   if (true)
      std::cout << "::::::::::::::::::::::::::: handle_drag_and_drop_string: " << uri_in << " " << std::endl;

   if (! tried_already) {
      // OK, was it an HTTP type string?
      if (url.length() > 9) {

         if (url.substr(0,7) == "http://" || url.substr(0,8) == "https://") {
            tried_already = true;
            int l = url.length();
            if (url[l-1] == '\n') {
               // std::cout << "extra \\n" << std::endl;
               url = url.substr(0, l-1);
            }

            l = url.length();
            int c = url[l-1];
            if (url[l-1] == '\r') {
               // std::cout << "extra \\r" << std::endl;
               url = url.substr(0, l-1);
            }

            if (true) { // OK, we made it (or had it)
               std::string url_file_name_file = url;

               std::string ext = coot::util::file_name_extension(url);

               if (ext == ".png") {
                  // special rule - convert the url of the png to that
                  // of an accession code.
                  if (url.find("/PDBimages/")) {
                     std::pair<std::string, std::string> s   = coot::util::split_string_on_last_slash(url);
                     std::pair<std::string, std::string> ss  = coot::util::split_string_on_last_slash(s.first);
                     std::pair<std::string, std::string> sss = coot::util::split_string_on_last_slash(ss.first);
                     tried_already = true;
                     handled = FALSE;
                     if (ss.second.length() == 2) {
                        if (sss.second.length() == 2) {
                           std::string code;
                           code += ss.second[0];
                           code += sss.second;
                           code += ss.second[1];
                           // get_coords_for_accession_code(code.c_str());
                           network_get_accession_code_entity(code, 0);
                        }
                     }
                  }

               } else {

                  // it was coords or mtz - we presume
                  std::string::size_type pos = url.find_last_of('/');
                  if (pos != std::string::npos) {
                     // normal path
                     url_file_name_file = url.substr(pos);
                  }
                  // use XDG BDP here.
                  std::string file_name = coot::util::append_dir_file("coot-download", url_file_name_file);
                  // return 0 on success
                  int status = coot_get_url(url.c_str(), file_name.c_str());
                  handled = handle_drag_and_drop_single_item(file_name);
               }
            }
         }
      }
   }

   if (! tried_already) {
      // was it a 4-letter (accession number)?
      int l = uri_in.length();
      if (l == 4) {
         get_coords_for_accession_code(uri_in.c_str());
         tried_already = true;
         handled = TRUE;
      }
   }

   if (! tried_already) {
      if (coot::file_exists(url)) {
	 handled = handle_drag_and_drop_single_item(url);
         tried_already = true;
      }
   }

   if (! tried_already) {
      if (uri.length() > 7) {
         if (uri.find("file:///") != std::string::npos) {
            std::string fn = uri.substr(7);
            std::string ext = coot::util::file_name_extension(fn);
            if (ext == ".cif") read_coordinates(fn);
            if (ext == ".pdb") read_coordinates(fn);
            if (ext == ".mtz") auto_read_make_and_draw_maps(fn.c_str());
         }
      }
   }
   return handled;
}

int handle_drag_and_drop_single_item(const std::string &file_name) {

   int handled = FALSE;
   std::cout << "Here in handle_drag_and_drop_single_item() with file_name "
             << file_name << std::endl;

   std::string ext = coot::util::file_name_extension(file_name);
   if (ext == ".cif") {
      // try as restraints file
      int n_bonds = read_cif_dictionary(file_name.c_str());
      if (n_bonds > 0) {
	 handled = TRUE;
      }
   }

   if (handled == FALSE) { 
      std::string ext_tmp = coot::util::file_name_extension(file_name);
      if (file_type_coords(file_name.c_str())) {
	 int imol = read_pdb(file_name);
	 if (is_valid_model_molecule(imol))
	    handled = TRUE;
	 else
	    std::cout << "INFO:: " << file_name << " was not a coordinates file" << std::endl;
      } else { 
	 if (ext == ".mtz") {
	    std::vector<int> imol_map = auto_read_make_and_draw_maps(file_name.c_str());
	    if (is_valid_map_molecule(imol_map.front()))
	       handled = TRUE;
	 }
      }
   }
   return handled;
}
