/* src/drag-and-drop.cc
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2013 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#if 0 // 20220602-PE FIXME
gboolean
on_gl_canvas_drag_drop(GtkWidget *widget,
		       GdkDragContext *context,
		       gint x, gint y,
		       guint time,
		       gpointer user_data) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   // 20220528-PE FIXME mouse
   return 0;
#else
   gboolean is_valid_drop_site = TRUE;
   // Request the data from the source.
   GList *targets = gdk_drag_context_list_targets(context);
   if (targets) {
      GdkAtom target_type = GDK_POINTER_TO_ATOM(g_list_nth_data(targets, TARGET_STRING));
      
      gtk_drag_get_data(widget, context,  
			target_type,    /* the target type we want (a string) */
			time);
   } else {
      std::cout << "ERROR:: null dnd context" << std::endl;
   } 
   return  is_valid_drop_site;
#endif
}
#endif

#if 0 // 20220602-PE FIXME
void
on_drag_data_received (GtkWidget *widget, 
		       GdkDragContext *context, 
		       gint x, gint y,
		       GtkSelectionData *selection_data,
		       guint target_type, 
		       guint time,
		       gpointer data) {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
   // 20220528-PE FIXME mouse
#else
   gboolean dnd_success = FALSE;
   gboolean delete_selection_data = FALSE;
   
   // Deal with what the source sent over
   gint len = gtk_selection_data_get_length(selection_data);
   if ((selection_data != NULL) && (len >= 0)) {
      std::string uri_string;
      if (target_type == TEXT_URL) {
         // we have an url to deal with
         uri_string = (gchar *) gtk_selection_data_get_text(selection_data);
         dnd_success = handle_drag_and_drop_string(uri_string);
      }
      else if (target_type == TEXT_URI) {
         // we have a text uri (file!?)
         gchar **uris;
         gint i = 0;
         gchar *res = 0;
         uris = g_uri_list_extract_uris((gchar*) gtk_selection_data_get_text(selection_data));
         if (uris) {
            while (uris[i] != 0) {
               res = g_filename_from_uri(uris[i], NULL, NULL);
               i++;
               if (res != NULL) {
                  // everything fine
                  dnd_success = handle_drag_and_drop_single_item((gchar *)res);
               } else {
                  // not a file (shouldnt necessary happen - urls are dealt above and 
                  // simple strings below
                  uri_string = (gchar *) gtk_selection_data_get_data(selection_data);
                  dnd_success = handle_drag_and_drop_string(uri_string);
               }
            }
         }
         g_free(res);
         g_strfreev(uris);
      }
      else if (target_type == TARGET_STRING) {
         // simple string could call an extra function here too
         uri_string = (gchar *) gtk_selection_data_get_text(selection_data);
         dnd_success = handle_drag_and_drop_string(uri_string);
      }
      delete_selection_data = TRUE;
      
   }
   gtk_drag_finish (context, dnd_success, delete_selection_data, time);
#endif
}
#endif


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

   // std::cout << ":::::::::::::::: handle_drag_and_drop_string(" << uri_in << ")" << std::endl;

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
	    
	    
	    int status = make_directory_maybe("coot-download");
	    if (status == 0) { // OK, we made it (or had it)
	       std::string url_file_name_file = url;

	       std::string ext = coot::util::file_name_extension(url);
               
	       if (ext == ".png") {
		  // special rule - convert the url of the png to that
		  // of an accession code.
		  if (url.find("/PDBimages/")) {
		     std::pair<std::string, std::string> s =
			coot::util::split_string_on_last_slash(url);
		     std::pair<std::string, std::string> ss =
			coot::util::split_string_on_last_slash(s.first);
		     std::pair<std::string, std::string> sss =
			coot::util::split_string_on_last_slash(ss.first);
		     tried_already = true;
		     handled = FALSE;
		     if (ss.second.length() == 2) {
			if (sss.second.length() == 2) {
			   std::string code;
			   code += ss.second[0];
			   code += sss.second;
			   code += ss.second[1];
			   get_coords_for_accession_code(code.c_str());
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
		  std::string file_name =
		     coot::util::append_dir_file("coot-download", url_file_name_file);
		  coot_get_url(url.c_str(), file_name.c_str());
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
      std::cout << "here at the end of handle_drag_and_drop_string() " << std::endl;
      if (coot::file_exists(url)) {
	 handled = handle_drag_and_drop_single_item(url);
      } 
   } 
   return handled;
}

int handle_drag_and_drop_single_item(const std::string &file_name) {

   int handled = FALSE;
   // std::cout << "handle_drag_and_drop_single_item() " << file_name << ":" << std::endl; 

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
	 int imol = read_pdb(file_name.c_str());
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
