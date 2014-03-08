
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

gboolean
on_gl_canvas_drag_drop(GtkWidget *widget,
		       GdkDragContext *context,
		       gint x, gint y,
		       guint time,
		       gpointer user_data) {

   gboolean is_valid_drop_site = TRUE;
   // Request the data from the source.
   if (context->targets) {
      GdkAtom target_type =
	 GDK_POINTER_TO_ATOM(g_list_nth_data(context->targets, TARGET_STRING));
      
      gtk_drag_get_data(widget, context,  
			target_type,    /* the target type we want (a string) */
			time);
   } else {
      std::cout << "ERROR:: null dnd context" << std::endl;
   } 
   return  is_valid_drop_site;
}

void
on_drag_data_received (GtkWidget *widget, 
		       GdkDragContext *context, 
		       gint x, gint y,
		       GtkSelectionData *selection_data, 
		       guint target_type, 
		       guint time,
		       gpointer data) {

   gboolean dnd_success = FALSE;
   gboolean delete_selection_data = FALSE;
   
   // Deal with what the source sent over
   if((selection_data != NULL) && (selection_data-> length >= 0)) {
      if (target_type == TARGET_STRING) {
	 std::string uri_string = (gchar*)selection_data-> data;
	 dnd_success = handle_drag_and_drop_string(uri_string);
     delete_selection_data = TRUE;
      } 
   }
   gtk_drag_finish (context, dnd_success, delete_selection_data, time);
}


//! \brief handle the string that get when a file or URL is dropped.
//
// uri string can be a concatenation of string with \ns between them
// (and a \n to end)
// 
int handle_drag_and_drop_string(const std::string &uri_in) {

   int handled = FALSE;
   bool tried_already = false;
   std::string uri = uri_in;
   // std::cout << "handle_drag_and_drop_string() :" << uri<< ":" << std::endl;
   // was it a file://xx?
   std::string url;
   std::string::size_type pos = uri.find_first_of('\n');
   while (pos != std::string::npos) {
      url = uri.substr(0, pos-1); // front part
      uri = uri.substr(pos+1); // back part
      // std::cout << "Now URI is " << uri.size() << ":" << uri << ":" << std::endl;
      // std::cout << "Now URL is " << url.size() << ":" << url << ":" << std::endl;
      if (url.length() > 7) {
	 if (url.substr(0,7) == "file://") {
	    // std::cout << "---:" << url << ": was a file " << std::endl;
	    std::string file;
#ifdef WINDOWS_MINGW
	    file = url.substr(8);
#else
	    file = url.substr(7);
#endif
	    handled = handle_drag_and_drop_single_item(file);
	    tried_already = true;
	    std::cout << "setting tried_already here 1" << std::endl;
	 }
      }
      pos = uri.find_first_of('\n');
      
   }

   if (! tried_already) {
      // OK, was it an HTTP type string?
      std::string url = uri_in;
      if (url.length() > 8) {

	 if (url.substr(0,7) == "http://") {
	    tried_already = true;
	    int l = url.length();
	    if (url[l-1] == '\n') { 
	       // std::cout << "extra \\n" << std::endl;
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
      if (file_type_coords(file_name.c_str())) {
	 int imol = read_pdb(file_name.c_str());
	 if (is_valid_model_molecule(imol))
	    handled = TRUE;
	 else
	    std::cout << "INFO:: " << file_name << " was not a coordinates file" << std::endl;
      } else { 
	 std::string ext = coot::util::file_name_extension(file_name);
	 if (ext == ".mtz") {
	    int imol_map = auto_read_make_and_draw_maps(file_name.c_str());
	    if (is_valid_map_molecule(imol_map))
	       handled = TRUE;
	 }
      }
   }
   return handled;
}
