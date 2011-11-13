
#ifdef HAVE_GOOCANVAS

#include <iostream>
#include "lbg-drag-and-drop.hh"

#include "lbg.hh"

extern "C" G_MODULE_EXPORT gboolean
on_lbg_drag_drop (GtkWidget *widget,
		  GdkDragContext *context,
		  gint x, gint y,
		  guint time,
		  gpointer user_data) {

   
   // gboolean is_valid_drop_site = TRUE;
   gboolean retval = FALSE;
   // Request the data from the source.
   if (context->targets) {
      GdkAtom target_type =
	 GDK_POINTER_TO_ATOM(g_list_nth_data(context->targets, TARGET_STRING));
      
      gtk_drag_get_data(widget, context,  
			target_type,    /* the target type we want (a string) */
			time);
      retval = TRUE;
   } else {
      std::cout << "ERROR:: null dnd context" << std::endl;
   }
   return retval;
}


extern "C" G_MODULE_EXPORT void
on_lbg_drag_data_received (GtkWidget *widget, 
			   GdkDragContext *context, 
			   gint x, gint y,
			   GtkSelectionData *selection_data, 
			   guint target_type, 
			   guint time,
			   gpointer user_data) {

   gboolean dnd_success = FALSE;
   gboolean delete_selection_data = FALSE;
   
   GtkWidget *canvas = GTK_WIDGET(user_data);
   lbg_info_t *l = static_cast<lbg_info_t *> (gtk_object_get_user_data(GTK_OBJECT(widget)));
   std::cout << "Here in on_lbg_drag_data_received() " << l << std::endl;
   if (l) {
      // Deal with what the source sent over
      if((selection_data != NULL) && (selection_data-> length >= 0)) {
	 if (target_type == TARGET_STRING) {
	    std::string uri_string = (gchar*)selection_data-> data;
	    dnd_success = l->handle_lbg_drag_and_drop_string(uri_string);
	 } 
      }
   }
   gtk_drag_finish (context, dnd_success, delete_selection_data, time);
}



int
lbg_info_t::handle_lbg_drag_and_drop_string(const std::string &uri_in) {

   int handled = FALSE;

   std::string uri = uri_in;
   
   std::cout << "handle this string :" << uri << ": " << std::endl;
   std::string::size_type pos = uri.find_first_of('\n');
   if (pos != std::string::npos) {
      // there was a carriage return, strip down the string
      uri = uri.substr(0, pos-1); // front part
   }
   handled = handle_lbg_drag_and_drop_single_item(uri);
   return handled;
}

int
lbg_info_t::handle_lbg_drag_and_drop_single_item(const std::string &uri) {

   int handled = FALSE;
   if (uri.length() > 7) {
      if (uri.substr(0,7)== "file://") {
	 std::cout << "---:" << uri << ": was a file:// string " << std::endl;
	 std::string file_name;
#ifdef WINDOWS_MINGW
	 file_name = uri.substr(8);
#else
	 file_name = uri.substr(7);
#endif

	 std::string ext = coot::util::file_name_extension(file_name);
	 if (ext == ".mdl" || ext == ".mol" || ext == ".mol2") { 
	    import_mol_from_file(file_name);
	 }
	 handled = TRUE;
      }
      if (uri.substr(0,7) == "http://") {
	 if (get_url_func_ptr_flag) { // we have a function get get urls?
	    int l = uri.length();
	    std::string uri_clean = uri;
	    if (uri[l-1] == '\n') 
	       uri_clean = uri.substr(0, l-1);
	    int status = coot::util::create_directory("coot-download"); // like make_directory_maybe()
	    if (status == 0) { // OK, we made it (or had it)
	       std::string::size_type pos = uri_clean.find_last_of('/');
	       if (pos != std::string::npos) {
		  // normal path
		  std::string url_file_name_file = uri_clean.substr(pos+1);
		  std::string ext = coot::util::file_name_extension(uri_clean);
		  if (ext == ".mol" || ext == ".sdf") {
		     std::string file_name =
			coot::util::append_dir_file("coot-download", url_file_name_file);
		     get_url_func_ptr(uri_clean.c_str(), file_name.c_str());
		     import_mol_from_file(file_name);
		  }

		  std::cout << "url_file_name_file :" << url_file_name_file << ":" << std::endl;
		  if (url_file_name_file == "image.png") {
		     if (uri_clean.find_last_of("/molecules/DB") != std::string::npos) {
			std::pair<std::string, std::string> s =
			   coot::util::split_string_on_last_slash(uri_clean);
			std::pair<std::string, std::string> ss =
			   coot::util::split_string_on_last_slash(s.first);
			if (ss.second.find("DB") != std::string::npos) {
			   std::string local_file = ss.second;
			   local_file += ".mol";
			   std::string drugbank_mol_url = "http://drugbank.ca/drugs/";
			   drugbank_mol_url += local_file;
			   std::string file_name =
			      coot::util::append_dir_file("coot-download", local_file);
			   get_url_func_ptr(drugbank_mol_url.c_str(), file_name.c_str());
			   import_mol_from_file(file_name);
			} 
		     }
		  }
	       }
	    }
	 }
      }
   }
   return handled;
}

#endif // HAVE_GOOCANVAS
