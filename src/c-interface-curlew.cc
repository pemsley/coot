
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON // for std::ptrdiff_t
#include <string>
#include <cstddef>

#ifdef HAVE_CXX11
   #if defined(__clang__)
      #define BUILD_CURLEW
   #else
      #if defined(__GNUC__)
         #if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) < 40805
            // no curlew
         #else
            #define BUILD_CURLEW
         #endif
      #endif
   #endif
#endif


#ifdef BUILD_CURLEW
#include "json.hpp" // clever stuff from Niels Lohmann
using json = nlohmann::json;
#endif


#include <gtk/gtk.h>

#include "support.h"
#include "interface.h"
#include "c-interface.h"
#include "cc-interface.hh" // for coot_get_url()
#include "coot-version.hh"

#include "graphics-info.h" // for extensions register // after json.hpp

#include "curlew.h"
#include "curlew.hh"


void remove_file_curlew_menu_item_maybe() {

#ifdef BUILD_CURLEW

   // OK, keep the menu item in.

#else
   GtkWidget *menubar = lookup_widget(graphics_info_t::statusbar, "menubar1");
   if (menubar) {
      // gtk_container_foreach(GTK_CONTAINER(menubar), my_delete_file_curlew_menu_item, menubar);

      GList *dlist_1 = gtk_container_children(GTK_CONTAINER(menubar));
      while (dlist_1) {
         GtkWidget *w = static_cast<GtkWidget *>(dlist_1->data);
         std::string l = gtk_menu_item_get_label(GTK_MENU_ITEM(w));
         // std::cout << "l: " << l << std::endl;
         if (l == "_File") {
            GtkWidget *w_submenu = gtk_menu_item_get_submenu(GTK_MENU_ITEM(w));
            GList *dlist_2 = gtk_container_children(GTK_CONTAINER(w_submenu));
            GtkWidget *curlew_menu = 0;
            while (dlist_2) {
               GtkWidget *w_inner = static_cast<GtkWidget *>(dlist_2->data);
               std::string l_inner = gtk_menu_item_get_label(GTK_MENU_ITEM(w_inner));
               // std::cout << "l_inner: " << l_inner << std::endl;
               if (l_inner == "Curlew")
                  curlew_menu = w_inner;
               dlist_2 = dlist_2->next;
            }
            if (curlew_menu) {
               gtk_container_remove(GTK_CONTAINER(w_submenu), curlew_menu);
            }
         }
         dlist_1 = dlist_1->next;
      }
   } else {
      std::cout << "WARNING:: remove_file_curlew_menu_item_maybe() ooops no menubar" << std::endl;
   }
#endif // BUILD_CURLEW   
}

// put this in a widget header (maybe its own?)
GtkWidget *make_and_add_curlew_extension_widget(GtkWidget *dialog,
						GtkWidget *item_hbox,
						int idx,
						const std::string &icon,
						const std::string &name,
						const std::string &description,
						const std::string &date,
						const std::string &version,
						const std::string &checksum,
						const std::string &file_name,
						const std::string &download_dir,
						const std::string &url_curlew_prefix);

void curlew() {

#ifdef BUILD_CURLEW

   GtkWidget *w = create_curlew_dialog();
   graphics_info_t g;

   GtkWidget *vbox = lookup_widget(w, "curlew_vbox_for_extensions");
   GtkWidget *install_button = lookup_widget(w, "curlew_install_button");
   if (vbox) {
      std::string download_dir = "coot-download"; // FIXME
      std::string dl_fn = download_dir + "/info.json";

      // not https, that transfers nothing
      // (probably a curl configuration thing)
      // 2019-07-31 https is the only way now
      //
      std::string url_prefix = "https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/";
      url_prefix += "extensions";

      std::string url_curlew_prefix = url_prefix + "/curlew";
      std::string json_url = url_curlew_prefix + "/info.json";

      int r = coot_get_url(json_url.c_str(), dl_fn.c_str());
      bool is_empty = true; // now check that it isn't
      struct stat buf;
      int istat = stat(dl_fn.c_str(), &buf);
      if (istat == 0) { // OK, it exists...
	 if (buf.st_size > 0) {
	    is_empty = false;
	 } else {
	    std::cout << "WARNING:: empty file " << dl_fn << std::endl;
	    std::cout << "          maybe your curl needs OpenSSL?" << std::endl;
	    std::string s = "WARNING:: empty file " + dl_fn;
	    add_status_bar_text(s.c_str());
	 }
      }

      if (is_empty) {
	 if (coot::file_exists(dl_fn)) {
	    std::fstream f(dl_fn);
	    if (! f) {
	       std::cout << "WARNING:: Missing/bad curlew info file " << dl_fn << std::endl;
            } else {

	       std::string s;
	       f.seekg(0, std::ios::end);
	       s.reserve(f.tellg());
	       f.seekg(0, std::ios::beg);

	       s.assign((std::istreambuf_iterator<char>(f)),
			std::istreambuf_iterator<char>());
	       unsigned int n_already_done = 0;

	       try {
		  json j = json::parse(s);
		  json ls = j["extensions"];
		  // std::cout << "found " << ls.size() << " extensions" << std::endl;
		  int n_extensions = ls.size();

		  for (std::size_t i=0; i<ls.size(); i++) {
		     json &item = ls[i];
		     std::string name;
		     std::string description;
		     std::string date;
		     std::string version;
		     std::string icon;
		     std::string file_name;
		     std::string checksum;
		     std::string expired_version; // which version of coot has this built in
		     // so that the extension is no longer needed
		     bool expired = false;
		     bool have_this_or_more_recent = false;

		     json::iterator it;
		     it = item.find(std::string("name"));
		     if (it != item.end()) { name = it.value(); }
		     it = item.find(std::string("description"));
		     if (it != item.end()) { description = it.value(); }
		     it = item.find(std::string("date"));
		     if (it != item.end()) { date = it.value(); }
		     it = item.find(std::string("icon"));
		     if (it != item.end()) { icon = it.value(); }
		     it = item.find(std::string("file-name"));
		     if (it != item.end()) { file_name = it.value(); }
		     it = item.find(std::string("version"));
		     if (it != item.end()) { version = it.value(); }
		     it = item.find(std::string("checksum"));
		     if (it != item.end()) { checksum = it.value(); }
		     it = item.find(std::string("expired_version"));
		     if (it != item.end()) { expired_version = it.value(); }

		     // set expired here
		     if (! expired_version.empty()) {
			std::string c = coot_version();
			if (c > expired_version) {
			   expired = true;
			}
		     }

		     // set "have more recent" (or same) here
		     std::string vv = g.get_version_for_extension(file_name);
		     if (! vv.empty())
			if (vv >= version)
			   have_this_or_more_recent = true;

		     if (have_this_or_more_recent)
			n_already_done++;

		     GtkWidget *hbox = make_and_add_curlew_extension_widget(w, vbox, i, icon,
									    name, description, date,
									    version, checksum, file_name,
									    download_dir, url_curlew_prefix);
		     if (expired || have_this_or_more_recent)
			gtk_widget_set_sensitive(hbox, FALSE);

		  }

		  if (install_button)
		     g_object_set_data(G_OBJECT(install_button), "n_extensions",
				       GINT_TO_POINTER(n_extensions));

	       }
	       catch(const nlohmann::detail::type_error &e) {
		  std::cout << "ERROR:: " << e.what() << std::endl;
	       }
	       catch(const nlohmann::detail::parse_error &e) {
		  std::cout << "ERROR:: " << e.what() << std::endl;
	       }

	       GtkWidget *done_label = lookup_widget(GTK_WIDGET(w), "curlew_already_installed_label");
	       if (done_label) {
		  if (n_already_done>0) {
		     std::string txt = coot::util::int_to_string(n_already_done);
		     txt += " extension";
		     if (n_already_done != 1) txt += "s";
		     txt += " already installed";
		     gtk_label_set_text(GTK_LABEL(done_label), txt.c_str());
		     gtk_widget_show(done_label);
		  } else {
		     gtk_widget_hide(done_label);
		  }
	       }
	    }
	 }
      } // we've done the "empty" message already
   }

   gtk_widget_show(w);

#else
   // well take out the menu item then!
   std::cout << "No curlew with old GCC" << std::endl;
#endif // BUILD_CURLEW
}


GtkWidget *make_and_add_curlew_extension_widget(GtkWidget *dialog,
						GtkWidget *vbox,
						int idx,
						const std::string &icon,
						const std::string &name,
						const std::string &description,
						const std::string &date,
						const std::string &version,
						const std::string &checksum,
						const std::string &file_name,
						const std::string &download_dir,
						const std::string &url_curlew_prefix) {

   GtkWidget *item_hbox = gtk_hbox_new(FALSE, 0);
   
   std::string item_hbox_name = "curlew_extension_hbox_";
   item_hbox_name += coot::util::int_to_string(idx);
   g_object_set_data_full(G_OBJECT(dialog),
			  item_hbox_name.c_str(),
			  item_hbox,
			  (GtkDestroyNotify) gtk_widget_unref);
   gtk_widget_ref(item_hbox);

   // --------------- Icon -----------------
   GtkWidget *icon_widget = 0;
   if (icon.size() > 0) {
      std::string icon_url = url_curlew_prefix + "/" + icon;
      std::string icon_fn  =
	 coot::util::append_dir_file(download_dir,
				     coot::util::file_name_non_directory(icon));
      // std::cout << "get " << icon_url << " to " << icon_fn << std::endl;
      coot_get_url(icon_url.c_str(), icon_fn.c_str());
      if (coot::file_exists(icon_fn)) {
	 GError *error = NULL;
	 GtkWidget *w = gtk_image_new_from_file(icon_fn.c_str());
	 if (w) {
	    icon_widget = w;
	 } else {
	    std::cout << "Null icon" << std::endl;
	 }
      } else {
	 icon_widget = gtk_label_new("  Icon");
	 gtk_misc_set_alignment (GTK_MISC(icon_widget), 0, 0.5);
      }
   } else {
      std::cout << "No icon in item " << std::endl;
      icon_widget = gtk_label_new("  ----");
   }
   gtk_widget_set_usize(icon_widget, 50, -1);

   // --------------- Description -----------------
   std::string rr = "<b>";
   rr += name;
   rr += "</b>\n";
   rr += description;
   GtkWidget *description_label = gtk_label_new(rr.c_str());
   gtk_label_set_use_markup(GTK_LABEL(description_label), TRUE);
   gtk_misc_set_alignment (GTK_MISC(description_label), 0, 0.5);
   gtk_widget_set_usize(description_label, 320, -1);
   // --------------- Version -----------------
   GtkWidget *version_label = gtk_label_new(version.c_str());
   // --------------- Date -----------------
   GtkWidget *date_label = gtk_label_new(date.c_str());
   // --------------- Select -----------------
   GtkWidget *selected_check_button = gtk_check_button_new();
   std::string cb_name = "curlew_selected_check_button_";
   cb_name += coot::util::int_to_string(idx);
   // --------------------------------------

   gtk_box_pack_start(GTK_BOX(item_hbox), icon_widget,           TRUE, TRUE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox), description_label,     TRUE, TRUE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox), version_label,         TRUE, TRUE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox),    date_label,         TRUE, TRUE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox), selected_check_button, TRUE, TRUE, 0);

   gtk_widget_show(icon_widget);
   gtk_widget_show(description_label);
   gtk_widget_show(version_label);
   gtk_widget_show(date_label);
   gtk_widget_show(selected_check_button);
   gtk_widget_show(item_hbox);

   gtk_box_pack_start(GTK_BOX(vbox), item_hbox, TRUE, TRUE, 0);

   g_object_set_data_full(G_OBJECT(dialog),
			  cb_name.c_str(),
			  selected_check_button,
			  (GtkDestroyNotify) gtk_widget_unref);

   char *file_name_copy = new char[file_name.size() +1];
   strcpy(file_name_copy, file_name.c_str());
   g_object_set_data(G_OBJECT(selected_check_button),
		     "file-name", (gpointer) file_name_copy);

   if (! checksum.empty()) {
      char *checksum_copy = new char[checksum.size() + 1];
      strcpy(checksum_copy, checksum.c_str());
      g_object_set_data(G_OBJECT(selected_check_button), "checksum",
			(gpointer) checksum_copy);
   }

   gtk_widget_ref(selected_check_button); // ref after set_data?

   return item_hbox;
}

/*! \brief register an extension */
void register_extension(const std::string &name, const std::string &version) {

   graphics_info_t g;
   g.register_extension(name, version);

}


std::string version_for_extension(const std::string &name) {

   graphics_info_t g;
   return g.get_version_for_extension(name);

}
