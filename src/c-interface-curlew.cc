
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

#include "c-interface-curlew.hh"

void remove_file_curlew_menu_item_maybe() {

#ifdef BUILD_CURLEW

   // OK, keep the menu item in.

#else
   std::cout << "No curlew - not compiled" << std::endl;
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
						const std::string &url_curlew_prefix,
                                                bool make_and_add_curlew_extension_widget);

void curlew() {

#ifdef BUILD_CURLEW

   GtkWidget *w = create_curlew_dialog();
   graphics_info_t g;

   GtkWidget *vbox = lookup_widget(w, "curlew_vbox_for_extensions");
   GtkWidget *install_selected_button = lookup_widget(w, "curlew_install_button");
   // new method does individual installs, not all at once
   gtk_widget_hide(install_selected_button);

   // install_button callback:  curlew_dialog_install_extensions()
   if (vbox) {
      std::string download_dir = "coot-download";
      download_dir = coot::get_directory(download_dir.c_str());
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

      // hack in a pre-downloaded file
      // dl_fn = "hack-info.json";

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

      if (! is_empty) {
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

                  unsigned int n_already_installed = 0;
                  unsigned int n_available = 0;
                  for (std::size_t i=0; i<ls.size(); i++) {
                     json &item = ls[i];
                     std::string file_name;
                     std::string version;
                     json::iterator it;
                     it = item.find(std::string("file-name"));
                     if (it != item.end()) { file_name = it.value(); }
                     it = item.find(std::string("version"));
                     if (it != item.end()) { version = it.value(); }
                     // set "have more recent" (or same) here
                     std::string vv = g.get_version_for_extension(file_name);
                     if (! vv.empty()) {
                        if (vv >= version)
                           n_already_installed++;
                        else
                           n_available++;
                     } else {
                        n_available++;
                     }
                  }

                  std::cout << "DEBUG:: n_already_installed: " << n_already_installed
                           << " n_available " << n_available << std::endl;

                  for(unsigned int iround=0; iround<2; iround++) {
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

                        if (iround == 0) // add these up only once
                           if (have_this_or_more_recent)
                              n_already_done++;

                        bool do_it = false;
                        if (iround == 0)
                           if (have_this_or_more_recent)
                              do_it = true;
                        if (iround == 1)
                           if (! have_this_or_more_recent)
                              do_it = true;

                        if (do_it)
                              GtkWidget *hbox = make_and_add_curlew_extension_widget(w, vbox, i, icon,
                                    name, description, date,
                                    version, checksum, file_name,
                                    download_dir, url_curlew_prefix,
                                    have_this_or_more_recent);

                     }
                  } // rounds
               } // try

               catch(const nlohmann::detail::type_error &e) {
                  std::cout << "ERROR:: " << e.what() << std::endl;
               }
               catch(const nlohmann::detail::parse_error &e) {
                  std::cout << "ERROR:: " << e.what() << std::endl;
               }

               GtkWidget *done_label = lookup_widget(GTK_WIDGET(w), "curlew_already_installed_label");
               if (done_label) {
                  if (n_already_done > 0) {
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

void curlew_install_extension(GtkWidget *w, gpointer data) {

   gchar *file_name_cstr = static_cast<gchar *> (g_object_get_data(G_OBJECT(w), "file-name"));
   gchar *checksum_cstr  = static_cast<gchar *> (g_object_get_data(G_OBJECT(w), "checksum"));
   if (file_name_cstr) {
      if (checksum_cstr) {
         std::string fn(file_name_cstr);
         std::string checksum(checksum_cstr);
         curlew_install_extension_file(fn, checksum);
      } else {
         std::cout << "Null thing in curlew_install_extension" << std::endl;
      }
   } else {
      std::cout << "Null thing in curlew_install_extension" << std::endl;
   }
}

// return uninstall status (true for done)
void curlew_uninstall_extension(GtkWidget *w, gpointer data) {

   gchar *file_name_cstr = static_cast<gchar *> (g_object_get_data(G_OBJECT(w), "file-name"));
   if (file_name_cstr) {
      std::string fn(file_name_cstr);
      bool status = curlew_uninstall_extension_file(fn);
      if (status) {
         // Hide the uninstall button, because once it's been uninstalled we can't install it
         // any more
         gtk_widget_hide(w);
         if (data) {
            GtkWidget *install_button = static_cast<GtkWidget *>(data);
            gtk_widget_show(install_button);
         }
      }
   } else {
      std::cout << "Null thing in curlew_uninstall_extension" << std::endl;
   }
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
						const std::string &url_curlew_prefix,
                                                bool have_this_or_more_recent) {

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
      std::string icon_fn  = coot::util::append_dir_file(download_dir,
                             coot::util::file_name_non_directory(icon));
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
   gtk_widget_set_size_request(description_label, 340, -1);
   // --------------- Version -----------------
   GtkWidget *version_label = gtk_label_new(version.c_str());
   gtk_widget_set_size_request(version_label, 40, -1);
   // --------------- Date -----------------
   GtkWidget *date_label = gtk_label_new(date.c_str());
   // --------------- Uninstall -----------------
   GtkWidget *uninstall_frame = gtk_frame_new(NULL);
   GtkWidget *uninstall_button = gtk_button_new();
   std::string ucb_name = "curlew_uninstall_button_";
   ucb_name += coot::util::int_to_string(idx);
   gtk_button_set_label(GTK_BUTTON(uninstall_button), "Uninstall");
   gtk_widget_set_size_request(uninstall_frame, 100, -1);
   gtk_frame_set_shadow_type(GTK_FRAME(uninstall_frame), GTK_SHADOW_NONE);
   // --------------- Install -----------------
   GtkWidget *install_frame = gtk_frame_new(NULL);
   GtkWidget *install_button = gtk_button_new();
   std::string icb_name = "curlew_install_button_";
   icb_name += coot::util::int_to_string(idx);
   gtk_button_set_label(GTK_BUTTON(install_button), "Install");
   gtk_widget_set_size_request(install_frame, 100, -1);
   gtk_frame_set_shadow_type(GTK_FRAME(install_frame), GTK_SHADOW_NONE);
   // --------------------------------------

   // question to self: does strcpy() set the ending \0?
   char *file_name_copy_1 = new char[file_name.size() +1];
   char *file_name_copy_2 = new char[file_name.size() +1];
   strcpy(file_name_copy_1, file_name.c_str());
   strcpy(file_name_copy_2, file_name.c_str());
   g_object_set_data(G_OBJECT(  install_button), "file-name", (gpointer) file_name_copy_1);
   g_object_set_data(G_OBJECT(uninstall_button), "file-name", (gpointer) file_name_copy_2);
   char *checksum_copy = new char[checksum.size() +1];
   strcpy(checksum_copy, checksum.c_str());
   g_object_set_data(G_OBJECT(  install_button), "checksum",  (gpointer) checksum_copy);

   GdkColor color_green;
   GdkColor color_blue;
   gdk_color_parse ("#aabbaa", &color_green);
   gdk_color_parse ("#aaaabb", &color_blue);
   gtk_widget_modify_bg(GTK_WIDGET(  install_button), GTK_STATE_NORMAL, &color_green);
   gtk_widget_modify_bg(GTK_WIDGET(uninstall_button), GTK_STATE_NORMAL, &color_blue);

   gtk_container_add(GTK_CONTAINER(  install_frame),   install_button);
   gtk_container_add(GTK_CONTAINER(uninstall_frame), uninstall_button);
   gtk_box_pack_start(GTK_BOX(item_hbox), icon_widget,       TRUE, FALSE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox), description_label, TRUE, FALSE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox), version_label,     FALSE, FALSE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox), date_label,        TRUE, FALSE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox), install_frame,     TRUE, FALSE, 0);
   gtk_box_pack_start(GTK_BOX(item_hbox), uninstall_frame,   TRUE, FALSE, 0);

   gtk_widget_show(icon_widget);
   gtk_widget_show(description_label);
   gtk_widget_show(version_label);
   gtk_widget_show(date_label);
   gtk_widget_show(install_frame);
   gtk_widget_show(uninstall_frame);
   gtk_widget_show(item_hbox);

   if (have_this_or_more_recent)
      gtk_widget_show(uninstall_button);
   else
      gtk_widget_show(install_button);

   gtk_box_pack_start(GTK_BOX(vbox), item_hbox, TRUE, TRUE, 0);

   g_signal_connect(  install_button, "clicked", G_CALLBACK(curlew_install_extension),   NULL);
   g_signal_connect(uninstall_button, "clicked", G_CALLBACK(curlew_uninstall_extension), install_button);

   g_object_set_data_full(G_OBJECT(dialog),
			  icb_name.c_str(),
			  install_button,
			  (GtkDestroyNotify) gtk_widget_unref);

   g_object_set_data_full(G_OBJECT(dialog),
			  ucb_name.c_str(),
			  uninstall_button,
			  (GtkDestroyNotify) gtk_widget_unref);

   if (! checksum.empty()) {
      char *checksum_copy = new char[checksum.size() + 1];
      strcpy(checksum_copy, checksum.c_str());
      g_object_set_data(G_OBJECT(install_button), "checksum",
                        (gpointer) checksum_copy);
   }

   gtk_widget_ref(install_button); // ref after set_data?
   gtk_widget_ref(uninstall_button); // ref after set_data?

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
