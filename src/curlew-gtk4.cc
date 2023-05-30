
#include <iostream>
#include <fstream>
#include "curlew-gtk4.hh"
#include "widget-from-builder.hh"
#include "utils/coot-utils.hh"
#include "graphics-info.h"

#include "json.hpp" // clever stuff from Niels Lohmann
using json = nlohmann::json;

// in c-interface.hh
int coot_get_url(const std::string &url, const std::string &file_name);
// in c-interface.h
void run_script(const char *filename);

// return a status and a failure message.
// return bool true on success.
std::pair<bool, std::string>
curlew_install_extension_file_gtk4(const std::string &script_here_file_name) {

   bool success = false;
   std::string failure_message;
   if (coot::file_exists_and_non_empty(script_here_file_name)) {
      std::string home_directory = coot::get_home_dir();
      if (!home_directory.empty()) {
         std::string file_name = coot::util::file_name_non_directory(script_here_file_name);
         std::string preferences_dir = coot::util::append_dir_dir(home_directory, ".coot");
         std::string preferences_file_name = coot::util::append_dir_file(preferences_dir, file_name);
         // std::cout << "debug:: attempting to copy \"" << script_here_file_name << "\" as \"" << preferences_file_name
         // << "\"" << std::endl;
         int status = coot::copy_file(script_here_file_name, preferences_file_name); // it returns a bool actually
         if (status == false) {
            // std::cout << "WARNING:: Copy file script failed: " << script_here_file_name << std::endl;
            FILE *fp = fopen(script_here_file_name.c_str(), "r");
            PyRun_SimpleFile(fp, script_here_file_name.c_str());
            fclose(fp);
            failure_message = "WARNING:: Copy file script failed: " + script_here_file_name;
         } else {
            // cool.
            // std::cout << "debug:: run_script() called on " << preferences_file_name << std::endl;
            FILE *fp = fopen(script_here_file_name.c_str(), "r");
            PyRun_SimpleFile(fp, preferences_file_name.c_str());
            fclose(fp);
            success = true;
         }
      }
   }
   return std::make_pair(success, failure_message);
}

// put this in coot utils
int coot_rename(const std::string &f1, const std::string &f2) {

   // return 0 on success

   // std::cout << "coot_rename \"" << f1 << "\" to \"" << f2 << "\"" << std::endl;

#ifndef WINDOWS_MINGW
   int status = rename(f1.c_str(), f2.c_str());
#else
   int status = coot::rename_win(f1.c_str(), f2.c_str());
#endif
   return status;
}

int
curlew_uninstall_extension_file_gtk4(const std::string &script_file_name) {

   std::string home_directory = coot::get_home_dir();
   std::string preferences_dir = coot::util::append_dir_dir(home_directory, ".coot");
   std::string preferences_file_name = coot::util::append_dir_file(preferences_dir, script_file_name);
   std::string renamed_file = preferences_file_name + "_uninstalled";
   int status = coot_rename(preferences_file_name, renamed_file);
   return status;
}

// in c-interface-gui.cc
std::pair<bool, std::string>checksums_match(const std::string &file_name, const std::string &checksum);

GtkWidget *
curlew_dialog() {

   auto get_curlew_url_prefix = [] () {
      // put this in utils somewhere
      std::string url_prefix;
#ifndef WINDOWS_MINGW
      url_prefix += "https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/";
#else
      url_prefix += "https://bernhardcl.github.io/coot/";
#endif

      std::string coot_version_dir_prefix = "curlew-extensions/gtk4/Coot-1";
      std::string scripts_dir_prefix = "scripts";
      std::string url = coot::util::append_dir_dir(url_prefix, coot_version_dir_prefix);
      return url;
   };

   auto clear_the_grid = [] (GtkWidget *grid) {
      GtkWidget *child = gtk_widget_get_first_child(GTK_WIDGET(grid));
      while (child) {
         GtkWidget *next = gtk_widget_get_next_sibling(child);
         gtk_grid_remove(GTK_GRID(grid), child);
         child = next; // for next round
      }
   };

   auto file_to_string = [] (const std::string &dl_fn) {
      std::string s;
      std::fstream f(dl_fn);
      f.seekg(0, std::ios::end);
      s.reserve(f.tellg());
      f.seekg(0, std::ios::beg);
      s.assign((std::istreambuf_iterator<char>(f)),
               std::istreambuf_iterator<char>());
      return s;
   };

   class extension_info_t {
   public:
      std::string name;
      std::string file_name;
      std::string version;
      std::string description;
      std::string icon;
      std::string date;
      std::string checksum;
      std::string installed_extension_version; // can be blank of course
      extension_info_t(std::string &name, std::string &file_name, std::string &version, std::string &description,
                       std::string &icon, std::string &date, std::string &checksum, const std::string &vv) :
         name(name), file_name(file_name), version(version), description(description), icon(icon),
         date(date), checksum(checksum), installed_extension_version(vv) {}
      std::string make_label() const {
         std::string s = "<b>" + name + "</b>\n";
         s += description;
         return s;
      }
      bool already_installed_p() const {
         bool status = false;
         if (! installed_extension_version.empty())
            if (installed_extension_version >= version)
               status = true;
         return status;
      }
   };

   std::string url_prefix = get_curlew_url_prefix();
   std::cout << "url_prefix: " << url_prefix << std::endl;
   std::string scripts_dir = coot::util::append_dir_file(url_prefix, "scripts");

   std::string download_dir = "coot-download";
   download_dir = coot::get_directory(download_dir.c_str());

   auto add_extension_to_grid = [download_dir, url_prefix, scripts_dir] (const extension_info_t &extension, GtkWidget *grid, int row) {

      bool already_installed = extension.already_installed_p();
      GtkWidget *icon_widget = 0;
      if (! extension.icon.empty()) {
         std::string icon_file_name_nd = coot::util::file_name_non_directory(extension.icon);
         std::string icon_file_name = coot::util::append_dir_file(download_dir, icon_file_name_nd);
         if (coot::file_exists_and_non_empty(icon_file_name)) {
         } else {
            std::string icon_dir = coot::util::append_dir_dir(url_prefix, "images");
            std::string icon_url = coot::util::append_dir_file(icon_dir, icon_file_name_nd);
            std::cout << "icon_url: " << icon_url << std::endl;
            coot_get_url(icon_url, icon_file_name);
         }
         if (coot::file_exists_and_non_empty(icon_file_name)) {
            icon_widget = gtk_image_new_from_file(icon_file_name.c_str());
            gtk_widget_set_size_request(icon_widget, 60, -1);
         } else {
            icon_widget = gtk_label_new(" "); // a placeholder
         }
      } else {
         icon_widget = gtk_label_new(" "); // a placeholder
      }

      GtkWidget *w_0 = icon_widget;
      GtkWidget *w_1 = gtk_label_new(extension.make_label().c_str());
      GtkWidget *w_2 = gtk_label_new(extension.version.c_str());
      GtkWidget *w_3 = gtk_label_new(extension.date.c_str());
      GtkWidget *w_i = gtk_button_new_with_label("Install");
      GtkWidget *w_u = gtk_button_new_with_label("Uninstall");
      gtk_label_set_use_markup(GTK_LABEL(w_1), TRUE);
      gtk_grid_attach(GTK_GRID(grid), w_0, 0, row, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), w_1, 1, row, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), w_2, 2, row, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), w_3, 3, row, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), w_i, 4, row, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), w_u, 5, row, 1, 1);
      gtk_widget_set_margin_start(w_0, 6);
      gtk_widget_set_margin_end(w_0, 6);
      gtk_widget_set_margin_top(w_0, 2);
      gtk_widget_set_margin_bottom(w_0, 2);
      gtk_widget_set_margin_start(w_3, 6);
      gtk_widget_set_margin_end(w_3, 6);
      gtk_widget_set_margin_top(w_i, 4);
      gtk_widget_set_margin_bottom(w_i, 4);
      gtk_widget_set_margin_start(w_i, 4);
      gtk_widget_set_margin_end(w_i, 4);
      gtk_widget_set_margin_top(w_u, 4);
      gtk_widget_set_margin_bottom(w_u, 4);
      gtk_widget_set_margin_start(w_u, 4);
      gtk_widget_set_margin_end(w_u, 4);

      if (already_installed)
         gtk_widget_set_visible(w_i, FALSE);
      else
         gtk_widget_set_visible(w_u, FALSE); // we can't uninstall it

      auto uninstall_callback_func = +[] (GtkWidget *self, gpointer user_data) {
         const std::string *file_name_p = static_cast<std::string *>(g_object_get_data(G_OBJECT(self), "file_name"));
         std::string file_name(*file_name_p);
         std::cout << "uninstall file_name " << file_name << std::endl;
         int status = curlew_uninstall_extension_file_gtk4(file_name);
         if (status == 0) { // OS success return value
            GtkWidget *install_button = GTK_WIDGET(g_object_get_data(G_OBJECT(self), "install-button"));
            gtk_widget_set_visible(self, FALSE);
            gtk_widget_set_visible(install_button, TRUE);
         } else {
            std::string m = "WARNING:: failed to uninstall " + file_name;
            graphics_info_t::info_dialog(m);
         }
      };

      auto install_callback_func = +[] (GtkWidget *self, gpointer user_data) {
         const std::string *file_name_p    = static_cast<std::string *>(g_object_get_data(G_OBJECT(self), "file_name"));
         const std::string *download_dir_p = static_cast<std::string *>(g_object_get_data(G_OBJECT(self), "download_dir"));
         const std::string *checksum_p     = static_cast<std::string *>(g_object_get_data(G_OBJECT(self), "checksum"));
         std::string file_name(*file_name_p);
         std::string download_dir(*download_dir_p);
         std::string checksum(*checksum_p);
         GtkWidget *uninstall_button = GTK_WIDGET(g_object_get_data(G_OBJECT(self), "uninstall-button"));
         // std::cout << "get this url_file_name: " << file_name << " to dir " << download_dir << std::endl;
         std::string file_name_here = coot::util::append_dir_file(download_dir, file_name);
         std::string *scripts_dir_p = static_cast<std::string *>(user_data);
         std::string scripts_dir = *scripts_dir_p;
         std::string url = coot::util::append_dir_file(scripts_dir, file_name);
         coot_get_url(url, file_name_here);
         if (coot::file_exists_and_non_empty(file_name_here)) {
            std::pair<bool, std::string> checksum_result = checksums_match(file_name_here, checksum);
            if (checksum_result.first) { 
               std::pair<bool, std::string> result = curlew_install_extension_file_gtk4(file_name_here);
               if (result.first) {
                  // hide the "Install" button
                  gtk_widget_set_visible(self, FALSE);
                  // show the "Uninstall" button
                  gtk_widget_set_visible(uninstall_button, TRUE);
               }
            } else {
               std::string m = std::string("WARNING:: checksums do not match ") + file_name_here
                  + std::string("\n") + checksum_result.second;
               graphics_info_t::info_dialog(m);
            }
         }
      };

      std::string *sp = new std::string(extension.file_name);
      std::string *cp = new std::string(extension.checksum);
      std::string *dp = new std::string(download_dir);
      std::string *scripts_dir_p = new std::string(scripts_dir);
      g_object_set_data(G_OBJECT(w_i), "file_name",    sp);
      g_object_set_data(G_OBJECT(w_u), "file_name",    sp); // same as above.. is this safe?
      g_object_set_data(G_OBJECT(w_i), "download_dir", dp);
      g_object_set_data(G_OBJECT(w_i), "checksum",     cp);
      g_object_set_data(G_OBJECT(w_i), "uninstall-button", w_u);
      g_object_set_data(G_OBJECT(w_u),   "install-button", w_i);
      g_signal_connect(G_OBJECT(w_i), "clicked", G_CALLBACK(  install_callback_func), scripts_dir_p);
      g_signal_connect(G_OBJECT(w_u), "clicked", G_CALLBACK(uninstall_callback_func), scripts_dir_p);
   };

   auto add_extensions_to_grid = [add_extension_to_grid] (const std::vector<extension_info_t> &extensions, GtkWidget *grid) {

      int row = 1; // starting row
      for (const auto &e : extensions) {
         add_extension_to_grid(e, grid, row);
         row++; // for next
      }
   };

   auto add_header_to_grid = [] (GtkWidget *grid) {
      GtkWidget *label_0 = gtk_label_new("   "); // icon
      GtkWidget *label_1 = gtk_label_new("Extension");
      GtkWidget *label_2 = gtk_label_new(" Version ");
      GtkWidget *label_3 = gtk_label_new("Date");
      GtkWidget *label_4 = gtk_label_new("                            ");
      GtkWidget *label_5 = gtk_label_new("                            ");
      gtk_grid_attach(GTK_GRID(grid), label_0, 0, 0, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), label_1, 1, 0, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), label_2, 2, 0, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), label_3, 3, 0, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), label_4, 4, 0, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), label_5, 5, 0, 1, 1);
   };

   GtkWidget *dialog = widget_from_builder("curlew_dialog");
   GtkWidget *grid   = widget_from_builder("curlew_grid");
   clear_the_grid(grid);
   add_header_to_grid(grid);
   gtk_widget_set_vexpand(grid, TRUE);

   std::string json_url   = coot::util::append_dir_file(url_prefix, "info/curlew-info.json");

   // std::cout << "debug:: here with json_url: " << json_url << std::endl;
   // std::cout << "debug:: here with scripts_dir: " << scripts_dir << std::endl;
   std::string dl_fn = coot::util::append_dir_file(download_dir, "curlew-info.json");
   int r = coot_get_url(json_url, dl_fn);

   if (coot::file_exists_and_non_empty(dl_fn)) {
      std::string s = file_to_string(dl_fn);
      std::vector<extension_info_t> extensions;
      try {
         graphics_info_t g;
         json j = json::parse(s);
         json ls = j["extensions"];
         for (std::size_t i=0; i<ls.size(); i++) {
            json &item = ls[i];
            std::string name;
            std::string file_name;
            std::string version;
            std::string description;
            std::string icon;
            std::string date;
            std::string checksum;
            json::iterator it;
            it = item.find(std::string("name"));
            if (it != item.end()) { name = it.value(); }
            it = item.find(std::string("file-name"));
            if (it != item.end()) { file_name = it.value(); }
            it = item.find(std::string("version"));
            if (it != item.end()) { version = it.value(); }
            it = item.find(std::string("description"));
            if (it != item.end()) { description = it.value(); }
            it = item.find(std::string("icon"));
            if (it != item.end()) { icon = it.value(); }
            it = item.find(std::string("date"));
            if (it != item.end()) { date = it.value(); }
            it = item.find(std::string("checksum"));
            if (it != item.end()) { checksum = it.value(); }
            // set "have more recent" (or same) here
            std::string vv = g.get_version_for_extension(file_name);

            if (false)
               std::cout << "debug:: " << " name: " << name << " file_name: " << file_name << " version: " << version
                         << " description: " << description << " icon: " << icon << " date: " << date
                         << " checksum: " << checksum << " vv: " << vv << std::endl;
            extension_info_t ei(name, file_name, version, description, icon, date, checksum, vv);
            extensions.push_back(ei);
         }
         add_extensions_to_grid(extensions, grid);
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
      catch(const nlohmann::detail::type_error &e) {
         std::cout << "ERROR:: " << e.what() << std::endl;
      }
      catch(const nlohmann::detail::parse_error &e) {
         std::cout << "ERROR:: " << e.what() << std::endl;
      }
   } else {
      std::cout << "WARNING:: file " << dl_fn << " does not exist or was empty" << std::endl;
   }

   return dialog;

}
