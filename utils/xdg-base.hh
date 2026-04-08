/* utils/xdg-base.hh
 *
 * Copyright 2024 by Paul Emsley
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


#ifndef UTILS_XDG_BASE_HH
#define UTILS_XDG_BASE_HH

// 20240902-PE Tweak from Charles
#ifdef _MSC_VER
#define WINDOWS
#endif

#ifdef WINDOWS
#else
#include <pwd.h>
#include <unistd.h>
#include <sys/types.h>
#endif

#include <vector>
#include <iostream> // remove this when done debugging
#include <cstdlib>
#include <string>
#include <filesystem>

class xdg_t {
   std::filesystem::path   data_home;
   std::filesystem::path  state_home;
   std::filesystem::path  cache_home;
   std::filesystem::path config_home;
   std::filesystem::path runtime_dir;
   std::string data_dirs;
   std::string config_dirs;
   std::string package_name; // Coot

   std::filesystem::path get_home_dir() {
      std::string home;
#ifdef WINDOWS
#else
      struct passwd *pw = getpwuid(getuid());
      const char *home_str = pw->pw_dir;
      if (home_str)
         home = std::string(home_str);
#endif
      return std::filesystem::path(home);
   }

   std::vector<std::filesystem::path> get_scripts_internal(std::filesystem::path path,
                                                           const std::string &extension) const {
      std::vector<std::filesystem::path> scripts;
      for (const auto &entry : std::filesystem::directory_iterator(path)) {
         if (entry.path().extension() == extension)
            scripts.push_back(entry);
      }
      return scripts;
   }

public:
   xdg_t() : package_name("Coot") {
      init();
   }
   explicit xdg_t(const std::string &pn) : package_name(pn) {
      init();
   }

   void init() {

      // config_home is for hand-crafted/user-specified configurations
      // cache_home is non-essential data (backups)
      // state_home is for scripts (auto) written on Exit
      // quick-save-as state script goes into state_home
      // curlew scripts go into config_home (they are not state files)
      // history scripts goes into STATE_HOME
      char *e;
      e = std::getenv("XDG_DATA_HOME");   if (e)   data_home = e;
      e = std::getenv("XDG_STATE_HOME");  if (e)  state_home = e;
      e = std::getenv("XDG_CACHE_HOME");  if (e)  cache_home = e;
      e = std::getenv("XDG_CONFIG_HOME"); if (e) config_home = e;
      e = std::getenv("XDG_RUNTIME_DIR"); if (e) runtime_dir = e;
      e = std::getenv("XDG_DATA_DIRS");   if (e)   data_dirs = e;
      e = std::getenv("XDG_CONFIG_DIRS"); if (e) config_dirs = e;
      if (data_home.empty()) {
         std::filesystem::path d = get_home_dir();
         d.append(".local");
         d.append("share");
         d.append(package_name);
         data_home = d;
      }
      if (config_home.empty()) {
         std::filesystem::path d = get_home_dir();
         d.append(".config");
         d.append(package_name);
         config_home = d;
      }
      if (state_home.empty()) {
         std::filesystem::path d = get_home_dir();
         d.append(".local");
         d.append("state");
         d.append(package_name);
         state_home = d;
      }
      if (cache_home.empty()) {
         std::filesystem::path d = get_home_dir();
         d.append(".cache");
         d.append(package_name);
         cache_home = d;
      }
      if (data_dirs.empty()) {
         data_dirs = "/usr/local/share:/usr/share";
      }
      if (config_dirs.empty()) {
         config_dirs = "/etc/xdg";
      }
   }
   std::filesystem::path get_state_home() const {
      try {
         if (!std::filesystem::is_directory(state_home))
            std::filesystem::create_directories(state_home);
      }
      catch (const std::filesystem::filesystem_error &fse) {
         std::cout << "ERROR:: " << fse.what() << std::endl;
      }
      return state_home;
   }
   std::filesystem::path get_data_home() const {
      try {
         if (!std::filesystem::is_directory(data_home))
            std::filesystem::create_directories(data_home);
      }
      catch (const std::filesystem::filesystem_error &fse) {
         std::cout << "ERROR:: " << fse.what() << std::endl;
      }
      return data_home;
   }
   std::filesystem::path get_config_home() const {
      try {
         if (!std::filesystem::is_directory(config_home))
            std::filesystem::create_directories(config_home);
      }
      catch (const std::filesystem::filesystem_error &fse) {
         std::cout << "ERROR:: " << fse.what() << std::endl;
      }
      return config_home;
   }
   std::filesystem::path get_cache_home() const {
      try {
         if (!std::filesystem::is_directory(cache_home))
            std::filesystem::create_directories(cache_home);
      }
      catch (const std::filesystem::filesystem_error &fse) {
         std::cout << "ERROR:: " << fse.what() << std::endl;
      }
      return cache_home;
   }
   std::filesystem::path get_runtime_dir() const {
      // 20250113-PE add a check here that runtime_dir is on a local filesystem
      // Using runtime_dir is non-trivial. It needs more work.
      // Bottom line is don't use it at the moment.
      try {
         if (!std::filesystem::is_directory(runtime_dir))
            std::filesystem::create_directories(runtime_dir);
      }
      catch (const std::filesystem::filesystem_error &fse) {
         std::cout << "ERROR:: " << fse.what() << std::endl;
      }
      return runtime_dir;
   }
   std::string get_data_dirs() const { return data_dirs; }
   std::string get_config_dirs() const { return config_dirs; }
   std::vector<std::filesystem::path> get_python_config_scripts() const {
      return get_scripts_internal(get_config_home(), ".py"); }
   std::vector<std::filesystem::path> get_scheme_config_scripts() const{
      return get_scripts_internal(get_config_home(), ".scm"); }

   std::pair<bool, std::filesystem::path> get_python_state_script() const {
      bool status = false;
      std::filesystem::path p;
      std::vector<std::filesystem::path> scripts = get_scripts_internal(get_state_home(), ".py");
      for (const auto &script : scripts) {
         if (script.filename() == "0-coot.state.py") {
            status = true;
            p = script;
            break;
         }
      }
      return std::make_pair(status, p);
   }
   std::pair<bool, std::filesystem::path> get_scheme_state_script() const{
      bool status = false;
      std::filesystem::path p;
      std::vector<std::filesystem::path> scripts = get_scripts_internal(get_state_home(), ".scm");
      for (const auto &script : scripts) {
         if (script.filename() == "0-coot.state.scm") {
            status = true;
            p = script;
            break;
         }
      }
      return std::make_pair(status, p);
   }
   std::filesystem::path join(const std::filesystem::path &p, const std::string &file_name) const {
      auto d = p / file_name;
      return d;
   }
   // utility function for coot
   std::string get_download_dir() const {
      std::filesystem::path c = get_cache_home();
      auto d = c / "coot-download";
      return d;
   }
};

#endif // UTILS_XDG_BASE_HH
