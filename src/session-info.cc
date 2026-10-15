
#include <iostream>
#include <fstream>
#include <filesystem>
#include <ctime>

#include "coot-utils/json.hpp"
#include "utils/xdg-base.hh"

#include "session-info.hh"

namespace {

   // Maximum number of entries kept in the recent-files list.
   const std::size_t recent_files_max = 15;

   std::filesystem::path session_json_path() {
      xdg_t xdg;
      // get_state_home() creates ~/.local/state/Coot if needed
      return xdg.get_state_home().append("session.json");
   }

   std::string canonical_path_string(const std::string &path) {
      // The file has just been read, so it should exist - but be defensive
      // and fall back to the absolute (or raw) path if canonicalisation fails.
      try {
         return std::filesystem::weakly_canonical(std::filesystem::absolute(path)).string();
      }
      catch (const std::exception &e) {
         try {
            return std::filesystem::absolute(path).string();
         }
         catch (const std::exception &e2) {
            return path;
         }
      }
   }

   std::string now_iso8601() {
      std::time_t t = std::time(nullptr);
      char buf[32];
      std::tm tm_local;
#if defined(_WIN32)
      localtime_s(&tm_local, &t);
#else
      localtime_r(&t, &tm_local);
#endif
      if (std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm_local) == 0)
         return std::string();
      return std::string(buf);
   }

   // Load the existing document, or return a fresh one. Never throws.
   nlohmann::json load_session_json(const std::filesystem::path &path) {
      nlohmann::json j;
      bool ok = false;
      if (std::filesystem::exists(path)) {
         try {
            std::ifstream f(path);
            if (f) {
               f >> j;
               if (j.is_object() && j.contains("recent_files") && j["recent_files"].is_array())
                  ok = true;
            }
         }
         catch (const std::exception &e) {
            // malformed file - start fresh rather than lose the feature
            std::cout << "WARNING:: session.json parse failed (" << e.what()
                      << ") - starting a fresh recent-files list" << std::endl;
         }
      }
      if (! ok) {
         j = nlohmann::json::object();
         j["schema_version"] = 1;
         j["recent_files"] = nlohmann::json::array();
      }
      return j;
   }

   // Write atomically: write a temp file in the same directory, then rename
   // over the target so a kill mid-write cannot corrupt session.json.
   void write_session_json(const std::filesystem::path &path, const nlohmann::json &j) {
      std::filesystem::path tmp = path;
      tmp += ".tmp";
      try {
         {
            std::ofstream f(tmp);
            if (! f) {
               std::cout << "WARNING:: cannot write " << tmp.string() << std::endl;
               return;
            }
            f << j.dump(2) << std::endl;
         }
         std::filesystem::rename(tmp, path);
      }
      catch (const std::exception &e) {
         std::cout << "WARNING:: failed to write session.json (" << e.what() << ")" << std::endl;
         std::error_code ec;
         std::filesystem::remove(tmp, ec); // best-effort cleanup, ignore error
      }
   }

   bool same_recent_entry(const nlohmann::json &a, const nlohmann::json &b) {
      // Coordinates and dictionaries de-dup on (type, path). Maps also
      // compare the column selection, so the same MTZ opened with different
      // columns (e.g. 2Fo-Fc vs Fo-Fc) keeps a distinct recent entry.
      if (a.value("type", "") != b.value("type", ""))
         return false;
      if (a.value("path", "") != b.value("path", ""))
         return false;
      if (a.value("type", "") == std::string("map")) {
         if (a.value("f_col", "")      != b.value("f_col", ""))      return false;
         if (a.value("phi_col", "")    != b.value("phi_col", ""))    return false;
         if (a.value("weight_col", "") != b.value("weight_col", "")) return false;
         if (a.value("is_diff_map", false) != b.value("is_diff_map", false)) return false;
      }
      if (a.value("type", "") == std::string("ccp4-map")) {
         if (a.value("is_diff_map", false) != b.value("is_diff_map", false)) return false;
      }
      return true;
   }

   // Append entry to the recent-files MRU (most-recent-first), de-duplicated
   // and capped. entry must already carry its "type" and "path".
   void note_recent_entry(nlohmann::json entry) {
      entry["last_used"] = now_iso8601();

      std::filesystem::path path = session_json_path();
      nlohmann::json j = load_session_json(path);

      nlohmann::json &recents = j["recent_files"];

      nlohmann::json rebuilt = nlohmann::json::array();
      rebuilt.push_back(entry);
      for (const auto &e : recents) {
         if (rebuilt.size() >= recent_files_max)
            break;
         if (! same_recent_entry(e, entry))
            rebuilt.push_back(e);
      }
      j["recent_files"] = rebuilt;

      write_session_json(path, j);
   }
}

std::vector<coot::session::recent_file_t>
coot::session::get_recent_files() {

   std::vector<recent_file_t> v;
   std::filesystem::path path = session_json_path();
   nlohmann::json j = load_session_json(path);
   const nlohmann::json &recents = j["recent_files"];
   for (const auto &e : recents) {
      recent_file_t rf;
      rf.type        = e.value("type", "");
      rf.path        = e.value("path", "");
      rf.f_col       = e.value("f_col", "");
      rf.phi_col     = e.value("phi_col", "");
      rf.weight_col  = e.value("weight_col", "");
      rf.use_weights = e.value("use_weights", false);
      rf.is_diff_map = e.value("is_diff_map", false);
      rf.last_used   = e.value("last_used", "");
      if (! rf.type.empty() && ! rf.path.empty())
         v.push_back(rf);
   }
   return v;
}

void
coot::session::note_recent_coordinates(const std::string &path) {

   nlohmann::json entry;
   entry["type"] = "coordinates";
   entry["path"] = canonical_path_string(path);
   note_recent_entry(entry);
}

void
coot::session::note_recent_map(const std::string &mtz_file_name,
                               const std::string &f_col,
                               const std::string &phi_col,
                               const std::string &weight_col,
                               bool use_weights,
                               bool is_diff_map) {

   nlohmann::json entry;
   entry["type"] = "map";
   entry["path"] = canonical_path_string(mtz_file_name);
   entry["f_col"] = f_col;
   entry["phi_col"] = phi_col;
   entry["weight_col"] = use_weights ? weight_col : std::string();
   entry["use_weights"] = use_weights;
   entry["is_diff_map"] = is_diff_map;
   note_recent_entry(entry);
}

void
coot::session::note_recent_ccp4_map(const std::string &path, bool is_diff_map) {

   nlohmann::json entry;
   entry["type"] = "ccp4-map";
   entry["path"] = canonical_path_string(path);
   entry["is_diff_map"] = is_diff_map;
   note_recent_entry(entry);
}

void
coot::session::note_recent_dictionary(const std::string &path) {

   nlohmann::json entry;
   entry["type"] = "dictionary";
   entry["path"] = canonical_path_string(path);
   note_recent_entry(entry);
}
