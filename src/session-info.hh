
#ifndef SESSION_INFO_HH
#define SESSION_INFO_HH

#include <string>
#include <vector>

// Recent-files store for the "Open Recent" GUI. Each time a file is read,
// an entry is appended (most-recent-first) to the recent-files list in
// ~/.local/state/Coot/session.json (XDG state home). The list is capped,
// de-duplicated by canonical path, and written atomically. For maps the
// MTZ column selection is recorded, so a map can later be reopened without
// re-prompting for columns.
//
// These functions are safe to call even if the write fails (e.g. a
// read-only state dir): they log and return rather than throw. The caller
// is expected to call only on a successful read, and to honour
// graphics_info_t::disable_state_script_writing at the call site.

namespace coot {
   namespace session {

      // A parsed recent-files entry. For coordinates and dictionaries only
      // type and path are meaningful; for maps the MTZ column selection is
      // filled in too.
      class recent_file_t {
      public:
         std::string type;        // "coordinates", "map" or "dictionary"
         std::string path;
         std::string f_col;
         std::string phi_col;
         std::string weight_col;
         bool use_weights = false;
         bool is_diff_map = false;
         std::string last_used;   // ISO-8601, informational
         recent_file_t() {}
         bool is_coordinates() const { return type == "coordinates"; }
         bool is_map()         const { return type == "map"; }       // from MTZ, with columns
         bool is_ccp4_map()    const { return type == "ccp4-map"; }   // a CCP4/MRC map file
         bool is_dictionary()  const { return type == "dictionary"; }
      };

      // The recent-files list, most-recent-first, as stored in session.json.
      // Returns an empty vector if the file is absent or unparseable.
      std::vector<recent_file_t> get_recent_files();

      void note_recent_coordinates(const std::string &path);

      void note_recent_map(const std::string &mtz_file_name,
                           const std::string &f_col,
                           const std::string &phi_col,
                           const std::string &weight_col,
                           bool use_weights,
                           bool is_diff_map);

      void note_recent_ccp4_map(const std::string &path, bool is_diff_map);

      void note_recent_dictionary(const std::string &path);

   }
}

#endif // SESSION_INFO_HH
