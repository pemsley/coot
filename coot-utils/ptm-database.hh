

#ifndef COOT_UTILS_PTM_DATABASE_HH
#define COOT_UTILS_PTM_DATABASE_HH

#include <string>
#include <map>
#include <vector>
#include <optional>

namespace coot {

   // A single PTM entry describing a modification
   struct ptm_entry_t {
      std::string modified_residue_name;  // e.g. "ACE", "SEP"
      std::string category;               // e.g. "Terminal acetylation", "Named protein modification"
      std::string type;                   // e.g. "Methylation", "Phosphorylation", "None"
      std::string position;               // e.g. "Amino-acid backbone", "Amino-acid side chain"
      std::string location;               // e.g. "N-terminal", "C-terminal", "Any position"

      ptm_entry_t() = default;
      ptm_entry_t(const std::string &mod_res, const std::string &cat,
                  const std::string &typ, const std::string &pos,
                  const std::string &loc)
         : modified_residue_name(mod_res), category(cat), type(typ),
           position(pos), location(loc) {}
   };

   // Database storing all PTM entries, indexed by parent residue name
   class ptm_database_t {
   public:
      ptm_database_t() = default;

      // Parse the PTM database from a JSON file
      // Returns true on success, false on failure
      bool read(const std::string &file_name);

      // Get all PTM entries for a given parent residue (e.g. "SER")
      // Returns empty vector if residue not found
      std::vector<ptm_entry_t> get_entries_for_residue(const std::string &parent_residue) const;

      // Look up a specific modification by parent residue and modified residue name
      // Returns nullopt if not found
      std::optional<ptm_entry_t> get_entry(const std::string &parent_residue,
                                           const std::string &modified_residue) const;

      // Get all parent residue names in the database
      std::vector<std::string> get_parent_residues() const;

      // Get total number of PTM entries
      std::size_t size() const;

      // Check if database is empty
      bool empty() const { return data_.empty(); }

   private:
      // Map from parent residue name (e.g. "SER") to map of
      // modified residue name (e.g. "SEP") to PTM entry
      std::map<std::string, std::map<std::string, ptm_entry_t>> data_;
   };

} // namespace coot

#endif // COOT_UTILS_PTM_DATABASE_HH
