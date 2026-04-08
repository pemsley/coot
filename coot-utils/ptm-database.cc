

#include <fstream>
#include <iostream>
#include "ptm-database.hh"
#include "json.hpp"

using json = nlohmann::json;

bool coot::ptm_database_t::read(const std::string &file_name) {

   std::ifstream f(file_name);
   if (!f.is_open()) {
      std::cerr << "ERROR: ptm_database_t::read() failed to open " << file_name << std::endl;
      return false;
   }

   try {
      json j = json::parse(f);

      // Structure: { "ALA": { "ACE": { "category": ..., "type": ..., ... }, ... }, ... }
      for (auto &[parent_residue, modifications] : j.items()) {
         std::map<std::string, ptm_entry_t> mod_map;
         for (auto &[mod_residue, props] : modifications.items()) {
            ptm_entry_t entry;
            entry.modified_residue_name = mod_residue;
            entry.category = props.value("category", "");
            entry.type     = props.value("type", "");
            entry.position = props.value("position", "");
            entry.location = props.value("location", "");
            mod_map[mod_residue] = entry;
         }
         data_[parent_residue] = mod_map;
      }
      return true;

   } catch (const json::parse_error &e) {
      std::cerr << "ERROR: ptm_database_t::read() JSON parse error: " << e.what() << std::endl;
      return false;
   } catch (const std::exception &e) {
      std::cerr << "ERROR: ptm_database_t::read() exception: " << e.what() << std::endl;
      return false;
   }
}

std::vector<coot::ptm_entry_t>
coot::ptm_database_t::get_entries_for_residue(const std::string &parent_residue) const {

   std::vector<ptm_entry_t> result;
   auto it = data_.find(parent_residue);
   if (it != data_.end()) {
      for (const auto &[mod_name, entry] : it->second) {
         result.push_back(entry);
      }
   }
   return result;
}

std::optional<coot::ptm_entry_t>
coot::ptm_database_t::get_entry(const std::string &parent_residue,
                          const std::string &modified_residue) const {

   auto parent_it = data_.find(parent_residue);
   if (parent_it != data_.end()) {
      auto mod_it = parent_it->second.find(modified_residue);
      if (mod_it != parent_it->second.end()) {
         return mod_it->second;
      }
   }
   return std::nullopt;
}

std::vector<std::string>
coot::ptm_database_t::get_parent_residues() const {

   std::vector<std::string> result;
   result.reserve(data_.size());
   for (const auto &[parent, mods] : data_) {
      result.push_back(parent);
   }
   return result;
}

std::size_t
coot::ptm_database_t::size() const {

   std::size_t count = 0;
   for (const auto &[parent, mods] : data_) {
      count += mods.size();
   }
   return count;
}

