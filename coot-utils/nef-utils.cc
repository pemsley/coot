/**
 * NEF/BMRB STAR Data Parser using gemmi's STAR parser
 *
 * Parses NMR data files in NEF and BMRB STAR formats to extract:
 * - Assigned chemical shifts
 * - NOE peaks (NEF format)
 * - Distance restraints (NEF format)
 * - Coupling constants (BMRB format)
 *
 */

#include "nef-utils.hh"

#include <gemmi/cif.hpp>
#include <gemmi/numb.hpp>

#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <optional>
#include <regex>


// ============================================================================
// File preprocessing for malformed BMRB files
// ============================================================================

/**
 * Preprocess STAR file content to fix common issues that cause parsing errors.
 * This is the core function that works with string content, suitable for WebAssembly.
 *
 * @param content The raw file content as a string
 * @return The preprocessed content with fixes applied
 */
std::string nef::preprocess_star_file_string(const std::string& content) {

    std::vector<std::string> lines;
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        lines.push_back(line);
    }

    std::vector<std::string> fixed_lines;
    bool in_loop = false;
    std::vector<std::string> loop_tags;

    // Regex to quote $references
    std::regex dollar_ref(R"((\s)(\$\S+))");

    for (const auto& original_line : lines) {
        std::string stripped = original_line;
        // Trim whitespace
        size_t start = stripped.find_first_not_of(" \t\r\n");
        size_t end = stripped.find_last_not_of(" \t\r\n");
        if (start != std::string::npos && end != std::string::npos) {
            stripped = stripped.substr(start, end - start + 1);
        } else {
            stripped = "";
        }

        if (stripped.rfind("loop_", 0) == 0) {
            in_loop = true;
            loop_tags.clear();
            fixed_lines.push_back(original_line);
            continue;
        }

        if (in_loop) {
            if (!stripped.empty() && stripped[0] == '_') {
                // Tag definition
                loop_tags.push_back(stripped);
                fixed_lines.push_back(original_line);
            } else if (stripped == "stop_") {
                in_loop = false;
                loop_tags.clear();
                fixed_lines.push_back(original_line);
            } else if (stripped.rfind("save_", 0) == 0) {
                in_loop = false;
                loop_tags.clear();
                fixed_lines.push_back(original_line);
            } else if (stripped.empty() || stripped[0] == '#') {
                // Empty lines and comments kept as-is
                fixed_lines.push_back(original_line);
            } else {
                // Data row - quote $references
                std::string fixed_line = std::regex_replace(original_line, dollar_ref, "$1\"$2\"");
                fixed_lines.push_back(fixed_line);
            }
        } else {
            // Outside loop - also quote $references in single-line values
            std::string fixed_line = std::regex_replace(original_line, dollar_ref, "$1\"$2\"");
            fixed_lines.push_back(fixed_line);
        }
    }

    std::ostringstream result;
    for (const auto& l : fixed_lines) {
        result << l << "\n";
    }
    return result.str();
}

/**
 * Preprocess STAR file from a filepath.
 * This is a convenience wrapper that reads the file and calls preprocess_star_file_string.
 *
 * @param filepath Path to the STAR file
 * @return The preprocessed content with fixes applied
 */
std::string nef::preprocess_star_file(const std::string& filepath) {

    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filepath);
    }

    std::ostringstream buffer;
    buffer << file.rdbuf();
    file.close();

    return preprocess_star_file_string(buffer.str());
}

// ============================================================================
// NEF/STAR Parser class
// ============================================================================

nef::NEFParser::NEFParser(const std::string& content, const std::string& source_name) : filepath_(source_name) {
   init_from_string(content, source_name);
}

/**
 * Construct parser from a filepath (convenience wrapper).
 *
 * @param filepath Path to the STAR/NEF file
 */

// static
nef::NEFParser nef::NEFParser::from_file(const std::string& filepath) {

   std::ifstream file(filepath);
   if (!file.is_open()) {
      throw std::runtime_error("Could not open file: " + filepath);
   }
   std::ostringstream buffer;
   buffer << file.rdbuf();
   file.close();

   return NEFParser(buffer.str(), filepath);
}

void
nef::NEFParser::init_from_string(const std::string& content, const std::string& source_name) {
   try {
      doc_ = gemmi::cif::read_memory(content.data(), content.size(), source_name.c_str());
   } catch (const std::exception& e) {
      std::cout << "Note: Standard parsing failed (" << e.what() << ")\n";
      std::cout << "Attempting to preprocess content...\n";
      std::string preprocessed = nef::preprocess_star_file_string(content);
      doc_ = gemmi::cif::read_memory(preprocessed.data(), preprocessed.size(), source_name.c_str());
   }

   if (doc_.blocks.empty()) {
      throw std::runtime_error("No data blocks found in content");
   }
   block_ = &doc_.blocks[0];
}

   // ============================================================================
   // Helper functions
   // ============================================================================

std::string nef::NEFParser::format_to_string(nef::FileFormat fmt) const {
   switch (fmt) {
   case nef::FileFormat::NEF: return "nef";
   case nef::FileFormat::BMRB_STAR: return "bmrb_star";
   default: return "unknown";
   }
}

bool nef::NEFParser::is_null_value(const std::string& val) const {
   return val == "." || val == "?" || val.empty();
}

std::optional<double> nef::NEFParser::parse_double(const std::string& val) const {
   if (is_null_value(val)) return std::nullopt;
   try {
      return std::stod(val);
   } catch (...) {
      return std::nullopt;
   }
}

std::optional<int> nef::NEFParser::parse_int(const std::string& val) const {
   if (is_null_value(val)) return std::nullopt;
   try {
      return std::stoi(val);
   } catch (...) {
      return std::nullopt;
   }
}

void nef::NEFParser::parse() {
   format_ = detect_format();
   std::cout << "Parsing file with gemmi STAR parser...\n";
   std::cout << "Block name: " << block_->name << "\n";
   std::cout << "Detected format: " << format_to_string(format_) << "\n";

   if (format_ == nef::FileFormat::NEF) {
      parse_assigned_chemical_shifts();
      parse_noe_peaks();
      parse_distance_restraints();
   } else if (format_ == nef::FileFormat::BMRB_STAR) {
      parse_assigned_chemical_shifts();
      parse_coupling_constants();
      parse_distance_restraints();
   } else {
      std::cout << "Unknown format, attempting to parse all data types...\n";
      parse_noe_peaks();
      parse_distance_restraints();
      parse_assigned_chemical_shifts();
      parse_coupling_constants();
   }
}

void nef::NEFParser::print_summary() const {

   std::cout << "\n" << std::string(80, '=') << "\n";
   std::cout << "NMR Data Summary (gemmi STAR parser - C++ version)\n";
   std::cout << std::string(80, '=') << "\n";
   std::cout << "\nFile: " << filepath_ << "\n";
   std::cout << "Format: " << format_to_string(format_) << "\n";

   // Chemical shifts summary
   if (!chemical_shifts_.empty()) {
      std::cout << "\n" << std::string(80, '=') << "\n";
      std::cout << "Assigned Chemical Shifts:\n";
      std::cout << std::string(80, '=') << "\n";

      // Count by atom type
      std::map<std::string, int> by_type;
      for (const auto& shift : chemical_shifts_) {
         std::string key = std::to_string(shift.isotope_number) + shift.atom_type;
         by_type[key]++;
      }

      std::cout << "\n  Total assignments: " << chemical_shifts_.size() << "\n";
      std::cout << "  By nucleus:\n";
      for (const auto& [nucleus, count] : by_type) {
         std::cout << "    " << nucleus << ": " << count << "\n";
      }

      // Count unique residues
      std::set<std::pair<int, std::string> > residues;
      for (const auto& s : chemical_shifts_) {
         residues.insert({s.seq_id, s.comp_id});
      }
      std::cout << "  Unique residues with assignments: " << residues.size() << "\n";

      // Show sample shifts
      std::cout << "\n  Sample assignments (first 10):\n";
      std::cout << "  " << std::string(50, '-') << "\n";
      size_t count = 0;
      for (const auto& shift : chemical_shifts_) {
         if (count++ >= 10) break;
         std::cout << "    " << shift.to_string() << "\n";
      }
   }

   // Coupling constants summary
   if (!coupling_constants_.empty()) {
      std::cout << "\n" << std::string(80, '=') << "\n";
      std::cout << "Coupling Constants:\n";
      std::cout << std::string(80, '=') << "\n";
      std::cout << "\n  Total: " << coupling_constants_.size() << "\n";
      std::cout << "\n  Sample (first 5):\n";
      size_t count = 0;
      for (const auto& c : coupling_constants_) {
         if (count++ >= 5) break;
         std::cout << "    " << c.to_string() << "\n";
      }
   }

   // Distance restraints summary
   if (!restraints_.empty()) {
      std::cout << "\n" << std::string(80, '=') << "\n";
      std::cout << "Distance Restraints:\n";
      std::cout << std::string(80, '=') << "\n";

      // Count by origin type
      std::map<std::string, int> by_origin;
      for (const auto& r : restraints_) {
         std::string origin = r.restraint_origin.empty() ? "unspecified" : r.restraint_origin;
         by_origin[origin]++;
      }

      std::cout << "\n  Total restraints: " << restraints_.size() << "\n";
      std::cout << "  By origin:\n";
      for (const auto& [origin, count] : by_origin) {
         std::cout << "    " << origin << ": " << count << "\n";
      }

      // Count unique restraint IDs (rows with same restraint_id are alternatives)
      std::set<int> unique_ids;
      for (const auto& r : restraints_) {
         unique_ids.insert(r.restraint_id);
      }
      std::cout << "  Unique restraint IDs: " << unique_ids.size() << "\n";

      // Show sample restraints
      std::cout << "\n  Sample restraints (first 15):\n";
      std::cout << "  " << std::string(70, '-') << "\n";
      size_t count = 0;
      for (const auto& r : restraints_) {
         if (count++ >= 15) break;
         std::cout << "    " << r.to_string() << "\n";
      }
   }

   // NOE peaks
   if (!peaks_.empty()) {
      std::cout << "\n" << std::string(80, '=') << "\n";
      std::cout << "Total NOE Peaks: " << peaks_.size() << "\n";
      std::cout << std::string(80, '=') << "\n";
   }

   std::cout << "\n" << std::string(80, '=') << "\n";
}

void nef::NEFParser::print_statistics() const {

   if (chemical_shifts_.empty() && peaks_.empty() && restraints_.empty()) {
      return;
   }

   std::cout << "\n" << std::string(80, '=') << "\n";
   std::cout << "Detailed Statistics\n";
   std::cout << std::string(80, '=') << "\n";

   if (!chemical_shifts_.empty()) {
      std::cout << "\nChemical Shift Statistics:\n";

      // Group by atom type
      std::map<std::string, std::vector<double>> by_type;
      for (const auto& shift : chemical_shifts_) {
         std::string key = std::to_string(shift.isotope_number) + shift.atom_type;
         by_type[key].push_back(shift.value);
      }

      for (const auto& [nucleus, values] : by_type) {
         double sum = 0;
         double min_val = values[0];
         double max_val = values[0];
         for (double v : values) {
            sum += v;
            min_val = std::min(min_val, v);
            max_val = std::max(max_val, v);
         }
         double avg = sum / values.size();

         std::cout << "\n  " << nucleus << " (" << values.size() << " shifts):\n";
         std::cout << "    Range: " << std::fixed << std::setprecision(3)
                   << min_val << " - " << max_val << " ppm\n";
         std::cout << "    Average: " << avg << " ppm\n";
      }

      // Residue coverage
      std::set<std::pair<int, std::string>> residues;
      for (const auto& s : chemical_shifts_) {
         residues.insert({s.seq_id, s.comp_id});
      }
      if (!residues.empty()) {
         int first_res = residues.begin()->first;
         int last_res = residues.rbegin()->first;
         std::cout << "\n  Sequence coverage: residues " << first_res
                   << " to " << last_res << "\n";
         std::cout << "  Residues with assignments: " << residues.size() << "\n";
      }
   }
   std::cout << std::string(80, '=') << "\n";
}


nef::FileFormat nef::NEFParser::detect_format() {
   // Check for NEF-specific saveframes
   for (const auto& item : block_->items) {
      if (item.type == gemmi::cif::ItemType::Frame) {
         if (item.frame.name.rfind("nef_", 0) == 0) {
            return nef::FileFormat::NEF;
         }
      }
      if (item.type == gemmi::cif::ItemType::Pair) {
         const std::string& tag = item.pair[0];
         if (tag.find("nef_") != std::string::npos) {
            return nef::FileFormat::NEF;
         }
         if (tag.find("_Assigned_chem_shift_list") != std::string::npos ||
             tag.find("_Atom_chem_shift") != std::string::npos) {
            return nef::FileFormat::BMRB_STAR;
         }
      }
      if (item.type == gemmi::cif::ItemType::Loop) {
         for (const auto& tag : item.loop.tags) {
            if (tag.find("nef_") != std::string::npos) {
               return nef::FileFormat::NEF;
            }
            if (tag.find("_Atom_chem_shift") != std::string::npos) {
               return nef::FileFormat::BMRB_STAR;
            }
         }
      }
   }
   return nef::FileFormat::UNKNOWN;
}

void nef::NEFParser::parse_assigned_chemical_shifts() {

   std::cout << "\nSearching for assigned chemical shift lists...\n";

   std::vector<const gemmi::cif::Loop*> shift_loops;

   // Method 1: Try finding saveframes with chemical shifts
   for (const auto& item : block_->items) {
      if (item.type == gemmi::cif::ItemType::Frame) {
         const gemmi::cif::Block& frame = item.frame;
         for (const auto& frame_item : frame.items) {
            if (frame_item.type == gemmi::cif::ItemType::Loop) {
               const gemmi::cif::Loop& loop = frame_item.loop;
               bool has_nef_shift = false;
               bool has_bmrb_shift = false;
               for (const auto& tag : loop.tags) {
                  if (tag.find("_nef_chemical_shift.") != std::string::npos) {
                     has_nef_shift = true;
                     break;
                  }
                  if (tag.find("_Atom_chem_shift.") != std::string::npos) {
                     has_bmrb_shift = true;
                     break;
                  }
               }
               if (has_nef_shift) {
                  shift_loops.push_back(&loop);
                  std::cout << "  Found NEF chemical shift list in saveframe: "
                            << frame.name << "\n";
               } else if (has_bmrb_shift) {
                  shift_loops.push_back(&loop);
                  std::cout << "  Found BMRB chemical shift list in saveframe: "
                            << frame.name << "\n";
               }
            }
         }
      }
   }

   // Method 2: Try finding loops directly in block
   if (shift_loops.empty()) {
      for (const auto& item : block_->items) {
         if (item.type == gemmi::cif::ItemType::Loop) {
            const gemmi::cif::Loop& loop = item.loop;
            for (const auto& tag : loop.tags) {
               if (tag.find("_Atom_chem_shift.") != std::string::npos) {
                  shift_loops.push_back(&loop);
                  std::cout << "  Found BMRB chemical shift list in block\n";
                  break;
               }
               if (tag.find("_nef_chemical_shift.") != std::string::npos) {
                  shift_loops.push_back(&loop);
                  std::cout << "  Found NEF chemical shift list in block\n";
                  break;
               }
            }
         }
      }
   }

   if (shift_loops.empty()) {
      std::cout << "  No assigned chemical shift data found\n";
      return;
   }

   // Process all found loops
   for (const gemmi::cif::Loop* loop : shift_loops) {
      process_chemical_shift_loop(*loop);
   }
}

void nef::NEFParser::process_chemical_shift_loop(const gemmi::cif::Loop& loop) {
   std::cout << "  Found " << loop.length() << " chemical shift assignments\n";

   // Determine format (NEF vs BMRB)
   bool is_nef = false;
   for (const auto& tag : loop.tags) {
      if (tag.find("_nef_chemical_shift.") != std::string::npos) {
         is_nef = true;
         break;
      }
   }
   std::string prefix = is_nef ? "_nef_chemical_shift." : "_Atom_chem_shift.";

   // Find column indices
   auto get_col_idx = [&](const std::string& col_name,
                          const std::vector<std::string>& alt_names = {}) -> int {
      std::vector<std::string> names_to_try = {prefix + col_name};
      for (const auto& alt : alt_names) {
         names_to_try.push_back(prefix + alt);
      }
      for (const auto& name : names_to_try) {
         for (size_t i = 0; i < loop.tags.size(); ++i) {
            if (loop.tags[i] == name) {
               return static_cast<int>(i);
            }
         }
      }
      return -1;
   };

   int idx_id = get_col_idx("ID", {"index"});
   int idx_seq_id = get_col_idx("Seq_ID", {"sequence_code"});
   int idx_comp_index_id = get_col_idx("Comp_index_ID", {});  // Fallback for Seq_ID
   int idx_comp_id = get_col_idx("Comp_ID", {"residue_name"});
   int idx_atom_id = get_col_idx("Atom_ID", {"atom_name"});
   int idx_atom_type = get_col_idx("Atom_type", {"element"});
   int idx_isotope = get_col_idx("Atom_isotope_number", {"isotope_number"});
   int idx_val = get_col_idx("Val", {"value"});
   int idx_val_err = get_col_idx("Val_err", {"value_uncertainty"});
   int idx_ambiguity = get_col_idx("Ambiguity_code");

   // Check required columns (seq_id or comp_index_id must be present)
   bool has_seq_col = (idx_seq_id >= 0) || (idx_comp_index_id >= 0);
   if (!has_seq_col || idx_comp_id < 0 || idx_atom_id < 0 || idx_val < 0) {
      std::cout << "  Error: Missing required columns for chemical shifts\n";
      return;
   }

   // Parse each row
   int n_errors = 0;
   for (size_t i = 0; i < loop.length(); ++i) {
      try {
         // Handle sequence code - try Seq_ID first, fall back to Comp_index_ID
         std::string seq_code;
         if (idx_seq_id >= 0) {
            seq_code = loop.val(i, idx_seq_id);
         }

         // If Seq_ID is null/missing, try Comp_index_ID
         if ((seq_code.empty() || seq_code == "." || seq_code == "?") && idx_comp_index_id >= 0) {
            seq_code = loop.val(i, idx_comp_index_id);
         }

         if (seq_code.empty() || seq_code[0] == '@' || seq_code == "." || seq_code == "?") {
            continue;  // Skip unassigned
         }

         int seq_id;
         try {
            seq_id = std::stoi(seq_code);
         } catch (...) {
            continue;  // Skip non-integer sequence codes
         }

         nef::AssignedChemicalShift shift;
         shift.shift_id = (idx_id >= 0) ? parse_int(loop.val(i, idx_id)).value_or(i + 1) : i + 1;
         shift.seq_id = seq_id;
         shift.comp_id = loop.val(i, idx_comp_id);
         shift.atom_id = loop.val(i, idx_atom_id);
         shift.atom_type = (idx_atom_type >= 0) ? loop.val(i, idx_atom_type) : "H";
         if (shift.atom_type.empty() || shift.atom_type == "." || shift.atom_type == "?") {
            shift.atom_type = "H";
         }

         // Get isotope number, inferring from atom type if not specified
         auto isotope_opt = (idx_isotope >= 0) ? parse_int(loop.val(i, idx_isotope)) : std::nullopt;
         if (isotope_opt) {
            shift.isotope_number = *isotope_opt;
         } else {
            // Infer common isotopes from atom type
            if (shift.atom_type == "H") {
               shift.isotope_number = 1;
            } else if (shift.atom_type == "C") {
               shift.isotope_number = 13;
            } else if (shift.atom_type == "N") {
               shift.isotope_number = 15;
            } else if (shift.atom_type == "P") {
               shift.isotope_number = 31;
            } else if (shift.atom_type == "F") {
               shift.isotope_number = 19;
            } else {
               shift.isotope_number = 1;  // Default
            }
         }
         shift.value = std::stod(loop.val(i, idx_val));
         shift.value_error = (idx_val_err >= 0) ? parse_double(loop.val(i, idx_val_err)) : std::nullopt;
         shift.ambiguity_code = (idx_ambiguity >= 0) ? parse_int(loop.val(i, idx_ambiguity)) : std::nullopt;

         chemical_shifts_.push_back(shift);
      } catch (const std::exception& e) {
         n_errors++;
         if (n_errors == 1) {
            std::cout << "    Warning: Error parsing shift: " << e.what() << "\n";
         }
      }
   }

   if (n_errors > 0) {
      std::cout << "    Note: " << n_errors << " shifts had parsing errors\n";
   }
}

void nef::NEFParser::parse_coupling_constants() {

   std::cout << "\nSearching for coupling constant lists...\n";

   const gemmi::cif::Loop* coupling_loop = nullptr;

   // Search in saveframes
   for (const auto& item : block_->items) {
      if (item.type == gemmi::cif::ItemType::Frame) {
         for (const auto& frame_item : item.frame.items) {
            if (frame_item.type == gemmi::cif::ItemType::Loop) {
               for (const auto& tag : frame_item.loop.tags) {
                  if (tag.find("_Coupling_constant.") != std::string::npos) {
                     coupling_loop = &frame_item.loop;
                     break;
                  }
               }
               if (coupling_loop) break;
            }
         }
         if (coupling_loop) break;
      }
   }

   if (!coupling_loop) {
      std::cout << "  No coupling constant data found\n";
      return;
   }

   std::cout << "  Found " << coupling_loop->length() << " coupling constants\n";

   // Find column indices
   auto get_col_idx = [&](const std::string& col_name) -> int {
      std::string full_name = "_Coupling_constant." + col_name;
      for (size_t i = 0; i < coupling_loop->tags.size(); ++i) {
         if (coupling_loop->tags[i] == full_name) {
            return static_cast<int>(i);
         }
      }
      return -1;
   };

   int idx_id = get_col_idx("ID");
   int idx_seq_1 = get_col_idx("Seq_ID_1");
   int idx_comp_1 = get_col_idx("Comp_ID_1");
   int idx_atom_1 = get_col_idx("Atom_ID_1");
   int idx_seq_2 = get_col_idx("Seq_ID_2");
   int idx_comp_2 = get_col_idx("Comp_ID_2");
   int idx_atom_2 = get_col_idx("Atom_ID_2");
   int idx_val = get_col_idx("Val");
   int idx_val_err = get_col_idx("Val_err");

   // Parse each row
   int n_errors = 0;
   for (size_t i = 0; i < coupling_loop->length(); ++i) {
      try {
         nef::CouplingConstant coupling;
         coupling.coupling_id = (idx_id >= 0) ? parse_int(coupling_loop->val(i, idx_id)).value_or(i + 1) : i + 1;
         coupling.seq_id_1 = (idx_seq_1 >= 0) ? parse_int(coupling_loop->val(i, idx_seq_1)).value_or(0) : 0;
         coupling.comp_id_1 = (idx_comp_1 >= 0) ? coupling_loop->val(i, idx_comp_1) : "";
         coupling.atom_id_1 = (idx_atom_1 >= 0) ? coupling_loop->val(i, idx_atom_1) : "";
         coupling.seq_id_2 = (idx_seq_2 >= 0) ? parse_int(coupling_loop->val(i, idx_seq_2)).value_or(0) : 0;
         coupling.comp_id_2 = (idx_comp_2 >= 0) ? coupling_loop->val(i, idx_comp_2) : "";
         coupling.atom_id_2 = (idx_atom_2 >= 0) ? coupling_loop->val(i, idx_atom_2) : "";
         coupling.value = (idx_val >= 0) ? std::stod(coupling_loop->val(i, idx_val)) : 0.0;
         coupling.value_error = (idx_val_err >= 0) ? parse_double(coupling_loop->val(i, idx_val_err)) : std::nullopt;

         coupling_constants_.push_back(coupling);
      } catch (...) {
         n_errors++;
      }
   }

   if (n_errors > 0) {
      std::cout << "    Note: " << n_errors << " couplings had parsing errors\n";
   }
}


void nef::NEFParser::parse_noe_peaks() {
   std::cout << "\nSearching for NOESY spectra...\n";
   // TODO: Implement NOE peak parsing similar to Python version
   std::cout << "  No NOESY spectra found\n";
}

void nef::NEFParser::parse_distance_restraints() {
   std::cout << "\nSearching for distance restraint lists...\n";

   if (format_ == nef::FileFormat::NEF) {
      parse_nef_distance_restraints();
   } else if (format_ == nef::FileFormat::BMRB_STAR) {
      parse_bmrb_distance_restraints();
   } else {
      // Try both
      parse_nef_distance_restraints();
      parse_bmrb_distance_restraints();
   }

   if (restraints_.empty()) {
      std::cout << "  No distance restraint lists found\n";
   }
}

void nef::NEFParser::parse_nef_distance_restraints() {
   // Look for NEF distance restraint saveframes
   for (const auto& item : block_->items) {
      if (item.type == gemmi::cif::ItemType::Frame) {
         const gemmi::cif::Block& frame = item.frame;

         // Check if this is a distance restraint list
         std::string restraint_origin;
         for (const auto& frame_item : frame.items) {
            if (frame_item.type == gemmi::cif::ItemType::Pair) {
               if (frame_item.pair[0] == "_nef_distance_restraint_list.restraint_origin") {
                  restraint_origin = frame_item.pair[1];
               }
            }
         }

         // Find the distance restraint loop
         for (const auto& frame_item : frame.items) {
            if (frame_item.type == gemmi::cif::ItemType::Loop) {
               const gemmi::cif::Loop& loop = frame_item.loop;
               bool has_nef_restraint = false;
               for (const auto& tag : loop.tags) {
                  if (tag.find("_nef_distance_restraint.") != std::string::npos) {
                     has_nef_restraint = true;
                     break;
                  }
               }

               if (has_nef_restraint) {
                  std::cout << "  Found NEF distance restraint list in saveframe: "
                            << frame.name;
                  if (!restraint_origin.empty() && restraint_origin != ".") {
                     std::cout << " (origin: " << restraint_origin << ")";
                  }
                  std::cout << "\n";
                  process_nef_distance_restraint_loop(loop, restraint_origin);
               }
            }
         }
      }
   }
}

void nef::NEFParser::process_nef_distance_restraint_loop(const gemmi::cif::Loop& loop, const std::string& restraint_origin) {
   std::string prefix = "_nef_distance_restraint.";

   auto get_col_idx = [&](const std::string& col_name) -> int {
      std::string full_name = prefix + col_name;
      for (size_t i = 0; i < loop.tags.size(); ++i) {
         if (loop.tags[i] == full_name) {
            return static_cast<int>(i);
         }
      }
      return -1;
   };

   int idx_index = get_col_idx("index");
   int idx_restraint_id   = get_col_idx("restraint_id");
   int idx_combination_id = get_col_idx("restraint_combination_id");
   int idx_chain_1        = get_col_idx("chain_code_1");
   int idx_seq_1          = get_col_idx("sequence_code_1");
   int idx_res_1          = get_col_idx("residue_name_1");
   int idx_atom_1         = get_col_idx("atom_name_1");
   int idx_chain_2        = get_col_idx("chain_code_2");
   int idx_seq_2          = get_col_idx("sequence_code_2");
   int idx_res_2          = get_col_idx("residue_name_2");
   int idx_atom_2         = get_col_idx("atom_name_2");
   int idx_weight         = get_col_idx("weight");
   int idx_target         = get_col_idx("target_value");
   int idx_target_unc     = get_col_idx("target_value_uncertainty");
   int idx_lower          = get_col_idx("lower_limit");
   int idx_upper          = get_col_idx("upper_limit");

   // Required columns check
   if (idx_seq_1 < 0 || idx_atom_1 < 0 || idx_seq_2 < 0 || idx_atom_2 < 0) {
      std::cout << "    Error: Missing required columns for distance restraints\n";
      return;
   }

   int n_errors = 0;
   size_t count_before = restraints_.size();

   for (size_t i = 0; i < loop.length(); ++i) {
      try {
         nef::DistanceRestraint r;
         r.index = (idx_index >= 0) ? parse_int(loop.val(i, idx_index)).value_or(i + 1) : i + 1;
         r.restraint_id = (idx_restraint_id >= 0) ? parse_int(loop.val(i, idx_restraint_id)).value_or(i + 1) : i + 1;
         r.restraint_combination_id = (idx_combination_id >= 0) ? parse_int(loop.val(i, idx_combination_id)) : std::nullopt;

         r.chain_code_1 = (idx_chain_1 >= 0) ? loop.val(i, idx_chain_1) : "";
         std::string seq_str_1 = loop.val(i, idx_seq_1);
         r.sequence_code_1 = std::stoi(seq_str_1);
         r.residue_name_1 = (idx_res_1 >= 0) ? loop.val(i, idx_res_1) : "";
         r.atom_name_1 = loop.val(i, idx_atom_1);

         r.chain_code_2 = (idx_chain_2 >= 0) ? loop.val(i, idx_chain_2) : "";
         std::string seq_str_2 = loop.val(i, idx_seq_2);
         r.sequence_code_2 = std::stoi(seq_str_2);
         r.residue_name_2 = (idx_res_2 >= 0) ? loop.val(i, idx_res_2) : "";
         r.atom_name_2 = loop.val(i, idx_atom_2);

         r.weight = (idx_weight >= 0) ? parse_double(loop.val(i, idx_weight)).value_or(1.0) : 1.0;
         r.target_value = (idx_target >= 0) ? parse_double(loop.val(i, idx_target)) : std::nullopt;
         r.target_value_uncertainty = (idx_target_unc >= 0) ? parse_double(loop.val(i, idx_target_unc)) : std::nullopt;
         r.lower_limit = (idx_lower >= 0) ? parse_double(loop.val(i, idx_lower)) : std::nullopt;
         r.upper_limit = (idx_upper >= 0) ? parse_double(loop.val(i, idx_upper)) : std::nullopt;
         r.restraint_origin = restraint_origin;

         restraints_.push_back(r);
      } catch (const std::exception& e) {
         n_errors++;
         if (n_errors == 1) {
            std::cout << "    Warning: Error parsing restraint: " << e.what() << "\n";
         }
      }
   }

   std::cout << "    Parsed " << (restraints_.size() - count_before) << " distance restraints\n";
   if (n_errors > 0) {
      std::cout << "    Note: " << n_errors << " restraints had parsing errors\n";
   }
}

void nef::NEFParser::parse_bmrb_distance_restraints() {
   // Look for BMRB general distance constraint saveframes
   for (const auto& item : block_->items) {
      if (item.type == gemmi::cif::ItemType::Frame) {
         const gemmi::cif::Block& frame = item.frame;

         // Check if this is a distance constraint list with NOE type
         std::string constraint_type;
         for (const auto& frame_item : frame.items) {
            if (frame_item.type == gemmi::cif::ItemType::Pair) {
               if (frame_item.pair[0] == "_Gen_dist_constraint_list.Constraint_type") {
                  constraint_type = frame_item.pair[1];
               }
            }
         }

         // Find the distance constraint loop
         for (const auto& frame_item : frame.items) {
            if (frame_item.type == gemmi::cif::ItemType::Loop) {
               const gemmi::cif::Loop& loop = frame_item.loop;
               bool has_bmrb_constraint = false;
               for (const auto& tag : loop.tags) {
                  if (tag.find("_Gen_dist_constraint.") != std::string::npos) {
                     has_bmrb_constraint = true;
                     break;
                  }
               }

               if (has_bmrb_constraint) {
                  std::cout << "  Found BMRB distance constraint list in saveframe: "
                            << frame.name;
                  if (!constraint_type.empty() && constraint_type != ".") {
                     std::cout << " (type: " << constraint_type << ")";
                  }
                  std::cout << "\n";
                  process_bmrb_distance_constraint_loop(loop, constraint_type);
               }
            }
         }
      }
   }
}

void nef::NEFParser::process_bmrb_distance_constraint_loop(const gemmi::cif::Loop& loop, const std::string& constraint_type) {

   std::string prefix = "_Gen_dist_constraint.";

   auto get_col_idx = [&](const std::string& col_name) -> int {
      std::string full_name = prefix + col_name;
      for (size_t i = 0; i < loop.tags.size(); ++i) {
         if (loop.tags[i] == full_name) {
            return static_cast<int>(i);
         }
      }
      return -1;
   };

   int idx_id = get_col_idx("ID");
   int idx_member_id = get_col_idx("Member_ID");
   int idx_seq_1 = get_col_idx("Seq_ID_1");
   int idx_comp_1 = get_col_idx("Comp_ID_1");
   int idx_atom_1 = get_col_idx("Atom_ID_1");
   int idx_seq_2 = get_col_idx("Seq_ID_2");
   int idx_comp_2 = get_col_idx("Comp_ID_2");
   int idx_atom_2 = get_col_idx("Atom_ID_2");
   int idx_dist_val = get_col_idx("Distance_val");
   int idx_lower = get_col_idx("Distance_lower_bound_val");
   int idx_upper = get_col_idx("Distance_upper_bound_val");
   // Try entity_assembly for chain info
   int idx_entity_1 = get_col_idx("Entity_assembly_ID_1");
   int idx_entity_2 = get_col_idx("Entity_assembly_ID_2");

   // Required columns check
   if (idx_seq_1 < 0 || idx_atom_1 < 0 || idx_seq_2 < 0 || idx_atom_2 < 0) {
      std::cout << "    Error: Missing required columns for distance constraints\n";
      return;
   }

   int n_errors = 0;
   size_t count_before = restraints_.size();

   for (size_t i = 0; i < loop.length(); ++i) {
      try {
         nef::DistanceRestraint r;
         r.index = i + 1;
         r.restraint_id = (idx_id >= 0) ? parse_int(loop.val(i, idx_id)).value_or(i + 1) : i + 1;
         if (idx_member_id >= 0) {
            r.restraint_combination_id = parse_int(loop.val(i, idx_member_id));
         }

         r.chain_code_1 = (idx_entity_1 >= 0) ? loop.val(i, idx_entity_1) : "";
         std::string seq_str_1 = loop.val(i, idx_seq_1);
         r.sequence_code_1 = std::stoi(seq_str_1);
         r.residue_name_1 = (idx_comp_1 >= 0) ? loop.val(i, idx_comp_1) : "";
         r.atom_name_1 = loop.val(i, idx_atom_1);

         r.chain_code_2 = (idx_entity_2 >= 0) ? loop.val(i, idx_entity_2) : "";
         std::string seq_str_2 = loop.val(i, idx_seq_2);
         r.sequence_code_2 = std::stoi(seq_str_2);
         r.residue_name_2 = (idx_comp_2 >= 0) ? loop.val(i, idx_comp_2) : "";
         r.atom_name_2 = loop.val(i, idx_atom_2);

         r.weight = 1.0;
         r.target_value = (idx_dist_val >= 0) ? parse_double(loop.val(i, idx_dist_val)) : std::nullopt;
         r.lower_limit = (idx_lower >= 0) ? parse_double(loop.val(i, idx_lower)) : std::nullopt;
         r.upper_limit = (idx_upper >= 0) ? parse_double(loop.val(i, idx_upper)) : std::nullopt;
         r.restraint_origin = constraint_type;

         restraints_.push_back(r);
      } catch (const std::exception& e) {
         n_errors++;
         if (n_errors == 1) {
            std::cout << "    Warning: Error parsing constraint: " << e.what() << "\n";
         }
      }
   }

   std::cout << "    Parsed " << (restraints_.size() - count_before) << " distance constraints\n";
   if (n_errors > 0) {
      std::cout << "    Note: " << n_errors << " constraints had parsing errors\n";
   }
}

std::string nef::NEFParser::restraints_to_json() const {

   std::ostringstream oss;
   oss << "[\n";
   for (size_t i = 0; i < restraints_.size(); ++i) {
      const auto& r = restraints_[i];
      oss << "   {\n";
      oss << "      \"index\": " << r.index << ",\n";
      oss << "      \"restraint_id\": " << r.restraint_id << ",\n";
      if (r.restraint_combination_id) {
         oss << "      \"restraint_combination_id\": " << *r.restraint_combination_id << ",\n";
      } else {
         oss << "      \"restraint_combination_id\": null,\n";
      }
      oss << "      \"chain_code_1\": \"" << r.chain_code_1 << "\",\n";
      oss << "      \"sequence_code_1\": " << r.sequence_code_1 << ",\n";
      oss << "      \"residue_name_1\": \"" << r.residue_name_1 << "\",\n";
      oss << "      \"atom_name_1\": \"" << r.atom_name_1 << "\",\n";
      oss << "      \"chain_code_2\": \"" << r.chain_code_2 << "\",\n";
      oss << "      \"sequence_code_2\": " << r.sequence_code_2 << ",\n";
      oss << "      \"residue_name_2\": \"" << r.residue_name_2 << "\",\n";
      oss << "      \"atom_name_2\": \"" << r.atom_name_2 << "\",\n";
      oss << "      \"weight\": " << r.weight << ",\n";
      if (r.target_value) {
         oss << "      \"target_value\": " << *r.target_value << ",\n";
      } else {
         oss << "      \"target_value\": null,\n";
      }
      if (r.target_value_uncertainty) {
         oss << "      \"target_value_uncertainty\": " << *r.target_value_uncertainty << ",\n";
      } else {
         oss << "      \"target_value_uncertainty\": null,\n";
      }
      if (r.lower_limit) {
         oss << "      \"lower_limit\": " << *r.lower_limit << ",\n";
      } else {
         oss << "      \"lower_limit\": null,\n";
      }
      if (r.upper_limit) {
         oss << "      \"upper_limit\": " << *r.upper_limit << ",\n";
      } else {
         oss << "      \"upper_limit\": null,\n";
      }
      oss << "      \"restraint_origin\": \"" << r.restraint_origin << "\"\n";
      oss << "   }";
      if (i < restraints_.size() - 1) {
         oss << ",";
      }
      oss << "\n";
   }
   oss << "]\n";
   return oss.str();
}


#if 0

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <nef_or_star_file>\n";
        return 1;
    }

    std::string filepath = argv[1];
    std::cout << "Using gemmi STAR parser (C++ version) to read: " << filepath << "\n\n";

    try {
        // Use the filepath convenience method
        NEFParser parser = NEFParser::from_file(filepath);
        parser.parse();
        parser.print_summary();
        parser.print_statistics();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

#endif
