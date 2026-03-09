
#ifndef COOT_UTILS_PARSE_GEMMI_NEF
#define COOT_UTILS_PARSE_GEMMI_NEF

#ifdef USE_GEMMI

#include <gemmi/cif.hpp>
#include <gemmi/numb.hpp>

#include <iostream>
#include <string>
#include <optional>
#include <iomanip>

// ============================================================================
// File preprocessing for malformed BMRB files
// ============================================================================

// ============================================================================
// Data structures
// ============================================================================

namespace nef {

   struct AssignedChemicalShift {
       int shift_id;
       int seq_id;
       std::string comp_id;      // Residue type (MET, LYS, etc.)
       std::string atom_id;      // Atom name (HA, HB2, C, CA, N, etc.)
       std::string atom_type;    // Element type (H, C, N, etc.)
       int isotope_number;       // Isotope (1 for 1H, 13 for 13C, 15 for 15N)
       double value;             // Chemical shift value in ppm
       std::optional<double> value_error;  // Uncertainty
       std::optional<int> ambiguity_code;  // Assignment ambiguity

       std::string to_string() const {
           std::ostringstream oss;
           oss << std::setw(3) << seq_id << " "
               << std::setw(3) << comp_id << " "
               << std::setw(4) << atom_id << " ("
               << std::setw(2) << isotope_number << atom_type << "): "
               << std::setw(8) << std::fixed << std::setprecision(3) << value;
           if (value_error) {
               oss << " ±" << *value_error;
           }
           return oss.str();
       }
   };

   struct CouplingConstant {
       int coupling_id;
       int seq_id_1;
       std::string comp_id_1;
       std::string atom_id_1;
       int seq_id_2;
       std::string comp_id_2;
       std::string atom_id_2;
       double value;
       std::optional<double> value_error;

       std::string to_string() const {
           std::ostringstream oss;
           oss << seq_id_1 << " " << comp_id_1 << " " << atom_id_1 << " - "
               << seq_id_2 << " " << comp_id_2 << " " << atom_id_2 << ": "
               << std::fixed << std::setprecision(2) << value << " Hz";
           return oss.str();
       }
   };

   struct NOEPeak {
       int peak_id;
       double volume;
       double height;
       double position_1, position_2, position_3;
       // Atom assignments (chain, seq, res, atom)
       std::tuple<std::string, std::string, std::string, std::string> atom_1;
       std::tuple<std::string, std::string, std::string, std::string> atom_2;
       std::tuple<std::string, std::string, std::string, std::string> atom_3;
       std::string spectrum_name;
   };

   struct DistanceRestraint {
       int index;
       int restraint_id;
       std::optional<int> restraint_combination_id;
       // Atom 1 identification
       std::string chain_code_1;
       int sequence_code_1;
       std::string residue_name_1;
       std::string atom_name_1;
       // Atom 2 identification
       std::string chain_code_2;
       int sequence_code_2;
       std::string residue_name_2;
       std::string atom_name_2;
       // Distance parameters
       double weight;
       std::optional<double> target_value;
       std::optional<double> target_value_uncertainty;
       std::optional<double> lower_limit;
       std::optional<double> upper_limit;
       // Origin (noe, hbond, etc.)
       std::string restraint_origin;

       std::string to_string() const {
           std::ostringstream oss;
           oss << std::setw(4) << restraint_id << ": "
               << chain_code_1 << "/" << std::setw(3) << sequence_code_1 << " "
               << std::setw(3) << residue_name_1 << " " << std::setw(5) << atom_name_1
               << " -- "
               << chain_code_2 << "/" << std::setw(3) << sequence_code_2 << " "
               << std::setw(3) << residue_name_2 << " " << std::setw(5) << atom_name_2;
           if (lower_limit || upper_limit) {
               oss << "  [";
               if (lower_limit) oss << std::fixed << std::setprecision(2) << *lower_limit;
               else oss << "?";
               oss << " - ";
               if (upper_limit) oss << std::fixed << std::setprecision(2) << *upper_limit;
               else oss << "?";
               oss << " A]";
           }
           if (target_value) {
               oss << " target=" << std::fixed << std::setprecision(2) << *target_value << " A";
           }
           return oss.str();
       }
   };

   enum class FileFormat {
       NEF,
       BMRB_STAR,
       UNKNOWN
   };

   std::string preprocess_star_file_string(const std::string& content);
   std::string preprocess_star_file(const std::string& filepath);



// ============================================================================
// NEF/STAR Parser class
// ============================================================================

   class NEFParser {

   public:
      /**
       * Construct parser from file content string (suitable for WebAssembly).
       *
       * @param content The STAR/NEF file content as a string
       * @param source_name A name to identify the source (for error messages)
       */
      NEFParser(const std::string& content, const std::string& source_name);

      /**
       * Construct parser from a filepath (convenience wrapper).
       *
       * @param filepath Path to the STAR/NEF file
       */
      static NEFParser from_file(const std::string& filepath);

   private:

      void init_from_string(const std::string& content, const std::string& source_name);

      // ============================================================================
      // Helper functions
      // ============================================================================

      std::string format_to_string(nef::FileFormat fmt) const;
      bool is_null_value(const std::string& val) const;
      std::optional<double> parse_double(const std::string& val) const;
      std::optional<int> parse_int(const std::string& val) const;

   public:

      void parse();
      void print_summary() const;
      void print_statistics() const;

      // Accessors
      const std::vector<nef::AssignedChemicalShift>& chemical_shifts() const { return chemical_shifts_; }
      const std::vector<nef::CouplingConstant>& coupling_constants() const { return coupling_constants_; }
      const std::vector<nef::NOEPeak>& peaks() const { return peaks_; }
      const std::vector<nef::DistanceRestraint>& restraints() const { return restraints_; }

      // JSON export
      std::string restraints_to_json() const;

   private:
      nef::FileFormat detect_format();
      void parse_assigned_chemical_shifts();
      void process_chemical_shift_loop(const gemmi::cif::Loop& loop);
      void parse_coupling_constants();
      void parse_noe_peaks();
      void parse_distance_restraints();
      void parse_nef_distance_restraints();
      void process_nef_distance_restraint_loop(const gemmi::cif::Loop& loop, const std::string& restraint_origin);
      void parse_bmrb_distance_restraints();
      void process_bmrb_distance_constraint_loop(const gemmi::cif::Loop& loop, const std::string& constraint_type);

      std::string filepath_;
      gemmi::cif::Document doc_;
      gemmi::cif::Block* block_;
      nef::FileFormat format_ = nef::FileFormat::UNKNOWN;

      std::vector<nef::AssignedChemicalShift> chemical_shifts_;
      std::vector<nef::CouplingConstant> coupling_constants_;
      std::vector<nef::NOEPeak> peaks_;
      std::vector<nef::DistanceRestraint> restraints_;
   };

}

#endif // USE_GEMMI
#endif // COOT_UTILS_PARSE_GEMMI_NEF
