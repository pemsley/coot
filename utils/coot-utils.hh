/* coot-utils/coot-utils.hh
 *
 * Copyright 2004, 2005, 2006 by The University of York
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef COOT_UTILS_HH
#define COOT_UTILS_HH

#include <string>
#include <vector>
#include <stdlib.h>

#include "colour-holder.hh"

namespace coot {

   std::string get_home_dir();

   // The user can set COOT_DATA_DIR (in fact this is the usual case
   // when using binaries) and that should over-ride the built-in
   // PKGDATADIR.
   //
   // Use this to find things in $prefix/share/coot
   std::string package_data_dir();

   std::string prefix_dir(); // base, i.e. $HOME/coot-build

   // Use this to find things in $prefix/share/RDKit
   std::string rdkit_package_data_dir();

   // use env var COOT_N_THREADS (or fallback) to get the number of threads
   unsigned int get_max_number_of_threads();
   // sets this:
   static unsigned int coot_n_threads;
   // using this:
   long get_number_of_threads_by_system_call();

   // if you can try to get the directoy dir in this directory.
   // if not, try to make it in this directory.
   // if not, try to find it in $HOME
   // if not try to make it in $HOME
   // if not, return the empty string
   std::string get_directory(const std::string &dir);

   // return true for success - Use to copy text files with new lines.
   bool copy_file(const std::string &from_file, const std::string &to_file);

   namespace sequence {

      class fasta {
      public:
	 std::string label;
	 std::string sequence;
	 fasta(const std::string &in_string);
	 bool is_fasta_aa(const std::string &a) const;
      };

      bool is_sequence_triplet(const std::string &s);
   }

   // return empty string on failure
   std::string suggest_new_comp_id(const std::string &comp_id_in);

   std::pair<std::string, std::string> get_userid_name_pair();

   namespace util {

      int round_up_by_hundreds(int num);
      std::string current_working_dir();
      std::string append_dir_dir (const std::string &s1, const std::string &dir);
      std::string append_dir_file(const std::string &s1, const std::string &file);

      // If cwd is a substring of f, then return the basename of f (i.e. cwd
      // stripped from f).  If cwd is not a substring of f, then return f;
      //
      std::string relativise_file_name(const std::string &f, const std::string &cwd);
      std::string absolutise_file_name(const std::string &file_name);
      std::string name_sans_extension(const std::string &f);
      std::string Upper(const std::string &s);
      std::string remove_leading_spaces(const std::string &s);
      std::string remove_whitespace(const std::string &s);
      std::string remove_trailing_whitespace(const std::string &s); // "ALA X  " -> "ALA X";
      std::string int_to_string(int i);
      std::string long_int_to_string(long int i);
      std::string float_to_string(float f);
      std::string float_to_string_using_dec_pl(float f, unsigned short int n_dec_pl);
      std::string float_to_unspaced_string_using_dec_pl(float f, unsigned short int n_dec_pl);
      // throw an runtime_error exception on unable to convert
      int string_to_int(const std::string &s);
      // throw an exception on unable to convert
      float string_to_float(const std::string &s);
      double string_to_double(const std::string &s);
      //
      std::pair<std::string, std::string> split_string_on_last_slash(const std::string &string_in);
      std::vector<std::string> split_string(const std::string &string_in,
					    const std::string &splitter);
      std::vector<std::string> split_string_no_blanks(const std::string &string_in,
						      const std::string &splitter=" ");
      std::vector<std::string> split_string_on_whitespace_no_blanks(const std::string &string_in);

      // can throw a std::runtime_error exception.  If this returns, it guarantees a useful result.
      std::pair<std::string, long> extract_number_string(const std::string &s_in);

      std::string plain_text_to_sequence(const std::string &s);
      std::string plain_text_to_pir(const std::string &title, const std::string &sequence, short int il);
      short int is_fasta_aa(const std::string &s); // single letter code
      std::string single_quote(const std::string &s, const std::string &quote_char="\"");
      // return 0 on success, something else on failure
      int create_directory(const std::string &dir_name);
      // use dds.erase(dds.size() - 1) to strip trailing "/" - see prefix_dir()
      std::string file_name_directory(const std::string &file_name);
      std::string file_name_extension(const std::string &file_name); // including the "."
      std::string file_name_non_directory(const std::string &file_name);
      bool extension_is_for_shelx_coords(const std::string &ext);
      bool extension_is_for_mdl_mol_or_mol2_coords(const std::string &ext);
      bool extension_is_for_coords(const std::string &ext);
      bool extension_is_for_auto_datasets(const std::string &ext);
      bool extension_is_for_scripts(const std::string &ext);
      // void template<T> swap(*T v1, *T v2);

      // is ALA, GLY, TRP, MET, MSE...? (RNA, DNA allowed too)
      bool is_standard_residue_name(const std::string &residue_name);
      // as above but only protein atom names allowed (and MSE).
      bool is_standard_amino_acid_name(const std::string &residue_name);
      // as above but only nucleotide names allowed.
      bool is_standard_nucleotide_name(const std::string &residue_name);

      // return a set of string that match the glob, with the directory name pre-appended
      std::vector<std::string> glob_files(const std::string &dir, const std::string &glob_pattern);

      std::string downcase(const std::string &s);
      std::string upcase(const std::string &s);
      std::string capitalise(const std::string &s); // capitalise first, downcase rest

      std::vector<std::pair<std::string, int> > atomic_number_atom_list();
      int atomic_number(const std::string &atom_name,
			const std::vector<std::pair<std::string, int> > &atom_list);

      // return a long int between 0 and RAND_MAX
      long int random();
      float random_f();
      std::string intelligent_debackslash(const std::string &s);
      std::string remove_trailing_slash(const std::string &s);

      // remove the first bit from long
      std::string remove_string(const std::string &long_string, const std::string &bit);

      int decode_keysym(const std::string &s);
      std::vector<std::pair<std::string, int> > key_sym_vec();

      bool is_number(char c);
      bool is_letter(char c);

      bool even_p(int ii);

      bool close_double_p(const double &d1, const double &d2, const double &diff_crit=0.005);

      inline bool sd_compare(const std::pair<std::string, double> &p1,
			     const std::pair<std::string, double> &p2) {
	 return p1.second < p2.second;
      }


   } // end of util name space

   bool is_member_p(const std::vector<std::string> &v, const std::string &a);
   bool is_member_p(const std::vector<int> &v, const int &a);
   void remove_member(std::vector<int> *v_p, const int &a);

   short int
   is_mmcif_filename(const std::string &filename);

   bool file_exists(const std::string &filename);

   bool file_is_empty(const std::string &filename);

   bool file_exists_and_non_empty(const std::string &file_name);

   bool file_exists_and_non_tiny(const std::string &file_name, unsigned int tiny_size_max = 200);

   bool is_directory_p(const std::string &filename);

   // Gauss Legendre Quadrature

   class gauss_legendre_t {
      void fill_weight_abscicca(int N);
      int N;
      std::vector<std::pair<double, double> > weight_abscissa_;
   public:
      std::pair<double,double> weight_abscissa(int idx) {
	 if (weight_abscissa_.size() == 0)
	    fill_weight_abscicca(N);
	 return weight_abscissa_[idx];
      }
      double weight(int idx) {
	 if (weight_abscissa_.size() == 0)
	    fill_weight_abscicca(N);
	 return weight_abscissa_[idx].first;
      }
      double abscissa(int idx) {
	 if (weight_abscissa_.size() == 0)
	    fill_weight_abscicca(N);
	 return weight_abscissa_[idx].second;
      }
      gauss_legendre_t();
   };



}

#endif // COOT_UTILS_HH
