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

namespace coot {

   // The user can set COOT_DATA_DIR (in fact this is the usual case
   // when using binaries) and that should over-ride the built-in
   // PKGDATADIR.  
   //
   // Use this to find things in $prefix/share/coot
   std::string package_data_dir();

   namespace sequence {

      class fasta {
      public:
	 std::string label;
	 std::string sequence;
	 fasta(const std::string &in_string);
	 bool is_fasta_aa(const std::string &a) const;
      };

   }

   std::pair<std::string, std::string> get_userid_name_pair();

   namespace util { 

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
      // throw an exception on unable to convert
      int string_to_int(const std::string &s);
      // throw an exception on unable to convert
      float string_to_float(const std::string &s);
      // 
      std::pair<std::string, std::string> split_string_on_last_slash(const std::string &string_in);
      std::vector<std::string> split_string(const std::string &string_in,
					    const std::string &splitter);
      std::vector<std::string> split_string_no_blanks(const std::string &string_in,
						      const std::string &splitter);

      std::string plain_text_to_sequence(const std::string &s);
      std::string plain_text_to_pir(const std::string &title, const std::string &sequence, short int il);
      short int is_fasta_aa(const std::string &s); // single letter code
      std::string single_quote(const std::string &s);
      // return 0 on success, something else on failure
      int create_directory(const std::string &dir_name);
      std::string file_name_directory(const std::string &file_name);
      std::string file_name_extension(const std::string &file_name);
      std::string file_name_non_directory(const std::string &file_name);
      bool extension_is_for_shelx_coords(const std::string &ext);
      bool extension_is_for_mdl_mol_or_mol2_coords(const std::string &ext);
      // void template<T> swap(*T v1, *T v2);

      // is ALA, GLY, TRP, MET, MSE...? (RNA, DNA allowed too)
      bool is_standard_residue_name(const std::string &residue_name);

      std::string downcase(const std::string &s);
      std::string upcase(const std::string &s);
      std::string capitalise(const std::string &s); // capitalise first, downcase rest

      std::vector<std::pair<std::string, int> > atomic_number_atom_list();
      int atomic_number(const std::string &atom_name, 
			const std::vector<std::pair<std::string, int> > &atom_list);

      // return a long int between 0 and RAND_MAX
      long int random();
      std::string intelligent_debackslash(const std::string &s);

      // remove the first bit from long
      std::string remove_string(const std::string &long_string, const std::string &bit);

      int decode_keysym(const std::string &s);
      std::vector<std::pair<std::string, int> > key_sym_vec();

      bool is_number(char c);
      bool is_letter(char c);

      bool even_p(int ii);
      

   } // end of util name space
   
   bool is_member_p(const std::vector<std::string> &v, const std::string &a);
   bool is_member_p(const std::vector<int> &v, const int &a);
   void remove_member(std::vector<int> *v_p, const int &a);

   short int
   is_mmcif_filename(const std::string &filename);

   bool file_exists(const std::string &filename);

   bool is_directory_p(const std::string &filename);

   class colour_holder {
   public:
      // values between 0 and 1.0
      float red;
      float green;
      float blue;
      colour_holder(const float &r, const float &g, const float &b) {
	 red = r;
	 green = g;
	 blue = b;
      }
      // needed because it's in a vector.
      colour_holder() { 
	 red = 0.5;
	 green = 0.5;
	 blue = 0.5;
      }
      colour_holder(const std::string &hex_colour_string);
      std::string hex() const;
      friend std::ostream& operator<< (std::ostream& s, const colour_holder &ch);
   };
   std::ostream& operator<< (std::ostream& s, const colour_holder &ch);

   // colour conversion
   std::vector<float> convert_hsv_to_rgb(const std::vector<float> &hsv);
   std::vector<float> convert_rgb_to_hsv(const std::vector<float> &in_vals);
   colour_holder hsv_to_colour(const std::vector<float> &hsv);

}

#endif // COOT_UTILS_HH

