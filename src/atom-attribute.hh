/* src/atom-attributes.cc
 *
 * Copyright 2010 by the University of Oxford
 * Copyright 2015 by Medical Research Council
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


namespace coot {
   class atom_attribute_setting_help_t {
   public:
     enum { UNSET, IS_FLOAT, IS_STRING, IS_INT};
      short int type;
      int i;
      float val;
      std::string s;
      atom_attribute_setting_help_t(const std::string &s_in) {
         s = s_in;
         type = IS_STRING;
         i = -1;
      }
      atom_attribute_setting_help_t(float v) {
         val = v;
         type = IS_FLOAT;
         i = -1;
      }
      atom_attribute_setting_help_t(int iin) {
         i = iin;
         type = IS_INT;
      }
      atom_attribute_setting_help_t() {
         type = UNSET;
         i = -1;
      }
   };

   class atom_attribute_setting_t {
   public:
     atom_spec_t atom_spec;
     std::string attribute_name;
     atom_attribute_setting_help_t attribute_value;
     atom_attribute_setting_t(const std::string &chain_id_in,
			      int resno_in,
			      const std::string &inscode_in,
			      const std::string &atom_name_in,
			      const std::string &alt_conf_in,
			      const std::string &attribute_name_in,
			      const atom_attribute_setting_help_t &att_val) {
       atom_spec = atom_spec_t(chain_id_in, resno_in, inscode_in, atom_name_in, alt_conf_in);
       attribute_name = attribute_name_in;
       attribute_value = att_val;
     }
   };

}
