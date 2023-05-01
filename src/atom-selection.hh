/* src/atom-selection.hh
 * 
 * Copyright 2011 by the University of Oxford
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
#ifndef ATOM_SELECTION_HH
#define ATOM_SELECTION_HH

#include <string>
#include <mmdb2/mmdb_manager.h>

namespace coot {
   
   class atom_selection_info_t { 
   public:
      enum { UNSET, BY_STRING, BY_ATTRIBUTES }; 
      int type;
      std::string chain_id;
      int resno_start;
      int resno_end;
      std::string ins_code;
      std::string altconf;
      bool alt_conf_is_set;
      // or:
      std::string atom_selection_str;
      atom_selection_info_t(const std::string &s) { 
	 atom_selection_str = s;
	 type = BY_STRING;
	 alt_conf_is_set = 0;
      }
      atom_selection_info_t(const std::string &chain_id_in, 
			    int resno_start_in, 
			    int resno_end_in,
			    const std::string &ins_code_in) { 
	 chain_id = chain_id_in;
	 resno_start = resno_start_in;
	 resno_end = resno_end_in;
	 ins_code = ins_code_in;
	 type = BY_ATTRIBUTES;
	 alt_conf_is_set = 0;
      }
      atom_selection_info_t(const std::string &chain_id_in, 
			    int resno_start_in, 
			    int resno_end_in,
			    const std::string &ins_code_in,
			    const std::string &alt_conf_in) { 
	 chain_id = chain_id_in;
	 resno_start = resno_start_in;
	 resno_end = resno_end_in;
	 ins_code = ins_code_in;
	 type = BY_ATTRIBUTES;
	 altconf = alt_conf_in;
	 alt_conf_is_set = 1;
      }
      atom_selection_info_t() { 
	 type = UNSET;
	 alt_conf_is_set = 0;
      }

      // Return the selection handle.  It is up to the caller to
      // dispose of the atom selection, with a DeleteSelection().
      int select_atoms(mmdb::Manager *mol) const;
      void using_altconf(const std::string &altconf_in) {
	 altconf = altconf_in;
	 alt_conf_is_set = 1;
      }
      std::string name() const;
      std::string mmdb_string() const;
   };

}


#endif // ATOM_SELECTION_HH

