/* src/c-interface-scm.cc
 * 
 * Copyright 2007, 2008, 2009 by The University of Oxford
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


#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#include <string>
#include <vector>

#include "utils/coot-utils.hh"
#include "c-interface-scm.hh"

#include "graphics-info.h"

#ifdef USE_GUILE

#include "guile-fixups.h"


// This is a common denominator really.  It does not depend on mmdb,
// but it can't be declared in c-interface.h because then we'd have to
// include c-interface.h which would cause (resolvable, I think, not
// checked) problems.
// 
// return a scm string, decode to c++ using scm_to_locale_string();
SCM display_scm(SCM o) {

   SCM dest = SCM_BOOL_F;
   SCM mess = scm_from_locale_string("object: ~s");
   return scm_simple_format(dest, mess, scm_list_1(o));
}

bool scm_is_undefined(SCM o) {

   // #t and #f don't have that bit set.  Not sure about others, but
   // undefined does.
   return (1024 & SCM_UNPACK(o)) > 0 ? 1 : 0;
}

// e.g. '("B" 41 "" " CA " "")
std::pair<bool, coot::atom_spec_t>
make_atom_spec(SCM spec) {

   bool good_spec = 0;
   coot::atom_spec_t as;
   SCM spec_length_scm = scm_length(spec);
   int spec_length = scm_to_int(spec_length_scm);

   if (spec_length == 5) {
      SCM  chain_id_scm = scm_list_ref(spec, scm_from_int(0));
      SCM     resno_scm = scm_list_ref(spec, scm_from_int(1));
      SCM  ins_code_scm = scm_list_ref(spec, scm_from_int(2));
      SCM atom_name_scm = scm_list_ref(spec, scm_from_int(3));
      SCM  alt_conf_scm = scm_list_ref(spec, scm_from_int(4));
      std::string chain_id = scm_to_locale_string(chain_id_scm);
      int resno = scm_to_int(resno_scm);
      std::string ins_code  = scm_to_locale_string(ins_code_scm);
      std::string atom_name = scm_to_locale_string(atom_name_scm);
      std::string alt_conf  = scm_to_locale_string(alt_conf_scm);
      as = coot::atom_spec_t(chain_id, resno, ins_code, atom_name, alt_conf);
      good_spec = 1;
   }
   return std::pair<bool, coot::atom_spec_t> (good_spec, as);
}


std::pair<bool, coot::residue_spec_t>
make_residue_spec(SCM spec) {
   bool good_spec = 0;
   coot::residue_spec_t rs("A", 1);
   SCM spec_length_scm = scm_length(spec);
   int spec_length = scm_to_int(spec_length_scm);
   // we can now allow specs that are of length 4.  specs of length
   // are created by het-groups (amongst other things) and have a
   // state in the first position, which we skip (using offset = 1).
   int offset = 0;
   if (spec_length == 4) offset = 1;
   if (spec_length >= 3) {
      SCM chain_id_scm = scm_list_ref(spec, scm_from_int(0+offset));
      SCM resno_scm    = scm_list_ref(spec, scm_from_int(1+offset));
      SCM ins_code_scm = scm_list_ref(spec, scm_from_int(2+offset));
      std::string chain_id = scm_to_locale_string(chain_id_scm);
      int resno = scm_to_int(resno_scm);
      std::string ins_code  = scm_to_locale_string(ins_code_scm);
      rs = coot::residue_spec_t(chain_id, resno, ins_code);
      good_spec = 1;
   }
   
   return std::pair<bool, coot::residue_spec_t>(good_spec, rs);
} 

int key_sym_code_scm(SCM s_scm) {

   int r = -1;
   SCM s_test = scm_string_p(s_scm);
   if (scm_is_true(s_test)) { 
      std::string s = scm_to_locale_string(s_scm);
      r = coot::util::decode_keysym(s);
   }
   return r;
}

// Convert a scheme list of strings to a clipper::Spacegroup.  return
// a Spacegroup that is_null() on failure.
// 
clipper::Spacegroup
scm_symop_strings_to_space_group(SCM symop_string_list) {

   clipper::Spacegroup sg;
   if (scm_is_true(scm_list_p(symop_string_list))) {
      SCM n_scm = scm_length(symop_string_list);
      int n = scm_to_int(n_scm);
      std::string sgo;
      for (int i=0; i<n; i++) {
	 SCM s = scm_list_ref(symop_string_list, scm_from_int(i));
	 std::string se = scm_to_locale_string(s);
	 sgo += se;
	 sgo += " ; ";
      }
      if (sgo.length() > 0) {
	 try {
	    sg.init(clipper::Spgr_descr(sgo, clipper::Spgr_descr::Symops));
	 } catch ( clipper::Message_base exc ) {
	    std::string mess = "Can't make spacegroup from ";
	    mess += sgo;
	    std::cout << "WARNING:: " << mess << std::endl;
	 }
      } 
   } else {
      std::cout << "WARNING:: " << scm_to_locale_string(display_scm(symop_string_list))
		<< " is not a list" << std::endl;
   } 
   return sg;
}

SCM
atom_spec_to_scm(const coot::atom_spec_t &spec) {

   graphics_info_t g;
   return g.atom_spec_to_scm(spec);
}


#endif // USE_GUILE

