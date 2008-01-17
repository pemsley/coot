/* src/c-interface-scm.cc
 * 
 * Copyright 2007 by The University of Oxford
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

#include <string>
#include <vector>

#include "c-interface-scm.hh"

#ifdef USE_GUILE

#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
// no fix up needed 
#else    
#include <guile/gh.h>
#define scm_to_int gh_scm2int
#define scm_to_locale_string SCM_STRING_CHARS
#define scm_to_double  gh_scm2double
#define  scm_is_true gh_scm2bool
#endif // SCM version


// This is a common denominator really.  It does not depend on mmdb,
// but it can't be declared in c-interface.h because then we'd have to
// include c-interface.h which would cause (resolvable, I think, not
// checked) problems.
// 
// return a scm string, decode to c++ using scm_to_locale_string();
SCM display_scm(SCM o) {

   SCM dest = SCM_BOOL_F;
   SCM mess = scm_makfrom0str("object: ~s\n");
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
      SCM  chain_id_scm = scm_list_ref(spec, SCM_MAKINUM(0));
      SCM     resno_scm = scm_list_ref(spec, SCM_MAKINUM(1));
      SCM  ins_code_scm = scm_list_ref(spec, SCM_MAKINUM(2));
      SCM atom_name_scm = scm_list_ref(spec, SCM_MAKINUM(3));
      SCM  alt_conf_scm = scm_list_ref(spec, SCM_MAKINUM(4));
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


#endif // USE_GUILE
