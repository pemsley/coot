/* src/c-interface-nucleotides.cc
 * 
 * Copyright 2008 The University of Oxford
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

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "mmdb_manager.h"
#include "mmdb-extras.h"

#include "graphics-info.h"

#include "c-interface.h"

#ifdef USE_GUILE
#include <guile/gh.h>
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
// no fix up needed 
#else    
#define scm_to_int gh_scm2int
#define scm_to_locale_string SCM_STRING_CHARS
#define scm_to_double  gh_scm2double
#define scm_is_true gh_scm2bool
#define scm_car SCM_CAR
#endif // SCM version
#endif // USE_GUILE

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
#ifdef USE_PYTHON
#include "Python.h"
#if (PY_MINOR_VERSION > 4) 
// no fixup needed 
#else
#define Py_ssize_t int
#endif
#endif // USE_PYTHON

#include "cc-interface.hh"

#ifdef USE_GUILE
SCM pucker_info_scm(int imol, SCM residue_spec_scm, int do_pukka_pucker_check) {

   std::string altconf = "";
   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_spec_from_scm(residue_spec_scm);
      CResidue *res_p = graphics_info_t::molecules[imol].get_residue(residue_spec);
      if (res_p) {
	 try {
	    coot::pucker_analysis_info_t pi(res_p, altconf);
	    if (do_pukka_pucker_check) {
	       CResidue *following_res =
		  graphics_info_t::molecules[imol].get_following_residue(residue_spec);
	       if (following_res) {
// 		  std::cout << "   DEBUG:: " << coot::residue_spec_t(following_res)
// 			    << " follows " << residue_spec << std::endl;
		  try {
		     double phosphate_d = pi.phosphate_distance(following_res);
		     r = SCM_EOL;
		     r = scm_cons(scm_double2num(pi.plane_distortion), r);
		     r = scm_cons(scm_double2num(pi.out_of_plane_distance), r);
		     r = scm_cons(scm_makfrom0str(pi.puckered_atom().c_str()), r);
		     r = scm_cons(scm_double2num(phosphate_d), r);
		  }
		  catch (std::runtime_error phos_mess) {
		     std::cout << " Failed to find Phosphate for "
			       << coot::residue_spec_t(following_res) << " " 
			       << phos_mess.what() << std::endl;
		  } 
		  
	       } else {
		  r = SCM_EOL;
		  // r = scm_cons(scm_double2num(pi.plane_distortion), r);
		  // r = scm_cons(scm_double2num(pi.out_of_plane_distance), r);
	       } 
	    } else { 
	       r = SCM_EOL;
	       r = scm_cons(scm_double2num(pi.plane_distortion), r);
	       r = scm_cons(scm_double2num(pi.out_of_plane_distance), r);
	       r = scm_cons(scm_makfrom0str(pi.puckered_atom().c_str()), r);
	    }
	 }
	 catch (std::runtime_error mess) {
	    std::cout << " failed to find pucker for " << residue_spec << " " 
		      << mess.what() << std::endl;
	 } 
      } 
   }
   return r;
}
#endif /* USE_GUILE */


