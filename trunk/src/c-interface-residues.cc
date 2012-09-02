/* src/c-interface-residues.cc
 * 
 * Copyright 2012 by The University of Oxford
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#ifdef USE_GUILE
#include <libguile.h>
#endif

#include <vector>
#include "named-rotamer-score.hh"

#include "cc-interface.hh"
#include "graphics-info.h"

#include "c-interface.h" // for is_valid_model_molecule()


std::vector<coot::named_rotamer_score>
score_rotamers(int imol, 
	       const char *chain_id, 
	       int res_no, 
	       const char *ins_code, 
	       const char *alt_conf, 
	       int imol_map, 
	       int clash_flag,
	       float lowest_probability) {

   std::vector<coot::named_rotamer_score> v;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
	 graphics_info_t g;
	 clipper::Xmap<float> xmap = g.molecules[imol_map].xmap_list[0];
	 v = graphics_info_t::molecules[imol].score_rotamers(chain_id, res_no, ins_code,
							     alt_conf,
							     clash_flag, lowest_probability,
							     xmap, *g.Geom_p());
      }
   }
   return v; 
} 
 

#ifdef USE_GUILE
SCM score_rotamers_scm(int imol, 
		       const char *chain_id, 
		       int res_no, 
		       const char *ins_code, 
		       const char *alt_conf, 
		       int imol_map, 
		       int clash_flag, 
		       float lowest_probability) {

   SCM r = SCM_EOL;
   std::vector<coot::named_rotamer_score> v =
      score_rotamers(imol, chain_id, res_no, ins_code, alt_conf,
		     imol_map, clash_flag, lowest_probability);
   for (unsigned int i=0; i<v.size(); i++) {
      SCM name_scm  = scm_from_locale_string(v[i].name.c_str());
      SCM prob_scm  = scm_double2num(v[i].rotamer_probability_score);
      SCM fit_scm   = scm_double2num(v[i].density_fit_score);
      SCM clash_scm = scm_double2num(v[i].clash_score);
      SCM atom_list_scm = SCM_EOL;
      for (unsigned int iat=0; iat<v[i].density_score_for_atoms.size(); iat++) {
	 SCM p1 = scm_from_locale_string(v[i].density_score_for_atoms[iat].first.c_str());
	 SCM p2 = scm_double2num(v[i].density_score_for_atoms[iat].second);
	 SCM atom_item = scm_list_2(p1,p2);
	 atom_list_scm = scm_cons(atom_item, atom_list_scm);
      }
      atom_list_scm = scm_reverse(atom_list_scm);
      SCM item = scm_list_5(name_scm, prob_scm, atom_list_scm, fit_scm, clash_scm);
      
      r = scm_cons(item, r);
   }
   r = scm_reverse(r);
   return r;
} 
#endif

