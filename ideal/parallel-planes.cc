/* ideal/parallel-planes.cc
 * 
 * Copyright 2013, 2015 by Medical Research Council
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

#include "utils/coot-utils.hh"
#include "simple-restraint.hh"
#include "parallel-planes.hh"


// try to set matches (if the line correctly parses as a stacking restraint).
// 
coot::parallel_planes_t::parallel_planes_t(const std::string &line) {

   matches = false;
   target_angle  = 0;
   sigma_angle = 5;
   sigma_distance = 0.2;
   sigma_combined_planes = 0.5;
   
   std::vector<std::string> words = coot::util::split_string_no_blanks(line, " ");
   std::vector<std::string> u(words.size());
   for (unsigned int i=0; i<words.size(); i++)
      u[i] = coot::util::upcase(words[i]);

   if (u.size() > 10) {
      if (u[0].length() > 3) {
	 if (u[0].substr(0,4) == "EXTE") {
	    if (u[1].length() > 3) {
	       if (u[1].substr(0,4) == "STAC") {
		  if (u[2].length() > 3) {
		     if (u[2].substr(0,4) == "PLAN") {
			if (u[3] == "1") {
			   if (u[4].length() > 3) {
			      if (u[4].substr(0,4) == "FIRS") {
				 if (u[5].length() > 3) {
				    if (u[5].substr(0,4) == "RESI") {
				       // u[6] is the residue-number string
				       // u[7] is INS
				       // u[8] is . (usually)
				       // u[9] is CHAI

				       if (u[7] == "INS") {
					  std::string ins_code = "";
					  if (u[8] != ".")
					     ins_code = u[8];
				    
					  if (u[9].length() > 3) {
					     if (u[9].substr(0,4) == "CHAI") {
						try {
						   plane_1_atoms.res_spec = residue_spec_t(u[10], util::string_to_int(u[6]), ins_code);
						} catch (const std::runtime_error &rte) {
						}
						if (u[11].length() > 3) {
						   if (u[11].substr(0,4) == "ATOM") {
						      if (u[12] == "{") {
							 int o=0;
							 for (unsigned int i=13; i<words.size(); i++) {
							    o++;
							    if (u[i] == "}")
							       break;
							    plane_1_atoms.atom_names.push_back(u[i]);
							 }
							 int o_2 = parse_2nd_plane(u, o);
							 if (o_2 > o) {
							    parse_dist_and_type(u, o_2);
							 }
						      }
						   }
						}
					     }
					  }
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (plane_1_atoms.size() > 3) {
      if (plane_2_atoms.size() > 3) {
	 try { 
	    matches = true;
	 }
	 catch (const std::runtime_error &rte) {
	 } 
      }
   }
}



int 
coot::parallel_planes_t::parse_2nd_plane(const std::vector<std::string> &u, int offset) {

   int o_1 = offset;
   int o_2 = o_1;

   if (u[13+o_1].length() > 3) {
      if (u[13+o_1].substr(0,4) == "PLAN") {
	 if (u[14+o_1] == "2") {
	    if (u[15+o_1].length() > 3) {
	       if (u[15+o_1].substr(0,4) == "FIRS") {
		  if (u[16+o_1].length() > 3) {
		     if (u[16+o_1].substr(0,4) == "RESI") {
			
			// u[17+o_1] is residue number string
			// u[18+o_1] is INS
			// u[19+o_1] is .
			// u[20+o_1] is CHAI

			if (u[18+o_1] == "INS") {
			   std::string ins_code = "";
			   if (u[19+o_1] != ".")
			      ins_code = u[19+o_1];
			
			   if (u[20+o_1].length() > 3) {
			      if (u[20+o_1].substr(0,4) == "CHAI") {

				 plane_2_atoms.res_spec = residue_spec_t(u[21+o_1], util::string_to_int(u[17+o_1]), ins_code);
			      
				 if (u[22+o_1].length() > 3) {
				    if (u[22+o_1].substr(0,4) == "ATOM") {
				    
				       if (u[23+o_1] == "{") {

					  for (unsigned int i=24+o_1; i<u.size(); i++) {
					     o_2++;
					     if (u[i] == "}")
						break;
					     plane_2_atoms.atom_names.push_back(u[i]);
					  }
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return o_2;
} 


int 
coot::parallel_planes_t::parse_dist_and_type(const std::vector<std::string> &u, int offset) {

   if (u.size() > 10) {
      for (unsigned int i=24+offset; i<u.size(); i++) {
	 if (u[i].length() > 3) {
	    if (u[i].substr(0,4) == "DIST") {
	       if (u.size() > i+1) {
		  std::string dist_str = u[i+1];
		  try {
		     distance = std::pair<bool, double> (true, util::string_to_float(dist_str));
		  }
		  catch (const std::runtime_error &rte) {
		  } 
		  
	       }
	    }
	 }
      }
   } 
   return offset;

}

std::ostream &
operator<<(std::ostream &s, coot::parallel_planes_t pp) {

   s << "pp-restr: " << pp.plane_1_atoms.res_spec << " " << pp.plane_2_atoms.res_spec;
   return s;

} 

