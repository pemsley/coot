
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
   sigma_combined_planes = 0.2;
   
   std::vector<std::string> words = coot::util::split_string_no_blanks(line, " ");
   std::vector<std::string> u(words.size());
   for (unsigned int i=0; i<words.size(); i++)
      u[i] = coot::util::upcase(words[i]);

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
				    if (u[7].length() > 3) {
				       if (u[7].substr(0,4) == "CHAI") {
					  try { 
					     plane_1_atoms.res_spec = residue_spec_t(u[8], util::string_to_int(u[6]));
					  } catch (const std::runtime_error &rte) {
					  }
					  if (u[9].length() > 3) {
					     if (u[9].substr(0,4) == "ATOM") {
						if (u[10] == "{") {
						   int o=0;
						   for (unsigned int i=11; i<words.size(); i++) {
						      o++;
						      if (u[i] == "}")
							 break;
						      // std::cout << "adding plan1 atom name " << u[i] << std::endl;
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
   
   if (u[11+o_1].length() > 3) {
      if (u[11+o_1].substr(0,4) == "PLAN") {
	 if (u[12+o_1] == "2") {
	    if (u[13+o_1].length() > 3) {
	       if (u[13+o_1].substr(0,4) == "FIRS") {
	    
		  if (u[14+o_1].length() > 3) {
		     if (u[14+o_1].substr(0,4) == "RESI") {
			
			if (u[16+o_1].length() > 3) {
			   if (u[16+o_1].substr(0,4) == "CHAI") {

			      plane_2_atoms.res_spec = residue_spec_t(u[17+o_1], util::string_to_int(u[15+o_1]));
			      
			      
			      if (u[18+o_1].length() > 3) {
				 if (u[18+o_1].substr(0,4) == "ATOM") {
				    
				    if (u[19+o_1] == "{") {

				       for (unsigned int i=20+o_1; i<u.size(); i++) {
					  o_2++;
					  if (u[i] == "}")
					     break;
					  // std::cout << "adding plan2 atom name " << u[i] << std::endl;
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
   return o_2;
} 


int 
coot::parallel_planes_t::parse_dist_and_type(const std::vector<std::string> &u, int offset) {

   if (u.size() > 10) {
      for (unsigned int i=20+offset; i<u.size(); i++) { 
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

