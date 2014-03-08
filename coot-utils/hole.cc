
#include <fstream>

#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-hole.hh"

coot::hole::hole(CMMDBManager *mol_in, 
		 const clipper::Coord_orth &from_pt_in,
		 const clipper::Coord_orth &to_pt_in,
		 const coot::protein_geometry &geom) {

   mol = mol_in;
   from_pt = from_pt_in;
   to_pt = to_pt_in;
   assign_vdw_radii(geom); 
   colour_map_multiplier = 1.0;
   colour_map_offset = 0.0;
}

void
coot::hole::assign_vdw_radii(const coot::protein_geometry &geom) {

   bool use_vdwH_flag = 0; // extended atoms

   std::map<std::pair<std::string, std::string>, double> cached_radii;
   std::map<std::pair<std::string, std::string>, double>::const_iterator it;

   double radius; 
   
   radius_handle = mol->RegisterUDReal(UDR_ATOM, "atom_radius");
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    std::string residue_name = residue_p->GetResName();
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);

	       std::string atom_name = at->name;

	       // try cache first
	       std::pair<std::string, std::string> p(atom_name, residue_name);

	       if (cached_radii.find(p) != cached_radii.end()) {
		  radius = it->second;
	       } else {
		  radius = geom.get_vdw_radius(atom_name, residue_name, use_vdwH_flag);
	       }
	       if (radius > 0) {
		  at->PutUDData(radius_handle, radius);
	       } else {
		  std::string ele = at->element;
		  // make a reasonable default
		  realtype radius = 1.7;
 		  if (ele == " N")
 		     radius = 1.55;
 		  if (ele == " 0")
 		     radius = 1.52;
 		  if (ele == " H")
 		     radius = 1.2;
		  at->PutUDData(radius_handle, radius); 
	       }
	    }
	 }
      }
   }
}

void
coot::hole::debug_atom_radii() const {

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    std::string residue_name = residue_p->GetResName();
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       realtype radius;
	       at->GetUDData(radius_handle, radius);
	       std::cout << "   " << coot::atom_spec_t(at) << " with radius " << radius<< std::endl;
	    }
	 }
      }
   }
} 

void 
coot::hole::make_atom_selection(int selhnd, const clipper::Coord_orth &pt,
				double radius_prev) const {

   // what protein atoms are within radius_prev of pt?

   mol->SelectAtoms(selhnd, 0, "*",
		    ANY_RES, "*", ANY_RES, "*",
		    "!HOH", // residue names
		    "*", // atom names
		    "H,C,O,N,S", // elements
		    "*", "*", "*", // altlocs, segments, charges
		    -1, -1, // occupancy (no limits)
		    pt.x(), pt.y(), pt.z(), radius_prev, SKEY_NEW);

   if (0) {
      PPCAtom atom_selection = 0;
      int n_selected_atoms;
      mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
      std::cout << "in make_atom_selection() selected " << n_selected_atoms
		<< " atoms in sphere selection around  " << pt.format() << std::endl;
   }
}


std::pair<std::vector<std::pair<clipper::Coord_orth, double> >, std::vector<coot::hole_surface_point_t> >
coot::hole::generate() {

   int nsteps = 200;
   double radius_prev = 7; // max prob size + 1.7 (max radius) + a bit.

   // position and sphere radius
   // 
   std::vector<std::pair<clipper::Coord_orth, double> > probe_path;
   
   clipper::Coord_orth diff = to_pt - from_pt;
   v_hat = clipper::Coord_orth(diff.unit());

   clipper::Coord_orth prev_point = from_pt; // initally.
   clipper::Coord_orth frag = 1/double(nsteps) * diff;
   
   for (unsigned int istep=0; istep<=nsteps; istep++) {
      clipper::Coord_orth l(double(istep) * frag);
      clipper::Coord_orth pt_linear(from_pt + l);

      // We are considering where to place the start point for the new
      // level (pt).  It should be
      // 1) just down a bit (diff) from the previous level
      // and/or
      // 2) a step along the vector from_pt -> to_pt.
      //
      // If/when the point previous new_pt has wandered off into the
      // unknown (because it was so high that it managed to escape
      // being in the channel) then what do we do? - because using 1)
      // from there would mean that we would miss then channel
      // completely.
      //
      // But following the 2) rule would mean that we would not be
      // able to follow kinky paths and that the precision of the
      // (start and) finish points would be important.
      //
      // So let's choose rule 1 except where that means we have gone
      // way off (5A?)
      // 
      clipper::Coord_orth pt = prev_point + frag;
      double dist_diff = clipper::Coord_orth::length(pt_linear, pt);

      if (0)
	 std::cout << "dist_diff: " << dist_diff << " " << pt_linear.format()
		   << " " << pt.format() << " frag: " << frag.format() << " l:"
		   << l.format() << " i_step: " << istep << std::endl;

      if (dist_diff > 7)
	 pt = pt_linear;
      

      int selhnd = mol->NewSelection();
      make_atom_selection(selhnd, pt, radius_prev);
      std::pair<clipper::Coord_orth, double> new_pt = optimize_point(pt, selhnd);
      probe_path.push_back(new_pt);
      double d = sqrt(l.lengthsq());
      if (0) 
	 std::cout << "istep: " << istep << " l: " << d
		   << " ss: " << new_pt.second << "        "
		   << new_pt.first.x() << " "
		   << new_pt.first.y() << " "
		   << new_pt.first.z() << " "
		   << std::endl;
      mol->DeleteSelection(selhnd);
      // next round
      prev_point = new_pt.first;
   }

   std::vector<coot::hole_surface_point_t> surface_points = get_surface_points(probe_path);
   // write_probe_path_using_spheres(surface_points, "probe-surface-points");

   return std::pair<std::vector<std::pair<clipper::Coord_orth, double> >, std::vector<coot::hole_surface_point_t> > (probe_path, surface_points);
}

// By moving about in the plane we try to maximise the distance
// between the point and the atoms around it.
// 
// Return the position and the max dist.
// 
std::pair<clipper::Coord_orth, double>
coot::hole::optimize_point(const clipper::Coord_orth &pt, int selhnd) {

   int max_moves_without_improvement = 80;
   int move_count = 0; // increment if we try to make a move that does not improve
   double biggest_allowed_distance_to_protein = 5.0;
   
   PPCAtom atom_selection = 0;
   int n_selected_atoms;
   mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
   // std::cout << "optimize_point(): n_selected_atoms: " << n_selected_atoms
   // << std::endl;
   
   clipper::Coord_orth current_pt = pt;

   double D_max = 0.1;
   float rmi = 1.0/float(RAND_MAX);

   double current_ss = sphere_size(pt, selhnd);
   while (move_count < max_moves_without_improvement) {
      double d1 = 2.0 * coot::util::random() * rmi -1.0;
      double d2 = 2.0 * coot::util::random() * rmi -1.0;
      double d3 = 2.0 * coot::util::random() * rmi -1.0;
      clipper::Coord_orth y_rand(d1, d2, d3);
      clipper::Coord_orth y_hat_rand(y_rand.unit());
      clipper::Coord_orth y_prime_rand =
	 y_hat_rand - (clipper::Coord_orth::dot(v_hat, y_hat_rand)) * v_hat;
      
      clipper::Coord_orth trial_pt = current_pt + D_max * y_prime_rand;

      // sphere_size() can return negative (no atoms selected).
      // 
      double ss = sphere_size(trial_pt, selhnd);
      if (ss > current_ss) {
	    current_pt = trial_pt;
	    current_ss = ss;
// 	    std::cout << "   updated current_ss to " << current_ss
// 		      << " for start point " << pt.format() << " and trial point "
// 		      << trial_pt.format() << "\n";
	    move_count = 0; // reset
	    if (ss > biggest_allowed_distance_to_protein) {
	       // std::cout << "breaking..." << std::endl;
	       break; // that's enough, the point has escaped presumably
	    }
      } else {
	 move_count++;
      }
   }
   return std::pair<clipper::Coord_orth, double> (current_pt, current_ss);
}


double
coot::hole::sphere_size(const clipper::Coord_orth &pt, int selhnd) const {

   double r = -1;
   PPCAtom atom_selection = 0;
   int n_selected_atoms;
   mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

   // std::cout << "sphere_size: n_selected_atoms: " << n_selected_atoms
   // << std::endl;
   
   double largest_possible_sphere = 99999;
   bool was_set = 0;
   realtype atom_vdw_radius;

   for (unsigned int iat=0; iat<n_selected_atoms; iat++) {
      clipper::Coord_orth atom_pos(atom_selection[iat]->x,
				   atom_selection[iat]->y,
				   atom_selection[iat]->z);
      double r_1 = clipper::Coord_orth::length(atom_pos, pt);
      atom_selection[iat]->GetUDData(radius_handle, atom_vdw_radius);
      double r = r_1 - atom_vdw_radius;
      if (r < largest_possible_sphere) {
	 largest_possible_sphere = r;
	 if (! was_set)
	    was_set = 1;
      }
   }
   
   if (was_set)
      return largest_possible_sphere;
   else
      return -1; // fail
}

// point on the wall of the probe path.
// 
void
coot::hole::write_probe_path(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path) const {

   std::string file_name = "probe-surface-points";

   unsigned int n_theta_points = 20;
   double inv_n_theta = 2 * M_PI/double(n_theta_points);

   std::ofstream render_stream(file_name.c_str());

   if (!render_stream) {
      std::cout << "failed to open " << file_name << std::endl;
   } else { 
  
      clipper::Coord_orth arb(1.1, 1.2, 1.3);
      clipper::Coord_orth plane_vect(clipper::Coord_orth::cross(v_hat, arb));
      clipper::Coord_orth unit_plane_vect(plane_vect.unit());

      for (unsigned int i=0; i<probe_path.size(); i++) {

	 unsigned int max_upstream_check_limit = i-5;
	 unsigned int max_downstream_check_limit = i+5;

	 if (max_upstream_check_limit<0)
	    max_upstream_check_limit = 0;
	 if (max_downstream_check_limit>=probe_path.size())
	     max_downstream_check_limit=probe_path.size()-1;

	 render_stream << probe_path[i].first.x() << " "
		       << probe_path[i].first.y() << " "
		       << probe_path[i].first.z() << " \"red\"\n";
	 
	 std::string colour = "blue";

	 if (probe_path[i].second < 2.1)
	    colour = "blue";
	 if (probe_path[i].second < 1.9)
	    colour = "cyan";
	 if (probe_path[i].second < 1.7)
	    colour = "green";
	 if (probe_path[i].second < 1.5)
	    colour = "greentint";
	 if (probe_path[i].second < 1.3)
	    colour = "sea";
	 if (probe_path[i].second < 1.1)
	    colour = "yellow";
	 if (probe_path[i].second < 0.9)
	    colour = "yelllowtint";
	 if (probe_path[i].second < 0.7)
	    colour = "orange";
	 if (probe_path[i].second < 0.5)
	    colour = "red";
	 if (probe_path[i].second < 0.3)
	    colour = "hotpink";
	 

	 for (unsigned int itheta=0; itheta<n_theta_points; itheta++) {
	    double theta = inv_n_theta * double(itheta);
	    if (coot::util::even_p(i))
	       theta += inv_n_theta * 0.5;

	    // now let's rotate unit_plane_vect around the circle.
	    clipper::Coord_orth pos(probe_path[i].second * unit_plane_vect);
	    clipper::Coord_orth circle_point =
	       coot::util::rotate_round_vector(v_hat,
					       pos,
					       clipper::Coord_orth(0,0,0),
					       theta);
	    clipper::Coord_orth surface_point = probe_path[i].first + circle_point;

	    render_stream << surface_point.x() << " "
			  << surface_point.y() << " "
			  << surface_point.z() << " "
			  << "\"" << colour << "\""
			  << "\n";
	 }
      }
   }

}

void
coot::hole::write_probe_path_using_spheres(const std::vector<coot::hole_surface_point_t> &surface_points,
					   const std::string &file_name) const {

   std::ofstream render_stream(file_name.c_str());
   if (!render_stream) {
      std::cout << "failed to open " << file_name << std::endl;
   } else {
      for (unsigned int i=0; i<surface_points.size(); i++) { 
	 render_stream << surface_points[i].position.format() << " "
		       << surface_points[i].normal.format() << " "
		       << surface_points[i].colour << "\n";
      }
   } 
}


std::vector<coot::hole_surface_point_t> 
coot::hole::get_surface_points(const std::vector<std::pair<clipper::Coord_orth, double> > &probe_path) const {

   std::vector<coot::hole_surface_point_t> surface_points; // return this.
   
   double dot_density = 0.35;
   // 
   double phi_step = 5.0 * (M_PI/180.0);
   double theta_step = 5.0 * (M_PI/180.0);

   for (unsigned int i=0; i<probe_path.size(); i++) {

      // int because we max_upstream_check_limit might go negative (briefly).
      // 
      int max_upstream_check_limit = i-10;
      int max_downstream_check_limit = i+10;
      const clipper::Coord_orth &this_centre = probe_path[i].first;
      const double &r = probe_path[i].second;

      if (max_upstream_check_limit<0)
	 max_upstream_check_limit = 0;
      if (max_downstream_check_limit>=probe_path.size())
	 max_downstream_check_limit=probe_path.size()-1;

      coot::colour_holder vertex_colour = probe_size_to_colour(probe_path[i].second);
      // std::cout << "         " << vertex_colour << std::endl;
      bool even = 1;

      // the begin and end point we don't want to draw the tops and bottoms of the "bulbs".
      // Leave them open.
      //
      bool check_ends = 0;
      double dot_prod_multiplier = 1;
      if (i == 0) {
	 check_ends = 1;
      }
      if (i == (probe_path.size() -1)) {
	 check_ends = 1;
	 dot_prod_multiplier = -1;
      }
	    
      for (double theta=0; theta<M_PI; theta+=theta_step) {
	 double phi_step_inner = phi_step + 0.1 * pow(theta-0.5*M_PI, 2);
	 double sin_theta = sin(theta);
	 double cos_theta = cos(theta);
	 for (double phi=0; phi<2*M_PI; phi+=phi_step_inner) {
	    if (even) { // even
	       clipper::Coord_orth origin_based_sphere_point(r*cos(phi)*sin_theta,
							     r*sin(phi)*sin_theta,
							     r*cos_theta);
	       clipper::Coord_orth surface_point = this_centre + origin_based_sphere_point;
	       clipper::Coord_orth normal(origin_based_sphere_point.unit());
	       bool reject = 0;
	       for (unsigned int iprev=max_upstream_check_limit; iprev<=max_downstream_check_limit; iprev++) {
		  if (iprev != i) {

		     clipper::Coord_orth diff = probe_path[iprev].first - surface_point;
		     double d_sqrd = diff.lengthsq();
		     double probe_prev_sqrd = probe_path[iprev].second * probe_path[iprev].second;
		     if (d_sqrd < probe_prev_sqrd) {
			reject = 1;
			break;
		     }
		  }
	       }
	       if (! reject) {

		  bool do_it = 1;
		     
		  // turn off do_it so that we have open cups at the end points.
		  // 
		  if (check_ends == 1) {
		     double dp = clipper::Coord_orth::dot(v_hat, origin_based_sphere_point);
		     if (dp * dot_prod_multiplier < 0)
			do_it = 0;
		  }
		     
		  if (do_it) { 
		     coot::hole_surface_point_t sp(surface_point, normal, vertex_colour);
		     surface_points.push_back(sp);
		  }
	       }
	    }
	    // next round
	    even = 1 - even;
	 }
      }
   }
   return surface_points;
} 
