
#include <iostream>
#include <string>

#include "coot-utils/coot-coord-utils.hh"
#include "quaternion.hh"

#include "canyon.hh"

coot::canyon::canyon(mmdb::Manager *mol_in) {
   mol = mol_in;
   selHnd = -1;
   // unset values for surface points
   for (unsigned int i=0; i<N_CANYON_STEPS; i++) { 
      for (unsigned int j=0; j<N_THETA_STEPS; j++) { 
	 surface_points[i][j].first = false;
      }
   }
}

void
coot::canyon::update_atom_selection(const clipper::Coord_orth &pt) {

   mmdb::realtype radius = 12;
   if (mol) {
      mol->DeleteSelection(selHnd);
      selHnd = mol->NewSelection();
      mol->SelectSphere(selHnd, mmdb::STYPE_ATOM, pt.x(), pt.y(), pt.z(), radius);
   }
}

void
coot::canyon::trace() {

   std::cout << "------------- calling get_initial_trace() " << std::endl;
   trace_info_t ti = get_initial_trace();

   std::cout << "------------- calling get_refined_trace() " << std::endl;
   ti = get_refined_trace(ti);
}

coot::trace_info_t
coot::canyon::get_initial_trace() {

   clipper::Coord_orth start_diff = start_path_point - start_point;
   clipper::Coord_orth start_path_uv(start_diff.unit());

   // end point and its diff go in the other direction to the start
   clipper::Coord_orth end_diff = end_path_point - end_point;
   clipper::Coord_orth end_path_uv(-end_diff.unit());
   

   clipper::Coord_orth current_point = start_point;
   double d_tot = sqrt((end_point - start_point).lengthsq());
   double path_step = d_tot/double(N_CANYON_STEPS);

   std::vector<clipper::Coord_orth> probe_path_points;

   for (unsigned int i=0; i<N_CANYON_STEPS; i++) {

      // every 5 time
      if (i%5 == 0) {
	 update_atom_selection(current_point);
      }

      // Adding i as part of the function that calculates pt, is not good because
      // current_point now moves to the centre of the canyon points.
      // 
      // clipper::Coord_orth pt(current_point + i * path_uv * path_step);

      double f = double(i)/double(N_CANYON_STEPS);
      clipper::Coord_orth path_uv = get_path_uv(start_path_uv, end_path_uv, f);
      
      clipper::Coord_orth pt(current_point + path_uv * path_step);
      
      trace_info_t dummy_trace_info;
      current_point = add_probe_points(i, f, pt, dummy_trace_info, false);

      probe_path_points.push_back(current_point);
   }

   // output_grid();
   
   clipper::Array2d<double> beta_vec_hats = polynomial_path_fit(probe_path_points);

   std::string file_name = "canyon-path-fit-initial.tab";
   write_fit_path(probe_path_points, beta_vec_hats, file_name);

   trace_info_t trace_info(probe_path_points, beta_vec_hats);
   return trace_info;
}

coot::trace_info_t
coot::canyon::get_refined_trace(const coot::trace_info_t &ti) {

   clipper::Coord_orth start_diff = start_path_point - start_point;
   clipper::Coord_orth start_path_uv(start_diff.unit());

   // end point and its diff go in the other direction to the start
   clipper::Coord_orth end_diff = end_path_point - end_point;
   clipper::Coord_orth end_path_uv(-end_diff.unit());
   

   clipper::Coord_orth current_point = start_point;
   double d_tot = sqrt((end_point - start_point).lengthsq());
   double path_step = d_tot/double(N_CANYON_STEPS);

   std::vector<clipper::Coord_orth> probe_path_points;

   for (unsigned int i=0; i<N_CANYON_STEPS; i++) {

      // every 5 time
      if (i%5 == 0) {
	 update_atom_selection(current_point);
      }

      // Adding i as part of the function that calculates pt, is not good because
      // current_point now moves to the centre of the canyon points.
      // 
      // clipper::Coord_orth pt(current_point + i * path_uv * path_step);

      double f = double(i)/double(N_CANYON_STEPS);
      clipper::Coord_orth path_uv = get_path_uv(ti, f);
      
      clipper::Coord_orth pt(current_point + path_uv * path_step);

      
      if (1) {   // debug uv using model path
	 clipper::Coord_orth tmp_pt(current_point + 2 * path_uv); 
	 std::cout << "line-uv "
		   << current_point.x() << " " 
		   << current_point.y() << " " 
		   << current_point.z() << "     "
		   << tmp_pt.x() << " " 
		   << tmp_pt.y() << " " 
		   << tmp_pt.z() << " "
		   << "#aabb22" << std::endl;
      }

      current_point = add_probe_points(i, f, pt, ti, true);

      probe_path_points.push_back(current_point);
   }

   output_grid();
   
   clipper::Array2d<double> beta_vec_hats = polynomial_path_fit(probe_path_points);

   std::string file_name = "canyon-path-fit-refined.tab";
   write_fit_path(probe_path_points, beta_vec_hats, file_name);

   trace_info_t trace_info(probe_path_points, beta_vec_hats);
   return trace_info;
}


void 
coot::canyon::write_fit_path(const std::vector<clipper::Coord_orth> &path_points,
			     const clipper::Array2d<double> &beta,
			     const std::string &file_name) const {

   std::ofstream f(file_name.c_str());
   for (unsigned int i=0; i<path_points.size(); i++) {
      double t = double(i)/double(path_points.size());
      f << path_points[i].x() << " "
	<< path_points[i].y() << " "
	<< path_points[i].z() << "     ";
      f << beta(0,0) + beta(0,1) * t + beta(0,2) * t * t << " "
	<< beta(1,0) + beta(1,1) * t + beta(1,2) * t * t << " "
	<< beta(2,0) + beta(2,1) * t + beta(2,2) * t * t << "\n";
   }
   f.close();
} 

clipper::Coord_orth
coot::canyon::get_path_uv(const trace_info_t &ti, double t) const {


   // we differentiate the beta vec hats
   // dx = \beta_1 + 2 \beta_2 * t
   // 
   clipper::Coord_orth r(ti.beta_vec_hats(0,1) + ti.beta_vec_hats(0,2) * 2 * t,
			 ti.beta_vec_hats(1,1) + ti.beta_vec_hats(1,2) * 2 * t,
			 ti.beta_vec_hats(2,1) + ti.beta_vec_hats(2,2) * 2 * t); 
   return clipper::Coord_orth(r.unit());
} 



clipper::Coord_orth
coot::canyon::get_path_uv(double frac) const {

   clipper::Coord_orth start_diff = start_path_point - start_point;
   clipper::Coord_orth start_path_uv(start_diff.unit());

   // end point and its diff go in the other direction to the start
   clipper::Coord_orth end_diff = end_path_point - end_point;
   clipper::Coord_orth end_path_uv(-end_diff.unit());
   clipper::Coord_orth r = get_path_uv(start_path_uv, end_path_uv, frac);

   return r;
}

clipper::Coord_orth
coot::canyon::get_path_uv(const clipper::Coord_orth &start_path_uv,
			  const clipper::Coord_orth &end_path_uv,
			  double frac) const {

   clipper::Coord_orth r = start_path_uv;
   quaternion q1(start_path_uv);
   quaternion q2(  end_path_uv);
   quaternion s = quaternion::slerp(q1, q2, frac);

   if (false) { 
      std::cout << "q1: " << q1 << " with 3-part length: " << q1.len3() << std::endl;
      std::cout << "q2: " << q2 << " with 3-part length: " << q2.len3() << std::endl;
      std::cout << " s: " << s  << " with 3-part length: " <<  s.len3() << std::endl;
   }

   r = s.first3();
   
   return r;
} 

// We need a signal that the get_path_uv() function should use the trace_info
// 
clipper::Coord_orth
coot::canyon::add_probe_points(int iround, double frac, clipper::Coord_orth &pt,
			       const trace_info_t &ti,
			       bool use_trace_info_for_path_uv) {

   clipper::Coord_orth tracked_centre = pt;  // initially

   if (true) 
      std::cout << "centre point " << "null" << " path-length " << "null" << " "
		<< pt.x() << " " 
		<< pt.y() << " " 
		<< pt.z() << " "
		<< "#804040" << "\n";   
   
   mmdb::PPAtom atom_selection = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(selHnd, atom_selection, n_selected_atoms);
   std::vector<clipper::Coord_orth> pt_interesting;

   // double delta_theta = 2 * M_PI / double(n_theta_values); // radians

   double path_max = 7; // how far from the current centre should we look for contacts?
                        // Maybe this should be a user-settable variable.

   clipper::Coord_orth arb(0.1, 0.2, 0.3);

   // path_uv needs to get updated - so that we track along the canyon?
   // 
   // clipper::Coord_orth diff = start_path_point - start_point;
   // clipper::Coord_orth path_uv(diff.unit());

   clipper::Coord_orth path_uv = get_path_uv(frac);
   if (use_trace_info_for_path_uv)
      path_uv = get_path_uv(ti, frac);
   std::cout << "path_uv " << path_uv.format() << " " << frac << std::endl;

   clipper::Coord_orth cpu(clipper::Coord_orth::cross(path_uv, arb).unit());

   double probe_radius = 1.4;
   
   double theta_max = 2*M_PI;
   for (unsigned int i_theta=0; i_theta<N_THETA_STEPS; i_theta++) {

      double t = double(i_theta) * 2 * M_PI / double(N_THETA_STEPS);

      // now rotate cpu by theta t around diff vector
      clipper::Coord_orth origin(0,0,0);
      //                                                      dir pos orig-shift angle
      clipper::Coord_orth rpu(coot::util::rotate_round_vector(path_uv, cpu, origin, t).unit());

      for (double ps=0; ps<path_max; ps += 0.2) {
	 clipper::Coord_orth test_pos(pt + rpu * ps);

	 bool bumped = false;
	 clipper::Coord_orth contact_pt(0, 0, 0); // gets set when bumped gets set.

	 int idx_bumped_atom = -1;
	 double bumped_atom_radius = 1.4; // gets set when idx_bumped_atom gets set
	 
	 for (unsigned int iat=0; iat<n_selected_atoms; iat++) { 
	    double atom_radius = 1.7; // look this up per atom
	    double d_crit_sq = (probe_radius + atom_radius) * (probe_radius + atom_radius);
	    double this_dist_sq =
	       ((test_pos.x() - atom_selection[iat]->x) *  (test_pos.x() - atom_selection[iat]->x) + 
		(test_pos.y() - atom_selection[iat]->y) *  (test_pos.y() - atom_selection[iat]->y) +
		(test_pos.z() - atom_selection[iat]->z) *  (test_pos.z() - atom_selection[iat]->z));
	    if (this_dist_sq < d_crit_sq) {
	       bumped = true;
	       if (idx_bumped_atom == -1) { 
		  idx_bumped_atom = iat;
		  bumped_atom_radius = atom_radius;
	       } else {
		  // is this atom closer?
		  mmdb::Atom *current_bump_atom = atom_selection[iat];
		  double d_sq =
		     (current_bump_atom->x - test_pos.x()) + (current_bump_atom->x - test_pos.x()) + 
		     (current_bump_atom->y - test_pos.y()) + (current_bump_atom->y - test_pos.y()) + 
		     (current_bump_atom->z - test_pos.z()) + (current_bump_atom->z - test_pos.z());
		  if (d_sq < this_dist_sq) {
		     idx_bumped_atom = iat;
		     bumped_atom_radius = atom_radius;
		  } 
	       }
	       // break;
	    }
	 }

	 if (bumped) {

	    double fm = probe_radius/(probe_radius + bumped_atom_radius);
	    double om_fm = 1.0 - fm;
	    mmdb::Atom *bump_atom = atom_selection[idx_bumped_atom];
	    clipper::Coord_orth av_pos((fm * bump_atom->x + om_fm * test_pos.x()),
				       (fm * bump_atom->y + om_fm * test_pos.y()),
				       (fm * bump_atom->z + om_fm * test_pos.z()));
	    
	    pt_interesting.push_back(test_pos);
 	    surface_points[iround][i_theta].first = true;
 	    surface_points[iround][i_theta].second = test_pos;
	    break;
	 }
      }
   }

   if (pt_interesting.size() > 4) {
      // get the average of pt_interesting and assign that to tracked_centre.

      // save this for debugging output
      clipper::Coord_orth inital_tracked_centre = tracked_centre;
      tracked_centre = coot::util::average_position(pt_interesting);
      std::cout << "updated tracked_centre from " << inital_tracked_centre.format()
		<< " to " << tracked_centre.format() << std::endl;

   }
   return tracked_centre;
}

// points_vec_index is either 0, 1 or 2 to estimate the x,y or z
// coordinates from the points.
//
// Here we use a convenient 3x3 container:
// 
clipper::Array2d<double>
coot::canyon::polynomial_path_fit(const std::vector<clipper::Coord_orth> &points_in) const {

   clipper::Array2d<double> r(3,3);

   std::vector<double> e_x = polynomial_path_fit(points_in, 0);
   std::vector<double> e_y = polynomial_path_fit(points_in, 1);
   std::vector<double> e_z = polynomial_path_fit(points_in, 2);

   r(0,0) = e_x[0]; r(0,1) = e_x[1]; r(0,2) = e_x[2];
   r(1,0) = e_y[0]; r(1,1) = e_y[1]; r(1,2) = e_y[2];
   r(2,0) = e_z[0]; r(2,1) = e_z[1]; r(2,2) = e_z[2];
   return r;
} 

// points_vec_index is either 0, 1 or 2 to estimate the x,y or z
// coordinates from the points.
// 
std::vector<double> 
coot::canyon::polynomial_path_fit(const std::vector<clipper::Coord_orth> &points_in,
				  unsigned int points_vec_index) const {

   std::vector<clipper::Coord_orth> points = points_in;

   // maybe beta \toparrow is \vec{\beta}
   // so maybe
   // \hat{\vec{beta}} is what we want (and return).

   // x_i = \beta_0 + \beta_1 t + \beta_2 t^2

   // beta \toparrow \hat is the vector of parameters that we wish to
   // evaluate

   // For n points:
   // 
   // X is
   //
   // (1 t_1 t_1^2)
   // (1 t_2 t_2^2)
   // (1 t_3 t_3^2)
   // (.. ..      )
   // (1 t_n t_n^2)
   //
   // \hat{\vec{beta}} is a 3 element vertical vector
   // 
   // (\beta_0)
   // (\beta_1)
   // (\beta_2)
   //
   // \hat{\vec{\beta}} = (X^TX)^{-1}X^T\vec{x}

   // First, let's make X.
   //
   t_matrix X(points.size(), 3);
   double inv_n = 1.0/double(points.size());
   for (unsigned int ii=0; ii<points.size(); ii++) { 
      double f = double(ii) * inv_n;
      X(ii,0) = 1.0;
      X(ii,1) = f;
      X(ii,2) = f*f;
   }

   clipper::Mat33<double> XTX = X.XTX();
   // std::cout << "XTX\n" << XTX.format() << std::endl;
   clipper::Mat33<double> XTX_inv = XTX.inverse();
   // std::cout << "XTX_inv\n" << XTX_inv.format() << std::endl;
   
   // multiply the x vector (1 column, many rows) by X^T (transpose)
   // (which has many columns, 3 rows).
   // 
   std::vector<double> XTx(3);
   for (unsigned int i=0; i<3; i++) {
      double sum = 0;
      for (unsigned int ii=0; ii<points.size(); ii++) {
	 // get the only x coord of the points initially
	 sum += points[ii][points_vec_index] * X(ii,i);
      }
      XTx[i] = sum;
   }

   if (0) { 
      std::cout << "XTx\n";
      for (unsigned int i=0; i<3; i++)
	 std::cout << "  " << XTx[i] << std::endl;
   }
   

   // now multiply XTx by XTX_inv
   //
   std::vector<double> beta_vec_hat(3);  // degree 3 polynomial at the moment
   
   for (unsigned int i=0; i<3; i++) {
      double sum = 0;
      for (unsigned int j=0; j<3; j++) {
	 sum += XTX_inv(i,j) * XTx[j];
      }
      beta_vec_hat[i] = sum;
   }

   std::cout << "beta vec hat" << std::endl;
   for (unsigned int i=0; i<3; i++) {
      std::cout << "    " << i << " " << beta_vec_hat[i] << std::endl;
   }

   return beta_vec_hat;
}



void
coot::canyon::output_grid() const {

   bool show_grid = true;
   double d_max_sq = 3.0; // don't draw triangles with lines that are longer than this.

   std::ofstream f("canyon-grid-lines.data");
   
   for (unsigned int i=0; i<(N_CANYON_STEPS-1); i++) { 
      for (unsigned int j=0; j<N_THETA_STEPS; j++) {
	 unsigned int next_theta = j + 1;
	 if (next_theta == N_THETA_STEPS)
	    next_theta = 0;

	 // easy case: all grid points exist:
	 if (surface_points[i][j].first) {
	    if (surface_points[i+1][j].first) {
	       if (surface_points[i][next_theta].first) {
		  if (surface_points[i+1][next_theta].first) {

		     if (clipper::Coord_orth(surface_points[i][j].second - surface_points[i][next_theta].second).lengthsq() < d_max_sq)
			f << "line "
			  << surface_points[i][j].second.x() << " " 
			  << surface_points[i][j].second.y() << " " 
			  << surface_points[i][j].second.z() << " " 
			  << surface_points[i][next_theta].second.x() << " " 
			  << surface_points[i][next_theta].second.y() << " " 
			  << surface_points[i][next_theta].second.z() << " " 
			  << "sea\n";

		     if (show_grid) { 

			if (clipper::Coord_orth(surface_points[i][j].second - surface_points[i+1][j].second).lengthsq() < d_max_sq)
			   f << "line "
			     << surface_points[i][j].second.x() << " " 
			     << surface_points[i][j].second.y() << " " 
			     << surface_points[i][j].second.z() << " " 
			     << surface_points[i+1][j].second.x() << " " 
			     << surface_points[i+1][j].second.y() << " " 
			     << surface_points[i+1][j].second.z() << " " 
			     << "sea\n";
		     
			if (clipper::Coord_orth(surface_points[i+1][j].second - surface_points[i+1][next_theta].second).lengthsq() < d_max_sq)
			   f << "line "
			     << surface_points[i+1][j].second.x() << " " 
			     << surface_points[i+1][j].second.y() << " " 
			     << surface_points[i+1][j].second.z() << " " 
			     << surface_points[i+1][next_theta].second.x() << " " 
			     << surface_points[i+1][next_theta].second.y() << " " 
			     << surface_points[i+1][next_theta].second.z() << " " 
			     << "sea\n";
			if (clipper::Coord_orth(surface_points[i][next_theta].second - surface_points[i+1][next_theta].second).lengthsq() < d_max_sq)
			   f << "line "
			     << surface_points[i][next_theta].second.x() << " " 
			     << surface_points[i][next_theta].second.y() << " " 
			     << surface_points[i][next_theta].second.z() << " " 
			     << surface_points[i+1][next_theta].second.x() << " " 
			     << surface_points[i+1][next_theta].second.y() << " " 
			     << surface_points[i+1][next_theta].second.z() << " " 
			     << "sea\n";
		     }
		  }
	       }
	    }
	 }
      }
   }
   f.close();
} 


// start:   (58.042  6.9541 22.4479)
// path pt: (55.635 10.6790 17.4240)


int main(int argc, char **argv) {

   int status = 0;

   if (argc > 1) {
      std::string pdb_file_name = argv[1];
      mmdb::Manager *mol = new mmdb::Manager;
      mmdb::ERROR_CODE ec = mol->ReadCoorFile(pdb_file_name.c_str());
      if (! ec) {
	 coot::canyon c(mol);

// test-canyon.pdb	 
// 	 clipper::Coord_orth s(58.042,  6.9541, 22.4479);
// 	 clipper::Coord_orth p(55.635, 10.6790, 17.4240);

// test DNA	 
 	 clipper::Coord_orth s(-14.73, -16.04, 124.89);
 	 clipper::Coord_orth p(-16.14, -14.06, 125.99);

	 // e is the end point of the canyon
	 // d is a point "into" the canyon from the end point
	 clipper::Coord_orth e(-25.72, -20.88, 137.03);
 	 clipper::Coord_orth d(-26.34, -20.14, 136.03);
	 
	 c.set_start_point(s);
	 c.set_start_path_point(p);

	 c.set_end_point(e);
	 c.set_end_path_point(d);
	 
	 c.trace();
      }
   }
   return status;
} 
