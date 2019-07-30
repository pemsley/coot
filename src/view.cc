/* src/view.cc
 * 
 * Copyright 2009 by the University of Oxford
 * Copyright 2016 by Medical Research Council
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
#include "Python.h"
#endif // USE_PYTHON
#include "view.hh"
#include "graphics-info.h"

coot::view_info_t
coot::view_info_t::interpolate(const coot::view_info_t &view1,
			       const coot::view_info_t &view2_in,
			       int n_steps) {
   coot::view_info_t view;
   graphics_info_t g;
   view_info_t view2(view2_in);

//    std::cout << "start quat interpolation: zooms: " << view1.zoom << " " << view2.zoom
// 	     << " and centres: "
// 	     << view1.rotation_centre << " to " << view2.rotation_centre << std::endl;

//     std::cout << "quaternion interpolation using " << n_steps << " steps"
// 	      << std::endl;

   float total_zoom_by = view2.zoom/view1.zoom;
   int smooth_scroll_state = graphics_info_t::smooth_scroll;
   graphics_info_t::smooth_scroll = 0;

   if (true) {
      double dp = dot_product(view1, view2);
      if (dot_product(view1, view2) < 0.0) {
	 view2.negate_quaternion();
	 dp = dot_product(view1, view2); // dp needs updating, else wierdness
      }
      double dd = (view1.quat_length()*view2.quat_length());
      double ff = dp/dd;
      if (ff > 1.0) ff = 1.0; // stabilize
      double omega = acos(ff);

      // std::cout << "here with dd " << dd << " ff " << ff << " omega: " << omega << std::endl;

      if (omega != 0.0) { 
	 // slerping
	 //
	 if (n_steps < 1)
	    n_steps = 1;
	 // double frac = double(1.0)/double(n_steps);
	 // for (double f=0; f<=1.0; f+=frac) {
         for (int istep=0; istep<=n_steps; istep++) {
            double f = static_cast<double>(istep) / static_cast<double>(n_steps);
	    double one_over_sin_omega = 1/sin(omega);
	    double frac1 = sin((1-f)*omega) * one_over_sin_omega;
	    double frac2 = sin(f*omega) * one_over_sin_omega;
	    for (int iq=0; iq<4; iq++)
	       g.quat[iq] = frac1*view1.quat[iq] + frac2*view2.quat[iq];
	    coot::Cartesian rct =
	       view1.rotation_centre + (view2.rotation_centre - view1.rotation_centre).by_scalar(f);

	    // I don't want to setRotationCentre() because that sets the old rotation centre
	    // (and we need the original version of that to go "back")
	    // g.setRotationCentre(rct);
	    g.setRotationCentreSimple(rct);
	    
	    g.zoom = view1.zoom + pow(f,0.5)*(view2.zoom-view1.zoom);
	    // std::cout << "f " << f << " sqrt(t) " << sqrt(f)
	    // << " zoom " << g.zoom << "   " << rct << std::endl;

	    graphics_info_t::graphics_draw();
	 }
	 
      } else {
	 // non slerping

	 // do animation if the views don't match
	 //
	 bool do_animation = ! (view1.matches_view(view2));

	 if (do_animation) {
	    for (int i=0; i<=n_steps; i++) {
	       double frac = double(i)/double(n_steps);
	       coot::Cartesian rct =
		  view1.rotation_centre + (view2.rotation_centre - view1.rotation_centre).by_scalar(frac);
	       for (int iq=0; iq<4; iq++)
		  g.quat[iq] = view1.quat[iq] + frac*(view2.quat[iq]-view1.quat[iq]);
	       g.setRotationCentre(rct);
	       g.zoom = view1.zoom + frac*(view2.zoom-view1.zoom);
	       graphics_info_t::graphics_draw();
	    }
	 }
      }
   }

   graphics_info_t::smooth_scroll = smooth_scroll_state;
   return view;
}

float
coot::view_info_t::quat_length() const {

   float d = 0.0;
   for (int i=0; i<4; i++) {
      d += quat[i]*quat[i];
   }
   return sqrt(d);
}


// static
float 
coot::view_info_t::dot_product(const coot::view_info_t &view1,
			       const coot::view_info_t &view2) {

   float d = 0.0;
   for (int i=0; i<4; i++) {
      d += view1.quat[i]*view2.quat[i];
   }
   return d;
}

std::ostream&
coot::operator<<(std::ostream &f, const coot::view_info_t &view) {

   // position quaternion zoom view-name
   //


   std::cout << "debug: in view output operator(): view_name is \"" << view.view_name << "\"" << std::endl;

#ifdef USE_GUILE
   if (! view.is_simple_spin_view_flag) { 
      f << "(add-view ";
      f << "(list ";
      f << "   ";
      f << view.rotation_centre.x();
      f << " ";
      f << view.rotation_centre.y();
      f << " ";
      f << view.rotation_centre.z();
      f << ")\n";

      f << "   (list ";
      f << view.quat[0]; 
      f << " ";
      f << view.quat[1]; 
      f << " ";
      f << view.quat[2]; 
      f << " ";
      f << view.quat[3];
      f << ")\n";
      
      f << "   ";
      f << view.zoom; 
      f << "\n";

      f << "   ";
      f << coot::util::single_quote(view.view_name);
   
      f << ")\n";
   } else {
      f << "(add-spin-view ";
      f << coot::util::single_quote(view.view_name);
      f << " ";
      f << view.n_spin_steps;
      f << " ";
      f << view.degrees_per_step * view.n_spin_steps;
      f << ")\n";
   }
#else
#ifdef USE_PYTHON
   if (! view.is_simple_spin_view_flag) { 
      f << "add_view(";
      f << "[";
      f << view.rotation_centre.x();
      f << ", ";
      f << view.rotation_centre.y();
      f << ", ";
      f << view.rotation_centre.z();
      f << "],\n";

      f << "   [";
      f << view.quat[0]; 
      f << ", ";
      f << view.quat[1]; 
      f << ", ";
      f << view.quat[2]; 
      f << ", ";
      f << view.quat[3];
      f << "],\n";
      
      f << "   ";
      f << view.zoom; 
      f << ",\n";

      f << "   ";
      f << coot::util::single_quote(view.view_name);
   
      f << ")\n";
   } else {
      f << "add_spin_view(";
      f << coot::util::single_quote(view.view_name);
      f << ", ";
      f << view.n_spin_steps;
      f << ", ";
      f << view.degrees_per_step * view.n_spin_steps;
      f << ")\n";
   }
#endif // USE_PYTHON
#endif // USE_GUILE
   return f;

}

   
bool
coot::view_info_t::matches_view (const coot::view_info_t &view) const {
   
   float frac = 0.01;
   bool matches = false;
   // maybe there's a better way
   float xfrac, yfrac, zfrac, q0frac, q1frac, q2frac, q3frac;
   xfrac = yfrac = zfrac = q0frac = q1frac = q2frac = q3frac = frac;
   if (rotation_centre.x() < 0) xfrac = -xfrac;
   if (rotation_centre.y() < 0) yfrac = -yfrac;
   if (rotation_centre.z() < 0) zfrac = -zfrac;
   if (quat[0] < 0) q0frac = -q0frac;
   if (quat[1] < 0) q1frac = -q1frac;
   if (quat[2] < 0) q2frac = -q2frac;
   if (quat[3] < 0) q3frac = -q3frac;

   if (zoom < view.zoom*(1+frac)) {
      if (zoom > view.zoom*(1-frac)) {
	 if (rotation_centre.x() < view.rotation_centre.x()*(1+xfrac)) {
	    if (rotation_centre.x() > view.rotation_centre.x()*(1-xfrac)) { 
	       if (rotation_centre.y() < view.rotation_centre.y()*(1+yfrac)) { 
		  if (rotation_centre.y() > view.rotation_centre.y()*(1-yfrac)) { 
		     if (rotation_centre.z() < view.rotation_centre.z()*(1+zfrac)) { 
			if (rotation_centre.z() > view.rotation_centre.z()*(1-zfrac)) {
			   if (quat[0] < view.quat[0]*(1+q0frac)) {
			      if (quat[0] > view.quat[0]*(1-q0frac) ){ 
				 if (quat[1] < view.quat[1]*(1+q1frac)) { 
				    if (quat[1] > view.quat[1]*(1-q1frac) ){ 
				       if (quat[2] < view.quat[2]*(1+q2frac)) { 
					  if (quat[2] > view.quat[2]*(1-q2frac) ){ 
					     if (quat[3] < view.quat[3]*(1+q3frac)) { 
						if (quat[3] > view.quat[3]*(1-q3frac) ){
						   matches = true;
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
   return matches;
}
