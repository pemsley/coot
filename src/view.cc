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
#include <glm/gtc/quaternion.hpp> // needed?
#include <glm/gtx/string_cast.hpp>

// It now doesn't seem right that this is merely a member function of view_info_t
// because it touches graphics_info_t so much. Needs reworking to be a function
// of graphics_info_t that takes 2 views.
coot::view_info_t
coot::view_info_t::interpolate(const coot::view_info_t &view1,
			       const coot::view_info_t &view2_in,
			       int n_steps) {
   coot::view_info_t view;
   graphics_info_t g;
   view_info_t view2(view2_in);

   if (false) {
      std::cout << "start quat interpolation: zooms: " << view1.zoom << " " << view2.zoom
                << " and centres: "
                << view1.rotation_centre << " to " << view2.rotation_centre << std::endl;
      std::cout << "quaternion interpolation using " << n_steps << " steps"
                << std::endl;
   }

   float total_zoom_by = view2.zoom/view1.zoom;
   int smooth_scroll_state = graphics_info_t::smooth_scroll;
   graphics_info_t::smooth_scroll = 0;

   {
      double dp = dot_product(view1, view2);
      if (dot_product(view1, view2) < 0.0) {
	 view2.negate_quaternion();
	 dp = dot_product(view1, view2); // dp needs updating, else wierdness
      }
      double dd = (view1.quat_length()*view2.quat_length());
      double ff = dp/dd;
      if (ff > 1.0) ff = 1.0; // stabilize
      double omega = acos(ff);

      graphics_info_t::reorienting_residue_start_view = view1;
      graphics_info_t::reorienting_residue_end_view   = view2;

      // std::cout << "here with dd " << dd << " ff " << ff << " omega: " << omega << std::endl;

      // do we want omega == 0 to stop animation?

      if (omega != 0.0) {
	 // slerping
	 //
	 if (n_steps < 1)
	    n_steps = 1;
	 // double frac = double(1.0)/double(n_steps);
	 // for (double f=0; f<=1.0; f+=frac) {


         // replace the direct "updating of the positions and orientation, and redraw"
         // with a callback animation function

         // Do this with sinusoidal acceleration (rotation and translation) for more yum-factor

         // return a gboolean
         auto animation_func = [] (GtkWidget *widget,
                                   GdkFrameClock *frame_clock,
                                   gpointer data) {

                                  coot::view_info_t &view1 = graphics_info_t::reorienting_residue_start_view;
                                  coot::view_info_t &view2 = graphics_info_t::reorienting_residue_end_view;

                                  double dp = dot_product(view1, view2);
                                  if (dot_product(view1, view2) < 0.0) {
                                     view2.negate_quaternion();
                                     dp = dot_product(view1, view2); // dp needs updating, else wierdness
                                  }
                                  double dd = (view1.quat_length()*view2.quat_length());
                                  double ff = dp/dd;
                                  if (ff > 1.0) ff = 1.0; // stabilize
                                  double omega = acos(ff);

                                  gboolean do_continue = G_SOURCE_REMOVE;
                                  float frac = 1.0;
                                  int n_steps = graphics_info_t::smooth_scroll_steps;
                                  n_steps = 50;
                                  int i_current_step = graphics_info_t::smooth_scroll_current_step;
                                  graphics_info_t g; // for rotation centre debugging.
                                  if (n_steps > 0)
                                     frac = 1.0/static_cast<float>(n_steps);
                                  coot::Cartesian this_step_delta = graphics_info_t::smooth_scroll_delta * frac;
                                  graphics_info_t::smooth_scroll_current_step += 1; // update now
                                  if (i_current_step < n_steps) {
                                     graphics_info_t::add_vector_to_rotation_centre(this_step_delta);

                                     // now the orientation
                                     double f = static_cast<double>(i_current_step) / static_cast<double>(n_steps);
                                     double one_over_sin_omega = 1.0/sin(omega);
                                     double frac1 = sin((1.0-f)*omega) * one_over_sin_omega;
                                     double frac2 = sin(f*omega) * one_over_sin_omega;

                                     // for (int iq=0; iq<4; iq++)
                                     // g.quat[iq] = frac1*view1.quat[iq] + frac2*view2.quat[iq];

                                     glm::quat quat_start (view1.quaternion);
                                     glm::quat quat_target(view2.quaternion);

                                     // now update glm_quat
                                     float ff = static_cast<float>(f);
                                     glm::quat mixed = glm::mix(quat_start, quat_target, ff);
                                     graphics_info_t::glm_quat = glm::normalize(mixed);

                                     if (false) {
                                        std::cout << "lambda animation_func: this_step "
                                                  << graphics_info_t::smooth_scroll_current_step
                                                  << " this_step_delta: " << this_step_delta
                                                  << " for ff " << ff
                                                  << " Rotation centre now "
                                                  << g.RotationCentre() << std::endl;
                                     }

                                     if (false) {
                                        std::cout << " start  " << glm::to_string(quat_start)
                                                  << " target " << glm::to_string(quat_target)
                                                  << " mixed  " << glm::to_string(mixed)
                                                  << std::endl;
                                     }
                                     graphics_info_t::graphics_draw(); // adds to the queue
                                     do_continue = G_SOURCE_CONTINUE;
                                  } else {
                                     do_continue = G_SOURCE_REMOVE;
                                  }
                                  return do_continue;
                               };
         gpointer user_data = 0;

         graphics_info_t::smooth_scroll_current_step = 0; // reset
         gtk_widget_add_tick_callback(graphics_info_t::glareas[0], animation_func, user_data, NULL);

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
               g.glm_quat = view1.quaternion;
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
   return glm::sqrt(glm::dot(quaternion, quaternion));
}


// static
float 
coot::view_info_t::dot_product(const coot::view_info_t &view1,
			       const coot::view_info_t &view2) {

   return glm::dot(view1.quaternion, view2.quaternion);
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

      f << "   ";
      f << glm::to_string(view.quaternion);
      f << ",\n";
      
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

   if (zoom < view.zoom*(1+frac)) {
      if (zoom > view.zoom*(1-frac)) {
	 if (rotation_centre.x() < view.rotation_centre.x()*(1+xfrac)) {
	    if (rotation_centre.x() > view.rotation_centre.x()*(1-xfrac)) { 
	       if (rotation_centre.y() < view.rotation_centre.y()*(1+yfrac)) { 
		  if (rotation_centre.y() > view.rotation_centre.y()*(1-yfrac)) { 
		     if (rotation_centre.z() < view.rotation_centre.z()*(1+zfrac)) { 
			if (rotation_centre.z() > view.rotation_centre.z()*(1-zfrac)) {
                           // test similar quaternions here (deleted float-based code)
                           matches = true;
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
