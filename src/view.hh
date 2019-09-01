/* src/graphics-info.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by the University of York
 * Copyright 2007, 2008, 2009 by the University of Oxford
 * Copyright 2016 by Medical Research Council
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

#ifndef VIEW_HH
#define VIEW_HH

#include "coords/Cartesian.h"

namespace coot {
   
   // This is for a Pymol-like feature to switch between views:
   // 
   class view_info_t {
   public:
      float zoom;
      coot::Cartesian rotation_centre;
      std::string view_name;
      std::string description;
      bool is_simple_spin_view_flag;
      bool is_action_view_flag;
      int n_spin_steps;
      float degrees_per_step;
      float quat[4];
      std::string action;
      view_info_t(float *quat_in, const coot::Cartesian &rot_centre_in,
		  float zoom_in, const std::string &view_name_in) {
	is_simple_spin_view_flag = 0;
	is_action_view_flag = 0;
	 zoom = zoom_in;
	 rotation_centre = rot_centre_in;
	 view_name = view_name_in;
	 for (int i=0; i<4; i++) 
	    quat[i] = quat_in[i];
      }
      view_info_t(const view_info_t &v_in) {
	 zoom = v_in.zoom;
	 rotation_centre = v_in.rotation_centre;
	 description = v_in.description;
	 is_simple_spin_view_flag = v_in.is_simple_spin_view_flag;
	 is_action_view_flag = v_in.is_action_view_flag;
	 n_spin_steps = v_in.n_spin_steps;
	 degrees_per_step = v_in.degrees_per_step;
	 action = v_in.action;
         view_name = v_in.view_name;
	 for (std::size_t i=0; i<4; i++)
	    quat[i] = v_in.quat[i];
      }
      view_info_t() {
      	is_simple_spin_view_flag = 0;
	is_action_view_flag = 0;
      }
      // a spin view 
      view_info_t(const std::string &view_name_in, int nsteps, float degrees_total) {
	is_simple_spin_view_flag = 1;
	is_action_view_flag = 0;
	view_name = view_name_in;
	n_spin_steps = nsteps;
	if (n_spin_steps > 0) 
	  degrees_per_step = degrees_total/n_spin_steps;
	else
	  degrees_per_step = 0.5;
      }
      // an action view
      view_info_t(const std::string &view_name_in, const std::string &funct) {
	is_action_view_flag = 1;
	view_name = view_name_in;
	action = funct;
      }
      bool matches_view (const coot::view_info_t &view) const;
      float quat_length() const;
      void add_description(const std::string &descr) { 
	description = descr;
      }
      // -q is the same rotation as q, but we don't want to go the long way round
      // when moving from quat_1 to quat_2, if the dot_product is negative,
      // we need to negate_quaternion().
      void negate_quaternion() {
	 for (unsigned int i=0; i<4; i++)
	    quat[i] = -quat[i];
      }
      static view_info_t interpolate(const view_info_t &view1,
				     const view_info_t &view2,
				     int n_steps);
      static float dot_product(const view_info_t &view1,
			       const view_info_t &view2);

      friend std::ostream& operator<<(std::ostream &stream, 
				      const view_info_t &view);
   };
   std::ostream& operator<<(std::ostream &stream, 
			    const view_info_t &view);
}

#endif // VIEW_HH
