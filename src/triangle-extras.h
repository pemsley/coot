/*
 * src/triangle-extras.h
 *
 * Copyright 2007 by University of York
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

// Triangle extras

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif


BEGIN_C_DECLS

// triangle things
int tri_count = 0; 

struct point3d {
   float x[3];
};

struct my_TRIANGLE { 
   struct point3d point[3]; 
};

#define MAX_TRIANGLES 1000000

struct my_TRIANGLE triangle_list[MAX_TRIANGLES];

void visible(int vis); 

static void idle(void); 

void inc_molecule_rot_angle(void);

END_C_DECLS

// Don't adjust this, add this information to graphics_info.
// and get rid of molecule_rot_t (or at least encapsulate 
// it in graphics_info_t.
//
class molecule_rot_t { 

 public:   
   static float x_axis_angle;
   static float y_axis_angle;
}; 
   
