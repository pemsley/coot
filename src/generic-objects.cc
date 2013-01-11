/* src/generic-objects.cc
 * 
 * Copyright 2005, 2006 The University of York
 * Copyright 2007, 2012 The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

#include <gtk/gtk.h>
#include <math.h>

#include "c-interface.h"
#include "cc-interface.hh"
#include "graphics-info.h"

/*  ----------------------------------------------------------------------- */
/*                  Generic Objects                                         */
/*  ----------------------------------------------------------------------- */
void to_generic_object_add_line(int object_number,
				const char *colour_name,
				int line_width,
				float from_x1, 
				float from_y1, 
				float from_z1, 
				float to_x2, 
				float to_y2, 
				float to_z2) {

   clipper::Coord_orth x1(from_x1, from_y1, from_z1);
   clipper::Coord_orth x2(to_x2, to_y2, to_z2);

   graphics_info_t g;
   std::pair<clipper::Coord_orth, clipper::Coord_orth> coords(x1, x2);

   std::string c(colour_name);
   coot::colour_holder colour =
      coot::generic_display_object_t::colour_values_from_colour_name(c);
   if (object_number >= 0) { 
      unsigned int object_number_u(object_number);
      if (object_number_u < g.generic_objects_p->size()) { 
	 (*g.generic_objects_p)[object_number].add_line(colour,
							c,
							line_width,
							coords);
	 
      } else {
	 std::cout << "BAD object_number in to_generic_object_add_line"
		   << " out of range high" << object_number << std::endl;
      }
   } else {
      std::cout << "BAD object_number (out of range low) in to_generic_object_add_line"
		<< object_number << std::endl;
   }
}

/*! \brief add line to generic object object_number */
void to_generic_object_add_dashed_line(int object_number, 
				       const char *colour,
				       int line_width,
				       float dash_density,
				       float from_x1, 
				       float from_y1, 
				       float from_z1, 
				       float to_x2, 
				       float to_y2, 
				       float to_z2) {

   clipper::Coord_orth p_start(from_x1, from_y1, from_z1);
   clipper::Coord_orth p_end(to_x2, to_y2, to_z2);
   float ll = clipper::Coord_orth::length(p_start, p_end);
   int n_dashes = int(dash_density * ll);
   bool visible = 1;
   
   for (unsigned int idash=0; idash<(n_dashes-1); idash++) {
      if (visible) { 
	 float fracs = float(idash)/float(n_dashes);
	 float fracn = float(idash+1)/float(n_dashes);
	 clipper::Coord_orth p1 = p_start + fracs * (p_end - p_start);
	 clipper::Coord_orth p2 = p_start + fracn * (p_end - p_start);
	 to_generic_object_add_line(object_number, colour, line_width,
				    p1.x(), p1.y(), p1.z(), 
				    p2.x(), p2.y(), p2.z());
      }
      visible = !visible;
   }
} 





void to_generic_object_add_point(int object_number, 
				 const char *colour_name,
				 int point_width,
				 float from_x1, 
				 float from_y1, 
				 float from_z1) {

   graphics_info_t g;
   clipper::Coord_orth x1(from_x1, from_y1, from_z1);
   std::string c(colour_name);
   coot::colour_holder colour =
      coot::generic_display_object_t::colour_values_from_colour_name(c);

//    std::cout << "debug:: colour input " << c << " gave colour "
// 	     << colour << std::endl;

   if (object_number >=0 && object_number < g.generic_objects_p->size()) { 

      (*g.generic_objects_p)[object_number].add_point(colour,
						      c,
						      point_width,
						      x1);
   } else {
      std::cout << "BAD object_number in to_generic_object_add_point: "
		<< object_number << std::endl;
   } 
}


/*! \brief add point to generic object object_number */
void to_generic_object_add_arc(int object_number, 
			       const char *colour_name,
			       float radius,
			       float radius_inner,
			       float from_angle,
			       float to_angle,
			       float start_point_x,
			       float start_point_y,
			       float start_point_z,
			       float start_dir_x,
			       float start_dir_y,
			       float start_dir_z,
			       float normal_x, 
			       float normal_y, 
			       float normal_z) {
   graphics_info_t g;
   if (object_number >=0 && object_number < g.generic_objects_p->size()) {
      coot::generic_display_object_t::arc_t arc(from_angle, to_angle,
						clipper::Coord_orth(start_point_x,
								    start_point_y,
								    start_point_z),
						clipper::Coord_orth(start_dir_x,
								    start_dir_y,
								    start_dir_z),
						clipper::Coord_orth(normal_x,
								    normal_y,
								    normal_z),
						radius, radius_inner);
      coot::colour_holder colour =
	 coot::generic_display_object_t::colour_values_from_colour_name(std::string(colour_name));
      // bleugh!  Make your colour holders consistent!  - use the utils version throughout!
      arc.col.col[0] = colour.red;
      arc.col.col[1] = colour.green;
      arc.col.col[2] = colour.blue;
      (*g.generic_objects_p)[object_number].arcs.push_back(arc);
   } else {
      std::cout << "BAD object_number in to_generic_object_add_arc: "
		<< object_number << std::endl;
   } 
}


void to_generic_object_add_display_list_handle(int object_number, int display_list_id) { 

   graphics_info_t g;
   if (object_number >=0 && object_number < g.generic_objects_p->size()) { 
      (*g.generic_objects_p)[object_number].GL_display_list_handles.push_back(display_list_id);
   } else {
      std::cout << "BAD object_number in to_generic_object_add_point: "
		<< object_number << std::endl;
   } 
   

}


void set_display_generic_object(int object_number, short int istate) {

   graphics_info_t g;
   if (object_number >=0  && object_number < g.generic_objects_p->size()) {
      (*g.generic_objects_p)[object_number].is_displayed_flag = istate;
   } else {
      std::cout << "BAD object_number in to_generic_object_add_point: "
		<< object_number << std::endl;
   }
   graphics_draw();
}

/*! \brief is generic display object displayed?

  @return 1 for yes, otherwise 0  */
int generic_object_is_displayed_p(int object_number) {

   int is_displayed = 0;
   graphics_info_t g;
   if (object_number >=0  && object_number < g.generic_objects_p->size()) {
      is_displayed = (*g.generic_objects_p)[object_number].is_displayed_flag;
   }
   return is_displayed;
}



int new_generic_object_number(const char *name) {

   graphics_info_t g;
   std::string n;
   if (name)
      n = std::string(name);
   return g.new_generic_object_number(n);

}


// return the index of the object with name name, if not, return -1;
// 
int generic_object_index(const char *name) {

   graphics_info_t g;
   int index = -1; 
   if (name) {
      std::string n(name);
      int nobjs = g.generic_objects_p->size();
      for (int iobj=0; iobj<nobjs; iobj++) {
	 if ((*g.generic_objects_p)[iobj].name == n) {
	    if (!(*g.generic_objects_p)[iobj].is_closed_flag) { 
	       index = iobj;
	       break;
	    }
	 }
      }
   }
   return index;
}


// OLD code passing back a const char (yeuch)
//
// /*! \brief what is the name of generic object number obj_number? 

// return 0 (NULL) #f  on obj_number not available */
// const char *generic_object_name(int obj_number) {

//    std::string s;
//    const char *r = 0;
//    graphics_info_t g;
//    int n_objs = g.generic_objects_p->size();
//    // funny way! (no good reason)
//    for (unsigned int i=0; i<n_objs; i++) {
//       if (i == obj_number)
// 	 r = (*g.generic_objects_p)[i].name.c_str();
//    }
//    return r;
// }



/*! \brief what is the name of generic object number obj_number? 
  return #f on obj_number not available */
#ifdef USE_GUILE
SCM generic_object_name_scm(int obj_number) {
   graphics_info_t g;
   int n_objs = g.generic_objects_p->size();
   SCM r = SCM_BOOL_F;
   for (int i=(n_objs-1); i>=0; i--) {
      if (i == obj_number) {
	 if (!(*g.generic_objects_p)[i].is_closed_flag) { 
	    r = scm_makfrom0str((*g.generic_objects_p)[i].name.c_str());
	 }
      }
   }
   return r;
}
#endif /* USE_GUILE */

#ifdef USE_PYTHON
PyObject *generic_object_name_py(int obj_number) {
   graphics_info_t g;
   int n_objs = g.generic_objects_p->size();
   PyObject *r;
   r = Py_False;
   for (int i=(n_objs-1); i>=0; i--) {
      if (i == obj_number) {
	 if (!(*g.generic_objects_p)[i].is_closed_flag) { 
	    r = PyString_FromString((*g.generic_objects_p)[i].name.c_str());
	    break;
	 }
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
} 
#endif /* USE_PYTHON */


/*! \brief clear out the lines and points from object_number, but keep
  it displayable (not closed). */
void generic_object_clear(int object_number) {

   graphics_info_t g;
   if (object_number >= 0) {
      if (object_number < int(g.generic_objects_p->size())) {
	 (*g.generic_objects_p)[object_number].clear();
      }
   } 
 

}



/*! \brief pass a filename that contains molprobity's probe output in XtalView 
format */
void handle_read_draw_probe_dots(const char *dots_file) {

   // std::cout << "handle reading dots file" << dots_file << std::endl;
   if (dots_file) {

      FILE* dots = fopen(dots_file, "r" );
      if ( dots == NULL ) { 
	 std::cout << "handle_read_draw_probe_dots  - Could not read: "
		   << dots_file << std::endl;
	 // fclose(dots);
      } else {

	 // clear up what probe contacts/overlaps we already have:
	 std::vector<std::string> deletable_names;
	 deletable_names.push_back("wide contact");
	 deletable_names.push_back("close contact");
	 deletable_names.push_back("small overlap");
	 deletable_names.push_back("bad overlap");
	 deletable_names.push_back("H-bonds");
	 int nobjs = graphics_info_t::generic_objects_p->size();
	 for (unsigned int i=0; i< nobjs; i++) {
	    for (unsigned int d=0; d<deletable_names.size(); d++) { 
	       if ((*graphics_info_t::generic_objects_p)[i].name == deletable_names[d]) {
		  close_generic_object(i); // empty it, really
	       }
	    }
	 }
	 
	 int n_lines = 0;
	 int n_points = 0;
	 std::string current_colour = "blue"; // should be reset.
	 std::string current_name   = "Unassigned";
	 int obj_no = number_of_generic_objects();
	 
	 char line[240];
	 char s[240];
	 char s1[240];
	 char s2[240];
	 char s3[240];
	 float x1, x2, x3, x4, x5, x6;
	 while ( fgets( line, 240, dots ) != NULL ) {
	    if (sscanf(line, "# %s %s %s", s1, s2, s3)) {
	       std::string st1, st2, st3;
	       st1 = s1;
	       st2 = s2;
	       st3 = s3;
	       // a comment line, what type of dots are we looking at?
	       // std::cout << "# line --- " << line << std::endl;
	       clipper::String contact_line(s);
	       std::vector<std::string> vs;
	       vs.push_back(st1);
	       vs.push_back(st2);
	       vs.push_back(st3);
	       int n_bits = 3;

	       std::pair<short int, std::string> p = is_interesting_dots_object_next_p(vs);
	       if (p.first) {
		  if (p.second != current_name) {
		     int maybe_old_object = generic_object_index(p.second.c_str());
		     if (maybe_old_object > -1) {
			obj_no = maybe_old_object;
		     } else { 
			obj_no = new_generic_object_number(p.second.c_str());
		     }
		     // non-member function usage, so that we don't do the redraw.
		     (*graphics_info_t::generic_objects_p)[obj_no].is_displayed_flag = 1;
		     (*graphics_info_t::generic_objects_p)[obj_no].is_closed_flag = 0;
		     current_name = p.second;
		  }
	       }
	       
	    } else {
	       if (sscanf(line, "%f %f %f %f %f %f %s", &x1, &x2, &x3, &x4, &x5, &x6, s)) {
		  current_colour = s;
		  float length2 = pow((x1-x4),2) + pow((x2-x5),2) + pow((x3-x6),2);
		  if (length2 > 0.1) {
		     n_lines++;
		     to_generic_object_add_line(obj_no, current_colour.c_str(), 3,
						x1, x2, x3, x4, x5, x6);
		  } else {
		     n_points++;
		     to_generic_object_add_point(obj_no, current_colour.c_str(), 2,
						 x1, x2, x3);
		  }
	       } else {
		  if (strlen(line) > 0) 
		     std::cout << ":" << line << ": failed to scan" << std::endl;
	       }
	    }
	 }
	 // std::cout << "INFO:: added " << n_lines << " lines and " << n_points << " points\n";
      }
      fclose(dots);
   }
}


void handle_read_draw_probe_dots_unformatted(const char *dots_file, int imol,
					     int show_clash_gui_flag) {

   int dot_size_scale_factor = 1; // user param?

   std::vector<coot::atom_spec_t> clash_atoms;
   bool hybrid_36_enabled_probe_flag = 1; // new style
   //bool hybrid_36_enabled_probe_flag = 0; // old style
   
   if (dots_file) {

      FILE* dots = fopen(dots_file, "r" );
      if ( dots == NULL ) { 
	 std::cout << "handle_read_draw_probe_dots  - Could not read: "
		   << dots_file << std::endl;
	 // fclose(dots);
      } else {

	 // clear up what probe contacts/overlaps we already have:
	 std::vector<std::string> deletable_names;
	 deletable_names.push_back("wide contact");
	 deletable_names.push_back("close contact");
	 deletable_names.push_back("small overlap");
	 deletable_names.push_back("bad overlap");
	 deletable_names.push_back("H-bonds");
	 int nobjs = graphics_info_t::generic_objects_p->size();
	 for (unsigned int i=0; i< nobjs; i++) {
	    for (unsigned int d=0; d<deletable_names.size(); d++) { 
	       if ((*graphics_info_t::generic_objects_p)[i].name == deletable_names[d]) {
		  close_generic_object(i); // empty it, really
	       }
	    }
	 }

	 int n_input_lines = 0; 
	 int n_lines = 0;
	 int n_points = 0;
	 std::string current_colour = "blue"; // should be reset.
	 int obj_no = number_of_generic_objects();

	 int imol_source, imol_target;
	 float gap1, factor3, gap2, factor4, factor5, factor6;
	 char line[240];
	 char c_type[20];
	 char atom_id1[30];
	 char atom_id2[30];
	 char contact_type1[200];
	 char contact_type2[3];
	 float x1, x2, x3, x4, x5, x6;
	 std::string current_contact_type = "none";
	 int dot_size = 2; // 3 for hydrogen bonds

	 // null string arrays
	 for (int i=0; i<20; i++) c_type[i] = 0;
	 for (int i=0; i<30; i++) atom_id1[i] = 0;
	 for (int i=0; i<30; i++) atom_id2[i] = 0;
	 std::string current_useful_name;

	 std::string scan_line = ":%d->%d:%2c:%15c:%15c:%f:%f:%f:%f:%f:%f:%f:%s";
	 if (hybrid_36_enabled_probe_flag) 
	    scan_line = ":%d->%d:%2c:%16c:%16c:%f:%f:%f:%f:%f:%f:%f:%s";	 

	 while ( fgets( line, 240, dots ) != NULL ) {
	    n_input_lines++; 
	    for (int i=0; i<3; i++) contact_type1[i] = 0;
	    for (int i=0; i<3; i++) contact_type2[i] = 0;


	    if (sscanf(line,
		       scan_line.c_str(),
		       &imol_source, &imol_target, c_type, atom_id1, atom_id2,
		       &gap1, &gap2, &x1, &x2, &x3, &factor3, &factor4,
		       contact_type1)) {

	       if (strlen(contact_type1) > 2) {
		  int colon_count = 0;
		  int offset = 0; 

		  // std::cout << "colon search |" << contact_type1 << std::endl;
		  for (int i=0; i<10; i++) {
		     if (contact_type1[i] == ':')
			colon_count++;
		     if (colon_count == 2) {
			offset = i + 1;
			break;
		     }
		  }

		  if (offset > 0) { 
		     // std::cout << "scanning |" << contact_type1+offset << std::endl;
		     
		     if (sscanf((contact_type1+offset), "%f:%f:%f:%f:%f", 
				&x4, &x5, &x6, &factor5, &factor6)) {
			
			std::string atom_id_1_str(atom_id1);
			std::string atom_id_2_str(atom_id2);
			std::string contact_type(c_type);
			
//  			std::cout << "\"" << contact_type << "\"..\"" << atom_id_1_str
// 				  << "\"..to..\""
//  				  << atom_id_2_str << "\" " << gap1 << " " << gap2
// 				  << " (" << x1 << "," << x2 << "," << x3 << ")"
//  				  << " (" << x4 << "," << x5 << "," << x6 << ")"
//  				  << "\n";
//  			std::cout << "\"" << contact_type << "\"..\"" << atom_id_1_str
//  				  << "\"..to..\""
//  				  << atom_id_2_str << "\" " << gap1 << " " << gap2
//  				  << " (" << x1 << "," << x2 << "," << x3 << ")"
//  				  << " (" << x4 << "," << x5 << "," << x6 << ")"
//  				  << "\n";
	    
			// assign the colour and dot size
			if (contact_type == "hb") { 
			   current_colour = "greentint";
			   dot_size = 3;
			} else {
			   dot_size = 2;
			   current_colour = "hotpink";
			   if (gap2> -0.4) current_colour = "red";
			   if (gap2> -0.3) current_colour = "orange";
			   if (gap2> -0.2) current_colour = "yellow";
			   if (gap2> -0.1) current_colour = "yellowtint";
			   if (gap2>  0.0) current_colour = "green";
			   if (gap2> 0.15) current_colour = "sea";
			   if (gap2> 0.25) current_colour = "sky";
			   if (gap2> 0.35) current_colour = "blue";
			}

			// do we need a new object?
			if (contact_type != current_contact_type) {

			   current_useful_name = probe_dots_short_contact_name_to_expanded_name(contact_type);
			   // do we have an object of that name already?
			   int maybe_old_object = generic_object_index(current_useful_name.c_str());
			   current_contact_type = contact_type;

			   if (maybe_old_object > -1) {
			      obj_no = maybe_old_object;
			   } else { 
			      obj_no = new_generic_object_number(current_useful_name.c_str());
			      // std::cout << "changing type to " << contact_type << std::endl;
			   }
			   // non-member function usage, so that we don't do the redraw.
			   (*graphics_info_t::generic_objects_p)[obj_no].is_displayed_flag = 1;
			   (*graphics_info_t::generic_objects_p)[obj_no].is_closed_flag = 0;
			}

			float length2 = pow((x1-x4),2) + pow((x2-x5),2) + pow((x3-x6),2);
			if (length2 > 0.04) {

			   // hydrogen bond are not drawn as spikes,
			   // they should be drawn as pillow surfaces.
			   if (contact_type == "hb") {
			      to_generic_object_add_point(obj_no, current_colour.c_str(),
							  dot_size * dot_size_scale_factor,
							  x1, x2, x3);
			      to_generic_object_add_point(obj_no, current_colour.c_str(),
							  dot_size * dot_size_scale_factor,
							  x4, x5, x6);
			   } else { 
			      n_lines++;
			      to_generic_object_add_line(obj_no, current_colour.c_str(), 3 * dot_size_scale_factor,
							 x1, x2, x3, x4, x5, x6);
			   }
			} else {
			   n_points++;
			   to_generic_object_add_point(obj_no, current_colour.c_str(), dot_size * dot_size_scale_factor,
						       x1, x2, x3);
			}

			if (length2 > 5) {
			   // a really long
			   std::cout << "DEBUG:: long line at line number "
				     << n_input_lines
				     << "  (" << x1 << "," << x2 << "," << x3 << ")"
				     << "  (" << x4 << "," << x5 << "," << x6 << ")"
				     << std::endl;
			   std::cout << line;
			}

			if (gap2 < -0.42) {
			   // it's a bad contact.  Add the 1st clash atom to the list of bad clashes.
		  
			   std::string chain_id = atom_id_1_str.substr(0,2);
			   std::string atom_name = atom_id_1_str.substr(11,4);
			   std::string insertion_code = atom_id_1_str.substr(6,1);
			   std::string altconf = atom_id_1_str.substr(15,1);
			   if (altconf == " ")
			      altconf = "";
		     
			   int resno = atoi(atom_id_1_str.substr(2,4).c_str());

			   // std::cout << "   makes atom " << chain_id << ":" << resno << ":"
			   // << insertion_code << ":" << atom_name << ":"
			   // << altconf << ":" << "\n";
			   
			   coot::atom_spec_t r(chain_id, resno, insertion_code, atom_name, altconf);

			   short int ifound = 0;
			   for (unsigned int ic=0; ic<clash_atoms.size(); ic++) {
			      if (resno == clash_atoms[ic].resno) {
				 if (insertion_code == clash_atoms[ic].insertion_code) {
				    if (altconf == clash_atoms[ic].alt_conf) {
				       if (chain_id == clash_atoms[ic].chain) {
				 
					  // we are intested in bad (low) gaps
					  if (gap2 < clash_atoms[ic].float_user_data)
					     clash_atoms[ic].float_user_data = gap2;
				 
					  ifound = 1;
					  break;
				       }
				    }
				 }
			      }
			   }
			   if (ifound == 0) {
			      r.int_user_data = imol;
			      r.float_user_data = gap2;
			      r.string_user_data = current_useful_name;
			      // don't put hydrogen bonds in the bad contact list
			      if (contact_type != "hb") 
				 clash_atoms.push_back(r);
			   }
			} 
		     }
		  }
	       }
	    }
	 }
	 if (show_clash_gui_flag) {
	    if (graphics_info_t::use_graphics_interface_flag) {
		if (show_clash_gui_flag == 2) {
#ifdef USE_PYGTK
		   graphics_info_t g;
		   std::vector<std::string> cmd_strings;
		   cmd_strings.push_back("run_with_gtk_threading");
		   cmd_strings.push_back("interesting_things_gui");
		   cmd_strings.push_back(single_quote("Molprobity Probe Clash Gaps"));
		   std::sort(clash_atoms.begin(), clash_atoms.end(), coot::compare_atom_specs_user_float);
		   std::string ls = coot::util::interesting_things_list_py(clash_atoms);
		   cmd_strings.push_back(ls);
		   std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
		   safe_python_command(s);
#else
#if defined USE_GUILE && !defined WINDOWS_MINGW
		   graphics_info_t g;
		   std::vector<std::string> cmd_strings;
		   cmd_strings.push_back("interesting-things-gui");
		   cmd_strings.push_back(single_quote("Molprobity Probe Clash Gaps"));
		   std::sort(clash_atoms.begin(), clash_atoms.end(), coot::compare_atom_specs_user_float);
		   std::string ls = coot::util::interesting_things_list(clash_atoms);
		   cmd_strings.push_back(ls);
		   std::string s = g.state_command(cmd_strings, coot::STATE_SCM);
		   safe_scheme_command(s);
#endif // GUILE
#endif // PYGTK
		} else {
#if defined USE_GUILE && !defined WINDOWS_MINGW
		   graphics_info_t g;
		   std::vector<std::string> cmd_strings;
		   cmd_strings.push_back("interesting-things-gui");
		   cmd_strings.push_back(single_quote("Molprobity Probe Clash Gaps"));
		   std::sort(clash_atoms.begin(), clash_atoms.end(), coot::compare_atom_specs_user_float);
		   std::string ls = coot::util::interesting_things_list(clash_atoms);
		   cmd_strings.push_back(ls);
		   std::string s = g.state_command(cmd_strings, coot::STATE_SCM);
		   safe_scheme_command(s);
#else
#ifdef USE_PYGTK
		   graphics_info_t g;
		   std::vector<std::string> cmd_strings;
		   cmd_strings.push_back("run_with_gtk_threading");
		   cmd_strings.push_back("interesting_things_gui");
		   cmd_strings.push_back(single_quote("Molprobity Probe Clash Gaps"));
		   std::sort(clash_atoms.begin(), clash_atoms.end(), coot::compare_atom_specs_user_float);
		   std::string ls = coot::util::interesting_things_list_py(clash_atoms);
		   cmd_strings.push_back(ls);
		   std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
		   safe_python_command(s);
#endif // PYGTK
#endif // GUILE
		}
	    }
	 }
         fclose(dots);
      }
   }
}
