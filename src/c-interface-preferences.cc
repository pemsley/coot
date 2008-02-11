/* src/main.cc
 * 
 * Copyright 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007 The University of Oxford
 * Author: Paul Emsley, Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_IRUSR S_IREAD
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define sleep Sleep
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER

#include <dirent.h>   // for extra scheme dir


#include "graphics-info.h"
#include "interface.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "coot-preferences.h"
#ifdef USE_PYTHON
#include "Python.h"
#endif // PYTHON

void preferences() {

   GtkWidget *w = create_preferences();
   gtk_widget_show(w);
}


void clear_preferences() {

   graphics_info_t::preferences_widget = NULL;

}

void set_mark_cis_peptides_as_bad(int istate) {

   graphics_info_t::mark_cis_peptides_as_bad_flag = istate;
}

int show_mark_cis_peptides_as_bad_state() {

   return graphics_info_t::mark_cis_peptides_as_bad_flag;

}

#if (GTK_MAJOR_VERSION > 1)
void show_hide_preferences_tabs(GtkToggleToolButton *toggletoolbutton, int preference_type) {

  GtkWidget *frame;

  std::vector<std::string> preferences_tabs;
  
  if (preference_type == COOT_GENERAL_PREFERENCES) {
    preferences_tabs = *graphics_info_t::preferences_general_tabs;
  }
  if (preference_type == COOT_BOND_PREFERENCES) {
    preferences_tabs = *graphics_info_t::preferences_bond_tabs;
  }
  if (preference_type == COOT_GEOMETRY_PREFERENCES) {
    preferences_tabs = *graphics_info_t::preferences_geometry_tabs;
  }
  if (preference_type == COOT_COLOUR_PREFERENCES) {
    preferences_tabs = *graphics_info_t::preferences_colour_tabs;
  }
  if (preference_type == COOT_MAP_PREFERENCES) {
    preferences_tabs = *graphics_info_t::preferences_map_tabs;
  }
  if (preference_type == COOT_OTHER_PREFERENCES) {
    preferences_tabs = *graphics_info_t::preferences_other_tabs;
  }

  for (int i=0; i<preferences_tabs.size(); i++) {
    frame = lookup_widget(GTK_WIDGET(toggletoolbutton), preferences_tabs[i].c_str());
    if (gtk_toggle_tool_button_get_active(toggletoolbutton)){
	  gtk_widget_show(frame);
	} else {
	  gtk_widget_hide(frame);
	}
  }

}
#endif // GTK_MAJOR_VERSION
 

GtkWidget *popup_window(const char *str) {

   GtkWidget *w = create_popup_info_window();
   GtkWidget *label = lookup_widget(w, "info_label");
   gtk_label_set_text(GTK_LABEL(label), str);
   return w;
}

void add_status_bar_text(const char *s) {

   graphics_info_t g;
   g.statusbar_text(std::string(s));
} 



/*  ----------------------------------------------------------------------- */
/*                  Other interface preferences                            */
/*  ----------------------------------------------------------------------- */

void set_model_fit_refine_dialog_stays_on_top(int istate) { 
   graphics_info_t::model_fit_refine_dialog_stays_on_top_flag = istate;
}

int model_fit_refine_dialog_stays_on_top_state() {

   return graphics_info_t::model_fit_refine_dialog_stays_on_top_flag;
} 


void save_accept_reject_dialog_window_position(GtkWidget *acc_rej_dialog) {

   // 20070801 crash reported by "Gajiwala, Ketan"

   // OK, we can reproduce a problem
   // Refine something
   // Close the window using WM delete window
   // Press return in Graphics window (globjects:key_press_event() GDK_Return case)
   // 
   // So, we need to set graphics_info_t::accept_reject_dialog to NULL
   // when we get a WM delete event on the Accept/Reject box
   
   if (acc_rej_dialog) { 
      gint upositionx, upositiony;
      if (acc_rej_dialog->window) {
	 gdk_window_get_root_origin (acc_rej_dialog->window, &upositionx, &upositiony);
	 graphics_info_t::accept_reject_dialog_x_position = upositionx;
	 graphics_info_t::accept_reject_dialog_y_position = upositiony;
      } else {
	 std::cout << "ERROR:: Trapped an error in save_accept_reject_dialog_window_position\n"
		   << "        Report to Central Control!\n"
		   << "        (What did you do to make this happen?)\n";
      }
   }
}

void set_accept_reject_dialog(GtkWidget *w) { /* used by callbacks to unset the widget.
						 (errr... it wasn't but it is now (as it
						 should be)). */

   graphics_info_t::accept_reject_dialog = w;
}

/* \brief set position of Model/Fit/Refine dialog */
void set_model_fit_refine_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::model_fit_refine_x_position = x_pos;
   graphics_info_t::model_fit_refine_y_position = y_pos;
}

/* \brief set position of Display Control dialog */
void set_display_control_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::display_manager_x_position = x_pos;
   graphics_info_t::display_manager_y_position = y_pos;
}

/* \brief set position of Go To Atom dialog */
void set_go_to_atom_window_position(int x_pos, int y_pos) {

   graphics_info_t::go_to_atom_window_x_position = x_pos;
   graphics_info_t::go_to_atom_window_y_position = y_pos;
}

/* \brief set position of Delete dialog */
void set_delete_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::delete_item_widget_x_position = x_pos;
   graphics_info_t::delete_item_widget_y_position = y_pos;
}

void set_rotate_translate_dialog_position(int x_pos, int y_pos) {

   graphics_info_t::rotate_translate_x_position = x_pos;
   graphics_info_t::rotate_translate_y_position = y_pos;
}

/*! \brief set position of the Accept/Reject dialog */
void set_accept_reject_dialog_position(int x_pos, int y_pos) {
   graphics_info_t::accept_reject_dialog_x_position = x_pos;
   graphics_info_t::accept_reject_dialog_y_position = y_pos;
}

/*! \brief set position of the Ramachadran Plot dialog */
void set_ramachandran_plot_dialog_position(int x_pos, int y_pos) {
   graphics_info_t::ramachandran_plot_x_position = x_pos;
   graphics_info_t::ramachandran_plot_y_position = y_pos;
}


/*  ------------------------------------------------------------------------ */
/*                     state (a graphics_info thing)                         */
/*  ------------------------------------------------------------------------ */
void set_save_state_file_name(const char *filename) {
   graphics_info_t::save_state_file_name = filename; 
}

const char *save_state_file_name_raw() {
   return graphics_info_t::save_state_file_name.c_str();
}


#ifdef USE_GUILE
SCM save_state_file_name_scm() {

//    char *f = (char *) malloc(graphics_info_t::save_state_file_name.length() +1);
//    strcpy(f, graphics_info_t::save_state_file_name.c_str());
//    return f;

   std::string f = graphics_info_t::save_state_file_name;
   return scm_makfrom0str(f.c_str());
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *save_state_file_name_py() {
   std::string f = graphics_info_t::save_state_file_name;
   return PyString_FromString(f.c_str());
}
#endif // USE_PYTHON


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
	 std::cout << "BAD object_number in to_generic_object_add_line" << std::endl;
      }
   } else {
      std::cout << "BAD object_number in to_generic_object_add_line" << std::endl;
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
   PyObject *r = Py_False;
   for (int i=(n_objs-1); i>=0; i--) {
      if (i == obj_number) {
	 if (!(*g.generic_objects_p)[i].is_closed_flag) { 
	    r = PyString_FromString((*g.generic_objects_p)[i].name.c_str());
	 }
      }
   }
   return r;
} 
#endif /* USE_PYTHON */




/*! \brief pass a filename that contains molprobity's probe output in XtalView 
format */
void handle_read_draw_probe_dots(const char *dots_file) {

   // std::cout << "handle reading dots file" << dots_file << std::endl;
   if (dots_file) {

      FILE* dots = fopen(dots_file, "r" );
      if ( dots == NULL ) { 
	 std::cout << "handle_read_draw_probe_dots  - Could not read: "
		   << dots_file << std::endl;
	 fclose(dots);
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
	 int n_points = 0;	 std::string current_colour = "blue"; // should be reset.
	 int obj_no = number_of_generic_objects();
	 std::string current_name = "Unassigned";
	 
	 char line[240];
	 char s[240];
	 char s1[240];
	 char s2[240];
	 char s3[240];
	 float x1, x2, x3, x4, x5, x6;
	 while ( fgets( line, 240, dots ) != NULL ) {
	    if (sscanf(line, "# %s %s %s", s1, s2, s3)) {
	       std::string st1, st2, st3;
	       if (s1)
		  st1 = s1;
	       if (s2)
		  st2 = s2;
	       if (s3)
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
		  if (s)
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
   } 
}


void handle_read_draw_probe_dots_unformatted(const char *dots_file, int imol,
					     int show_clash_gui_flag) {

   std::vector<coot::atom_spec_t> clash_atoms;
   
   if (dots_file) {

      FILE* dots = fopen(dots_file, "r" );
      if ( dots == NULL ) { 
	 std::cout << "handle_read_draw_probe_dots  - Could not read: "
		   << dots_file << std::endl;
	 fclose(dots);
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

	 while ( fgets( line, 240, dots ) != NULL ) {
	    n_input_lines++; 
	    for (int i=0; i<3; i++) contact_type1[i] = 0;
	    for (int i=0; i<3; i++) contact_type2[i] = 0;


	    if (sscanf(line,
		       ":%d->%d:%2c:%15c:%15c:%f:%f:%f:%f:%f:%f:%f:%s",
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
			
// 			std::cout << "\"" << contact_type << "\"..\"" << atom_id_1_str
// 				  << "\"..to..\""
// 				  << atom_id_2_str << "\" " << gap1 << " " << gap2
// 				  << " (" << x1 << "," << x2 << "," << x3 << ")"
// 				  << " (" << x4 << "," << x5 << "," << x6 << ")"
// 				  << "\n";

	    
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
			   n_lines++;
			   to_generic_object_add_line(obj_no, current_colour.c_str(), 3,
						      x1, x2, x3, x4, x5, x6);
			} else {
			   n_points++;
			   to_generic_object_add_point(obj_no, current_colour.c_str(), dot_size,
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
		  
			   std::string chain_id = atom_id_1_str.substr(0,1);
			   std::string atom_name = atom_id_1_str.substr(10,4);
			   std::string insertion_code = atom_id_1_str.substr(5,1);
			   std::string altconf = atom_id_1_str.substr(14,1);
			   if (altconf == " ")
			      altconf = "";
		     
			   int resno = atoi(atom_id_1_str.substr(1,4).c_str());

			   // 		  std::cout << "   makes atom " << chain_id << ":" << resno << ":"
			   // 			    << insertion_code << ":" << atom_name << ":" << altconf << ":" << "\n";
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
#ifdef USE_GUILE
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
// BL says:: we need some python code here too!!!
#ifdef USE_PYGTK
           graphics_info_t g;
           std::vector<std::string> cmd_strings;
           cmd_strings.push_back("interesting_things_gui");
           cmd_strings.push_back(single_quote("Molprobity Probe Clash Gaps"));
           std::sort(clash_atoms.begin(), clash_atoms.end(), coot::compare_atom_specs_user_float);
           std::string ls = coot::util::interesting_things_list(clash_atoms);
           cmd_strings.push_back(ls);
           std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
           safe_python_command(s);
#endif // PYGTK
#endif // GUILE
	 }
      }
   }
}

#ifdef USE_PYTHON
void handle_read_draw_probe_dots_unformatted_py(const char *dots_file, int imol,
					     int show_clash_gui_flag) {

   std::vector<coot::atom_spec_t> clash_atoms;
   
   if (dots_file) {

      FILE* dots = fopen(dots_file, "r" );
      if ( dots == NULL ) { 
	 std::cout << "handle_read_draw_probe_dots  - Could not read: "
		   << dots_file << std::endl;
	 fclose(dots);
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

	 while ( fgets( line, 240, dots ) != NULL ) {
	    n_input_lines++; 
	    for (int i=0; i<3; i++) contact_type1[i] = 0;
	    for (int i=0; i<3; i++) contact_type2[i] = 0;


	    if (sscanf(line,
		       ":%d->%d:%2c:%15c:%15c:%f:%f:%f:%f:%f:%f:%f:%s",
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
			
// 			std::cout << "\"" << contact_type << "\"..\"" << atom_id_1_str
// 				  << "\"..to..\""
// 				  << atom_id_2_str << "\" " << gap1 << " " << gap2
// 				  << " (" << x1 << "," << x2 << "," << x3 << ")"
// 				  << " (" << x4 << "," << x5 << "," << x6 << ")"
// 				  << "\n";

	    
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
			   n_lines++;
			   to_generic_object_add_line(obj_no, current_colour.c_str(), 3,
						      x1, x2, x3, x4, x5, x6);
			} else {
			   n_points++;
			   to_generic_object_add_point(obj_no, current_colour.c_str(), dot_size,
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
		  
			   std::string chain_id = atom_id_1_str.substr(0,1);
			   std::string atom_name = atom_id_1_str.substr(10,4);
			   std::string insertion_code = atom_id_1_str.substr(5,1);
			   std::string altconf = atom_id_1_str.substr(14,1);
			   if (altconf == " ")
			      altconf = "";
		     
			   int resno = atoi(atom_id_1_str.substr(1,4).c_str());

			   // 		  std::cout << "   makes atom " << chain_id << ":" << resno << ":"
			   // 			    << insertion_code << ":" << atom_name << ":" << altconf << ":" << "\n";
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
#ifdef USE_PYGTK
           graphics_info_t g;
           std::vector<std::string> cmd_strings;
           cmd_strings.push_back("interesting_things_gui");
           cmd_strings.push_back(single_quote("Molprobity Probe Clash Gaps"));
           std::sort(clash_atoms.begin(), clash_atoms.end(), coot::compare_atom_specs_user_float);
           std::string ls = coot::util::interesting_things_list(clash_atoms);
           cmd_strings.push_back(ls);
           std::string s = g.state_command(cmd_strings, coot::STATE_PYTHON);
           safe_python_command(s);
#endif // PYGTK
	 }
      }
   }
}
#endif // PYTHON

char *unmangle_hydrogen_name(const char *pdb_hydrogen_name) {

   std::string atom_name(pdb_hydrogen_name);
   std::string new_atom_name = atom_name;

   if (atom_name.length() == 4) { 
      if (atom_name[0] == '1' ||
	  atom_name[0] == '2' ||
	  atom_name[0] == '3' ||
	  atom_name[0] == '4' ||
	  atom_name[0] == '*') {
	 // switch it.
	 if (atom_name[3] == ' ') {
	    new_atom_name = " "; // lead with a space, testing.
	    new_atom_name += atom_name.substr(1,2) + atom_name[0];
	 } else { 
	    new_atom_name = atom_name.substr(1,3) + atom_name[0];
	 }
      }
   } else { 
      if (atom_name[3] != ' ') { 
	 if (atom_name[3] == ' ') {
	    new_atom_name = atom_name.substr(1,2) + atom_name[0];
	    new_atom_name += ' ';
	 }
	 if (atom_name[2] == ' ') {
	    new_atom_name = atom_name.substr(1,1) + atom_name[0];
	 new_atom_name += ' ';
	 new_atom_name += ' ';
	 }
      } else {
	 // atom_name length is 3 presumably
	 new_atom_name = ' ';
	 new_atom_name += atom_name.substr(1,2) + atom_name[0];
      }
   }
   
   char *s = new char[5];
   for (int i=0; i<5; i++) s[i] = 0;
   strncpy(s, new_atom_name.c_str(), 4);

//    std::cout << "mangle debug:: :" << pdb_hydrogen_name << ": to :" << s << ":" << std::endl;
   return s;
}



short int do_probe_dots_on_rotamers_and_chis_state() {
   return graphics_info_t::do_probe_dots_on_rotamers_and_chis_flag;
} 

void set_do_probe_dots_on_rotamers_and_chis(short int state) {
   graphics_info_t::do_probe_dots_on_rotamers_and_chis_flag = state;
}

void set_do_probe_dots_post_refine(short int state) {
   graphics_info_t::do_probe_dots_post_refine_flag = state;
} 

short int do_probe_dots_post_refine_state() {
   return graphics_info_t::do_probe_dots_post_refine_flag;
}




/*! \brief return the number of generic display objects */
int number_of_generic_objects() {

   graphics_info_t g;
   return g.generic_objects_p->size();
}


void generic_objects_gui_wrapper() {

   std::vector<std::string> cmd;
   cmd.push_back("generic-objects-gui");
   graphics_info_t g;

#ifdef USE_GUILE

   std::string s = g.state_command(cmd, coot::STATE_SCM);
   safe_scheme_command(s);


#else
#ifdef USE_PYGTK

   std::string s = g.state_command(cmd, coot::STATE_PYTHON);
   safe_python_command(s);

#endif // USE_PYGTK

#endif // USE_GUILE

} 

/*! \brief print to the console the name and display status of the
  generic display objects */
void generic_object_info() {

   graphics_info_t g;
   unsigned int n_obs = g.generic_objects_p->size();
   std::cout << "There are " << n_obs << " generic objects\n";
   if (n_obs) {
      for (int i=0; i<n_obs; i++) {
	 std::string display_str(":Displayed:");
	 if ((*g.generic_objects_p)[i].is_displayed_flag == 0)
	    display_str = ":Not Displayed:";
	 std::string closed_str(":Closed:");
	 if ((*g.generic_objects_p)[i].is_closed_flag == 0)
	    closed_str = ":Not Closed:";
	 std::cout << " # " << i << " \"" << (*g.generic_objects_p)[i].name << "\" "
		   << display_str << " " << closed_str << std::endl;
      }
   } else {
      std::cout << "No Generic Display Objects" << std::endl;
   } 
} 

// generic object obj_no has things to display?
// Return 0 or 1.
// 
short int generic_object_has_objects_p(int object_number) {

   short int r = 0;
   graphics_info_t g;
   if ((object_number >=0) && (object_number < g.generic_objects_p->size())) {
      if ((*g.generic_objects_p)[object_number].lines_set.size() > 0)
	 r = 1;
      if ((*g.generic_objects_p)[object_number].points_set.size() > 0)
	 r = 1;
   } else {
      std::cout << "WARNING:: object_number in generic_objects_p "
		<< object_number << std::endl;
   } 

   return r;

} 

std::pair<short int, std::string>
is_interesting_dots_object_next_p(const std::vector<std::string> &vs) {

   std::pair<short int, std::string> r(0, "");

   if (vs.size() == 3) {
//       std::cout << "Looking at bits:  \n  "; 
//       for (unsigned int i=0; i<3; i++) { 
// 	 std::cout << ":" << vs[i] << ": ";
//       }
//       std::cout << "\n"; 
      if ((vs[1] == "wide") && (vs[2] == "contact)")) {
	 r.first = 1;
	 r.second = "wide contact";
      }
      if ((vs[1] == "close") && (vs[2] == "contact)")) {
	 r.first = 1;
	 r.second = "close contact";
      }
      if ((vs[1] == "small") && (vs[2] == "overlap)")) {
	 r.first = 1;
	 r.second = "small overlap";
      }
      if ((vs[1] == "bad") && (vs[2] == "overlap)")) {
	 r.first = 1;
	 r.second = "bad overlap";
      }
      if (vs[1] == "H-bonds)") { 
	 r.first = 1;
	 r.second = "H-bonds";
      }
   }
   return r;
}

std::string probe_dots_short_contact_name_to_expanded_name(const std::string &short_name) {

   std::vector<std::pair<std::string, std::string> > names;
   names.push_back(std::pair<std::string, std::string>("wc", "wide contact"));
   names.push_back(std::pair<std::string, std::string>("cc", "close contact"));
   names.push_back(std::pair<std::string, std::string>("so", "small overlap"));
   names.push_back(std::pair<std::string, std::string>("bo", "bad overlap"));
   names.push_back(std::pair<std::string, std::string>("hb", "H-bonds"));

   std::string r = "unknown";
   for (int i=0; i<5; i++) {
      if (names[i].first == short_name) {
	 r = names[i].second;
	 break;
      }
   } 
   return r;
} 



/*! \brief close generic object, clear the lines/points etc, not
  available for buttons/displaying etc */
void close_generic_object(int object_number) {

   graphics_info_t g;
   if (object_number >=0) {
      if (object_number < int(g.generic_objects_p->size())) {
	 (*g.generic_objects_p)[object_number].close_yourself();
      }
   }
} 

/*! \brief has the generic object been closed? 

   @return 1 for yes, 0 othersize
*/
short int is_closed_generic_object_p(int object_number) {

   short int state = 0;
   graphics_info_t g;
   if (object_number >=0) { 
      if (object_number < int(g.generic_objects_p->size())) {
	 state = (*g.generic_objects_p)[object_number].is_closed_flag;
      }
   }
   return state;
} 


// This is tedious and irritating to parse in C++.
// 
// a const when a member function
// 
// Note that the filenames have a trailing "/".
// 
std::vector<std::pair<std::string, std::string> > 
parse_ccp4i_defs(const std::string &filename) {

   std::vector<std::pair<std::string, std::string> > v;

   // put the current directory in, whether or not we can find the
   // ccp4 project dir
   char *pwd = getenv("PWD");
   if (pwd) {
      v.push_back(std::pair<std::string, std::string> (std::string(" - Current Dir - "),
						       std::string(pwd) + "/"));
   }
   
   struct stat buf;
   int stat_status = stat(filename.c_str(), &buf);
   if (stat_status != 0) {
     // silently return nothing if we can't find the file.
     return v;
   } 

   std::ifstream cin(filename.c_str());

   // Let's also add ccp4_scratch to the list if the environment
   // variable is declared and if directory exists
   char *scratch = getenv("CCP4_SCR");
   if (scratch) {
      struct stat buf;
      int istat_scratch = stat(scratch, &buf);
      if (istat_scratch == 0) {
	 if (S_ISDIR(buf.st_mode)) {
	    v.push_back(std::pair<std::string, std::string>(std::string("CCP4_SCR"),
							    std::string(scratch) + "/"));
	 }
      }
   }

   if (! cin) {
      std::cout << "WARNING:: failed to open " << filename << std::endl;
   } else {
      // std::string s;
      char s[1000];
      std::vector <coot::alias_path_t> alias;
      std::vector <coot::alias_path_t> path;
      std::string::size_type ipath;
      std::string::size_type ialias;
      int index;
      int icomma;
      short int path_coming = 0;
      short int alias_coming = 0;
      bool alias_flag = 0;
      while (! cin.eof()) {
	 cin >> s;
	 std::string ss(s);
	 // std::cout << "parsing:" << ss << std::endl;
	 if (path_coming == 2) {
	    path.push_back(coot::alias_path_t(index, ss, alias_flag));
	    path_coming = 0;
	 }
	 if (alias_coming == 2) {
	    alias.push_back(coot::alias_path_t(index, ss, alias_flag));
	    alias_coming = 0;
	 }
	 if ( path_coming == 1)  path_coming++;
	 if (alias_coming == 1) alias_coming++;
	 ipath  = ss.find("PROJECT_PATH,");
	 ialias = ss.find("PROJECT_ALIAS,");
	 if (ipath != std::string::npos) {
	    // std::cout << "DEBUG::  found a project path..." << std::endl;
	    path_coming = 1;
	    alias_flag = 0;
	    icomma = ss.find_last_of(",");
	    // std::cout << icomma << " " << ss.length() << std::endl;
	    if ( (icomma+1) < int(ss.length())) {
	       index = atoi(ss.substr(icomma+1, ss.length()).c_str());
	       // std::cout << "index: " << index << std::endl;
	    }
	 }
	 if (ialias != std::string::npos) {
	    alias_coming = 1;
	    icomma = ss.find_last_of(",");
	    if ( (icomma+1) < int(ss.length())) {
	       index = atoi(ss.substr(icomma+1, ss.length()).c_str());
	    }
	 }

	 // Things called ALIASES at at the CCP4 top level are
	 // actually speciified by DEF_DIR_PATH and DEF_DIR_ALIAS 
	 // in the same way that PROJECT_ALIAS and PROJECT_PATH work.
	 //

	 ipath  = ss.find("DEF_DIR_PATH,");
	 ialias = ss.find("DEF_DIR_ALIAS,");
	 if (ipath != std::string::npos) {
	    // std::cout << "DEBUG::  found an ALIAS path..." << ss << std::endl;
	    path_coming = 1;
	    alias_flag = 1;
	    icomma = ss.find_last_of(",");
	    // std::cout << icomma << " " << ss.length() << std::endl;
	    if ( (icomma+1) < int(ss.length())) {
	       index = atoi(ss.substr(icomma+1, ss.length()).c_str());
	       // std::cout << "index: " << index << std::endl;
	    }
	 }
	 if (ialias != std::string::npos) {
	    alias_coming = 1;
	    // std::cout << "DEBUG::  found an ALIAS name..." << ss << std::endl;
	    icomma = ss.find_last_of(",");
	    if ( (icomma+1) < int(ss.length())) {
	       index = atoi(ss.substr(icomma+1, ss.length()).c_str());
	    }
	 }
	 
      }

//       std::cout << "----------- path pairs: ------------" << std::endl;
//       for (int i=0; i<path.size(); i++) {
//   	 std::cout << path[i].index << "  " << path[i].s << " " << path[i].flag << std::endl;
//       }
//       std::cout << "----------- alias pairs: ------------" << std::endl;
//       for (int i=0; i<alias.size(); i++)
//   	 std::cout << alias[i].index << "  " << alias[i].s << " " << alias[i].flag << std::endl;
//       std::cout << "-------------------------------------" << std::endl;
      

      std::string alias_str;
      std::string path_str;
      for (unsigned int j=0; j<alias.size(); j++) {
	 for (unsigned int i=0; i<path.size(); i++) {
	    if (path[i].index == alias[j].index) {
	       if (path[i].flag == alias[j].flag) {
		  // check for "" "" pair here.
		  alias_str = alias[j].s;
		  path_str  = path[i].s;
		  // if the file is a directory, we need to put a "/" at the
		  // end so that went we set that filename in the fileselection
		  // widget, we go into the directory, rather than being in the
		  // directory above with the tail as the selected file.
		  //
		  struct stat buf;
		  int status = stat(path_str.c_str(), &buf);
	       
		  // valgrind says that buf.st_mode is uninitialised here
		  // strangely.  Perhaps we should first test for status?
		  // Yes - that was it.  I was using S_ISDIR() on a file
		  // that didn't exist.  Now we skip if the file does not
		  // exist or is not a directory.

		  // std::cout << "stating "<< path_str << std::endl;

		  if (status == 0) { 
		     if (S_ISDIR(buf.st_mode)) {
			path_str += "/";

			if (alias_str == "\"\"") {
			   alias_str = "";
			   path_str  = "";
			}
			v.push_back(std::pair<std::string, std::string> (alias_str, path_str));
		     }
		     // } else { 
		     // // This is too boring to see every time we open a file selection
		     // std::cout << "INFO:: directory for a CCP4i project: " 
		     // << path_str << " was not found\n";
		  }
	       }
	    }
	 }
      }
   }
   return v;
}

std::string
ccp4_project_directory(const std::string &ccp4_project_name) {

   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();
   std::vector<std::pair<std::string, std::string> > v = 
      parse_ccp4i_defs(ccp4_defs_file_name);
   std::string r = "";
   for (unsigned int i=0; i<v.size(); i++) {
      if (v[i].first == ccp4_project_name) {
	 r = v[i].second;
	 break;
      }
   }
   return r;
}


/* movies */
void set_movie_file_name_prefix(const char *file_name) {
   graphics_info_t::movie_file_prefix = file_name;
}

void set_movie_frame_number(int frame_number) {
   graphics_info_t::movie_frame_number = frame_number;
}

#ifdef USE_GUILE
SCM movie_file_name_prefix() {
   SCM r = scm_makfrom0str(graphics_info_t::movie_file_prefix.c_str());
   return r;
}
#endif
#ifdef USE_PYTHON
PyObject * movie_file_name_prefix_py() {
   PyObject *r;
   r = PyString_FromString(graphics_info_t::movie_file_prefix.c_str());
   return r;
}
#endif // PYTHON

int movie_frame_number() {
   return graphics_info_t::movie_frame_number;
}

void set_make_movie_mode(int make_movie_flag) {
   graphics_info_t::make_movie_flag = make_movie_flag;
}


#ifdef USE_GUILE
void try_load_scheme_extras_dir() {

   char *s = getenv("COOT_SCHEME_EXTRAS_DIR");
   if (s) {
      struct stat buf;
      int status = stat(s, &buf);
      if (status != 0) {
	 std::cout << "WARNING:: no directory " << s << std::endl;
      } else {
	 if (S_ISDIR(buf.st_mode)) {

	    DIR *lib_dir = opendir(s);
	    if (lib_dir == NULL) {
	       std::cout << "An ERROR occured on opening the directory "
			 << s << std::endl;
	    } else {

	       struct dirent *dir_ent;

	       // loop until the end of the filelist (readdir returns NULL)
	       // 
	       while (1) {
		  dir_ent = readdir(lib_dir);
		  if (dir_ent == NULL) {
		     break;
		  } else {
		     std::string sub_part(std::string(dir_ent->d_name));
		     struct stat buf2;
		     std::string fp = s;
		     fp += "/";
		     fp += sub_part;
		     int status2 = stat(fp.c_str(), &buf2);
		     if (status2 != 0) {
			std::cout << "WARNING:: no file " << sub_part << std::endl;
		     } else {
			if (S_ISREG(buf2.st_mode)) {
			   if (coot::util::file_name_extension(sub_part) == ".scm") {
			      std::cout << "loading extra: " << fp << std::endl;
			      scm_c_primitive_load(fp.c_str()); 
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
#endif // USE_GUILE

void set_button_label_for_external_refinement(const char *button_label) {
   graphics_info_t::external_refinement_program_button_label = button_label;
}
