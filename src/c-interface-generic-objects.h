/* src/c-interface-generic-objects/.cc
 * 
 * Copyright 2011 by the University of Oxford
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#include <clipper/core/coords.h>
#include <string>
#include "utils/colour-holder.hh"
#include "generic-display-objects-c.h"


/*! \file
  \brief Coot Scripting Interface - Generic Objects interface
*/


/*  ----------------------------------------------------------------------- */
/*                  Generic Objects                                         */
/*  ----------------------------------------------------------------------- */
/*! \name Generic Objects */
/* \{ */

/*! \brief create a new generic object with name objname

   return the index of the object */
int new_generic_object_number(const std::string &obj_name);

/*! \brief create a new generic object with name objname
           and attach it to the given molecule

bool is_valid_generic_display_object_number(int obj);
  @return the index of the object */
int new_generic_object_number_for_molecule(const std::string &obj_name, int imol);

/*! \brief add line to generic object object_number */
void to_generic_object_add_line(int object_number,
                                const char *colour,
                                int line_width,
                                float to_x1,
                                float to_y1,
                                float to_z1,
                                float to_x2,
                                float to_y2,
                                float to_z2);

#ifdef USE_PYTHON
/*! \brief add multiple lines to generic object object_number

c.f. to_generic_object_add_points()
each list item is a [colour, width, x1, y1, z1, x2, y2, z2] */
void to_generic_object_add_lines(int object_number, PyObject *line_info_list_py);
#endif

void to_generic_object_add_cylinder(int object_number,
                                    const char *colour,
                                    float line_radius,
                                    int n_slices, // 4, 8, 16
                                    float from_x,
                                    float from_y,
                                    float from_z,
                                    float to_x,
                                    float to_y,
                                    float to_z,
                                    bool cap_start,
                                    bool cap_end);

/*! \brief add a dashed line to generic object object_number

dash_density is number of dashes per Angstrom.*/
void to_generic_object_add_dashed_line(int object_number,
				       const char *colour,
				       int line_width,
				       float dash_density,
				       float from_x1,
				       float from_y1,
				       float from_z1,
				       float to_x2,
				       float to_y2,
				       float to_z2);

/*! \brief add point to generic object object_number */
void to_generic_object_add_point(int object_number,
				 const char *colour,
				 int point_width,
				 float from_x1,
				 float from_y1,
				 float from_z1);

#ifdef USE_PYTHON
// point_info_list_py is a list of [colour, point_width, x, y, z]
void to_generic_object_add_points(int object_number, PyObject *point_info_list_py);
#endif


#ifndef SWIG
void to_generic_object_add_point_internal(int object_number,
				 const std::string &colour_name, // needed for indexing objects by colour
				 const coot::colour_holder &colour,
				 int point_width,
				 const clipper::Coord_orth &pt);
#endif // SWIG

void from_generic_object_remove_last_item(int object_number);

/*! \brief add point to generic object object_number */
void to_generic_object_add_arc(int object_number,
			       const char *colour,
			       float radius,
			       float radius_inner,
			       float angle_delta,
			       float start_point_x,
			       float start_point_y,
			       float start_point_z,
			       float start_dir_x,
			       float start_dir_y,
			       float start_dir_z,
			       float normal_x1,
			       float normal_y1,
			       float normal_z1);

void to_generic_object_add_torus(int object_number,
                                 const char *colour_name,
                                 float radius,
                                 float radius_inner,
                                 float centre_point_x,
                                 float centre_point_y,
                                 float centre_point_z,
                                 float normal_x,
                                 float normal_y,
                                 float normal_z);

void
to_generic_object_add_arrow(int object_number,
                            const char *colour_name,
                            float stem_radius,
                            float from_x1,
                            float from_y1,
                            float from_z1,
                            float to_x2,
                            float to_y2,
                            float to_z2);

void to_generic_object_add_dodecahedron(int object_number,
					const char *colour,
					float radius,
					float x,
					float y,
					float z);

void to_generic_object_add_pentakis_dodecahedron(int object_number,
						 const char *colour,
						 float stellation_factor,
						 float radius_factor,
						 float x,
						 float y,
						 float z);

#ifdef USE_PYTHON
void to_generic_object_add_mesh(int object_number, PyObject *mesh_py);
#endif

void to_generic_object_attach_translation_gizmo(int object_number);

void generic_object_mesh_calculate_normals(int object_number);

/*! \brief add a display list handle generic object */
void to_generic_object_add_display_list_handle(int object_number, int display_list_id); 

/*! \brief set the display status of object number object_number, 

  when they are created, by default objects are not displayed, so we
  generally need this function.  */
void set_display_generic_object(int object_number, short int istate);

/*! \brief set the display status of object number object_number, 

  set the state of a generic object to be drawn, but no redraw. 
  Use when enabling multiple generic objects.
*/
void set_display_generic_object_simple(int object_number, short int istate);


/*! \brief is generic display object displayed?

  @return 1 for yes, otherwise 0  */
int generic_object_is_displayed_p(int object_number);

/*! \brief return the index of the object with name name, if not, return -1; */
int generic_object_index(const std::string &name);

/*! \brief what is the name of generic object number obj_number? 

 @return 0 (NULL) (scheme False)  on obj_number not available */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM generic_object_name_scm(int obj_number);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *generic_object_name_py(unsigned int obj_number);
PyObject *get_generic_object_py(unsigned int obj_number);
#endif /* USE_PYTHON */
#endif /*  __cplusplus */

/*! \brief return the number of generic display objects */
int number_of_generic_objects();

/*! \brief print to the console the name and display status of the
  generic display objects */
void generic_object_info(); 

#ifdef USE_PYTHON
/*! \brief get generic display objects */
PyObject *get_generic_object_info(int obj_number);
#endif /* USE_PYTHON */

/*! \brief does generic display object number obj_no have things to
  display? (predicate name)

@return 0 for no things, 1 for things. */
short int generic_object_has_objects_p(int obj_no); 

/*! \brief close generic object
 *
 * clear the lines/points etc, not available for buttons/displaying etc
 *
 * param obj_number the object_number of the generic object to close
 * */
void close_generic_object(int object_number);

/*! \brief has the generic object been closed? 

   @return 1 for yes, 0 othersize
*/
short int is_closed_generic_object_p(int object_number);

/*! \brief clear out the lines and points from object_number, but keep
  it displayable (not closed). */
void generic_object_clear(int object_number);

/*! \brief attach the generic object to a particular molecule 

one might do this if the generic object is specific to a molecule.
 */
void attach_generic_object_to_molecule(int obj_number, int imol);

// This no longer maeks sense.
// void set_display_generic_objects_as_solid(int state);


/* \} */
