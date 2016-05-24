
/*  ----------------------------------------------------------------------- */
/*                  Generic Objects                                         */
/*  ----------------------------------------------------------------------- */
/*! \name Generic Objects */
/* \{ */

/*! \brief create a new generic object with name objname and return the index 
   of the object */
int new_generic_object_number(const char *objname);

/*! \brief add line to generic object object_number */
void to_generic_object_add_line(int object_number, 
				const char *colour,
				int line_width,
				float from_x1, 
				float from_y1, 
				float from_z1, 
				float to_x2, 
				float to_y2, 
				float to_z2);

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

/*! \brief add point to generic object object_number */
void to_generic_object_add_arc(int object_number, 
			       const char *colour,
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
			       float normal_x1, 
			       float normal_y1, 
			       float normal_z1);



/*! \brief add a display list handle generic object */
void to_generic_object_add_display_list_handle(int object_number, int display_list_id); 

/*! \brief set the display status of object number object_number, 

  when they are created, by default objects are not displayed, so we
  generally need this function.  */
void set_display_generic_object(int object_number, short int istate);

/*! \brief display (1) or undisplay (0) all generic display objects */
void set_display_all_generic_objects(int state);


/*! \brief is generic display object displayed?

  @return 1 for yes, otherwise 0  */
int generic_object_is_displayed_p(int object_number);

/*! \brief return the index of the object with name name, if not, return -1; */
int generic_object_index(const char *name);

/*! \brief what is the name of generic object number obj_number? 

 @return 0 (NULL) (scheme False)  on obj_number not available */
#ifdef __cplusplus
#ifdef USE_GUILE
SCM generic_object_name_scm(int obj_number);
#endif /* USE_GUILE */
#ifdef USE_PYTHON
PyObject *generic_object_name_py(int obj_number);
#endif /* USE_PYTHON */
#endif /*  __cplusplus */

/*! \brief return the number of generic display objects */
int number_of_generic_objects();

/*! \brief print to the console the name and display status of the
  generic display objects */
void generic_object_info(); 

/*! \brief does generic display object number obj_no have things to
  display? (predicate name)

@return 0 for no things, 1 for things. */
short int generic_object_has_objects_p(int obj_no); 

/*! \brief close generic object, clear the lines/points etc, not
  available for buttons/displaying etc */
void close_generic_object(int object_number);

/*! \brief has the generic object been closed? 

   @return 1 for yes, 0 othersize
*/
short int is_closed_generic_object_p(int object_number);

void close_all_generic_objects();

/*! \brief clear out the lines and points from object_number, but keep
  it displayable (not closed). */
void generic_object_clear(int object_number);

/*! \brief a kludgey thing, so that the generic objects gui can be
  called from a callback.  */
void generic_objects_gui_wrapper();

/*! \brief attach the generic object to a particular molecule 

one might do this if the generic object is specific to a molecule.
 */
void attach_generic_object_to_molecule(int obj_number, int imol);

void set_display_generic_objects_as_solid(int state);


/* \} */
