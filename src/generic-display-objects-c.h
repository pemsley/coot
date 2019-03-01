
#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS extern
#define END_C_DECLS     
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS


/*! \file
  \brief Coot Scripting Interface - Generic objects
*/

void clear_generic_objects_dialog_pointer();

/*! \brief display (1) or undisplay (0) all generic display objects */
void set_display_all_generic_objects(int state);

/*! \brief a kludgey thing, so that the generic objects gui can be
  called from a callback.  */
void generic_objects_gui_wrapper();

/*! \brief close all generic display objects  */
void close_all_generic_objects();

END_C_DECLS
