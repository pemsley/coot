
#ifndef GTK_WIDGET_CONVERSION_UTILS_H
#define GTK_WIDGET_CONVERSION_UTILS_H

#include <gtk/gtk.h>

#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS
#define END_C_DECLS     
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS

struct entry_info_t {
  short int float_is_set;
  short int string_is_set;
  int val;
  float val_as_float;
  const char *string;
}; 

struct entry_info_t coot_entry_to_val(GtkEntry *entry);

END_C_DECLS

#endif /* GTK_WIDGET_CONVERSION_UTILS_H */
