
#include <stdexcept>
#include <iostream>

#include <gtk/gtk.h>
#include "gtk-widget-conversion-utils.h"

#include "utils/coot-utils.hh"

struct entry_info_t coot_entry_to_val(GtkEntry *entry) { 

  struct entry_info_t ei;
  const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
  ei.float_is_set = 0;
  ei.val_as_float = 0;
  ei.val = 0;

  if (text) {
     ei.string_is_set = 1;
     ei.string = text;
     try {
	ei.val_as_float = coot::util::string_to_float(text);
     }
     catch (const std::runtime_error &rte) {
	std::cout << rte.what() << std::endl;
     } 
  } 

  return ei;

}

