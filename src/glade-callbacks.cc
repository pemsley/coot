
#include "Python.h"

#include <iostream>
#include <gtk/gtk.h>

#include "c-interface.h"

extern "C" G_MODULE_EXPORT
void tt_fn(GtkMenuItem *menuitem, gpointer user_data) {
   std::cout << "Load tutorial model and data" << std::endl;
   load_tutorial_model_and_data();
}


extern "C" G_MODULE_EXPORT
void on___glade_unnamed_94_activate(GtkMenuItem *menuitem, gpointer user_data) {
   std::cout << "Load tutorial model and data" << std::endl;
   load_tutorial_model_and_data();
}
