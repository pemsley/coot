#ifndef COOT_APPLICATION_HH
#define COOT_APPLICATION_HH

#include <gtk/gtk.h>

void
application_activate(GtkApplication *app,
                     gpointer        user_data);

// return status
int start_using_application(int argc, char **argv);


#endif // COOT_APPLICATION_HH
