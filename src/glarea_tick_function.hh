#ifndef GLAREA_TICK_FUNCTION_HH
#define GLAREA_TICK_FUNCTION_HH

// also in draw-2.hh

#include <gtk/gtk.h>

gboolean
glarea_tick_func(GtkWidget *widget,
                 GdkFrameClock *frame_clock,
                 gpointer data);

#endif
