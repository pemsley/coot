#ifndef EVENT_CONTROLLER_CALLBACKS_HH
#define EVENT_CONTROLLER_CALLBACKS_HH

#include "gtk/gtk.h"

// ------------------- 20220429-PE lovely new controller code ---------------------

void on_glarea_drag_begin(GtkGestureDrag *gesture, double x, double y, GtkWidget *area);
void on_glarea_drag_update(GtkGestureDrag *gesture, double delta_x, double delta_y, GtkWidget *area);
void on_glarea_drag_end(GtkGestureDrag *gesture, double x, double y, GtkWidget *area);
gboolean
on_key_controller_key_pressed(GtkEventControllerKey *controller,
                              guint                  keyval,
                              guint                  keycode,
                              guint                  modifiers,
                              GtkButton             *button);
void
on_key_controller_key_released(GtkEventControllerKey *controller,
                               guint                  keyval,
                               guint                  keycode,
                               guint                  modifiers,
                               GtkButton             *button);



#endif // EVENT_CONTROLLER_CALLBACKS_HH
