/*
 * src/event-controller-callbacks.hh
 *
 * Copyright 2022 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
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
