/* src/c-interface.cc
 * 
 * Copyright 2013 by Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

 
// Load the head if it hasn't been included.
#ifdef USE_PYTHON
#ifndef PYTHONH
#define PYTHONH
#include <Python.h>
#endif
#endif

#include "graphics-info.h"

#include "c-interface.h"


/*  ----------------------------------------------------------------------- */
/*                  ribbons and spheroids display                           */
/*  ----------------------------------------------------------------------- */

void
set_show_ribbons(int imol, short int state) {

  graphics_info_t::molecules[imol].cootribbons = state;
  if (state == 1) {
    graphics_info_t::molecules[imol].make_ribbons();
  }
  graphics_draw();

}

void
set_show_aniso_spheroids(int imol, short int state) {

  graphics_info_t::molecules[imol].cootanisospheroids = state;
  if (state == 1) {
    graphics_info_t::molecules[imol].make_aniso_spheroids();
  }
  graphics_draw();

}


void
set_ribbon_param_int(int imol, const char *name, int value) {

  graphics_info_t::molecules[imol].set_ribbon_param(name, value);
}

void
set_ribbon_param_float(int imol, const char *name, float value) {

  graphics_info_t::molecules[imol].set_ribbon_param(name, value);
}

void
set_global_ccp4mg_param_int(int imol, const char *name, int value) {

  graphics_info_t::molecules[imol].set_global_ccp4mg_param(name, value);
}

void
set_global_ccp4mg_param_float(int imol, const char *name, float value) {

  graphics_info_t::molecules[imol].set_global_ccp4mg_param(name, value);
}
