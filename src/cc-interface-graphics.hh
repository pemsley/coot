/*
 * src/cc-interface-graphics.hh
 *
 * Copyright 2023 by Medical Research Council
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
#ifndef CC_INTERFACE_GRAPHICS_HH
#define CC_INTERFACE_GRAPHICS_HH

#include <string>

void show_coot_points_frame();

/*! \brief Set the view
 *
 * the view is a JSON string of a dict of the orientation quaternion, the zoom and the rotation centre.
 *
 * @param view_as_json the view parameter as a JSON string
 * */
void set_view_from_json(const std::string &view_as_json);

/*! \brief get the view
 *
 * the view is a JSON string of a dict of the orientation quaternion, the zoom and the rotation centre.
 *
 * */
std::string get_view_as_json();

#endif // CC_INTERFACE_GRAPHICS_HH
