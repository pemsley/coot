/*
 * coot-utils/shiftfield.h
 *
 * Copyright 2018 Kevin Cowtan & University of York all rights reserved
 * Author: Kevin Cowtan
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


#ifndef SHIFTFIELD
#define SHIFTFIELD

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

class Shift_field_refine{

 public:

  static bool shift_field_coord(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
    clipper::Xmap<float>& x1map,      clipper::Xmap<float>& x2map,      clipper::Xmap<float>& x3map,
    float rad, int filter );
  static bool shift_field_coord_const(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
    clipper::Xmap<float>& x1map,      clipper::Xmap<float>& x2map,      clipper::Xmap<float>& x3map,
    float rad, int filter );

  static bool shift_field_u_iso(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
    clipper::Xmap<float>& x1map,
    float rad, int filter );
  static bool shift_field_u_iso_const(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
    clipper::Xmap<float>& x1map,
    float rad, int filter );

  static bool shift_field_u_aniso(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
		clipper::Xmap<float>& x1map, clipper::Xmap<float>& x2map, clipper::Xmap<float>& x3map,
    clipper::Xmap<float>& x4map, clipper::Xmap<float>& x5map, clipper::Xmap<float>& x6map, 
    float rad, clipper::Resolution& reso, int filter );
  static bool shift_field_u_aniso_const(
    const clipper::Xmap<float>& cmap, const clipper::Xmap<float>& dmap, const clipper::Xmap<float>& mask,
		clipper::Xmap<float>& x1map, clipper::Xmap<float>& x2map, clipper::Xmap<float>& x3map,
    clipper::Xmap<float>& x4map, clipper::Xmap<float>& x5map, clipper::Xmap<float>& x6map, 
    float rad, clipper::Resolution& reso, int filter );

};

#endif

