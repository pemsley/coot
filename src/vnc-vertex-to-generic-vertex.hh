/*
 * src/vnc-vertex-to-generic-vertex.hh
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
#ifndef VNC_VERTEX_TO_GENERIC_VERTEX_HH
#define VNC_VERTEX_TO_GENERIC_VERTEX_HH

#include "coot-utils/vertex.hh"
#include "generic-vertex.hh"

s_generic_vertex vnc_vertex_to_generic_vertex(const coot::api::vnc_vertex &v);
s_generic_vertex vn_vertex_to_generic_vertex(const coot::api::vn_vertex &v);

#endif // VNC_VERTEX_TO_GENERIC_VERTEX_HH
