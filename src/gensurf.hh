/*
 * src/gensurf.hh
 *
 * Copyright 2019 by Medical Research Council
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

#ifndef COOT_SRC_GENSURF_HH
#define COOT_SRC_GENSURF_HH

#include <clipper/core/xmap.h>
#include "coords/Cartesian.hh"
#include "density-contour/density-contour-triangles.hh"

void gensurf_and_add_vecs_threaded_workpackage(const clipper::Xmap<float> *xmap_p,
					       float contour_level, float dy_radius,
					       coot::Cartesian centre,
					       int iream_start, int n_reams,
					       int isample_step,
					       bool is_em_map,
                                               bool use_vertex_gradients_for_map_normals_flag,
					       std::vector<coot::density_contour_triangles_container_t> *draw_vector_sets_p);

// std::vector<std::pair<const coot::CartesianPair *,int> > draw_vector_sets;


#endif // COOT_SRC_GENSURF_HH
