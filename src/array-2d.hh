/*
 * src/array-2d.hh
 *
 * Copyright 2020 by Medical Research Council
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

#ifndef ARRAY_2D_HH
#define ARRAY_2D_HH


class array_2d {
public:
   std::vector<float> v;
   unsigned int nx;
   array_2d(unsigned int nx, unsigned int ny) : nx(nx) {
      v.resize(nx * ny);
   }
   void set(const unsigned int &i, const unsigned int &j, const float &f) {
      v[j * nx + i] = f;
   }
   float get(const unsigned int &i, const unsigned int &j) const {
      return v[j * nx + i];
   }

};

class coord_array_2d {
public:
   std::vector<std::pair<clipper::Coord_orth, float> > v;
   unsigned int nx;
   coord_array_2d(unsigned int nx, unsigned int ny) : nx(nx) {
      v.resize(nx * ny);
   }
   void set(const unsigned int &i, const unsigned int &j, const clipper::Coord_orth &co, const float &f) {
      v[j * nx + i] = std::pair<clipper::Coord_orth, float> (co, f);
   }
   std::pair<clipper::Coord_orth, float>  get(const unsigned int &i, const unsigned int &j) const {
      return v[j * nx + i];
   }
   clipper::Coord_orth get_co(const unsigned int &i, const unsigned int &j) const {
      return v[j * nx + i].first;
   }
   float get_f(const unsigned int &i, const unsigned int &j) const {
      return v[j * nx + i].second;
   }

};

#endif // ARRAY_2D

