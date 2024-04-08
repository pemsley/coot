/*
 * cootilus/cootilus-build.h
 *
 * Copyright 2012 by Kevin Cowtan
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


class Coot_nucleic_acid_build {
 public:
  //! Constructor: takes filename for nucleic acid library file
  Coot_nucleic_acid_build( std::string filename );
  //! Build or extend a model
  bool build( mmdb::Manager* mmdb, const clipper::Xmap<float>& xmap, const clipper::Coord_orth& centre, double radius ) const;
 private:
  std::string filename_;
};

