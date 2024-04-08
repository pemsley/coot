/*
 * cootilus/nautilus-join.h
 *
 * Copyright 2011 by Kevin Cowtan
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


#ifndef NAUTILUS_JOIN_H
#define NAUTILUS_JOIN_H


#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


class NucleicAcidJoin {
 public:
  static clipper::MiniMol join( const clipper::MiniMol& mol );
 private:
  class Node { public: float score; std::vector<int> ptrs; };
  static std::vector<int> best_chain( std::vector<Node>& nodes );
};


#endif
