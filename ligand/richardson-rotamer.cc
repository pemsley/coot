/*
 * ligand/richardson-rotamer.cc
 *
 * Copyright 2007 by University of York
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

#include "richardson-rotamer.hh"

std::vector<float>
coot::richardson_rotamer::probabilities() const {

   std::string rt = Residue_Type();
   if (rt == "MSE")
      rt = "MET";

   std::vector<simple_rotamer> rots = get_rotamers(rt, Probability_limit());

   std::vector<float> p(rots.size());
   for(unsigned int i=0; i<rots.size(); i++)
      p[i] = rots[i].Probability_rich();
   
   return p;
}
