/*
 * src/named-rotamer-score.hh
 *
 * Copyright 2012 by University of York
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
#ifndef NAMED_ROTAMER_SCORE_HH
#define NAMED_ROTAMER_SCORE_HH

#include <string>
#include <vector>

namespace coot {

   class named_rotamer_score {
   public:
      std::string name;
      float clash_score;
      float rotamer_probability_score;
      std::vector<std::pair<std::string, float> > density_score_for_atoms;
      float density_fit_score;
      named_rotamer_score(const std::string &name_in,
			  float rotamer_probability_score_in,
			  float clash_score_in,
			  const std::vector<std::pair<std::string, float> > &atom_density_in,
			  float density_fit_score_in) {
	 name = name_in;
	 clash_score = clash_score_in;
	 density_fit_score = density_fit_score_in;
	 rotamer_probability_score = rotamer_probability_score_in;
	 density_score_for_atoms = atom_density_in;
      } 
   };
} 

#endif // NAMED_ROTAMER_SCORE_HH
