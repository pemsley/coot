/* ligand/ligand.cc 
 * 
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Author Paul Emsley
 * Copyright 2008, 2009 by The University of Oxford
 * Copyright 2012, 2015, 2016 by Medical Research Council
 * Author Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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
 * 02110-1301, USA.
 */


#include "ligand-score-card.hh"

// void
// coot::ligand_score_card::make_cached_scores(bool score_for_pro) {
//    double s = get_score(score_for_pro);
//    cached_score = std::pair<bool, double> (true, s);
// }


// return the score.
// score_for_pro is a optional argument default false
double
coot::ligand_score_card::get_score() const {

   // maybe use a cache for the scoring?

   return atom_point_score;

}
