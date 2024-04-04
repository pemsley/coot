/*
 * src/data-pair-remover.hh
 *
 * Copyright 2015 by Medical Research Council
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


class data_pair_remover {
   public:
   double r_lim;
   data_pair_remover(const double &r_in) { r_lim = r_in; }
   bool operator()(const std::pair<double,double> &rp) const {
      double r1 = coot::util::random()/float(RAND_MAX);
      // std::cout << "   " << r1 << " " << r_lim << "\n";
      return (r1 > r_lim);
   }
};
