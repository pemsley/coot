/*
 * analysis/kolmogorov.hh
 *
 * Copyright 2012 by The Medical Research Council
 * Author: Robert Nicholls
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

#ifndef INCLUDE_KOLMOGOROV_HH
#define INCLUDE_KOLMOGOROV_HH

#include <vector>

namespace nicholls {

   double get_KS(const std::vector<double> &v1, const std::vector<double> &v2);

   // Assumption - data are positive.
   // Kullback-Liebler 
   std::pair<double, double> get_KL(const std::vector<double> &v1, const std::vector<double> &v2);

}

#endif // INCLUDE_KOLMOGOROV_HH
