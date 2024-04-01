/*
 * coot-utils/exp-fit.hh
 *
 * Copyright 2022 by Medical Research Council
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
#ifndef EXP_FIT_HH
#define EXP_FIT_HH

#include <cmath>
#include <vector>

namespace coot {

   class exponential_fit_with_offset {
   public:

      // for y = a + bexp(cx)
      double a; // the offset
      double b; // scale fact
      double c; // scale for x

      // A_data must be in increasing order of x (up to caller to sort that out)
      explicit exponential_fit_with_offset(const std::vector<std::pair<double, double> > &A_data);
      double at(const double &x) const {
         double v = a + b * std::exp(c * x);
         return v;
      }
      double average_deviation(const std::vector<std::pair<double, double> > &A_data) const {
         if (! A_data.empty()) {
            double sum_delta = 0.0;
            for (unsigned int i=0; i<A_data.size(); i++) {
               const double &x = A_data[i].first;
               const double &y = A_data[i].second;
               double y_hat = at(x);
               double delta = std::fabs(y_hat - y);
               sum_delta += delta;
            }
            double av_delta = sum_delta/static_cast<double>(A_data.size());
            return av_delta;
         } else {
            return 0.0;
         }
      }
   };
}

#endif // EXP_FIT_HH
