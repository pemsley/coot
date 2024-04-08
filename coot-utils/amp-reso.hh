/*
 * coot-utils/amp-reso.hh
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

#ifndef AMP_RESO_HH
#define AMP_RESO_HH

namespace coot {

   class amplitude_vs_resolution_point {
   public:
      amplitude_vs_resolution_point() {
	 sum = 0.0; count = 0; resolution_recip_sum = 0.0; finished = false; }
      double sum;
      double average;
      unsigned int count;
      double resolution_recip_sum;
      double resolution_recip;
      bool finished;
      void add(const float &f_in, const float &inv_res_sq_in) {
	 sum += f_in;
	 resolution_recip_sum += inv_res_sq_in;
	 count += 1;
      }
      void finish() {
	 if (count > 0) {
	    average = sum/static_cast<float>(count);
	    resolution_recip = resolution_recip_sum/static_cast<float>(count);
	 }
	 finished = true;
      }
      double get_average_fsqrd() const {
	 if (finished) {
	    return average;
	 } else {
	    std::cout << "amplitude_vs_resolution_point() Not finihsed " << std::endl;
	    return 0.0;
	 }
      }
      double get_invresolsq() const {
	 if (finished) {
	    return resolution_recip;
	 } else {
	    std::cout << "amplitude_vs_resolution_point() Not finihsed " << std::endl;
	    return 0.0;
	 }
      }
   };
}

#endif // AMP_RESO_HH

