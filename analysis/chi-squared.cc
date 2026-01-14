/* analysis/chi-squared.cc
 *
 * Copyright 2025,  by Global Phasing Ltd.
 * Author: ClAuS Flensburg, Clemens Vonrhein, Gerard Bricogne
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

#include "chi-squared.hh"

float gphl::prob_to_radius(float prob) {
   // Return the scale factor for the confidence level at "prob" percent
   // of a Chi^2 distribution with three degrees of freedom.
   // Note that the 50% level correspond to a radius of 1.5382

  static const float confidence_scale[99] = {
            0.3389, 0.4299, 0.4951, 0.5478, 0.5931, 0.6334, 0.6699, 0.7035, 0.7349,
    0.7645, 0.7925, 0.8191, 0.8447, 0.8694, 0.8932, 0.9162, 0.9387, 0.9605, 0.9818,
    1.0026, 1.0230, 1.0430, 1.0627, 1.0821, 1.1012, 1.1200, 1.1386, 1.1570, 1.1751,
    1.1932, 1.2110, 1.2288, 1.2464, 1.2639, 1.2812, 1.2985, 1.3158, 1.3330, 1.3501,
    1.3672, 1.3842, 1.4013, 1.4183, 1.4353, 1.4524, 1.4695, 1.4866, 1.5037, 1.5209,
    1.5382, 1.5555, 1.5729, 1.5904, 1.6080, 1.6257, 1.6436, 1.6615, 1.6797, 1.6980,
    1.7164, 1.7351, 1.7539, 1.7731, 1.7923, 1.8119, 1.8318, 1.8519, 1.8724, 1.8932,
    1.9144, 1.9360, 1.9580, 1.9804, 2.0034, 2.0269, 2.0510, 2.0757, 2.1012, 2.1274,
    2.1544, 2.1824, 2.2114, 2.2416, 2.2730, 2.3059, 2.3404, 2.3767, 2.4152, 2.4563,
    2.5003, 2.5478, 2.5997, 2.6571, 2.7216, 2.7955, 2.8829, 2.9912, 3.1364, 3.3682
  };
  int iprob = int(prob + 0.5);
  if (iprob <  1) iprob =  1;
  if (iprob > 99) iprob = 99;
  return confidence_scale[iprob-1];
}
