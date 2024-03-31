/*
 * analysis/b-factor-histogram.hh
 *
 * Copyright 2009 by University of Oxford
 * Author: Paul Emsley
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


#include <vector>
#include <mmdb2/mmdb_manager.h>

namespace coot {

   class b_factor_histogram {
      int n_bins;
      int n_atoms;
      float b_max;
      void init();
      int get_n_atoms() const;
      int get_n_bins() const;
      int b_to_bin(const float &b) const;
      std::vector<std::vector<float> > b_vector;
      float alpha_estimate;
      float beta_estimate;
      float shift_estimate;
      double ig(const double &x) const;
      double Gamma(const double &b) const;
      void optimize_estimates();
      std::vector<double> select_from_model() const;
   public:
      b_factor_histogram(mmdb::Manager *mol);
      b_factor_histogram(mmdb::Manager *mol, int atom_selection_handle);
      void model();
      std::vector<std::pair<double, double> > get_data() const;
      std::vector<std::pair<double, double> > get_model() const; // call model() first
   };
}

