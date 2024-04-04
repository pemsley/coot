/*
 * src/ligand-check.hh
 *
 * Copyright 2014 by Medical Research Council
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

#ifndef LIGAND_CHECK_HH
#define LIGAND_CHECK_HH

#include "probe-clash-score.hh"

namespace coot { 

   class ligand_report_percentiles_t {
   public:
      // double mogul_percentile;
      double dictionary_geometry_percentile;
      double density_correlation_percentile;
      double probe_clash_percentile;
      ligand_report_percentiles_t(double a, double b, double c) {
	 density_correlation_percentile = a;
	 // mogul_percentile = b;
         dictionary_geometry_percentile = b;
	 probe_clash_percentile = c;
      }
      std::vector<double> density_correlation_vec;
      std::vector<double> mogul_z_worst_vec;
      std::vector<int> bad_contacts_vec;
   };

   class ligand_report_absolute_t {
   public:
      // double mogul_z_score;
      double dictionary_geometry_distortion_max;
      double density_correlation;
      probe_clash_score_t pcs;
      ligand_report_absolute_t(double a, double b, probe_clash_score_t p) {
	 density_correlation = a;
	 //   mogul_z_score = b;
	 dictionary_geometry_distortion_max = b;
	 pcs = p;
      }
      ligand_report_percentiles_t make_percentiles() const;
   };
   int ligands_db_sql_callback(void *NotUsed, int argc, char **argv, char **azColName);
}


#endif // LIGAND_CHECK_HH
