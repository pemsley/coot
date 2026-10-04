/* coot-utils/servalcat-refine-progress.hh
 *
 * Copyright 2026 by Medical Research Council
 * Author: Paul Emsley
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

#ifndef COOT_UTILS_SERVALCAT_REFINE_PROGRESS_HH
#define COOT_UTILS_SERVALCAT_REFINE_PROGRESS_HH

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <cstddef>

namespace coot {

   //! The parsed per-cycle statistics from a servalcat refinement.
   //!
   //! servalcat writes (and updates as the refinement runs) a statistics file
   //! named "<prefix>_stats.json", a JSON array with one object per cycle. This
   //! class holds that data in a form convenient for plotting - and it does so
   //! generically: every scalar is stored keyed by its servalcat name, so the
   //! same container serves both x-ray refinement (Rwork, Rfree, CCFwork, ...)
   //! and cryo-EM refinement (FSCaverage, ...). The consumer (e.g. the graphs in
   //! the GUI) chooses which keys to read.
   //!
   //! This class has no GTK/Cairo dependency - it is pure data and parsing.
   class servalcat_refine_progress_t {
   public:

      //! one resolution shell of a cycle's "binned" table
      class binned_shell_t {
      public:
         double d_min;
         double d_max;
         //! the per-shell values, e.g. "Rwork", "Rfree", "CCFwork", "Cmpl" ...
         std::map<std::string, double> values;
         binned_shell_t() : d_min(0.0), d_max(0.0) {}
         //! the resolution (A) at the centre of the shell (harmonic mean of the limits)
         double d_mid() const;
      };

      //! the statistics for one refinement cycle
      class cycle_t {
      public:
         int n_cycle; //!< servalcat "Ncyc"
         //! data.summary scalars, e.g. "Rwork", "Rfree", "-LL", "FOM" ...
         std::map<std::string, double> summary;
         //! data.binned, one entry per resolution shell
         std::vector<binned_shell_t> binned;
         //! geom.summary: outer key ("r.m.s.Z", "r.m.s.d.", ...) -> (restraint type -> value)
         std::map<std::string, std::map<std::string, double> > geom;
         cycle_t() : n_cycle(-1) {}
      };

      servalcat_refine_progress_t() {}

      //! the cycles, in the order they appear in the file
      std::vector<cycle_t> cycles;

      //! Parse a stats-json string. On success the cycles are replaced and true
      //! is returned. On a JSON error (e.g. a half-written file caught
      //! mid-update during polling) false is returned and the existing cycles
      //! are left unchanged, so a transient read does not wipe the graph.
      bool parse_json_string(const std::string &json_string);

      //! As parse_json_string(), but reads from a file. Returns false if the
      //! file cannot be read or parsed.
      bool parse_json_file(const std::string &file_name);

      bool empty() const { return cycles.empty(); }
      std::size_t size() const { return cycles.size(); }

      //! index of the final (highest-numbered) cycle, or -1 if there are none
      int last_cycle_index() const;

      // ---- convenience accessors for drawing ---------------------------------

      //! A summary value across cycles, as (cycle_number, value) pairs, skipping
      //! cycles that lack the key. e.g. summary_series("Rfree").
      std::vector<std::pair<int, double> >
      summary_series(const std::string &key) const;

      //! A geometry restraint value across cycles, e.g.
      //! geom_series("r.m.s.Z", "Bond distances, non H").
      std::vector<std::pair<int, double> >
      geom_series(const std::string &outer_key, const std::string &restraint_type) const;

      //! The binned values of one column at a given cycle, as (1/d^2, value)
      //! pairs ordered by increasing 1/d^2 (low to high resolution), e.g.
      //! binned_series(last_cycle_index(), "Rfree").
      std::vector<std::pair<double, double> >
      binned_series(std::size_t cycle_index, const std::string &column) const;

      // ---- introspection (so the GUI need not hardcode key names) ------------

      //! the union of summary keys present across all cycles
      std::vector<std::string> summary_keys() const;
      //! the union of binned-column names present across all cycles
      std::vector<std::string> binned_columns() const;
      //! the union of geom outer keys ("r.m.s.Z", ...) present across all cycles
      std::vector<std::string> geom_outer_keys() const;
      //! the union of geom restraint-type names present across all cycles
      std::vector<std::string> geom_restraint_types() const;

   };

}

#endif // COOT_UTILS_SERVALCAT_REFINE_PROGRESS_HH
