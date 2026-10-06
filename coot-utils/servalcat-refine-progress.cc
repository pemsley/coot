/* coot-utils/servalcat-refine-progress.cc
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

#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>

#include "json.hpp"
using json = nlohmann::json;

#include "servalcat-refine-progress.hh"

double
coot::servalcat_refine_progress_t::binned_shell_t::d_mid() const {

   // harmonic mean of the resolution limits (matches servalcat-tracker.cc)
   if (d_min <= 0.0 || d_max <= 0.0) return 0.0;
   return 2.0 / (1.0/d_min + 1.0/d_max);
}

// Parse one cycle object into a cycle_t. A missing/odd section is skipped
// rather than throwing, so that a slightly different schema (e.g. cryo-EM
// vs x-ray) still yields what data it can.
//
namespace {

   coot::servalcat_refine_progress_t::cycle_t
   parse_one_cycle(const json &j_cycle) {

      coot::servalcat_refine_progress_t::cycle_t cycle;

      if (j_cycle.contains("Ncyc") && j_cycle["Ncyc"].is_number())
         cycle.n_cycle = j_cycle["Ncyc"].get<int>();

      if (j_cycle.contains("data")) {
         const json &j_data = j_cycle["data"];

         // data.summary: a flat object of scalars
         if (j_data.contains("summary") && j_data["summary"].is_object()) {
            for (auto it = j_data["summary"].begin(); it != j_data["summary"].end(); ++it)
               if (it.value().is_number())
                  cycle.summary[it.key()] = it.value().get<double>();
         }

         // data.binned: an array of resolution shells
         if (j_data.contains("binned") && j_data["binned"].is_array()) {
            for (const auto &j_shell : j_data["binned"]) {
               if (! j_shell.is_object()) continue;
               coot::servalcat_refine_progress_t::binned_shell_t shell;
               for (auto it = j_shell.begin(); it != j_shell.end(); ++it) {
                  if (! it.value().is_number()) continue;
                  double v = it.value().get<double>();
                  if (it.key() == "d_min") shell.d_min = v;
                  else if (it.key() == "d_max") shell.d_max = v;
                  else shell.values[it.key()] = v;
               }
               cycle.binned.push_back(shell);
            }
         }
      }

      // geom.summary: outer key -> (restraint type -> value)
      if (j_cycle.contains("geom")) {
         const json &j_geom = j_cycle["geom"];
         if (j_geom.contains("summary") && j_geom["summary"].is_object()) {
            const json &j_gs = j_geom["summary"];
            for (auto it_outer = j_gs.begin(); it_outer != j_gs.end(); ++it_outer) {
               if (! it_outer.value().is_object()) continue;
               std::map<std::string, double> &inner = cycle.geom[it_outer.key()];
               for (auto it = it_outer.value().begin(); it != it_outer.value().end(); ++it)
                  if (it.value().is_number())
                     inner[it.key()] = it.value().get<double>();
            }
         }
      }

      return cycle;
   }
}

bool
coot::servalcat_refine_progress_t::parse_json_string(const std::string &json_string) {

   if (json_string.empty()) return false;

   try {
      json j = json::parse(json_string);
      if (! j.is_array()) {
         std::cout << "WARNING:: servalcat_refine_progress_t: top-level JSON is not an array"
                   << std::endl;
         return false;
      }
      std::vector<cycle_t> new_cycles;
      new_cycles.reserve(j.size());
      for (const auto &j_cycle : j)
         if (j_cycle.is_object())
            new_cycles.push_back(parse_one_cycle(j_cycle));

      cycles.swap(new_cycles); // only replace on a clean parse
      return true;
   }
   catch (const json::exception &e) {
      // a half-written file caught mid-update is the expected cause here - keep
      // the current cycles and let the caller try again on the next poll.
      std::cout << "WARNING:: servalcat_refine_progress_t: JSON parse error " << e.what()
                << std::endl;
      return false;
   }
}

bool
coot::servalcat_refine_progress_t::parse_json_file(const std::string &file_name) {

   std::ifstream f(file_name);
   if (! f) {
      std::cout << "WARNING:: servalcat_refine_progress_t: cannot open " << file_name << std::endl;
      return false;
   }
   std::stringstream ss;
   ss << f.rdbuf();
   return parse_json_string(ss.str());
}

int
coot::servalcat_refine_progress_t::last_cycle_index() const {

   if (cycles.empty()) return -1;
   // the cycle with the highest n_cycle (fall back to the last element)
   int best_index = static_cast<int>(cycles.size()) - 1;
   int best_n = cycles[best_index].n_cycle;
   for (std::size_t i=0; i<cycles.size(); i++) {
      if (cycles[i].n_cycle > best_n) {
         best_n = cycles[i].n_cycle;
         best_index = static_cast<int>(i);
      }
   }
   return best_index;
}

std::vector<std::pair<int, double> >
coot::servalcat_refine_progress_t::summary_series(const std::string &key) const {

   std::vector<std::pair<int, double> > v;
   for (const auto &cycle : cycles) {
      std::map<std::string, double>::const_iterator it = cycle.summary.find(key);
      if (it != cycle.summary.end())
         v.push_back(std::make_pair(cycle.n_cycle, it->second));
   }
   return v;
}

std::vector<std::pair<int, double> >
coot::servalcat_refine_progress_t::geom_series(const std::string &outer_key,
                                               const std::string &restraint_type) const {

   std::vector<std::pair<int, double> > v;
   for (const auto &cycle : cycles) {
      std::map<std::string, std::map<std::string, double> >::const_iterator it_outer =
         cycle.geom.find(outer_key);
      if (it_outer != cycle.geom.end()) {
         std::map<std::string, double>::const_iterator it = it_outer->second.find(restraint_type);
         if (it != it_outer->second.end())
            v.push_back(std::make_pair(cycle.n_cycle, it->second));
      }
   }
   return v;
}

std::vector<std::pair<double, double> >
coot::servalcat_refine_progress_t::binned_series(std::size_t cycle_index,
                                                 const std::string &column) const {

   std::vector<std::pair<double, double> > v;
   if (cycle_index >= cycles.size()) return v;
   const cycle_t &cycle = cycles[cycle_index];
   for (const auto &shell : cycle.binned) {
      std::map<std::string, double>::const_iterator it = shell.values.find(column);
      if (it != shell.values.end()) {
         double d_mid = shell.d_mid();
         if (d_mid > 0.0) {
            double inv_d_sq = 1.0 / (d_mid * d_mid);
            v.push_back(std::make_pair(inv_d_sq, it->second));
         }
      }
   }
   std::sort(v.begin(), v.end(),
             [] (const std::pair<double, double> &a, const std::pair<double, double> &b) {
                return a.first < b.first;
             });
   return v;
}

// ---- introspection ---------------------------------------------------------

std::vector<std::string>
coot::servalcat_refine_progress_t::summary_keys() const {
   std::set<std::string> s;
   for (const auto &cycle : cycles)
      for (const auto &kv : cycle.summary)
         s.insert(kv.first);
   return std::vector<std::string>(s.begin(), s.end());
}

std::vector<std::string>
coot::servalcat_refine_progress_t::binned_columns() const {
   std::set<std::string> s;
   for (const auto &cycle : cycles)
      for (const auto &shell : cycle.binned)
         for (const auto &kv : shell.values)
            s.insert(kv.first);
   return std::vector<std::string>(s.begin(), s.end());
}

std::vector<std::string>
coot::servalcat_refine_progress_t::geom_outer_keys() const {
   std::set<std::string> s;
   for (const auto &cycle : cycles)
      for (const auto &kv : cycle.geom)
         s.insert(kv.first);
   return std::vector<std::string>(s.begin(), s.end());
}

std::vector<std::string>
coot::servalcat_refine_progress_t::geom_restraint_types() const {
   std::set<std::string> s;
   for (const auto &cycle : cycles)
      for (const auto &outer : cycle.geom)
         for (const auto &kv : outer.second)
            s.insert(kv.first);
   return std::vector<std::string>(s.begin(), s.end());
}
