/*
 * ideal/model-bond-deltas.hh
 * 
 * Copyright 2017 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#ifndef MODEL_BOND_DELTAS_HH
#define MODEL_BOND_DELTAS_HH

#include <map>
#include <string>
#include <vector>

#include "geometry/protein-geometry.hh"

namespace coot {

   class model_bond_deltas {

      void fill_dictionaries(std::map<std::string, dictionary_residue_restraints_t> *dictionaries) const;
      void resolve(const std::map<std::string, dictionary_residue_restraints_t> &dictionaries);
      int imol;
      mmdb::Manager *mol;
      protein_geometry *geom_p;
   public:
      model_bond_deltas(mmdb::Manager *mol_in, int imol_in, protein_geometry *geom_p_in);
      model_bond_deltas() {
	 mol = 0;
	 imol = protein_geometry::IMOL_ENC_UNSET;
	 geom_p = 0;
      }
      void resolve();
      class xyz_deltas_t {
      public:
	 std::vector<double> x[3];
	 std::vector<double> deltas; // not resolved along axes
	 unsigned int n;
	 xyz_deltas_t() {n = 0;}
	 void add(double delta_in, const clipper::Coord_orth &d) {
	    deltas.push_back(delta_in);
	    x[0].push_back(d.x());
	    x[1].push_back(d.y());
	    x[2].push_back(d.z());
	    ++n;
	 }
	 unsigned int size() const { return n;}
      };
      xyz_deltas_t xyzd;
      void add(double delta_in, const clipper::Coord_orth &d) {
	 xyzd.add(delta_in, d);
      }
      unsigned int size() const { return xyzd.size(); }
   };
}

#endif // MODEL_BOND_DELTAS_HH
