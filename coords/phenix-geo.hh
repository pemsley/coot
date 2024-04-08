/* canyon/phenix-geo.hh
 * 
 * Copyright 2014 by Medical Research Council
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
 * Lesser General Public License for more details.
 * 
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#ifndef PHENIX_GEO_HH
#define PHENIX_GEO_HH

#include <vector>
#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class phenix_geo_bond {
   public:
      atom_spec_t atom_1;
      atom_spec_t atom_2;
      double ideal;
      double model;
      double delta;
      double sigma;
      double weight;
      double residual;
      bool geom_set_flag;
      phenix_geo_bond(const atom_spec_t &a1,
		      const atom_spec_t &a2) {
	 geom_set_flag = false;
	 atom_1 = a1;
	 atom_2 = a2;
      }
      void set_geom(double ideal_in, double model_in, double delta_in,
		    double sigma_in, double weight_in, double residual_in) {
	 ideal  = ideal_in;
	 model  = model_in;
	 delta  = delta_in;
	 sigma  = sigma_in;
	 weight = weight_in;
	 residual = residual_in;
      }
#ifndef SWIG
      friend std::ostream &operator<<(std::ostream &s, phenix_geo_bond bg);
#endif
   };
#ifndef SWIG
   std::ostream &operator<<(std::ostream &s, phenix_geo_bond bg);
#endif

   class phenix_geo_bonds {
      atom_spec_t parse_line_for_atom_spec(const std::string &l) const;
      
   public:
      std::vector<phenix_geo_bond> bonds;
      phenix_geo_bonds() {}
      phenix_geo_bonds(const std::string &file_name);
      unsigned int size() const { return bonds.size(); }
#ifndef SWIG
      const phenix_geo_bond &operator[](const unsigned int &idx) const { return bonds[idx]; }
#endif
      void add_bond(const phenix_geo_bond &b) {
	 bonds.push_back(b);
      } 
   };
}

#endif // PHENIX_GEO_HH
