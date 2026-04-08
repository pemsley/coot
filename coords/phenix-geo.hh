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

   namespace phenix_geo {

      atom_spec_t parse_line_for_atom_spec(const std::string &l);

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
                         const atom_spec_t &a2) : atom_1(a1), atom_2(a2), geom_set_flag(false) {}
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

      class phenix_geo_angle {
         public:
         atom_spec_t atom_1;
         atom_spec_t atom_2;
         atom_spec_t atom_3;
         double ideal;
         double model;
         double delta;
         double sigma;
         double weight;
         double residual;
         bool geom_set_flag;
         phenix_geo_angle(const atom_spec_t &a1,
                          const atom_spec_t &a2,
                          const atom_spec_t &a3) : atom_1(a1), atom_2(a2), atom_3(a3), geom_set_flag(false) {}
         void set_geom(double ideal_in, double model_in, double delta_in,
                       double sigma_in, double weight_in, double residual_in) {
           ideal  = ideal_in;
           model  = model_in;
           delta  = delta_in;
           sigma  = sigma_in;
           weight = weight_in;
           residual = residual_in;
         }
      };

      class phenix_geo_bonds {

      public:
         std::vector<phenix_geo_bond> bonds;
         phenix_geo_bonds() {}
         unsigned int size() const { return bonds.size(); }
#ifndef SWIG
         const phenix_geo_bond &operator[](const unsigned int &idx) const { return bonds[idx]; }
#endif
         void add_bond(const phenix_geo_bond &b) {
            bonds.push_back(b);
         } 
      };

      class phenix_geo_angles {
      public:
         phenix_geo_angles() {}
         std::vector<phenix_geo_angle> angles;
         unsigned int size() const { return angles.size(); }
         void add_angle(const phenix_geo_angle &a) { angles.push_back(a); }
      };

      class phenix_geo_dihedral {
      public:
         atom_spec_t atom_1;
         atom_spec_t atom_2;
         atom_spec_t atom_3;
         atom_spec_t atom_4;
         double ideal;
         double model;
         double delta;
         int periodicity; // sinusoidal: 1,2,3... harmonic: 0
         double sigma;
         double weight;
         double residual;
         bool is_sinusoidal; // true for sinusoidal, false for harmonic
         bool geom_set_flag;
         phenix_geo_dihedral(const atom_spec_t &a1,
                             const atom_spec_t &a2,
                             const atom_spec_t &a3,
                             const atom_spec_t &a4) :
            atom_1(a1), atom_2(a2), atom_3(a3), atom_4(a4),
            periodicity(0), is_sinusoidal(true), geom_set_flag(false) {}
         void set_geom(double ideal_in, double model_in, double delta_in,
                       int periodicity_in, double sigma_in, double weight_in,
                       double residual_in, bool is_sinusoidal_in) {
            ideal = ideal_in;
            model = model_in;
            delta = delta_in;
            periodicity = periodicity_in;
            sigma = sigma_in;
            weight = weight_in;
            residual = residual_in;
            is_sinusoidal = is_sinusoidal_in;
            geom_set_flag = true;
         }
      };

      class phenix_geo_dihedrals {
      public:
         phenix_geo_dihedrals() {}
         std::vector<phenix_geo_dihedral> dihedrals;
         unsigned int size() const { return dihedrals.size(); }
         void add_dihedral(const phenix_geo_dihedral &d) { dihedrals.push_back(d); }
      };

      class phenix_geo_chiral {
      public:
         atom_spec_t atom_centre;
         atom_spec_t atom_1;
         atom_spec_t atom_2;
         atom_spec_t atom_3;
         bool both_signs;
         double ideal;
         double model;
         double delta;
         double sigma;
         double weight;
         double residual;
         bool geom_set_flag;
         phenix_geo_chiral(const atom_spec_t &centre,
                           const atom_spec_t &a1,
                           const atom_spec_t &a2,
                           const atom_spec_t &a3) :
            atom_centre(centre), atom_1(a1), atom_2(a2), atom_3(a3),
            both_signs(false), geom_set_flag(false) {}
         void set_geom(bool both_signs_in, double ideal_in, double model_in,
                       double delta_in, double sigma_in, double weight_in,
                       double residual_in) {
            both_signs = both_signs_in;
            ideal = ideal_in;
            model = model_in;
            delta = delta_in;
            sigma = sigma_in;
            weight = weight_in;
            residual = residual_in;
            geom_set_flag = true;
         }
      };

      class phenix_geo_chirals {
      public:
         phenix_geo_chirals() {}
         std::vector<phenix_geo_chiral> chirals;
         unsigned int size() const { return chirals.size(); }
         void add_chiral(const phenix_geo_chiral &c) { chirals.push_back(c); }
      };

      class phenix_geo_nonbonded {
      public:
         atom_spec_t atom_1;
         atom_spec_t atom_2;
         double model;
         double vdw;
         bool geom_set_flag;
         phenix_geo_nonbonded(const atom_spec_t &a1,
                              const atom_spec_t &a2) :
            atom_1(a1), atom_2(a2), geom_set_flag(false) {}
         void set_geom(double model_in, double vdw_in) {
            model = model_in;
            vdw = vdw_in;
            geom_set_flag = true;
         }
      };

      class phenix_geo_nonbondeds {
      public:
         phenix_geo_nonbondeds() {}
         std::vector<phenix_geo_nonbonded> nonbondeds;
         unsigned int size() const { return nonbondeds.size(); }
         void add_nonbonded(const phenix_geo_nonbonded &nb) { nonbondeds.push_back(nb); }
      };

      class phenix_geometry {
         public:
         phenix_geometry() {}
         phenix_geo_bonds geo_bonds;
         phenix_geo_angles geo_angles;
         phenix_geo_dihedrals geo_dihedrals;
         phenix_geo_chirals geo_chirals;
         phenix_geo_nonbondeds geo_nonbondeds;
         void parse(const std::string &file_name);
      };
   }
}

#endif // PHENIX_GEO_HH
