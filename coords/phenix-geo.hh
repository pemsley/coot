
#ifndef PHENIX_GEO_HH
#define PHENIX_GEO_HH

#include <vector>
#include "coot-utils/residue-and-atom-specs.hh"

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
   std::ostream &operator<<(std::ostream &s, phenix_geo_bond bg);

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
