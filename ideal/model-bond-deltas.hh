
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
	 unsigned int n;
	 xyz_deltas_t() {n = 0;}
	 void add(const clipper::Coord_orth &d) {
	    x[0].push_back(d.x());
	    x[1].push_back(d.y());
	    x[2].push_back(d.z());
	    ++n;
	 }
	 unsigned int size() const { return n;}
      };
      xyz_deltas_t xyz;
      void add(const clipper::Coord_orth &d) {
	 xyz.add(d);
      }
      unsigned int size() const { return xyz.size(); }
   };
}

#endif // MODEL_BOND_DELTAS_HH
