
#ifndef COOT_UTILS_CFC_HH
#define COOT_UTILS_CFC_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <vector>
#include <string>
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include <GraphMol/MolChemicalFeatures/MolChemicalFeatureFactory.h>

namespace cfc {

   class input_info_t {
   public:
      mmdb::Manager *mol;
      int imol;
      std::string res_name;
      input_info_t(mmdb::Manager *mol, int imol, const std::string &rn) : mol(mol), imol(imol), res_name(rn) {}
   };

   class typed_cluster_t {
   public:
      std::string family;
      std::string type;
      int idx;
      RDGeom::Point3D pos;
      std::vector<RDGeom::Point3D> contributing_points;
      std::vector<std::pair<int, coot::residue_spec_t> > imols_with_specs;
      typed_cluster_t(const std::string &f, const std::string &t, int idx) :
         family(f), type(t), idx(idx) {}
      void add_imol(int imol, const coot::residue_spec_t &rs) {
         imols_with_specs.push_back(std::make_pair(imol, rs)); }
      bool imol_is_part_of_cluster(int imol_in) const {
         bool state = false;
         for (const auto &item : imols_with_specs) {
            if (item.first == imol_in) {
               state = true;
               break;
            }
         }
         return state;
      }
      std::string make_name() const {
         return family + " " + type + " " + std::to_string(idx);
      }
   };

   class water_info_t {
   public:
      int imol;
      coot::residue_spec_t residue_spec;
      RDGeom::Point3D pos;
      water_info_t(int imol_in, const coot::residue_spec_t &rs, const RDGeom::Point3D &p) :
         imol(imol_in), residue_spec(rs), pos(p) {}
      water_info_t() : imol(-1) {}
   };

   std::pair<std::vector<typed_cluster_t>, std::vector<std::vector<water_info_t> > >
   chemical_feature_clustering(const std::vector<input_info_t> &mol_infos,
                               const coot::protein_geometry &geom);

}

#endif // MAKE_ENHANCED_LIGAND_TOOLS
#endif // COOT_UTILS_CFC_HH
