
#include "lidia-core/rdkit-interface.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
// #include "generic-display-object.hh"

#include "coot-utils/simple-mesh.hh"

namespace chemical_features { 

   // not for public access 
   // void show(int imol, const RDKit::ROMol &rdkm, std::string name);

   // choose iconf = 0 if unsure.
   std::vector<coot::simple_mesh_t> generate_meshes(int imol, const RDKit::ROMol &rdkm, int iconf, const std::string &name);

   std::vector<coot::simple_mesh_t> generate_meshes(int imol, mmdb::Residue *residue_p, const coot::protein_geometry &geom);

   std::pair<bool, clipper::Coord_orth> get_normal_info(RDKit::MolChemicalFeature *feat,
							const RDKit::ROMol &mol,
							const RDKit::Conformer &conf);
   std::pair<bool, clipper::Coord_orth> get_normal_info_aromatic(RDKit::MolChemicalFeature *feat,
								 const RDKit::Conformer &conf);
   std::pair<bool, clipper::Coord_orth> get_normal_info_donor(RDKit::MolChemicalFeature *feat,
							      const RDKit::ROMol &mol,
							      const RDKit::Conformer &conf);
   // return null on inability to make factory.
   RDKit::MolChemicalFeatureFactory *get_feature_factory();
   
}
