#ifndef SDF_INTERFACE_FOR_EXPORT_HH
#define SDF_INTERFACE_FOR_EXPORT_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include "lidia-core/rdkit-interface.hh"
#include <GraphMol/MolChemicalFeatures/MolChemicalFeature.h>
#include "coot-utils/simple-mesh.hh"

namespace chemical_features {

   std::vector<coot::simple_mesh_t> generate_meshes(int imol, mmdb::Residue *residue_p, const coot::protein_geometry &geom);

   // choose iconf = 0 if unsure.
   std::vector<coot::simple_mesh_t> generate_meshes(int imol, const RDKit::ROMol &rdkm, int iconf, const std::string &name);
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS

#endif // SDF_INTERFACE_FOR_EXPORT_HH
