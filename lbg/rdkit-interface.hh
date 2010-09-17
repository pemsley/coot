
#ifndef ENTERPRISE_HH
#define ENTERPRISE_HH

#include "use-rdkit.hh"

#include "mmdb_manager.h"
#include "protein-geometry.hh"
#include "lbg-molfile.hh"

namespace coot { 

   // can throw an runtime_error exception (residue not in dictionary)
   // 
   RDKit::RWMol rdkit_mol(CResidue *residue_p, const protein_geometry &geom);
   int add_2d_conformer(RDKit::ROMol *rdkmol_in); // tweak rdkmol_in
   RDKit::Bond::BondType convert_bond_type(const std::string &t);

   lig_build::molfile_molecule_t make_molfile_molecule(const RDKit::ROMol &rdkm, int iconf);
   // lig_build::molfile_molecule_t make_molfile_molecule(const RDKit::RWMol &rdkm);
   lig_build::bond_t::bond_type_t convert_bond_type(const RDKit::Bond::BondType &type);

} 

#endif // ENTERPRISE_HH
