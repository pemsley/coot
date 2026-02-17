

#include "molecules-container.hh"
#include "utils/base64-encode-decode.hh"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/UFF/Builder.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <ForceField/ForceField.h>

// Prevents preprocessor substitution of `VERSION` in `MolPickler.h`
#ifndef RD_MOLPICKLE_H
#ifdef VERSION
#define __COOT_VERSION_VALUE VERSION
#undef VERSION
#endif
#include <GraphMol/MolPickler.h>
#ifdef __COOT_VERSION_VALUE
#define VERSION __COOT_VERSION_VALUE
#undef __COOT_VERSION_VALUE
#endif
#endif // RD_MOLPICKLE_H

#include "lidia-core/rdkit-interface.hh"
// #define LIBCOOTAPI_BUILD // means no python - this is set with cmake now.
#include "pyrogen/restraints.hh"
#endif

// C++ port of pyrogen/pyrogen.py pad_atom_name() and add_atom_names()
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
namespace {

   // Pad atom name to PDB 4-char format
   // 1-char elements: " X# " or " X## "
   // 2-char elements: "XX# " or "XX##"
   std::string pad_atom_name(const std::string &name, const std::string &element) {
      std::string padded = name;
      if (element.length() == 1) {
         if (name.length() == 2)
            padded = " " + name + " ";
         else if (name.length() == 3)
            padded = " " + name;
         else if (name.length() == 4)
            padded = name; // already full width
      } else if (element.length() == 2) {
         if (name.length() == 2)
            padded = name + "  ";
         else if (name.length() == 3)
            padded = name + " ";
         else if (name.length() == 4)
            padded = name;
      }
      return padded;
   }

   // Assign unique PDB-style atom names to all atoms in the molecule
   void add_atom_names(RDKit::RWMol &mol) {
      std::map<unsigned int, int> nz; // element atomic number -> count
      std::vector<std::string> atom_names;
      for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
         RDKit::Atom *atom = mol.getAtomWithIdx(i);
         std::string existing_name;
         try {
            atom->getProp("name", existing_name);
            atom_names.push_back(existing_name);
         } catch (KeyErrorException &) {
            unsigned int z = atom->getAtomicNum();
            if (nz.find(z) != nz.end())
               nz[z]++;
            else
               nz[z] = 1;
            std::string ele = atom->getSymbol();
            // Make element uppercase (matches Python .upper())
            for (auto &c : ele) c = std::toupper(c);

            int inc = 0;
            std::string name = ele + std::to_string(nz[z] + inc);
            std::string p_name = pad_atom_name(name, ele);
            while (std::find(atom_names.begin(), atom_names.end(), p_name) != atom_names.end()) {
               inc++;
               name = ele + std::to_string(nz[z] + inc);
               p_name = pad_atom_name(name, ele);
            }
            atom->setProp("name", p_name);
            atom_names.push_back(p_name);
         }
      }
   }

}
#endif


int
molecules_container_t::pyrogen_from_SMILES(const std::string &smiles_string,
                                           const std::string &compound_id) {

   int imol = -1;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   try {
      // 1. Parse SMILES
      RDKit::RWMol *mol_raw = RDKit::SmilesToMol(smiles_string);
      if (!mol_raw) {
         std::cout << "WARNING:: pyrogen_from_smiles(): failed to parse SMILES: " << smiles_string << std::endl;
         return -1;
      }
      std::unique_ptr<RDKit::RWMol> mol(mol_raw);

      // 2. Add hydrogens
      RDKit::MolOps::addHs(*mol);

      // 3. Generate 3D conformer
      int conf_id = RDKit::DGeomHelpers::EmbedMolecule(*mol);
      if (conf_id < 0) {
         std::cout << "WARNING:: pyrogen_from_smiles(): EmbedMolecule() failed for " << smiles_string << std::endl;
         return -1;
      }

      // 4. UFF force field minimize
      ForceFields::ForceField *ff = RDKit::UFF::constructForceField(*mol);
      if (ff) {
         ff->initialize();
         ff->minimize();
         delete ff;
      }

      // 5. Gasteiger charges
      RDKit::computeGasteigerCharges(*mol);

      // 6. Add atom names (C++ port of pyrogen add_atom_names())
      add_atom_names(*mol);

      // 7. Assign CCP4 energy types
      coot::set_energy_lib_atom_types(mol.get());

      // 8. Build dictionary
      auto dict_pair = coot::mmcif_dict_from_mol_using_energy_lib(compound_id, compound_id,
                                                                   static_cast<const RDKit::ROMol &>(*mol),
                                                                   true, true);
      if (!dict_pair.first) {
         std::cout << "WARNING:: pyrogen_from_smiles(): mmcif_dict_from_mol_using_energy_lib() failed" << std::endl;
         return -1;
      }

      // 9. Store in geometry
      int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
      geom.replace_monomer_restraints(compound_id, imol_enc, dict_pair.second);

      // 10-11. Create molecule from dictionary
      mmdb::Manager *mmdb_mol = geom.mol_from_dictionary(compound_id, imol_enc, true);
      if (mmdb_mol) {
         imol = molecules.size();
         std::string name = compound_id + "_from_smiles";
         atom_selection_container_t asc = make_asc(mmdb_mol);
         coot::molecule_t m(asc, imol, name);
         molecules.push_back(m);
      } else {
         std::cout << "WARNING:: pyrogen_from_smiles(): mol_from_dictionary() returned null for "
                   << compound_id << std::endl;
      }

   } catch (const std::exception &e) {
      std::cout << "WARNING:: pyrogen_from_smiles() exception: " << e.what() << std::endl;
   }

#else
   std::cout << "WARNING:: pyrogen_from_smiles(): not built with MAKE_ENHANCED_LIGAND_TOOLS" << std::endl;
#endif

   return imol;
}

int molecules_container_t::pyrogen_from_ccd_file(const std::string &ccd_file_name) {

   int imol = -1;

   // 1. Parse CIF dictionary
   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   coot::read_refmac_mon_lib_info_t r = geom.init_refmac_mon_lib(ccd_file_name, cif_dictionary_read_number,
                                                                 imol_enc);
   cif_dictionary_read_number++;

   if (!r.success || r.comp_id.empty()) {
      std::cout << "WARNING:: pyrogen_from_ccd_file(): failed to read " << ccd_file_name << std::endl;
      return -1;
   }

   // 2-4. Create molecule from dictionary
   mmdb::Manager *mmdb_mol = geom.mol_from_dictionary(r.comp_id, imol_enc, true);
   if (mmdb_mol) {
      imol = molecules.size();
      std::string name = r.comp_id + "_from_ccd";
      atom_selection_container_t asc = make_asc(mmdb_mol);
      coot::molecule_t m(asc, imol, name);
      molecules.push_back(m);
   } else {
      std::cout << "WARNING:: pyrogen_from_ccd_file(): mol_from_dictionary() returned null for "
                << r.comp_id << std::endl;
   }

   return imol;
}

//! this is the interface to use from the molecule sketcher (say) where the
//! calling function has an RDKit Mol
//!
//! @param rdkit_mol_pickled_string the rdkit mol as a picked string
//! @param compound_id is the compound_id that should be assigned to the new dictionary
//!        and molecule
//! @return the new molecule index or -1 on failure
int molecules_container_t::pyrogen_from_rdkit_mol_pickle_base64(const std::string &rdkit_mol_pickle_base64_string,
                                                                const std::string &compound_id) {

   int imol = -1;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   try {
      // 1. Decode base64 and unpickle
      std::string pickle_string = moorhen_base64::base64_decode(rdkit_mol_pickle_base64_string);
      std::unique_ptr<RDKit::RWMol> mol = std::make_unique<RDKit::RWMol>();
      RDKit::MolPickler::molFromPickle(pickle_string, mol.get());

      if (!mol || mol->getNumAtoms() == 0) {
         std::cout << "WARNING:: pyrogen_from_rdkit_mol_pickle_base64(): failed to unpickle molecule" << std::endl;
         return -1;
      }

      // 2. Add hydrogens if not already present
      bool has_hydrogens = false;
      for (unsigned int i = 0; i < mol->getNumAtoms(); i++) {
         if (mol->getAtomWithIdx(i)->getAtomicNum() == 1) {
            has_hydrogens = true;
            break;
         }
      }
      if (!has_hydrogens) {
         RDKit::MolOps::addHs(*mol);
      }

      // 3. Generate 3D conformer if not present
      if (mol->getNumConformers() == 0) {
         int conf_id = RDKit::DGeomHelpers::EmbedMolecule(*mol);
         if (conf_id < 0) {
            std::cout << "WARNING:: pyrogen_from_rdkit_mol_pickle_base64(): EmbedMolecule() failed" << std::endl;
            return -1;
         }

         // 4. UFF force field minimize (only if we just embedded)
         ForceFields::ForceField *ff = RDKit::UFF::constructForceField(*mol);
         if (ff) {
            ff->initialize();
            ff->minimize();
            delete ff;
         }
      }

      // 5. Gasteiger charges
      RDKit::computeGasteigerCharges(*mol);

      // 6. Add atom names (C++ port of pyrogen add_atom_names())
      add_atom_names(*mol);

      // 7. Assign CCP4 energy types
      coot::set_energy_lib_atom_types(mol.get());

      // 8. Build dictionary
      auto dict_pair = coot::mmcif_dict_from_mol_using_energy_lib(compound_id, compound_id,
                                                                   static_cast<const RDKit::ROMol &>(*mol),
                                                                   true, true);
      if (!dict_pair.first) {
         std::cout << "WARNING:: pyrogen_from_rdkit_mol_pickle_base64(): mmcif_dict_from_mol_using_energy_lib() failed" << std::endl;
         return -1;
      }

      // 9. Store in geometry
      int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
      geom.replace_monomer_restraints(compound_id, imol_enc, dict_pair.second);

      // 10-11. Create molecule from dictionary
      mmdb::Manager *mmdb_mol = geom.mol_from_dictionary(compound_id, imol_enc, true);
      if (mmdb_mol) {
         imol = molecules.size();
         std::string name = compound_id + "_from_pickle";
         atom_selection_container_t asc = make_asc(mmdb_mol);
         coot::molecule_t m(asc, imol, name);
         molecules.push_back(m);
      } else {
         std::cout << "WARNING:: pyrogen_from_rdkit_mol_pickle_base64(): mol_from_dictionary() returned null for "
                   << compound_id << std::endl;
      }

   } catch (const std::exception &e) {
      std::cout << "WARNING:: pyrogen_from_rdkit_mol_pickle_base64() exception: " << e.what() << std::endl;
   }

#else
   std::cout << "WARNING:: pyrogen_from_rdkit_mol_pickle_base64(): not built with MAKE_ENHANCED_LIGAND_TOOLS" << std::endl;
#endif

   return imol;
}

