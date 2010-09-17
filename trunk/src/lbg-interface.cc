
#ifdef MAKE_ENTERPRISE_TOOLS

#include "enterprise.hh"
#include "lbg-interface.hh"
#include "graphics-info.h"
#include "lbg.hh"

#include "c-interface.h"

void residue_to_ligand_builder(int imol, const char *chain_id, int resno, const char *ins_code) {

   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      CResidue *residue_p = graphics_info_t::molecules[imol].get_residue(chain_id, resno, ins_code);
      if (residue_p) {
	 try {
	    RDKit::RWMol rdkm = coot::rdkit_mol(residue_p, *g.Geom_p());
	    RDKit::ROMol *rdk_mol_with_no_Hs_ro = RDKit::MolOps::removeHs(rdkm);
	    RDKit::RWMol rdk_mol_with_no_Hs = *rdk_mol_with_no_Hs_ro;

	    // clear out any cached properties
	    rdk_mol_with_no_Hs.clearComputedProps();
	    // clean up things like nitro groups
 	    RDKit::MolOps::cleanUp(rdk_mol_with_no_Hs);
 	    // update computed properties on atoms and bonds:
 	    rdk_mol_with_no_Hs.updatePropertyCache();
 	    RDKit::MolOps::Kekulize(rdk_mol_with_no_Hs);
 	    RDKit::MolOps::assignRadicals(rdk_mol_with_no_Hs);
	    
	    // then do aromaticity perception
	    // RDKit::MolOps::setAromaticity(rdkm);
    
	    // set conjugation
	    RDKit::MolOps::setConjugation(rdk_mol_with_no_Hs);
	       
	    // set hybridization
	    RDKit::MolOps::setHybridization(rdk_mol_with_no_Hs); // non-linear ester bonds, yay!

	    // remove bogus chirality specs:
	    RDKit::MolOps::cleanupChirality(rdk_mol_with_no_Hs);


	    if (0) {  // debugging
	       // what are the bond types after sanitization/kekulization?
	       std::cout << "after sanitization/kekulization:: " << std::endl;
	       int n_bonds = rdkm.getNumBonds();
	       for (unsigned int ib=0; ib<n_bonds; ib++) {
		  const RDKit::Bond *bond_p = rdk_mol_with_no_Hs.getBondWithIdx(ib);
		  std::cout << ib << "bond:  " << *bond_p << std::endl;
	       }
	    }

	    CMMDBManager *mol = coot::util::create_mmdbmanager_from_residue(g.molecules[imol].atom_sel.mol,
									    residue_p);
	    int mol_2d_depict_conformer = coot::add_2d_conformer(&rdk_mol_with_no_Hs);
	    lig_build::molfile_molecule_t m =
	       coot::make_molfile_molecule(rdk_mol_with_no_Hs, mol_2d_depict_conformer);
	    if (!mol) {
	       std::cout << "ERROR:: failed to make mol for lbg" << std::endl;
	    } else { 
	       lbg(m, mol, "");
	    }
	 }
	 catch (std::runtime_error coot_error) {
	    std::cout << coot_error.what() << std::endl;
	    std::string m = "Residue type ";
	    m += residue_p->GetResName();
	    m += " not found in dictionary.";
	    add_status_bar_text(m.c_str());
	 }
	 catch (std::exception rdkit_error) {
	    std::cout << rdkit_error.what() << std::endl;
	 }
      } 
   } 
} 


#endif // MAKE_ENTERPRISE_TOOLS
