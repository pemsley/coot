

//  ---------------------------------------------
//           20230109-PE this needs to be restored
//  ---------------------------------------------


// #include "c-interface.h" // for set_display_generic_object()
//   #include "c-interface-widgets.hh" // for add_generic_display_object
// #include "c-interface-gtk-widgets.h" // for move_molecule_to_screen_centre_internal()

//! \brief
//! import a molecule from a smiles string
//!
//! RDKit is used to interpret the SMILES string
//!
//! no dictionary is generated
//!
//! @return a molecule number or -1 on failure
int import_rdkit_mol_from_smiles(const std::string &smiles_str, const std::string &comp_id) {

   int imol = -1;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   try {
      RDKit::RWMol *m = RDKit::SmilesToMol(smiles_str);
      if (m) {
	 bool explicit_only = false;
	 bool add_coords = true;
	 RDKit::MolOps::addHs(*m, explicit_only, add_coords);
	 int iconf = RDKit::DGeomHelpers::EmbedMolecule(*m);
	 if(iconf < 0) {
	    std::cout << "WARNING:: RDKit::embedding failed." << std::endl;
	 } else {
	    double vdwThresh=10.0;
	    int confId = -1;
	    bool ignoreInterfragInteractions=true; // sensible?
	    ForceFields::ForceField *ff =
	       RDKit::UFF::constructForceField(*m,
					       vdwThresh, confId,
					       ignoreInterfragInteractions);
	    ff->initialize();
	    int maxIters = 500;
	    int res=ff->minimize(maxIters);
	    delete ff;

	    mmdb::Residue *residue_p = coot::make_residue(*m, iconf, comp_id);
	    if (residue_p) {
	       mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
	       if (mol) {
		  graphics_info_t g;
		  imol = g.create_molecule();
		  std::string label = "Imported ";
		  label += comp_id;
		  g.molecules[imol].install_model(imol, mol, g.Geom_p(), label, 1, false, false);
		  move_molecule_to_screen_centre_internal(imol);
	       }
	       delete residue_p;
	    }
	 }
      } else {
	 std::cout << "WARNING:: BAD SMILES " << smiles_str << std::endl;
	 std::string s = "WARNING:: Bad SMILES: " + smiles_str;
	 info_dialog(s.c_str());
      }
   }
   catch (const std::exception &e) {
      std::cout << "WARNING:: exception caught " << e.what() << std::endl;
   }
   
#endif // MAKE_ENHANCED_LIGAND_TOOLS

   return imol;
}

