
#include "Python.h"
#include "graphics-info.h"
#include "molecule-class-info.h"

// make and add to the scene
void
molecule_class_info_t::make_molecularrepresentationinstance() {

#ifdef USE_MOLECULES_TO_TRIANGLES
#ifdef HAVE_CXX11
   if (atom_sel.mol) {

      auto my_mol = std::make_shared<MyMolecule> (atom_sel.mol);

      // (.. (seletionText, name), "colour")
      // sheet colour: 130 190 170 -> #82BEAA

//       auto init_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("/*/*/*.*/*:*", "NoSecondary"), "chocolate");
//       auto none_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_None", "SSE_None"), "chocolate");
//       auto dark_helix_cr  = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Helix", "SSE_Helix"), "#3A2314");
//       auto stand_cr       = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Strand", "SSE_Strand"), "#82ce9a");
//       auto brown_dna_cr   = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("NUCLEICACIDS", "NucleicAcids"), "brown");

      auto init_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("/*/*/*.*/*:*", "NoSecondary"), "grey");
      auto none_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_None", "SSE_None"), "grey");
      auto dark_helix_cr  = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Helix", "SSE_Helix"), "forestgreen");
      auto stand_cr       = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Strand", "SSE_Strand"), "firebrick");
      auto brown_dna_cr   = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("NUCLEICACIDS", "NucleicAcids"), "brown");

      // in a molecularrepresentation, set the following:
      // radiusOneHelix=0.95;
      // radiusTwoHelix=0.15;
      // radiusTwoStrand=0.3;
      

//       SolidColorRule green_mol_scr;
//       FCXXCoord green(0.4, 0.7, 0.4);
//       green_mol_scr.setColor(green);
//       std::shared_ptr<SolidColorRule> green_mol_cr = std::make_shared<SolidColorRule> (green_mol_scr);
//       std::shared_ptr<ColorScheme> green_mol_cs = std::make_shared<ColorScheme>(green_mol_cr);

      auto ss_cs = ColorScheme::colorBySecondaryScheme();

      auto        ramp_cs = ColorScheme::colorRampChainsScheme();
      auto ribbon_ramp_cs = ColorScheme::colorRampChainsScheme();
      auto ele_cs         = ColorScheme::colorByElementScheme();
      auto chains_cs      = ColorScheme::colorChainsScheme();
      auto handles = chains_cs->prepareForMMDB(atom_sel.mol);
      auto handles_2 = ramp_cs->prepareForMMDB(atom_sel.mol);

      ss_cs->addRule(init_cr);
      ss_cs->addRule(none_cr);
      ss_cs->addRule(dark_helix_cr);
      ss_cs->addRule(stand_cr);
      ramp_cs->addRule(brown_dna_cr);

      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    std::cout << "here with chain " << chain_p << std::endl;
	    std::pair<bool, std::pair<int, int> > residue_min_max = coot::util::min_max_residues_in_polymer_chain(chain_p);
	    if (residue_min_max.first) {
	       std::cout << "min: max for chain " << chain_p->GetChainID() << " "
			 << residue_min_max.second.first << " "
			 << residue_min_max.second.second << " "
			 << std::endl;
	       AtomPropertyRampColorRule apcrr;
	       apcrr.setStartValue(residue_min_max.second.first);
	       apcrr.setEndValue(residue_min_max.second.second);
	       auto apcrr_p = std::make_shared<AtomPropertyRampColorRule> (apcrr);
	       ribbon_ramp_cs->addRule(apcrr_p);
	       std::string selection_str = "//";
	       std::string chain_id = chain_p->GetChainID();
	       selection_str += chain_id;

	       FCXXCoord green(0.4, 0.7, 0.4);
	       auto green_cr = SolidColorRule::colorRuleForSelectionStringAndColor("//", green);
	       auto green_chains_cs      = ColorScheme::colorChainsScheme();
	       green_chains_cs->addRule(green_cr);

	       // VdWSurface and AccessibleSurface don't seem to work at the moment
	       //
	       // auto molrepinst_1 = MolecularRepresentationInstance::create(my_mol, ribbon_ramp_cs, selection_str, "MolecularSurface");
	       // auto molrepinst_1 = MolecularRepresentationInstance::create(my_mol, ramp_cs, selection_str, "Ribbon");
	       // auto molrepinst_2 = MolecularRepresentationInstance::create(my_mol, chains_cs, mmdb_chain, "DishyBases");
	       // auto molrepinst_2 = MolecularRepresentationInstance::create(my_mol, chains_cs, mmdb_chain, "Calpha");

	       auto molrepinst_1 = MolecularRepresentationInstance::create(my_mol, ele_cs, "//A/22-24/CA", "MolecularSurface");
	       graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(molrepinst_1);
	       molrepinsts.push_back(molrepinst_1);

	    } else {
	       std::cout << "No min/max found for chain " << chain_p << " " << chain_p->GetChainID() << std::endl;
	    }
	 }
      }
   }
#endif // HAVE_CXX11
#endif // USE_MOLECULES_TO_TRIANGLES
}

void
molecule_class_info_t::set_mol_triangles_is_displayed(int state) {
   
#ifdef USE_MOLECULES_TO_TRIANGLES
   if (molrepinsts.size()) {
      if (state) {
	 for (auto mri : molrepinsts)
	    graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(mri);
      } else {
	 for (auto mri : molrepinsts)
	    graphics_info_t::mol_tri_scene_setup->removeRepresentationInstance(mri);
      }
   }
#endif // USE_MOLECULES_TO_TRIANGLES

}

#ifdef USE_MOLECULES_TO_TRIANGLES
int
molecule_class_info_t::add_molecular_representation(const std::string &atom_selection,
						    const std::string &colour_scheme,
						    const std::string &style) {

   int status = 0;
   return status;
}
#endif // USE_MOLECULES_TO_TRIANGLES
