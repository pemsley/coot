
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
      auto green_helix_cr = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Helix", "SSE_Helix"), "lightgreen");
      auto red_stand_cr = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Strand", "SSE_Strand"), "indianred");
      auto brown_dna_cr = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("NUCLEICACIDS", "NucleicAcids"), "lightseagreen");

      auto cs_ss = ColorScheme::colorBySecondaryScheme();

      auto        ramp_cs = ColorScheme::colorRampChainsScheme();
      auto ribbon_ramp_cs = ColorScheme::colorRampChainsScheme();
      auto ele_cs         = ColorScheme::colorByElementScheme();
      auto chains_cs      = ColorScheme::colorChainsScheme();
      auto handles = chains_cs->prepareForMMDB(atom_sel.mol);
      auto handles_2 = ramp_cs->prepareForMMDB(atom_sel.mol);

      // cs_ss->addRule(green_helix_cr);
      // cs_ss->addRule(red_stand_cr);
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
	       std::string mmdb_chain = "//";
	       std::string chain_id = chain_p->GetChainID();
	       mmdb_chain += chain_id;
	       auto molrepinst_1 = MolecularRepresentationInstance::create(my_mol,     cs_ss, mmdb_chain, "MolecularSurface");
	       // auto molrepinst_2 = MolecularRepresentationInstance::create(my_mol, chains_cs, mmdb_chain, "DishyBases");
	       auto molrepinst_2 = MolecularRepresentationInstance::create(my_mol, chains_cs, mmdb_chain, "Calpha");

	       if (chain_id == "C") {
		  auto molrepinst_c = MolecularRepresentationInstance::create(my_mol, cs_ss, mmdb_chain, "VdWSurface");
		  graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(molrepinst_c);
		  molrepinsts.push_back(molrepinst_c);
	       } else {
		  graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(molrepinst_1);
		  graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(molrepinst_2);
	       }
	       molrepinsts.push_back(molrepinst_1);
	       molrepinsts.push_back(molrepinst_2);
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
