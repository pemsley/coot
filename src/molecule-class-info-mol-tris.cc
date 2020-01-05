

#ifdef USE_PYTHON
#include "Python.h"
#endif

#include "graphics-info.h"
#include "molecule-class-info.h"

// make and add to the scene
int
molecule_class_info_t::make_molecularrepresentationinstance(const std::string &atom_selection,
   const std::string &colour_scheme,
   const std::string &style) {

      std::cout << "-------------------- make_molecularrepresentationinstance() start " << std::endl;

      int status = -1;

      #ifdef USE_MOLECULES_TO_TRIANGLES
      if (atom_sel.mol) {

         auto my_mol = std::make_shared<MyMolecule> (atom_sel.mol);

         // (.. (seletionText, name), "colour")
         // sheet colour: 130 190 170 -> #82BEAA

         // auto init_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("/*/*/*.*/*:*", "NoSecondary"), "chocolate");
         // auto none_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_None", "SSE_None"), "chocolate");
         auto init_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("/*/*/*.*/*:*", "NoSecondary"), "#707070");
         auto none_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_None", "SSE_None"), "#707070");
         auto dark_helix_cr  = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Helix", "SSE_Helix"), "#3A2314");
         auto stand_cr       = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Strand", "SSE_Strand"), "#82ce9a");
         auto brown_dna_cr   = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("NUCLEICACIDS", "NucleicAcids"), "brown");

         if (false) {
            auto init_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("/*/*/*.*/*:*", "NoSecondary"), "grey");
            auto none_cr        = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_None", "SSE_None"), "grey");
            auto dark_helix_cr  = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Helix", "SSE_Helix"), "forestgreen");
            auto stand_cr       = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("SSE_Strand", "SSE_Strand"), "firebrick");
            auto brown_dna_cr   = SolidColorRule::colorRuleForSelectionAndName(std::make_shared<CompoundSelection>("NUCLEICACIDS", "NucleicAcids"), "brown");
         }

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
         // ramp_cs->addRule(brown_dna_cr);

         if (false) {
            std::cout << "Pre-check chains " << std::endl;
            for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
               mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
               if (model_p) {
                  int n_chains = model_p->GetNumberOfChains();
                  for (int ichain=0; ichain<n_chains; ichain++) {
                     mmdb::Chain *chain_p = model_p->GetChain(ichain);
                     int nres = chain_p->GetNumberOfResidues();
                     std::cout << "pre-check chain " << chain_p  << " which has chain id "
                     << chain_p->GetChainID() << " "
                     << chain_p->GetNumberOfResidues() << " residues " << std::endl;
                  }
               }
            }
         }

         std::cout << "----------------- make_molecularrepresentationinstance() model loop " << std::endl;

         for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
            if (model_p) {

               auto this_cs = ele_cs;
               auto this_atom_selection = atom_selection;

               // do need need to add/override things?

               if (colour_scheme == "colorRampChainsScheme") {

                  std::vector<mmdb::Chain *> chain_selection =
                  coot::util::chains_in_atom_selection(atom_sel.mol, imod, atom_selection);
                  int n_chains = model_p->GetNumberOfChains();

                  // run through the chains in the atom selection.
                  for (unsigned int ichain=0; ichain<chain_selection.size(); ichain++) {
                     mmdb::Chain *chain_p = chain_selection[ichain];

                     std::cout << "here with chain " << chain_p << " which has chain id " << chain_p->GetChainID() << std::endl;
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
                        // if we are using a ribbon_ramp_cs, then we need to apply that colour scheme only to the selected chain
                        std::string chain_atom_selection = "/" + coot::util::int_to_string(imod) + "/" + chain_p->GetChainID();

                        FCXXCoord green(0.4, 0.7, 0.4);
                        auto green_cr = SolidColorRule::colorRuleForSelectionStringAndColor(atom_selection, green);
                        auto green_chains_cs      = ColorScheme::colorChainsScheme();
                        green_chains_cs->addRule(green_cr);

                        this_atom_selection = chain_atom_selection;

                        // VdWSurface and AccessibleSurface don't seem to work at the moment
                        //
                        // auto molrepinst_1 = MolecularRepresentationInstance::create(my_mol, ribbon_ramp_cs, selection_str, "MolecularSurface");
                        // auto molrepinst_1 = MolecularRepresentationInstance::create(my_mol, ramp_cs, selection_str, "Ribbon");
                        // auto molrepinst_2 = MolecularRepresentationInstance::create(my_mol, chains_cs, mmdb_chain, "DishyBases");
                        // auto molrepinst_2 = MolecularRepresentationInstance::create(my_mol, chains_cs, mmdb_chain, "Calpha");
                        // auto molrepinst_1 = MolecularRepresentationInstance::create(my_mol, ss_cs, selection_str, style);

                     }
                  }
               }

               if (colour_scheme == "colorChainsScheme")
               this_cs = chains_cs;
               if (colour_scheme == "colorRampChainsScheme") {
                  this_cs = ribbon_ramp_cs;
               }
               if (colour_scheme == "colorByElementScheme") {
                  this_cs = ele_cs;
               }
               if (colour_scheme == "colorBySecondaryScheme") {
                  this_cs = ss_cs;
               }

               std::cout << "###### calling create() with this_atom_selection " << atom_selection << std::endl;
               GLenum err = glGetError(); if (true) std::cout << "in " << __FUNCTION__ << " A with err " << err << std::endl;
               auto molrepinst_1 = MolecularRepresentationInstance::create(my_mol, this_cs, atom_selection, style);
               err = glGetError(); if (true) std::cout << "in " << __FUNCTION__ << " B with err " << err << std::endl;
               if (!graphics_info_t::mol_tri_scene_setup) {
                  std::cout << "ERROR:: null mol_tri_scene_setup" << std::endl;
               } else {
                  graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(molrepinst_1);
                  err = glGetError(); if (true) std::cout << "in " << __FUNCTION__ << " C with err " << err << std::endl;
                  // auto molrepinst_2 = MolecularRepresentationInstance::create(my_mol, ss_cs, "//A", "MolecularSurface");
                  // graphics_info_t::mol_tri_scene_setup->addRepresentationInstance(molrepinst_2);

                  molrepinsts.push_back(molrepinst_1);
                  status = molrepinsts.size() -1; // index of back
               }
            }
         }
      }

      #endif // USE_MOLECULES_TO_TRIANGLES

      return status;
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

   // return index in molrepinsts
   int status = make_molecularrepresentationinstance(atom_selection, colour_scheme, style);
   return status;
}

#endif // USE_MOLECULES_TO_TRIANGLES


#ifdef USE_MOLECULES_TO_TRIANGLES
void
molecule_class_info_t::remove_molecular_representation(int idx) {

   if (idx >= 0) {

      // this will shuffle the indices of the other molecule representations, hmm...
      // molrepinsts.erase();
      if (molrepinsts.size() > 0) {
    std::vector<std::shared_ptr<MolecularRepresentationInstance> >::const_iterator it = molrepinsts.end();
    it --;
    molrepinsts.erase(it);
    std::cout << "erased - now molrepinsts size " << molrepinsts.size() << std::endl;
      }
   }
}

#endif // USE_MOLECULES_TO_TRIANGLES
