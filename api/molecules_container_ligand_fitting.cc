
#include "ligand/wligand.hh"

#include "molecules_container.hh"

#include <clipper/ccp4/ccp4_map_io.h> // debugging mapout


   //! Ligand Fit
   //! @return a vector or the best fitting ligands to this blob.
   //! I am not yet clear what extra cut-offs and flags need to be added here.
std::vector<int>
molecules_container_t::fit_ligand_right_here(int imol_protein, int imol_map, int imol_ligand, float x, float y, float z,
                                             float n_rmsd,
                                             bool use_conformers, unsigned int n_conformers) {

   auto get_first_residue_name = [] (mmdb::Manager *mol) {
      std::string res_name;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  res_name = residue_p->GetResName();
                  break;
               }
            }
            if (! res_name.empty()) break;
         }
      }
      return res_name;
   };

   std::vector<int> mol_list;

   if (is_valid_model_molecule(imol_protein)) {
      if (is_valid_model_molecule(imol_ligand)) {
         if (is_valid_map_molecule(imol_map)) {
            clipper::Coord_orth pt(x,y,z);

            coot::wligand wlig;
            wlig.set_verbose_reporting();
            wlig.set_debug_wiggly_ligands();
            
            try {
               coot::minimol::molecule mmol(molecules[imol_ligand].atom_sel.mol);
               std::string res_name = get_first_residue_name(molecules[imol_ligand].atom_sel.mol);
#ifdef EMSCRIPTEN
               int n_threads = 3;
#else
               int n_threads = coot::get_max_number_of_threads();
               n_threads = n_threads - 2;
               if (n_threads < 1) n_threads = 1;
#endif
               ctpl::thread_pool thread_pool(n_threads);
               if (use_conformers) {
                  bool optim_geom = true;
                  bool fill_return_vec = false; // would give input molecules (not conformers)
                  wlig.install_simple_wiggly_ligands(&geom, mmol, imol_ligand,
                                                     n_conformers, optim_geom,
                                                     fill_return_vec, &thread_pool, n_threads);
               } else {

                  mmdb::Manager *ligand_mol = molecules[imol_ligand].atom_sel.mol;
                  wlig.install_ligand(ligand_mol);
               }

               clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
               wlig.import_map_from(xmap);
               short int mask_waters_flag = true;
               wlig.set_map_atom_mask_radius(2.0);  // Angstroms
               mmdb::Manager *protein_mol = molecules[imol_protein].atom_sel.mol;
               wlig.mask_map(protein_mol, mask_waters_flag);

               wlig.cluster_from_point(pt, n_rmsd);
               wlig.fit_ligands_to_clusters(1); // just this cluster.
                  
               // now add in the solution ligands:
               int n_final_ligands = wlig.n_clusters_final();

               if (n_final_ligands == 1) {
                  unsigned int iclust = 0;
                  unsigned int isol   = 0;
                  coot::minimol::molecule m = wlig.get_solution(isol, iclust);

                  if (false)
                     std::cout << "########## fit_ligand_right_here(): m has " << m.count_atoms() << " atoms " << std::endl;
                  mmdb::Manager *ligand_mol = m.pcmmdbmanager();
                  
                  coot::hetify_residues_as_needed(ligand_mol);
                  atom_selection_container_t asc = make_asc(ligand_mol);
                  int imol_in_hope = molecules.size();
                  std::string name = "Fitted ligand " + res_name;
                  coot::molecule_t mm(asc, imol_in_hope, name);
                  molecules.push_back(mm);
                  mol_list.push_back(imol_in_hope);

                  if (false) {
                     std::string file_name = "debug.map";
                     clipper::Xmap<float> masked_map = wlig.masked_map();
                     clipper::CCP4MAPfile mapout;
                     mapout.open_write(file_name);
                     mapout.export_xmap(xmap);
                     mapout.close_write();
                  }

               }
            }
            catch (const std::runtime_error &mess) {
               std::cout << "ERROR:: in flexible ligand definition.\n";
               std::cout << mess.what() << std::endl;
            }

         }
      }
   }
   return mol_list;
}


//! "Jiggle-Fit Ligand"
//! if n_trials is 0, then a sensible default value will be used.
//! if translation_scale_factor is negative then a sensible default value will be used.
//! @return a value less than -99.9 on failure to fit.
float
molecules_container_t::fit_to_map_by_random_jiggle(int imol, const coot::residue_spec_t &res_spec, int n_trials, float translation_scale_factor) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (n_trials == 0) n_trials = 100;
      if (translation_scale_factor < 0.0) translation_scale_factor = 1.0;
      if (is_valid_map_molecule(imol_refinement_map)) {
         clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
         float rmsd = molecules[imol_refinement_map].get_map_rmsd_approx();
         float score = molecules[imol].fit_to_map_by_random_jiggle(res_spec, xmap, rmsd, n_trials, translation_scale_factor);
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;


}

//! "Jiggle-Fit Ligand" with different interface
//! @return a value less than -99.9 on failure to fit.
float
molecules_container_t::fit_to_map_by_random_jiggle_using_cid(int imol, const std::string &cid, int n_trials, float translation_scale_factor) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t spec = residue_cid_to_residue_spec(imol, cid);
      if (is_valid_map_molecule(imol_refinement_map)) {
         clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
         float rmsd = molecules[imol_refinement_map].get_map_rmsd_approx();
         status = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, translation_scale_factor);
      } else {
         std::cout << "ERROR:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_refinement_map << std::endl;
      }
   } else {
      std::cout << "ERROR:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#endif

#include "lidia-core/svg-molecule.hh"


//! This is a ligand function, not really a ligand-fitting function.
//!
//! But more importantly than that, it doesn't work yet.
std::string
molecules_container_t::get_svg_for_residue_type(int imol, const std::string &comp_id,
                                                bool dark_bg_flag) {

   std::string s = "Needs-to-be-compiled-with-the-RDKit";

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   std::map<std::string, std::string>::const_iterator it = ligand_svg_store.find(comp_id);
   if (it != ligand_svg_store.end()) {

      return it->second;

   } else {

      std::pair<bool, coot::dictionary_residue_restraints_t> mr = geom.get_monomer_restraints(comp_id, imol);
      if (mr.first) {
         try {
            const coot::dictionary_residue_restraints_t &restraints = mr.second;
            svg_molecule_t svg;
            RDKit::RWMol mol = coot::rdkit_mol(restraints);
            // bool undelocalize_flag = true;
            // used undelocalize_flag in RDKit::RWMol mol_rw = coot::rdkit_mol(r, rest, "", undelocalize_flag);
            // RDKit::RWMol mol_rw = coot::rdkit_mol_sanitized(r, imol, geom);
            RDKit::MolOps::removeHs(mol);
            RDKit::MolOps::Kekulize(mol);
            int iconf = RDDepict::compute2DCoords(mol, NULL, true);
            RDKit::Conformer &conf = mol.getConformer(iconf);
            RDKit::WedgeMolBonds(mol, &conf);
            svg.import_rdkit_mol(&mol, iconf);
            s = svg.render_to_svg_string(dark_bg_flag);
            ligand_svg_store[comp_id] = s;
         }
         catch (const Invar::Invariant &e) {
            std::cout << "error " << e.what() << std::endl;
         }

      } else {
         s = std::string("No dictionary for ") + comp_id;
      }
   }

#endif // MAKE_ENHANCED_LIGAND_TOOLS

   return s;
}
