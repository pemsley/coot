
#include <clipper/ccp4/ccp4_map_io.h> // debugging mapout
#include <stdexcept>

#include "utils/base64-encode-decode.hh"
#include "ligand/wligand.hh"

#include "molecules-container.hh"


// Give this ex-lambda function a home?
std::string get_first_residue_name(mmdb::Manager *mol) {

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
}

#if 0
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
// Another ligand but non-ligand-fitting function:
RDKit::RWMol molecules_container_t::get_rdkit_mol(const std::string &residue_name, int imol_enc) {

   RDKit::RWMol m;

   try {
      std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
         geom.get_monomer_restraints(residue_name, imol_enc);
      if (r_p.first) {
         const auto &restraints = r_p.second;
         m = coot::rdkit_mol(restraints);
         std::string prop_string = "ligand-from-dictionary-" + residue_name + "-" + std::to_string(imol_enc);
         m.setProp("moorhen-id", prop_string);
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
   }
   return m;
}
#endif
#endif


#ifdef MAKE_ENHANCED_LIGAND_TOOLS

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

#endif //RD_MOLPICKLE_H


std::string
molecules_container_t::get_rdkit_mol_pickle_base64(const std::string &residue_name, int imol_enc) {

   RDKIT_GRAPHMOL_EXPORT RDKit::MolPickler mp;
   std::string pickle_string;
   try {
      std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
         geom.get_monomer_restraints(residue_name, imol_enc);
      if (r_p.first) {
         const auto &restraints = r_p.second;
         RDKit::RWMol mol = coot::rdkit_mol(restraints);
         if (mol.getNumAtoms() > 0) {
            mp.pickleMol(mol, pickle_string);
            return moorhen_base64::base64_encode((const unsigned char*)pickle_string.c_str(), pickle_string.size());
         }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "DEBUG:: get_rdkit_mol_pickle_base64 " << rte.what() << std::endl;
   }
   return pickle_string;
}
#endif


//! Ligand Fit
//! @return a vector or the best fitting ligands to this blob.
//! I am not yet clear what extra cut-offs and flags need to be added here.
std::vector<int>
molecules_container_t::fit_ligand_right_here(int imol_protein, int imol_map, int imol_ligand, float x, float y, float z,
                                             float n_rmsd,
                                             bool use_conformers, unsigned int n_conformers) {


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

#ifdef LETS_USE_THE_THREAD_POOL

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
#else

               if (use_conformers) {
                  // bool optim_geom = true;
                  bool optim_geom = false;
                  for (unsigned int i_conf=0; i_conf<n_conformers; i_conf++) {
                     int imol_enc = coot::protein_geometry::IMOL_ENC_ANY; // pass this
                     wlig.install_simple_wiggly_ligand(&geom, mmol, imol_ligand, i_conf, optim_geom);
                  }

               } else {

                  mmdb::Manager *ligand_mol = molecules[imol_ligand].atom_sel.mol;
                  wlig.install_ligand(ligand_mol);
               }
#endif

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
                  int n_atoms = m.count_atoms();

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


//! Ligand Fitting
//!
//! @return a vector of indices of molecules for the best fitting ligands to this blob.
std::vector<molecules_container_t::fit_ligand_info_t>
molecules_container_t::fit_ligand(int imol_protein, int imol_map, int imol_ligand,
                                  float n_rmsd, bool use_conformers, unsigned int n_conformers) {

   std::vector<fit_ligand_info_t> mol_list;

   if (is_valid_model_molecule(imol_protein)) {
      if (is_valid_model_molecule(imol_ligand)) {
         if (is_valid_map_molecule(imol_map)) {

            coot::wligand wlig;
            // wlig.set_verbose_reporting();
            // wlig.set_debug_wiggly_ligands();

            // 20231121-PE No thread pool. Add it here if needed
            try {

               coot::minimol::molecule mmol(molecules[imol_ligand].atom_sel.mol);
               std::string res_name = get_first_residue_name(molecules[imol_ligand].atom_sel.mol);
               if (use_conformers) {
                  // bool optim_geom = true;
                  bool optim_geom = false;
                  for (unsigned int i_conf=0; i_conf<n_conformers; i_conf++) {
                     // std::cout << "installing ligand " << i_conf << std::endl;
                     wlig.install_simple_wiggly_ligand(&geom, mmol, imol_ligand, i_conf, optim_geom);
                  }
               } else {
                  mmdb::Manager *ligand_mol = molecules[imol_ligand].atom_sel.mol;
                  wlig.install_ligand(ligand_mol);
               }
               clipper::Xmap<float> &xmap = molecules[imol_map].xmap;

               {
                  std::string map_file_name_1 = "molecules-container-fit-ligand-A.map";
                  clipper::CCP4MAPfile mapout;
                  mapout.open_write(map_file_name_1);
                  mapout.export_xmap(xmap);
                  mapout.close_write();
               }

               wlig.import_map_from(xmap);

               {
                  std::string map_file_name_1 = "molecules-container-fit-ligand-B.map";
                  clipper::CCP4MAPfile mapout;
                  mapout.open_write(map_file_name_1);
                  mapout.export_xmap(xmap);
                  mapout.close_write();
               }

               {
                  wlig.output_map("wlig.output_map.map");
               }

               short int mask_waters_flag = true;
               wlig.set_map_atom_mask_radius(2.0);  // Angstroms
               mmdb::Manager *protein_mol = molecules[imol_protein].atom_sel.mol;
               wlig.mask_map(protein_mol, mask_waters_flag);

               float ligand_acceptable_fit_fraction = 0.85; // was 0.75
               int find_ligand_n_top_ligands = 10;

               wlig.find_clusters(n_rmsd);  // trashes the xmap
               wlig.set_acceptable_fit_fraction(ligand_acceptable_fit_fraction);
               wlig.fit_ligands_to_clusters(find_ligand_n_top_ligands); // 10 clusters

               // now add in the solution ligands: 20231121-PE (What did I mean by this?)
               int n_clusters = wlig.n_clusters_final();

               coot::minimol::molecule m;
               for (int iclust=0; iclust<n_clusters; iclust++) {

                  if (wlig.n_ligands_for_cluster(iclust) > 0) {

                     float frac_lim = 0.9;
                     float correl_frac_lim = 0.9;
                     bool find_ligand_multiple_solutions_per_cluster_flag = true;
                     float cv = wlig.get_cluster_volume(iclust);

                     // nino-mode
                     unsigned int n_ligands_for_cluster = wlig.n_ligands_for_cluster(iclust, frac_lim);
                     wlig.score_and_resort_using_correlation(iclust, n_ligands_for_cluster);

                     if (find_ligand_multiple_solutions_per_cluster_flag == false) {
                        n_ligands_for_cluster = 1;
                        correl_frac_lim = 0.975;
                     }

                     if (n_ligands_for_cluster > 12) n_ligands_for_cluster = 12; // arbitrary limit of max 12 solutions per cluster
                     float tolerance = 20.0;
                     // limit_solutions should be run only after a post-correlation sort.
                     //
                     wlig.limit_solutions(iclust, correl_frac_lim, n_ligands_for_cluster, tolerance, true);

                     for (unsigned int isol=0; isol<n_ligands_for_cluster; isol++) {
                        m = wlig.get_solution(isol, iclust);
                        if (! m.is_empty()) {
                           // std::cout << "------------------ found a solution " << isol << " for iclust " << iclust << std::endl;
                           coot::minimol::molecule m = wlig.get_solution(isol, iclust);
                           mmdb::Manager *ligand_mol = m.pcmmdbmanager();
                           coot::hetify_residues_as_needed(ligand_mol);
                           atom_selection_container_t asc = make_asc(ligand_mol);
                           int imol_in_hope = molecules.size();
                           std::string name = "Fitted ligand " + res_name;
                           coot::molecule_t mm(asc, imol_in_hope, name);
                           molecules.push_back(mm);
                           fit_ligand_info_t fli(imol_in_hope, iclust, isol);
                           // set fli fitting score and cluster volume
                           fli.cluster_volume = cv;
                           coot::ligand_score_card lsc = wlig.get_solution_score(iclust, isol);
                           fli.fitting_score = lsc.correlation.second;
                           mol_list.push_back(fli);
                        }
                     }
                  }
               }
            }
            catch (const std::runtime_error &e) {
               std::cout << "WARNING::" << e.what() << std::endl;
            }
         }
      }
   }

   return mol_list;
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
         float tf_1 = translation_scale_factor;
         float tf_2 = 0.5 * translation_scale_factor;
         float tf_3 = 0.25 * translation_scale_factor;
         float tf_4 = 0.18 * translation_scale_factor;
         status = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, tf_1);
         status = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, tf_2);
         status = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, tf_3);
         status = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, tf_4);
      } else {
         std::cout << "ERROR:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_refinement_map << std::endl;
      }
   } else {
      std::cout << "ERROR:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}

float
molecules_container_t::fit_to_map_by_random_jiggle_with_blur_using_cid(int imol, int imol_map, const std::string &cid, float b_factor,
                                                                       int n_trials, float translation_scale_factor) {

   float r = -999.0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         // make blurred copy of imol_map's map and use that for fitting
         clipper::Xmap<float> xmap = molecules[imol_map].xmap;
         coot::util::sharpen_blur_map(&xmap, b_factor);
         auto p = coot::util::mean_and_variance(xmap);
         float rmsd = std::sqrt(p.second);
         float tf_1 = translation_scale_factor;
         float tf_2 = 0.5 * translation_scale_factor;
         float tf_3 = 0.25 * translation_scale_factor;
         float tf_4 = 0.18 * translation_scale_factor;
         r = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, tf_1);
         r = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, tf_2);
         r = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, tf_3);
         r = molecules[imol].fit_to_map_by_random_jiggle_using_atom_selection(cid, xmap, rmsd, n_trials, tf_4);
      } else {
         std::cout << "WARNING:: " << imol_map << " is not a valid map"<< std::endl;
      }
   } else {
      std::cout << "WARNING:: " << imol_map << " is not a valid model"<< std::endl;
   }
   return r;
}

// move this into coot-utils, I think - maybe fast-eigens get_fast_eigenvalues_for_residue_atoms()
#include "coot-utils/fast-eigens.hh"
std::vector<double>
get_eigenvalues(mmdb::Residue *residue_p) {

   // c.f. coot::distortion_score_plane_internal()
   std::vector<double> v;
   std::vector<double> x, y, z;
   if (residue_p) {
      int n_atoms = residue_p->GetNumberOfAtoms();
      for (int iat=0; iat<n_atoms; iat++) {
         mmdb::Atom *at = residue_p->GetAtom(iat);
         if (! at->isTer()) {
            x.push_back(at->x);
            y.push_back(at->y);
            z.push_back(at->z);
         }
      }
      if (! x.empty()) {
         coot::stats::single stats_x(x);
         coot::stats::single stats_y(y);
         coot::stats::single stats_z(z);
         double x_mean = stats_x.mean();
         double y_mean = stats_y.mean();
         double z_mean = stats_z.mean();
         clipper::Matrix<double> mat(3,3);
         for (int iat=0; iat<n_atoms; iat++) {
            mmdb::Atom *at = residue_p->GetAtom(iat);
            if (! at->isTer()) {
               mat(0,0) += (double(at->x) - x_mean) * (double(at->x) - x_mean);
               mat(1,1) += (double(at->y) - y_mean) * (double(at->y) - y_mean);
               mat(2,2) += (double(at->z) - z_mean) * (double(at->z) - z_mean);
               mat(0,1) += (double(at->x) - x_mean) * (double(at->y) - y_mean);
               mat(0,2) += (double(at->x) - x_mean) * (double(at->z) - z_mean);
               mat(1,2) += (double(at->y) - y_mean) * (double(at->z) - z_mean);
            }
         }
         mat(1,0) = mat(0,1);
         mat(2,0) = mat(0,2);
         mat(2,1) = mat(1,2);
         std::tuple<double, double, double> eigens = coot::fast_eigens(mat, false);
         // std::cout << "eigens: " << std::get<0>(eigens) << " " << std::get<1>(eigens) << " " << std::get<2>(eigens) << std::endl;
         v.push_back(std::get<0>(eigens));
         v.push_back(std::get<1>(eigens));
         v.push_back(std::get<2>(eigens));
      }
   }
   return v;
}

std::vector<double>
molecules_container_t::get_eigenvalues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   std::vector<double> v;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
      mmdb::Residue *r = molecules[imol].get_residue(residue_spec);
      if (r) {
         v = ::get_eigenvalues(r);
      } else {
         std::cout << "WARNING:: get_eigenvalues(): No residue " << chain_id << " " << res_no
                   << " in molecule " << imol << std::endl;
      }
   }
   return v;
}


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include "lidia-core/svg-molecule.hh"
#include "svg-store-key.hh"
#endif

//! This is a ligand function, not really a ligand-fitting function.
//!
std::string
molecules_container_t::get_svg_for_residue_type(int imol, const std::string &comp_id,
                                                bool use_rdkit_svg,
                                                const std::string &background_type) {

   std::string s = "Needs-to-be-compiled-with-the-RDKit"; // gets changed to the correct string (hopefully)

   bool dark_bg_flag = false;
   // "dark-mode" representation
   if (background_type == "light-bonds/transparent-bg") dark_bg_flag = true;
   if (background_type == "light-bonds/opaque-bg")      dark_bg_flag = true;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

   svg_store_key_t key(imol, comp_id);
   std::map<svg_store_key_t, std::string>::const_iterator it = ligand_svg_store.find(key);
   if (it != ligand_svg_store.end()) {

      return it->second;

   } else {

      std::pair<bool, coot::dictionary_residue_restraints_t> mr = geom.get_monomer_restraints(comp_id, imol);
      if (mr.first) {

         try {
            const coot::dictionary_residue_restraints_t &restraints = mr.second;
            std::pair<int, RDKit::RWMol> mol_pair = coot::rdkit_mol_with_2d_depiction(restraints);
            int conformer_id = mol_pair.first;

            if (conformer_id == -1) { // there was no depiction in the restraints

               RDKit::RWMol mol = coot::rdkit_mol(restraints);
               // bool undelocalize_flag = true;
               // used undelocalize_flag in RDKit::RWMol mol_rw = coot::rdkit_mol(r, rest, "", undelocalize_flag);
               // RDKit::RWMol mol_rw = coot::rdkit_mol_sanitized(r, imol, geom);
               RDKit::MolOps::removeHs(mol);
               RDKit::MolOps::Kekulize(mol);
               int iconf = RDDepict::compute2DCoords(mol, NULL, true);
               RDKit::Conformer &conf = mol.getConformer(iconf);
               RDKit::WedgeMolBonds(mol, &conf);
               if (use_rdkit_svg) {
                  std::stringstream ss;
                  RDKit::MolDraw2DSVG svg_drawer_ss(400, 400, ss);
                  svg_drawer_ss.drawMolecule(mol);
                  svg_drawer_ss.finishDrawing();
                  s = ss.str();
               } else {
                  svg_molecule_t svg;
                  svg.import_rdkit_mol(&mol, iconf);
                  double sf = 400.0;
                  bool add_bg = true;
                  if (background_type == "light-bonds/transparent-bg") add_bg = false;
                  if (background_type == "dark-bonds/transparent-bg")  add_bg = false;
                  s = svg.render_to_svg_string(sf, dark_bg_flag, add_bg);
                  ligand_svg_store[key] = s;
               }
            } else {

               // std::cout << "Use the depiction from the restraints ---------------------------- " << std::endl;

               // there was a depiction in the restraints - use that.
               auto &rdkit_mol(mol_pair.second);

               if (use_rdkit_svg) {

                  if (false) { // this resets the wedges, which we don't want to do
                     RDKit::Conformer &conf = rdkit_mol.getConformer(conformer_id);
                     RDKit::WedgeMolBonds(rdkit_mol, &conf);
                  }
                  std::stringstream ss;
                  RDKit::MolDraw2DSVG svg_drawer_ss(400, 400, ss);
                  // svg_drawer_ss.drawMolecule(rdkit_mol);
                  svg_drawer_ss.drawMolecule(rdkit_mol, nullptr, nullptr, nullptr, conformer_id);
                  std::string dt = svg_drawer_ss.getDrawingText();
                  svg_drawer_ss.finishDrawing();
                  s = ss.str();

               } else {

                  RDKit::Conformer &conf = rdkit_mol.getConformer(conformer_id);
                  RDKit::WedgeMolBonds(rdkit_mol, &conf);
                  svg_molecule_t svg;
                  svg.import_rdkit_mol(&rdkit_mol, conformer_id);
                  double sf = 400.0;
                  bool add_background_rect = true;
                  if (background_type == "light-bonds/transparent-bg") add_background_rect = false;
                  if (background_type == "dark-bonds/transparent-bg")  add_background_rect = false;
                  s = svg.render_to_svg_string(sf, dark_bg_flag, add_background_rect);
               }
               ligand_svg_store[key] = s;
            }
         }
         catch (const Invar::Invariant &e) {
            std::cout << "error " << e.what() << std::endl;
         }
         catch (const std::runtime_error &e) {
            std::cout << "error " << e.what() << std::endl;
         }
         catch (...) {
            std::cout << "error something" << std::endl;
         }

      } else {
         s = std::string("No dictionary for ") + comp_id;
      }
   }

#endif // MAKE_ENHANCED_LIGAND_TOOLS

   return s;
}


//! Fit ligands
//! ``multi_ligand_molecule_number_list`` is a colon-separated list of molecules, *e.g.* "2:3:4"
std::vector<molecules_container_t::fit_ligand_info_t>
molecules_container_t::fit_ligand_multi_ligand(int imol_protein, int imol_map, const std::string &multi_ligand_molecule_number_list,
                                               float n_rmsd, bool use_conformers, unsigned int n_conformers) {

   std::vector<fit_ligand_info_t> v;
   // now decode multi_ligand_molecule_number_list
   std::vector<int> ligand_imols;

   return v;
}
