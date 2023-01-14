
#include <iomanip>

#include <sys/types.h> // for stating
#include <sys/stat.h>

#include "molecules_container.hh"
#include "ideal/pepflip.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"

#include "coords/Bond_lines.h"

#include "coot-utils/oct.hh"

// statics
std::atomic<bool> molecules_container_t::restraints_lock(false);
std::atomic<bool> molecules_container_t::on_going_updating_map_lock(false);
ctpl::thread_pool molecules_container_t::static_thread_pool(4); // or so
std::string molecules_container_t::restraints_locking_function_name; // I don't know why this needs to be static
std::vector<atom_pull_info_t> molecules_container_t::atom_pulls;
// 20221018-PE not sure that this needs to be static.
clipper::Xmap<float> *molecules_container_t::dummy_xmap = new clipper::Xmap<float>;
bool molecules_container_t::make_backups_flag = true;

bool
molecules_container_t::is_valid_model_molecule(int imol) const {
   bool status = false;
   if (imol >= 0) {
      int ms = molecules.size();
      if (imol < ms) {
         status = molecules[imol].is_valid_model_molecule();
      }
   }
   return status;
}

bool
molecules_container_t::is_valid_map_molecule(int imol) const {

   bool status = false;
   if (imol >= 0) {
      int ms = molecules.size();
      if (imol < ms) {
         status = molecules[imol].is_valid_map_molecule();
      }
   }
   return status;
}

int
molecules_container_t::close_molecule(int imol) {

   int status = 0;
   int ms = molecules.size(); // type change
   if (imol < ms) {
      if (imol >= 0) {
         molecules[imol].close_yourself();
         status = 1;
      }
   }
   return status;
}

std::string
molecules_container_t::get_molecule_name(int imol) const {

   int ms = molecules.size();
   if (imol < ms)
      if (imol >= 0)
         return molecules[imol].get_name();
   return std::string("");
}


void
molecules_container_t::display_molecule_names_table() const {

   for (unsigned int imol=0; imol<molecules.size(); imol++) {
      std::cout << imol << " " << std::setw(40) << molecules[imol].get_name() << std::endl;
   }
}

void
molecules_container_t::set_draw_missing_residue_loops(bool state) {
   draw_missing_residue_loops_flag = state;
}

void
molecules_container_t::read_standard_residues() {

   std::string standard_env_dir = "COOT_STANDARD_RESIDUES";

   const char *filename = getenv(standard_env_dir.c_str());
   if (! filename) {

      std::string dir = coot::package_data_dir();
      std::string standard_file_name = coot::util::append_dir_file(dir, "standard-residues.pdb");

      struct stat buf;
      int status = stat(standard_file_name.c_str(), &buf);
      if (status != 0) { // standard-residues file was not found in
	                 // default location either...
	 std::cout << "WARNING: Can't find standard residues file in the "
		   << "default location \n";
	 std::cout << "         and environment variable for standard residues ";
	 std::cout << standard_env_dir << "\n";
	 std::cout << "         is not set.";
	 std::cout << " Mutations will not be possible\n";
	 // mark as not read then:
	 standard_residues_asc.read_success = 0;
	 standard_residues_asc.n_selected_atoms = 0;
	 // std::cout << "DEBUG:: standard_residues_asc marked as empty" << std::endl;
      } else {
	 // stat success:
	 standard_residues_asc = get_atom_selection(standard_file_name, true, false, false);
      }
   } else { 
      standard_residues_asc = get_atom_selection(filename, true, false, false);
   }
}



//! get the active atom
std::pair<int, std::string>
molecules_container_t::get_active_atom(float x, float y, float z, const std::string &displayed_model_molecules_list) const {

   auto atom_to_cid = [] (mmdb::Atom *at) {
      if (at) {
         std::string s = "/";
         s += std::to_string(at->GetModelNum());
         s += "/";
         s += std::to_string(at->GetSeqNum());
         s += std::string(at->GetInsCode());
         s += "/";
         s += std::string(at->GetAtomName());
         s += std::string(at->altLoc);
         return s;
      } else {
         return std::string("");
      }
   };

   int imol_closest = -1;
   std::string cid;
   std::vector<std::string> number_strings = coot::util::split_string(displayed_model_molecules_list, ":");
   std::vector<int> mols;
   for (const auto &item : number_strings) {
      int idx = coot::util::string_to_int(item);
      if (is_valid_model_molecule(idx))
         mols.push_back(idx);
   }

   float best_distance_sqrd = 99999999999999999.9;
   int best_imol = -1;
   mmdb::Atom *best_atom = 0;
   coot::Cartesian screen_centre(x,y,z);
   for (unsigned int ii=0; ii<mols.size(); ii++) {
      int imol = mols[ii];
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      if (mol) {
         for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int n_res = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<n_res; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        int n_atoms = residue_p->GetNumberOfAtoms();
                        for (int iat=0; iat<n_atoms; iat++) {
                           mmdb::Atom *at = residue_p->GetAtom(iat);
                           if (! at->isTer()) {
                              coot::Cartesian atom_pos(at->x, at->y, at->z);
                              float dd = coot::Cartesian::lengthsq(screen_centre, atom_pos);
                              if (dd < best_distance_sqrd) {
                                 best_distance_sqrd = dd;
                                 best_imol = imol;
                                 best_atom = at;
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (best_atom) {
      imol_closest = best_imol;
      cid = atom_to_cid(best_atom);
   }
   return std::make_pair(imol_closest, cid);
}


void
molecules_container_t::update_updating_maps(int imol) {

   // should I return the stats?

   if (updating_maps_info.imol_model == imol) {
      if (updating_maps_info.maps_need_an_update) {
         if (is_valid_model_molecule(imol)) {
            if (is_valid_map_molecule(updating_maps_info.imol_2fofc)) {
               if (is_valid_map_molecule(updating_maps_info.imol_fofc)) {
                  coot::util::sfcalc_genmap_stats_t stats =
                     sfcalc_genmaps_using_bulk_solvent(imol,
                                                       updating_maps_info.imol_2fofc,
                                                       updating_maps_info.imol_fofc,
                                                       updating_maps_info.imol_with_data_info_attached);
                  updating_maps_info.maps_need_an_update = false;
               }
            }
         }
      } else {
         // 20221106-PE add debugging for now
         std::cout << "in updating_maps_info() maps_need_an_update is false" << std::endl;
      }
   }
}

void
molecules_container_t::set_updating_maps_need_an_update(int imol) {

   if (updating_maps_info.imol_model == imol)
      updating_maps_info.maps_need_an_update = true;

}

coot::atom_spec_t
molecules_container_t::atom_cid_to_atom_spec(int imol, const std::string &cid) const {

   coot::atom_spec_t spec;
   if (is_valid_model_molecule(imol)) {
      auto p = molecules[imol].cid_to_atom_spec(cid);
      if (p.first) {
         spec = p.second;
      } else {
         std::cout << "WARNING:: molecule_class_info_t::atom_cid_to_atom_spec() no matching atom " << cid << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return spec;
}


coot::residue_spec_t
molecules_container_t::residue_cid_to_residue_spec(int imol, const std::string &cid) const   {

   coot::residue_spec_t spec;

   if (is_valid_model_molecule(imol)) {
      auto p = molecules[imol].cid_to_residue_spec(cid);
      if (p.first) {
         spec = p.second;
      } else {
         std::cout << "WARNING:: molecule_class_info_t::residue_cid_to_residue_spec() no matching residue " << cid << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return spec;
}


void
molecules_container_t::geometry_init_standard() {
   geom.init_standard();
}


int
molecules_container_t::undo(int imol) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].undo();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::redo(int imol) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].redo();
         } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}




int
molecules_container_t::flip_peptide(int imol, const coot::atom_spec_t &as, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      result = molecules[imol].flip_peptide(as, alt_conf);
      set_updating_maps_need_an_update(imol);
   }
   return result;
}

int
molecules_container_t::flip_peptide_using_cid(int imol, const std::string &atom_cid, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      auto &m = molecules[imol];
      std::pair<bool, coot::atom_spec_t> as = m.cid_to_atom_spec(atom_cid);
      if (as.first) {
         const auto &atom_spec = as.second;
         result = molecules[imol].flip_peptide(atom_spec, alt_conf); // N check in here
         set_updating_maps_need_an_update(imol);
      }
   }
   return result;
}


int
molecules_container_t::read_pdb(const std::string &file_name) {

   int status = -1;
   atom_selection_container_t asc = get_atom_selection(file_name);
   if (asc.read_success) {

      if (false) { // for debugging atom selections in the future
         int udd_mess_indexing_1 = asc.mol->RegisterUDInteger(mmdb::UDR_ATOM, "xx");
         int udd_mess_indexing_2 = asc.mol->RegisterUDInteger(mmdb::UDR_ATOM, "xy");
         int udd_mess_indexing_3 = asc.mol->RegisterUDInteger(mmdb::UDR_ATOM, "xz");

         std::cout << "read_pdb(): udd_mess_indexing: "
                   << udd_mess_indexing_1 << " " << udd_mess_indexing_2 << " " << udd_mess_indexing_3 << std::endl;
      }

      // 20221011-PE this constructor doesn't (at the moment) call make_bonds(). I
      // don't know if that is a good idea.
      int imol = molecules.size();
      coot::molecule_t m = coot::molecule_t(asc, imol, file_name);
      // m.make_bonds(&geom, &rot_prob_tables); // where does this go? Here or as a container function?
      molecules.push_back(m);
      status = imol;
   } else {
      std::cout << "debug:: in read_pdb() asc.read_success was " << asc.read_success
                << " for " << file_name << std::endl;
   }
   return status;
}

int
molecules_container_t::import_cif_dictionary(const std::string &cif_file_name, int imol_enc) {

   coot::read_refmac_mon_lib_info_t r = geom.init_refmac_mon_lib(cif_file_name, cif_dictionary_read_number, imol_enc);
   cif_dictionary_read_number++;

   std::cout << "debug:: import_cif_dictionary() cif_file_name(): success " << r.success << " with "
             << r.n_atoms << " atoms " << r.n_bonds << " bonds " << r.n_links << " links and momoner index "
             << r.monomer_idx << std::endl;

   return r.success;

}


int
molecules_container_t::get_monomer_from_dictionary(const std::string &comp_id,
                                                   int imol_enc,
                                                   bool idealised_flag) {

   int istat = -1; // unfound molecule

   // int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   mmdb::Manager *mol = geom.mol_from_dictionary(comp_id, imol_enc, idealised_flag);
   if (mol) {
      int imol = molecules.size();
      std::string name = comp_id;
      name += "_from_dict";
      // graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      // move_molecule_to_screen_centre_internal(imol);
      atom_selection_container_t asc = make_asc(mol);
      coot::molecule_t m = coot::molecule_t(asc, imol, name);
      molecules.push_back(m);
      istat = imol;
   } else {
      std::cout << "WARNING:: Null mol from mol_from_dictionary() with comp_id " << comp_id << " "
		<< idealised_flag << std::endl;
   }
   return istat;
}


int
molecules_container_t::get_monomer(const std::string &comp_id) {

   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   int imol = get_monomer_from_dictionary(comp_id, imol_enc, true); // idealized
   return imol;
}

int
molecules_container_t::get_monomer_and_position_at(const std::string &comp_id, int imol_enc, float x, float y, float z) {

   int imol = get_monomer_from_dictionary(comp_id, imol_enc, true);
   if (is_valid_model_molecule(imol)) {
      move_molecule_to_new_centre(imol, x, y, z);
   }
   return imol;
}


// 20221030-PE nice to have one day
// int
// molecules_container_t::get_monomer_molecule_by_network(const std::string &text) {

//    int imol = -1;
//    return imol;
// }


int
molecules_container_t::write_coordinates(int imol, const std::string &file_name) const {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].write_coordinates(file_name);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


int
molecules_container_t::read_mtz(const std::string &file_name,
                                const std::string &f, const std::string &phi, const std::string &weight,
                                bool use_weight, bool is_a_difference_map) {

   int imol = -1; // currently unset
   int imol_in_hope = molecules.size();

   std::string name_in = file_name + std::string(" ") + std::string(f) + std::string(" ") + std::string(phi);
   coot::molecule_t m(name_in, imol_in_hope);
   bool status = coot::util::map_fill_from_mtz(&m.xmap, file_name, f, phi, weight, use_weight, map_sampling_rate);
   if (is_a_difference_map)
      m.set_map_is_difference_map(true);
   if (status) {
      molecules.push_back(m);
      imol = imol_in_hope;
      std::cout << "DEBUG:: in read_mtz() " << file_name << " imol map: " << imol
                << " diff-map-status: " << is_a_difference_map << std::endl;
   }
   return imol;
}

#include "clipper-ccp4-map-file-wrapper.hh"
#include "coot-utils/slurp-map.hh"

int
molecules_container_t::read_ccp4_map(const std::string &file_name, bool is_a_difference_map) {

   int imol = -1; // currently unset
   int imol_in_hope = molecules.size();
   bool done = false;

   if (coot::util::is_basic_em_map_file(file_name)) {

      // fill xmap
      bool check_only = false;
      coot::molecule_t m(file_name, imol_in_hope);
      clipper::Xmap<float> &xmap = m.xmap;
      done = coot::util::slurp_fill_xmap_from_map_file(file_name, &xmap, check_only);
      try {
         // what's the point of this now?
         clipper_map_file_wrapper file;
         file.open_read(file_name);
         // set_is_em_map(file); // sets is_em_map_cached_flag
         // em = is_em_map_cached_flag;
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << file_name << std::endl;
         // bad_read = true;
      }

      if (done) {
         molecules.push_back(m);
         imol = imol_in_hope;
      }
   }

   if (! done) {
         std::cout << "INFO:: attempting to read CCP4 map: " << file_name << std::endl;
         // clipper::CCP4MAPfile file;
         clipper_map_file_wrapper w_file;
         try {
            w_file.open_read(file_name);

            // em = set_is_em_map(file);

            bool use_xmap = true; // not an nxmap
            if (true) {
               clipper::Grid_sampling fgs = w_file.grid_sampling();
               clipper::Cell fcell = w_file.cell();
               double vol = fcell.volume();
               if (vol < 1.0) {
                  std::cout << "WARNING:: read_ccp4_map(): non-sane unit cell volume " << vol << " - skip read"
                            << std::endl;
                  // bad_read = true;
               } else {
                  try {
                     clipper::CCP4MAPfile file;
                     file.open_read(file_name);
                     clipper::Xmap<float> xmap;
                     file.import_xmap(xmap);
                     if (xmap.is_null()) {
                        std::cout << "ERROR:: failed to read the map" << file_name << std::endl;
                     } else {
                        std::string name = file_name;
                        coot::molecule_t m(name, imol_in_hope);
                        m.xmap = xmap;
                        if (is_a_difference_map)
                           m.set_map_is_difference_map(true);
                        molecules.push_back(m); // oof.
                        imol = imol_in_hope;
                     }
                  }
                  catch (const clipper::Message_generic &exc) {
                     std::cout << "WARNING:: failed to read " << file_name
                               << " Bad ASU (inconsistant gridding?)." << std::endl;
                     // bad_read = true;
                  }
               }
            }
         } catch (const clipper::Message_base &exc) {
            std::cout << "WARNING:: failed to open " << file_name << std::endl;
            // bad_read = true;
         }
   }
   return imol;
}



coot::validation_information_t
molecules_container_t::density_fit_analysis(int imol_model, int imol_map) const {

   coot::validation_information_t r;
   r.name = "Density fit analysis";
#ifdef EMSCRIPTEN
   r.type = "DENSITY";
#else
   r.type = coot::DENSITY;
#endif
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
         // fill these
         mmdb::PResidue *SelResidues = 0;
         int nSelResidues = 0;

         auto atom_sel = molecules[imol_model].atom_sel;
         int selHnd = atom_sel.mol->NewSelection(); // yes, it's deleted.
         int imod = 1; // multiple models don't work on validation graphs

         atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                              "*", // chain_id
                              mmdb::ANY_RES, "*",
                              mmdb::ANY_RES, "*",
                              "*",  // residue name
                              "*",  // Residue must contain this atom name?
                              "*",  // Residue must contain this Element?
                              "*",  // altLocs
                              mmdb::SKEY_NEW // selection key
                              );
         atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

         for (int ir=0; ir<nSelResidues; ir++) {
            mmdb::Residue *residue_p = SelResidues[ir];
            coot::residue_spec_t res_spec(residue_p);
            mmdb::PAtom *residue_atoms=0;
            int n_residue_atoms;
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            double residue_density_score =
               coot::util::map_score(residue_atoms, n_residue_atoms, molecules[imol_map].xmap, 1);
            std::string l = res_spec.label();
            std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
            const std::string &chain_id = res_spec.chain_id;
            int this_resno = res_spec.res_no;
            coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
            coot::residue_validation_information_t rvi(res_spec, atom_spec, residue_density_score, l);
            r.add_residue_validation_information(rvi, chain_id);
         }
         atom_sel.mol->DeleteSelection(selHnd);
      }
   }
   r.set_min_max();
   return r;
}

//! density correlation validation information
coot::validation_information_t
molecules_container_t::density_correlation_analysis(int imol_model, int imol_map) const {

   coot::validation_information_t r;
   r.name = "Density correlation analysis";
#ifdef EMSCRIPTEN
   r.type = "CORRELATION";
#else
   r.type = coot::CORRELATION;
#endif
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

         mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;
         const clipper::Xmap<float> &xmap = molecules.at(imol_map).xmap;

         unsigned short int atom_mask_mode = 0;
         float atom_radius = 2.0;

         std::vector<coot::residue_spec_t> residue_specs;
         std::vector<mmdb::Residue *> residues = coot::util::residues_in_molecule(mol);
         for (unsigned int i=0; i<residues.size(); i++)
            residue_specs.push_back(coot::residue_spec_t(residues[i]));

         std::vector<std::pair<coot::residue_spec_t, float> > correlations =
            coot::util::map_to_model_correlation_per_residue(mol,
                                                             residue_specs,
                                                             atom_mask_mode,
                                                             atom_radius, // for masking
                                                             xmap);

         std::vector<std::pair<coot::residue_spec_t, float> >::const_iterator it;
         for (it=correlations.begin(); it!=correlations.end(); ++it) {
            const auto &r_spec(it->first);
            const auto &correl(it->second);

            std::string atom_name = " CA ";
            coot::atom_spec_t atom_spec(r_spec.chain_id, r_spec.res_no, r_spec.ins_code, atom_name, "");
            std::string label = "Correl: ";
            coot::residue_validation_information_t rvi(r_spec, atom_spec, correl, label);
            r.add_residue_validation_information(rvi, r_spec.chain_id);
         }

      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_model << std::endl;
   }
   r.set_min_max();
   return r;
}


//! rotamer validation information
coot::validation_information_t
molecules_container_t::rotamer_analysis(int imol_model) const {

   coot::validation_information_t r;
   r.name = "Rotamer analysis";
#ifdef EMSCRIPTEN
   r.type = "PROBABILITY";
#else
   r.type = coot::PROBABILITY;
#endif

   if (is_valid_model_molecule(imol_model)) {

      mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;

      // fill these
      mmdb::PResidue *SelResidues = 0;
      int nSelResidues = 0;

      int selHnd = mol->NewSelection(); // yes, it's deleted.
      int imod = 1; // multiple models don't work on validation graphs

      mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                           "*", // chain_id
                           mmdb::ANY_RES, "*",
                           mmdb::ANY_RES, "*",
                           "*",  // residue name
                           "*",  // Residue must contain this atom name?
                           "*",  // Residue must contain this Element?
                           "*",  // altLocs
                           mmdb::SKEY_NEW // selection key
                           );
      mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
         coot::residue_spec_t res_spec(residue_p);
         mmdb::PAtom *residue_atoms=0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

         // double residue_density_score = coot::util::map_score(residue_atoms, n_residue_atoms, xmap, 1);

         if (n_residue_atoms > 5) {

            std::string res_name = residue_p->GetResName();
            if (true) {

               coot::rotamer rot(residue_p);
               coot::rotamer_probability_info_t rpi = rot.probability_of_this_rotamer();
               double prob = rpi.probability;

               std::string l = res_spec.label();
               std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
               const std::string &chain_id = res_spec.chain_id;
               int this_resno = res_spec.res_no;
               coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
               coot::residue_validation_information_t rvi(res_spec, atom_spec, prob, l);
               r.add_residue_validation_information(rvi, chain_id);
            }
         }
      }
      mol->DeleteSelection(selHnd);
   }
   r.set_min_max();
   return r;
}

double
molecules_container_t::phi_psi_probability(const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) const {

      const clipper::Ramachandran *rama = &rc.rama;

      if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
      if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

      // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
      // if (phi_psi.is_pre_pro())
      // if (phi_psi.residue_name() != "GLY")
      // rama = &rc.rama_pre_pro;

      double rama_prob = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                           clipper::Util::d2rad(phi_psi.psi()));
      return rama_prob;
}

//! ramachandran validation information (formatted for a graph, not 3d)
coot::validation_information_t
molecules_container_t::ramachandran_analysis(int imol_model) const {

   coot::validation_information_t vi;
   vi.name = "Ramachandran plot Probability";
#ifdef EMSCRIPTEN
   vi.type = "PROBABILITY";
#else
   vi.type = coot::PROBABILITY;
#endif
   std::vector<coot::phi_psi_prob_t> rv = ramachandran_validation(imol_model);

   for (unsigned int i=0; i<rv.size(); i++) {
      std::string chain_id = rv[i].phi_psi.chain_id;
      coot::residue_spec_t residue_spec(rv[i].phi_psi.chain_id, rv[i].phi_psi.residue_number, rv[i].phi_psi.ins_code);
      double pr = rv[i].probability;
      std::string label = rv[i].phi_psi.chain_id + std::string(" ") + std::to_string(rv[i].phi_psi.residue_number);
      if (! rv[i].phi_psi.ins_code.empty())
         label += std::string(" ") + rv[i].phi_psi.ins_code;
      coot::atom_spec_t atom_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " CA ", "");
      coot::residue_validation_information_t rvi(residue_spec, atom_spec, pr, label);
      if (false)
         std::cout << "         " << residue_spec << " " << rv[i].phi_psi.phi() << " " << rv[i].phi_psi.psi()
                   << " pr " << pr << " " << std::endl;
      vi.add_residue_validation_information(rvi, chain_id);
   }
   vi.set_min_max();
   return vi;
}

//! peptide omega validation information
//! @returns a `validation_information_t`
coot::validation_information_t
molecules_container_t::peptide_omega_analysis(int imol) const {

   coot::validation_information_t vi;
   vi.name = "Peptide Omega Deviation";
#ifdef EMSCRIPTEN
   vi.type = "TORSION_ANGLE";
#else
   vi.type = coot::TORSION_ANGLE;
#endif

   if (is_valid_model_molecule(imol)) {

      bool mark_cis_peptides_as_bad_flag = false;
      bool m = mark_cis_peptides_as_bad_flag;
      std::vector<std::string> chain_ids = molecules[imol].chains_in_model();
      for (const auto &chain_id : chain_ids) {
         coot::chain_validation_information_t cvi(chain_id);
         coot::omega_distortion_info_container_t odi = molecules.at(imol).peptide_omega_analysis(geom, chain_id, m);
         for (const auto &od : odi.omega_distortions) {
            // oops - we have forgotten about the insertion code.
            coot::residue_spec_t res_spec(chain_id, od.resno, "");
            coot::atom_spec_t atom_spec(chain_id, od.resno, "", " CA ", "");
            std::string label = od.info_string;
            coot::residue_validation_information_t rvi(res_spec, atom_spec, od.distortion, label);
            cvi.add_residue_validation_information(rvi);
         }
         vi.cviv.push_back(cvi);
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return vi;
}


// #include "vertex.hh" // neeeded?

coot::simple_mesh_t
molecules_container_t::test_origin_cube() const {

   coot::simple_mesh_t mesh;

   std::vector<coot::api::vnc_vertex> vertices;
   std::vector<g_triangle> triangles;

   glm::vec4 c(0.5, 0.2, 0.5, 1.0); // colour

   // bottom
   coot::api::vnc_vertex v0(glm::vec3(0, 0, 0), glm::vec3(0,0,-1), c); vertices.push_back(v0);
   coot::api::vnc_vertex v1(glm::vec3(1, 0, 0), glm::vec3(0,0,-1), c); vertices.push_back(v1);
   coot::api::vnc_vertex v2(glm::vec3(0, 1, 0), glm::vec3(0,0,-1), c); vertices.push_back(v2);
   coot::api::vnc_vertex v3(glm::vec3(1, 1, 0), glm::vec3(0,0,-1), c); vertices.push_back(v3);

   // top
   coot::api::vnc_vertex v4(glm::vec3(0, 0, 1), glm::vec3(0,0,1), c); vertices.push_back(v4);
   coot::api::vnc_vertex v5(glm::vec3(1, 0, 1), glm::vec3(0,0,1), c); vertices.push_back(v5);
   coot::api::vnc_vertex v6(glm::vec3(0, 1, 1), glm::vec3(0,0,1), c); vertices.push_back(v6);
   coot::api::vnc_vertex v7(glm::vec3(1, 1, 1), glm::vec3(0,0,1), c); vertices.push_back(v7);

   // left
   coot::api::vnc_vertex v8 (glm::vec3(0, 0, 0), glm::vec3(-1,0,0), c); vertices.push_back(v8);
   coot::api::vnc_vertex v9 (glm::vec3(0, 1, 0), glm::vec3(-1,0,0), c); vertices.push_back(v9);
   coot::api::vnc_vertex v10(glm::vec3(0, 0, 1), glm::vec3(-1,0,0), c); vertices.push_back(v10);
   coot::api::vnc_vertex v11(glm::vec3(0, 1, 1), glm::vec3(-1,0,0), c); vertices.push_back(v11);

   // right
   coot::api::vnc_vertex v12(glm::vec3(1, 0, 0), glm::vec3(1,0,0), c); vertices.push_back(v12);
   coot::api::vnc_vertex v13(glm::vec3(1, 1, 0), glm::vec3(1,0,0), c); vertices.push_back(v13);
   coot::api::vnc_vertex v14(glm::vec3(1, 0, 1), glm::vec3(1,0,0), c); vertices.push_back(v14);
   coot::api::vnc_vertex v15(glm::vec3(1, 1, 1), glm::vec3(1,0,0), c); vertices.push_back(v15);

   // front
   coot::api::vnc_vertex v16(glm::vec3(0, 0, 0), glm::vec3(0,-1,0), c); vertices.push_back(v16);
   coot::api::vnc_vertex v17(glm::vec3(1, 0, 0), glm::vec3(0,-1,0), c); vertices.push_back(v17);
   coot::api::vnc_vertex v18(glm::vec3(0, 0, 1), glm::vec3(0,-1,0), c); vertices.push_back(v18);
   coot::api::vnc_vertex v19(glm::vec3(1, 0, 1), glm::vec3(0,-1,0), c); vertices.push_back(v19);

   // back
   coot::api::vnc_vertex v20(glm::vec3(0, 1, 0), glm::vec3(0,1,0), c); vertices.push_back(v20);
   coot::api::vnc_vertex v21(glm::vec3(1, 1, 0), glm::vec3(0,1,0), c); vertices.push_back(v21);
   coot::api::vnc_vertex v22(glm::vec3(0, 1, 1), glm::vec3(0,1,0), c); vertices.push_back(v22);
   coot::api::vnc_vertex v23(glm::vec3(1, 1, 1), glm::vec3(0,1,0), c); vertices.push_back(v23);

   triangles.push_back(g_triangle( 0, 1, 2));
   triangles.push_back(g_triangle( 1, 3, 2));
   triangles.push_back(g_triangle( 4, 5, 6));
   triangles.push_back(g_triangle( 5, 7, 6));
   triangles.push_back(g_triangle( 8, 9,10));
   triangles.push_back(g_triangle( 9,11,10));
   triangles.push_back(g_triangle(12,13,14));
   triangles.push_back(g_triangle(13,15,14));
   triangles.push_back(g_triangle(16,17,18));
   triangles.push_back(g_triangle(17,19,18));
   triangles.push_back(g_triangle(20,21,22));
   triangles.push_back(g_triangle(21,23,22));

   for (auto &vertex : vertices) {
      vertex.pos *= 10.0;
      vertex.pos += glm::vec3(-5.0, -5.0, -5.0);
   }

   coot::simple_mesh_t m(vertices, triangles);
   // m.translate(glm::vec3(-0.5, -0.5, -0.5));
   return m;
}

std::vector<coot::phi_psi_prob_t>
molecules_container_t::ramachandran_validation(int imol) const {

   std::vector<coot::phi_psi_prob_t> v;
   if (is_valid_model_molecule(imol))
      v = molecules[imol].ramachandran_validation(ramachandrans_container);
   return v;
}

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

coot::simple_mesh_t
molecules_container_t::get_ramachandran_validation_markup_mesh(int imol) const {

   // this function should be pushed into the coot::molecule_t class
   // (which means that the mesh will be copied)

   unsigned int num_subdivisions = 2;  // pass this
   float rama_ball_radius = 0.5;

   auto prob_raw_to_colour_rotation = [] (float prob) {
                                         if (prob > 0.5) prob = 0.5; // 0.4 and 2.5 f(for q) might be better (not tested)
                                         // good probabilities have q = 0
                                         // bad probabilities have q 0.66
                                         double q = (1.0 - 2.0 * prob);
                                         q = pow(q, 20.0);
                                         return q;
   };

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
   };

   auto phi_psi_probability = [] (const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) {

      const clipper::Ramachandran *rama = &rc.rama;

      if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
      if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

      // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
      // if (phi_psi.is_pre_pro())
      // if (phi_psi.residue_name() != "GLY")
      // rama = &rc.rama_pre_pro;

      double rama_prob = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                           clipper::Util::d2rad(phi_psi.psi()));
      return rama_prob;
   };

   auto test_ramachandran_probabilities = [] (const ramachandrans_container_t &rc) {

      std::vector<const clipper::Ramachandran *> ramas = { &rc.rama, &rc.rama_gly, &rc.rama_pro, &rc.rama_non_gly_pro };

      for (unsigned int ir=0; ir<ramas.size(); ir++) {
         for (unsigned int i=0; i<10; i++) {
            for (unsigned int j=0; j<10; j++) {
               double phi = static_cast<double>(i * 36.0) - 180.0;
               double psi = static_cast<double>(j * 36.0) - 180.0;
               double p = rc.rama.probability(phi, psi);
               std::cout << ir << "   "
                         << std::setw(10) << phi << " " << std::setw(10) << psi << " "
                         << std::setw(10) << p << std::endl;
            }
         }
      }
   };

   // test_ramachandran_probabilities(ramachandrans_container); // don't use rama_pre_pro without CLIPPER_HAS_TOP8000

   coot::simple_mesh_t mesh;

   // 20221126-PE Calm down the ultra-bright rama balls:
   float sober_factor = 0.75;

   if (is_valid_model_molecule(imol)) {

      std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaball = tessellate_octasphere(num_subdivisions);

      std::vector<coot::phi_psi_prob_t> ramachandran_goodness_spots = ramachandran_validation(imol);
      // now convert positions into meshes of balls
      int n_ramachandran_goodness_spots = ramachandran_goodness_spots.size();
      for (int i=0; i<n_ramachandran_goodness_spots; i++) {
         const coot::Cartesian &position = ramachandran_goodness_spots[i].position;
         // std::cout << "goodness spot " << i << " position " << position << std::endl;
         const coot::phi_psi_prob_t &phi_psi = ramachandran_goodness_spots[i];
         double prob_raw = phi_psi.probability;
         double q = prob_raw_to_colour_rotation(prob_raw);
         coot::colour_holder col = coot::colour_holder(q, 0.0, 1.0, false, std::string(""));
         glm::vec3 ball_position = cartesian_to_glm(position);
         unsigned int idx_base = mesh.vertices.size();
         unsigned int idx_tri_base = mesh.triangles.size();
         for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
            glm::vec4 col_v4(sober_factor * col.red, sober_factor * col.green, sober_factor * col.blue, 1.0f);
            const glm::vec3 &vertex_position = octaball.first[ibv];
            coot::api::vnc_vertex vertex(ball_position + rama_ball_radius * vertex_position, vertex_position, col_v4);
            mesh.vertices.push_back(vertex);
         }
         std::vector<g_triangle> octaball_triangles = octaball.second;
         mesh.triangles.insert(mesh.triangles.end(), octaball_triangles.begin(), octaball_triangles.end());

         for (unsigned int k=idx_tri_base; k<mesh.triangles.size(); k++)
            mesh.triangles[k].rebase(idx_base);
      }
   }
   return mesh;
}


mmdb::Atom *
molecules_container_t::get_atom(int imol, const coot::atom_spec_t &atom_spec) const {

   mmdb::Atom *r = nullptr;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_atom(atom_spec);
   }
   return r;
}

mmdb::Residue *
molecules_container_t::get_residue(int imol, const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *r = nullptr;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_residue(residue_spec);
   }
   return r;
}

// returns either the specified atom or null if not found
mmdb::Atom *
molecules_container_t::get_atom_using_cid(int imol, const std::string &cid) const {

   mmdb::Atom *at = nullptr;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::atom_spec_t> p = molecules[imol].cid_to_atom_spec(cid);
      if (p.first)
         at = molecules[imol].get_atom(p.second);
   }
   return at;
}

// returns either the specified residue or null if not found
mmdb::Residue *
molecules_container_t::get_residue_using_cid(int imol, const std::string &cid) const {
   mmdb::Residue *residue_p = nullptr;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::residue_spec_t> p = molecules[imol].cid_to_residue_spec(cid);
      if (p.first)
         residue_p = molecules[imol].get_residue(p.second);
   }
   return residue_p;
}


int
molecules_container_t::move_molecule_to_new_centre(int imol, float x, float y, float z) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::Cartesian new_centre(x,y,z);
      status = molecules[imol].move_molecule_to_new_centre(new_centre);
      set_updating_maps_need_an_update(imol); // weird thing to do usually
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

//! get the atom position - don't use this in emscript
std::pair<bool, coot::Cartesian>
molecules_container_t::get_atom_position(int imol, coot::atom_spec_t &atom_spec) {

   mmdb::Atom *at = get_atom(imol, atom_spec);
   if (at) {
      return std::pair<bool, coot::Cartesian> (true, coot::Cartesian(at->x, at->y, at->z));
   } else {
      return std::pair<bool, coot::Cartesian> (false, coot::Cartesian(0,0,0));
   }

}



coot::Cartesian
molecules_container_t::get_molecule_centre(int imol) const {

   coot::Cartesian c;
   if (is_valid_model_molecule(imol)) {
      c = molecules[imol].get_molecule_centre();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return c;
}


int
molecules_container_t::write_map(int imol, const std::string &file_name) const {

   int status= 0;
   if (is_valid_map_molecule(imol)) {
      status = molecules[imol].write_map(file_name);
   }
   return status;

}

// Mode is "COLOUR-BY-CHAIN-AND-DICTIONARY" or "CA+LIGANDS"
coot::simple_mesh_t
molecules_container_t::get_bonds_mesh(int imol, const std::string &mode,
                                      bool against_a_dark_background,
                                      float bonds_width, float atom_radius_to_bond_width_ratio,
                                      int smoothness_factor) {

   bool draw_hydrogen_atoms_flag = true; // pass this

   auto tp_0 = std::chrono::high_resolution_clock::now();

   coot::simple_mesh_t sm;
   if (is_valid_model_molecule(imol)) {
      sm = molecules[imol].get_bonds_mesh(mode, &geom, against_a_dark_background, bonds_width, atom_radius_to_bond_width_ratio,
                                          smoothness_factor, draw_hydrogen_atoms_flag, draw_missing_residue_loops_flag);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   if (show_timings) {
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "---------- timings: for get_bonds_mesh(): : " << d10 << " milliseconds " << std::endl;
   }

   return sm;
}


coot::simple_mesh_t
molecules_container_t::get_map_contours_mesh(int imol, double position_x, double position_y, double position_z,
                                             float radius, float contour_level) {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   coot::simple_mesh_t mesh;
   try {
      if (is_valid_map_molecule(imol)) {
         clipper::Coord_orth position(position_x, position_y, position_z);

         if (updating_maps_info.maps_need_an_update)
            update_updating_maps(updating_maps_info.imol_model);

         mesh = molecules[imol].get_map_contours_mesh(position, radius, contour_level);
      }
   }
   catch (...) {
      std::cout << "An error occured in " << __FUNCTION__<< "() - this should not happen " << std::endl;
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   if (show_timings) {
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "---------- timings: for get_map_contours_mesh(): : " << d10 << " milliseconds " << std::endl;
   }
   return mesh;
}


// get the rotamer dodecs for the model
coot::simple_mesh_t
molecules_container_t::get_rotamer_dodecs(int imol) {
   coot::simple_mesh_t m;
   if (is_valid_model_molecule(imol)) {
      return molecules[imol].get_rotamer_dodecs(&geom, &rot_prob_tables);
   } else {
      std::cout << "WARNING:: in " << __FUNCTION__ << "() imol " << imol << " was not a valid model molecule " << std::endl;
   }
   return m;
}

//! get the rotamer dodecs for the model, not const because it regenerates the bonds.
//! @return an `instanced_mesh_t`
coot::instanced_mesh_t
molecules_container_t::get_rotamer_dodecs_instanced(int imol) {

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      im = molecules[imol].get_rotamer_dodecs_instanced(&geom, &rot_prob_tables);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return im;
}



int
molecules_container_t::auto_fit_rotamer(int imol,
                                        const std::string &chain_id, int res_no, const std::string &ins_code,
                                        const std::string &alt_conf,
                                        int imol_map) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
         std::cout << "debug:: mc::auto_fit_rotamer() calling the coot_molecule version with "
                   << chain_id << " " << res_no << std::endl;
         status = molecules[imol].auto_fit_rotamer(chain_id, res_no, ins_code, alt_conf, xmap, geom);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "debug:: mc::auto_fit_rotamer() not a valid map index " << imol_map << std::endl;
      }
   } else {
      std::cout << "debug:: mc::auto_fit_rotamer() not a valid model molecule " << imol << std::endl;
   }
   return status;
}


int
molecules_container_t::delete_atom(int imol,
                                   const std::string &chain_id, int res_no, const std::string &ins_code,
                                   const std::string &atom_name, const std::string &alt_conf) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec(chain_id, res_no, ins_code, atom_name, alt_conf);
      status = molecules[imol].delete_atom(atom_spec);
      set_updating_maps_need_an_update(imol);
   }
   return status;
}

int
molecules_container_t::delete_atom_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      status = molecules[imol].delete_atom(atom_spec);
      set_updating_maps_need_an_update(imol);
   }
   return status;
}



int
molecules_container_t::delete_residue(int imol,
                                      const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
      status = molecules[imol].delete_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
   }
   return status;
}


int
molecules_container_t::delete_residue_using_cid(int imol, const std::string &residue_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_cid_to_residue_spec(imol, residue_cid);
      status = molecules[imol].delete_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
   }
   return status;
}

int
molecules_container_t::delete_residue_atoms_using_cid(int imol, const std::string &atom_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].delete_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
   }
   return status;
}

int
molecules_container_t::delete_residue_atoms_with_alt_conf(int imol, const std::string &chain_id,
                                                          int res_no, const std::string &ins_code,
                                                          const std::string &alt_conf) {
  int status = 0;
  if (is_valid_model_molecule(imol)) {
     std::string atom_cid = std::string("//") + chain_id + std::string("/") + std::to_string(res_no) + ins_code;
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].delete_residue(residue_spec);
      set_updating_maps_need_an_update(imol);
  }
  return status;
}



int
molecules_container_t::delete_chain_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].delete_chain_using_atom_cid(cid);
      set_updating_maps_need_an_update(imol);
   }
   return status;
}


//! delete the atoms specified in the CID selection
//! @return 1 on successful deletion, return 0 on failure to delete.
int
molecules_container_t::delete_literal_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].delete_literal_using_cid(cid);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}



//where scope in ["ATOM","WATER", "RESIDUE","CHAIN","MOLECULE"]
int
molecules_container_t::delete_using_cid(int imol, const std::string &cid, const std::string &scope) {

   int status = 0;
   if (scope == "ATOM")
      status = delete_atom_using_cid(imol, cid);
   if (scope == "RESIDUE")
      status = delete_residue_atoms_using_cid(imol, cid);
   if (scope == "CHAIN")
      status = delete_chain_using_cid(imol, cid);
   if (scope == "LITERAL")
      status = delete_literal_using_cid(imol, cid);
   if (scope == "MOLECULE")
      status = close_molecule(imol);
   return status;
}

// Old API
// int
// molecules_container_t::load_dictionary_file(const std::string &monomer_cif_file_name) {

//    int status = 0;

//    int read_number = 44;
//    geom.init_refmac_mon_lib(monomer_cif_file_name, read_number);
//    return status;
// }

std::vector<std::string>
molecules_container_t::non_standard_residue_types_in_model(int imol) const {
   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].non_standard_residue_types_in_model();
   }
   return v;
}

float
molecules_container_t::get_map_rmsd_approx(int imol) const {
   float rmsd = -1.1;
   if (is_valid_map_molecule(imol)) {
      rmsd = molecules[imol].get_map_rmsd_approx();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return rmsd;
}


std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::difference_map_peaks(int imol_map, int imol_protein, float n_rmsd) const {

   std::vector<coot::molecule_t::interesting_place_t> v;
   if (is_valid_model_molecule(imol_protein)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *m = get_mol(imol_protein);
         v = molecules[imol_map].difference_map_peaks(m, n_rmsd);
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_protein << std::endl;
   }
   return v;
}



// return a useful message if the addition did not work
std::pair<int, std::string>
molecules_container_t::add_terminal_residue_directly(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   std::string new_res_type = "ALA";
   int status = 0;
   std::string message;

   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map)) {
         clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
         coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
         std::pair<int, std::string> m = molecules[imol].add_terminal_residue_directly(residue_spec, new_res_type,
                                                                                       geom, xmap);
         status  = m.first;
         message = m.second;
         if (! message.empty())
            std::cout << "WARNING:: add_terminal_residue_directly(): " << message << std::endl;
         // write_coordinates(imol, "post-add-terminal-residue.pdb");
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_refinement_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return std::make_pair(status, message);
}

// 20221023-PE return an int for now so that I can write the binding
int
molecules_container_t::add_terminal_residue_directly_using_cid(int imol, const std::string &cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      if (! atom_spec.empty()) {
         auto p = add_terminal_residue_directly(imol, atom_spec.chain_id, atom_spec.res_no, atom_spec.ins_code);
         status = p.first;
      }
   }
   return status;
}


// reset the rail_points (calls reset_the_rail_points()), updates the maps (using internal/clipper SFC)
// so, update your contour lines meshes after calling this function.
int
molecules_container_t::connect_updating_maps(int imol_model, int imol_with_data_info_attached, int imol_map_2fofc, int imol_map_fofc) {

   int status = 0;

   rail_point_history.clear();
   updating_maps_info.imol_model = imol_model;
   updating_maps_info.imol_2fofc = imol_map_2fofc;
   updating_maps_info.imol_fofc  = imol_map_fofc;
   updating_maps_info.imol_with_data_info_attached = imol_with_data_info_attached;
   imol_difference_map = imol_map_fofc;

   // Let's force a sfcalc_genmap here.
   update_updating_maps(imol_model);

   return status;
}

void
molecules_container_t::associate_data_mtz_file_with_map(int imol, const std::string &data_mtz_file_name,
                                                        const std::string &f_col, const std::string &sigf_col,
                                                        const std::string &free_r_col) {

   if (is_valid_map_molecule(imol) || is_valid_model_molecule(imol)) {
      // 20221018-PE if free_r_col is not valid then Coot will (currently) crash on the structure factor calculation
      molecules[imol].associate_data_mtz_file_with_map(data_mtz_file_name, f_col, sigf_col, free_r_col);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid molecule " << imol << std::endl;
   }
}

/*! \brief Calculate structure factors from the model and update the given difference
           map accordingly */

// copied from:
// void
// graphics_info_t::sfcalc_genmap(int imol_model,
//                                int imol_map_with_data_attached,
//                                int imol_updating_difference_map) {
void
molecules_container_t::sfcalc_genmap(int imol_model,
                                     int imol_map_with_data_attached,
                                     int imol_updating_difference_map) {

   // I am keen for this function to be fast - so that it can be used with cryo-EM structures
   //
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_with_data_attached)) {
         if (true) {
            if (is_valid_map_molecule(imol_updating_difference_map)) {
               if (molecules[imol_updating_difference_map].is_difference_map_p()) {
                  clipper::Xmap<float> *xmap_p = &molecules[imol_updating_difference_map].xmap;
                  try {
                     if (! on_going_updating_map_lock) {
                        on_going_updating_map_lock = true;
                        molecules[imol_map_with_data_attached].fill_fobs_sigfobs();
                        const clipper::HKL_data<clipper::data32::F_sigF> *fobs_data =
                           molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
                        const clipper::HKL_data<clipper::data32::Flag> *free_flag =
                           molecules[imol_map_with_data_attached].get_original_rfree_flags();
                        if (fobs_data && free_flag) {
                           molecules[imol_model].sfcalc_genmap(*fobs_data, *free_flag, xmap_p);
                        } else {
                           std::cout << "sfcalc_genmap() either fobs_data or free_flag were not set " << std::endl;
                        }
                        on_going_updating_map_lock = false;
                     } else {
                        std::cout << "DEBUG:: on_going_updating_map_lock was set! - aborting map update." << std::endl;
                     }
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << rte.what() << std::endl;
                  }
               } else {
                  std::cout << "sfcalc_genmap() not a valid difference map " << imol_updating_difference_map << std::endl;
               }
            } else {
               std::cout << "sfcalc_genmap() not a valid map (diff) " << imol_updating_difference_map << std::endl;
            }
         }
      } else {
         std::cout << "sfcalc_genmap() not a valid map " << imol_map_with_data_attached << std::endl;
      }
   } else {
      std::cout << "sfcalc_genmap() not a valid model " << imol_model << std::endl;
   }
}


coot::util::sfcalc_genmap_stats_t
molecules_container_t::sfcalc_genmaps_using_bulk_solvent(int imol_model,
                                                         int imol_map_2fofc,  // this map should have the data attached.
                                                         int imol_map_fofc,
                                                         int imol_with_data_info_attached) {
   coot::util::sfcalc_genmap_stats_t stats;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_2fofc)) {
         if (is_valid_map_molecule(imol_map_fofc)) {
            if (molecules[imol_map_fofc].is_difference_map_p()) {
               try {
                  if (! on_going_updating_map_lock) {
                     on_going_updating_map_lock = true;
                     molecules[imol_with_data_info_attached].fill_fobs_sigfobs();

                     // 20210815-PE used to be const reference (get_original_fobs_sigfobs() function changed too)
                     // const clipper::HKL_data<clipper::data32::F_sigF> &fobs_data = molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
                     // const clipper::HKL_data<clipper::data32::Flag> &free_flag = molecules[imol_map_with_data_attached].get_original_rfree_flags();
                     // now the full object (40us for RNAse test).
                     // 20210815-PE OK, the const reference was not the problem. But we will leave it as it is now, for now.
                     //
                     clipper::HKL_data<clipper::data32::F_sigF> *fobs_data_p = molecules[imol_with_data_info_attached].get_original_fobs_sigfobs();
                     clipper::HKL_data<clipper::data32::Flag>   *free_flag_p = molecules[imol_with_data_info_attached].get_original_rfree_flags();

                     if (fobs_data_p && free_flag_p) {

                        if (true) { // sanity check data

                           const clipper::HKL_info &hkls_check = fobs_data_p->base_hkl_info();
                           const clipper::Spacegroup &spgr_check = hkls_check.spacegroup();
                           const clipper::Cell &cell_check = fobs_data_p->base_cell();
                           const clipper::HKL_sampling &sampling_check = fobs_data_p->hkl_sampling();

                           std::cout << "DEBUG:: in sfcalc_genmaps_using_bulk_solvent() imol_map_with_data_attached "
                                     << imol_map_2fofc << std::endl;

                           std::cout << "DEBUG:: Sanity check in graphics_info_t:sfcalc_genmaps_using_bulk_solvent(): HKL_info: "
                                     << "base_cell: " << cell_check.format() << " "
                                     << "spacegroup: " << spgr_check.symbol_xhm() << " "
                                     << "sampling is null: " << sampling_check.is_null() << " "
                                     << "resolution: " << hkls_check.resolution().limit() << " "
                                     << "invsqreslim: " << hkls_check.resolution().invresolsq_limit() << " "
                                     << "num_reflections: " << hkls_check.num_reflections()
                                     << std::endl;
                        }

                        clipper::Xmap<float> &xmap_2fofc = molecules[imol_map_2fofc].xmap;
                        clipper::Xmap<float> &xmap_fofc  = molecules[imol_map_fofc].xmap;
                        stats = molecules[imol_model].sfcalc_genmaps_using_bulk_solvent(*fobs_data_p, *free_flag_p, &xmap_2fofc, &xmap_fofc);

                     } else {
                        std::cout << "ERROR:: null data pointer in graphics_info_t::sfcalc_genmaps_using_bulk_solvent() " << std::endl;
                     }
                     on_going_updating_map_lock = false;
                  }
               }
               catch (const std::runtime_error &rte) {
                  std::cout << rte.what() << std::endl;
               }
            }
         }
      }
   }
   return stats;
}

int
molecules_container_t::rail_points_total() const { // the sum of all the rail ponts accumulated
   return rail_points_t::total(rail_point_history);
}

int
molecules_container_t::calculate_new_rail_points() {

   float rmsd = get_map_rmsd_approx(imol_difference_map);
   if (! rail_point_history.empty()) {
      const rail_points_t &prev = rail_point_history.back();
      rail_points_t new_points(rmsd, prev);
      rail_point_history.push_back(new_points);
      return new_points.map_rail_points_delta;
   } else {
      rail_points_t prev = rail_points_t(rmsd);
      rail_points_t new_points(rmsd, prev);
      rail_point_history.push_back(new_points);
      return new_points.map_rail_points_delta;
   }
}


// static
void
molecules_container_t::thread_for_refinement_loop_threaded() {

   // I think that there is a race condition here
   // check_and_warn_inverted_chirals_and_cis_peptides()
   // get called several times when the refine loop ends
   // (with success?).

   bool use_graphics_interface_flag = false;
   bool refinement_immediate_replacement_flag = true;

#if 0 // 20221018-PE this might not be the right thing

   if (restraints_lock) {
      if (false)
         std::cout << "debug:: thread_for_refinement_loop_threaded() restraints locked by "
                   << restraints_locking_function_name << std::endl;
      return;
   } else {

      if (use_graphics_interface_flag) {

         if (!refinement_immediate_replacement_flag) {

            // if there's not a refinement redraw function already running start up a new one.
            if (threaded_refinement_redraw_timeout_fn_id == -1) {
               GSourceFunc cb = GSourceFunc(regenerate_intermediate_atoms_bonds_timeout_function_and_draw);
	       // int id = gtk_timeout_add(15, cb, NULL);

               int timeout_ms = 15;
               timeout_ms = 30; // 20220503-PE try this value
	       int id = g_timeout_add(timeout_ms, cb, NULL);
               threaded_refinement_redraw_timeout_fn_id = id;
            }
         }
      }

      continue_threaded_refinement_loop = true;
      std::thread r(refinement_loop_threaded);
      r.detach();
   }
#endif

}


int
molecules_container_t::refine_direct(int imol, std::vector<mmdb::Residue *> rv, const std::string &alt_loc) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map)) {
         const clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
         molecules[imol].refine_direct(rv, alt_loc, xmap, map_weight, geom, refinement_is_quiet);
      }
   }
   return status;
}

int
molecules_container_t::refine_residues_using_atom_cid(int imol, const std::string &cid, const std::string &mode) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map)) {
         coot::atom_spec_t spec = atom_cid_to_atom_spec(imol, cid);
         status = refine_residues(imol, spec.chain_id, spec.res_no, spec.ins_code, spec.alt_conf, mode);
      } else {
         std::cout << "Not a valid map molecule " << imol_refinement_map << std::endl;
      }
   } else {
      std::cout << "Not a valid model molecule " << imol << std::endl;
   }
   return status;
}



int
molecules_container_t::refine_residues(int imol, const std::string &chain_id, int res_no, const std::string &ins_code,
                                       const std::string &alt_conf, const std::string &mode) {
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec(chain_id, res_no, ins_code);
      std::vector<mmdb::Residue *> rv = molecules[imol].select_residues(residue_spec, mode);
      if (! rv.empty()) {
         status = refine_direct(imol, rv, alt_conf);
         set_updating_maps_need_an_update(imol);
      } else {
         std::cout << "WARNING:: in refine_residues() - empty residues." << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::refine_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      std::vector<mmdb::Residue *> rv = molecules[imol].select_residues(chain_id, res_no_start, res_no_end);
      if (! rv.empty()) {
         std::string alt_conf = "";
         status = refine_direct(imol, rv, alt_conf);
      } else {
         std::cout << "WARNING:: in refine_residues() - empty residues." << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}



coot::refinement_results_t
molecules_container_t::refine_residues_vec(int imol,
                                           const std::vector<mmdb::Residue *> &residues,
                                           const std::string &alt_conf,
                                           mmdb::Manager *mol) {
   bool use_map_flag = true;
   if (false)
      std::cout << "INFO:: refine_residues_vec() with altconf \""
		<< alt_conf << "\"" << std::endl;

   coot::refinement_results_t rr = generate_molecule_and_refine(imol, residues, alt_conf, mol, use_map_flag);
   return rr;
}

// return -1 on failure to find a residue for insertion index
//
int
molecules_container_t::find_serial_number_for_insert(int seqnum_new,
                                                     const std::string &ins_code_for_new,
                                                     mmdb::Chain *chain_p) const {

   int iserial_no = -1;
   if (chain_p) {
      int current_diff = 999999;
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) { // ires is a serial number
	 mmdb::Residue *residue = chain_p->GetResidue(ires);

	 // we are looking for the smallest negative diff:
	 //
	 int diff = residue->GetSeqNum() - seqnum_new;
	 if ( (diff > 0) && (diff < current_diff) ) {
	    iserial_no = ires;
	    current_diff = diff;
	 } else {
	    if (diff == 0) {
	       std::string ins_code_this = residue->GetInsCode();
	       if (ins_code_this > ins_code_for_new) {
		  iserial_no = ires;
		  break;
	       }
	    }
	 }
      }
   }
   return iserial_no;
}


#include "coords/mmdb-extras.h"

std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> >
molecules_container_t::create_mmdbmanager_from_res_vector(const std::vector<mmdb::Residue *> &residues,
                                                          int imol,
                                                          mmdb::Manager *mol_in,
                                                          std::string alt_conf) {

   // returned entities
   mmdb::Manager *new_mol = 0;
   std::vector<mmdb::Residue *> rv; // gets checked

   float dist_crit = 5.0;
   bool debug = false;

   if (debug) {
      std::cout << "############ starting create_mmdbmanager_from_res_vector() with these "
		<< " residues " << std::endl;
      for (std::size_t ii=0; ii<residues.size(); ii++)
	 std::cout << "   " << coot::residue_spec_t(residues[ii])  << std::endl;
      int udd_atom_index_handle = mol_in->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
      std::cout << "############ udd for atom index from seeding molecule " << udd_atom_index_handle
		<< std::endl;
      for (std::size_t ii=0; ii<residues.size(); ii++) {
	 mmdb::Residue *residue_p = residues[ii];
	 mmdb::Atom **residue_atoms = 0;
	 int n_residue_atoms;
	 residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    mmdb::Atom *at = residue_atoms[iat];
	    int idx = -1;
	    at->GetUDData(udd_atom_index_handle, idx);
	    std::cout << "#### input residue atom " << coot::atom_spec_t(at) << " had udd index "
		      << idx << std::endl;
	 }
      }
   }

   int n_flanker = 0; // a info/debugging counter

   if (residues.size() > 0) {

      // Also add the index of the reference residue (the one in molecules[imol].atom_selection.mol)
      // to the molecule that we are construction here. So that we can properly link
      // the residues in restraints_container (there we rather need to know the references indices,
      // not the indices from the fragment molecule)
      //

      std::pair<bool,std::string> use_alt_conf(false, "");
      if (! alt_conf.empty())
	 use_alt_conf = std::pair<bool, std::string> (true, alt_conf);

      std::cout << "----------------- in create_mmdbmanager_from_res_vector() alt_conf is "
                << "\"" << alt_conf << "\"" << std::endl;
      std::cout << "----------------- in create_mmdbmanager_from_res_vector() use_alt_conf is "
                << use_alt_conf.first << "\"" << use_alt_conf.second << "\"" << std::endl;

      std::pair<bool, mmdb::Manager *> n_mol_1 =
	 coot::util::create_mmdbmanager_from_residue_vector(residues, mol_in, use_alt_conf);

      // check that first is sane, so indent all this lot (when it works)

      if (n_mol_1.first) {

	 int index_from_reference_residue_handle =
	    n_mol_1.second->GetUDDHandle(mmdb::UDR_RESIDUE, "index from reference residue");

	 if (false) { // debug
	    int imod = 1;
	    mmdb::Model *model_p = n_mol_1.second->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     int n_atoms = residue_p->GetNumberOfAtoms();
		     for (int iat=0; iat<n_atoms; iat++) {
			mmdb::Atom *at = residue_p->GetAtom(iat);
			int idx = -1;
			at->GetUDData(index_from_reference_residue_handle, idx);
			std::cout << "   create_mmdbmanager_from_residue_vector() returns this mol atom "
				  << iat << " " << coot::atom_spec_t(at) << " with idx " << idx << std::endl;
		     }
		  }
	       }
	    }
	 }

	 new_mol = n_mol_1.second;
	 mmdb::Model *model_p = new_mol->GetModel(1);

	 // how many (free) residues were added to that model? (add them to rv)
	 //
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       rv.push_back(residue_p);
	    }
	 }

	 if (false) {
	    for (std::size_t ir=0; ir<rv.size(); ir++) {
	       mmdb::Residue *r = rv[ir];
	       std::cout << "Moving Residue " << coot::residue_spec_t(r) << std::endl;
	       mmdb::Atom **residue_atoms = 0;
	       int n_residue_atoms;
	       r->GetAtomTable(residue_atoms, n_residue_atoms);
	       for (int iat=0; iat<n_residue_atoms; iat++) {
		  mmdb::Atom *at = residue_atoms[iat];
		  std::cout << "    " << coot::atom_spec_t(at) << std::endl;
	       }
	    }
	 }

	 short int whole_res_flag = 0;
	 int atom_index_udd_handle = molecules[imol].atom_sel.UDDAtomIndexHandle;

	 // Now the flanking residues:
	 //
	 std::vector<mmdb::Residue *> flankers_in_reference_mol;
	 flankers_in_reference_mol.reserve(32); // say

	 // find the residues that are close to the residues of
	 // residues that are not part of residues
	 //
	 // We don't have quite the function that we need in coot-utils,
	 // so we need to munge residues in to local_residues:
	 std::vector<std::pair<bool, mmdb::Residue *> > local_residues;
	 local_residues.resize(residues.size());
	 for (std::size_t ires=0; ires<residues.size(); ires++)
	    local_residues[ires] = std::pair<bool, mmdb::Residue *>(false, residues[ires]);
	 std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnr =
	    coot::residues_near_residues(local_residues, mol_in, dist_crit);
	 // now fill @var{flankers_in_reference_mol} from rnr, avoiding residues
	 // already in @var{residues}.
	 std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
	 for (it=rnr.begin(); it!=rnr.end(); ++it) {
	    const std::set<mmdb::Residue *> &s = it->second;
	    std::set<mmdb::Residue *>::const_iterator its;
	    for (its=s.begin(); its!=s.end(); ++its) {
	       mmdb::Residue *tres = *its;
	       if (std::find(residues.begin(), residues.end(), tres) == residues.end())
		  if (std::find(flankers_in_reference_mol.begin(), flankers_in_reference_mol.end(), tres) == flankers_in_reference_mol.end())
		     flankers_in_reference_mol.push_back(tres);
	    }
	 }

	 // So we have a vector of residues that were flankers in the
	 // reference molecule, we need to add copies of those to
	 // new_mol (making sure that they go into the correct chain).
	 //
	 if (false) { // debug
	    std::cout << "debug:: ############ Found " << flankers_in_reference_mol.size()
		      << " flanking residues" << std::endl;

	    for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++)
	       std::cout << "     #### flankers_in_reference_mol: " << ires << " "
			 << coot::residue_spec_t(flankers_in_reference_mol[ires]) << std::endl;
	 }


	 for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++) {
	    mmdb::Residue *r;

	    std::string ref_res_chain_id = flankers_in_reference_mol[ires]->GetChainID();

	    mmdb::Chain *chain_p = NULL;
	    int n_new_mol_chains = model_p->GetNumberOfChains();
	    for (int ich=0; ich<n_new_mol_chains; ich++) {
	       if (ref_res_chain_id == model_p->GetChain(ich)->GetChainID()) {
		  chain_p = model_p->GetChain(ich);
		  break;
	       }
	    }

	    if (! chain_p) {
	       // Add a new one then.
	       chain_p = new mmdb::Chain;
	       chain_p->SetChainID(ref_res_chain_id.c_str());
	       model_p->AddChain(chain_p);
	    }

	    if (false)
	       std::cout << "debug:: flankers_in_reference_mol " << ires << " "
			 << coot::residue_spec_t(flankers_in_reference_mol[ires]) << " "
			 << "had index " << flankers_in_reference_mol[ires]->index
			 << std::endl;

            // get rid of this function at some stage
            bool embed_in_chain = false;
	    r = coot::deep_copy_this_residue_old_style(flankers_in_reference_mol[ires],
					               alt_conf, whole_res_flag,
					               atom_index_udd_handle, embed_in_chain);

	    if (r) {

	       r->PutUDData(index_from_reference_residue_handle, flankers_in_reference_mol[ires]->index);

	       // copy over the atom indices. UDDAtomIndexHandle in mol_n becomes UDDOldAtomIndexHandle
	       // indices in the returned molecule

	       int sni = find_serial_number_for_insert(r->GetSeqNum(), r->GetInsCode(), chain_p);

	       if (false) { // debug
		  mmdb::Atom **residue_atoms = 0;
		  int n_residue_atoms;
		  std::cout << "Flanker Residue " << coot::residue_spec_t(r) << std::endl;
		  r->GetAtomTable(residue_atoms, n_residue_atoms);
		  for (int iat=0; iat<n_residue_atoms; iat++) {
		     mmdb::Atom *at = residue_atoms[iat];
		     std::cout << "    " << coot::atom_spec_t(at) << std::endl;
		  }
	       }

	       if (sni == -1)
		  chain_p->AddResidue(r); // at the end
	       else
		  chain_p->InsResidue(r, sni);
	       r->seqNum = flankers_in_reference_mol[ires]->GetSeqNum();
	       r->SetResName(flankers_in_reference_mol[ires]->GetResName());
	       n_flanker++;

	       if (false)
		  std::cout << "debug:: create_mmdbmanager_from_residue_vector() inserted/added flanker "
			    << coot::residue_spec_t(r) << std::endl;

	    }
	 }

	 // super-critical for correct peptide bonding in refinement!
	 //
	 coot::util::pdbcleanup_serial_residue_numbers(new_mol);

	 if (debug) {
	    int imod = 1;
	    mmdb::Model *model_p = new_mol->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     std::cout << "create_mmdb..  ^^^ " << coot::residue_spec_t(residue_p) << " "
			       << residue_p << " index " << residue_p->index
			       << std::endl;
		  }
	       }
	    }
	 }

	 if (debug)
	    std::cout << "DEBUG:: in create_mmdbmanager_from_res_vector: " << rv.size()
		      << " free residues and " << n_flanker << " flankers" << std::endl;
      }
   }

   return std::pair <mmdb::Manager *, std::vector<mmdb::Residue *> > (new_mol, rv);
}



std::string
molecules_container_t::adjust_refinement_residue_name(const std::string &resname) const {

   std::string r = resname;
   if (resname == "UNK") r = "ALA"; // hack for KC/buccaneer.
   if (resname.length() > 2)
      if (resname[2] == ' ')
	 r = resname.substr(0,2);
   return r;
}


// Return 0 (first) if any of the residues don't have a dictionary
// entry and a list of the residue type that don't have restraints.
//
std::pair<int, std::vector<std::string> >
molecules_container_t::check_dictionary_for_residue_restraints(int imol, mmdb::PResidue *SelResidues, int nSelResidues) {

   int status;
   bool status_OK = 1; // pass, by default
   std::vector<std::string> res_name_vec;

   for (int ires=0; ires<nSelResidues; ires++) {
      std::string resn(SelResidues[ires]->GetResName());
      std::string resname = adjust_refinement_residue_name(resn);
      status = geom.have_dictionary_for_residue_type(resname, imol, cif_dictionary_read_number);
      cif_dictionary_read_number++;
      if (! status) {
	 status_OK = 0;
	 res_name_vec.push_back(resname);
      }

      if (0)
	 std::cout << "DEBUG:: have_dictionary_for_residues() on residue "
		   << ires << " of " << nSelResidues << ", "
		   << resname << " returns "
		   << status << std::endl;
      cif_dictionary_read_number++;
   }
   return std::pair<int, std::vector<std::string> > (status_OK, res_name_vec);
}

std::pair<int, std::vector<std::string> >
molecules_container_t::check_dictionary_for_residue_restraints(int imol, const std::vector<mmdb::Residue *> &residues) {

   std::vector<std::string> res_name_vec;
   std::pair<int, std::vector<std::string> > r(0, res_name_vec);
   for (unsigned int i=0; i<residues.size(); i++) {
      std::string resname = adjust_refinement_residue_name(residues[i]->GetResName());
      int status = geom.have_dictionary_for_residue_type(resname, imol, cif_dictionary_read_number);
      if (! status) {
	 r.first = 0;
	 r.second.push_back(resname);
      }
      cif_dictionary_read_number++; // not sure why this is needed.
   }
   return r;
}


std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > >
molecules_container_t::make_rotamer_torsions(const std::vector<std::pair<bool, mmdb::Residue *> > &local_residues) const {

   std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > v;
   for (unsigned int i=0; i<local_residues.size(); i++) {
      if (! local_residues[i].first) {
         mmdb::Residue *residue_p = local_residues[i].second;
         std::string rn(residue_p->GetResName());
         if (coot::util::is_standard_amino_acid_name(rn)) {
            std::string alt_conf; // run through them all, ideally.
            coot::rotamer rot(residue_p, alt_conf, 1);
            coot::closest_rotamer_info_t cri = rot.get_closest_rotamer(rn);
            if (cri.residue_chi_angles.size() > 0) {
               std::vector<coot::dict_torsion_restraint_t> dictionary_vec;
               std::vector<std::vector<std::string> > rotamer_atom_names = rot.rotamer_atoms(rn);

               if (cri.residue_chi_angles.size() != rotamer_atom_names.size()) {

                  std::cout << "-------------- mismatch for " << coot::residue_spec_t(residue_p) << " "
                            << cri.residue_chi_angles.size() << " "  << rotamer_atom_names.size()
                            << " ---------------" << std::endl;
               } else {

                  for (unsigned int ichi=0; ichi<cri.residue_chi_angles.size(); ichi++) {
                     // we have to convert chi angles to atom names
                     double esd = 3.0; // 20210315-PE was 10.0. I want them tighter than that.
                     int per = 1;
                     std::string id = "chi " + coot::util::int_to_string(cri.residue_chi_angles[ichi].first);
                     const std::string &atom_name_1 = rotamer_atom_names[ichi][0];
                     const std::string &atom_name_2 = rotamer_atom_names[ichi][1];
                     const std::string &atom_name_3 = rotamer_atom_names[ichi][2];
                     const std::string &atom_name_4 = rotamer_atom_names[ichi][3];
                     double torsion = cri.residue_chi_angles[ichi].second;
                     coot::dict_torsion_restraint_t dr(id, atom_name_1, atom_name_2, atom_name_3, atom_name_4,
                                                       torsion, esd, per);
                     dictionary_vec.push_back(dr);
                  }

                  if (dictionary_vec.size() > 0) {
                     std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > p(residue_p, dictionary_vec);
                     v.push_back(p);
                  }
               }
            }
         }
      }
   }
   return v;
}



atom_selection_container_t
molecules_container_t::make_moving_atoms_asc(mmdb::Manager *residues_mol,
                                             const std::vector<mmdb::Residue *> &residues) const {

   // This also rebonds the imol_moving_atoms molecule

   atom_selection_container_t local_moving_atoms_asc;
   local_moving_atoms_asc.UDDAtomIndexHandle = -1;
   local_moving_atoms_asc.UDDOldAtomIndexHandle = residues_mol->GetUDDHandle(mmdb::UDR_ATOM, "old atom index");

   int SelHnd = residues_mol->NewSelection();

   for (unsigned int ir=0; ir<residues.size(); ir++) {
      const char *chain_id = residues[ir]->GetChainID();
      const char *inscode = residues[ir]->GetInsCode();
      int resno = residues[ir]->GetSeqNum();
      residues_mol->Select(SelHnd, mmdb::STYPE_ATOM,
                           0, chain_id,
                           resno, // starting resno, an int
                           inscode, // any insertion code
                           resno, // ending resno
                           inscode, // ending insertion code
                           "*", // any residue name
                           "*", // atom name
                           "*", // elements
                           "*",  // alt loc.
                           mmdb::SKEY_OR);
   }

   local_moving_atoms_asc.mol = residues_mol;
   local_moving_atoms_asc.SelectionHandle = SelHnd;
   residues_mol->GetSelIndex(local_moving_atoms_asc.SelectionHandle,
                             local_moving_atoms_asc.atom_selection,
                             local_moving_atoms_asc.n_selected_atoms);


   if (true) {
      std::cout << "returning an atom selection for all moving atoms "
                << local_moving_atoms_asc.n_selected_atoms << " atoms "
                << std::endl;
   }

   // This new block added so that we don't draw atoms in the "static" molecule when we have the
   // corresponding atoms in the moving atoms.
   //
#if 0 // 20221018-PE there is no drawing at the momment
   const atom_selection_container_t &imol_asc = molecules[imol_moving_atoms].atom_sel;
   std::set<int> atom_set = coot::atom_indices_in_other_molecule(imol_asc, local_moving_atoms_asc);

   if (false) { // debug atoms in other molecule
      std::set<int>::const_iterator it;
      for(it=atom_set.begin(); it!=atom_set.end(); it++) {
         int idx = *it;
         mmdb::Atom *at = imol_asc.atom_selection[idx];
         coot::atom_spec_t as(at);
         std::cout << " this is a moving atom: " << idx << " " << as << std::endl;
      }
   }

   if (false) { // debug old atom index
      for (int i=0; i<local_moving_atoms_asc.n_selected_atoms; i++) {
         mmdb::Atom *at = local_moving_atoms_asc.atom_selection[i];
         coot::atom_spec_t as(at);
         int idx = -1;
         at->GetUDData(local_moving_atoms_asc.UDDOldAtomIndexHandle, idx);
         std::cout << "DEBUG:: in make_moving_atoms_asc " << as << " idx " << idx << std::endl;
      }
   }
   // now rebond molecule imol without bonds to atoms in atom_set
   if (atom_set.size() > 0) {
      if (regenerate_bonds_needs_make_bonds_type_checked_flag) {
         molecules[imol_moving_atoms].make_bonds_type_checked(atom_set, __FUNCTION__);
      }
   }
#endif

   return local_moving_atoms_asc;
}

// static
void
molecules_container_t::all_atom_pulls_off() {
   for (std::size_t i=0; i<atom_pulls.size(); i++)
      atom_pulls[i].off();
   atom_pulls.clear();
}


// return the state of having found restraints.
bool
molecules_container_t::make_last_restraints(const std::vector<std::pair<bool,mmdb::Residue *> > &local_residues,
                                      const std::vector<mmdb::Link> &links,
                                      const coot::protein_geometry &geom,
                                      mmdb::Manager *mol_for_residue_selection,
                                      const std::vector<coot::atom_spec_t> &fixed_atom_specs,
                                      coot::restraint_usage_Flags flags,
                                      bool use_map_flag,
                                      const clipper::Xmap<float> *xmap_p) {

   bool do_torsion_restraints = true; // make this a data member
   double torsion_restraints_weight = 10.0;
   bool convert_dictionary_planes_to_improper_dihedrals_flag = false;
   double geometry_vs_map_weight = 25.5;
   bool do_trans_peptide_restraints = true;
   double rama_plot_restraints_weight = 20.0;
   bool do_rama_restraints = false;
   bool make_auto_h_bond_restraints_flag = false;
   coot::pseudo_restraint_bond_type pseudo_bonds_type = coot::NO_PSEUDO_BONDS;
   bool use_harmonic_approximation_for_NBCs = false;
   double pull_restraint_neighbour_displacement_max_radius = 1.0;
   double lennard_jones_epsilon = 1.0;
   int restraints_rama_type = 1;
   bool do_rotamer_restraints = false;
   double geman_mcclure_alpha = 0.1;
   bool do_numerical_gradients =  false;
   bool draw_gl_ramachandran_plot_flag = false;
   bool use_graphics_interface_flag = false;


   if (last_restraints) {
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "    ERROR:: A: last_restraints not cleared up " << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
   }

   if (false) { // these are the passed residues, nothing more.
      std::cout << "debug:: on construction of restraints_container_t local_residues: "
                << std::endl;
      for (std::size_t jj=0; jj<local_residues.size(); jj++) {
         std::cout << "   " << coot::residue_spec_t(local_residues[jj].second)
                   << " is fixed: " << local_residues[jj].first << std::endl;
      }
   }

   // moving_atoms_extra_restraints_representation.clear();
   continue_threaded_refinement_loop = true; // no longer set in refinement_loop_threaded()

   // the refinment of torsion seems a bit confused? If it's in flags, why does it need an flag
   // of its own? I suspect that it doesn't. For now I will keep it (as it was).
   //
   bool do_residue_internal_torsions = false;
   if (do_torsion_restraints) {
      do_residue_internal_torsions = 1;
   }

   last_restraints = new
      coot::restraints_container_t(local_residues,
                                   links,
                                   geom,
                                   mol_for_residue_selection,
                                   fixed_atom_specs, xmap_p);

   std::cout << "debug:: on creation last_restraints is " << last_restraints << std::endl;

   last_restraints->set_torsion_restraints_weight(torsion_restraints_weight);

   if (convert_dictionary_planes_to_improper_dihedrals_flag) {
      last_restraints->set_convert_plane_restraints_to_improper_dihedral_restraints(true);
   }

   // This seems not to work yet.
   // last_restraints->set_dist_crit_for_bonded_pairs(9.0);

   if (use_map_flag)
      last_restraints->add_map(geometry_vs_map_weight);

   unsigned int n_threads = coot::get_max_number_of_threads();
   if (n_threads > 0)
      last_restraints->thread_pool(&static_thread_pool, n_threads);

   all_atom_pulls_off();
   particles_have_been_shown_already_for_this_round_flag = false;

   // elsewhere do this:
   // gtk_widget_remove_tick_callback(glareas[0], wait_for_hooray_refinement_tick_id);

   // moving_atoms_visited_residues.clear(); // this is used for HUD label colour

   int n_restraints = last_restraints->make_restraints(imol_moving_atoms,
                                                       geom, flags,
                                                       do_residue_internal_torsions,
                                                       do_trans_peptide_restraints,
                                                       rama_plot_restraints_weight,
                                                       do_rama_restraints,
                                                       true, true, make_auto_h_bond_restraints_flag,
                                                       pseudo_bonds_type);
                                                       // link and flank args default true

   if (use_harmonic_approximation_for_NBCs) {
      std::cout << "INFO:: using soft harmonic restraints for NBC" << std::endl;
      last_restraints->set_use_harmonic_approximations_for_nbcs(true);
   }

   if (pull_restraint_neighbour_displacement_max_radius > 1.99) {
      last_restraints->set_use_proportional_editing(true);
      last_restraints->pull_restraint_neighbour_displacement_max_radius =
         pull_restraint_neighbour_displacement_max_radius;
   }

   last_restraints->set_geman_mcclure_alpha(geman_mcclure_alpha);
   last_restraints->set_lennard_jones_epsilon(lennard_jones_epsilon);
   last_restraints->set_rama_type(restraints_rama_type);
   last_restraints->set_rama_plot_weight(rama_plot_restraints_weight); // >2? danger of non-convergence
                                                                       // if planar peptide restraints are used
   // Oh, I see... it's not just the non-Bonded contacts of the hydrogens.
   // It's the planes, chiral and angles too. Possibly bonds too.
   // How about marking non-H atoms in restraints that contain H atoms as
   // "invisible"? i.e. non-H atoms are not influenced by the positions of the
   // Hydrogen atoms (but Hydrogen atoms *are* influenced by the positions of the
   // non-Hydrogen atoms). This seems like a lot of work. Might be easier to turn
   // off angle restraints for H-X-X (but not H-X-H) in the first instance, that
   // should go most of the way to what "invisible" atoms would do, I imagine.
   // is_H_non_bonded_contact should be renamed to is_H_turn_offable_restraint
   // or something.
   //
   // last_restraints->set_apply_H_non_bonded_contacts(false);

   if (do_rotamer_restraints) {
      std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > rotamer_torsions = make_rotamer_torsions(local_residues);
      std::cout << "debug:: calling add_or_replace_torsion_restraints_with_closest_rotamer_restraints() from make_last_restraints() " << std::endl;
      last_restraints->add_or_replace_torsion_restraints_with_closest_rotamer_restraints(rotamer_torsions);
   }

   if (molecules[imol_moving_atoms].extra_restraints.has_restraints()) {
      std::cout << "debug:: calling add_extra_restraints() from make_last_restraints() " << std::endl;
      last_restraints->add_extra_restraints(imol_moving_atoms, "user-defined from make_last_restraints()",
                                            molecules[imol_moving_atoms].extra_restraints, geom);
   }

   if (do_numerical_gradients)
      last_restraints->set_do_numerical_gradients();

   bool found_restraints_flag = false;

   if (last_restraints->size() > 0) {

      last_restraints->analyze_for_bad_restraints();
      thread_for_refinement_loop_threaded();
      found_restraints_flag = true;
      // rr.found_restraints_flag = true;
      draw_gl_ramachandran_plot_flag = true;

      // are you looking for conditionally_wait_for_refinement_to_finish() ?

      if (refinement_immediate_replacement_flag) {
         // wait until refinement finishes
         while (restraints_lock) {
            std::this_thread::sleep_for(std::chrono::milliseconds(7));
            std::cout << "INFO:: make_last_restraints() [immediate] restraints locked by "
                      << restraints_locking_function_name << std::endl;
         }
      }

   } else {
      continue_threaded_refinement_loop = false;
      if (use_graphics_interface_flag) {
         // GtkWidget *widget = create_no_restraints_info_dialog();
         // GtkWidget *widget = widget_from_builder("no_restraints_info_dialog");
         // gtk_widget_show(widget);
      }
   }

   return found_restraints_flag;
}


// simple mmdb::Residue * interface to refinement.  20081216
coot::refinement_results_t
molecules_container_t::generate_molecule_and_refine(int imol,  // needed for UDD Atom handle transfer
                                                    const std::vector<mmdb::Residue *> &residues_in,
                                                    const std::string &alt_conf,
                                                    mmdb::Manager *mol,
                                                    bool use_map_flag) {

   // 20221018-PE make a function in the class
   auto set_refinement_flags = [] () {
      return coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
   };
   int cif_dictionary_read_number = 44; // make this a class member

   bool do_torsion_restraints = true;
   bool do_rama_restraints = false; // or true?
   bool moving_atoms_have_hydrogens_displayed = false;


   coot::refinement_results_t rr(0, GSL_CONTINUE, "");

   if (is_valid_map_molecule(imol_refinement_map) || (! use_map_flag)) {
      // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      coot::restraint_usage_Flags flags = set_refinement_flags();
      bool do_residue_internal_torsions = false;
      if (do_torsion_restraints) {
         do_residue_internal_torsions = 1;
         flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
      }

      if (do_rama_restraints)
         // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
         flags = coot::ALL_RESTRAINTS;

      std::vector<coot::atom_spec_t> fixed_atom_specs = molecules[imol].get_fixed_atoms();

      // refinement goes a bit wonky if there are multiple occurrances of the same residue
      // in input residue vector, so let's filter out duplicates here
      //
      std::vector<mmdb::Residue *> residues;
      std::set<mmdb::Residue *> residues_set;
      std::set<mmdb::Residue *>::const_iterator it;
      for (std::size_t i=0; i<residues_in.size(); i++)
         residues_set.insert(residues_in[i]);
      residues.reserve(residues_set.size());
      for(it=residues_set.begin(); it!=residues_set.end(); ++it)
         residues.push_back(*it);

      // OK, so the passed residues are the residues in the graphics_info_t::molecules[imol]
      // molecule.  We need to do 2 things:
      //
      // convert the mmdb::Residue *s of the passed residues to the mmdb::Residue *s of residues mol
      //
      // and
      //
      // in create_mmdbmanager_from_res_vector() make sure that that contains the flanking atoms.
      // The create_mmdbmanager_from_res_vector() from this class is used, not coot::util
      //
      // The flanking atoms are fixed the passed residues are not fixed.
      // Keep a clear head.

      std::vector<std::string> residue_types = coot::util::residue_types_in_residue_vec(residues);
      // use try_dynamic_add()
      bool have_restraints = geom.have_dictionary_for_residue_types(residue_types, imol, cif_dictionary_read_number);
      cif_dictionary_read_number += residue_types.size();

      if (have_restraints) {

         std::string residues_alt_conf = alt_conf;
         imol_moving_atoms = imol;
         std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> > residues_mol_and_res_vec =
            create_mmdbmanager_from_res_vector(residues, imol, mol, residues_alt_conf);

         if (true) { // debug
            mmdb::Manager *residues_mol = residues_mol_and_res_vec.first;
            int imod = 1;
            mmdb::Model *model_p = residues_mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  std::cout << "DEBUG:: in generate_molecule_and_refine() residues_mol_and_res_vec mol: chain: "
                            << chain_p->GetChainID() << std::endl;
                  int nres = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     std::cout << "DEBUG:: in generate_molecule_and_refine() residues_mol_and_res_vec mol:   residue "
                               << coot::residue_spec_t(residue_p) << " residue "
                               << residue_p << " chain " << residue_p->chain << " index "
                               << residue_p->index << std::endl;
                  }
               }
            }
         }

         // We only want to act on these new residues and molecule, if
         // there is something there.
         //
         if (residues_mol_and_res_vec.first) {

            // Now we want to do an atom name check.  This stops exploding residues.
            //
            bool check_hydrogens_too_flag = false;
            std::pair<bool, std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > >
               icheck_atoms = geom.atoms_match_dictionary(imol, residues, check_hydrogens_too_flag, false);

            if (! icheck_atoms.first) {

               std::cout << "WARNING:: non-matching atoms! " << std::endl;

            } else {

               moving_atoms_have_hydrogens_displayed = true;
               if (! molecules[imol].hydrogen_atom_should_be_drawn())
                  moving_atoms_have_hydrogens_displayed = false;

               atom_selection_container_t local_moving_atoms_asc =
                  make_moving_atoms_asc(residues_mol_and_res_vec.first, residues);

               // 20221018-PE make_moving_atoms_graphics_object(imol, local_moving_atoms_asc); not today!

               int n_cis = coot::util::count_cis_peptides(local_moving_atoms_asc.mol);
               // moving_atoms_n_cis_peptides = n_cis; // 20221018-PE not today

               std::vector<std::pair<bool,mmdb::Residue *> > local_residues;  // not fixed.
               for (unsigned int i=0; i<residues_mol_and_res_vec.second.size(); i++)
                  local_residues.push_back(std::pair<bool, mmdb::Residue *>(0, residues_mol_and_res_vec.second[i]));

               moving_atoms_asc_type = NEW_COORDS_REPLACE;

               int imol_for_map = imol_refinement_map;
               clipper::Xmap<float> *xmap_p = dummy_xmap;

               if (is_valid_map_molecule(imol_for_map))
                  xmap_p = &molecules[imol_for_map].xmap;

               bool found_restraints_flag = make_last_restraints(local_residues,
                                                                 local_moving_atoms_asc.links,
                                                                 geom,
                                                                 residues_mol_and_res_vec.first,
                                                                 fixed_atom_specs,
                                                                 flags, use_map_flag, xmap_p);

               if (last_restraints) {
                  // 20220423-PE I can't do this here because setup_minimize() has not been called yet
                  // rr = last_restraints->get_refinement_results();
               }
               rr.found_restraints_flag = found_restraints_flag;

            }
         }
      } else {
         // we didn't have restraints for everything.
         //
         std::pair<int, std::vector<std::string> > icheck =
            check_dictionary_for_residue_restraints(imol, residues);
         if (icheck.first == 0) {
            std::cout << "WARNING:: <some info here about missing residue types> " << std::endl;
         }
      }
   }
   return rr;
}


int
molecules_container_t::mutate(int imol, const std::string &cid, const std::string &new_residue_type) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].mutate(residue_spec, new_residue_type);
      set_updating_maps_need_an_update(imol);
      // qstd::cout << "mutate status " << status << std::endl;
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

#include "coot-utils/blob-line.hh"

std::pair<bool, clipper::Coord_orth>
molecules_container_t::go_to_blob(float x1, float y1, float z1, float x2, float y2, float z2, float contour_level) {

   std::pair<bool, clipper::Coord_orth> p;

   clipper::Coord_orth p1(x1,y1,z1);
   clipper::Coord_orth p2(x2,y2,z2);

   // iterate through all the maps (another day)

   if (is_valid_map_molecule(imol_refinement_map)) {
      const clipper::Xmap<float> &xmap = molecules[imol_refinement_map].xmap;
      std::pair<bool, clipper::Coord_orth> pp = coot::find_peak_along_line_favour_front(p1, p2, contour_level, xmap);
      p = pp;
   }
   return p;
}


int
molecules_container_t::side_chain_180(int imol, const std::string &atom_cid) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t residue_spec(atom_spec);
      status = molecules[imol].side_chain_180(residue_spec, atom_spec.alt_conf, &geom);
      set_updating_maps_need_an_update(imol); // won't change much
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}

std::string
molecules_container_t::jed_flip(int imol, const std::string &atom_cid, bool invert_selection) {

   std::string message;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = atom_cid_to_atom_spec(imol, atom_cid);
      coot::residue_spec_t res_spec(atom_spec);
      std::string atom_name = atom_spec.atom_name;
      std::string alt_conf = "";
      message = molecules[imol].jed_flip(res_spec, atom_name, alt_conf, invert_selection, &geom);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return message;
}


void
molecules_container_t::coot_all_atom_contact_dots_instanced(mmdb::Manager *mol, int imol) {

   // 20221025-PE fill me later.

}

#include "ligand/ligand.hh"

int
molecules_container_t::add_waters(int imol_model, int imol_map) {

   int n_waters_added = 0;
   int ligand_water_n_cycles = 3;
   float ligand_water_to_protein_distance_lim_max = 3.2;
   float ligand_water_to_protein_distance_lim_min = 2.4;
   float ligand_water_variance_limit = 0.12;
   float sigma_cut_off = 1.0; // what is this?

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
	 coot::ligand lig;
	 int n_cycles = ligand_water_n_cycles; // 3 by default

	 // n_cycles = 1; // for debugging.

	 short int mask_waters_flag; // treat waters like other atoms?
	 // mask_waters_flag = g.find_ligand_mask_waters_flag;
	 mask_waters_flag = 1; // when looking for waters we should not
	 // ignore the waters that already exist.
	 // short int do_flood_flag = 0;    // don't flood fill the map with waters for now.

	 lig.import_map_from(molecules[imol_map].xmap,
			     molecules[imol_map].get_map_rmsd_approx());
	 // lig.set_masked_map_value(-2.0); // sigma level of masked map gets distorted
	 lig.set_map_atom_mask_radius(1.9); // Angstroms
	 lig.set_water_to_protein_distance_limits(ligand_water_to_protein_distance_lim_max,
						  ligand_water_to_protein_distance_lim_min);
	 lig.set_variance_limit(ligand_water_variance_limit);
	 lig.mask_map(molecules[imol_model].atom_sel.mol, mask_waters_flag);
	 // lig.output_map("masked-for-waters.map");
	 std::cout << "debug:: add_waters(): using n-sigma cut off " << sigma_cut_off << std::endl;

	 lig.water_fit(sigma_cut_off, n_cycles);

	 coot::minimol::molecule water_mol = lig.water_mol();
         molecules[imol_model].insert_waters_into_molecule(water_mol);
         n_waters_added = water_mol.count_atoms();
      }
   }
   return n_waters_added;
}


int
molecules_container_t::delete_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;

   // 20221025-PE Fill me later
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, res_no, ins_code);
      molecules[imol].delete_side_chain(res_spec);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::fill_side_chain(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   int status = 0;

   // 20221025-PE Fill me later
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, res_no, ins_code);
      // molecules[imol].fill_side_chain(res_spec);
      std::cout << "------------- fill fill_side_chain() here!" << std::endl;
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}


std::vector<std::string>
molecules_container_t::get_chains_in_model(int imol) const {

   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].chains_in_model();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

std::vector<std::pair<coot::residue_spec_t, std::string> >
molecules_container_t::get_single_letter_codes_for_chain(int imol, const std::string &chain_id) const {

   std::vector<std::pair<coot::residue_spec_t, std::string> > v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_single_letter_codes_for_chain(chain_id);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


std::vector<std::string>
molecules_container_t::get_residue_names_with_no_dictionary(int imol) const {

   std::vector<std::string> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_residue_names_with_no_dictionary(geom);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}



int
molecules_container_t::apply_transformation_to_atom_selection(int imol, const std::string &atoms_selection_cid,
                                                              int n_atoms, // for validation of the atom selection
                                                              float m00, float m01, float m02,
                                                              float m10, float m11, float m12,
                                                              float m20, float m21, float m22,
                                                              float c0, float c1, float c2, // the centre of the rotation
                                                              float t0, float t1, float t2) { // translation

   int n_atoms_moved = 0;
   if (is_valid_model_molecule(imol)) {
      clipper::Coord_orth rotation_centre(c0, c1, c2);
      clipper::Coord_orth t(t0, t1, t2);
      clipper::Mat33<double> m(m00, m01, m02, m10, m11, m12, m20, m21, m22);
      clipper::RTop_orth rtop_orth(m, t);
      n_atoms_moved = molecules[imol].apply_transformation_to_atom_selection(atoms_selection_cid, n_atoms, rotation_centre, rtop_orth);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return n_atoms_moved;

}


int
molecules_container_t::new_positions_for_residue_atoms(int imol, const std::string &residue_cid, std::vector<coot::molecule_t::moved_atom_t> &moved_atoms) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].new_positions_for_residue_atoms(residue_cid, moved_atoms);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;
}

int
molecules_container_t::new_positions_for_atoms_in_residues(int imol, const std::vector<coot::molecule_t::moved_residue_t> &moved_residues) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = molecules[imol].new_positions_for_atoms_in_residues(moved_residues);
      set_updating_maps_need_an_update(imol);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return status;

}

// put this in a new file molecules_container_validation.cc

#include "coot-utils/pepflip-using-difference-map.hh"


std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::pepflips_using_difference_map(int imol_coords, int imol_difference_map, float n_sigma) const {

   auto mmdb_to_clipper = [] (mmdb::Atom *at) {
      return clipper::Coord_orth(at->x, at->y, at->z);
   };

   std::vector<coot::molecule_t::interesting_place_t> v;

   if (is_valid_model_molecule(imol_coords)) {
      if (is_valid_map_molecule(imol_difference_map)) {
         if (molecules[imol_difference_map].is_difference_map_p()) {
            const clipper::Xmap<float> &diff_xmap = molecules[imol_difference_map].xmap;
            mmdb::Manager *mol = get_mol(imol_coords);
            coot::pepflip_using_difference_map pf(mol, diff_xmap);
            std::vector<coot::residue_spec_t> flips = pf.get_suggested_flips(n_sigma);
            for (std::size_t i=0; i<flips.size(); i++) {
               const auto &res_spec = flips[i];
               mmdb::Residue *residue_this_p = get_residue(imol_coords, res_spec);
               if (residue_this_p) {
                  coot::residue_spec_t res_spec_next =  res_spec.next();
                  mmdb::Residue *residue_next_p = get_residue(imol_coords, res_spec);
                  if (residue_next_p) {
                     std::string feature_type = "Difference Map Suggest Pepflip";
                     std::string label = "Flip this Bad boy: " + res_spec.format();
                     mmdb::Atom *at_1 = residue_this_p->GetAtom(" CA ");
                     mmdb::Atom *at_2 = residue_next_p->GetAtom(" CA ");
                     if (at_1 && at_2) {
                        clipper::Coord_orth pt_1 = mmdb_to_clipper(at_1);
                        clipper::Coord_orth pt_2 = mmdb_to_clipper(at_2);
                        clipper::Coord_orth pos = 0.5 * (pt_1 + pt_2);
                        float f = static_cast<float>(i)/static_cast<float>(flips.size());
                        float badness = 20.0 + 50.0 * (1.0 - f);
                        coot::molecule_t::interesting_place_t ip(feature_type, res_spec, pos, label);
                        ip.set_badness_value(badness);
                        v.push_back(ip);
                     }
                  }
               }
            }
         }
      }
   }
   std::cout << "DEBUG:: pepflips_using_difference_map() returns " << v.size() << " flips" << std::endl;
   return v;

}

//! @return a vector of residue specifiers for the ligand residues - the residue name is encoded
//! in the `string_user_data` data item of the residue specifier
std::vector<coot::residue_spec_t>
molecules_container_t::get_non_standard_residues_in_molecule(int imol) const {

   std::vector<coot::residue_spec_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_non_standard_residues_in_molecule();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


coot::simple_mesh_t
molecules_container_t::get_molecular_representation_mesh(int imol, const std::string &cid, const std::string &colour_scheme,
                                                         const std::string &style) const {

   coot::simple_mesh_t mesh;
   int status = 0;
   if (is_valid_model_molecule(imol)) {
      mesh = molecules[imol].get_molecular_representation_mesh(cid, colour_scheme, style);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return mesh;
}


//! get interesting places (does not work yet)
//! @return a vector of `validation_information_t`
std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::get_interesting_places(int imol, const std::string &mode) const {

   std::vector<coot::molecule_t::interesting_place_t> v;
   std::cout << "Nothing here yet" << std::endl;

   return v;
}


//! get Gaussian surface representation
coot::simple_mesh_t
molecules_container_t::get_gaussian_surface(int imol) const {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      mesh = molecules[imol].get_gaussian_surface();
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return mesh;

}


//! get chemical feaatures for the given residue
coot::simple_mesh_t
molecules_container_t::get_chemical_features_mesh(int imol, const std::string &cid) const {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      mesh = molecules[imol].get_chemical_features_mesh(cid, geom);
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return mesh;
}
