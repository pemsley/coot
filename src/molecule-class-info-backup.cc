
// #include "geometry/residue-and-atom-specs.hh"
#include <utility>
#include <cmath>
#include <math.h> // for fabsf
#include "molecule-class-info.h"

#include "utils/logging.hh"
extern logging logger;

/*! \brief Make a backup for a model molecule
 *
 * @param imol the model molecule index
 * @description a description that goes along with this back point
 */
int molecule_class_info_t::make_backup_checkpoint(const std::string &description) {

   make_backup(description);
   std::cout << "DEBUG:: mci make_backup_checkpoint() returns " << history_index << std::endl;
   return history_index;
}

/*! \brief Restore molecule from backup
 *
 * restore model @p imol to checkpoint backup @p backup_index
 *
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 */
int molecule_class_info_t::restore_to_backup_checkpoint(int backup_index) {

   std::string cwd = coot::util::current_working_dir();
   bool status = restore_from_backup(backup_index, cwd);
   if (status)
      return backup_index;
   else
      return -1;
}

/*! \brief Compare current model to backup
 *
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 * @return a list of residue specs for residues that have
 *         at least one atom in a different place.
 *   the first says is the backup_index was valid.
 */
std::pair<bool, std::vector<coot::residue_spec_t> > molecule_class_info_t::compare_current_model_to_backup(int backup_index) {

   bool status = false;
   std::vector<coot::residue_spec_t> rs;

   auto find_residue_differences = [] (mmdb::Manager *mol_1, mmdb::Manager *mol_2) {

      std::vector<mmdb::Residue *> moved_residues; // get converted to specs

      int imod = 1;
      mmdb::Model *model_1_p = mol_1->GetModel(imod);
      mmdb::Model *model_2_p = mol_2->GetModel(imod);
      if (model_1_p) {
         int n_chains_1 = model_1_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains_1; ichain++) {
            mmdb::Chain *chain_1_p = model_1_p->GetChain(ichain);
            std::string chain_id_1 = chain_1_p->GetChainID();
            int n_chains_2 = model_2_p->GetNumberOfChains();
            for (int jchain=0; jchain<n_chains_2; jchain++) {
               mmdb::Chain *chain_2_p = model_2_p->GetChain(ichain);
               std::string chain_id_2 = chain_2_p->GetChainID();
               if (chain_id_1 == chain_id_2) {
                  int n_res_1 = chain_1_p->GetNumberOfResidues();
                  for (int ires=0; ires<n_res_1; ires++) {
                     mmdb::Residue *residue_1_p = chain_1_p->GetResidue(ires);
                     if (residue_1_p) {
                        int res_no_1 = residue_1_p->GetSeqNum();
                        const char *ins_code_1 = residue_1_p->GetInsCode();
                        // eek!
                        int n_res_2 = chain_2_p->GetNumberOfResidues();
                        for (int jres=0; jres<n_res_2; jres++) {
                           mmdb::Residue *residue_2_p = chain_2_p->GetResidue(jres);
                           if (residue_2_p) {
                              int res_no_2 = residue_2_p->GetSeqNum();
                              const char *ins_code_2 = residue_1_p->GetInsCode();
                              if (res_no_1 == res_no_2) {
                                 if (strcmp(ins_code_1, ins_code_2) == 0) {
                                    bool diff_found = false;
                                    int n_atoms_1 = residue_1_p->GetNumberOfAtoms();
                                    for (int iat=0; iat<n_atoms_1; iat++) {
                                       mmdb::Atom *at_1 = residue_1_p->GetAtom(iat);
                                       if (! at_1->isTer()) {
                                          const char *atom_name_1 = at_1->GetAtomName();
                                          const char *atom_alt_conf_1 = at_1->altLoc;
                                          int n_atoms_2 = residue_2_p->GetNumberOfAtoms();
                                          for (int jat=0; jat<n_atoms_2; jat++) {
                                             mmdb::Atom *at_2 = residue_2_p->GetAtom(jat);
                                             if (! at_2->isTer()) {
                                                const char *atom_name_2 = at_2->GetAtomName();
                                                const char *atom_alt_conf_2 = at_2->altLoc;

                                                if (strcmp(atom_name_1, atom_name_2) == 0) {
                                                   if (strcmp(atom_alt_conf_1, atom_alt_conf_2) == 0) {
                                                      float delta_x = at_1->x - at_2->x;
                                                      float delta_y = at_1->y - at_2->y;
                                                      float delta_z = at_1->z - at_2->z;
                                                      if (fabsf(delta_x) > 0.01) {
                                                         if (fabsf(delta_y) > 0.01) {
                                                            if (fabsf(delta_z) > 0.01) {
                                                               if (false)
                                                                  std::cout << "DEBUG:: atom " << coot::atom_spec_t(at_1) << " " << coot::atom_spec_t(at_2) << " "
                                                                            << delta_x << " " << delta_y << " " << delta_z << std::endl;
                                                               if (std::find(moved_residues.begin(), moved_residues.end(), residue_1_p) == moved_residues.end()) {
                                                                  moved_residues.push_back(residue_1_p);
                                                               }
                                                               diff_found = true;
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                             if (diff_found) break;
                                          }
                                       }
                                       if (diff_found) break;
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
      }
      std::vector<coot::residue_spec_t> res_specs(moved_residues.size());
      for (unsigned int i=0; i<moved_residues.size(); i++) {
         res_specs[i] = coot::residue_spec_t(moved_residues[i]);
      }
      return res_specs;
   };

   if (atom_sel.mol) {
      if (int(history_filename_vec.size()) > backup_index) {
         // there should be a log function for unsigned int, int.
         logger.log(log_t::INFO, logging::function_name_t(__FUNCTION__),
                    {logging::ltw(history_filename_vec.size()), logging::ltw(backup_index)});

         if (backup_index < int(history_filename_vec.size()) && backup_index >= 0) {
            std::string filename = history_filename_vec[backup_index].backup_file_name;

            bool use_gemmi = false; // for now
            atom_selection_container_t asc = get_atom_selection(filename, use_gemmi);
            if (asc.read_success) {
               rs = find_residue_differences(asc.mol, atom_sel.mol);
               status = true;
               asc.clear_up();
            }
         }
      } else {
         logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__),
                    {logging::ltw(std::string("bad backup_index")),
                     logging::ltw(history_filename_vec.size()), logging::ltw(backup_index)});
      }
   }
   return std::make_pair(status, rs);
}

/*! \brief Get backup info
 * 
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 * @return a Python list of the given description (str)
 *         and a timestamp (str).
 */
coot::backup_file_info_t molecule_class_info_t::get_backup_info(int backup_index) {
   
   coot::backup_file_info_t bfi;
   if (backup_index < int(history_filename_vec.size()) && backup_index >= 0) {
      coot::backup_file_info_t bfi = history_filename_vec[backup_index];
   }
   return bfi;
}


void molecule_class_info_t::print_backup_history_info() const {

   std::cout << "DEBUG:: =========== history has " << history_filename_vec.size() << " entries ================"
             << std::endl;
   for (unsigned int i=0; i<history_filename_vec.size(); i++) {
      const auto &history = history_filename_vec[i];
      std::cout << "DEBUG:: " << i << " " << history.imol << " " << history.backup_file_name << " "
                << history.name << " description: " << history.description << " time: " << history.get_timespec_string() << std::endl;
   }



}
