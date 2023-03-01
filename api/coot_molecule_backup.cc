
#include <sys/stat.h>

#include "coot_molecule.hh"
#include "molecules_container.hh" // for the 

#include "coords/mmdb.hh" // for write_atom_selection_file()

std::string
coot::molecule_t::make_backup() {

   if (false) {
      std::cout << "start make_backup() for molecule " << imol_no << std::endl;
      std::cout << "start make_backup() for molecule with " << atom_sel.n_selected_atoms << " atoms " << std::endl;
   }

   if (molecules_container_t::make_backups_flag == false)
      return std::string("No Backups");

   std::string info_message;

   auto make_maybe_backup_dir = [] (const std::string &backup_dir) {
      return util::create_directory(backup_dir);
   };

   bool backup_this_molecule = true; // if needed, give user control later
   bool backup_compress_files_flag = false;

   if (backup_this_molecule) {
      std::string backup_dir("coot-backup");

      //shall we use the environment variable instead?
      char *env_var = getenv("COOT_BACKUP_DIR");

#ifdef EMSCRIPTEN

#else

      if (env_var) {
         struct stat buf;

         // we better debackslash the directory (for windows)
         std::string tmp_dir = env_var;
         tmp_dir = coot::util::intelligent_debackslash(tmp_dir);
         int err = stat(tmp_dir.c_str(), &buf);

         if (!err) {
            if (! S_ISDIR(buf.st_mode)) {
               env_var = NULL;
            }
         } else {
            env_var = NULL;
         }
      }
#endif

      if (env_var)
         backup_dir = env_var;

      if (atom_sel.mol) {
         int dirstat = make_maybe_backup_dir(backup_dir);

         if (dirstat != 0) {
            // fallback to making a directory in $HOME
            std::string home_dir = coot::get_home_dir();
            if (! home_dir.empty()) {
               backup_dir = coot::util::append_dir_dir(home_dir, "coot-backup");
               dirstat = make_maybe_backup_dir(backup_dir);
               if (dirstat != 0) {
                  std::cout << "WARNING:: backup directory "<< backup_dir
                            << " failure to exist or create" << std::endl;
               } else {
                  std::cout << "INFO using backup directory " << backup_dir << std::endl;
               }
            } else {
               std::cout << "WARNING:: backup directory "<< backup_dir
                         << " failure to exist or create" << std::endl;
            }
         }

         if (dirstat == 0) {
            // all is hunkey-dorey.  Directory exists.

            std::string backup_file_name = get_save_molecule_filename(backup_dir);
            // std::cout << "INFO:: make_backup() backup file name " << backup_file_name << std::endl;

            mmdb::byte gz;
            if (backup_compress_files_flag) {
               gz = mmdb::io::GZM_ENFORCE;
            } else {
               gz = mmdb::io::GZM_NONE;
            }

            // Writing out a modified binary mmdb like this results in the
            // file being unreadable (crash in mmdb read).
            //
            int istat;
            if (! is_from_shelx_ins_flag) {
               bool write_as_cif = false;
               if (coot::is_mmcif_filename(name))
                  write_as_cif = true;

               istat = write_atom_selection_file(atom_sel, backup_file_name, write_as_cif, gz);

               if (true) { // 20221021-PE this should not be needed
                  struct stat buf;
                  int err = stat(backup_file_name.c_str(), &buf);
                  if (err == 0) {
                     // std::cout << "DEBUG:: in make_backup() " << backup_file_name << " confirmed as existing" << std::endl;
                  } else {
                     std::cout << "DEBUG:: in make_backup() " << backup_file_name << " does not exist!" << std::endl;
                  }
               }

               // WriteMMDBF returns 0 on success, else mmdb:Error_CantOpenFile (15)
               if (istat) {
                  std::string warn;
                  warn = "WARNING:: WritePDBASCII failed! Return status ";
                  warn += istat;
                  // g.info_dialog_and_text(warn);
                  info_message = warn;
               }
            } else {
               std::pair<int, std::string> p = write_shelx_ins_file(backup_file_name);
               istat = p.first;
            }

            save_history_file_name(backup_file_name);
            // 20221016-PE old history counting - now use save_info.
            // if (history_index == max_history_index)
            // max_history_index++;
            // history_index++;
         }
      } else {
         std::cout << "WARNING:: BACKUP:: Ooops - no atoms to backup for this empty molecule"
                   << std::endl;
      }
   } else {
      // Occasionally useful but mostly tedious...
      // std::cout << "INFO:: backups turned off on this molecule"
      // << std::endl;
   }
   return info_message;
}

