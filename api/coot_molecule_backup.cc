
#include <sys/stat.h>

#include "utils/coot-utils.hh"
#include "coot_molecule.hh"
#include "molecules_container.hh" // for the 

#include "coords/mmdb.hh" // for write_atom_selection_file()

std::string
coot::molecule_t::make_backup(const std::string &modification_info_string) {

   std::string info_message; // non-empty on failture
   info_message = modification_info.make_backup(atom_sel.mol, modification_info_string);
   return info_message;
}


int
coot::molecule_t::undo() {

   int status = 0;
   mmdb::Manager *mol_new = modification_info.undo();
   if (! mol_new) {
      std::cout << "ERROR:: undo failed" << std::endl;
   } else {
      atom_sel.clear_up();
      atom_sel = make_asc(mol_new);
   }
   return status;
}

int
coot::molecule_t::redo() {

   int status = 0;
   mmdb::Manager *mol_new = modification_info.redo();
   if (! mol_new) {
      std::cout << "ERROR:: undo failed" << std::endl;
   } else {
      atom_sel.clear_up();
      atom_sel = make_asc(mol_new);
   }
   return status;

}



// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
//                    modification_info_t
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------


std::string
coot::molecule_t::modification_info_t::get_backup_file_name_from_index(int idx) const {

   std::string s;

   auto get_extension = [this] () {
      std::string ext = ".pdb";
      if (this->is_mmcif_flag)
         ext = ".cif";
      return ext;
   };

   std::string fn = mol_name + "-" + get_index_string(idx) + get_extension();
   if (backup_dir.empty()) {
      s = fn;
   } else {
      util::create_directory(backup_dir); // maybe
      s = util::append_dir_file(backup_dir, fn);
   }
   return s;
}

bool
coot::molecule_t::modification_info_t::have_unsaved_changes() const {

   return true;

}

std::string
coot::molecule_t::modification_info_t::make_backup(mmdb::Manager *mol, const std::string &modification_info_string) {

   std::string message;
   int index = save_info.size();
   std::string fn = get_backup_file_name_from_index(index);

   if (is_mmcif_flag) {

      // from write_atom_selection_file():

      // WriteCIFASCII() seems to duplicate the atoms (maybe related to aniso?)
      // So let's copy the molecule and throw away the copy, that way we don't
      // duplicate the atoms in the original molecule.
      
      mmdb::Manager *mol_copy  = new mmdb::Manager;
      mol_copy->Copy(mol, mmdb::MMDBFCM_All);
      int ierr = mol_copy->WriteCIFASCII(fn.c_str());
      delete mol_copy;

      if (ierr != mmdb::Error_NoError) {
         std::cout << "get the error message " << fn << std::endl;
      }
      save_info.push_back(save_info_t(fn, modification_info_string));
      modification_index = save_info.size();

   } else {

      int ierr = mol->WritePDBASCII(fn.c_str());

      if (ierr != mmdb::Error_NoError) {
         std::cout << "get the error message " << fn << std::endl;
      }
   }

   return message;

}

mmdb::Manager *
coot::molecule_t::modification_info_t::save_info_t::get_mol() {

   mmdb::Manager *MMDBManager = nullptr;
   MMDBManager = new mmdb::Manager;
   MMDBManager->SetFlag ( mmdb::MMDBF_IgnoreBlankLines |
                          mmdb::MMDBF_IgnoreDuplSeqNum |
                          mmdb::MMDBF_IgnoreNonCoorPDBErrors |
                          mmdb::MMDBF_IgnoreHash |
                          mmdb::MMDBF_IgnoreRemarks);
   mmdb::ERROR_CODE err = MMDBManager->ReadCoorFile(file_name.c_str());
   if (err != mmdb::Error_NoError) {
      int  error_count;
      char error_buf[500];
      std::cout << "ERROR::" << file_name << " " << mmdb::GetErrorDescription(err) << std::endl;
      MMDBManager->GetInputBuffer(error_buf, error_count);
      if (error_count >= 0)
         std::cout << "ERROR:: LINE #" << error_count << "\n     " << error_buf << std::endl;
   }
   return MMDBManager;
}


mmdb::Manager *
coot::molecule_t::modification_info_t::undo() {

   mmdb::Manager *MMDBManager = nullptr;
   std::cout << "undo" << std::endl;
   int idx = modification_index - 1;
   if (idx >= 0) {
      if (idx < int(save_info.size())) {
         MMDBManager = save_info[idx].get_mol();
      }
   }
   return MMDBManager;
}


mmdb::Manager *
coot::molecule_t::modification_info_t::redo() {

   mmdb::Manager *MMDBManager = nullptr;
   std::cout << "redo" << std::endl;
   int idx = modification_index + 1;
   if (idx >= 0) {
      if (idx < int(save_info.size())) {
         MMDBManager = save_info[idx].get_mol();
      }
   }
   return MMDBManager;
}

