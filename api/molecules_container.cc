
#include "molecules_container.hh"
#include "ideal/pepflip.hh"
#include "coot-utils/coot-map-utils.hh"

bool
molecules_container_t::is_valid_model_molecule(int imol) {
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
molecules_container_t::is_valid_map_molecule(int imol) {
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
molecules_container_t::flipPeptide(int imol, const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      result = molecules[imol].flipPeptide(rs, alt_conf);
   }
   return result;
}

int
molecules_container_t::flipPeptide(int imol, const std::string &cid, const std::string &alt_conf) {

   int result = 0;
   if (is_valid_model_molecule(imol)) {
      auto &m = molecules[imol];
      std::pair<bool, coot::residue_spec_t> rs = m.cid_to_residue_spec(cid);
      if (rs.first)
         result = molecules[imol].flipPeptide(rs.second, alt_conf);
   }
   return result;
}


int
molecules_container_t::read_pdb(const std::string &file_name) {

   int status = -1;
   atom_selection_container_t asc = get_atom_selection(file_name);
   if (asc.read_success) {
      molecules.push_back(coot::molecule_t(asc));
      status = molecules.size() -1;
   }
   return status;
}

int
molecules_container_t::read_mtz(const std::string &file_name,
                                const std::string &f, const std::string &phi, const std::string &weight,
                                bool use_weight, bool is_a_difference_map) {

   int imol = -1;

   coot::molecule_t m;
   bool status = coot::util::map_fill_from_mtz(&m.xmap, file_name, f, phi, weight, use_weight, is_a_difference_map, 0, 0);
   if (status) {
      molecules.push_back(m);
      imol = molecules.size() -1;
   }
   return imol;
}


coot::validation_information_t
molecules_container_t::density_fit_analysis(int imol_model, int imol_map) {

   coot::validation_information_t r;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_model_molecule(imol_map)) {
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
            r.add_residue_valiation_informtion(rvi, chain_id);
         }
         atom_sel.mol->DeleteSelection(selHnd);
      }
   }
   return r;
}
