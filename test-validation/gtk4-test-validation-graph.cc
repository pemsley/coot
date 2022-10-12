
// use: ./gtk4-test-validation-graph tutorial-modern.pdb rnasa-1.8-all_refmac1.mtz
// there files are in the data directory.

#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/coot-map-utils.hh"
#include "validation-information.hh"

bool
read_mtz(const std::string &file_name,
         const std::string &f, const std::string &phi, const std::string &weight,
         bool use_weight, bool is_a_difference_map,
         clipper::Xmap<float> *xmap_p) {

   bool status = coot::util::map_fill_from_mtz(xmap_p, file_name, f, phi, weight, use_weight, is_a_difference_map);
   return status;
}



coot::validation_information_t
density_fit_analysis(const std::string &pdb_file_name, const std::string &mtz_file_name) {

   coot::validation_information_t r;

   // fill these
   mmdb::PResidue *SelResidues = 0;
   int nSelResidues = 0;

   auto atom_sel = get_atom_selection(pdb_file_name, true, false, false);
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

   clipper::Xmap<float> xmap;
   bool status = read_mtz(mtz_file_name, "FWT", "PHWT", "W", 0, 0, &xmap);

   if (! status) {

      std::cout << "Bad mtz file read " << mtz_file_name << std::endl;

   } else {

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
         coot::residue_spec_t res_spec(residue_p);
         mmdb::PAtom *residue_atoms=0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         double residue_density_score =
            coot::util::map_score(residue_atoms, n_residue_atoms, xmap, 1);
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
   return r;
}


int main(int argc, char **argv) {

   if (argc > 2) {
      std::string pdb_file_name = argv[1];
      std::string mtz_file_name = argv[2];
      coot::validation_information_t vi = density_fit_analysis(pdb_file_name, mtz_file_name);

      // now do something (i.e. make a pretty interactive graph) with vi.

      for (const auto &cvi : vi.cviv) {
         std::cout << "Chain " << cvi.chain_id << std::endl;
         for (const auto &ri : cvi.rviv) {
            std::cout << " Residue " << ri.residue_spec << " " << ri.distortion << std::endl;
         }
      }
   } else {
      std::cout << "Two commandline args needed.\n";
      return 2;
   }

   return 0;
}
