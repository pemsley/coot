
#include "side-chain.hh"
#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-coord-extras.hh"

#include "chi-angles.hh"

#include "compat/coot-sysdep.h"

void
coot::do_180_degree_side_chain_flip(const coot::residue_spec_t &spec,
                                    const std::string &alt_conf,
                                    mmdb::Manager *mol,
                                    coot::protein_geometry *geom_p) {


   int nSelResidues = 0;
   mmdb::PResidue *SelResidues = NULL;
   int selnd = mol->NewSelection();
   mol->Select(selnd, mmdb::STYPE_RESIDUE, 0,
               spec.chain_id.c_str(),
               spec.res_no, spec.ins_code.c_str(),
               spec.res_no, spec.ins_code.c_str(),
               "*", "*", "*", "*", mmdb::SKEY_NEW);
   mol->GetSelIndex(selnd, SelResidues, nSelResidues);
   if (nSelResidues > 0 ) {
      mmdb::Residue *residue = SelResidues[0];
      std::string resname = residue->GetResName();

      int nth_chi = -1; // unset

      // if (resname == "ARG") nth_chi = 4;
      if (resname == "ARG") nth_chi = 5;
      if (resname == "ASP") nth_chi = 2;
      if (resname == "ASN") nth_chi = 2;
      if (resname == "CYS") nth_chi = 1;
      if (resname == "GLN") nth_chi = 3;
      if (resname == "GLU") nth_chi = 3;
      if (resname == "PHE") nth_chi = 2;
      if (resname == "HIS") nth_chi = 2;
      if (resname == "SER") nth_chi = 1;
      if (resname == "THR") nth_chi = 1;
      if (resname == "VAL") nth_chi = 1;
      if (resname == "TRP") nth_chi = 2;
      if (resname == "TYR") nth_chi = 2;

      mmdb::Atom **residue_atoms = NULL;
      int nResidueAtoms;
      residue->GetAtomTable(residue_atoms, nResidueAtoms);

      if (nth_chi != -1) {
         mmdb::Residue *residue_copy =
            coot::util::deep_copy_this_residue_add_chain(residue, alt_conf, 0, 0);

         // Which atoms have we got in residue_copy?
         int n_atom_residue_copy;
         mmdb::PAtom *residue_atoms_copy = 0;
         residue_copy->GetAtomTable(residue_atoms_copy, n_atom_residue_copy);
         //             for (int iat=0; iat<n_atom_residue_copy; iat++)
         //                std::cout << residue_atoms_copy[iat] << std::endl;

         // check that the N comes before the CA and the CA comes before
         // the CB (if it has one).
         if (coot::util::is_standard_amino_acid_name(resname)) {
            bool needs_reordering = false;
            int idx_N  = -1;
            int idx_CA = -1;
            int idx_CB = -1;
            for (int iat=0; iat<n_atom_residue_copy; iat++) {
               mmdb::Atom *at = residue_atoms_copy[iat];
               std::string at_name(at->GetAtomName());
               if (at_name == " N  ") {  // PDBv3 FIXME
                  idx_N = iat;
               }
               if (at_name == " CA ") {  // PDBv3 FIXME
                  idx_CA = iat;
               }
               if (at_name == " CB ") {  // PDBv3 FIXME
                  idx_CB = iat;
               }
            }
            if (idx_N != -1) {
               if (idx_CA != -1) {
                  if (idx_N > idx_CA) {
                     needs_reordering = true;
                  }
               }
            }
            if (idx_CB != -1) {
               if (idx_CA != -1) {
                  if (idx_CA > idx_CB) {
                     needs_reordering = true;
                  }
               }
            }
            if (needs_reordering) {
               coot::put_amino_acid_residue_atom_in_standard_order(residue_copy);
            }
         }

         chi_angles chi_ang(residue_copy, 0);
         std::vector<std::vector<int> > contact_indices(n_atom_residue_copy);
         bool add_reverse_contacts = false;
         contact_indices = coot::util::get_contact_indices_from_restraints(residue_copy, geom_p, 1,
                                                                           add_reverse_contacts);
         double diff = 180.0;
         std::pair<short int, float> istat = chi_ang.change_by(nth_chi, diff, contact_indices);

         if (istat.first) { // failure
            std::cout << "Failure to flip" << std::endl;
         } else {

            // OK, we need transfer the coordinates of the
            // altconfed atoms of residue_copy to residue:
            //
            for (int iatc=0; iatc<n_atom_residue_copy; iatc++) {
               // std::cout << residue_atoms_copy[iat] << std::endl;
               std::string atom_copy_altconf = residue_atoms_copy[iatc]->altLoc;
               if (atom_copy_altconf == alt_conf) {
                  // we need to find this atom in residue
                  std::string atom_copy_name = residue_atoms_copy[iatc]->name;
                  for (int iato=0; iato<nResidueAtoms; iato++) {
                     std::string orig_atom_altconf = residue_atoms[iato]->altLoc;
                     std::string orig_atom_name    = residue_atoms[iato]->name;
                     if (orig_atom_name == atom_copy_name) {
                        if (atom_copy_altconf == orig_atom_altconf) {
                           //                               std::cout << "DEBUG:: copying coords from "
                           //                                         << residue_atoms_copy[iatc] << std::endl;
                           residue_atoms[iato]->x = residue_atoms_copy[iatc]->x;
                           residue_atoms[iato]->y = residue_atoms_copy[iatc]->y;
                           residue_atoms[iato]->z = residue_atoms_copy[iatc]->z;
                        }
                     }
                  }
               }
            }
            // Now let's get rid of residue_copy:
            delete residue_copy;
            residue_copy = 0;
         }
      }
   }
 }
