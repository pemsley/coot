
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot_molecule.hh"
#include "ideal/pepflip.hh"
#include "rama-plot-phi-psi.hh"

bool
coot::molecule_t::is_valid_model_molecule() const {

   bool status = false;
   if (atom_sel.mol)
      status = true;
   return status;

}

bool
coot::molecule_t::is_valid_map_molecule() const {

   bool status = false;
   if (! xmap.is_null()) {
      status = true;
   }
   return status;
}

std::pair<bool, coot::residue_spec_t>
coot::molecule_t::cid_to_residue_spec(const std::string &cid) {

   bool status = false;
   coot::residue_spec_t rs;
   if (atom_sel.mol) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::Residue **SelResidues;
      int nSelResidues = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_CHAIN, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      atom_sel.mol->DeleteSelection(selHnd);
      if (nSelResidues > 0) {
         mmdb::Residue *residue_p = SelResidues[0];
         coot::residue_spec_t rs_inner(residue_p);
         rs = rs_inner;
         status = true;
      }
   }
   return std::make_pair(status, rs);
}

void
coot::molecule_t::make_backup() {

   // nothing done yet.

}

bool
coot::molecule_t::moving_atom_matches(mmdb::Atom *at, int this_mol_index_maybe) const {

   bool matches = false;
   if (atom_sel.n_selected_atoms > 0) {
      if (this_mol_index_maybe >= atom_sel.n_selected_atoms) {
         return false;
      } else {
         std::string atom_name_mov = at->name;
         std::string ins_code_mov  = at->GetInsCode();
         std::string alt_conf_mov  = at->altLoc;
         std::string chain_id_mov  = at->GetChainID();
         int resno_mov = at->GetSeqNum();

         std::string atom_name_ref = atom_sel.atom_selection[this_mol_index_maybe]->name;
         std::string ins_code_ref  = atom_sel.atom_selection[this_mol_index_maybe]->GetInsCode();
         std::string alt_conf_ref  = atom_sel.atom_selection[this_mol_index_maybe]->altLoc;
         std::string chain_id_ref  = atom_sel.atom_selection[this_mol_index_maybe]->GetChainID();
         int resno_ref = atom_sel.atom_selection[this_mol_index_maybe]->GetSeqNum();

         if (atom_name_ref == atom_name_mov) {
            if (ins_code_ref == ins_code_mov) {
               if (resno_ref == resno_mov) {
                  if (alt_conf_ref == alt_conf_mov) {
                     if (chain_id_mov == chain_id_ref) { // 20170612 extra condition added, Oliver Clarke bug
                        matches = true;
                     }
                  }
               }
            }
         }
      }
   }
   return matches;
}


int
coot::molecule_t::full_atom_spec_to_atom_index(const coot::atom_spec_t &atom_spec) const {

   return full_atom_spec_to_atom_index(atom_spec.chain_id,
                                       atom_spec.res_no,
                                       atom_spec.ins_code,
                                       atom_spec.atom_name,
                                       atom_spec.alt_conf);

}

// return -1 on no atom found.
int
coot::molecule_t::full_atom_spec_to_atom_index(const std::string &chain,
                                               int resno,
                                               const std::string &insertion_code,
                                               const std::string &atom_name,
                                               const std::string &alt_conf) const {

   int iatom_index = -1;

   // some protection for null molecule.
   if (! atom_sel.mol) {
      std::cout << "ERROR:: null molecule for molecule number "
                << imol_no << " pointer: " << atom_sel.mol
                << " (in full_atom_spec_to_atom_index)" << std::endl;
      return -1;
   }

   int selHnd = atom_sel.mol->NewSelection();
   int idx = 0;

   atom_sel.mol->SelectAtoms(selHnd, 0, chain.c_str(),
                            resno, insertion_code.c_str(), // start, insertion code
                            resno, insertion_code.c_str(), // end, insertion code
                            "*", // residue name
                            atom_name.c_str(),
                            "*", // elements
                            alt_conf.c_str()); // alt locs

   int nSelAtoms;
   mmdb::PPAtom local_SelAtom;
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (false)
      std::cout << "DEBUG:: full_atom_spec_to_atom_index() for :" << chain << ": "
                << resno << " :" << insertion_code << ": :"
                << atom_name << ": :" << alt_conf << ": finds " << nSelAtoms <<  " atoms\n";

   if (nSelAtoms == 0) {

      std::cout << "WARNING:: full_atom_spec_to_atom_index() Could not find "
                << "\"" << atom_name << "\"," << "\"" << alt_conf  << "\"" << "/"
                << resno << insertion_code << "/" << chain << " in this molecule: ("
                <<  imol_no << ") " << name << std::endl;

      int selHnd2 = atom_sel.mol->NewSelection(); // d

      atom_sel.mol->SelectAtoms(selHnd2, 0,
                                chain.c_str(),
                                resno, "*", // start, insertion code
                                resno, "*", // end, insertion code
                                "*", // residue name
                                "*", // atom name
                                "*", // elements
                                "*"); // alt locs

      atom_sel.mol->GetSelIndex(selHnd2, local_SelAtom, nSelAtoms);

      if (false) { // debugging.
         std::cout << "There were " << nSelAtoms << " atoms in that residue:\n";
         std::cout << "debgu:: full_atom_spec_to_atom_index() resno " << resno
                   << " (cf MinInt4) " << mmdb::MinInt4 << "\n";
         if (resno == mmdb::MinInt4) {
            std::cout << "      residue with resno MinInt4\n";
         } else {
            for (int i=0; i<nSelAtoms; i++) {
               std::cout << "      " << local_SelAtom[i] << "\n";
            }
         }
      }

      atom_sel.mol->DeleteSelection(selHnd2);

   } else {

      if (nSelAtoms != 1) {
         // the wildcard atom selection case "*HO2"
         short int found = 0;
         for (int i=0; i<nSelAtoms; i++) {
            if (std::string(local_SelAtom[i]->GetChainID()) == chain) {
               if (local_SelAtom[i]->residue->seqNum == resno) {
                  if (std::string(local_SelAtom[i]->GetInsCode()) == insertion_code) {
                     if (std::string(local_SelAtom[i]->name) == atom_name) {
                        if (std::string(local_SelAtom[i]->altLoc) == alt_conf) {
                           found = 0;
                           idx = i;
                           break;
                        }
                     }
                  }
               }
            }
         }
      }

      int iatom_index_udd = -1;
      int ic;
      if (local_SelAtom[idx]->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == mmdb::UDDATA_Ok) {
         iatom_index_udd = ic;
      }
      iatom_index = iatom_index_udd;
   }
   atom_sel.mol->DeleteSelection(selHnd); // Oh dear, this should have
                                          // been in place for years
                                          // (shouldn't it?) 20071121
   return iatom_index;
}




// helper function for above function
bool
coot::molecule_t::movable_atom(mmdb::Atom *mol_atom, bool replace_coords_with_zero_occ_flag) const {

   bool m = true;

   if ((mol_atom->occupancy < 0.0001) &&
       (mol_atom->occupancy > -0.0001))
      if (replace_coords_with_zero_occ_flag == 0)
         m = 0; // zero occupancy and "dont move zero occ atoms is set"
   return m;
}


// We just added a new atom to a residue, now we need to adjust the
// occupancy of the other atoms (so that we don't get residues with
// atoms whose occupancy is greater than 1.0 (Care for SHELX molecule?)).
// at doesn't have to be in residue.
//
// Perhaps this can be a utility function?
//
void
coot::molecule_t::adjust_occupancy_other_residue_atoms(mmdb::Atom *at,
                                                       mmdb::Residue *residue,
                                                       short int force_sum_1_flag) {

   if (!is_from_shelx_ins_flag) {
      int nResidueAtoms;
      mmdb::PPAtom ResidueAtoms = 0;
      residue->GetAtomTable(ResidueAtoms, nResidueAtoms);
      float new_atom_occ = at->occupancy;
      std::string new_atom_name(at->name);
      std::string new_atom_altconf(at->altLoc);
      std::vector<mmdb::Atom *> same_name_atoms;
      float sum_occ = 0;
      for (int i=0; i<nResidueAtoms; i++) {
         std::string this_atom_name(ResidueAtoms[i]->name);
         std::string this_atom_altloc(ResidueAtoms[i]->altLoc);
         if (this_atom_name == new_atom_name) {
            if (this_atom_altloc != new_atom_altconf) {
               same_name_atoms.push_back(ResidueAtoms[i]);
               sum_occ += ResidueAtoms[i]->occupancy;
            }
         }
      }

      //
      if (sum_occ > 0.01) {
         if (same_name_atoms.size() > 0) {
            float other_atom_occ_sum = 0.0;
            for (unsigned int i=0; i<same_name_atoms.size(); i++)
               other_atom_occ_sum += same_name_atoms[i]->occupancy;

            float remainder = 1.0 - new_atom_occ;
            float f = remainder/other_atom_occ_sum;
            for (unsigned int i=0; i<same_name_atoms.size(); i++) {
               if (0) // debug
                  std::cout << "debug " << same_name_atoms[i]
                            << " mulitplying occ " << same_name_atoms[i]->occupancy
                            << " by " << remainder << "/" << other_atom_occ_sum << "\n";
               same_name_atoms[i]->occupancy *= f;
            }
         }
      }
   }
}


// Put the regularization results back into the molecule:
//
//// Recall that regularized_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
void
coot::molecule_t::replace_coords(const atom_selection_container_t &asc,
                                 bool change_altconf_occs_flag,
                                 bool replace_coords_with_zero_occ_flag) {


   int n_atom = 0;
   int tmp_index;
   bool debug = false;
   float add_alt_conf_new_atoms_occupancy = 0.5; // was a static in graphics_info_t

   make_backup();

   // debug::
   if (debug) {
      std::cout << "DEBUG:: --------------- replace_coords replacing "
                << asc.n_selected_atoms << " atoms " << std::endl;
      for (int i=0; i<asc.n_selected_atoms; i++) {
         mmdb::Atom *atom = asc.atom_selection[i];
         bool is_ter_state = atom->isTer();
         std::cout << "DEBUG:: in replace_coords, intermediate atom: chain-id :"
                   << atom->residue->GetChainID() <<  ": "
                   << atom->residue->seqNum << " inscode \""
                   << atom->GetInsCode() << "\" name \""
                   << atom->name << "\" altloc \""
                   << atom->altLoc << "\" occupancy: "
                   << atom->occupancy << " :"
                   << " TER state: " << is_ter_state << std::endl;
      }
   }

   // For each atom in the new set of atoms:
   //
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idx = -1;
      mmdb::Atom *atom = asc.atom_selection[i];
      if (! atom->isTer()) {

         if (debug) { // debug
            std::cout << "considering replacement for selected atom " << coot::atom_spec_t(atom) << std::endl;

            //
            // idx = atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
            // atom->residue->seqNum, std::string(atom->name));

         }

         if (asc.UDDOldAtomIndexHandle >= 0) { // OK for fast atom indexing

            if (debug)
               std::cout << "... OK for fast atom indexing, asc.UDDOldAtomIndexHandle: "
                         << asc.UDDOldAtomIndexHandle
                         << " for atom " << coot::atom_spec_t(atom) << std::endl;

            if (atom->GetUDData(asc.UDDOldAtomIndexHandle, tmp_index) == mmdb::UDDATA_Ok) {

               if (debug)
                  std::cout << "OK, good GetUDData() for atom " << coot::atom_spec_t(atom) << std::endl;
               if (tmp_index >= 0) {
                  if (moving_atom_matches(atom, tmp_index)) {
                     // std::cout << "      DEBUG:: successfully found old atom index" << std::endl;
                     idx = tmp_index;
                  } else {
                     // std::cout << "DEBUG:: atom index mismatch" << std::endl;
                     idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                                        atom->residue->seqNum,
                                                        std::string(atom->GetInsCode()),
                                                        std::string(atom->name),
                                                        std::string(atom->altLoc));
                     // std::cout << "DEBUG:: full_atom_spec_to_atom_index gives index: " << idx << std::endl;
                  }
               } else {
                  // This shouldn't happen.
                  std::cout << "Good Handle, bad index found for old atom: specing" << std::endl;
                  idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                                     atom->residue->seqNum,
                                                     std::string(atom->GetInsCode()),
                                                     std::string(atom->name),
                                                     std::string(atom->altLoc));
               }
            } else {

               std::cout << "ERROR:: non-bad handle (" << asc.UDDOldAtomIndexHandle
                         <<  "), but bad GetUDData() for atom " << coot::atom_spec_t(atom) << std::endl;

            }
         } else {

            if (debug)
               std::cout << "DEBUG:: asc.UDDOldAtomIndexHandle is "
                         << asc.UDDOldAtomIndexHandle << " using full atom spec to atom index..."
                         << std::endl;

            idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                               atom->residue->seqNum,
                                               std::string(atom->GetInsCode()),
                                               std::string(atom->name),
                                               std::string(atom->altLoc));
            if (idx == -1) {
               std::cout << "DEBUG:: idx: " << idx << "\n";
               std::cout << "ERROR:: failed to find atom in molecule: chain-id :"
                         << std::string(atom->residue->GetChainID()) <<  ": res_no "
                         << atom->residue->seqNum << " inscode :"
                         << std::string(atom->GetInsCode()) << ": name :"
                         << std::string(atom->name) << ": altloc :"
                         << std::string(atom->altLoc) << ":" << std::endl;
            }
         }

         if (change_altconf_occs_flag) {
            if (idx >= 0) {
               n_atom++;
               mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
               float atom_occ = atom->occupancy;
               // if this is a shelx molecule, then we don't change
               // occupancies this way.  We do it by changing the FVAR
               if (is_from_shelx_ins_flag) {
                  atom_occ = mol_atom->occupancy;

                  // OK, one more go.  We have an occupancy of 31 or -31
                  // say.  Now, the alt conf atoms has been immmediately
                  // added with the old occupancy for the actual FVAR number
                  // - this happens before we get to twiddle the occupancy
                  // slider.  So here we have to find out the index of the
                  // replaced atom and set it's fvar to whatever the slider
                  // value had been set to.

                  int fvar_number = coot::ShelxIns::shelx_occ_to_fvar(atom_occ);
                  if (fvar_number > 1) {
                     //                std::cout << "DEBUG:: replace_coords: setting fvar number "
                     //                          <<  fvar_number << " (generated from occ " << atom_occ << ") to "
                     //                          << graphics_info_t::add_alt_conf_new_atoms_occupancy << std::endl;
                     shelxins.set_fvar(fvar_number, add_alt_conf_new_atoms_occupancy);
                  }

                  if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag))
                     mol_atom->SetCoordinates(atom->x,
                                              atom->y,
                                              atom->z,
                                              atom_occ,
                                              mol_atom->tempFactor);
               } else {
                  if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag))
                     mol_atom->SetCoordinates(atom->x,
                                              atom->y,
                                              atom->z,
                                              atom_occ,
                                              mol_atom->tempFactor);
               }

               // similarly we adjust occupancy if this is not a shelx molecule
               if (! is_from_shelx_ins_flag) {
                  adjust_occupancy_other_residue_atoms(mol_atom, mol_atom->residue, 0);
               }
               // std::cout << atom << " coords replace " << idx << " " << mol_atom << std::endl;
            } else {
               std::cout << "ERROR:: bad atom index in replace_coords replacing atom: "
                         << atom << std::endl;
            }
         } else {

            // "don't change alt confs" mode

            if (idx != -1 ) {
               mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
               if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag)) {
                  if (debug) {
                     coot::Cartesian old_pos(mol_atom->x, mol_atom->y, mol_atom->z);
                     coot::Cartesian new_pos(atom->x, atom->y, atom->z);
                     double d = (new_pos - old_pos).amplitude();
                     std::cout << "    changing coords for atom with idx " << idx << " "
                               << coot::atom_spec_t(mol_atom) << std::endl;
                     std::cout << "   " << old_pos << " " << new_pos << " moved-by " << d << std::endl;
                  }
                  mol_atom->SetCoordinates(atom->x,
                                           atom->y,
                                           atom->z,
                                           mol_atom->occupancy,
                                           mol_atom->tempFactor);
                  n_atom++;
               }
            } else {
               std::cout << "WARNING:: bad atom idx -1" << std::endl;
            }
         }
      }
   }
   std::cout << "INFO:: replace_coords: " << n_atom << " atoms updated." << std::endl;
   // have_unsaved_changes_flag = 1;
   save_info.new_modification();

   if (show_symmetry) {  // internal
      update_symmetry();
   }

   // make_bonds_type_checked(__FUNCTION__);

}


int coot::molecule_t::flipPeptide(const coot::residue_spec_t &rs, const std::string &alt_conf) {

   int result = coot::pepflip(atom_sel.mol, rs.chain_id, rs.res_no, rs.ins_code, alt_conf);
   return result;

}

// public - because currently making bonds is not done on molecule construction
void
coot::molecule_t::make_bonds(coot::protein_geometry *geom, coot::rotamer_probability_tables *rot_prob_tables_p) {

   bonds_box_type = coot::COLOUR_BY_CHAIN_BONDS;
   make_bonds_type_checked(geom, rot_prob_tables_p, __FUNCTION__);

   std::cout << "debug:: in molecule_t::make_bonds() " << bonds_box.n_bonds() << " bonds " << bonds_box.n_atoms() << " atoms "
             << std::endl;
}



// private
void
coot::molecule_t::makebonds(coot::protein_geometry *geom, coot::rotamer_probability_tables *rotamer_tables_p, std::set<int> &no_bonds_to_these_atoms) {

   bool force_rebond = true;
   bool do_rotamer_markup = true; // pass this
   make_colour_by_chain_bonds(geom, no_bonds_to_these_atoms, true, false, do_rotamer_markup, rotamer_tables_p, force_rebond);

}


void
coot::molecule_t::make_colour_by_chain_bonds(coot::protein_geometry *geom,
                                             const std::set<int> &no_bonds_to_these_atoms,
                                             bool change_c_only_flag,
                                             bool goodsell_mode,
                                             bool do_rota_markup,
                                             coot::rotamer_probability_tables *tables_p,
                                             bool force_rebonding) {

   // We don't want to rebond if we don't have to (i.e the mode requested is the current mode)
   // so check the previous value of bonds_box_type so that we can know if it can be skipped.

   bool draw_hydrogens_flag = true; // pass this
   bool draw_missing_loops_flag = true; // pass this

   Bond_lines_container bonds(geom, no_bonds_to_these_atoms, draw_hydrogens_flag);

   std::cout << "debug:: in make_colour_by_chain_bonds() adding tables_p " << tables_p << std::endl;

   bonds.add_rotamer_tables(tables_p);
   bonds.do_colour_by_chain_bonds(atom_sel, false, imol_no, draw_hydrogens_flag,
                                  draw_missing_loops_flag,
                                  change_c_only_flag, goodsell_mode, do_rota_markup);

   // std::cout << "------------------- calling make_graphical_bonds_no_thinning() " << std::endl;

   bonds_box = bonds.make_graphical_bonds_no_thinning(); // make_graphical_bonds() is pretty
                                                         // stupid when it comes to thining.

   // bonds_box = bonds.make_graphical_bonds(); // make_graphical_bonds() is pretty
                                                // stupid when it comes to thining.

   // testing previous values of bonds_box_type
   if (bonds_box_type != coot::COLOUR_BY_CHAIN_BONDS)
      force_rebonding = true;

   if (goodsell_mode)
      if (bonds_box_type != coot::COLOUR_BY_CHAIN_GOODSELL)
         force_rebonding = true;

   bonds_box_type = coot::COLOUR_BY_CHAIN_BONDS;

   if (goodsell_mode)
      bonds_box_type = coot::COLOUR_BY_CHAIN_GOODSELL;

   // 20221011-PE Hmm... is this needed in this API? I don't think so
   //
   // if (force_rebonding)
   //    make_glsl_bonds_type_checked(__FUNCTION__);

}

void
coot::molecule_t::make_ca_bonds() {

}


void
coot::molecule_t::make_bonds_type_checked(coot::protein_geometry *geom_p, coot::rotamer_probability_tables *rotamer_probability_tables_p, const char *caller) {

   bool draw_missing_loops_flag = false; // pass this
   bool rotate_colour_map_on_read_pdb_c_only_flag = true; // pass this or make class data item

   bool debug = false;

   // Note caller can be 0 (e.g. with clang) - so be aware of that when debugging.

   std::string caller_s("NULL");
   if (caller) caller_s = std::string(caller);

   bool is_intermediate_atoms_molecule = false; // 20221005-PE IMPORT-HACK
   if (debug)
      std::cout << "debug:: plain make_bonds_type_checked() --------start--------- called by "
                << caller_s << "() with is_intermediate_atoms_molecule: " << is_intermediate_atoms_molecule
                << std::endl;
   if (debug)
      std::cout << "--------- make_bonds_type_checked() called with bonds_box_type "
                << bonds_box_type << " vs "
                << "NORMAL_BONDS " << coot::NORMAL_BONDS << " "
                << "BONDS_NO_HYDROGENS " << coot::BONDS_NO_HYDROGENS << " "
                << "COLOUR_BY_CHAIN_BONDS " << coot::COLOUR_BY_CHAIN_BONDS << " "
                << "COLOUR_BY_MOLECULE_BONDS " << coot::COLOUR_BY_MOLECULE_BONDS << " "
                << "CA_BONDS " << coot::CA_BONDS << " "
                << "CA_BONDS_PLUS_LIGANDS " << coot::CA_BONDS_PLUS_LIGANDS << " "
                << "COLOUR_BY_USER_DEFINED_COLOURS___BONDS " << coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS << " "
                << std::endl;

   // Delete this in due course
   // graphics_info_t g; // urgh!  (But the best solution?)

   bool force_rebonding = true; // if we get here, this must be true (?)

   // coot::protein_geometry *geom_p = g.Geom_p();

   std::set<int> dummy;

   if (bonds_box_type == coot::NORMAL_BONDS)
      makebonds(geom_p, nullptr, dummy);

   if (bonds_box_type == coot::BONDS_NO_HYDROGENS)
      makebonds(geom_p, nullptr, dummy);
   if (bonds_box_type == coot::CA_BONDS)
      make_ca_bonds();
   if (bonds_box_type == coot::COLOUR_BY_CHAIN_BONDS || bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL) {
      // Baah, we have to use the static in graphics_info_t here as it
      // is not a per-molecule property.
      std::set<int> s;
      bool goodsell_mode = false;
      if (bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL)
         goodsell_mode = true;
      bool do_rota_markup = true;

      make_colour_by_chain_bonds(geom_p, s, rotate_colour_map_on_read_pdb_c_only_flag, goodsell_mode, do_rota_markup, rotamer_probability_tables_p, force_rebonding);
   }

#if 0 // not implemenented yet
   if (bonds_box_type == coot::COLOUR_BY_MOLECULE_BONDS)
      make_colour_by_molecule_bonds(force_rebonding);
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS)
      make_ca_plus_ligands_bonds(g.Geom_p());
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS)
      make_ca_plus_ligands_and_sidechains_bonds(g.Geom_p());
   if (bonds_box_type == coot::BONDS_NO_WATERS)
      bonds_no_waters_representation();
   if (bonds_box_type == coot::BONDS_SEC_STRUCT_COLOUR)
      bonds_sec_struct_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR)
      ca_plus_ligands_sec_struct_representation(g.Geom_p());
   if (bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS)
      ca_plus_ligands_rainbow_representation(g.Geom_p());
   if (bonds_box_type == coot::COLOUR_BY_OCCUPANCY_BONDS)
      occupancy_representation();
   if (bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS)
      b_factor_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR)
      b_factor_representation_as_cas();
   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS)
      user_defined_colours_representation(g.Geom_p(), true, g.draw_missing_loops_flag); // hack,
                                                             // because we need to remeber somehow
                                                             // if this was called with all-atom or CA-only.
                                                             // See c-interface.cc
                                                             // graphics_to_user_defined_atom_colours_representation()
                                                             // Perhaps we need two functions
                                                             // user_defined_colours_representation_all()
                                                             // user_defined_colours_representation_Calpha() [+ ligands]
   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS)
      user_defined_colours_representation(geom_p, false, draw_missing_loops_flag); // hack,

#endif


#if 0 // 20221005-PE not sure what these are
   // all these will need to be changed or removed
   update_additional_representations(glci, g.Geom_p());
   update_fixed_atom_positions();
   update_ghosts();
   update_extra_restraints_representation();
#endif

   if (debug) {
      std::cout << "debug:: -------------- make_bonds_type_checked() done " << std::endl;
   }
}

std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> >
coot::molecule_t::ramachandran_validation() const {

   auto get_HA_unit_vector = [] (mmdb::Residue *r) {
      bool status = false;
      coot::Cartesian dir;
      mmdb::Atom *CA = r->GetAtom(" CA ");
      mmdb::Atom *C  = r->GetAtom(" C  ");
      mmdb::Atom *N  = r->GetAtom(" N  ");
      mmdb::Atom *CB = r->GetAtom(" CB ");

      if (CA && C && N && CB) {
         coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
         coot::Cartesian  c_pos( C->x,  C->y,  C->z);
         coot::Cartesian  n_pos( N->x,  N->y,  N->z);
         coot::Cartesian cb_pos(CB->x, CB->y, CB->z);
         coot::Cartesian dir_1 = ca_pos - c_pos;
         coot::Cartesian dir_2 = ca_pos - n_pos;
         coot::Cartesian dir_3 = ca_pos - cb_pos;
         coot::Cartesian r = dir_1 + dir_2 + dir_3;
         dir = r.unit();
         status = true;
      } else {
         if (CA && C && N) {
            coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
            coot::Cartesian  c_pos( C->x,  C->y,  C->z);
            coot::Cartesian  n_pos( N->x,  N->y,  N->z);
            coot::Cartesian dir_1 = ca_pos - c_pos;
            coot::Cartesian dir_2 = ca_pos - n_pos;
            coot::Cartesian r = dir_1 + dir_2;
            dir = r.unit();
            status = true;
         }
      }
      return std::make_pair(status, dir);
   };

   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> > v;

   float rama_ball_pos_offset_scale = 0.6;

   rama_plot::phi_psis_for_model_t ppm(atom_sel.mol);
   // This: std::map<coot::residue_spec_t, phi_psi_t> ppm.phi_psi  is now filled

   std::map<residue_spec_t, rama_plot::phi_psi_t>::const_iterator it;
   for (it=ppm.phi_psi.begin(); it!=ppm.phi_psi.end(); ++it) {
      const auto &phi_psi(it->second);
      mmdb::Residue *rp = phi_psi.residue_prev;
      mmdb::Residue *rt = phi_psi.residue_this;
      mmdb::Residue *rn = phi_psi.residue_next;
      if (rp && rt && rn) {
         mmdb::Atom *at = rt->GetAtom(" CA "); // 20221006-PE alt-confs another day
         if (at) {
            coot::Cartesian pos(at->x, at->y, at->z);
            std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(rt);
            coot::Cartesian offset(0,0,rama_ball_pos_offset_scale);
            if (hav.first) offset = hav.second * rama_ball_pos_offset_scale;
            coot::util::phi_psi_t cupp(rp, rt, rn);
            auto p = std::make_pair(pos + offset, cupp);
            v.push_back(p);
         }
      }
   }

   return v;
}


// returns either the specified atom or null if not found
mmdb::Atom *
coot::molecule_t::get_atom(const coot::atom_spec_t &atom_spec) const {

   mmdb::Atom *at = coot::util::get_atom(atom_spec, atom_sel.mol);
   return at;
}


#include "utils/dodec.hh"

// returns either the specified residue or null if not found
mmdb::Residue *
coot::molecule_t::get_residue(const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *r = coot::util::get_residue(residue_spec, atom_sel.mol);
   return r;

}

coot::simple_mesh_t
coot::molecule_t::get_rotamer_dodecs(coot::protein_geometry *geom_p,
                                     coot::rotamer_probability_tables *rpt) {

   // THis function is an API version of:
   //
   // void
   // Mesh::make_graphical_bonds_rotamer_dodecs(const graphical_bonds_container &gbc,
   // const glm::vec3 &screen_up_dir)

   simple_mesh_t m;

   // use bonds_box

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
      return glm::vec3(c.x(), c.y(), c.z()); };

   auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
      return glm::vec3(c.x(), c.y(), c.z()); };

   auto clipper_to_cartesian = [] (const clipper::Coord_orth &c) {
      return Cartesian(c.x(), c.y(), c.z()); };

   auto colour_holder_to_glm = [] (const coot::colour_holder &ch) {
                                  return glm::vec4(ch.red, ch.green, ch.blue, 1.0f);
                               };

   auto get_HA_unit_vector = [] (mmdb::Residue *r) {
      bool status = false;
      coot::Cartesian dir;
      mmdb::Atom *CA = r->GetAtom(" CA ");
      mmdb::Atom *C  = r->GetAtom(" C  ");
      mmdb::Atom *N  = r->GetAtom(" N  ");
      mmdb::Atom *CB = r->GetAtom(" CB ");

      if (CA && C && N && CB) {
         coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
         coot::Cartesian  c_pos( C->x,  C->y,  C->z);
         coot::Cartesian  n_pos( N->x,  N->y,  N->z);
         coot::Cartesian cb_pos(CB->x, CB->y, CB->z);
         coot::Cartesian dir_1 = ca_pos - c_pos;
         coot::Cartesian dir_2 = ca_pos - n_pos;
         coot::Cartesian dir_3 = ca_pos - cb_pos;
         coot::Cartesian r = dir_1 + dir_2 + dir_3;
         dir = r.unit();
         status = true;
      } else {
         if (CA && C && N) {
            coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
            coot::Cartesian  c_pos( C->x,  C->y,  C->z);
            coot::Cartesian  n_pos( N->x,  N->y,  N->z);
            coot::Cartesian dir_1 = ca_pos - c_pos;
            coot::Cartesian dir_2 = ca_pos - n_pos;
            coot::Cartesian r = dir_1 + dir_2;
            dir = r.unit();
            status = true;
         }
      }
      return std::make_pair(status, dir);
   };

   std::set<int> dummy;
   bool do_rota_markup = true;
   make_colour_by_chain_bonds(geom_p, dummy, true, false, do_rota_markup, rpt, true);

   if (bonds_box.n_rotamer_markups > 0) {

      auto &vertices = m.vertices;
      auto &triangles = m.triangles;

      glm::vec4 col(0.6, 0.2, 0.8, 1.0); // starting colour
      dodec d;
      std::vector<clipper::Coord_orth> coords = d.coords();
      std::vector<glm::vec3> dodec_postions(coords.size());
      for (unsigned int i=0; i<coords.size(); i++)
         dodec_postions[i] = clipper_to_glm(coords[i]);

      std::vector<coot::api::vnc_vertex> dodec_vertices;
      std::vector<g_triangle> dodec_triangles;
      dodec_triangles.reserve(36);

      for (unsigned int iface=0; iface<12; iface++) {

         std::vector<coot::api::vnc_vertex> face_verts;
         std::vector<g_triangle> face_triangles;
         face_triangles.reserve(3);

         std::vector<unsigned int> indices_for_face = d.face(iface);
         glm::vec3 ns(0,0,0);
         for (unsigned int j=0; j<5; j++)
            ns += dodec_postions[indices_for_face[j]];
         glm::vec3 normal = glm::normalize(ns);

         for (unsigned int j=0; j<5; j++) {
            glm::vec3 &pos = dodec_postions[indices_for_face[j]];
            coot::api::vnc_vertex v(0.5f * pos, normal, col);
            face_verts.push_back(v);
         }

         face_triangles.push_back(g_triangle(0,1,2));
         face_triangles.push_back(g_triangle(0,2,3));
         face_triangles.push_back(g_triangle(0,3,4));

         unsigned int idx_base = dodec_vertices.size();
         unsigned int idx_tri_base = dodec_triangles.size();
         dodec_vertices.insert(dodec_vertices.end(), face_verts.begin(), face_verts.end());
         dodec_triangles.insert(dodec_triangles.end(), face_triangles.begin(), face_triangles.end());
         for (unsigned int jj=idx_tri_base; jj<dodec_triangles.size(); jj++)
            dodec_triangles[jj].rebase(idx_base);
      }

      // now there is a dodec at the origin, dodec_vertices and dodec_triangle

      // let's make copies of that and move them around to the residues

      double rama_ball_pos_offset_scale = 1.2; // may need tweaking
      for (int i=0; i<bonds_box.n_rotamer_markups; i++) {
         const rotamer_markup_container_t &rm = bonds_box.rotamer_markups[i];
         const residue_spec_t &residue_spec = rm.spec;
         mmdb::Residue *residue_p = get_residue(residue_spec);
         Cartesian offset(0,0,rama_ball_pos_offset_scale);
         if (residue_p) {
            std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(residue_p);
            if (hav.first) offset = hav.second; // * rama_ball_pos_offset_scale;
         }

         glm::vec3 atom_pos = cartesian_to_glm(rm.pos) + cartesian_to_glm(offset);
         auto this_dodec_colour = colour_holder_to_glm(rm.col);

         std::vector<coot::api::vnc_vertex> this_dodec_vertices = dodec_vertices; // at the origin to start

         // now move it.
         for (unsigned int j=0; j<this_dodec_vertices.size(); j++) {
            auto &vertex = this_dodec_vertices[j];
            vertex.pos  += atom_pos;
            vertex.color = this_dodec_colour;
            if (false)
               std::cout << "atom_pos " << glm::to_string(vertex.pos)
                         << " rama_markup_col " << rm.col
                         << " color " << glm::to_string(vertex.color) << std::endl;
         }
         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         vertices.insert(vertices.end(), this_dodec_vertices.begin(), this_dodec_vertices.end());
         triangles.insert(triangles.end(), dodec_triangles.begin(), dodec_triangles.end());
         for (unsigned int jj=idx_tri_base; jj<triangles.size(); jj++)
            triangles[jj].rebase(idx_base);
      }
   }
   return m;
}

#include "ligand/backrub-rotamer.hh"

std::pair<bool,float>
coot::molecule_t::backrub_rotamer(const std::string &chain_id, int res_no,
                                  const std::string &ins_code, const std::string &alt_conf,
                                  const clipper::Xmap<float> &xmap_in,
                                  const coot::protein_geometry &pg) {

   bool status = false;
   float score = -1;
   bool refinement_move_atoms_with_zero_occupancy_flag = true; // pass this?

   std::cout << "debug:: molecule_t::backrub_rotamer() starts " << chain_id << " " << res_no << std::endl;

   residue_spec_t res_spec(chain_id, res_no, ins_code);
   mmdb::Residue *res = get_residue(res_spec);
   if (! res) {
      std::cout << "   WARNING:: residue in molecule :" << chain_id << ": "
                << res_no << " inscode :" << ins_code << ": altconf :"
                << alt_conf << ":" << std::endl;
   } else {
      std::string monomer_type = res->GetResName();
      std::pair<short int, coot::dictionary_residue_restraints_t> p =
         pg.get_monomer_restraints(monomer_type, imol_no);
      coot::dictionary_residue_restraints_t restraints = p.second;

      if (p.first) {
         try {

            mmdb::Residue *prev_res = coot::util::previous_residue(res);
            mmdb::Residue *next_res = coot::util::next_residue(res);
            mmdb::Manager *mol = atom_sel.mol;
            coot::backrub br(chain_id, res, prev_res, next_res, alt_conf, mol,
                             &xmap_in); // use a pointer for the map
            std::pair<coot::minimol::molecule,float> m = br.search(restraints);
            std::vector<coot::atom_spec_t> baddie_waters = br.waters_for_deletion();
            score = m.second;
            status = true;
            atom_selection_container_t fragment_asc = make_asc(m.first.pcmmdbmanager());
            replace_coords(fragment_asc, 0, refinement_move_atoms_with_zero_occupancy_flag);
            if (baddie_waters.size())
               delete_atoms(baddie_waters);
         }
         catch (const std::runtime_error &rte) {
            std::cout << "WARNING:: thrown " << rte.what() << std::endl;
         }
      } else {
         std::cout << " No restraints found for " << monomer_type << std::endl;
      }
   }
   return std::pair<bool,float> (status, score);
}


int
coot::molecule_t::auto_fit_rotamer(const std::string &chain_id, int res_no, const std::string &ins_code,
                                   const std::string &alt_conf,
                                   const clipper::Xmap<float> &xmap, const coot::protein_geometry &geom) {

   int status = 0;

   // let's just use backrub rotamer for now.
   std::pair<bool, float> r = backrub_rotamer(chain_id, res_no, ins_code, alt_conf, xmap, geom);
   status = r.first;

   return status;
}


int
coot::molecule_t::delete_atoms(const std::vector<coot::atom_spec_t> &atom_specs) {

   short int was_deleted = 0;
   int n_deleleted_atoms = 0;

   if (atom_sel.n_selected_atoms > 0) {
      if (atom_specs.size() > 0)
         make_backup();
      for (unsigned int i=0; i<atom_specs.size(); i++) {
         int SelHnd = atom_sel.mol->NewSelection();
         // how about a function that calls this:
         // select_atomspec_atoms(atom_sel.mol, atom_specs[i])
         // or  a member function of an atom_spec_t:
         //    atom_specs[i].select_atoms(mol)
         //
         mmdb::PAtom *atoms = NULL;
         int n_atoms;
         atom_sel.mol->SelectAtoms(SelHnd, 0, atom_specs[i].chain_id.c_str(),
                                   atom_specs[i].res_no, atom_specs[i].ins_code.c_str(),
                                   atom_specs[i].res_no, atom_specs[i].ins_code.c_str(),
                                   "*",
                                   atom_specs[i].atom_name.c_str(),
                                   "*",
                                   atom_specs[i].alt_conf.c_str()
                                   );
         atom_sel.mol->GetSelIndex(SelHnd, atoms, n_atoms);
         if (n_atoms) {
            delete atoms[0];
            atoms[0] = NULL;
            n_deleleted_atoms++;
            was_deleted = 1;
         }
         atom_sel.mol->DeleteSelection(SelHnd);
      }
   }

   // potentially
   if (was_deleted) {
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      // make_bonds_type_checked(__FUNCTION__);
      // have_unsaved_changes_flag = 1;
      save_info.new_modification();
      // unlikely to be necessary:
      trim_atom_label_table();

      // update_symmetry();
   }

   return n_deleleted_atoms;
}

void
coot::molecule_t::update_symmetry() {

   // 20221013-PE no symmetry yet
}

void
coot::molecule_t::trim_atom_label_table() {

   // 20221013-PE no atom labels yet
}

void
coot::molecule_t::delete_ghost_selections() {

   // 20221013-PE no ghosts at the moment
}

int
coot::molecule_t::delete_atom(coot::atom_spec_t &atom_spec) {

   int was_deleted = 0;
   mmdb::Chain *chain;
   mmdb::Residue *residue_of_deleted_atom = NULL;

   std::string chain_id = atom_spec.chain_id;
   int resno            = atom_spec.res_no;
   std::string ins_code = atom_spec.ins_code;
   std::string atname   = atom_spec.atom_name;
   std::string altconf  = atom_spec.alt_conf;

   // run over chains of the existing mol
   int nchains = atom_sel.mol->GetNumberOfChains(1);
   for (int ichain=0; ichain<nchains; ichain++) {

      chain = atom_sel.mol->GetChain(1,ichain);
      std::string mol_chain_id(chain->GetChainID());

      // Note, if in the PDB file, the chain id is not set to
      // something, A, B, C etc, then the chain id from mmdb is ""
      // (not " "!)

//       std::cout << "debug:: delete_atom comparing chain_ids :"
//                 << chain_id << ": vs :" << mol_chain_id << ":"
//                 << std::endl;

      if (chain_id == mol_chain_id) {

         int nres = chain->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::PResidue res = chain->GetResidue(ires);
            std::string ins_code_local = res->GetInsCode();

            if (res) {
//                std::cout << "debug:: delete_atom: comparing residue "
//                          << res->GetSeqNum() << " :" << ins_code_local
//                          << ":" << " to " << resno << " :" << ins_code << ":"
//                          << std::endl;
               if (res->GetSeqNum() == resno) {
                  if (ins_code_local == ins_code) {

                     // so we have a matching residue:
                     // std::cout << "debug:: delete_atom: we have a matching residue "
                     // << resno << " :" << ins_code << ":" << std::endl;

                     mmdb::PPAtom residue_atoms;
                     int nResidueAtoms;
                     std::string mol_atom_name;
                     res->GetAtomTable(residue_atoms, nResidueAtoms);
                     for (int iat=0; iat<nResidueAtoms; iat++) {

                        mol_atom_name = residue_atoms[iat]->name;
                        if (atname == mol_atom_name) {

                           if (std::string(residue_atoms[iat]->altLoc) == altconf) {

                              make_backup();
                              atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                              delete_ghost_selections();
                              res->DeleteAtom(iat);
                              was_deleted = 1;
                              residue_of_deleted_atom = res;
                              break;
                           }
                        }
                     }
                  }
               }
            }
            if (was_deleted)
               break;
         }
      }
      if (was_deleted)
         break;
   }

   // potentially
   if (was_deleted) {
      atom_sel.mol->FinishStructEdit();


      // Now reset/recalculate the occupancy and the altLoc of the
      // remaining atoms in the residue with the same atom name.
      //
      mmdb::PPAtom atoms = NULL;
      int n_atoms;
      mmdb::Atom *at = 0;
      int n_matching_name = 0;
      residue_of_deleted_atom->GetAtomTable(atoms, n_atoms);
      for (int iat=0; iat<n_atoms; iat++) {
         std::string res_atom_name = atoms[iat]->name;
         if (res_atom_name == atname) {
            at = atoms[iat];
            n_matching_name++;
         }
      }
      if (n_matching_name == 1) { // one atom of this name left in the residue, so
                                    // remove its altconf string
         if (at) {
            strncpy(at->altLoc, "", 2);
            // set the occupancy to 1.0 of the remaining atom if it was not zero.
            if (at->occupancy > 0.009)
               at->occupancy = 1.0;
         }
      }

      atom_sel = make_asc(atom_sel.mol);
      // make_bonds_type_checked(__FUNCTION__);
      // have_unsaved_changes_flag = 1;
      save_info.new_modification();
      // unlikely to be necessary:
      trim_atom_label_table();
      update_symmetry();
   }
   return was_deleted;
}


int
coot::molecule_t::delete_residue(coot::residue_spec_t &residue_spec) {

   int status = 0;

   return status;
}
