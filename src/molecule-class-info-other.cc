/* src/molecule-class-info-other.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by The University of Oxford
 * Copyright 2007, 2008, 2009 The University of Oxford
 * Copyright 2013, 2014, 2015 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include <stdlib.h>

#if !defined WINDOWS_MINGW && !defined _MSC_VER
#  include <glob.h>
#  include <unistd.h>
#else
#   ifdef _MSC_VER
#     include <windows.h>
#   else
#     include <unistd.h>
#     include <glob.h>
#   endif
#endif

#include <string.h>  // strncpy
#include <sys/types.h>  // for stating
#include <sys/stat.h>

#include <iostream>
#include <vector>

#include "compat/coot-sysdep.h"

#include "clipper/core/xmap.h"
#include "clipper/cns/cns_hkl_io.h"
#include "clipper/minimol/minimol_io.h"

#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
#undef V_UNKNOWN
#define V_UNKNOWNA V_UNKNOWN
#endif

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"

#include "graphics-info.h"
#include "xmap-utils.h"
#include "coot-utils/xmap-stats.hh"

#include "molecule-class-info.h"

#ifdef USE_DUNBRACK_ROTAMERS
#include "ligand/dunbrack.hh"
#else
#include "ligand/richardson-rotamer.hh"
#endif

#include "geometry/mol-utils.hh"
#include "ligand/ligand.hh"
#include "utils/coot-utils.hh"
#include "coot-utils/lsq-improve.hh"
#include "coot-utils/coot-trim.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-coord-utils.hh" // check_dictionary_for_residue
#include "coot-utils/coot-map-heavy.hh"   // situation's heavy... [simplex]
#include "ideal/pepflip.hh"
#include "ligand/backrub-rotamer.hh"
#include "ligand/torsion-general.hh"

#include "api/coot-molecule.hh"  // pick up RESIDUE_NUMBER_UNSET (it used to be in molecule-class-info.h)
                                 // I don't think that this is a good organization

#include "coot-nomenclature.hh"

// #include "GL/glu.h"
// #include "GL/glut.h"

#include "rotamer-search-modes.hh"


// ---------------------------------------------------------------------------------------
//                           molecule class
// ---------------------------------------------------------------------------------------
void
molecule_class_info_t::debug_selection() const {

   // sigh - debugging
   //

   int SelHnd = atom_sel.mol->NewSelection();
   mmdb::PPAtom atom = NULL;
   int n_atoms;

   atom_sel.mol->SelectAtoms(SelHnd,
                    0,
                    "A",
                    888, "*",
                    890, "*",
                    "*", // rnames
                    "*", // anames
                    "*", // elements
                    "*"  // altLocs
                    );

   atom_sel.mol->GetSelIndex(SelHnd, atom, n_atoms);
   if (n_atoms == 0) {
      std::cout << "debug_selection: no atoms selected" << std::endl;
   } else {
      std::cout << "debug_selection: selected atoms" << std::endl;
      for(int i=0; i<n_atoms; i++) {
         std::cout << atom[i] << std::endl;
      }
      std::cout << "----------- " << std::endl;
   }
}

// make this a bool
short int
molecule_class_info_t::molecule_is_all_c_alphas() const {

   short int is_ca = 1;

   int n_atoms = atom_sel.n_selected_atoms;
   if (n_atoms == 0) {
      is_ca = 0;
   } else {
      for (int i=0; i<n_atoms; i++) {
         std::string name_string(atom_sel.atom_selection[i]->name);
         if (name_string != " CA " ) {
            is_ca = 0;
            break;
         }
      }
   }
   return is_ca;
}

void
molecule_class_info_t::bond_representation(const coot::protein_geometry *geom_p,
                                           bool force_rebonding) {

   bool do_rebond = true;
   if (draw_hydrogens_flag && bonds_box_type == coot::NORMAL_BONDS)
      do_rebond = false;
   if (!draw_hydrogens_flag && bonds_box_type == coot::BONDS_NO_HYDROGENS)
      do_rebond = false;

   if  (force_rebonding)
      do_rebond = true;

   if (bonds_box_type != coot::NORMAL_BONDS)
      do_rebond = true;

   if (do_rebond) {
      std::set<int> dummy;
      makebonds(geom_p, dummy);
   }
}

void
molecule_class_info_t::ca_representation(bool force_rebonding) {

   if (bonds_box_type != coot::CA_BONDS || force_rebonding) {

      // std::cout << "calling make_ca_bonds() in molecule_class_info_t" << std::endl;
      bonds_box.clear_up();
      make_ca_bonds(2.4, 4.7);
      bonds_box_type = coot::CA_BONDS;
   }
}

void
molecule_class_info_t::ca_plus_ligands_representation(coot::protein_geometry *geom, bool force_rebonding) {

   if (bonds_box_type != coot::CA_BONDS_PLUS_LIGANDS || force_rebonding) {
      bonds_box.clear_up();
      make_ca_plus_ligands_bonds(geom);
      bonds_box_type = coot::CA_BONDS_PLUS_LIGANDS;
   }
}

void
molecule_class_info_t::ca_plus_ligands_and_sidechains_representation(coot::protein_geometry *geom) {

   bool force_rebonding = true; // for now
   if (force_rebonding) {
      bonds_box.clear_up();
      make_ca_plus_ligands_and_sidechains_bonds(geom);
      bonds_box_type = coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS;
   }
}

void
molecule_class_info_t::bonds_no_waters_representation() {

   if (bonds_box_type != coot::BONDS_NO_WATERS) {
      bonds_box.clear_up();
      Bond_lines_container bonds;
      bonds.do_normal_bonds_no_water(atom_sel, imol_no, 0.01, 1.9);
      bonds_box = bonds.make_graphical_bonds();
      bonds_box_type = coot::BONDS_NO_WATERS;
      make_glsl_bonds_type_checked(__FUNCTION__);
   }
}

void
molecule_class_info_t::bonds_sec_struct_representation() {

   if (bonds_box_type != coot::BONDS_SEC_STRUCT_COLOUR) {
      std::set<int> no_bonds_to_these_atom_indices;
      // I am guessing that this is the constructor that I wanted?
      // It looks like it was not getting the right one before.
      Bond_lines_container bonds(graphics_info_t::Geom_p(), no_bonds_to_these_atom_indices, draw_hydrogens_flag);
      bonds.do_colour_sec_struct_bonds(atom_sel, imol_no, 0.01, 1.9);
      bonds_box = bonds.make_graphical_bonds_no_thinning();
      bonds_box_type = coot::BONDS_SEC_STRUCT_COLOUR;
      make_glsl_bonds_type_checked(__FUNCTION__);
   }
}


void
molecule_class_info_t::ca_plus_ligands_sec_struct_representation(coot::protein_geometry *pg) {

   Bond_lines_container bonds;
   bonds.do_Ca_plus_ligands_colour_sec_struct_bonds(atom_sel, imol_no, pg, 2.4, 4.7,
                                                    draw_hydrogens_flag,
                                                    graphics_info_t::draw_missing_loops_flag);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR;
   make_glsl_bonds_type_checked(__FUNCTION__);

}

void
molecule_class_info_t::ca_plus_ligands_rainbow_representation(coot::protein_geometry *pg) {

    //
    Bond_lines_container bonds;
    bonds.do_Ca_plus_ligands_bonds(atom_sel, imol_no, pg,
                                   2.4, 4.7,
                                   graphics_info_t::draw_missing_loops_flag,
                                   coot::COLOUR_BY_RAINBOW,
                                   draw_hydrogens_flag); // not COLOUR_BY_RAINBOW_BONDS
    bonds_box = bonds.make_graphical_bonds_no_thinning();
    bonds_box_type = coot::COLOUR_BY_RAINBOW_BONDS;
    make_glsl_bonds_type_checked(__FUNCTION__);

}

void
molecule_class_info_t::b_factor_representation() {

   Bond_lines_container::bond_representation_type bond_type = Bond_lines_container::COLOUR_BY_B_FACTOR;

   Bond_lines_container bonds(atom_sel, imol_no, graphics_info_t::Geom_p(), bond_type);
   bonds_box = bonds.make_graphical_bonds_no_thinning();
   bonds_box_type = coot::COLOUR_BY_B_FACTOR_BONDS;
   make_glsl_bonds_type_checked(__FUNCTION__);
}

void
molecule_class_info_t::b_factor_representation_as_cas() {

   // std::cout << "************************************************ b_factor_representation_as_cas() " << std::endl;

   Bond_lines_container::bond_representation_type bond_type = Bond_lines_container::COLOUR_BY_B_FACTOR;
   Bond_lines_container bonds;
   bonds.do_Ca_plus_ligands_bonds(atom_sel, imol_no, NULL, 2.4, 4.7,
                                  graphics_info_t::draw_missing_loops_flag,
                                  Bond_lines_container::COLOUR_BY_B_FACTOR,
                                  draw_hydrogens_flag); // pass a dictionary
   bonds_box = bonds.make_graphical_bonds_no_thinning();
   bonds_box_type = coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR;
   make_glsl_bonds_type_checked(__FUNCTION__);
}

void
molecule_class_info_t::occupancy_representation() {

   Bond_lines_container::bond_representation_type bond_type =
      Bond_lines_container::COLOUR_BY_OCCUPANCY;

   // 20241130-PE Constructor I - it was using Constructor A until now.
   Bond_lines_container bonds(atom_sel, imol_no, graphics_info_t::Geom_p(), bond_type);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::COLOUR_BY_OCCUPANCY_BONDS;
   make_glsl_bonds_type_checked(__FUNCTION__);
}


int
molecule_class_info_t::set_atom_attribute(std::string chain_id, int resno, std::string ins_code,
                                          std::string atom_name, std::string alt_conf,
                                          std::string attribute_name, float val) {

   int istate = 0;
   if (atom_sel.n_selected_atoms > 0) {

      int SelectionHandle = atom_sel.mol->NewSelection();
      atom_sel.mol->SelectAtoms(SelectionHandle, 0,
                                (char *) chain_id.c_str(),
                                resno, (char *) ins_code.c_str(),
                                resno, (char *) ins_code.c_str(),
                                "*",  // res type
                                (char *) atom_name.c_str(),
                                "*",  // any element
                                (char *) alt_conf.c_str(), mmdb::SKEY_NEW);
      int nSelAtoms;
      mmdb::PPAtom SelAtoms = NULL;
      atom_sel.mol->GetSelIndex(SelectionHandle, SelAtoms, nSelAtoms);
      if (nSelAtoms > 0) {
         mmdb::Atom *at = SelAtoms[0];
         if (attribute_name == "x")
            at->x = val;
         if (attribute_name == "y")
            at->y = val;
         if (attribute_name == "z")
            at->z = val;
         if (attribute_name == "B")
            at->tempFactor = val;
         if (attribute_name == "b")
            at->tempFactor = val;
         if (attribute_name == "occ")
            at->occupancy = val;
      }
      atom_sel.mol->DeleteSelection(SelectionHandle);
   }
   have_unsaved_changes_flag = 1;
   atom_sel.mol->FinishStructEdit();
   make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
   return istate;
}

int molecule_class_info_t::swap_atom_alt_conf(std::string chain_id, int res_no, std::string ins_code,
                                              std::string atom_name, std::string alt_conf) {

  int istate = 0;
   if (atom_sel.n_selected_atoms > 0) {
      int SelectionHandle = atom_sel.mol->NewSelection();
      atom_sel.mol->SelectAtoms(SelectionHandle, 0,
                                chain_id.c_str(),
                                res_no, ins_code.c_str(),
                                res_no, ins_code.c_str(),
                                "*",
                                atom_name.c_str(),
                                "*",
                                "*");
      int nSelAtoms;
      mmdb::PPAtom SelAtoms;
      atom_sel.mol->GetSelIndex(SelectionHandle, SelAtoms, nSelAtoms);
      if (nSelAtoms > 1) {
         mmdb::Atom *at_0 = SelAtoms[0];
         mmdb::Atom *at_1 = SelAtoms[1];
         std::string alt_conf_0 = at_0->altLoc;
         std::string alt_conf_1 = at_1->altLoc;
         if (alt_conf_0.length() <= 8){
            if (alt_conf_1.length() <= 8){
               strncpy(at_0->altLoc, alt_conf_1.c_str(), alt_conf_1.length()+1);
               strncpy(at_1->altLoc, alt_conf_0.c_str(), alt_conf_0.length()+1);
            }
         }
      }
   }
   have_unsaved_changes_flag = 1;
   atom_sel.mol->FinishStructEdit();
   make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
   return istate;
 }


int
molecule_class_info_t::set_atom_string_attribute(std::string chain_id, int resno, std::string ins_code,
                                                 std::string atom_name, std::string alt_conf,
                                                 std::string attribute_name, std::string val_str) {

   int istate = 0;
   if (atom_sel.n_selected_atoms > 0) {
      int SelectionHandle = atom_sel.mol->NewSelection();
      atom_sel.mol->SelectAtoms(SelectionHandle, 0,
                                chain_id.c_str(),
                                resno, ins_code.c_str(),
                                resno, ins_code.c_str(),
                                "*",
                                atom_name.c_str(),
                                "*",
                                alt_conf.c_str());
      int nSelAtoms;
      mmdb::PPAtom SelAtoms;
      atom_sel.mol->GetSelIndex(SelectionHandle, SelAtoms, nSelAtoms);
      if (nSelAtoms > 0) {
         mmdb::Atom *at = SelAtoms[0];
         if (attribute_name == "atom-name")
            at->SetAtomName((char *)val_str.c_str());
         if (attribute_name == "alt-conf") {
            strncpy(at->altLoc, val_str.c_str(), 2);
         }
         if (attribute_name == "element") {
            at->SetElementName(val_str.c_str());
         }
         if (attribute_name == "segid") {
            strncpy(at->segID, val_str.c_str(), 4);
         }
      }
      have_unsaved_changes_flag = 1;
      atom_sel.mol->FinishStructEdit();
      make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
   }
   return istate;
}

int
molecule_class_info_t::set_atom_attributes(const std::vector<coot::atom_attribute_setting_t> &v) {

   int istate = 0;
   if (has_model()) {
      if (v.size() > 0) {
         make_backup();
         for (unsigned int iv=0; iv<v.size(); iv++) {
            int SelectionHandle = atom_sel.mol->NewSelection();
            atom_sel.mol->SelectAtoms(SelectionHandle, 0,
                                      v[iv].atom_spec.chain_id.c_str(),
                                      v[iv].atom_spec.res_no, v[iv].atom_spec.ins_code.c_str(),
                                      v[iv].atom_spec.res_no, v[iv].atom_spec.ins_code.c_str(),
                                      "*",
                                      v[iv].atom_spec.atom_name.c_str(),
                                      "*",
                                      v[iv].atom_spec.alt_conf.c_str());
            int nSelAtoms;
            mmdb::PPAtom SelAtoms;
            atom_sel.mol->GetSelIndex(SelectionHandle, SelAtoms, nSelAtoms);

            if (nSelAtoms > 0) {
               mmdb::Atom *at = SelAtoms[0];
//                std::cout << "DEBUG:: "
//                          << v[iv].attribute_value.type << " " << v[iv].attribute_name << " :"
//                          << v[iv].attribute_value.s << ": " << v[iv].attribute_value.val
//                          << std::endl;
               if (v[iv].attribute_value.type == coot::atom_attribute_setting_help_t::IS_STRING) {
                  if (v[iv].attribute_name == "atom-name")
                     at->SetAtomName(v[iv].attribute_value.s.c_str());
                  if (v[iv].attribute_name == "alt-conf") {
                     strncpy(at->altLoc, v[iv].attribute_value.s.c_str(), 2);
                  }
                  if (v[iv].attribute_name == "element") {
                     at->SetElementName(v[iv].attribute_value.s.c_str());
                  }
                  if (v[iv].attribute_name == "segid") {
                     for (int ii=0; ii<20; ii++) // 20 is magic number - for char[20] SegID
                        at->segID[ii] = '\n';
                     strncpy(at->segID, v[iv].attribute_value.s.c_str(), 19);
                  }
               }
               if (v[iv].attribute_value.type == coot::atom_attribute_setting_help_t::IS_FLOAT) {
                  if (v[iv].attribute_name == "x")
                     at->x = v[iv].attribute_value.val;
                  if (v[iv].attribute_name == "y")
                     at->y = v[iv].attribute_value.val;
                  if (v[iv].attribute_name == "z")
                     at->z = v[iv].attribute_value.val;
                  if (v[iv].attribute_name == "b")
                     at->tempFactor = v[iv].attribute_value.val;
                  if (v[iv].attribute_name == "B")
                     at->tempFactor = v[iv].attribute_value.val;
                  if (v[iv].attribute_name == "occ")
                     at->occupancy = v[iv].attribute_value.val;
               }
            }
            atom_sel.mol->DeleteSelection(SelectionHandle);
         }
         have_unsaved_changes_flag = 1;
         atom_sel.mol->FinishStructEdit();
         make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
      }
   }
   return istate;
}

void
molecule_class_info_t::set_residue_name(std::string chain_id, int res_no, std::string ins_code,
                                        std::string new_name) {

   make_backup();
   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         std::string mol_chain_id = chain_p->GetChainID();
         if (mol_chain_id == chain_id) {
            int nres = chain_p->GetNumberOfResidues();
            mmdb::Residue *residue_p;
            for (int ires=0; ires<nres; ires++) {
               residue_p = chain_p->GetResidue(ires);
               if (res_no == residue_p->GetSeqNum()) {
                  if (ins_code == residue_p->GetInsCode()) {
                     residue_p->SetResName(new_name.c_str());
                  }
               }
            }
         }
      }
   }
   have_unsaved_changes_flag = 1;
}



// -----------------------------------------------------------------------------
//                     pepflip
// -----------------------------------------------------------------------------

// This is the function called by graphics-info-define (handling the atom click)
//
void
molecule_class_info_t::pepflip(int atom_index) {

   // std::cout << "-------------------- pepflip(atom_index) ---------"  << std::endl;
   const char *chain_id = atom_sel.atom_selection[atom_index]->residue->GetChainID();
   int resno = atom_sel.atom_selection[atom_index]->residue->seqNum;
   std::string atom_name = atom_sel.atom_selection[atom_index]->name;
   std::string inscode =  atom_sel.atom_selection[atom_index]->GetInsCode();
   std::string altconf =  atom_sel.atom_selection[atom_index]->altLoc;

   std::cout << "INFO:: flipping " << resno << " " << altconf << " "
             << chain_id << std::endl;
   if (atom_name == " N  ")
      resno--;
   if (atom_name == " H  ")
      resno--;

   pepflip_residue(chain_id, resno, inscode, altconf);

}
// model_refine_dialog_pepflip_button

int
molecule_class_info_t::pepflip_residue(const std::string &chain_id,
                                       int resno,
                                       const std::string &ins_code,
                                       const std::string &alt_conf) {

   // std::cout << "-------------------- pepflip_residue() ---------"  << std::endl;

   make_backup(); // must do it here, no intermediate.
   int iresult = coot::pepflip(atom_sel.mol, chain_id, resno, ins_code, alt_conf);
   if (iresult == 1) {
      std::cout << "INFO:: flipped " << resno << " " << chain_id << " success" << std::endl;
      make_bonds_type_checked(__FUNCTION__);
      have_unsaved_changes_flag = 1;
   } else {
      std::cout << "pepflip failed " << chain_id << " " << resno << std::endl;
   }
   return iresult;
}

void
molecule_class_info_t::pepflip(const coot::atom_spec_t &spec) {

   // std::cout << "-------------------- pepflip(spec) ---------"  << std::endl;
   // std::cout << "debug:: spec: " << spec << std::endl;
   std::string alt_conf = spec.alt_conf;
   int res_no = spec.res_no;
   if (spec.atom_name == " N  ")
      res_no--;
   if (spec.atom_name == " H  ")
      res_no--;

   pepflip_residue(spec.chain_id, res_no, spec.ins_code, alt_conf);
}

graphical_bonds_container
molecule_class_info_t::make_environment_bonds_box(const coot::residue_spec_t &residue_spec,
                                                  coot::protein_geometry *protein_geom_p) const {

   graphics_info_t g;
   graphical_bonds_container bonds_box_env;
   mmdb::Residue *residue_p = coot::util::get_residue(residue_spec, atom_sel.mol);

   if (residue_p) {

      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
      if (nResidueAtoms == 0) {
         std::cout << " something broken in atom residue selection in ";
         std::cout << "make_environment_bonds_box: got " << nResidueAtoms
                   << " atoms " << std::endl;
      } else {

         bool residue_is_water_flag = false;
         bool draw_bonds_to_hydrogens_flag = draw_hydrogens_flag; // class var
         std::string residue_name = residue_p->GetResName();
         if (residue_name == "HOH" || residue_name == "WAT")
            residue_is_water_flag = 1;
         Bond_lines_container bonds(atom_sel,residue_atoms, nResidueAtoms,
                                    protein_geom_p,
                                    residue_is_water_flag,
                                    draw_bonds_to_hydrogens_flag,
                                    g.environment_min_distance,
                                    g.environment_max_distance);
         bonds_box_env = bonds.make_graphical_bonds();
      }
   } else {
      std::cout << "ERROR:: NULL residue_p in make_environment_bonds_box() " << std::endl;
   }
   return bonds_box_env;
}

// pass min and max dist - don't use graphics_info_t!
graphical_bonds_container
molecule_class_info_t::make_environment_bonds_box(int atom_index,
                                                  coot::protein_geometry *protein_geom_p) const {

   graphics_info_t g;
   graphical_bonds_container bonds_box_env;

   if ((atom_index >= atom_sel.n_selected_atoms) || (atom_index < 0)) {

      std::cout << "ERROR:: trapped an atom index problem in make_environment_bonds_box()!!!\n"
                << "        Tell Paul - he wants to know...." << std::endl;
      std::cout << "ERROR:: " << atom_index << " " << atom_sel.n_selected_atoms << std::endl;
   } else {

      mmdb::Atom *point_atom_p = atom_sel.atom_selection[atom_index];

      mmdb::PPResidue SelResidues;
      int nSelResdues;

      int ires = point_atom_p->GetSeqNum();
      const char *chain_id = point_atom_p->GetChainID();
      mmdb::Residue *residue_p = point_atom_p->residue;

      if (residue_p) {
         coot::residue_spec_t residue_spec(residue_p);
         return make_environment_bonds_box(residue_spec, protein_geom_p);
      }
   }
   return bonds_box_env;
}

graphical_bonds_container
molecule_class_info_t::make_symmetry_environment_bonds_box(int atom_index,
                                                           coot::protein_geometry *protein_geom_p) const {
   graphical_bonds_container env_bonds_box;

   // std::cout << ":: entering make_symmetry_environment_bonds_box" << std::endl;
   if (atom_sel.atom_selection != NULL) {

      graphics_info_t g;
      mmdb::PPResidue SelResidues;
      int nSelResdues;

      // First select all the atoms in this residue:
      //
      if ((atom_index >= atom_sel.n_selected_atoms) || (atom_index < 0)) {

         std::cout << "ERROR:: trapped an atom index problem in make_symmetry_environment_bonds_box()!!!\n"
                   << "        Tell Paul - he wants to know...." << std::endl;
         std::cout << "ERROR:: " << atom_index << " " << atom_sel.n_selected_atoms << std::endl;
      } else {
         // Happy (normal) path
         mmdb::PAtom point_atom_p = atom_sel.atom_selection[atom_index];
         int ires = point_atom_p->GetSeqNum();
         char *chain_id = point_atom_p->GetChainID();

         int selHnd = atom_sel.mol->NewSelection();
         atom_sel.mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1,
                               chain_id, // chains
                               ires,"*", // starting res
                               ires,"*", // ending res
                               "*",  // residue name
                               "*",  // Residue must contain this atom name?
                               "*",  // Residue must contain this Element?
                               "*",  // altLocs
                               mmdb::SKEY_NEW // selection key
                               );
         atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResdues);

         if (nSelResdues != 1) {
            std::cout << " something broken in residue selection in ";
            std::cout << "make_environment_bonds_box: got " << nSelResdues
                      << " residues " << std::endl;
         } else {
            mmdb::PPAtom residue_atoms;
            int nResidueAtoms;
            SelResidues[0]->GetAtomTable(residue_atoms, nResidueAtoms);
            if (nResidueAtoms == 0) {
               std::cout << " something broken in atom residue selection in ";
               std::cout << "make_environment_bonds_box: got " << nResidueAtoms
                         << " atoms " << std::endl;
            } else {

               short int do_symmetry = 1;
               // std::cout << "... calling Bond_lines_container constructor" << std::endl;
               bool draw_bonds_to_hydrogens_flag = draw_hydrogens_flag; // class var
               Bond_lines_container bonds(atom_sel, residue_atoms, nResidueAtoms,
                                          g.environment_max_distance,
                                          g.environment_min_distance,
                                          draw_bonds_to_hydrogens_flag,
                                          do_symmetry);
               env_bonds_box = bonds.make_graphical_bonds();
            }
         }
         atom_sel.mol->DeleteSelection(selHnd);
      }
   }
   return env_bonds_box;
}

// return "N', "C" or "not-terminal-residue"
std::string
molecule_class_info_t::get_term_type_old(int atom_index) {

   std::string term_type;

   char *chainid = atom_sel.atom_selection[atom_index]->GetChainID();
   int ires_atom = atom_sel.atom_selection[atom_index]->GetSeqNum();
   mmdb::PChain chain = atom_sel.mol->GetChain(1,chainid);
   int nres = chain->GetNumberOfResidues();
   int lowest_res_no = 99999;
   int highest_res_no = -99999;
   for (int ires=0; ires<nres; ires++) {
      mmdb::PResidue res = chain->GetResidue(ires);
      if (res) { // could have been deleted (NULL)
         if (res->GetSeqNum() > highest_res_no) {
            highest_res_no = res->GetSeqNum();
         }
         if (res->GetSeqNum() < lowest_res_no) {
            lowest_res_no = res->GetSeqNum();
         }
      }
   }
   if (ires_atom == lowest_res_no) {
      term_type = "N";
   } else {
      if (ires_atom == highest_res_no) {
         term_type = "C";
      } else {
         term_type = "not-terminal-residue";
      }
   }
   return term_type;
}

// The atom_index is the atom index of the clicked atom.
//
// Initially, this routine tested for real terminii.
//
// The one day EJD asked me how I would build a few missing residues
// (in a loop or so) that arp-warp missed.  I said "add terminal
// residue".  Then I tried it and it failed, of course because that
// residue was not a real terminus.  On reflection, I don't think that
// it should fail in this situation.
//
// So, I don't want to test for a real terminus.
//
// Let's see if this residue has a residue on one side of it, but not
// the other.  If so, return N or C depending on whether the other
// residue is upstream or not.
//
// Return not-terminal-residue if 0 neighbours, return "M" (for mid)
// for both neighbours present (used by
// graphics_info_t::execute_add_terminal_residue()).  Realise that
// this is a bit of a kludge, because usually, the terminal type
// refers to the residue that we clicked on not the missing residue.
//
// In the "M" case, we refer to the missing residue.
//
// "singleton" is a possibile terminal type - for cases where this
// residue does not have neighbours.
//
// Note that this ignores altlocs and insertion codes. It should do
// altlocs at least.
//
std::string
molecule_class_info_t::get_term_type(int atom_index) const {

   std::string term_type = "not-terminal-residue"; // returned thing

   if (atom_index < 0) return "";
   if (atom_index >= atom_sel.n_selected_atoms) return "";

   int ires_atom = atom_sel.atom_selection[atom_index]->GetSeqNum();
   mmdb::PChain chain = atom_sel.atom_selection[atom_index]->GetChain();
   int nres = chain->GetNumberOfResidues();

   // including tests needed for single missing residue:
   short int has_up_neighb = 0;
   short int has_down_neighb = 0;
   short int has_up_up_neighb = 0;
   short int has_down_down_neighb = 0;

   // Check for neighbouring residues to the clicked atom. Don't count
   // waters as neighbours.
   //
   for (int ires=0; ires<nres; ires++) {
      mmdb::PResidue res = chain->GetResidue(ires);
      if (res) { // could have been deleted (NULL)
         if (res->GetSeqNum() == (ires_atom + 1))
            has_up_neighb = 1;
         if (res->GetSeqNum() == (ires_atom + 2))
            has_up_up_neighb = 1;
         if (res->GetSeqNum() == (ires_atom - 1))
            has_down_neighb = 1;
         if (res->GetSeqNum() == (ires_atom - 2))
            has_down_down_neighb = 1;
      }
   }

   if ( (has_up_neighb + has_down_neighb) == 1 ) {
      if (has_up_neighb)
         term_type = "N";
      if (has_down_neighb)
         term_type = "C";
   }

   if ((has_up_neighb == 0) && (has_down_neighb == 0))
      term_type = "singleton";

   // Now test for missing single residue, "M" (mid):
   //
   if ( (!has_up_neighb) && has_up_up_neighb)
      term_type = "MC"; // missing middle res, treat as C term

   if ( (!has_down_neighb) && has_down_down_neighb)
      term_type = "MN"; // missing middle res, treat as N term

   // std::cout << "DEBUG:: get_term_type Returning residue type " << term_type << std::endl;

   return term_type;
}

std::string
molecule_class_info_t::get_term_type(mmdb::Atom *atom) const {

   std::string term_type = "not-terminal-residue"; // returned thing

   int ires_atom      = atom->GetSeqNum();
   mmdb::Chain *chain = atom->GetChain();
   int nres = chain->GetNumberOfResidues();

   // including tests needed for single missing residue:
   short int has_up_neighb = 0;
   short int has_down_neighb = 0;
   short int has_up_up_neighb = 0;
   short int has_down_down_neighb = 0;

   // Check for neighbouring residues to the clicked atom. Don't count
   // waters as neighbours.
   //
   for (int ires=0; ires<nres; ires++) {
      mmdb::PResidue res = chain->GetResidue(ires);
      if (res) { // could have been deleted (NULL)
         if (res->GetSeqNum() == (ires_atom + 1))
            has_up_neighb = 1;
         if (res->GetSeqNum() == (ires_atom + 2))
            has_up_up_neighb = 1;
         if (res->GetSeqNum() == (ires_atom - 1))
            has_down_neighb = 1;
         if (res->GetSeqNum() == (ires_atom - 2))
            has_down_down_neighb = 1;
      }
   }

   if ( (has_up_neighb + has_down_neighb) == 1 ) {
      if (has_up_neighb)
         term_type = "N";
      if (has_down_neighb)
         term_type = "C";
   }

   if ((has_up_neighb == 0) && (has_down_neighb == 0))
      term_type = "singleton";

   // Now test for missing single residue, "M" (mid):
   //
   if ( (!has_up_neighb) && has_up_up_neighb)
      term_type = "MC"; // missing middle res, treat as C term

   if ( (!has_down_neighb) && has_down_down_neighb)
      term_type = "MN"; // missing middle res, treat as N term

   // std::cout << "DEBUG:: get_term_type Returning residue type " << term_type << std::endl;

   return term_type;
}

// Replace the atoms in this molecule by those in the given atom selection.
int
molecule_class_info_t::replace_fragment(atom_selection_container_t asc) {

   if (! asc.mol) return 0;

   bool move_zero_occ = 1;

   // replace an atom if you can, otherwise create a new atom (and a new residue and chain if needed)

   make_backup();
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idx = -1;
      mmdb::Atom *at = asc.atom_selection[i];

      if (! at->isTer()) {
         // can we find the atom with fast indexing?
         if (asc.UDDOldAtomIndexHandle >= 0) {
            // OK for fast atom indexing
            int ref_index = -1;
            if (at->GetUDData(asc.UDDOldAtomIndexHandle, ref_index) == mmdb::UDDATA_Ok) {
               if (ref_index >= 0) {
                  if (moving_atom_matches(at, ref_index)) {
                     idx = ref_index; // yay.
                  }
               }
            }
         }

         if (idx == -1) {
            idx = full_atom_spec_to_atom_index(coot::atom_spec_t(at));
         }

         if (idx != -1) {
            mmdb::Atom *ref_atom = atom_sel.atom_selection[idx];
            ref_atom->x = at->x;
            ref_atom->y = at->y;
            ref_atom->z = at->z;

         } else {

            // add the atom
            mmdb::Chain *chain_p = get_chain(at->GetChainID());
            mmdb::Residue *residue_p = get_residue(coot::residue_spec_t(coot::atom_spec_t(at)));

            if (! chain_p) {
               int imod = 1;
               mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
               if (model_p) {
                  mmdb::Chain *chain_p = new mmdb::Chain;
                  chain_p->SetChainID(at->GetChainID());
                  residue_p = new mmdb::Residue;
                  residue_p->seqNum = at->GetSeqNum();
                  residue_p->SetResName(at->residue->GetResName());
                  chain_p->AddResidue(residue_p);
                  model_p->AddChain(chain_p);
                  atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                  atom_sel.mol->FinishStructEdit();
               }
            } else {
               if (residue_p) {
                  // std::cout << "   ======= found the residue " << std::endl;
               } else {
                  residue_p = new mmdb::Residue;
                  residue_p->SetResID(at->residue->GetResName(), at->residue->seqNum, at->residue->insCode);
                  int res_no = at->GetSeqNum();
                  std::string ins_code(at->GetInsCode());
                  std::pair<int, mmdb::Residue *> sn =
                     find_serial_number_for_insert(res_no, ins_code, chain_p->GetChainID());

                  if (sn.first != -1) { // normal insert

                     int n_residues_before = chain_p->GetNumberOfResidues();
                     int n_chain_residues = chain_p->InsResidue(residue_p, sn.first);
                     mmdb::Residue *res_after_p = get_residue(coot::residue_spec_t(coot::atom_spec_t(at)));

                  } else {
                     chain_p->AddResidue(residue_p);
                     atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);// needed?
                     atom_sel.mol->FinishStructEdit();
                  }
               }
            }

            if (residue_p) {
               mmdb::Atom *at_copy(at);
               residue_p->AddAtom(at_copy);
               // residue_p->TrimAtomTable();
            }
         }
      }
   }

   atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);// needed?
   coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
   atom_sel.mol->FinishStructEdit();

   atom_sel = make_asc(atom_sel.mol);
   have_unsaved_changes_flag = 1;
   if (show_symmetry)
      update_symmetry();
   make_bonds_type_checked(__FUNCTION__);
   return 1;
}


// Delete residue: the whole residue.  We decide outside if this
// function or delete_residue_with_altconf should be used.
//
// Return the success status (0: failure to delete, 1 is deleted)
//
// regenerate the graphical bonds box if necessary.
//
// model number is either a specific model number of mmdb::MinInt4, meaning:
// any/all model(s).
//
// if we have delete_zone mode then we don't want to update the ghosts or the gui
// or make backups
short int
molecule_class_info_t::delete_residue(int model_number,
                                      const std::string &chain_id, int resno,
                                      const std::string &ins_code) {

   short int was_deleted = 0;
   mmdb::Chain *chain;

   // run over chains of the existing mol
   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) {

      if (false)
         std::cout << "debug:: delete_residue() comparing imod: "
                   << imod << " and model_number "
                   << model_number << std::endl;

      if ((imod == model_number) || (model_number == mmdb::MinInt4)) {

         int nchains = atom_sel.mol->GetNumberOfChains(imod);
         for (int ichain=0; ichain<nchains; ichain++) {

            chain = atom_sel.mol->GetChain(imod,ichain);
            std::string mol_chain_id(chain->GetChainID());

            if (chain_id == mol_chain_id) {

               // std::cout << "debug:: matching chain_ids on  " << chain_id << std::endl;

               int nres = chain->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *res = chain->GetResidue(ires);
                  if (res) {
                     if (res->GetSeqNum() == resno) {

                        // so we have a matching residue:
                        int iseqno = res->GetSeqNum();
                        mmdb::pstr inscode = res->GetInsCode();
                        std::string inscodestr(inscode);
                        if (ins_code == inscodestr) {
                           make_backup();
                           atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                           delete_ghost_selections();
                           chain->DeleteResidue(iseqno, inscode);
                           was_deleted = 1;
                           res = NULL;
                           break;
                        }
                     }
                  }
               }
            }
            if (was_deleted)
               break;
         }
      }
   }

   if (was_deleted) {

      // we can't do this after the modification: it has to be done before
      // atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);

      atom_sel.atom_selection = NULL;
      coot::residue_spec_t spec(model_number, chain_id, resno, ins_code);
      delete_any_link_containing_residue(spec);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
      //
      trim_atom_label_table();
      update_symmetry();
   }
   return was_deleted;
}

// wraps above
short int
molecule_class_info_t::delete_residue(const coot::residue_spec_t &spec) {

   short int status = 0;
   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imodel=1; imodel<=n_models; imodel++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imodel);
      if (model_p)
         status = delete_residue(imodel, spec.chain_id, spec.res_no, spec.ins_code);
   }
   return status;

}

short int
molecule_class_info_t::delete_residue_hydrogens(const std::string &chain_id, int resno,
                                                const std::string &ins_code,
                                                const std::string &altloc) {

   short int was_deleted = 0;
   mmdb::Chain *chain;

   // run over chains of the existing mol
   int nchains = atom_sel.mol->GetNumberOfChains(1);
   for (int ichain=0; ichain<nchains; ichain++) {

      chain = atom_sel.mol->GetChain(1,ichain);
      std::string mol_chain_id(chain->GetChainID());

      if (chain_id == mol_chain_id) {

         int nres = chain->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::PResidue res = chain->GetResidue(ires);
            if (res) {
               if (res->GetSeqNum() == resno) {

                  std::string local_ins_code = res->GetInsCode();
                  if (local_ins_code == ins_code) {

                     // Only make a backup, delete the atom selection
                     // (and proceed with the attempted deletion) if
                     // the residue has hydrogens.  And this should
                     // fix the esoteric crash that we get when
                     // deleting waters (because running this function
                     // meant that the atom selection was going out of
                     // date by the time it got to the
                     // delete_atom_by_atom_index() function.
                     //
                     bool residue_has_hydrogens = 0;
                     mmdb::PPAtom residue_atoms;
                     int nResidueAtoms;
                     res->GetAtomTable(residue_atoms, nResidueAtoms);
                     for (int iat=0; iat<nResidueAtoms; iat++) {
                        std::string ele(residue_atoms[iat]->element);
                        if (ele == " H") {
                           residue_has_hydrogens = 1;
                           break;
                        }
                        if (ele == " D") {
                           residue_has_hydrogens = 1;
                           break;
                        }
                     }

                     if (residue_has_hydrogens) {

                        // so we have a matching residue:
                        make_backup();
                        atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                        delete_ghost_selections();
                        was_deleted = 1; // This is moved here because we
                        // don't want to DeleteSelection on
                        // multiple atom delete, so we set
                        // it before the delete actually
                        // happens - it does mean that in
                        // the pathological case
                        // was_deleted is set but no atoms
                        // get deleted - which is harmless
                        // (I think).

                        // getatom table here and check the elements in the lst

                        if (nResidueAtoms == 0) {
                           std::cout << "WARNING:: No atoms in residue (strange!)\n";
                        } else {
                           for (int i=0; i<nResidueAtoms; i++) {
                              std::string element(residue_atoms[i]->element);
                              if (element == " H" || element == " D") {
                                 res->DeleteAtom(i);
                                 // was_deleted = 1; // done upstairs
                              }
                           }
                           if (was_deleted)
                              res->TrimAtomTable();
                        }
                     }
                  }
               }
            }
         }
      }
      if (was_deleted)
         break;
   }

   // potentially
   if (was_deleted) {
      atom_sel.atom_selection = NULL;
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
      //
      trim_atom_label_table();
      update_symmetry();
   }
   return was_deleted;
}


// Delete only the atoms of the residue that have the same altconf (as
// the selected atom).  If the selected atom has altconf "", you
// should call simply delete_residue().
//
// Return 1 if at least one atom was deleted, else 0.
//
short int
molecule_class_info_t::delete_residue_with_full_spec(int imodel,
                                                     const std::string &chain_id,
                                                     int resno,
                                                     const std::string &ins_code,
                                                     const std::string &altconf) {

   short int was_deleted = 0;
   mmdb::Chain *chain;
   mmdb::Residue *residue_for_deletion = NULL;
   std::vector<std::pair<std::string, float> > deleted_atom;
//    std::cout << "DEBUG:: start of delete-residue-with-altconf n-atoms: "
//              << atom_sel.n_selected_atoms << std::endl;

   // run over chains of the existing mol

   int n_models = atom_sel.mol->GetNumberOfModels();

   for (int imod=1; imod<=n_models; imod++) {

      if ((imod == imodel) || (imodel == mmdb::MinInt4)) {

         int nchains = atom_sel.mol->GetNumberOfChains(imod);
         for (int ichain=0; ichain<nchains; ichain++) {

            chain = atom_sel.mol->GetChain(imod,ichain);
            std::string mol_chain_id(chain->GetChainID());

            if (chain_id == mol_chain_id) {

               int nres = chain->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::PResidue res = chain->GetResidue(ires);
                  if (res) {
                     if (res->GetSeqNum() == resno) {

                        mmdb::pstr inscode_res = res->GetInsCode();
                        std::string inscode_res_str(inscode_res);

                        if (inscode_res_str == ins_code) {

                           // so we have a matching residue:
                           residue_for_deletion = res;

                           // delete the specific atoms of the residue:
                           mmdb::PPAtom atoms;
                           int n_atoms;
                           bool have_deletable_atom = 0;
                           res->GetAtomTable(atoms, n_atoms);
                           for (int i=0; i<n_atoms; i++) {
                              if (std::string(atoms[i]->altLoc) == altconf) {
                                 have_deletable_atom = 1;
                              }
                           }

                           if (have_deletable_atom) {
                              make_backup();
                              atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                              delete_ghost_selections();

                              for (int i=0; i<n_atoms; i++) {
                                 if (std::string(atoms[i]->altLoc) == altconf) {
                                    std::pair<std::string, float> p(atoms[i]->name,
                                                                    atoms[i]->occupancy);
                                    deleted_atom.push_back(p);
                                    res->DeleteAtom(i);
                                    was_deleted = 1;
                                 }
                                 // delete pointer if all atoms were deleted
                                 // should probably delete res too!?
                                 unsigned int ui_n_atoms = n_atoms;
                                 if (deleted_atom.size() == ui_n_atoms) {
                                    residue_for_deletion = NULL;
                                 }
                              }
                           }
                        }
                        break;
                     }
                  }
               }
            }
         }
      }
      // jump out (only) if something was deleted and we are not
      // deleting from all chains.
      if (was_deleted && (imodel != mmdb::MinInt4))
         break;
   }

   // potentially (usually, I imagine)
   if (was_deleted) {
      atom_sel.atom_selection = NULL;
      atom_sel.mol->FinishStructEdit();
      trim_atom_label_table();
      unalt_conf_residue_atoms(residue_for_deletion);
      atom_sel.mol->FinishStructEdit();

      // Now, add the occpancy of the delete atom to the remaining
      // atoms of that residue with that atom name (if there are any)
      //
      if (residue_for_deletion) {
         if (deleted_atom.size() > 0) {
            if (! is_from_shelx_ins_flag) {
               mmdb::PPAtom atoms = NULL;
               int n_atoms;
               residue_for_deletion->GetAtomTable(atoms, n_atoms);
               for (unsigned int idat=0; idat<deleted_atom.size(); idat++) {
                  std::vector<mmdb::Atom *> same_name;
                  for (int iresat=0; iresat<n_atoms; iresat++) {
                     if (deleted_atom[idat].first == std::string(atoms[iresat]->name)) {
                        same_name.push_back(atoms[iresat]);
                     }
                  }
                  if (same_name.size() > 0) {
                     float extra_occ = deleted_atom[idat].second / float(same_name.size());
                     for (unsigned int isn=0; isn<same_name.size(); isn++) {
                        // std::cout << "Adding " << extra_occ << " to "
                        // << same_name[isn]->occupancy
                        // << " for atom " << same_name[isn]->name << std::endl;
                        same_name[isn]->occupancy += extra_occ;
                     }
                  }
               }
            } else {
               // was a shelx molecule, if we are down to one alt conf
               // now, then we set the occupancy to 11.0, otherwise
               // leave it as it is.
               mmdb::PPAtom atoms = NULL;
               int n_atoms;
               residue_for_deletion->GetAtomTable(atoms, n_atoms);
               for (unsigned int idat=0; idat<deleted_atom.size(); idat++) {
                  std::vector<mmdb::Atom *> same_name;
                  for (int iresat=0; iresat<n_atoms; iresat++) {
                     if (deleted_atom[idat].first == std::string(atoms[iresat]->name)) {
                        same_name.push_back(atoms[iresat]);
                     }
                  }
                  if (same_name.size() == 1) {
                     same_name[0]->occupancy = 11.0;
                  }
               }
            }
         }
      }

      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
      update_symmetry();
   }
//    std::cout << "DEBUG:: End of delete-residue-with-altconf n-atoms: "
//              << atom_sel.n_selected_atoms << std::endl;
   return was_deleted;
}

// Return 1 if at least one atom was deleted, else 0.
//
short int
molecule_class_info_t::delete_residues(const std::vector<coot::residue_spec_t> &specs) {

   bool something_deleted = false;
   mmdb::Manager *mol = atom_sel.mol;

   // make a backup if there are residues to delete
   for (unsigned int i=0; i<specs.size(); i++) {
      const coot::residue_spec_t &spec = specs[i];
      mmdb::Residue *residue_p = get_residue(spec);
      if (residue_p) {
         make_backup();
         break;
      }
   }

   for (unsigned int i=0; i<specs.size(); i++) {
      const coot::residue_spec_t &spec = specs[i];
      mmdb::Residue *residue_p = get_residue(spec);
      if (residue_p) {
         mmdb::Chain *chain_p = residue_p->GetChain();
         if (chain_p) {
            bool a_link_was_deleted = coot::util::delete_residue_references_in_header_info(residue_p, mol);
            delete residue_p;
            something_deleted = true;
         }
      }
   }

   if (something_deleted) {
      atom_sel.atom_selection = NULL;
      atom_sel.mol->FinishStructEdit();
      trim_atom_label_table();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
      update_symmetry();
   }
   return something_deleted;
}

short int
molecule_class_info_t::delete_residue_sidechain(const coot::residue_spec_t &rs) {

   return delete_residue_sidechain(rs.chain_id, rs.res_no, rs.ins_code);
}



short int
molecule_class_info_t::delete_residue_sidechain(const std::string &chain_id,
                                                int resno,
                                                const std::string &inscode) {


   short int was_deleted = 0;
   mmdb::Chain *chain;
   mmdb::Residue *residue_for_deletion = NULL;

   // run over chains of the existing mol
   int nchains = atom_sel.mol->GetNumberOfChains(1);
   for (int ichain=0; ichain<nchains; ichain++) {

      chain = atom_sel.mol->GetChain(1,ichain);
      std::string mol_chain_id(chain->GetChainID());

      if (chain_id == mol_chain_id) {

         int nres = chain->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::PResidue res = chain->GetResidue(ires);
            if (res) {
               if (res->GetSeqNum() == resno) {

                  mmdb::pstr inscode_res = res->GetInsCode();
                  std::string inscode_res_str(inscode_res);

                  if (inscode_res_str == inscode) {

                     // so we have a matching residue:
                     residue_for_deletion = res;
                     make_backup();
                     atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                     delete_ghost_selections();
                     was_deleted = 1; // need to set this here because
                                      // we need to regenerate the
                                      // atom selection.
                     // delete the specific atoms of the residue:
                     mmdb::PPAtom atoms = 0;
                     int n_atoms = 0;
                     res->GetAtomTable(atoms, n_atoms);
                     for (int i=0; i<n_atoms; i++) {
                        if (! (coot::is_main_chain_or_cb_p(atoms[i]))) {
                           res->DeleteAtom(i);
                        }
                     }
                     if (was_deleted)
                        res->TrimAtomTable();
                  }
               }
            }
         }
      }
   }
   // potentially (usually, I imagine)
   if (was_deleted) {
      atom_sel.atom_selection = NULL;
      atom_sel.mol->FinishStructEdit();
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel = make_asc(atom_sel.mol);
      trim_atom_label_table();
      unalt_conf_residue_atoms(residue_for_deletion);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
   }
   return was_deleted;
}



int
molecule_class_info_t::delete_zone(const coot::residue_spec_t &res1,
                                   const coot::residue_spec_t &res2) {

   int first_res = res1.res_no;
   int last_res  = res2.res_no;
   // std::string alt_conf = res1.altconf;  // FIXME
   std::string alt_conf = "";
   std::string inscode = ""; // FIXME, more cleverness required.
                             // A range of values?

   if (first_res > last_res) {
      int tmp = last_res;
      last_res = first_res;
      first_res = tmp;
   }

   make_backup();
   // temporarily turn off backups when we delete this range:
   //
   int tmp_backup_this_molecule = backup_this_molecule;
   backup_this_molecule = 0;

   std::cout << "DEBUG:: in delete_zone() we have model numbers "
             << res1.model_number << " and "
             << res2.model_number << std::endl;

   // delete any/all residues in range, unless the both models in the
   // residue specs were set to something interesting (a specific
   // model number) and are the same as each other.
   //
   int model_number_ANY = mmdb::MinInt4;
   if (res1.model_number != mmdb::MinInt4)
      if (res2.model_number != mmdb::MinInt4)
         if (res1.model_number == res2.model_number)
            model_number_ANY = res1.model_number;

   bool was_deleted = false;

   std::vector<coot::residue_spec_t> deleted_residue_specs;
   std::vector<mmdb::Residue *> deleted_residues;

   // run over chains of the existing mol
   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) {

      int nchains = atom_sel.mol->GetNumberOfChains(imod);
      for (int ichain=0; ichain<nchains; ichain++) {

         mmdb::Chain *chain_p = atom_sel.mol->GetChain(imod, ichain);
         if (chain_p) {
            std::string mol_chain_id(chain_p->GetChainID());

            if (res1.chain_id == mol_chain_id) {

               int nres = chain_p->GetNumberOfResidues();
               // for (int ires=0; ires<nres; ires++) {
               for (int ires=nres-1; ires>=0; ires--) {
                  mmdb::Residue *res = chain_p->GetResidue(ires);
                  if (res) {
                     int res_no = res->GetSeqNum();
                     if (res_no >= first_res) {
                        if (res_no <= last_res) {
                           was_deleted = true;
                           deleted_residue_specs.push_back(coot::residue_spec_t(res));
                           deleted_residues.push_back(res);
                        }
                     }
                  }
               }
               // delete multiple residues like this, rather than chain_p->DeleteResidue(ires);
               for (unsigned int i=0; i<deleted_residues.size(); i++)
                  delete deleted_residues[i];
            }
         }
      }
   }
   backup_this_molecule = tmp_backup_this_molecule; // restore state

   if (was_deleted) {

      std::cout << "INFO... deleting links..." << std::endl;
      for (unsigned int ispec=0; ispec<deleted_residue_specs.size(); ispec++) {
         const coot::residue_spec_t &spec = deleted_residue_specs[ispec];
         delete_any_link_containing_residue(spec);
      }
      atom_sel.atom_selection = NULL;
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
      trim_atom_label_table();
      update_symmetry();
   }

   return 0;
}


void
molecule_class_info_t::unalt_conf_residue_atoms(mmdb::Residue *residue_p) {

   // We had a crash here (delete-residue-range), which is why there
   // is now more protection than you might imagine 20051231.

   if (residue_p) {
      mmdb::PPAtom atoms;
      int n_atoms;
      residue_p->GetAtomTable(atoms, n_atoms);
      std::cout << "There are " << n_atoms << " atoms in "
                << residue_p->GetChainID() << " " << residue_p->GetSeqNum()
                << std::endl;
      for (int i=0; i<n_atoms; i++) {
         std::string this_atom_name(atoms[i]->name);
         int n_match = 0;
         for (int j=0; j<n_atoms; j++) {
            if (atoms[j]) {
               std::string inner_name(atoms[j]->name);
               if (inner_name == this_atom_name) {
                  n_match++;
               }
            } else {
               std::cout << "ERROR:: null atom in unalt_conf" << std::endl;
            }
         }
         if (n_match == 1) {
	    std::string alt_conf(atoms[i]->altLoc);
            if (! alt_conf.empty()) {
	       std::string new_alt_conf("");
	       // force it down the atom's throat :) c.f. insert_coords_change_altconf
	       strncpy(atoms[i]->altLoc, new_alt_conf.c_str(), 2);
            }
         }
      }
   }
}

int
molecule_class_info_t::delete_water(const coot::atom_spec_t &atom_spec) {

   bool status = false;
   coot::residue_spec_t res_spec(atom_spec);
   mmdb::Residue *residue_p = get_residue(res_spec);
   if (residue_p) {
      std::string type = residue_p->GetResName();
      if (type == "HOH")
         status = delete_residue(res_spec);
   }
   return status;
}

bool
molecule_class_info_t::delete_atom(const coot::atom_spec_t &spec) {

   return delete_atom(spec.chain_id, spec.res_no, spec.ins_code, spec.atom_name, spec.alt_conf);
}

bool
molecule_class_info_t::delete_atom(const std::string &chain_id,
                                   int resno,
                                   const std::string &ins_code,
                                   const std::string &atname,
                                   const std::string &altconf) {

   short int was_deleted = 0;
   mmdb::Chain *chain;
   mmdb::Residue *residue_of_deleted_atom = NULL;

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
         strncpy(at->altLoc, "", 2);
         // set the occupancy to 1.0 of the remaining atom if it was not zero.
         if (at)
            if (at->occupancy > 0.009)
               at->occupancy = 1.0;
      }


      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked(__FUNCTION__);
      have_unsaved_changes_flag = 1;
      // unlikely to be necessary:
      trim_atom_label_table();
      update_symmetry();
   }
   return was_deleted;
}


int
molecule_class_info_t::delete_atoms(const std::vector<coot::atom_spec_t> &atom_specs) {

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
      make_bonds_type_checked(__FUNCTION__);
      have_unsaved_changes_flag = 1;
      // unlikely to be necessary:
      trim_atom_label_table();
      update_symmetry();
   }

   return n_deleleted_atoms;
}

int
molecule_class_info_t::delete_hydrogens(){  // return status of atoms deleted (0 -> none deleted).

   int n_deleted = 0;

   // we make a big list and then delete them all at once.  Deleting
   // them in place (inside the residue block crashes).
   //
   std::vector<mmdb::Atom *> atoms_to_be_deleted;

   if (molecule_has_hydrogens()) {
      for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            mmdb::Residue *residue_p;
            mmdb::Atom *at;
            for (int ires=0; ires<nres; ires++) {
               residue_p = chain_p->GetResidue(ires);
               mmdb::PPAtom residue_atoms = 0;
               int n_residue_atoms;
               residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  at = residue_atoms[iat];
                  std::string ele(at->element);
                  if (ele == " H")
                     atoms_to_be_deleted.push_back(at);
                  if (ele == " D")
                     atoms_to_be_deleted.push_back(at);
               }
            }
         }
      }

      if (atoms_to_be_deleted.size() > 0) {

         make_backup();
         for (unsigned int iat=0; iat<atoms_to_be_deleted.size(); iat++) {
            delete atoms_to_be_deleted[iat];
            atoms_to_be_deleted[iat] = NULL;
         }

         atom_sel.mol->FinishStructEdit();
         atom_sel = make_asc(atom_sel.mol);
         make_bonds_type_checked(__FUNCTION__);
         have_unsaved_changes_flag = 1;
         // unlikely to be necessary:
         trim_atom_label_table();
         update_symmetry();
      }
   }
   return atoms_to_be_deleted.size();
}

int
molecule_class_info_t::delete_waters() {

   std::vector<mmdb::Atom *> waters_to_be_deleted;
  for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         mmdb::Residue *residue_p;
         mmdb::Atom *at;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            std::string res_name(residue_p->GetResName());
            if (res_name == "HOH") {
               mmdb::PPAtom residue_atoms = 0;
               int n_residue_atoms;
               residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  at = residue_atoms[iat];
                  waters_to_be_deleted.push_back(at);
               }
            }
         }
      }
   }
   for (unsigned int iat=0; iat<waters_to_be_deleted.size(); iat++) {
      delete waters_to_be_deleted[iat];
      waters_to_be_deleted[iat] = NULL;
   }

   if (waters_to_be_deleted.size() > 0) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
   return waters_to_be_deleted.size();
}



// best fit rotamer stuff.
float
molecule_class_info_t::auto_fit_best_rotamer(int rotamer_search_mode,
                                             int resno,
                                             const std::string &altloc,
                                             const std::string &insertion_code,
                                             const std::string &chain_id,
                                             int imol_map,
                                             int clash_flag, float lowest_prob,
                                             const coot::protein_geometry &pg) {

   // 20090714 We decide here if we go into auto_fit_best_rotamer
   // (conventional mode with rigid body fitting) or backrub rotamers
   //
   float backrub_reso_limit = 2.9; // resolutions worse than this go
                                   // into backrub mode if automatic
                                   // is selected.

   bool do_backrub = false;
   if (rotamer_search_mode == ROTAMERSEARCHLOWRES)
      do_backrub = true;

   if (rotamer_search_mode == ROTAMERSEARCHAUTOMATIC) {
      if (graphics_info_t::is_valid_map_molecule(imol_map)) {
         float r = graphics_info_t::molecules[imol_map].data_resolution();
         if (r > backrub_reso_limit)
            do_backrub = 1;
      }
   }

   if (do_backrub) {
      std::pair<bool,float> br_score = backrub_rotamer(chain_id, resno, insertion_code, altloc, pg);
      if (br_score.first)
         return br_score.second;
      else
         return auto_fit_best_rotamer(resno, altloc, insertion_code, chain_id, imol_map,
                                      clash_flag, lowest_prob, pg);
   } else {
      return auto_fit_best_rotamer(resno, altloc, insertion_code, chain_id, imol_map,
                                   clash_flag, lowest_prob, pg);
   }
}


// best fit rotamer
float
molecule_class_info_t::auto_fit_best_rotamer(int resno,
                                             const std::string &altloc,
                                             const std::string &insertion_code,
                                             const std::string &chain_id,
                                             int imol_map,
                                             int clash_flag, float lowest_prob,
                                             const coot::protein_geometry &pg) {

   // First we will get the mmdb::Residue for the residue that we are
   // trying to fit.
   //
   // For each rotamer of that residue: (probabilities.size() is the
   // number of rotamers for a residue of that type)
   //
   //    Rigid body refine rotamer
   //
   //    score by density fit to map imol_map
   //
   //

   // First check that imol_map has a map.
   bool have_map_flag = true;
   if (imol_map < 0)
      have_map_flag = false;
   if (! graphics_info_t::is_valid_map_molecule(imol_map))
      have_map_flag = false;

   float f = -99.9;
   float clash_score_limit = 20.0; // Rotamers must have a clash score
                                   // less than this (if clashes are
                                   // tested).


   mmdb::Residue *res = get_residue(std::string(chain_id), resno, std::string(insertion_code));

   if (res) {
      if (false) {
         std::cout << " ==== fitting residue " << res->GetSeqNum() << res->GetInsCode()
         << " of chain " << chain_id;
         if (have_map_flag)
         std::cout << " to map number " << imol_map << " ======" << std::endl;
         else
            std::cout << " without a map =====" << std::endl;
      }

      if (coot::util::residue_has_hydrogens_p(res))
         clash_score_limit = 500; // be more generous... lots of hydrogen contacts

      // std::cout << "DEBUG:: found residue" << std::endl;
      std::string res_type(res->name);
      // use the old coords function - update this at some stage
      bool embed_in_new_chain_flag = false;
      mmdb::Residue *copied_res = coot::deep_copy_this_residue_old_style(res, altloc, 0, atom_sel.UDDAtomIndexHandle, embed_in_new_chain_flag);

      if (!copied_res) {
         std::cout << "WARNING:: residue copied - no atoms" << std::endl;
      } else {
         coot::richardson_rotamer d(copied_res, altloc, atom_sel.mol, lowest_prob, 0);
         std::vector<float> probabilities = d.probabilities();
         if (probabilities.size() == 0) {
            std::cout << "WARNING:: no rotamers probabilities for residue type "
                      << res_type << std::endl;
         } else {
            mmdb::Residue *rotamer_res;
            double best_score = -99.9;
            std::vector<mmdb::Atom *> clashing_waters_for_best_score;
            coot::minimol::molecule best_rotamer_mol;

            std::string monomer_type = res->GetResName();
            std::pair<short int, coot::dictionary_residue_restraints_t> p =
               pg.get_monomer_restraints(monomer_type, imol_no);

            std::pair<float, std::vector<mmdb::Atom *> > cs;
            std::vector<mmdb::Atom *> clashing_waters;

            if (p.first) {
               coot::dictionary_residue_restraints_t rest = p.second;

               if (have_map_flag) {
                  for (unsigned int i=0; i<probabilities.size(); i++) {
                     // std::cout << "--- Rotamer number " << i << " ------"  << std::endl;
                     rotamer_res = d.GetResidue(rest, i); // does a deep copy, needs deleting

                     // first make a minimol molecule for the residue so that we
                     // can install it into lig.
                     //
                     coot::minimol::residue residue_res(rotamer_res);
                     coot::minimol::molecule residue_mol;
                     coot::minimol::fragment frag;

                     int ifrag = residue_mol.fragment_for_chain(chain_id);
                     try {
                        residue_mol[ifrag].addresidue(residue_res, 0);

                        // now do "ligand" fitting setup and run:
                        coot::ligand lig;
                        lig.import_map_from(graphics_info_t::molecules[imol_map].xmap);
                        lig.set_acceptable_fit_fraction(0.5);  // at least half of the atoms
                        // have to be fitted into
                        // positive density, otherwise
                        // the fit failed, and we leave
                        // the atom positions as they
                        // are (presumably in an even
                        // worse position?)
                        lig.install_ligand(residue_mol);
                        lig.set_dont_write_solutions();
                        lig.find_centre_by_ligand(0);
                        lig.set_dont_test_rotations();
                        lig.fit_ligands_to_clusters(1);
                        // so we have the solution from lig, what was its score?
                        //
                        unsigned int iclust = 0;
                        unsigned int isol   = 0;
                        coot::ligand_score_card score_card = lig.get_solution_score(0, isol);
                        coot::minimol::molecule moved_mol  = lig.get_solution(isol, iclust);

                        float clash_score = 0.0;
                        if (clash_flag) {
                           bool score_H_atoms = false;
                           bool water_interaction_mode = 0; // ignore waters in scoring, return close waters for deletion
                           cs = get_clash_score(moved_mol, score_H_atoms, water_interaction_mode); // clash on atom_sel.mol
                           clash_score = cs.first;
                        }
                        if (clash_score < clash_score_limit) { // This value may need to be exported
                                                               // to the user interface.
                           if (score_card.get_score() > best_score) {
                              best_score = score_card.get_score();
                              // 20081120 best_rotamer_mol loses the insertion
                              // code for the residue.  Must fix.
                              best_rotamer_mol = moved_mol;
                              clashing_waters_for_best_score = cs.second;
                           }
                        }
                     }

                     catch (const std::runtime_error &rte) {
                        std::cout << "ERROR:: auto_fit_best_rotamer() " << rte.what() << std::endl;
                     }

                     if (rotamer_res) {
                        // implicitly delete rotamer_res too
                        delete rotamer_res->chain;
                     }
                     rotamer_res = NULL;

                  }
               } else {

                  // --- we don't have a map --- like KD suggests:

                  float clash_score;
                  best_score = -99.9; // clash score are better when lower,
                                      // but we should stay consistent with
                                      // the map based scoring which is the
                                      // other way round.
                  for (unsigned int i=0; i<probabilities.size(); i++) {
                     // std::cout << "--- Rotamer number " << i << " ------"  << std::endl;
                     // std::cout << "Getting rotamered residue... " << std::endl;
                     rotamer_res = d.GetResidue(rest, i); // does a deep copy, needs deleting
                     // std::cout << "Got rotamered residue... " << std::endl;
                     coot::minimol::residue  residue_res(rotamer_res);
                     coot::minimol::molecule residue_mol;
                     int ifrag = residue_mol.fragment_for_chain(chain_id);
                     try {
                        residue_mol[ifrag].addresidue(residue_res, 0);
                     }
                     catch (const std::runtime_error &rte) {
                        std::cout << "ERROR:: auto_fit_best_rotamer() 2 " << rte.what() << std::endl;
                     }
                     coot::minimol::molecule moved_mol = residue_mol;
                     // std::cout << "Getting clash score... " << std::endl;
                     bool score_H_atoms = false;
                     int water_interaction_mode = 1; // 2019 default - include waters in clash score
                     water_interaction_mode = 0; // don't clash waters
                     cs = get_clash_score(moved_mol, score_H_atoms, water_interaction_mode); // clash on atom_sel.mol
                     clash_score = -cs.first;
                     // std::cout << "INFO:: clash score: " << clash_score << "\n";
                     if (clash_score > best_score) {
                        best_score = clash_score;
                        best_rotamer_mol = moved_mol;
                        clashing_waters_for_best_score = cs.second;
                     }
                     if (rotamer_res) {
                        // implicitly delete rotamer_res too
                        delete rotamer_res->chain;
                     }
                     rotamer_res = NULL;
                  }
               }
            }

            // now we have tested each rotamer:
            if (best_score > -99.9) {
               // replace current atoms by the atoms in best_rotamer_mol
               //
               // we create a mol and then an asc and then use
               // replace_coords method:
               //
               mmdb::PManager mol = best_rotamer_mol.pcmmdbmanager();
               atom_selection_container_t asc = make_asc(mol);
               bool move_zero_occ = true;
               replace_coords(asc, 1, move_zero_occ); // fix other alt conf occ
               f = best_score;
               std::vector<coot::atom_spec_t> baddie_waters;
               if (clashing_waters_for_best_score.size() > 0) {
                  for (unsigned int ii=0; ii<clashing_waters_for_best_score.size(); ii++) {
                     baddie_waters.push_back(coot::atom_spec_t(clashing_waters_for_best_score[ii]));
                  }
                  delete_atoms(baddie_waters);
               }
            }

            // ALAs and GLY's get to here without entering the loops
            if (probabilities.size() == 0) { // ALAs and GLYs
               f = 0.0;
            }
         }
         delete copied_res->chain; // implicitly delete copied_res too.
      } // non-null copied res

   } else {
      std::cout << "WARNING:: residue not found in molecule" << std::endl;
   }

   return f;
}

// interface from atom picking:
float
molecule_class_info_t::auto_fit_best_rotamer(int rotamer_search_mode,
                                             int atom_index, int imol_map, int clash_flag,
                                             float lowest_probability,
                                             const coot::protein_geometry &pg) {

   int resno = atom_sel.atom_selection[atom_index]->GetSeqNum();
   std::string insertion_code = atom_sel.atom_selection[atom_index]->GetInsCode();
   std::string chain_id = atom_sel.atom_selection[atom_index]->GetChainID();
   std::string altloc(atom_sel.atom_selection[atom_index]->altLoc);

   return auto_fit_best_rotamer(rotamer_search_mode,
                                resno, altloc, insertion_code, chain_id, imol_map,
                                clash_flag, lowest_probability, pg);

}

std::vector<coot::named_rotamer_score>
molecule_class_info_t::score_rotamers(const std::string &chain_id,
                                      int res_no,
                                      const std::string &ins_code,
                                      const std::string &alt_conf,
                                      int clash_flag,
                                      float lowest_probability,
                                      const clipper::Xmap<float> &xmap_in,
                                      const coot::protein_geometry &pg) {

   std::vector<coot::named_rotamer_score> v;
   mmdb::Residue *res = get_residue(std::string(chain_id), res_no, std::string(ins_code));
   if (res) {
      std::string res_type(res->name);
      bool embed_in_new_chain_flag = false; // I presume.
      mmdb::Residue *copied_res = coot::deep_copy_this_residue_old_style(res, alt_conf, 0, atom_sel.UDDAtomIndexHandle, embed_in_new_chain_flag);
      if (!copied_res) {
         std::cout << "WARNING:: residue copied - no atoms" << std::endl;
      } else {
         coot::richardson_rotamer d(copied_res, alt_conf, atom_sel.mol, lowest_probability, 0);
         std::vector<float> probabilities = d.probabilities();
         if (probabilities.size() == 0) {
            std::cout << "WARNING:: no rotamers probabilities for residue type "
                      << res_type << std::endl;
         } else {
            std::pair<short int, coot::dictionary_residue_restraints_t> p =
               pg.get_monomer_restraints(res_type, imol_no);

            if (p.first) {
               coot::dictionary_residue_restraints_t rest = p.second;
               for (unsigned int i=0; i<probabilities.size(); i++) {
                  mmdb::Residue *rotamer_res = d.GetResidue(rest, i); // does a deep copy, deleted
                  std::string rotamer_name = d.rotamer_name(i);
                  float clash_score = -1; // unset

                  coot::rotamer rotamer_for_residue(rotamer_res);
                  coot::rotamer_probability_info_t rpi =
                     rotamer_for_residue.probability_of_this_rotamer();

                  coot::minimol::residue residue_res(rotamer_res);
                  if (clash_flag) {
                     // to get the clash score, we need a minimol molecule
                     coot::minimol::fragment frag(chain_id);
                     frag.addresidue(residue_res, 0);
                     coot::minimol::molecule mm(frag);
                     bool score_H_atoms = false;
                     int water_interaction_mode = 1; // 2019 default
                     water_interaction_mode = 0; // don't clash waters
                     std::pair<float, std::vector<mmdb::Atom *> > cs = get_clash_score(mm, score_H_atoms, water_interaction_mode);
                     clash_score = cs.first;
                  }

                  std::vector<std::pair<std::string, float> > atom_densities =
                     coot::util::score_atoms(residue_res, xmap_in);
                  float rot_prob = rpi.probability;
                  float total_atom_density_score = 0.0;
                  for (unsigned int iat=0; iat<atom_densities.size(); iat++)
                     if (! coot::is_main_chain_or_cb_p(atom_densities[iat].first))
                         total_atom_density_score += atom_densities[iat].second;

                  coot::named_rotamer_score rot(rotamer_name,
                                                rot_prob,
                                                clash_score,
                                                atom_densities,
                                                total_atom_density_score);
                  v.push_back(rot);
                  delete rotamer_res->chain; // strange (perhaps?) - deletes rotamer_res too.
               }
            }
         }
      }
   }
   return v;
}


// internal.
std::pair<bool,float>
molecule_class_info_t::backrub_rotamer(const std::string &chain_id, int res_no,
                                       const std::string &ins_code, const std::string &alt_conf,
                                       const coot::protein_geometry &pg) {

   bool status = 0;
   float score = -1;
   graphics_info_t g;
   int imol_map = g.Imol_Refinement_Map();
   if (imol_map >= 0) {
      if (imol_map < int(graphics_info_t::molecules.size())) {
         if (graphics_info_t::molecules[imol_map].has_xmap() ||
             graphics_info_t::molecules[imol_map].has_nxmap()) {
            mmdb::Residue *res = get_residue(chain_id, res_no, ins_code);
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

                     make_backup();
                     mmdb::Residue *prev_res = coot::util::previous_residue(res);
                     mmdb::Residue *next_res = coot::util::next_residue(res);
                     mmdb::Manager *mol = atom_sel.mol;
                     coot::backrub br(chain_id, res, prev_res, next_res, alt_conf, mol,
                                      &g.molecules[imol_map].xmap); // use a pointer for the map
                     std::pair<coot::minimol::molecule,float> m = br.search(restraints);
                     std::vector<coot::atom_spec_t> baddie_waters = br.waters_for_deletion();
                     score = m.second;
                     status = 1;
                     atom_selection_container_t fragment_asc = make_asc(m.first.pcmmdbmanager());
                     bool mzo = g.refinement_move_atoms_with_zero_occupancy_flag;
                     replace_coords(fragment_asc, 0, mzo);
                     // std::cout << "Debug:: waters for deletion size " << baddie_waters.size() << std::endl;
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
         } else {
            std::cout << "   WARNING:: " << imol_map << " is not a valid map molecule"
                      << std::endl;
         }
      }
   }
   return std::pair<bool,float>(status,score);
}

void
molecule_class_info_t::backrub_rotamer_residue_range(const std::string &chain_id, int resno_start, int resno_end,
                                                     const coot::protein_geometry &pg) {

   for (int resno=resno_start; resno<=resno_end; resno++)
      backrub_rotamer(chain_id, resno, "", "", pg);
}



int
molecule_class_info_t::set_residue_to_rotamer_number(coot::residue_spec_t res_spec,
                                                     const std::string &alt_conf_in,
                                                     int rotamer_number,
                                                     const coot::protein_geometry &pg) {

   int i_done = 0;
   mmdb::Residue *res = get_residue(res_spec);
   if (res) {

      make_backup();
#ifdef USE_DUNBRACK_ROTAMERS
      coot::dunbrack d(res, atom_sel.mol, 0.01, 0);
#else
      std::string alt_conf = ""; // fixme?
      coot::richardson_rotamer d(res, alt_conf, atom_sel.mol, 0.01, 0);
#endif // USE_DUNBRACK_ROTAMERS
      std::string monomer_type = res->GetResName();
      std::pair<short int, coot::dictionary_residue_restraints_t> p =
         pg.get_monomer_restraints(monomer_type, imol_no);

      if (p.first) {
         coot::dictionary_residue_restraints_t rest = p.second;
         mmdb::Residue *moving_res = d.GetResidue(rest, rotamer_number);
         if (moving_res) {
            i_done = set_residue_to_rotamer_move_atoms(res, moving_res);
            delete moving_res; // or moving_res->chain?
         }
      }
   }
   if (! i_done)
      std::cout << "WARNING:: set to rotamer number failed" << std::endl;
   return i_done;
}

int
molecule_class_info_t::set_residue_to_rotamer_name(coot::residue_spec_t res_spec,
                                                   const std::string &alt_conf_in,
                                                   const std::string &rotamer_name,
                                                   const coot::protein_geometry &pg) {

   int status = 0;
   mmdb::Residue *res = get_residue(res_spec);
   if (res) {
      make_backup();
      coot::richardson_rotamer rr(res, alt_conf_in, atom_sel.mol, 0.01, 0);
      std::string monomer_type = res->GetResName();
      std::pair<short int, coot::dictionary_residue_restraints_t> p =
         pg.get_monomer_restraints(monomer_type, imol_no);

      if (p.first) {
         coot::dictionary_residue_restraints_t rest = p.second;
         mmdb::Residue *moving_res = rr.GetResidue(rest, rotamer_name);
         if (moving_res) {
            status = set_residue_to_rotamer_move_atoms(res, moving_res);
            delete moving_res; // or delete moving_res->chain ??
         }
      }
   }
   return status;
}

int
molecule_class_info_t::set_residue_to_rotamer_move_atoms(mmdb::Residue *res, mmdb::Residue *moving_res) {

   int i_done = 0;
   int n_ref_atoms;
   mmdb::PPAtom ref_residue_atoms = 0;
   int n_mov_atoms;
   mmdb::PPAtom mov_residue_atoms= 0;

   res->GetAtomTable(ref_residue_atoms, n_ref_atoms);
   moving_res->GetAtomTable(mov_residue_atoms, n_mov_atoms);

   int n_atoms = 0;
   for (int imov=0; imov<n_mov_atoms; imov++) {
      std::string atom_name_mov(mov_residue_atoms[imov]->name);
      std::string alt_loc_mov(mov_residue_atoms[imov]->altLoc);
      for (int iref=0; iref<n_ref_atoms; iref++) {
         std::string atom_name_ref(ref_residue_atoms[iref]->name);
         std::string alt_loc_ref(ref_residue_atoms[iref]->altLoc);
         if (atom_name_mov == atom_name_ref) {
            if (alt_loc_mov == alt_loc_ref) {
               ref_residue_atoms[iref]->x = mov_residue_atoms[imov]->x;
               ref_residue_atoms[iref]->y = mov_residue_atoms[imov]->y;
               ref_residue_atoms[iref]->z = mov_residue_atoms[imov]->z;
               n_atoms++;
               i_done = 1;
            }
         }
      }
   }

   if (i_done) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
   return i_done;
}

// Return a copy of the pointer to the chain (only).  Return NULL
// on chain with given chain ID not found.
//
mmdb::Chain *
molecule_class_info_t::get_chain(const std::string &chain_id) const {

  mmdb::Chain *r = NULL;

  mmdb::Manager *mol = atom_sel.mol;

  if (mol) {
     int imod = 1;
     mmdb::Model *model_p = mol->GetModel(imod);
     mmdb::Chain *chain_p;
     int n_chains = model_p->GetNumberOfChains();
     for (int ichain=0; ichain<n_chains; ichain++) {
        chain_p = model_p->GetChain(ichain);
        std::string mol_chain_id = chain_p->GetChainID();
        if (chain_id == mol_chain_id) {
           r = chain_p;
           break;
        }
     }
  }
  return r;
}



// Return NULL on residue not found in this molecule.
//
mmdb::Residue *
molecule_class_info_t::get_residue(const std::string &chain_id,
                                   int resno,
                                   const std::string &insertion_code) const {

   mmdb::Residue *res = NULL;
   if (atom_sel.n_selected_atoms > 0) {
      res = coot::util::get_residue(chain_id, resno, insertion_code, atom_sel.mol);
   }
   return res;
}

mmdb::Residue *
molecule_class_info_t::get_residue(const coot::residue_spec_t &residue_spec) const {

   mmdb::Residue *res = get_residue(residue_spec.chain_id,
                                    residue_spec.res_no,
                                    residue_spec.ins_code);
   return res;
}

std::string
molecule_class_info_t::get_residue_name(const coot::residue_spec_t &rs) const {

   std::string rn;
   mmdb::Residue *r = get_residue(rs);
   if (r) {
      rn = r->GetResName();
   }
   return rn;
}

// Useful when we know that the molecule is just one residue
mmdb::Residue *
molecule_class_info_t::get_first_residue() {

   mmdb::Residue *res = 0;

   if (atom_sel.mol) {
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         mmdb::Residue *residue_p;
         mmdb::Atom *at;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            int n_atoms = residue_p->GetNumberOfAtoms();
            if (n_atoms) {
               res = residue_p;
               break;
            }
         }
      }
   }
   return res;
}



// Can return NULL.
mmdb::Residue *
molecule_class_info_t::get_following_residue(const coot::residue_spec_t &rs) const {

   mmdb::Residue *res = NULL;
   if (atom_sel.n_selected_atoms > 0) {
      res = coot::util::get_following_residue(rs, atom_sel.mol);
      if (res) {
         int new_res_res_number = res->GetSeqNum();
         if (new_res_res_number < rs.res_no) {
            res = NULL;
         }
      }
   }
   return res;
}



mmdb::Residue *
molecule_class_info_t::residue_from_external(int resno, const std::string &insertion_code,
                                             const std::string &chain_id) const {

   return get_residue(chain_id, resno, insertion_code);

}

mmdb::Atom *
molecule_class_info_t::get_atom(const std::string &go_to_residue_string,
                                const coot::atom_spec_t &active_atom_spec,
                                const coot::Cartesian &pt) const {

   mmdb::Atom *at = NULL;
   coot::goto_residue_string_info_t si(go_to_residue_string, atom_sel.mol);
   if (si.chain_id_is_set) {
      if (si.res_no_is_set) {
         mmdb::Residue *res_p = get_residue(si.chain_id, si.res_no, "");
         if (res_p) {
            mmdb::Atom *int_at = intelligent_this_residue_mmdb_atom(res_p);
            if (int_at) {
               at = int_at;
            }
         }
      } else {
         // the closest atom in the si.chain_id
         coot::at_dist_info_t cl_at = closest_atom(pt, 1, si.chain_id, 1);
         if (cl_at.atom)
            at = cl_at.atom;
      }
   } else {
      // use the chain_id from the active_atom_spec.chain then
      if (si.res_no_is_set) {
         mmdb::Residue *res_p = get_residue(active_atom_spec.chain_id, si.res_no, "");
         if (res_p) {
            mmdb::Atom *int_at = intelligent_this_residue_mmdb_atom(res_p);
            if (int_at) {
               at = int_at;
            }
         }
      }
   }
   return at;
}

// simple version of above, don't be intelligent.
//
mmdb::Atom *
molecule_class_info_t::get_atom(const coot::atom_spec_t &atom_spec) const {

   mmdb::Residue *res = get_residue(atom_spec.chain_id, atom_spec.res_no, atom_spec.ins_code);
   mmdb::Atom *at = NULL;

   auto get_atom_from_residue = [res, atom_spec] (const std::string &test_atom_name) {
                                   mmdb::Atom *at = 0;
                                   mmdb::PPAtom residue_atoms = 0;
                                   int nResidueAtoms = 0;
                                   res->GetAtomTable(residue_atoms, nResidueAtoms);
                                   for (int iat=0; iat<nResidueAtoms; iat++) {
                                      mmdb::Atom *test_at = residue_atoms[iat];
                                      std::string at_name(test_at->name);
                                      if (test_atom_name == at_name) {
                                         std::string at_alt_conf(test_at->altLoc);
                                         if (atom_spec.alt_conf == at_alt_conf) {
                                            at = test_at;
                                            break;
                                         }
                                      }
                                   }
                                   return at;
                                };

   if (res) {
      mmdb::PPAtom residue_atoms = 0;
      int nResidueAtoms = 0;
      res->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int iat=0; iat<nResidueAtoms; iat++) {
         mmdb::Atom *test_at = residue_atoms[iat];
         std::string at_name = test_at->name;
         std::string at_alt_conf = test_at->altLoc;
         if (atom_spec.atom_name == at_name) {
            if (atom_spec.alt_conf == at_alt_conf) {
               at = test_at;
               break;
            }
         }
         const std::size_t asnl = atom_spec.atom_name.length();
         if (asnl != 4) {
            // perhaps we were give an atom name with no spaces?
            if (asnl == 1) {
               std::string test_atom_name = " " + atom_spec.atom_name + "  ";
               at = get_atom_from_residue(test_atom_name);
               if (! at) {
                  test_atom_name = atom_spec.atom_name + "   ";
                  at = get_atom_from_residue(test_atom_name);
               }
            }
            if (asnl == 2) {
               std::string test_atom_name = " " + atom_spec.atom_name + " ";
               at = get_atom_from_residue(test_atom_name);
               if (! at) {
                  test_atom_name = atom_spec.atom_name + "  ";
                  at = get_atom_from_residue(test_atom_name);
               }
            }
            if (asnl == 3) {
               std::string test_atom_name = " " + atom_spec.atom_name;
               at = get_atom_from_residue(test_atom_name);               
            }
         }
      }
   }
   return at;
}

mmdb::Atom *
molecule_class_info_t::get_atom(int idx) const {

   mmdb::Atom *r = NULL;
   if (idx < atom_sel.n_selected_atoms)
      r = atom_sel.atom_selection[idx];
   return r;
}

mmdb::Atom *
molecule_class_info_t::get_atom(const pick_info &pi) const {

   mmdb::Atom *at = 0;
   if (pi.success == GL_TRUE)
      if (pi.atom_index < atom_sel.n_selected_atoms)
         at = atom_sel.atom_selection[pi.atom_index];
   return at;

}



// This should check that if "a" is typed, then set "a" as the
// chain_id if it exists, else convert to "A" (if that exists).
//
// When there is a space, then the first string is the chain id and second string is the res_no.
//
coot::goto_residue_string_info_t::goto_residue_string_info_t(const std::string &goto_residue_string,
                                                             mmdb::Manager *mol) {

   res_no_is_set    = false;
   chain_id_is_set  = false;
   res_no = mmdb::MinInt4;
   chain_id = "";

   std::vector<std::string> bits = coot::util::split_string_no_blanks(goto_residue_string);


   if (bits.size() == 1) {

      std::vector<std::string> chain_ids;
      if (mol) {
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<nchains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            chain_ids.push_back(chain_p->GetChainID());
         }
      }

      std::string::size_type l = goto_residue_string.length();
      std::string number_string = "";
      std::string chain_id_string = "";
      for (std::string::size_type i=0; i<l; i++) {
         if (coot::util::is_number(goto_residue_string[i])) {
            number_string += goto_residue_string[i];
            res_no_is_set = true;
         }
         if (coot::util::is_letter(goto_residue_string[i])) {
            chain_id_string += goto_residue_string[i];
            chain_id_is_set = true;
         }
      }

      if (res_no_is_set) {
         res_no = atoi(number_string.c_str());
      }

      if (chain_id_is_set) {
         chain_id = chain_id_string;
      }
   } else {

      // The Chris Oubridge case, the first set of characters might actually be a chain id.

      if (bits.size() == 2) {

         if (mol) {
            int imod = 1;
            mmdb::Model *model_p = mol->GetModel(imod);
            mmdb::Chain *chain_p;
            // run over chains of the existing mol
            int nchains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<nchains; ichain++) {
               chain_p = model_p->GetChain(ichain);
               std::string this_chain_id = chain_p->GetChainID();
               if (this_chain_id == bits[0]) {
                  try {
                     res_no = coot::util::string_to_int(bits[1]);
                     res_no_is_set = true;
                     chain_id = this_chain_id;
                     chain_id_is_set = true;
                  }
                  catch (const std::runtime_error &rte) {

                  }
               }
            }
         }
      }
   }
}



// pair.first is the status, 0 for bad, 1 for good.
//
std::pair<short int, mmdb::Atom *>
molecule_class_info_t::baton_build_delete_last_residue() {

   // this is the baton molecule, we don't get here if we aren't.
   //
   // How do we find the last residue?  Simply it is the last chain in
   // the list of chains and the last residue in the list of residues.
   //
   std::pair<short int, mmdb::Atom *> r;
   r.first = 0; // status

   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      mmdb::Chain *chain_p;
      int n_chains = atom_sel.mol->GetNumberOfChains(1);
      chain_p = model_p->GetChain(n_chains-1);

      int nres = chain_p->GetNumberOfResidues();
      mmdb::Residue *residue_p = chain_p->GetResidue(nres-1); // is a serial number

      int iseqno = residue_p->GetSeqNum();
      mmdb::pstr inscode = residue_p->GetInsCode();
      chain_p->DeleteResidue(iseqno, inscode);

      // Kevin's Baton Build crash: Argh.  I have no idea why.
      // atom_sel.SelectionHandle seems good and is generated by an
      // atom selection in make asc (when the latest Ca had been
      // placed).
      //
//       if (atom_sel.SelectionHandle > 0) {
//          atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle); // prefered.
//       }

      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_ca_bonds(2.4, 4.7);

      if (atom_sel.n_selected_atoms > 0 ) { // atom sel was modified.
         residue_p = chain_p->GetResidue(nres-2);
         mmdb::Atom *atom_p = residue_p->GetAtom(" CA ");
         if (atom_p) {
            r.first = 1; // OK
            r.second = atom_p;
         }
      }
   }
   return r;
}

std::pair<short int, clipper::Coord_grid>
molecule_class_info_t::search_for_skeleton_near(const coot::Cartesian &pos) const {

   // std::pair<short int, clipper::Coord_grid> r;
   clipper::Coord_orth co(pos.x(), pos.y(), pos.z());

   coot::CalphaBuild cab;
   std::pair<short int, clipper::Coord_grid> r =
      cab.search_for_skeleton_near(co, xskel_cowtan, skeleton_treenodemap);

   return r;
}

short int
molecule_class_info_t::add_OXT_to_residue(int reso, const std::string &insertion_code,
                                          const std::string &chain_id,
                                          coot::protein_geometry *geom_p) {

   mmdb::Residue *residue = get_residue(chain_id, reso, insertion_code);
   return add_OXT_to_residue(residue, geom_p);  // checks for null residue

}

short int
molecule_class_info_t::add_OXT_to_residue(mmdb::Residue *residue,
                                          coot::protein_geometry *geom_p) {

   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   short int istatus = 0; // fail
   if (!residue) {
      std::cout << "WARNING: NULL residue, no atom added." << std::endl;
   } else {
      residue->GetAtomTable(residue_atoms, nResidueAtoms);
      if (nResidueAtoms == 0) {
         std::cout << "WARNING: no atoms in this residue" << std::endl;
      } else {

         bool has_already = residue_has_oxt_p(residue);
         if (has_already) {
            std::cout << "WARNING:: This residue already has an OXT - aborting\n";
            return 0;
         }

         mmdb::Atom *n_atom = NULL;
         mmdb::Atom *c_atom = NULL;
         mmdb::Atom *ca_atom = NULL;
         mmdb::Atom *o_atom = NULL;

         clipper::Coord_orth co;
         clipper::Coord_orth c_atom_co;
         clipper::Coord_orth o_atom_co;
         clipper::Coord_orth n_atom_co;
         clipper::Coord_orth ca_atom_co;
         int n_found_atoms = 0;

         mmdb::Atom *atom;
         atom = residue->GetAtom(" N  ");
         if (atom) {
            n_found_atoms++;
            n_atom = atom;
            n_atom_co = clipper::Coord_orth(atom->x, atom->y, atom->z);
         }
         atom = residue->GetAtom(" O  ");
         if (atom) {
            n_found_atoms++;
            o_atom = atom;
            o_atom_co = clipper::Coord_orth(atom->x, atom->y, atom->z);
         }
         atom = residue->GetAtom(" CA ");
         if (atom) {
            n_found_atoms++;
            ca_atom = atom;
            ca_atom_co = clipper::Coord_orth(atom->x, atom->y, atom->z);
         }
         atom = residue->GetAtom(" C  ");
         if (atom) {
            n_found_atoms++;
            c_atom = atom;
            c_atom_co = clipper::Coord_orth(atom->x, atom->y, atom->z);
         }

         if (! (n_atom && c_atom && ca_atom && o_atom) ) {
            std::cout << "WARNING:: Not all reference atoms found in residue."
                      << std::endl;
            std::cout << "          No atom fitted." << std::endl;
            std::string m("WARNING:: Not all reference atoms found in residue\n");
            m += "          No OXT atom fitted.";
            GtkWidget *w = graphics_info_t::wrapped_nothing_bad_dialog(m);
            gtk_widget_set_visible(w, TRUE);
         } else {
            make_backup();
            double torsion = clipper::Coord_orth::torsion(n_atom_co, ca_atom_co, c_atom_co, o_atom_co);
            double angle = clipper::Util::d2rad(120.8);
            clipper::Coord_orth new_oxt(n_atom_co, ca_atom_co, c_atom_co,
                                        1.231, angle, torsion + M_PI);
            mmdb::Atom *new_oxt_atom = new mmdb::Atom;
            new_oxt_atom->SetCoordinates(new_oxt.x(),
                                         new_oxt.y(),
                                         new_oxt.z(), 1.0,
                                         graphics_info_t::default_new_atoms_b_factor);
            new_oxt_atom->SetAtomName(" OXT");
            new_oxt_atom->SetElementName(" O");
            if (coot::util::residue_has_hetatms(residue))
               new_oxt_atom->Het = true;
            residue->AddAtom(new_oxt_atom);

            atom_sel.mol->FinishStructEdit();
            atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);

            // Now handle the TER atom.
            //
            residue->GetAtomTable(residue_atoms, nResidueAtoms); // reset the pointers after FinishStructEdit()

            mmdb::Atom *ter_atom = NULL;
            int ter_index = -1;
            for (int iat=0; iat<nResidueAtoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               if (at) {
                  if (at->isTer()) {
                     ter_atom = at;
                     ter_index = iat;
                  }
               } else {
                  std::cout << "ERROR:: trapped null atom in add_OXT_to_residue() " << std::endl;
               }
            }
            if (ter_atom) {
               residue->DeleteAtom(ter_index);
               // create a new TER atom
               mmdb::Atom *new_ter_atom = new mmdb::Atom;
               new_ter_atom->Copy(new_oxt_atom);
               new_ter_atom->MakeTer();
               residue->AddAtom(new_ter_atom);
            }
            atom_sel.mol->FinishStructEdit();

            atom_sel = make_asc(atom_sel.mol);
            have_unsaved_changes_flag = 1;
            std::set<int> dummy;
            makebonds(geom_p, dummy); // not type checked, so that we can see the atom.
            istatus = 1;
            std::cout << "Added OXT at " << new_oxt_atom << std::endl;
         }
      }
   }
   return istatus;
}

// used by above.  Don't add if returns true.
bool
molecule_class_info_t::residue_has_oxt_p(mmdb::Residue *residue) const {

   bool r = 0;

   if (residue) {
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      residue->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
         std::string atname(residue_atoms[i]->name);
         if (atname == " OXT") {
            r = 1;
            break;
         }
      }
   }
   return r;
}


// return -1 on a problem.
int
molecule_class_info_t::residue_serial_number(const std::string &chain_id,
                                             int resno, const std::string &insertion_code) const {

   int iserial = -1;
   mmdb::Residue *res = get_residue(chain_id, resno, insertion_code);
   if (res) {
      // std::cout << "DEBUG:: residue_serial_number residue " << resno << " found " << std::endl;
      iserial = res->index;
      if (iserial == -1) {
         coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
         iserial = res->index;
      }
      if (iserial == -1) {
         std::cout << "WARNING:: residue_number_serial() returns -1 for " << chain_id << " " << resno
                   << " \"" << insertion_code << "\"" << std::endl;
      }
   } else {
      std::cout << "WARNING:: residue" << resno << " " << insertion_code
                << " " << chain_id << " not found" << std::endl;
   }
   return iserial;
}

std::string
molecule_class_info_t::res_name_from_serial_number(std::string chain_id, unsigned int serial_number) const {

   std::string r;
   if (atom_sel.mol) {
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
         mmdb::Chain *chain_p = atom_sel.mol->GetChain(1,ichain);
         std::string mol_chain_id(chain_p->GetChainID());
         if (mol_chain_id == std::string(chain_id)) {
            unsigned int nres = chain_p->GetNumberOfResidues();
            if (serial_number < nres) {
               int ch_n_res;
               mmdb::PResidue *residues;
               chain_p->GetResidueTable(residues, ch_n_res);
               mmdb::Residue *this_res = residues[serial_number];
               r = this_res->GetResName();
            }
         }
      }
   }
   return r;
}



void
molecule_class_info_t::apply_atom_edit(const coot::select_atom_info &sai) {

   // first select the atom
   mmdb::PPAtom SelAtoms = NULL;
   int nSelAtoms;
   int SelHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(SelHnd, 0, (char *) sai.chain_id.c_str(),
                             sai.residue_number, (char *) sai.insertion_code.c_str(),
                             sai.residue_number, (char *) sai.insertion_code.c_str(),
                             "*", // residue name
                             (char *) sai.atom_name.c_str(),
                             "*", // elements
                             (char *) sai.altconf.c_str()); // alt locs

   atom_sel.mol->GetSelIndex(SelHnd, SelAtoms, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "Sorry. Could not find " << sai.atom_name << ","
                << sai.altconf << "/"
                << sai.residue_number << sai.insertion_code << "/" << sai.chain_id
                << " in this molecule: ("
                <<  imol_no << ") " << name_ << std::endl;
   } else {

      if (nSelAtoms > 1) {
         std::cout << "Unexepected condition in apply_atom_edit: many atoms! "
                   << nSelAtoms << std::endl;
      } else {
         mmdb::Atom *atom_p = SelAtoms[0];

         if (sai.has_b_factor_edit())
            atom_p->tempFactor = sai.b_factor;
         if (sai.has_occ_edit())
            atom_p->occupancy = sai.occ;
      }
   }
}

void
molecule_class_info_t::apply_atom_edits(const std::vector<coot::select_atom_info> &saiv) {

   std::cout << "DEBUG:: in mci::apply_atom_edits() " << saiv.size() << std::endl;

   bool made_edit = false;
   make_backup();

   for (unsigned int i=0; i<saiv.size(); i++) {
      // std::cout << "mci::apply_atom_edits() " << i << std::endl;
      mmdb::Atom *at = saiv[i].get_atom(atom_sel.mol);
      if (at) {
         // std::cout << "mci::apply_atom_edits() B " << i << std::endl;
         if (saiv[i].has_b_factor_edit()) {
            // std::cout << "mci::apply_atom_edits() c " << i << std::endl;
            at->tempFactor = saiv[i].b_factor;
            made_edit = 1;
         }
         if (saiv[i].has_occ_edit()) {
            // std::cout << "mci::apply_atom_edits() d " << i << std::endl;
            at->occupancy = saiv[i].occ;
            made_edit = 1;
         }
         if (saiv[i].has_altloc_edit()) {
            // std::cout << "mci::apply_atom_edits() e " << i << std::endl;
            // mmmdb limit is char altloc[20];
            // strncpy(at->altLoc, saiv[i].altloc_new.c_str(), 2);
            // strncpy() writes n bytes no matter what - that is not what I want.
            if (saiv[i].altloc_new.length() < 20)
               strcpy(at->altLoc, saiv[i].altloc_new.c_str());
            made_edit = 1;
         }
      }
   }

   if (made_edit) {
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
}



// For edit phi psi, tinker with the passed (moving atoms) asc
// Return status, -1 on error (unable to change).
short int
molecule_class_info_t::residue_edit_phi_psi(atom_selection_container_t residue_asc,
                                            int atom_index, double phi, double psi) {

   // Don't forget we do this kind of setting positions by torsion value in
   // ligand/residue_by_phi_psi (but recall that we use minimol there).
   //

   mmdb::Residue *this_res = atom_sel.atom_selection[atom_index]->residue;

   // get the previous C and next N atoms.  If we can't find them then
   // we don't do any changes to the atom positions and return -1.

   // First get the residues:
   mmdb::Residue *next_res;
   mmdb::Residue *prev_res;

   int selHnd1 = atom_sel.mol->NewSelection();
   int selHnd2 = atom_sel.mol->NewSelection();
   mmdb::PPResidue SelResidues;
   int nSelResdues;

   atom_sel.mol->Select (selHnd1, mmdb::STYPE_RESIDUE, 0,
                         this_res->GetChainID(),
                         this_res->GetSeqNum() - 1, this_res->GetInsCode(),
                         this_res->GetSeqNum() - 1, this_res->GetInsCode(),
                         "*",
                         "*",
                         "*",
                         "*",
                         mmdb::SKEY_NEW
                         );
   atom_sel.mol->GetSelIndex(selHnd1, SelResidues, nSelResdues);
   if (nSelResdues < 1) {
      std::cout << "Can't find previous residue" << std::endl;
      return -1;
   } else {
      prev_res = SelResidues[0];
   }
   atom_sel.mol->Select (selHnd2, mmdb::STYPE_RESIDUE, 1,
                         this_res->GetChainID(),
                         this_res->GetSeqNum() + 1, this_res->GetInsCode(),
                         this_res->GetSeqNum() + 1, this_res->GetInsCode(),
                         "*",
                         "*",
                         "*",
                         "*",
                         mmdb::SKEY_NEW
                         );
   atom_sel.mol->GetSelIndex(selHnd2, SelResidues, nSelResdues);
   if (nSelResdues < 1) {
      std::cout << "Can't find next residue" << std::endl;
      return -1;
   } else {
      next_res = SelResidues[0];
   }

   // So we have both residues if we have got here:
   //
   // Now to get the relavent atoms (recall that we need Next Ca
   // because the O of this residue lies in the peptide plane of the
   // next residue (and that plane is defined by (amongst other
   // things) the Ca position of the next residue).

   mmdb::Atom *next_N  = next_res->GetAtom(" N  ");
   mmdb::Atom *next_Ca = next_res->GetAtom(" CA ");
   mmdb::Atom *prev_C  = prev_res->GetAtom(" C  ");

   if (next_N == NULL) {
      std::cout << "Can't find N in previous residue\n";
      return -1;
   }
   if (next_Ca == NULL) {
      std::cout << "Can't find CA in next residue\n";
      return -1;
   }
   if (prev_C == NULL) {
      std::cout << "Can't find C in next residue\n";
      return -1;
   }

   // So we have next_N and prev_C if we have got here.

   // delete the selections
   //
   atom_sel.mol->DeleteSelection(selHnd1);
   atom_sel.mol->DeleteSelection(selHnd2);


   //
   mmdb::Atom *this_C  = this_res->GetAtom(" C  ");
   // mmdb::Atom *this_O  = this_res->GetAtom(" O  ");
   mmdb::Atom *this_N  = this_res->GetAtom(" N  ");
   mmdb::Atom *this_Ca = this_res->GetAtom(" CA ");

   if (this_C == NULL) {
      std::cout << "Can't find C in this residue\n";
      return -1;
   }
   if (this_C == NULL) {
      std::cout << "Can't find O in this residue\n";
      return -1;
   }
   if (this_N == NULL) {
      std::cout << "Can't find N in this residue\n";
      return -1;
   }
   if (this_Ca == NULL) {
      std::cout << "Can't find Ca in this residue\n";
      return -1;
   }

   // OK, we we have all the atoms needed to calculate phi and psi if
   // we have got here.
   //
   // Except, we don't want to calculate them, we want to set our
   // atoms.

   // We want to set C (then O) by phi and N by psi.
   //
   // Convert to Coord_orth's

   clipper::Coord_orth cthis_C_orig  = to_coord_orth(this_C );
   // clipper::Coord_orth cthis_N  = to_coord_orth(this_N );
   clipper::Coord_orth cthis_Ca = to_coord_orth(this_Ca);
   clipper::Coord_orth cprev_C  = to_coord_orth(prev_C );
   clipper::Coord_orth cnext_N  = to_coord_orth(next_N );
   clipper::Coord_orth cnext_Ca = to_coord_orth(next_Ca);


   double angle, torsion;

   // we build up from the C terminal end (note that the new positions
   // for C and N are intertependent (because phi and psi depend on
   // the positions of the backbone atoms in this residue).
   //
   //

   // N
   angle   = clipper::Util::d2rad(111.200); // N-Ca-C
   torsion = clipper::Util::d2rad(psi);
   clipper::Coord_orth cthis_N(cnext_N, cthis_C_orig, cthis_Ca,
                               1.458, angle,torsion);

   angle   = clipper::Util::d2rad(111.200);  // N-CA-C
   torsion = clipper::Util::d2rad(phi);
   clipper::Coord_orth cthis_C(cprev_C, cthis_N, cthis_Ca,
                               1.525, angle, torsion);  // CA-C bond

   // O
   angle   = clipper::Util::d2rad(123.0); // N-C-O
   torsion = clipper::Util::d2rad(0.0);
   clipper::Coord_orth cthis_O(cnext_Ca, cnext_N, cthis_C,
                               1.231, angle, torsion);


   // Now apply those changes to the residue_asc (the (white) moving atoms):
   //
   int selHnd_moving = residue_asc.mol->NewSelection();

   residue_asc.mol->Select (selHnd_moving, mmdb::STYPE_RESIDUE, 0,
                            this_res->GetChainID(),
                            this_res->GetSeqNum(), this_res->GetInsCode(),
                            this_res->GetSeqNum(), this_res->GetInsCode(),
                            "*",
                            "*",
                            "*",
                            "*",
                            mmdb::SKEY_NEW
                            );
   residue_asc.mol->GetSelIndex(selHnd_moving, SelResidues, nSelResdues);

   if (nSelResdues != 1) {
      std::cout<< "Can't find this moving residue (bizarrely enough)" << std::endl;
      return -1;
   } else {

      mmdb::Atom *moving_N = SelResidues[0]->GetAtom(" N  ");
      mmdb::Atom *moving_C = SelResidues[0]->GetAtom(" C  ");
      mmdb::Atom *moving_O = SelResidues[0]->GetAtom(" O  ");


      moving_N->x = cthis_N.x();
      moving_N->y = cthis_N.y();
      moving_N->z = cthis_N.z();

      moving_C->x = cthis_C.x();
      moving_C->y = cthis_C.y();
      moving_C->z = cthis_C.z();

      moving_O->x = cthis_O.x();
      moving_O->y = cthis_O.y();
      moving_O->z = cthis_O.z();

      residue_asc.mol->DeleteSelection(selHnd_moving);
   }
   return 1;
}

clipper::Coord_orth
molecule_class_info_t::to_coord_orth(mmdb::Atom *atom) const {

   return clipper::Coord_orth(atom->x, atom->y, atom->z);
}

// Return the residue thats for moving_atoms_asc as a molecule.
//
atom_selection_container_t
molecule_class_info_t::edit_residue_pull_residue(int atom_index,
                                                 short int whole_residue_flag) {

   mmdb::Residue *res = NULL;
   // short int found_res = 0;
   atom_selection_container_t r;
   r.n_selected_atoms = 0; // added 20050402 for Wall compilation
   r.atom_selection = 0;
   r.mol = 0;
   r.read_success = 0;

   res = atom_sel.atom_selection[atom_index]->residue;
   std::string altconf(atom_sel.atom_selection[atom_index]->altLoc);

   if (res) {

      bool embed_in_new_chain_flag = false;
      mmdb::Residue *ret_res = coot::deep_copy_this_residue_old_style(res, altconf,
                                                       whole_residue_flag,
                                                       atom_sel.UDDAtomIndexHandle,
                                                       embed_in_new_chain_flag);
      if (ret_res) {

         mmdb::Manager *MMDBManager = new mmdb::Manager;
         mmdb::Model *model_p = new mmdb::Model;
         mmdb::Chain *chain_p = new mmdb::Chain;

         chain_p->SetChainID(res->GetChainID());

         model_p->AddChain(chain_p);
         MMDBManager->AddModel(model_p);
         chain_p->AddResidue(ret_res);

         atom_selection_container_t r = make_asc(MMDBManager);

         return r;
      }
   }
   return r;
}


// Return pair.first < -200 for error.
//
std::pair<double, double>
molecule_class_info_t::get_phi_psi(int atom_index) const {

   std::pair<double, double> pp;
   pp.first = -999.9; // magic number for bad return values


   if (atom_sel.n_selected_atoms > 0) {
      if (atom_index < atom_sel.n_selected_atoms) {
         mmdb::Residue *res = atom_sel.atom_selection[atom_index]->residue;
         int this_resno = res->GetSeqNum();

         mmdb::PPResidue SelResidues;
         int nSelResdues;

         int selHnd = atom_sel.mol->NewSelection();
         atom_sel.mol->Select ( selHnd,mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
                                res->GetChainID(), // Chain(s)
                                res->GetSeqNum() - 1, "*",  // starting res
                                res->GetSeqNum() + 1, "*",  // ending res
                                "*",  // residue name
                                "*",  // Residue must contain this atom name?
                                "*",  // Residue must contain this Element?
                                "*",  // altLocs
                                mmdb::SKEY_NEW // selection key
                                );
         atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResdues);

         if (nSelResdues == 3) {

            clipper::Coord_orth prev_C;
            clipper::Coord_orth this_N;
            clipper::Coord_orth this_Ca;
            clipper::Coord_orth this_C;
            clipper::Coord_orth next_N;
            int n_atoms = 0;

            for (int ires=0; ires<3; ires++) {
               if (SelResidues[ires]->GetSeqNum() == (this_resno-1)) {
                  // previous residue
                  mmdb::Atom *atom = SelResidues[ires]->GetAtom(" C  ");
                  if (atom) {
                     n_atoms++;
                     prev_C = clipper::Coord_orth(atom->x, atom->y, atom->z);
                  }
               }
               if (SelResidues[ires]->GetSeqNum() == this_resno) {
                  // this residue
                  mmdb::Atom *atom;

                  atom = SelResidues[ires]->GetAtom(" N  ");
                  if (atom) {
                     n_atoms++;
                     this_N = clipper::Coord_orth(atom->x, atom->y, atom->z);
                  }
                  atom = SelResidues[ires]->GetAtom(" CA ");
                  if (atom) {
                     n_atoms++;
                     this_Ca = clipper::Coord_orth(atom->x, atom->y, atom->z);
                  }
                  atom = SelResidues[ires]->GetAtom(" C  ");
                  if (atom) {
                     n_atoms++;
                     this_C = clipper::Coord_orth(atom->x, atom->y, atom->z);
                  }
               }
               if (SelResidues[ires]->GetSeqNum() == this_resno +1) {
                  // next residue
                  mmdb::Atom *atom = SelResidues[ires]->GetAtom(" N  ");
                  if (atom) {
                     n_atoms++;
                     next_N = clipper::Coord_orth(atom->x, atom->y, atom->z);
                  }
               }
            }

            if (n_atoms == 5) {
               double phi = clipper::Util::rad2d(clipper::Coord_orth::torsion(prev_C, this_N, this_Ca, this_C));
               double psi = clipper::Util::rad2d(clipper::Coord_orth::torsion(this_N, this_Ca, this_C, next_N));
               pp.first = phi;
               pp.second = psi;
            } else {
               std::cout << "WARNING found " << n_atoms << " atoms (not 5) "
                         << "Can't get phi psi" << std::endl;
            }
         } else {
            std::cout << "WARNING: found " << nSelResdues << " residues (not 3) "
                      << "Can't get phi psi" << std::endl;
         }
      }
   }

   return pp;
}


// short int Have_modifications_p() const { return history_index > 0 ? 1 : 0;}

bool
molecule_class_info_t::Have_redoable_modifications_p() const {

   bool r = false;
//     std::cout << "DEBUG:: redoable? history_index: " << history_index
//               << " max_history_index: "
//               << max_history_index << " " << name_ << std::endl;

   if (history_index < max_history_index) {
      // When there are 3 backups made and we are viewing molecule 2,
      // we don't want to restore from history_filename_vec[3]:
      //
      if ( int(history_filename_vec.size()) > (history_index + 1)) {
         r = 1;
      }
   } else {
      r = false;
   }
   return r;
}

int
molecule_class_info_t::get_history_index() const {

   return history_index;
}

// Do not be mislead this is only a flag, *not* the number of chi
// angles for this residue
int
molecule_class_info_t::N_chis(int atom_index) {

   mmdb::Residue *res_p = atom_sel.atom_selection[atom_index]->residue;
   int r;

   std::string resname(res_p->GetResName());
   graphics_info_t g;
   // we want to ask g->Geom_p() if it has torsions for residues of
   // this type.

   if ( (resname == "GLY") || (resname == "ALA") ) {
      r = 0;
   } else {
      if (g.Geom_p()->have_dictionary_for_residue_type(resname, imol_no,
                                                       graphics_info_t::cif_dictionary_read_number)) {
         std::vector <coot::dict_torsion_restraint_t> v =
            g.Geom_p()->get_monomer_torsions_from_geometry(resname, 0);

         if (v.size() > 0)
            r = v.size();
         else
            r = 0;
      } else {
         r = 0;
      }

   }

   return r;

}

// How are clashes scored?  Hmmm... well, I think no clashes at all should have a score
// of 0.0 (c.f. auto_best_fit_rotamer()).  I think a badly crashing residue should have a
// score of around 1000.  A single 2.0A crash will have a score of 16.7 and a 1.0A crash
// 66.7.
//
std::pair<float, std::vector<mmdb::Atom *> >
molecule_class_info_t::get_clash_score(const coot::minimol::molecule &a_rotamer, bool score_hydrogen_atoms_flag,
                                       int water_interaction_mode) const {

   float score = 0;
   std::vector<mmdb::Atom *> clashing_waters;
   float dist_crit = 2.7;

   // First, where is the middle of the rotamer residue atoms and what
   // is the mean and maximum distance of coordinates from that point?

   // double std_dev_residue_pos;

   std::pair<double, clipper::Coord_orth> rotamer_info = get_minimol_pos(a_rotamer);
   double max_dev_residue_pos = rotamer_info.first;
   clipper::Coord_orth mean_residue_pos = rotamer_info.second;
   if (rotamer_info.first < 0.0) {
      // there were no atoms then
      std::cout << "ERROR: clash score: there are no atoms in the residue" << std::endl;
   } else {

      // So now we know the centre of the residue and the maximum distance of one of its
      // atoms from the centre.  Now let's run over the atoms of the atom_sel.
      //
      // When we find a distance between the middle of the residue and an atom_sel atom
      // that is less than (max_dev_residue_pos + distance), then we have found
      // potentially clashing atoms, so check for a clash of that atom_sel atom with all
      // the atoms of the residue.
      double d;
      double d_atom;
      float badness;
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
         mmdb::Atom *at = atom_sel.atom_selection[i];
         clipper::Coord_orth atom_sel_atom_pos = coot::co(atom_sel.atom_selection[i]);
         std::string res_name(atom_sel.atom_selection[i]->residue->GetResName());
         // either ignore waters or accumulate them for deletion if too close
         if (res_name == "HOH") {
            if (water_interaction_mode == 0) {
               for (unsigned int ifrag=0; ifrag<a_rotamer.fragments.size(); ifrag++) {
                  for (int ires=a_rotamer[ifrag].min_res_no(); ires<=a_rotamer[ifrag].max_residue_number(); ires++) {
                     for (auto it=a_rotamer[ifrag][ires].atoms.begin(); it != a_rotamer[ifrag][ires].atoms.end(); it++) {
                        if (score_hydrogen_atoms_flag || it->element != " H") {
                           double dd = (it->pos - atom_sel_atom_pos).lengthsq();
                           if (dd < 2.6 * 2.6) {
                              clashing_waters.push_back(at);
                           }
                        }
                     }
                  }
               }
            }
         } else {
            d = clipper::Coord_orth::length(atom_sel_atom_pos, mean_residue_pos);
            if (d < (max_dev_residue_pos + dist_crit)) {
               for (unsigned int ifrag=0; ifrag<a_rotamer.fragments.size(); ifrag++) {
                  for (int ires=a_rotamer[ifrag].min_res_no(); ires<=a_rotamer[ifrag].max_residue_number(); ires++) {
                     std::string residue_name = a_rotamer[ifrag][ires].name;
                     bool is_standard_aa = false;
                     if (coot::util::is_standard_residue_name(residue_name))
                        is_standard_aa = true;
                     for (unsigned int iat=0; iat<a_rotamer[ifrag][ires].n_atoms(); iat++) {
                        if (score_hydrogen_atoms_flag || a_rotamer[ifrag][ires][iat].element != " H") {
                           d_atom = clipper::Coord_orth::length(a_rotamer[ifrag][ires][iat].pos, atom_sel_atom_pos);
                           if (d_atom < dist_crit) {
                              int atom_sel_atom_resno = atom_sel.atom_selection[i]->GetSeqNum();
                              std::string atom_sel_atom_chain(atom_sel.atom_selection[i]->GetChainID());

                              // we should test residue specs here
                              if (! ((ires == atom_sel_atom_resno) &&
                                     (a_rotamer[ifrag].fragment_id == atom_sel_atom_chain)) ) {

                                 if ( (!is_standard_aa) ||
                                      (is_standard_aa && ! coot::is_main_chain_p(a_rotamer[ifrag][ires][iat].name))) {
                                    badness = 100.0 * (1.0/d_atom - 1.0/dist_crit);
                                    if (badness > 100.0)
                                       badness = 100.0;
                                    //               std::cout << "DEBUG:: adding clash badness " << badness
                                    //                         << " for atom "
                                    //                  << a_rotamer[ifrag][ires][iat].name << " with d_atom = "
                                    //                  << d_atom << " to sel atom "
                                    //                         << atom_sel.atom_selection[i] << std::endl;
                                    score += badness;
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
   std::pair<float, std::vector<mmdb::Atom *> > p(score, clashing_waters);
   return p;
}


// Return a negative value in the pair.first if there were no atoms in a_rotamer.
// Also, print an error message because (AFAICS) it should never happen.
//
std::pair<double, clipper::Coord_orth>
molecule_class_info_t::get_minimol_pos(const coot::minimol::molecule &a_rotamer) const {

   return a_rotamer.get_pos();
}

// We have name_.  We use that to find backup files in coot-backup
// dir.  If there is one and it is more recent that the file name_
// then restore for it and return 1.  Else return 0.
coot::backup_file_info
molecule_class_info_t::recent_backup_file_info() const {

   coot::backup_file_info info;

#if !defined(_MSC_VER)
   if (has_model()) {
      std::string t_name_glob = name_;
      // convert "/" to "_"
      // and in Windows the ":" to "_" as well
      int slen = t_name_glob.length();
      for (int i=0; i<slen; i++)
#ifdef WINDOWS_MINGW
         if (t_name_glob[i] == '/' || t_name_glob[i] == ':')
            t_name_glob[i] = '_';
#else
         if (t_name_glob[i] == '/')
            t_name_glob[i] = '_';
#endif // MINGW

      // Let's make a string that we can glob:
      // "coot-backup/thing.pdb*.pdb.gz"

      // c.f. make_backup():
      char *es = getenv("COOT_BACKUP_DIR");
      std::string backup_name_glob = "coot-backup/";
      // very first check if COOT_BACKUP_DIR is defined
      if (es) {
        // first we shall check if es, i.e. COOT_BACKUP_DIR actually exists
        struct stat buf;
        int err = stat(es, &buf);
        if (!err) {
          if (! S_ISDIR(buf.st_mode)) {
            es = NULL;
          }
        } else {
          es = NULL;
        }
      }
      if (es) {
         backup_name_glob = es;
         // on windows we somehow need to add an /
#ifdef WINDOWS_MINGW
         backup_name_glob += "/";
#endif // MINGW
      }
      backup_name_glob += t_name_glob;

      // First only the ones withwout gz
      backup_name_glob += "*.pdb";

      glob_t myglob;
      int flags = 0;
      glob(backup_name_glob.c_str(), flags, 0, &myglob);
      // And finally the ones with gz
      backup_name_glob += ".gz";
      flags = GLOB_APPEND;
      glob(backup_name_glob.c_str(), flags, 0, &myglob);
      size_t count;

      char **p;
      std::vector<std::string> v;
      for (p = myglob.gl_pathv, count = myglob.gl_pathc; count; p++, count--) {
         std::string f(*p);
         v.push_back(f);
      }
      globfree(&myglob);

      if (v.size() > 0) {

         struct stat buf;

         // get the time of name_ file
         int status = stat(name_.c_str(), &buf);
         if (status == 0) {

            time_t name_mtime;
            name_mtime = buf.st_mtime;

            // now stat for dated back files:
            time_t mtime;
            time_t mtime_youngest = 1000000000;
            short int set_something = 0;
            std::string backup_filename;

            for (unsigned int i=0; i<v.size(); i++) {
               status = stat(v[i].c_str(),&buf);
               if (status == 0) {
                  mtime = buf.st_mtime;
                  if (mtime > mtime_youngest) {
                     set_something = 1;
                     mtime_youngest = mtime;
                     backup_filename = v[i];
                  }
               }
            }
            if (set_something) {
               if (mtime_youngest > name_mtime) {
//                   std::cout << "Restoring from a recent backup "
//                             << backup_filename << std::endl;
                  info.name = name_;
                  info.backup_file_name = backup_filename;
                  info.status = 1; // There is a file.
               }
            }
         }
      }
   }
#endif
   return info;
}


short int
molecule_class_info_t::execute_restore_from_recent_backup(std::string backup_file_name,
                                                          std::string cwd) {

   // consider passing this:
   bool convert_flag = graphics_info_t::convert_to_v2_atom_names_flag;
   bool allow_duplseqnum = graphics_info_t::allow_duplseqnum;

   // std::cout << "Recovering from file: " << backup_file_name << std::endl;

   std::string save_name = name_;
   int save_imol = imol_no;
   // similarly, it messes with the save_state_command_strings_, we
   // don't want that either:
   std::vector<std::string> save_save_state = save_state_command_strings_;
   short int is_undo_or_redo = 1;
   short int reset_rotation_centre_flag = 0;
   handle_read_draw_molecule(imol_no, backup_file_name,
                             cwd,
                             graphics_info_t::Geom_p(),
                             reset_rotation_centre_flag,
                             is_undo_or_redo,
                             allow_duplseqnum,
                             convert_flag,
                             bond_width, Bonds_box_type(),
                             false);
   save_state_command_strings_ = save_save_state;
   imol_no = save_imol;
   name_ = save_name;
   have_unsaved_changes_flag = 1;

   return 0;
}


void
molecule_class_info_t::convert_rgb_to_hsv_in_place(const float *rgb, float *hsv) const {

   // convert to hsv
   float maxc = -1.0;
   float minc = 9.0;

   for (int i=0; i<3; i++) {
      if (maxc < rgb[i]) maxc = rgb[i];
      if (minc > rgb[i]) minc = rgb[i];
   }
   hsv[2] = maxc;

   if (minc == maxc) {
      hsv[0] = 0.0;
      hsv[1] = 0.0;
      hsv[2] = maxc;
   } else {

      float range = maxc - minc;
      hsv[1] = range/maxc;
      float rc = (maxc - rgb[0]) / range;
      float gc = (maxc - rgb[1]) / range;
      float bc = (maxc - rgb[2]) / range;
      if (rgb[0] == maxc){
         hsv[0] = bc-gc;
      } else {
         if (rgb[1]==maxc) {
            hsv[0] = 2.0+rc-bc;
         } else {
            hsv[0] = 4.0 + gc-rc;
         }
      }
      hsv[0] = hsv[0]/6.0- floorf(hsv[0]/6.0);
   }
}

void
molecule_class_info_t::convert_hsv_to_rgb_in_place(const float* hsv, float *rgb) const {


   if (hsv[1] == 0.0) {
      rgb[0] = hsv[2];
      rgb[1] = hsv[2];
      rgb[2] = hsv[2];
   } else {
      float fi = floorf(hsv[0]*6.0);
      float f  = (hsv[0]*6.0) - fi;
      float p = hsv[2]*(1.0 - hsv[1]);
      float q = hsv[2]*(1.0 - hsv[1]*f);
      float t = hsv[2]*(1.0 - hsv[1]*(1.0-f));

      int i = int(fi);
      switch (i) {

      case 0:
      case 6:
         rgb[0] = hsv[2];
         rgb[1] = t;
         rgb[2] = p;
         break;

      case 1:
         rgb[0] = q;
         rgb[1] = hsv[2];
         rgb[2] = p;
         break;

      case 2:
         rgb[0] = p;
         rgb[1] = hsv[2];
         rgb[2] = t;
         break;

      case 3:
         rgb[0] = p;
         rgb[1] = q;
         rgb[2] = hsv[2];
         break;

      case 4:
         rgb[0] = t;
         rgb[1] = p;
         rgb[2] = hsv[2];
         break;

      case 5:
         rgb[0] = hsv[2];
         rgb[1] = p;
         rgb[2] = q;
         break;
      }
   }

}


// for widget label:
std::string
molecule_class_info_t::cell_text_with_embeded_newline() const {

   std::string s;

   // NXMAP-FIXME.

   if (has_xmap()) {

      s = "   ";
      s += graphics_info_t::float_to_string(xmap.cell().descr().a());
      s += " ";
      s += graphics_info_t::float_to_string(xmap.cell().descr().b());
      s += " ";
      s += graphics_info_t::float_to_string(xmap.cell().descr().c());
      s += "\n   ";
      s += graphics_info_t::float_to_string(clipper::Util::rad2d(xmap.cell().descr().alpha()));
      s += " ";
      s += graphics_info_t::float_to_string(clipper::Util::rad2d(xmap.cell().descr().beta()));
      s += " ";
      s += graphics_info_t::float_to_string(clipper::Util::rad2d(xmap.cell().descr().gamma()));
   }
   return s;
}


// http://ngfnblast.gbf.de/docs/fasta.html
//
// For those programs that use amino acid query sequences (BLASTP and
// TBLASTN), the accepted amino acid codes are:
//
//     A  alanine                         P  proline
//     B  aspartate or asparagine         Q  glutamine
//     C  cystine                         R  arginine
//     D  aspartate                       S  serine
//     E  glutamate                       T  threonine
//     F  phenylalanine                   U  selenocysteine
//     G  glycine                         V  valine
//     H  histidine                       W  tryptophane
//     I  isoleucine                      Y  tyrosine
//     K  lysine                          Z  glutamate or glutamine
//     L  leucine                         X  any
//     M  methionine                      *  translation stop
//     N  asparagine                      -  gap of indeterminate length

// sequence
void
molecule_class_info_t::assign_fasta_sequence(const std::string &chain_id, const std::string &seq_in) {

   std::cout << "INFO:: assign_fasta_sequence " << chain_id << "\n" << seq_in << std::endl;

   // format "> name\n <sequence>", we ignore everything that is not a
   // letter after the newline.

   // sequence is member data.  Let's fill it.

   std::string seq;

   int nchars = seq_in.length();
   short int found_greater = 0;
   short int found_newline = 0;
   std::string t;

   for (int i=0; i<nchars; i++) {

      // std::cout << "checking character: " << seq_in[i] << std::endl;

      if (found_newline && found_greater) {
         t = std::toupper(seq_in[i]);
         if (is_fasta_aa(t)) {
            //             std::cout << "adding character: " << seq_in[i] << std::endl;
            seq += t;
         }
      }
      if (seq_in[i] == '>') {
         //          std::cout << "DEBUG:: " << seq_in[i] << " is > (greater than)\n";
         found_greater = 1;
      }
      if (seq_in[i] == '\n') {
         if (found_greater) {
            //             std::cout << "DEBUG:: " << seq_in[i] << " is carriage return\n";
            found_newline = 1;
         }
      }
   }

   if (seq.length() > 0) {
      std::cout << "debug:: assign_fasta_sequence(): storing sequence: " << seq << " for chain id: " << chain_id
                << " in molecule number " << imol_no << std::endl;
      std::cout << "debug:: pushing back onto input_sequence" << std::endl;
      input_sequence.push_back(std::pair<std::string, std::string> (chain_id,seq));
      std::cout << "debug:: input_sequence size " << input_sequence.size() << std::endl;
   } else {
      std::cout << "WARNING:: assign_fasta_sequence(): no sequence found or improper fasta sequence format\n";
   }
}

// add the sequence the file (read depending on file name) to input_sequence vector (chain-id is blank
// as it could apply to any chain)
void
molecule_class_info_t::associate_sequence_from_file(const std::string &seq_file_name) {

   auto file_to_string = [] (const std::string &file_name) {
                            std::string s;
                            std::string line;
                            std::ifstream f(file_name.c_str());
                            if (!f) {
                               std::cout << "WARNING:: Failed to open " << file_name << std::endl;
                            } else {
                               while (std::getline(f, line)) {
                                  s += line;
                                  s += "\n";
                               }
                            }
                            return s;
                         };

   std::string extension = coot::util::file_name_extension(seq_file_name);
   std::string chain_id;
   if (coot::file_exists(seq_file_name)) {
      std::string seq = file_to_string(seq_file_name);
      if (extension == ".pir") {
         assign_pir_sequence(chain_id, seq);
      } else {
         assign_fasta_sequence(chain_id, seq);
      }
   } else {
      std::cout << "WARNING:: file does not exist: " << seq_file_name << std::endl;
   }
}

// this is not assigning the sequence! This is adding a PIR file for a particular chain id!
void
molecule_class_info_t::assign_pir_sequence(const std::string &chain_id, const std::string &seq_in) {

   // format "> ID;database-id\ntext-descr\n <sequence>", we ignore everything that is not a
   // letter after the newline.

   // sequence is member data.  Let's fill it.
   // std::cout << "in assign_pir_sequence\n";

   std::string seq;

   int nchars = seq_in.length();
   short int found_greater = 0;
   short int found_newline = 0;
   short int found_textdescr = 0;
   std::string t;

   for (int i=0; i<nchars; i++) {

      // std::cout << "checking character: " << seq_in[i] << std::endl;

      if (found_newline && found_greater && found_textdescr) {
         t = std::toupper(seq_in[i]);
         if (is_pir_aa(t)) {
            //             std::cout << "adding character: " << seq_in[i] << std::endl;
            seq += t;
            if (t == "*") // end of sequence
               break; // the for loop
         }
      }
      if (seq_in[i] == '>') {
         //          std::cout << "DEBUG:: " << seq_in[i] << " is > (greater than)\n";
         found_greater = 1;
      }
      if (seq_in[i] == '\n') {
         if (found_newline) {
            //             std::cout << "DEBUG:: " << seq_in[i] << " is carriage return\n";
            found_textdescr = 1;
         }
         if (found_greater) {
            //             std::cout << "DEBUG:: " << seq_in[i] << " is carriage return\n";
            found_newline = 1;
         }
      }
   }

   if (seq.length() > 0) {
      std::cout << "debug:: assign_pir_sequence():: storing sequence: "
                << seq << " for chain id: " << chain_id
                << " in molecule number " << imol_no << std::endl;
      // replace the sequence assigned to chain_id is possible
      bool found_chain = false;
      for (unsigned int ich=0; ich<input_sequence.size(); ich++) {
         if (input_sequence[ich].first == chain_id) {
            found_chain = true;
            input_sequence[ich].second = seq;
            break;
         }
      }

      if (! found_chain) {
         std::cout << "debug:: assign_pir_sequence() input_sequence pushing back for chain-id: \"" << chain_id << "\" seq: "
                   << seq << " imol " << imol_no << std::endl;
         input_sequence.push_back(std::pair<std::string, std::string> (chain_id, seq));
      }
   } else {
      std::cout << "WARNING:: assign_pir_sequence() no sequence found or improper pir sequence format\n";
   }

}

std::vector<std::string>
molecule_class_info_t::get_chain_ids() const {

   std::vector<std::string> v;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::string chain_id(chain_p->GetChainID());
         v.push_back(chain_id);
      }
   }
   return v;
}


// we let clipper assign the sequence for us from any sequence file
void
molecule_class_info_t::assign_sequence_from_file(const std::string &filename) {

   if (! atom_sel.mol) return;

   if (! coot::file_exists(filename)) {
      std::cout << "Sequence file not found: " << filename << std::endl;
   } else {

      clipper::SEQfile seq_file;
      clipper::MMoleculeSequence molecule_sequence;
      seq_file.read_file(filename);
      seq_file.import_molecule_sequence(molecule_sequence);

      std::vector<std::string> chain_ids = get_chain_ids();
      input_sequence.clear();
      for (unsigned int i=0; i<chain_ids.size(); i++) {
         const std::string &chain_id = chain_ids[i];

         int selHnd = atom_sel.mol->NewSelection(); // d
         mmdb::PResidue *SelResidues = NULL;
         int nSelResidues = 0;
         float wgap   = -3.0;
         float wspace = -0.4;

         atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
                              chain_id.c_str(),
                              mmdb::ANY_RES, "*",
                              mmdb::ANY_RES, "*",
                              "*",  // residue name
                              "*",  // Residue must contain this atom name?
                              "*",  // Residue must contain this Element?
                              "*",  // altLocs
                              mmdb::SKEY_NEW // selection key
                              );
         atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
         if (nSelResidues > 0) {
            float current_best_alignment_score = -1.0;
            std::string current_best_sequence;
            for (int j=0; j<molecule_sequence.size(); j++) {
               std::string target = molecule_sequence[j].sequence();
               coot::chain_mutation_info_container_t alignment =
                  align_on_chain(chain_id, SelResidues, nSelResidues,
                                 target, wgap, wspace, false, false);

               std::cout << "chain_id " << chain_id
                         << " alignment_score " << alignment.alignment_score.first
                         << " " << alignment.alignment_score.second
                         << " n-alignment-mutations " << alignment.mutations.size()
                         << " with " << nSelResidues << " residues in chain" << std::endl;

               if (alignment.alignment_score.first) {
                  if (alignment.alignment_score.second > 2.0 * 0.7 * nSelResidues) {
                     if (alignment.alignment_score.second > current_best_alignment_score) {
                        current_best_alignment_score = alignment.alignment_score.second;
                        current_best_sequence = target;
                     }
                  }
               }
            }
            if (! current_best_sequence.empty()) {
               std::pair<std::string, std::string> new_seq_info(chain_id, current_best_sequence);
               input_sequence.push_back(new_seq_info);
            }
         }
         atom_sel.mol->DeleteSelection(selHnd);
      }
   }

   if (true) {
      std::cout << "Now we have these sequences: " << std::endl;
      for (unsigned int i=0; i<input_sequence.size(); i++) {
         const std::string chain_id = input_sequence[i].first;
         const std::string seq = input_sequence[i].second;
         std::cout << "chain " << chain_id << "  " << seq << std::endl;
      }
   }

}

// to assign a sequence from a simple string
// Apply to NCS-related chains too if present.
//
void
molecule_class_info_t::assign_sequence_to_NCS_related_chains_from_string(const std::string &chain_id,
                                                                         const std::string &seq_in) {

   // std::cout << "in assign_sequence_from_string\n";

   std::string seq = seq_in;
   if (seq.length() > 0) {
      input_sequence.push_back(std::pair<std::string, std::string> (chain_id, seq));

      std::vector<std::string> ncs_related_chain_ids;
      for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
         if (ncs_ghosts[ighost].chain_id == chain_id)
            // add the target/master chain-id if it is not alreay in ncs_related_chain_ids
            if (std::find(ncs_related_chain_ids.begin(),
                          ncs_related_chain_ids.end(), ncs_ghosts[ighost].target_chain_id) ==
                ncs_related_chain_ids.end())
               ncs_related_chain_ids.push_back(ncs_ghosts[ighost].target_chain_id);
         if (ncs_ghosts[ighost].target_chain_id == chain_id)
            // add the ghost chain-id if it is not alreay in ncs_related_chain_ids
            if (std::find(ncs_related_chain_ids.begin(),
                          ncs_related_chain_ids.end(), ncs_ghosts[ighost].chain_id) ==
                ncs_related_chain_ids.end())
               ncs_related_chain_ids.push_back(ncs_ghosts[ighost].chain_id);
      }

      if (ncs_related_chain_ids.size())
         for (unsigned int in=0; in<ncs_related_chain_ids.size(); in++)
            assign_sequence_from_string_simple(ncs_related_chain_ids[in], seq);

   } else {
      std::cout << "WARNING:: assign_sequence_to_NCS_related_chains_from_string() no sequence found or improper string\n";
   }
}

void
molecule_class_info_t::assign_sequence_from_string_simple(const std::string &chain_id, const std::string &seq_in) {

   std::string seq = seq_in;
   if (seq.length() > 0) {
      std::cout << "debug:: assign_sequence_from_string_simple() storing sequence: "
                << seq << " for chain id: " << chain_id
                << " in molecule number " << imol_no << std::endl;
      input_sequence.push_back(std::pair<std::string, std::string> (chain_id, seq));
   }
}


// Delete all the associated sequences from the molecule
void
molecule_class_info_t::delete_all_sequences_from_molecule() {

  std::vector<std::pair<std::string, std::string> > seq =
    graphics_info_t::molecules[imol_no].sequence_info();
  input_sequence.clear();
  seq = graphics_info_t::molecules[imol_no].sequence_info();
}

// Delete the to chain_id associated sequence from the molecule
void
molecule_class_info_t::delete_sequence_by_chain_id(const std::string &chain_id_in) {

  std::vector<std::pair<std::string, std::string> > seq =
    graphics_info_t::molecules[imol_no].sequence_info();
  std::vector<std::pair<std::string, std::string> >::iterator iter;
  for (iter = input_sequence.begin(); iter != input_sequence.end(); iter++) {
    std::string chain_id = (*(iter)).first;
    if (chain_id == chain_id_in) {
      // delete the seq for this chain
      input_sequence.erase(iter);
      break;
    }
  }
  seq = graphics_info_t::molecules[imol_no].sequence_info();
}


// Return a flag that tells us if we did indeed find a proper next
// residue.  Return the 3-letter-code in second.
//
std::pair<bool, std::string>
molecule_class_info_t::residue_type_next_residue_by_alignment(const coot::residue_spec_t &clicked_residue,
                                                              mmdb::Chain *clicked_residue_chain_p,
                                                              short int is_n_term_addition,
                                                              mmdb::realtype alignment_wgap,
                                                              mmdb::realtype alignment_wspace) const {

   std::pair<bool, std::string> p(0, "");

   if (input_sequence.size() > 0) {
      std::string chain_id = clicked_residue.chain_id;
      for (unsigned int ich=0; ich<input_sequence.size(); ich++) {

         if (input_sequence[ich].first == chain_id) {

            if (input_sequence[ich].second.length() > 0) {
               std::vector<mmdb::PResidue> frag_residues =
                  coot::util::get_residues_in_fragment(clicked_residue_chain_p, clicked_residue);
               // copy from vector to array
               mmdb::PResidue *SelResidues = new mmdb::PResidue[frag_residues.size()];
               for (unsigned int ires=0; ires<frag_residues.size(); ires++)
                  SelResidues[ires] = frag_residues[ires];

               coot::chain_mutation_info_container_t a =
                  align_on_chain(chain_id, SelResidues, frag_residues.size(),
                                 input_sequence[ich].second,
                                 alignment_wgap, alignment_wspace);

               if ((a.insertions.size() +
                    a.mutations.size() +
                    a.deletions.size()) > (input_sequence[ich].second.length()/5)) {
                  std::cout << "WARNING:: Too many mutations, "
                            << "can't make sense of aligment "
                            << a.insertions.size() << " "
                            << a.mutations.size() << " "
                            << a.deletions.size() << " "
                            << input_sequence[ich].second.length()
                            << std::endl;
               } else {
                  // proceed
                  std::cout << a.alignedS << std::endl;
                  std::cout << a.alignedT << std::endl;

                  // where is clicked_residue in SelResidues?
                  bool found = 0;
                  int frag_seqnum;
                  for (unsigned int ires=0; ires<input_sequence[ich].second.length(); ires++) {
                     if (SelResidues[ires]->GetSeqNum() == clicked_residue.res_no) {
                        if (clicked_residue.chain_id == SelResidues[ires]->GetChainID()) {
                           // found clicked_residue
                           found = 1;
                           frag_seqnum = ires;
//                            std::cout << "DEBUG:: found frag_seqnum: " << frag_seqnum
//                                      << " " << SelResidues[ires]->GetSeqNum() << " "
//                                      << std::endl;
                           break;
                        }
                     }
                  }
                  if (found) {
                     int added_res_resno;
                     if (is_n_term_addition) {
                        added_res_resno = frag_seqnum - 1;
                     } else {
                        added_res_resno = frag_seqnum + 1;
                     }

                     if (int(a.alignedT.length()) > added_res_resno) {
                        if (added_res_resno >= 0) {
                           char code = a.alignedT[added_res_resno];
                           std::cout << " code: " << code << std::endl;
                           std::string res =
                              coot::util::single_letter_to_3_letter_code(code);
                           p = std::pair<bool, std::string>(1, res);
                           for (int off=5; off>=0; off--) {
                               char c = a.alignedT[added_res_resno-off];
                              std::cout << c;
                           }
                           std::cout << std::endl;
                        }
                     }
                  }
               }
               delete [] SelResidues;
            }
            break;
         }
      }
   }
   return p;
}



bool
molecule_class_info_t::is_fasta_aa(const std::string &a) const {

   short int r = 0;

   if (a == "A" || a == "G" ) {
      r = 1;
   } else {
      if (a == "B"
          || a == "C" || a == "D" || a == "E" || a == "F" || a == "H" || a == "I"
          || a == "K" || a == "L" || a == "M" || a == "N" || a == "P" || a == "Q"
          || a == "R" || a == "S" || a == "T" || a == "U" || a == "V" || a == "W"
          || a == "Y" || a == "Z" || a == "X" || a == "*" || a == "-") {
         r = 1;
      }
   }
   return r;
}

bool
molecule_class_info_t::is_pir_aa(const std::string &a) const {

   bool r = false;

   // 20130820 Allow U for RNA sequence assignment.

   if (a == "A" || a == "G" ) {
      r = 1;
   } else {
      if (   a == "C" || a == "D" || a == "E" || a == "F" || a == "H" || a == "I"
          || a == "K" || a == "L" || a == "M" || a == "N" || a == "P" || a == "Q"
          || a == "R" || a == "S" || a == "T" ||             a == "V" || a == "W"
          || a == "Y" || a == "Z" || a == "X"
          || a == "U"
             ) {
         r = 1;
      }
   }
   return r;
}

// render option (other functions)
coot::ray_trace_molecule_info
molecule_class_info_t::fill_raster_model_info(bool against_a_dark_background) {

   // 20150803-PE FIXME - pass Geom_p().
   graphics_info_t g;
   coot::protein_geometry *geom_p = g.Geom_p();

   coot::ray_trace_molecule_info rtmi;
   if (has_model()) {
      if (draw_it) {
         int restore_bonds = 0;
         if (g.raster3d_water_sphere_flag && bonds_box_type == coot::NORMAL_BONDS) {
            // remove waters
            bonds_no_waters_representation();
            restore_bonds = 1;
         }
         rtmi.bond_lines.resize(bonds_box.num_colours);
         for (int i=0; i<bonds_box.num_colours; i++) {
            set_bond_colour_by_mol_no(i, against_a_dark_background); //sets bond_colour_internal
            double thickness = g.raster3d_bond_thickness;
            if (bonds_box.bonds_[i].thin_lines_flag) thickness *= 0.5;
            for (int j=0; j<bonds_box.bonds_[i].num_lines; j++) {
               coot::ray_trace_molecule_info::bond_t b(bonds_box.bonds_[i].pair_list[j].positions.getStart(),
                                                       bonds_box.bonds_[i].pair_list[j].positions.getFinish(),
                                                       thickness);
               rtmi.bond_lines[i].bonds.push_back(b);
            }
            coot::colour_t c;
            c.col.resize(3);
            c.col[0] = bond_colour_internal[0];
            c.col[1] = bond_colour_internal[1];
            c.col[2] = bond_colour_internal[2];
            rtmi.bond_lines[i].colour = c;
         }
         // restore bond_box_type
         if (restore_bonds) {
            bonds_box_type = coot::NORMAL_BONDS;
            std::set<int> dummy;
            makebonds(geom_p, dummy);
         }

         std::cout << " There are " << bonds_box.n_atom_centres_
                   << " atom centres in this bonds_box\n";

         for (int i=0; i<bonds_box.n_atom_centres_; i++) {
            coot::colour_t c;
            c.col.resize(3);
            //sets bond_colour_internal
            set_bond_colour_by_mol_no(bonds_box.atom_centres_colour_[i], against_a_dark_background);
            c.col[0] = bond_colour_internal[0];
            c.col[1] = bond_colour_internal[1];
            c.col[2] = bond_colour_internal[2];
            // std::cout << " bonds_box for atoms " << i << " col " << c << std::endl;
            // here is the place to add tiny rastered hydrogen balls.
            // rtmi.atom.push_back(std::pair<coot::Cartesian, coot::colour_t>
            // (bonds_box.atom_centres_[i].second, c));
            double r = g.raster3d_atom_radius;
            // std::cout << "comparing colours " << bonds_box.atom_centres_colour_[i] << " vs "
            //           << HYDROGEN_GREY_BOND << std::endl;
            if (bonds_box.atom_centres_colour_[i] == HYDROGEN_GREY_BOND)
               r *= 0.5;

            coot::ray_trace_molecule_info::ball_t b(bonds_box.atom_centres_[i].position, c, r);
            rtmi.balls.push_back(b);
         }
         rtmi.molecule_name = name_;
         rtmi.molecule_number = imol_no;
      }
   }
   std::cout << "DEBUG:: Done fill_raster_model_info for "
             << imol_no << std::endl;
   return rtmi;
}


// Pass lev as +1, -1
//
coot::ray_trace_molecule_info
molecule_class_info_t::fill_raster_map_info(short int lev) const {

   coot::ray_trace_molecule_info rtmi;
   if (has_xmap()) {   // NXMAP-FIXME

      rtmi.bond_colour.resize(3);
      if (draw_it_for_map) {
         // if (draw_it_for_map_standard_lines) // let's see something
                                                // at least, if we are displaying a volume surface....
	if (1) {
	   if (lev == 1) {
	      if (! draw_vector_sets.empty()) {

	         rtmi.density_colour.col.resize(3);
                 rtmi.density_colour.col[0] = map_colour.red;
                 rtmi.density_colour.col[1] = map_colour.green;
                 rtmi.density_colour.col[2] = map_colour.blue;

                 for (std::size_t i=0; i<draw_vector_sets.size(); i++) {
                    for (unsigned int j=0; j<draw_vector_sets[i].point_indices.size(); j++) {
                       const clipper::Coord_orth &pt_1(draw_vector_sets[i].points[draw_vector_sets[i].point_indices[j].pointID[0]]);
                       const clipper::Coord_orth &pt_2(draw_vector_sets[i].points[draw_vector_sets[i].point_indices[j].pointID[1]]);
                       const clipper::Coord_orth &pt_3(draw_vector_sets[i].points[draw_vector_sets[i].point_indices[j].pointID[2]]);
                       // I have to do a bit of jiggery pokery here - I don't know why.
                       std::pair<coot::Cartesian, coot::Cartesian> p2;
                       p2.first  = coot::Cartesian(pt_1);
                       p2.second = coot::Cartesian(pt_2);
                       rtmi.density_lines.push_back(p2);
                       p2.first  = coot::Cartesian(pt_1);
                       p2.second = coot::Cartesian(pt_3);
                       rtmi.density_lines.push_back(p2);
                       p2.first  = coot::Cartesian(pt_2);
                       p2.second = coot::Cartesian(pt_3);
                       rtmi.density_lines.push_back(p2);
                    }
                 }
              }
           } else {
              if (! draw_diff_map_vector_sets.empty()) {

                 rtmi.density_colour.col.resize(3);
                 rtmi.density_colour.col[0] = map_colour_negative_level.red;
                 rtmi.density_colour.col[1] = map_colour_negative_level.green;
                 rtmi.density_colour.col[2] = map_colour_negative_level.blue;

                 for (std::size_t i=0; i<draw_diff_map_vector_sets.size(); i++) {
                    for (unsigned int j=0; j<draw_vector_sets[i].point_indices.size(); j++) {
                       const clipper::Coord_orth &pt_1(draw_diff_map_vector_sets[i].points[draw_diff_map_vector_sets[i].point_indices[j].pointID[0]]);
                       const clipper::Coord_orth &pt_2(draw_diff_map_vector_sets[i].points[draw_diff_map_vector_sets[i].point_indices[j].pointID[1]]);
                       const clipper::Coord_orth &pt_3(draw_diff_map_vector_sets[i].points[draw_diff_map_vector_sets[i].point_indices[j].pointID[2]]);
                       std::pair<coot::Cartesian, coot::Cartesian> p2;
                       p2.first  = coot::Cartesian(pt_1);
                       p2.second = coot::Cartesian(pt_2);
                       rtmi.density_lines.push_back(p2);
                       p2.first  = coot::Cartesian(pt_1);
                       p2.second = coot::Cartesian(pt_3);
                       rtmi.density_lines.push_back(p2);
                       p2.first  = coot::Cartesian(pt_2);
                       p2.second = coot::Cartesian(pt_3);
                       rtmi.density_lines.push_back(p2);
                    }
                 }
              }
           }
        }
      }

      if (fc_skeleton_draw_on == 1) {

         rtmi.bones_colour.col.resize(3);
         for (int i=0; i<3; i++)
            rtmi.bones_colour.col[i] = graphics_info_t::skeleton_colour[i];
         for (int l=0; l<fc_skel_box.num_colours; l++) {
            for (int j=0; j<fc_skel_box.bonds_[l].num_lines; j++) {
               std::pair<coot::Cartesian, coot::Cartesian>
                  p(fc_skel_box.bonds_[l].pair_list[j].positions.getStart(),
                    fc_skel_box.bonds_[l].pair_list[j].positions.getFinish());
               rtmi.bone_lines.push_back(p);
            }
         }
      }
      rtmi.molecule_name = name_;
      rtmi.molecule_number = imol_no;

   }
   return rtmi;
}

// For additional restraints.
coot::ray_trace_molecule_info
molecule_class_info_t::fill_raster_additional_info() const {

   coot::ray_trace_molecule_info rti;

   if (draw_it) {
      if (draw_it_for_extra_restraints) {

         double thickness = 0.02; // 0.012 for silkworm
         if (extra_restraints_representation_for_bonds_go_to_CA)
            thickness = 0.08;

         for (unsigned int ib=0; ib<extra_restraints_representation.bonds.size(); ib++) {
            const coot::extra_restraints_representation_t::extra_bond_restraints_respresentation_t &res =
               extra_restraints_representation.bonds[ib];

            // red if actual distance is greater than target
            //
            double d_sqd = (res.second - res.first).clipper::Coord_orth::lengthsq();

            if (res.esd > 0) {
               double b = (res.target_dist*res.target_dist - d_sqd)/res.esd * 0.005;
               if (b >  0.4999) b =  0.4999;
               if (b < -0.4999) b = -0.4999;
               double b_green = b;
               if (b > 0) b_green *= 0.2;
               coot::colour_t c(0.5-b, 0.5+b_green*0.9, 0.5+b);
               coot::Cartesian p1(res.first);
               coot::Cartesian p2(res.second);
               rti.add_extra_representation_line(p1, p2, c, thickness);
            }
         }

         coot::colour_t c(0.7, 0.7, 0.1);
         thickness = 0.04;
         for (unsigned int ipp=0; ipp<extra_restraints_representation.parallel_planes.size(); ipp++) {
            const coot::extra_restraints_representation_t::extra_parallel_planes_restraints_representation_t &ppr = extra_restraints_representation.parallel_planes[ipp];
            coot::Cartesian p1(ppr.ring_centre);
            coot::Cartesian p2(ppr.plane_projection_point);
            rti.add_extra_representation_line(p1, p2, c, thickness);
            // add a small ball at the projection point end:
            coot::ray_trace_molecule_info::ball_t b(p2 ,c, 0.08);
            rti.balls.push_back(b);


            // now the rings:
            //
            clipper::Coord_orth arb(0.2, 0.8, 0.1);
            clipper::Coord_orth cr(clipper::Coord_orth::cross(ppr.normal, arb).unit());
            clipper::Coord_orth first_pt = ppr.ring_centre + ppr.ring_radius * cr;
            clipper::Coord_orth first_pt_pp = ppr.plane_projection_point + ppr.pp_radius * cr;
            // std::cout << i << " r.plane_projection_point: " << r.plane_projection_point.format() << std::endl;

            unsigned int n_steps = 128;
            double step_frac = 1/double(n_steps);
            clipper::Coord_orth pt_1;
            clipper::Coord_orth pt_2;
            for (unsigned int istep=0; istep<n_steps; istep++) {
               double angle_1 = step_frac * 2.0 * M_PI * istep;
               double angle_2 = step_frac * 2.0 * M_PI * (istep + 1);
               pt_1 = coot::util::rotate_around_vector(ppr.normal, first_pt, ppr.ring_centre, angle_1);
               pt_2 = coot::util::rotate_around_vector(ppr.normal, first_pt, ppr.ring_centre, angle_2);
               coot::Cartesian p1(pt_1);
               coot::Cartesian p2(pt_2);
               rti.add_extra_representation_line(p1, p2, c, thickness);
            }

            n_steps = 64;
            step_frac = 1/double(n_steps);
            for (unsigned int istep=0; istep<n_steps; istep++) {
               double angle_1 = step_frac * 2.0 * M_PI * istep;
               double angle_2 = step_frac * 2.0 * M_PI * (istep + 1);
               pt_1 = coot::util::rotate_around_vector(ppr.normal, first_pt_pp, ppr.plane_projection_point, angle_1);
               pt_2 = coot::util::rotate_around_vector(ppr.normal, first_pt_pp, ppr.plane_projection_point, angle_2);
               coot::Cartesian p1(pt_1);
               coot::Cartesian p2(pt_2);
               rti.add_extra_representation_line(p1, p2, c, thickness);
            }

         }
      }
   }
   return rti;
}



int
molecule_class_info_t::trim_by_map(const clipper::Xmap<float> &xmap_in,
                                   float map_level, short int delete_or_zero_occ_flag) {


   // "Trim off the bits we don't need":
   // 20050610: I laugh at that.

   short int waters_only_flag = 0;
   int i = coot::util::trim_molecule_by_map(atom_sel.mol, xmap_in, map_level,
                                            delete_or_zero_occ_flag,
                                            waters_only_flag);

   std::cout << "INFO:: " << i << " atoms were trimmed\n";
   if (i > 0) {
      make_backup();
      update_molecule_after_additions(); // sets have_unsaved_changes_flag
   }
   return i;
}



std::vector<coot::atom_spec_t>
molecule_class_info_t::check_waters_by_difference_map(const clipper::Xmap<float> &xmap,
                                                      float sigma_level) const {

   std::vector<coot::atom_spec_t> v;

   std::vector<std::pair<coot::util::density_stats_info_t, coot::atom_spec_t> > dsi;

   std::pair<coot::util::density_stats_info_t, coot::atom_spec_t> pair;
   coot::atom_spec_t at_spec;

   for (int i=0; i<atom_sel.n_selected_atoms; i++) {
      std::string resname = atom_sel.atom_selection[i]->residue->name;
      if (! atom_sel.atom_selection[i]->isTer()) {
         if (resname == "WAT" || resname == "HOH") {
            std::string ele(atom_sel.atom_selection[i]->element);
            if (ele == " O") {
               clipper::Coord_orth p(atom_sel.atom_selection[i]->x,
                                     atom_sel.atom_selection[i]->y,
                                     atom_sel.atom_selection[i]->z);
               coot::atom_spec_t at_spec(atom_sel.atom_selection[i]->GetChainID(),
                                         atom_sel.atom_selection[i]->GetSeqNum(),
                                         atom_sel.atom_selection[i]->GetInsCode(),
                                         atom_sel.atom_selection[i]->GetAtomName(),
                                         atom_sel.atom_selection[i]->altLoc);
               pair = std::pair<coot::util::density_stats_info_t, coot::atom_spec_t>(coot::util::density_around_point(p, xmap, 1.5), at_spec);
               dsi.push_back(pair);
            }
         }
      }
   }

   float sum_v = 0.0;
   float sum_v_sq = 0.0;
   if (dsi.size()) {
      for (unsigned int id=0; id<dsi.size(); id++) {
         sum_v    += dsi[id].first.sum_sq;
         sum_v_sq += dsi[id].first.sum_sq * dsi[id].first.sum_sq;
      }

      float v_mean = sum_v/float(dsi.size());
      float v_variance = sum_v_sq/float(dsi.size()) - v_mean*v_mean;
      for (unsigned int id=0; id<dsi.size(); id++) {
         std::cout << "debug:: for atom " << dsi[id].second << " comparing "
                   << dsi[id].first.sum_sq << " - "
                   << v_mean << ")/sqrt(" << v_variance << ") = "
                   << (dsi[id].first.sum_sq - v_mean)/sqrt(v_variance)
                   << " > " << sigma_level
                   << std::endl;
         if ( (dsi[id].first.sum_sq - v_mean)/sqrt(v_variance) > sigma_level) {
            std::cout << "debug::   pushing back " << dsi[id].second << std::endl;
            v.push_back(dsi[id].second);
         }
      }
   }
   return v;
}


// Return a list of bad chiral volumes for this molecule:
//
// Return also a flag for the status of this test, were there any
// residues for we we didn't find restraints?  The flag is the number
// of residue names in first part of the returned pair.
//
//
std::pair<std::vector<std::string> , std::vector <coot::atom_spec_t> >
molecule_class_info_t::bad_chiral_volumes() const {

   return inverted_chiral_volumes();
}

std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t> >
molecule_class_info_t::inverted_chiral_volumes() const {

   std::vector <coot::atom_spec_t> v;
   std::vector<std::string> unknown_types_vec;
   std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t> > pair(unknown_types_vec, v);

   if (atom_sel.n_selected_atoms > 0) {
      // grr Geom_p() is not static
      graphics_info_t g;
      pair = coot::inverted_chiral_volumes(imol_no, atom_sel.mol, g.Geom_p(),
                                           graphics_info_t::cif_dictionary_read_number);
   }

   return pair;
}

std::pair<std::vector<std::string>, std::vector<std::pair<coot::atom_spec_t, double> > >
molecule_class_info_t::distorted_chiral_volumes(double chiral_volume_limit_for_outlier) const {

   std::pair<std::vector<std::string> , std::vector<std::pair<coot::atom_spec_t, double> > > p =
      coot::distorted_chiral_volumes(imol_no, atom_sel.mol, graphics_info_t::Geom_p(),
                                     graphics_info_t::cif_dictionary_read_number,
                                     chiral_volume_limit_for_outlier);

   return p;
}


float
molecule_class_info_t::score_residue_range_fit_to_map(int resno1, int resno2,
                                                      std::string altloc,
                                                      std::string chain_id,
                                                      int imol_for_map) {

   float f = 0;

   // Make an atom selection:
   // Convert atoms to coord_orths,
   // Use a coot util function to return the score.
   //

   int selHnd = atom_sel.mol->NewSelection();
   atom_sel.mol->SelectAtoms(selHnd, 0,  (char *) chain_id.c_str(),
                             resno1, "*",
                             resno2, "*",
                             "*", // residue name
                             "*", // atom name
                             "*", // element name
                             (char *) altloc.c_str());
   mmdb::PPAtom local_SelAtom = NULL;
   int nSelAtoms;
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "WARNING:: No atoms selected in "
                << "score_residue_range_fit_to_map\n";
   } else {
      f = coot::util::map_score(local_SelAtom, nSelAtoms,
                                graphics_info_t::molecules[imol_for_map].xmap,
                                0 // not score by atom type
                                );
      std::cout << "score for residue range " << resno1 << " " << resno2
                << " chain " <<  chain_id
                << ": " << f << std::endl;
   }
   atom_sel.mol->DeleteSelection(selHnd);
   return f;
}


void
molecule_class_info_t::fit_residue_range_to_map_by_simplex(int resno1, int resno2,
                                                           std::string altloc,
                                                           std::string chain_id,
                                                           int imol_for_map) {

   int selHnd = atom_sel.mol->NewSelection();
   atom_sel.mol->SelectAtoms(selHnd, 0,  chain_id.c_str(),
                             resno1, "*",
                             resno2, "*",
                             "*", // residue name
                             "*", // atom name
                             "*", // ?
                             altloc.c_str());
   mmdb::PPAtom local_SelAtom = NULL;
   int nSelAtoms;
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "WARNING:: No atoms selected in "
                << "score_residue_range_fit_to_map\n";
   } else {
      make_backup();
      coot::util::fit_to_map_by_simplex_rigid(local_SelAtom, nSelAtoms,
                                              graphics_info_t::molecules[imol_for_map].xmap);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }

   atom_sel.mol->DeleteSelection(selHnd);
}


// Return a pair, the first if the split was done correctly or not and
// the second is the string of the new alt conf.
//
// If this is not a shelxl molecule, don't add change the remaining
// atoms to negative occupancy, set them to zero.
//
std::pair<bool,std::string>
molecule_class_info_t::split_residue(int atom_index, int alt_conf_split_type) {

   std::pair<bool,std::string> pr(0, "");
   if (atom_index < atom_sel.n_selected_atoms) {
      int do_intermediate_atoms = 0;
      mmdb::Residue *res = atom_sel.atom_selection[atom_index]->residue;
      std::vector<std::string> residue_alt_confs =
         get_residue_alt_confs(res);
      std::string altconf(atom_sel.atom_selection[atom_index]->altLoc);
      int atom_index_udd = atom_sel.UDDAtomIndexHandle;

      int udd_afix_handle = -1; // don't use value
      if (is_from_shelx_ins_flag)
         udd_afix_handle = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");

      std::pair<mmdb::Residue *, atom_selection_container_t> p =
         coot::deep_copy_this_residue_and_make_asc(atom_sel.mol, res, altconf, 1,
                                                   atom_index_udd, udd_afix_handle);
      mmdb::Residue *res_copy = p.first;
      atom_selection_container_t residue_mol = p.second;

      // if is from SHELXL
      if (is_from_shelx_ins_flag)
         // if splitting after/with CA
         if (alt_conf_split_type == 0)
            // a member function?
            residue_mol = filter_atom_selection_container_CA_sidechain_only(residue_mol);

//       std::cout << "----------- state of residue mol: " << std::endl;
//       for (int i=0; i<residue_mol.n_selected_atoms; i++) {
//          std::cout << "residue mol " << residue_mol.atom_selection[i] << std::endl;
//       }

//       int udd_afix_handle_inter = residue_mol.mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");
//       std::cout << "DEBUG:: split_residue got udd_afix_handle_inter : "
//                 << udd_afix_handle_inter << std::endl;
//       for (int i=0; i<residue_mol.n_selected_atoms; i++) {
//          int afix_number = -1;
//          if (residue_mol.atom_selection[i]->GetUDData(udd_afix_handle_inter, afix_number) == mmdb::UDDATA_Ok)
//             std::cout << residue_mol.atom_selection[i] << " has afix number " << afix_number
//                       << std::endl;
//          else
//             std::cout << residue_mol.atom_selection[i]
//                       << " Failed get udd afix number in split_residue"
//                       << std::endl;
//       }


      // if we are splitting a water, then testing for alt_conf_split_type == 0
      // doen't make sense, so let's also skip that if we have a water.
      //
      short int is_water_flag = 0;
      std::string residue_type = res->GetResName();
      if (residue_type == "WAT" || residue_type == "HOH" || residue_type == "DUM")
         is_water_flag = 1;

      // We can't do a rotamer fit if we don't have C and N atoms:
      //
      // Also, we don't want to do a rotamer fit if the user has
      // turned this option off:
      //
      if (!graphics_info_t::show_alt_conf_intermediate_atoms_flag &&
          have_atoms_for_rotamer(res_copy)) {

         // We want to delete atoms if we chose alt_conf_split_type.. thing
         if (alt_conf_split_type == 0) { // partial split (not whole residue)
            mmdb::PPAtom residue_atoms = NULL;
            int nResidueAtoms;
            res_copy->GetAtomTable(residue_atoms, nResidueAtoms);
            std::string mol_atom_name;
            for (int iat=0; iat<nResidueAtoms; iat++) {
               mol_atom_name = std::string(residue_atoms[iat]->name);
               if (mol_atom_name == " N  " ||
                   mol_atom_name == " C  " ||
                   mol_atom_name == " H  " ||
                   mol_atom_name == " HA " || // CA hydrogen
                   mol_atom_name == " O  ") {
//                   std::cout << "DEBUG:: split residue deleting  "
//                             << mol_atom_name << std::endl;
                  res_copy->DeleteAtom(iat);
               } else {
                  // std::cout << "split residue accepting " << mol_atom_name
                  // << std::endl;
               }
            }
            res_copy->TrimAtomTable();
            // delete [] residue_atoms; // not with GetAtomTable(), fool.
         }

         // just go ahead and stuff in the atoms to the
         // molecule, no user intervention required.
         //
         // std::cout << "Calling  ------------ split_residue_then_rotamer path -------------\n";
         short int use_residue_mol_flag = 0;
         if (is_from_shelx_ins_flag)
            use_residue_mol_flag = 1;

         split_residue_then_rotamer(res_copy, altconf, residue_alt_confs, residue_mol,
                                    use_residue_mol_flag);

      } else {
         // we don't have all the atoms in the residue to do a
         // rotamer, therefore we must fall back to showing the
         // intermediate atoms.
         do_intermediate_atoms = 1;
         if (alt_conf_split_type == 0) { // partial split (not whole residue)

            if (! is_water_flag) {
               // so we need to delete some atoms from this residue:
               //
               // but we only want to do this if we are representing
               // intermediate atoms.  If we are going on to rotamers,
               // we need all the atoms (so that we can do a fit using
               // the main chain atoms of the rotamer).
               //
               //
               mmdb::PPAtom residue_atoms;
               int nResidueAtoms;
               res_copy->GetAtomTable(residue_atoms, nResidueAtoms);
               std::string mol_atom_name;
               for (int iat=0; iat<nResidueAtoms; iat++) {
                  mol_atom_name = std::string(residue_atoms[iat]->name);
                  if (mol_atom_name == " N  " ||
                      mol_atom_name == " C  " ||
                      mol_atom_name == " H  " ||
                      mol_atom_name == " HA " ||
                      mol_atom_name == " O  ") {
                     // std::cout << "split residue deleting  " << mol_atom_name
                     // << std::endl;
                     res_copy->DeleteAtom(iat);
                  } else {
                     // std::cout << "split residue accepting " << mol_atom_name
                     // << std::endl;
                  }
               }
            }
         } else {
            std::cout << "split_residue split type " << alt_conf_split_type
                      << " no deleting atoms  of this residue\n";
         }

         res_copy->TrimAtomTable();
         atom_selection_container_t asc_dummy;
         pr = split_residue_internal(res_copy, altconf, residue_alt_confs, p.second, 1);
      }
   } else {
      std::cout << "WARNING:: split_residue: bad atom index.\n";
   }
   std::cout << "split_residue(int atom_index, int alt_conf_split_type) returns "
             << pr.first << " " << pr.second << std::endl;
   return pr;
}

// change asc.
atom_selection_container_t
molecule_class_info_t::filter_atom_selection_container_CA_sidechain_only(atom_selection_container_t asc) const {

   std::string mol_atom_name;
   for (int iat=0; iat<asc.n_selected_atoms; iat++) {
      mol_atom_name = std::string(asc.atom_selection[iat]->name);
      if (mol_atom_name == " N  " ||
          mol_atom_name == " C  " ||
          mol_atom_name == " H  " ||
          // mol_atom_name == " HA " || // CA hydrogen we don't want
                                        // to delete it if we don't delete the CA.
          mol_atom_name == " H0 " || // shelx N hydrogen,
          mol_atom_name == " O  ") {
         mmdb::Residue *r = asc.atom_selection[iat]->residue;
         r->DeleteAtom(iat);
         // std::cout << "Filter out atom " << asc.atom_selection[iat] << std::endl;
      } else {
         // std::cout << "side chain atom " << asc.atom_selection[iat] << std::endl;
      }
   }
   asc.mol->FinishStructEdit();

   atom_selection_container_t ret_asc = make_asc(asc.mol);
   //    std::cout << "----------- immediate state of residue mol: " << std::endl;
   //    for (int i=0; i<ret_asc.n_selected_atoms; i++) {
   //       std::cout << "rest asc mol " << ret_asc.atom_selection[i] << std::endl;
   //    }

   return ret_asc;
}



std::vector<std::string>
molecule_class_info_t::get_residue_alt_confs(mmdb::Residue *res) const {

   std::vector<std::string> v = coot::util::get_residue_alt_confs(res);
   return v;
}

// What's in the residue     What we clicked   Old Coordinates   New Coordinates
//      ""                        ""                "" -> "A"       "B"
//    "A" "B"                     "A"              no change        "C"
//    "A" "B"                     "B"              no change        "C"
//    "" "A" "B"                  ""               [1]              "C"
//    "" "A" "B"                  "A"              [1]              "C"
//    "" "A" "B"                  "B"              [1]              "C"
//
// [1] depends on the split:
//     whole residue split: "" -> "A" , "A" and "B" remain the same
//     partial split:       no change
//
std::string
molecule_class_info_t::make_new_alt_conf(const std::vector<std::string> &residue_alt_confs,
                                         int alt_conf_split_type_in) const {

   std::string v("");
   std::vector<std::string> m;
   m.push_back("B");
   m.push_back("C");
   m.push_back("D");

   short int got;
   for (unsigned int im=0; im<m.size(); im++) {
      got = 0;
      for (unsigned int ir=0; ir<residue_alt_confs.size(); ir++) {
         if (m[im] == residue_alt_confs[ir]) {
            got = 1;
            break;
         }
      }
      if (got == 0) {
         v = m[im];
         break;
      }
   }

   return v;
}


// This is a molecule-class-info function.
//
// It calls do_accept_reject_dialog().  That is bad and should be cleaned up.
//
std::pair<bool, std::string>
molecule_class_info_t::split_residue_internal(mmdb::Residue *residue, const std::string &altconf,
                                              const std::vector<std::string> &all_residue_altconfs,
                                              atom_selection_container_t residue_mol,
                                              short int use_residue_mol_flag) {

   std::pair<bool,std::string> p(0,"");
   std::string ch(residue->GetChainID());

   // mmdb::Residue **SelResidues = new mmdb::Residue *; // memory leak
   mmdb::Residue **SelResidues = &residue; // just one

   // std::cout << "==================== in split_residue_internal ============= " << std::endl;

   atom_selection_container_t asc;
   if (!use_residue_mol_flag) {
      mmdb::Manager *mov_mol = create_mmdbmanager_from_res_selection(SelResidues,    // class function
                                                                    1, 0, 0, altconf,
                                                                    ch, 1);

      asc = make_asc(mov_mol);
   } else {
      asc = residue_mol;
      int udd_afix_handle = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");
//       std::cout << "DEBUG:: split_residue_internal got udd_afix_handle : "
//                 << udd_afix_handle << std::endl;
      for (int i=0; i<asc.n_selected_atoms; i++) {
         int afix_number = -1;
         if (asc.atom_selection[i]->GetUDData(udd_afix_handle, afix_number) == mmdb::UDDATA_Ok)
            std::cout << asc.atom_selection[i] << " has afix number " << afix_number << std::endl;
//          else
//             std::cout << asc.atom_selection[i]
//                       << " Failed get udd afix number in split residue_internal"
//                       << std::endl;
      }
   }

   std::string new_alt_conf = make_new_alt_conf(all_residue_altconfs,
                                                graphics_info_t::alt_conf_split_type);
   //std::cout << "DEBUG:: new_alt_conf " << new_alt_conf << " from ";
   // for (int ii=0; ii<all_residue_altconfs.size(); ii++) {
   // std::cout << "  :" << all_residue_altconfs[ii] << ":  ";
   // }
   std::cout << std::endl;
   p.first = 1;
   p.second = new_alt_conf;

   mmdb::Atom *at;
   for (int i=0; i<asc.n_selected_atoms; i++) {
      at = asc.atom_selection[i];
      at->x += 0.02;
      at->y += 0.2;
      at->z += 0.03;
      // apply the new_alt_conf:
      // std::cout << "applying alt_conf: " << new_alt_conf << std::endl;
      strncpy(at->altLoc, new_alt_conf.c_str(), 2);

      // set the occupancy:
      at->occupancy = graphics_info_t::add_alt_conf_new_atoms_occupancy;
      adjust_occupancy_other_residue_atoms(at, at->residue, 0);

//       at->SetAtomName(at->GetIndex(),
//                       at->serNum,
//                       at->name,
//                       new_alt_conf.c_str(),
//                       at->segID,
//                       at->element,
//                       at->charge);
   }

   // Rotamer?
   //
   // do_rotamers(0, imol);  no.  what is imol?
   //
   graphics_info_t g;

//    g.make_moving_atoms_graphics_object(asc);
//    g.imol_moving_atoms = imol_no;
//    g.moving_atoms_asc_type = coot::NEW_COORDS_INSERT_CHANGE_ALTCONF;

   g.set_moving_atoms(asc, imol_no, coot::NEW_COORDS_INSERT_CHANGE_ALTCONF);

// debugging
//    int nbonds = 0;
//    for (int i=0; i<regularize_object_bonds_box.num_colours; i++)
//       nbonds += regularize_object_bonds_box.bonds_[i].num_lines;
//    std::cout << "Post new bonds we have " << nbonds << " bonds lines\n";

   if (! graphics_info_t::show_alt_conf_intermediate_atoms_flag) {
      if (graphics_info_t::use_graphics_interface_flag)
         do_accept_reject_dialog("Alt Conf Split", coot::refinement_results_t());
   }
   return p;
}

// We don't create an intermediate atom.
// We create a molecule and do an accept moving atoms equivalent on them.
//
// What is residue here?  residue is a pure copy of the clicked on
// residue, including all altconfs.
//
// altconf is the altconf of the clicked atom
// all_altconfs are all the altconfs in that residue (used so that we
// can find a new altconf for the new atoms).
//
void
molecule_class_info_t::split_residue_then_rotamer(mmdb::Residue *residue, const std::string &altconf,
                                                  const std::vector<std::string> &all_residue_altconfs,
                                                  atom_selection_container_t residue_mol_asc,
                                                  short int use_residue_mol_flag) {


   mmdb::PResidue *SelResidues = nullptr;
   std::string ch(residue->GetChainID());

   // alt_conf_split_type is a graphics_info_t static data member
   //
   std::string new_altconf = make_new_alt_conf(all_residue_altconfs,
                                               graphics_info_t::alt_conf_split_type);

   // Move the atoms of residue a bit here?
   atom_selection_container_t mov_mol_asc;

//    std::cout << "DEBUG:: in split_residue_then_rotamer use_residue_mol_flag: " << use_residue_mol_flag
//              << std::endl;

   if (use_residue_mol_flag) {
      //       std::cout << "DEBUG:: in split_residue_then_rotamer shelxl path " << std::endl;
      mov_mol_asc = residue_mol_asc;
      int udd_afix_handle = residue_mol_asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");
      for (int i=0; i<residue_mol_asc.n_selected_atoms; i++) {
         int afix_number = -1;
         if (residue_mol_asc.atom_selection[i]->GetUDData(udd_afix_handle, afix_number) == mmdb::UDDATA_Ok)
            std::cout << residue_mol_asc.atom_selection[i] << " has afix number " << afix_number
                      << std::endl;
//          else
//             std::cout << residue_mol_asc.atom_selection[i]
//                       << " Failed get udd afix number in split residue_internal"
//                       << std::endl;
      }

   } else {
      // std::cout << "DEBUG:: in split_residue_then_rotamer normal path " << std::endl;
      SelResidues = &residue; // just one
      mmdb::Manager *mov_mol = create_mmdbmanager_from_res_selection(SelResidues,
                                                                    1, 0, 0, altconf,
                                                                    ch, 1);
      mov_mol_asc = make_asc(mov_mol);
   }

   for (int i=0; i<mov_mol_asc.n_selected_atoms; i++) {
      mmdb::Atom *at = mov_mol_asc.atom_selection[i];
      at->x += 0.1;
      strncpy(at->altLoc, new_altconf.c_str(), 2);
      // set the occupancy:
      at->occupancy = graphics_info_t::add_alt_conf_new_atoms_occupancy;
   }

   std::string at_name;
   if (mov_mol_asc.n_selected_atoms > 0) {
      at_name = mov_mol_asc.atom_selection[0]->name;
   }

   insert_coords_change_altconf(mov_mol_asc);

   // now find the atom index of an atom with the new alt conf:

   int resno = residue->GetSeqNum();
   std::string chain_id = residue->GetChainID();
   std::string ins_code = residue->GetInsCode();

   int atom_index = full_atom_spec_to_atom_index(chain_id,
                                                 resno,
                                                 ins_code,
                                                 at_name,
                                                 new_altconf);

   if (atom_index >= 0) {

      // std::cout << "DEBUG:: in split_residue_then_rotamer() " << std::endl;
      graphics_info_t g;
      // Use the other do_rotamers() function, because we can pass the atom, not the atom index
      // g.do_rotamers(atom_index, imol_no); // this imol, obviously.
      mmdb:: Atom *aa = atom_sel.atom_selection[atom_index];
      g.do_rotamers(imol_no, aa); // this imol, obviously.
   } else {
      std::cout << "ERROR bad atom index in split_residue_then_rotamer: "
                << atom_index << std::endl;
   }
}

// We just added a new atom to a residue, now we need to adjust the
// occupancy of the other atoms (so that we don't get residues with
// atoms whose occupancy is greater than 1.0 (Care for SHELX molecule?)).
// at doesn't have to be in residue.
//
// Perhaps this can be a utility function?
//
void
molecule_class_info_t::adjust_occupancy_other_residue_atoms(mmdb::Atom *at,
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


short int
molecule_class_info_t::have_atoms_for_rotamer(mmdb::Residue *res) const {

   short int ihave = 0;  // initially not
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   int n_mainchain = 0;
   res->GetAtomTable(residue_atoms, nResidueAtoms);
   short int have_c = 0;
   short int have_ca = 0;
   short int have_n = 0;
   for (int iat=0; iat<nResidueAtoms; iat++) {
      std::string at_name(residue_atoms[iat]->name);
      if (at_name == " C  ") {
         n_mainchain++;
         have_c = 1;
      }
      if (at_name == " CA ") {
         n_mainchain++;
         have_ca = 1;
      }
      if (at_name == " N  ") {
         n_mainchain++;
         have_n = 1;
      }

   }
   if ((n_mainchain > 2) && have_c && have_ca && have_n)
      ihave = 1;

   return ihave;
}


// The flanking residues (if any) are in the residue selection (SelResidues).
// The flags are not needed now we have made adjustments in the calling
// function.
//
// create_mmdbmanager_from_res_selection must make adjustments
//
mmdb::Manager *
molecule_class_info_t::create_mmdbmanager_from_res_selection(mmdb::PResidue *SelResidues,
                                                             int nSelResidues,
                                                             int have_flanking_residue_at_start,
                                                             int have_flanking_residue_at_end,
                                                             const std::string &altconf,
                                                             const std::string &chain_id_1,
                                                             short int residue_from_alt_conf_split_flag) {

   int start_offset = 0;
   int end_offset = 0;

//    if (have_flanking_residue_at_start)
//       start_offset = -1;
//    if (have_flanking_residue_at_end)
//       end_offset = +1;

   mmdb::Manager *residues_mol = new mmdb::Manager;
   mmdb::Model *model = new mmdb::Model;
   mmdb::Chain *chain = new mmdb::Chain;
   bool whole_res_flag = 0; // not all alt confs, only this one ("A") and "".

   // For the active residue range (i.e. not the flanking residues) we only want
   // to refine the atoms that have the alt conf the same as the picked atom
   // (and that is altconf, passed here).
   //
   // However, for *flanking residues* it's different.  Say we are refining a
   // non-split residue with alt conf "": Say that residue has a flanking
   // residue that is completely split, into A and B.  In that case we want
   // either "" or "A" for the flanking atoms.
   //
   // And say we want to refine the A alt conf of a completely split residue
   // that has a flanking neighbour that is completely unsplit (""), we want
   // atoms that are either "A" or "".
   //
   // So let's try setting whole_res_flag to 1 for flanking residues.

   mmdb::Residue *r;
   for (int ires=start_offset; ires<(nSelResidues + end_offset); ires++) {

      if ( (ires == 0) || (ires == nSelResidues -1) ) {
         if (! residue_from_alt_conf_split_flag)
            whole_res_flag = 1;
      } else {
         whole_res_flag = 0;
      }

      r = coot::util::deep_copy_this_residue_add_chain(SelResidues[ires], altconf,
                                                       whole_res_flag, 0);
      chain->AddResidue(r);
      r->seqNum = SelResidues[ires]->GetSeqNum();
      r->SetResName(SelResidues[ires]->GetResName());
   }
   chain->SetChainID(chain_id_1.c_str());
   model->AddChain(chain);
   residues_mol->AddModel(model);
   residues_mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   residues_mol->FinishStructEdit();

   return residues_mol;
}


// merge molecules
// Return +1 as status of pair if we did indeed do a merge
// Return a list of new chain ids as the second.
//
// Note, very often add_molecules will only be of size 1.
//
// Recall that we will often be merging ligands into (otherwise quite
// complete) proteins.  In that case, in which we add a single residue
// it would be The Right Thing to Do if we could find a chain that
// consisted only of the same residue type as is the new ligand and
// add it to that chain. If that is the case, we should return a spec
// for the residue, not just the chain id.
//
// If we can't find a chain like that - or the new molecule contains
// more than one residue, we simply add new chains (with new chain
// ids) to this molecule.
//
// Question to self: how do I deal with different models?
//
std::pair<int, std::vector<merge_molecule_results_info_t> >
molecule_class_info_t::merge_molecules(const std::vector<atom_selection_container_t> &add_molecules) {

   int istat = 0;
   make_backup(); // could be more clever, by doing this only when needed.
   std::vector<merge_molecule_results_info_t> resulting_merge_info;
   std::pair<bool, coot::residue_spec_t> done_merge_ligand_to_near_chain;
   done_merge_ligand_to_near_chain.first = false;

   std::vector<std::string> this_model_chains = coot::util::chains_in_molecule(atom_sel.mol);

   for (unsigned int imol=0; imol<add_molecules.size(); imol++) {
      mmdb::Manager *adding_mol = add_molecules[imol].mol;
      if (add_molecules[imol].n_selected_atoms > 0) {
         int nresidues = coot::util::number_of_residues_in_molecule(add_molecules[imol].mol);

         // We need to set multi_residue_add_flag appropriately.  We
         // unset it if the molecule to be added has only one residue
         //
         bool multi_residue_add_flag = true;

         std::vector<std::string> adding_model_chains
            = coot::util::chains_in_molecule(add_molecules[imol].mol);

         if (nresidues == 1) {

            // pass this?
            const coot::residue_spec_t &spec = graphics_info_t::merge_molecules_ligand_spec;
            bool done_add_specific = merge_molecules_just_one_residue_at_given_spec(add_molecules[imol], spec);
            bool done_homogeneous_addition_flag = false;

            // by "homogeneous" I mean is there a chain of residues of the same type as that we are adding
            // e.g. an SO4 to a chain of SO4s? This don't happen (much?) these days.
            // These days, with one residue, we expect to run merge_ligand_to_near_chain()

            if (! done_add_specific)
               done_homogeneous_addition_flag = merge_molecules_just_one_residue_homogeneous(add_molecules[imol]);

            if (done_add_specific)
               multi_residue_add_flag = false;
            else
               multi_residue_add_flag = ! done_homogeneous_addition_flag;

            if (! done_homogeneous_addition_flag) {

               if (! done_add_specific)
                  done_merge_ligand_to_near_chain = merge_ligand_to_near_chain(adding_mol);

               if (done_merge_ligand_to_near_chain.first) {
                  merge_molecule_results_info_t mmr;
                  mmr.is_chain = false;
                  mmr.spec = done_merge_ligand_to_near_chain.second;
                  resulting_merge_info.push_back(mmr);
                  istat = 1;
               } else {

                  if (done_add_specific) {
                     // JED ligand addition
                     merge_molecule_results_info_t mmr;
                     mmr.is_chain = false;
                     mmr.spec = spec;
                     // std::cout << "---- JED case pushing back mmr " << mmr.spec << std::endl;
                     resulting_merge_info.push_back(mmr);
                     istat = 1;
                  }
               }

               // set multi_residue_add_flag if that was not a successful merge
               if (! done_add_specific)
                  multi_residue_add_flag = ! done_merge_ligand_to_near_chain.first;

            }
         }

         // Now that multi_residue_add_flag has been set properly, we use it...

         if (multi_residue_add_flag) {
            // return state
            std::pair<bool, std::vector<std::string> > add_state = try_add_by_consolidation(adding_mol);

            // some mild hacking, but we need to return a proper state and the added chains
            istat = 0;
            for (unsigned int i=0; i<add_state.second.size(); i++) {
               merge_molecule_results_info_t mmr;
               mmr.is_chain = true;
               mmr.chain_id = add_state.second[i];
               resulting_merge_info.push_back(mmr);
            }
            if (add_state.first) {
               update_molecule_after_additions();
               if (graphics_info_t::show_symmetry == 1)
                  update_symmetry();
               multi_residue_add_flag = false; // we've added everything for this mol.
               istat = add_state.first;
            }
         }

         // this should happen rarely these days...
         //
         if (multi_residue_add_flag) {

            std::vector<std::string> mapped_chains =
               map_chains_to_new_chains(adding_model_chains, this_model_chains);

            std::cout << "INFO:: Merge From chains: " << std::endl;
            for (unsigned int ich=0; ich<adding_model_chains.size(); ich++)
               std::cout << " :" << adding_model_chains[ich] << ":";
            std::cout << std::endl;
            std::cout << "INFO:: Merge To chains: " << std::endl;
            for (unsigned int ich=0; ich<mapped_chains.size(); ich++)
               std::cout << " :" << mapped_chains[ich] << ":";
            std::cout << std::endl;

            if (mapped_chains.size() != adding_model_chains.size()) {
               // can't continue with merging - no chains available.
               std::cout << "can't continue with merging - no chains available." << std::endl;
            } else {
               // fine, continue

               // Add the chains of the new molecule to this atom_sel, chain by chain.
               int i_add_model = 1;
               int i_this_model = 1;

               mmdb::Model *model_p = add_molecules[imol].mol->GetModel(i_add_model);
               mmdb::Model *this_model_p = atom_sel.mol->GetModel(i_this_model);

               int n_add_chains = model_p->GetNumberOfChains();

               for (int iaddchain=0; iaddchain<n_add_chains; iaddchain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(iaddchain);
                  mmdb::Chain *copy_chain_p = new mmdb::Chain;
                  copy_chain_p->Copy(chain_p);
                  copy_chain_p->SetChainID(mapped_chains[iaddchain].c_str());
                  this_model_p->AddChain(copy_chain_p);
                  this_model_chains.push_back(mapped_chains[iaddchain].c_str());
                  merge_molecule_results_info_t mmr;
                  mmr.is_chain = true;
                  mmr.chain_id = mapped_chains[iaddchain];
                  resulting_merge_info.push_back(mmr);
               }

               if (n_add_chains > 0) {
                  atom_sel.mol->FinishStructEdit();
                  update_molecule_after_additions();
                  if (graphics_info_t::show_symmetry == 1)
                     update_symmetry();
               }
               istat = 1;
            }
         }
      }
   }

   fill_ghost_info(true, 0.7);

   // std::cout << "------- resulting_merge_info has size " << resulting_merge_info.size() << std::endl;
   if (!resulting_merge_info.empty())
      std::cout << "INFO:: in merge_molecules(): resulting_merge_info[0] " << resulting_merge_info[0].spec << std::endl;

   return std::pair<int, std::vector<merge_molecule_results_info_t> > (istat, resulting_merge_info);
}

// return the multi_residue_add_flag
// done_homogeneous_addition_flag
//
bool
molecule_class_info_t::merge_molecules_just_one_residue_homogeneous(atom_selection_container_t molecule_to_add) {

   // If there is a chain that has only residues of the same
   // type as is the (single) residue in the new adding
   // molecule then we add by residue addition to chain
   // rather than add by chain (to molecule).

   // Are there chains in this model that only consist of
   // residues of type adding_model_chains[0]?
   //
   bool done_homogeneous_addition_flag = false;

   bool has_single_residue_type_chain_flag = false;

   int i_this_model = 1;

   mmdb::Chain *add_residue_to_this_chain = NULL;

   mmdb::Model *this_model_p = atom_sel.mol->GetModel(i_this_model);

   int n_this_mol_chains = this_model_p->GetNumberOfChains();

   for (int ithischain=0; ithischain<n_this_mol_chains; ithischain++) {
      mmdb::Chain *this_chain_p = this_model_p->GetChain(ithischain);
      std::vector<std::string> r = coot::util::residue_types_in_chain(this_chain_p);

      if (r.size() == 1) {
         std::string adding_model_resname(molecule_to_add.atom_selection[0]->residue->GetResName());
         if (r[0] == adding_model_resname) {
            // poly-ala helices (say) should not go into concatenated residues in the same chain
            if (adding_model_resname != "ALA") {
               add_residue_to_this_chain = this_chain_p;
               has_single_residue_type_chain_flag = true;
               break;
            }
         }
      }
   }

   if (has_single_residue_type_chain_flag) {
      if (molecule_to_add.n_selected_atoms > 0) {
         mmdb::Residue *add_model_residue = molecule_to_add.atom_selection[0]->residue;
         copy_and_add_residue_to_chain(add_residue_to_this_chain, add_model_residue);
         done_homogeneous_addition_flag = true;
         atom_sel.mol->FinishStructEdit();
         update_molecule_after_additions();
         if (graphics_info_t::show_symmetry == 1)
            update_symmetry();
      }
   }
   return done_homogeneous_addition_flag;
}

bool
molecule_class_info_t::merge_molecules_just_one_residue_at_given_spec(atom_selection_container_t molecule_to_add,
                                                                      coot::residue_spec_t target_spec) {

   bool status = false;

   if (! target_spec.empty()) {
      mmdb::Residue *residue_p = get_residue(target_spec);
      if (! residue_p) {
         // more checks: does molecule_to_add have only one residue?
         int i_model = 1;

         int n_res = coot::util::number_of_residues_in_molecule(molecule_to_add.mol);

         if (n_res == 1) {
            mmdb::Model *this_model_p = atom_sel.mol->GetModel(i_model);
            mmdb::Chain *this_chain_p = this_model_p->GetChain(target_spec.chain_id.c_str());
            if (! this_chain_p) {
               this_chain_p = new mmdb::Chain;
               this_chain_p->SetChainID(target_spec.chain_id.c_str());
               this_model_p->AddChain(this_chain_p);
            } else {
               std::cout << "INFO:: merge_molecules_just_one_residue_at_given_spec() "
                         << " this chain not found in molecule (good)" << std::endl;
            }
            mmdb::Residue *r = coot::util::get_first_residue(molecule_to_add.mol);
            if (r) {
               make_backup();
               mmdb::Residue *new_residue_p = copy_and_add_residue_to_chain(this_chain_p, r);
               new_residue_p->seqNum = target_spec.res_no;
               status = true;
            }
         } else {
            if (true) // debug
               std::cout << "debug:: merge_molecules_just_one_residue_at_given_spec() oops "
                         << " n_res is " << n_res << std::endl;
         }
      } else {
         std::cout << "WARNING:: merge_molecules_just_one_residue_at_given_spec() residue already exists "
                   << "in molecule " << target_spec << std::endl;
      }
   } else {
      if (false) // debug
         std::cout << "merge_molecules_just_one_residue_at_given_spec() null residue spec" << std::endl;
   }

   if (status) {
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
      if (graphics_info_t::show_symmetry == 1)
         update_symmetry();
   }

   if (false) // debug
      std::cout << "merge_molecules_just_one_residue_at_given_spec() returns " << status << std::endl;

   return status;
}


// There is (or should be) only one residue in mol
std::pair<bool, coot::residue_spec_t>
molecule_class_info_t::merge_ligand_to_near_chain(mmdb::Manager *mol) {

   bool done_merge = false;
   coot::residue_spec_t res_spec;

   mmdb::Residue *adding_residue_p = 0;

   { // set adding_residue
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            if (chain_p) {
               int nres = chain_p->GetNumberOfResidues();
               if (nres > 0)
                  adding_residue_p = chain_p->GetResidue(0);
            }
         }
      }
   }

   if (adding_residue_p) {

      // OK, what atoms in this molecule are close to adding_residue?
      // First, we need a vector of atoms in adding_residue
      // Then a set of atoms that are close to those positions
      // Then find the chain that has most atoms in that set
      // Then find a suitable residue number
      // Then add the residue

      std::vector<clipper::Coord_orth> ligand_atom_positions;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      adding_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
         mmdb::Atom *atom = residue_atoms[i];
         if (!atom->isTer()) {
            clipper::Coord_orth pt = coot::co(atom);
            ligand_atom_positions.push_back(pt);
         }
      }
      std::set<mmdb::Atom *> near_atoms; // atoms in the protein (this molecule)
      double dist_crit = 4.2;
      double dist_crit_sqrd = dist_crit * dist_crit;

      // this is a slow position by position check. It can be speeded up, if needed,
      // by fiddling with residue and molecules (i.e. copying and merging) and using
      // SelectAtoms().
      //
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (chain_p) {
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (at) {
                           if (! at->isTer()) {
                              clipper::Coord_orth this_at_co = coot::co(at);
                              for (std::size_t j=0; j<ligand_atom_positions.size(); j++) {
                                 const clipper::Coord_orth &pos = ligand_atom_positions[j];
                                 double this_dist_sqrd = (this_at_co-pos).lengthsq();
                                 if (this_dist_sqrd < dist_crit_sqrd) {
                                    near_atoms.insert(at);
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
         if (near_atoms.size() > 1) {
            // make a map of the number of residue in each chain that are close to the ligand atoms
            std::map<mmdb::Chain *, int> chain_map;
            std::set<mmdb::Atom *>::const_iterator it;
            for (it=near_atoms.begin(); it!=near_atoms.end(); it++) {
               mmdb::Chain *this_chain = (*it)->GetChain();
               if (chain_map.find(this_chain) == chain_map.end()) {
                  chain_map[this_chain] = 1;
               } else {
                  chain_map[this_chain]++;
               }
            }

            // find the "best" chain for this added residue/ligand
            mmdb::Chain *max_atom_chain = 0;
            int n_atoms_max = 0;
            std::map<mmdb::Chain *, int>::const_iterator chain_it;
            for (chain_it=chain_map.begin(); chain_it!=chain_map.end(); chain_it++) {
               int n = chain_it->second;
               if (n > n_atoms_max) {
                  n_atoms_max = n;
                  max_atom_chain = chain_it->first;
               }
            }

            if (n_atoms_max > 0) {
               if (max_atom_chain) {
                   // does not backup or finalize or update
                   bool new_res_no_by_hundreds = true;
                   mmdb::Residue *new_residue_p =
                      copy_and_add_residue_to_chain(max_atom_chain, adding_residue_p, new_res_no_by_hundreds);
                   done_merge = true;
                   res_spec = coot::residue_spec_t(new_residue_p);
               }
            }
         }
      }
   }

   if (done_merge) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
      if (graphics_info_t::show_symmetry == 1)
         update_symmetry();
   }
   return std::pair<bool, coot::residue_spec_t> (done_merge, res_spec);
}


// return status and vector of resulting chain ids.
//
std::pair<bool, std::vector<std::string> >
molecule_class_info_t::try_add_by_consolidation(mmdb::Manager *adding_mol) {

   bool status = false;
   std::vector<std::string> chain_ids;

   // 20180104 Don't merge molecules made of ALA.

   // for this molecule molecule, make a map of chains that have one
   // residue type.  Could well be empty (or perhaps consist of just a
   // water chain)
   //
   std::map<std::string, std::pair<int, mmdb::Chain *> > single_res_type_map;
   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         std::vector<std::string> residue_types;
         mmdb::Residue *residue_p;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            std::string res_name(residue_p->GetResName());
            if (std::find(residue_types.begin(), residue_types.end(), res_name) == residue_types.end())
               residue_types.push_back(res_name);
            if (residue_types.size() > 1)
               break;
         }
         if (residue_types.size() == 1)
            single_res_type_map[residue_types[0]] = std::pair<int, mmdb::Chain *> (imod, chain_p);
      }
   }

   for(int imod = 1; imod<=adding_mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = adding_mol->GetModel(imod);
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         bool done_this_chain = false;
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         std::vector<std::string> residue_types;
         mmdb::Residue *residue_p;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            std::string res_name(residue_p->GetResName());
            if (std::find(residue_types.begin(), residue_types.end(), res_name) == residue_types.end())
               residue_types.push_back(res_name);
         }

         if (residue_types.size() == 1) {
            if (residue_types[0] != "ALA") {
               std::map<std::string, std::pair<int, mmdb::Chain *> >::const_iterator it =
                  single_res_type_map.find(residue_types[0]);
               if (it != single_res_type_map.end()) {
                  if (it->second.first == imod) {

                     // We got a match, now add all of adding_mol chain_p
                     // to this molecule's chain

                     // BL says:: we check in copy_and_add_chain_residues_to_chain if there
                     // are overlapping waters. Alternativley we could do it here already.

                     copy_and_add_chain_residues_to_chain(chain_p, it->second.second);
                     done_this_chain = true;
                     std::string cid = it->second.second->GetChainID();
                     if (std::find(chain_ids.begin(), chain_ids.end(), cid) == chain_ids.end())
                        chain_ids.push_back(cid);
                  }
               }
            }
         }

         if (! done_this_chain) {
            // copy whole chain to a new chain
            mmdb::Model *this_model_p = atom_sel.mol->GetModel(imod);
            if (this_model_p) {
               std::string current_chain_id = chain_p->GetChainID();
               std::string new_chain_id = suggest_new_chain_id(current_chain_id);
               mmdb::Chain *copy_chain_p = new mmdb::Chain;
               copy_chain_p->Copy(chain_p);
               copy_chain_p->SetChainID(new_chain_id.c_str());
               this_model_p->AddChain(copy_chain_p);
               if (std::find(chain_ids.begin(), chain_ids.end(), new_chain_id) == chain_ids.end())
                  chain_ids.push_back(new_chain_id);
            }
         }
         atom_sel.mol->FinishStructEdit();
         status = true;
      }
   }
   return std::pair<bool, std::vector<std::string> > (status, chain_ids);
}

// Copy residues of new_chain into this_model_chain
void
molecule_class_info_t::copy_and_add_chain_residues_to_chain(mmdb::Chain *new_chain, mmdb::Chain *this_molecule_chain) {

   // remove TER record of current last residue (if it has a TER).
   remove_TER_on_last_residue(this_molecule_chain);

   int nres = new_chain->GetNumberOfResidues();
   for (int ires=0; ires<nres; ires++) {
      copy_and_add_residue_to_chain(this_molecule_chain, new_chain->GetResidue(ires));
   }
}


// Merge molecules helper function.
//
// What we want to do is map chains ids in the new (adding) molecules
// to chain ids that are unused in thid model:
//
// e.g. adding:  B,C,D
//      this:    A,B,C
//
//      return:  D,E,F
//
std::vector<std::string>
molecule_class_info_t::map_chains_to_new_chains(const std::vector<std::string> &adding_model_chains,
                                                const std::vector<std::string> &this_model_chains) const {

   std::vector<std::string> rv;
   // we dont want ! [] or *, they are mmdb specials.
   // 20200110-PE no more specials at all
   // std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz#$%^&@?/~|-+=(){}:;.,'");
   std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

   // first remove from r, all the chain that already exist in this molecule.
   for (unsigned int i=0; i<this_model_chains.size(); i++) {
      // but only do that if the chain id is of length 1 (single character)
      if (this_model_chains[i].length() == 1) {
         std::string::size_type found_this_model_chain = r.find(this_model_chains[i]);
         if (found_this_model_chain != std::string::npos) {
            // there was a match
            r = coot::util::remove_string(r, this_model_chains[i]);
         } else {
            // else there was not a match, this chain id does not
            // exist in r (surprisingly).
         }
      }
   }

   for (unsigned int i=0; i<adding_model_chains.size(); i++) {

      std::string t = "A";
      std::cout << "finding new chain id for chain id :" << adding_model_chains[i] << ": "
                << i << "/" << adding_model_chains.size() << std::endl;

      if (r.length() > 0) {
         t = r[0];
         r = r.substr(1); // return r starting at position 1;
      } else {
         t = "A";
      }
      rv.push_back(t);
   }
   return rv;
}

// return "" on failure
std::string
molecule_class_info_t::suggest_new_chain_id(const std::string &current_chain_id) const {

   // current_chain_id is the chain_id in the molecule that we are adding to this one.

   std::string new_chain_id;

   // 20200110-PE no more specials at all
   // std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz#$%^&@?/~|-+=(){}:;.,'");
   std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

   std::vector<std::string> existing;
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   // do we have current_chain_id in our list of chains already?  If not, then simply
   // choose current_chain_id.
   bool found_it = false;
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      std::string chid = chain_p->GetChainID();
      if (chid == current_chain_id) {
         found_it = true;
         break;
      }
   }
   if (! found_it)
      new_chain_id = current_chain_id; // all done!

   // how about a multichar post-fix? We only want to do that if current_chain_id
   // what multichar to begin with (from a pdbx file):
   // (Wolfram Tempel)
   //
   if (new_chain_id.empty()) {
      if (current_chain_id.length() > 1) {
         std::string trial_chain_id = current_chain_id + "2";
         if (trial_chain_id.length() < (10-1)) { // magic mmdb chain id max length (mmdb_defs.h)
                                               // (we need a space for the terminal NULL too).
            bool found_it = false;
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               std::string chid = chain_p->GetChainID();
               if (chid == trial_chain_id) {
                  found_it = true;
                  break;
               }
            }
            if (! found_it)
               new_chain_id = trial_chain_id; // all done!
         }
      }
   }

   if (new_chain_id.empty()) { // not set yet
      for (int ichain=0; ichain<n_chains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         existing.push_back(chain_p->GetChainID());
      }
      unsigned int l = r.length();
      std::vector<std::string> candidates(l);
      for (unsigned int i=0; i<l; i++)
         candidates[i] = r[i];

      for (unsigned int i=0; i<existing.size(); i++)
         candidates.erase(std::remove(candidates.begin(), candidates.end(), existing[i]), candidates.end());

      if (candidates.size())
         new_chain_id = candidates[0];
   }
   return new_chain_id;
}

// This doesn't do a backup or finalise model.
mmdb::Residue *
molecule_class_info_t::copy_and_add_residue_to_chain(mmdb::Chain *this_model_chain,
                                                     mmdb::Residue *add_model_residue,
                                                     bool new_resno_by_hundreds_flag) {

   mmdb::Residue *res_copied = NULL;
   if (add_model_residue) {
      bool whole_res_flag = true;
      int udd_atom_index_handle = 1; // does this matter?
      bool add_this = true;
      // check for overlapping water (could be generalised for same residue type?!
      std::vector<mmdb::Residue *> close_residues;
      close_residues = coot::residues_near_residue(add_model_residue, atom_sel.mol, 0.05);
      for (unsigned int i=0; i<close_residues.size(); i++) {
         if (close_residues[i]->isSolvent() && add_model_residue->isSolvent()) {
            add_this = false;
            std::cout<<"INFO:: not adding water because of overlap\n"<<std::endl;
            break;
         }
      }
      if (add_this) {

         /* No - this does an implicit embed-in-chain - that is not what we want
         mmdb::Residue *residue_copy = coot::deep_copy_this_residue(add_model_residue,
                                                                    "",
                                                                    whole_res_flag,
                                                                    udd_atom_index_handle);
         */
         mmdb::Residue *residue_copy = coot::util::deep_copy_this_residue(add_model_residue);

         if (residue_copy) {
            std::pair<short int, int> res_info =
               next_residue_number_in_chain(this_model_chain, new_resno_by_hundreds_flag);
            int new_res_resno = 9999;
            if (res_info.first)
               new_res_resno = res_info.second;
            residue_copy->seqNum = new_res_resno; // try changing the seqNum before AddResidue().
            this_model_chain->AddResidue(residue_copy);
            res_copied = residue_copy;
         }
      }
   }
   return res_copied;
}


int
molecule_class_info_t::renumber_residue_range(const std::string &chain_id,
                                              int start_resno, int last_resno, int offset) {

   int status = 0;

   // PDBCleanup(mmdb::PDBCLEAN_SERIAL) doesn't move the residue to
   // the end of the chain when we change the seqNum. Boo.
   // So, let's make a vector of residues that we will add
   // (insert) after the original residues have been deleted.
   //
   // BL says:: I dont belive it, so lets try...
   // OK, seems to be the case, so just add a SortResidues and
   // Bob is your uncle.
   // BL says (30/6/15):: should check if we get overlapping residues,
   // if so then dont do it.

   if (start_resno > last_resno) {
      int tmp = start_resno;
      start_resno = last_resno;
      last_resno = tmp;
   }

   // check existing residues (except moving ones) in range we want to move to
   int residue_exists = 0;
   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(i_chain);
	    std::string mol_chain(chain_p->GetChainID());
	    if (mol_chain == chain_id) {
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::Residue *residue_p;
	       for (int ires=0; ires<nres; ires++) { // ires is a serial number
		  residue_p = chain_p->GetResidue(ires);
		  int res_no = residue_p->seqNum;
		  char *ins_code = residue_p->GetInsCode();
		  if (res_no >= start_resno) {
		     if (res_no <= last_resno) {
			int new_res_no = res_no + offset;
			// moving range, so check for overlap in non-moving range
			if ((new_res_no < start_resno) || (new_res_no > last_resno)) {
			   residue_exists = does_residue_exist_p(chain_p->GetChainID(), new_res_no, ins_code);
			   if (residue_exists)
			      break;
			}
		     }
		  }
               }
            }
         }
      }
   }


   if (!residue_exists) {
      if (atom_sel.n_selected_atoms > 0) {
         mmdb::Model *model_p = atom_sel.mol->GetModel(1);
	 if (model_p) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int i_chain=0; i_chain<n_chains; i_chain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(i_chain);
	       std::string mol_chain(chain_p->GetChainID());
	       if (mol_chain == chain_id) {

		  make_backup();
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) { // ires is a serial number
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     if (residue_p->seqNum >= start_resno) {
			if (residue_p->seqNum <= last_resno) {
			   coot::residue_spec_t old_res_spec(residue_p);
			   coot::residue_spec_t new_res_spec(residue_p); // adjustment needed
			   new_res_spec.res_no += offset;

			   residue_p->seqNum += offset;
			   status = 1; // found one residue at least.

			   update_any_link_containing_residue(old_res_spec, new_res_spec);
			}
		     }
		  }
	       }
	       if (status == 1)
		  chain_p->SortResidues();
	    }
	 }
      }
      if (status == 1) {
         have_unsaved_changes_flag = 1;
         atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
         atom_sel.mol->FinishStructEdit();
         update_molecule_after_additions();
      }
   } else {
      std::cout << "WARNING:: the new residue range overlaps with original one. "
                << "Please change the range. Nothing has been done." << std::endl;
   }
   return status;

}


int
molecule_class_info_t::renumber_residue_range_old(const std::string &chain_id,
                                              int start_resno, int last_resno, int offset) {

   int status = 0;

   // PDBCleanup(mmdb::PDBCLEAN_SERIAL) doesn't move the residue to
   // the end of the chain when we change the seqNum. Boo.
   // So, let's make a vector of residues that we will add
   // (insert) after the original residues have been deleted.
   //
   std::vector<mmdb::Residue *> renumbered_residues;
   std::vector<mmdb::Residue *> residues_to_be_deleted;
   mmdb::Chain *chain_p_active = NULL; // the renumbered residues are in this chain


   if (start_resno > last_resno) {
      int tmp = start_resno;
      start_resno = last_resno;
      last_resno = tmp;
   }

   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(i_chain);
	    std::string mol_chain(chain_p->GetChainID());
	    if (mol_chain == chain_id) {

	       chain_p_active = chain_p;
	       make_backup();
	       int nres = chain_p->GetNumberOfResidues();
	       for (int ires=0; ires<nres; ires++) { // ires is a serial number
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		  if (residue_p->seqNum >= start_resno) {
		     if (residue_p->seqNum <= last_resno) {

			mmdb::Residue *residue_copy = coot::util::deep_copy_this_residue(residue_p);
			renumbered_residues.push_back(residue_copy);
			residues_to_be_deleted.push_back(residue_p);

			residue_copy->seqNum += offset;
			status = true; // found one residue at least.

		     }
		  }
	       }
	    }
         }
      }
   }
   if (status) {
      have_unsaved_changes_flag = 1;

      for (unsigned int ires=0; ires<residues_to_be_deleted.size(); ires++) {
         delete residues_to_be_deleted[ires];
         residues_to_be_deleted[ires] = NULL;
      }

      for (unsigned int ires=0; ires<renumbered_residues.size(); ires++) {

         int rr_seq_num = renumbered_residues[ires]->GetSeqNum();
         std::string rr_ins_code = renumbered_residues[ires]->GetInsCode();

         mmdb::Residue *res_p = renumbered_residues[ires]; // short-hand

         // What is the current serial number of the residue that
         // should be immediately after residue to be inserted?
         int iser = -1; // unset

         int chain_n_residues = chain_p_active->GetNumberOfResidues();
         for (int ichres=0; ichres<chain_n_residues; ichres++) {
            mmdb::Residue *ch_res = chain_p_active->GetResidue(ichres);
            int ch_res_seq_num = ch_res->GetSeqNum();
            std::string ch_res_ins_code = ch_res->GetInsCode();
            if (ch_res_seq_num == rr_seq_num) {
               if (ch_res_ins_code > rr_ins_code) {
                  iser = ichres;
                  break;
               }
            }
            if (ch_res_seq_num > rr_seq_num) {
               iser = ichres;
               break;
            }
         }

         // std::cout << "iser: " << iser << " for new residue seqnum " << rr_seq_num << std::endl;

         // now use iser if it was set.
         if (iser >= 0) {
            chain_p_active->InsResidue(res_p, iser);
         } else {
            // there was no residue that was that should be
            // immediately after the residue to be inserted.
            // Therefore, insert the residue at the end -
            // i.e. AddResidue.
            // std::cout << "Adding Residue to active chain " << std::endl;
            chain_p_active->AddResidue(res_p);
         }
      }

      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL);
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
   }
   return status;
}

int
molecule_class_info_t::change_residue_number(const std::string &chain_id,
                                             int current_resno,
                                             const std::string &current_inscode_str,
                                             int new_resno,
                                             const std::string &new_inscode_str) {


   if (false) // debug
      std::cout << "debug:: change_residue_number() called with chain_id " << chain_id
                << " current_resno, ins " << current_resno << "\"" << current_inscode_str << "\""
                << " new_resno, ins " << new_resno << "\"" << new_inscode_str << "\""
                << std::endl;

   int done_it = 0;
   mmdb::Residue *current_residue_p = get_residue(chain_id, current_resno, current_inscode_str);
   if (current_residue_p) {
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *this_chain_p = model_p->GetChain(ichain);
            if (this_chain_p == current_residue_p->GetChain()) {

               make_backup();

               std::string alt_conf = ""; // loop over alt_confs if not lazy
               bool whole_res_flag = true;
               int atom_index_udd_handle = atom_sel.UDDAtomIndexHandle;

               // mmdb-extras.h
               mmdb::Residue *res_copy = coot::deep_copy_this_residue_old_style(current_residue_p, alt_conf,
                                                                      whole_res_flag, atom_index_udd_handle,
                                                                      false);

               res_copy->seqNum = new_resno;
               strncpy(res_copy->insCode, new_inscode_str.c_str(), 9);
               std::pair<int, mmdb::Residue *> sn = find_serial_number_for_insert(new_resno,
                                                                                  new_inscode_str,
                                                                                  chain_id);

               // 20180529-PE If you are copying this, don't forget the make_asc at the end. It took
               // me 6 hours to find that that was the problem.

               if (sn.first != -1) { // normal insert

                  this_chain_p->InsResidue(res_copy, sn.first);
                  this_chain_p->TrimResidueTable(); // probably not needed
                  int result = atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_INDEX);
                  if (result != 0) {
                     std::cout << "WARNING:: change_residue_number() PDBCleanup failed " << std::endl;
                  }
                  atom_sel.mol->FinishStructEdit();

                  delete current_residue_p;

                  atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL);
                  atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_INDEX);
                  atom_sel.mol->FinishStructEdit();

               } else {

                  // do an add
                  this_chain_p->AddResidue(res_copy);
                  atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                  atom_sel.mol->FinishStructEdit();
                  delete current_residue_p;
                  atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                  atom_sel.mol->FinishStructEdit();

               }

               done_it = true;
               have_unsaved_changes_flag = true;
               atom_sel = make_asc(atom_sel.mol);
               coot::residue_spec_t old_spec(chain_id, current_resno, current_inscode_str);
               coot::residue_spec_t new_spec(chain_id, new_resno, new_inscode_str);
               update_any_link_containing_residue(old_spec, new_spec);
               make_bonds_type_checked(__FUNCTION__);
            }
         }
      }
   }
   return done_it;
}



// for add OXT
std::pair<bool, int>
molecule_class_info_t::last_residue_in_chain(const std::string &chain_id) const {

   std::pair<short int, int> p(false,0);
   int biggest_resno = -99999;

   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(i_chain);
	    std::string mol_chain(chain_p->GetChainID());
	    if (mol_chain == chain_id) {
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::Residue *residue_p;
	       for (int ires=0; ires<nres; ires++) { // ires is a serial number
		  residue_p = chain_p->GetResidue(ires);
		  if (residue_p->GetSeqNum() > biggest_resno) {
		     biggest_resno = residue_p->GetSeqNum();
		     p.first = true;
		  }
	       }
	    }
	 }
      }
   }
   p.second = biggest_resno;
   return p;
}

std::pair<bool, int>
molecule_class_info_t::last_protein_residue_in_chain(const std::string &chain_id) const {

   std::pair<short int, int> p(false,0);
   int biggest_resno = -99999;

   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(i_chain);
	    std::string mol_chain(chain_p->GetChainID());
	    if (mol_chain == chain_id) {
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::Residue *residue_p;
	       for (int ires=0; ires<nres; ires++) { // ires is a serial number
		  residue_p = chain_p->GetResidue(ires);
		  if (residue_p) {
		     std::string rn = residue_p->GetResName();
		     if (coot::util::is_standard_amino_acid_name(rn)) {
			if (residue_p->GetSeqNum() > biggest_resno) {
			   biggest_resno = residue_p->GetSeqNum();
			   p.first = true;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   p.second = biggest_resno;
   return p;

}



std::pair<bool, int>
molecule_class_info_t::first_residue_in_chain(const std::string &chain_id) const {

   std::pair<bool, int> p(false,0);
   int smallest_resno = 999999;

   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(i_chain);
	    std::string mol_chain(chain_p->GetChainID());
	    if (mol_chain == chain_id) {
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::Residue *residue_p;
	       for (int ires=0; ires<nres; ires++) { // ires is a serial number
		  residue_p = chain_p->GetResidue(ires);
		  if (residue_p->GetSeqNum() < smallest_resno) {
		     smallest_resno = residue_p->GetSeqNum();
		     p.first = true;
		  }
	       }
	    }
	 }
      }
   }
   p.second = smallest_resno;
   return p;
}


// return NULL on no last residue.
mmdb::Residue *
molecule_class_info_t::last_residue_in_chain(mmdb::Chain *chain_p) const {

   mmdb::Residue *res = NULL;
   int biggest_resno = -99999;

   int n_residues = chain_p->GetNumberOfResidues();
   for (int i=0; i<n_residues; i++) {
      mmdb::Residue *r = chain_p->GetResidue(i);
      if (r->GetSeqNum() >= biggest_resno) {
         biggest_resno = r->GetSeqNum();
         res = r;
      }
   }
   return res;
}



// validation

void
molecule_class_info_t::find_deviant_geometry(float strictness) {

   if (atom_sel.n_selected_atoms > 0) {
      std::vector<coot::atom_spec_t> fixed_atom_specs;
      short int have_flanking_residue_at_end = 0;
      short int have_flanking_residue_at_start = 0;
      // int resno_1, resno_2;

      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(i_chain);
	    std::string chain_id(chain_p->GetChainID());

	    std::pair<short int, int> resno_1 = first_residue_in_chain(chain_id);
	    std::pair<short int, int> resno_2 =  last_residue_in_chain(chain_id);

	    if (! (resno_1.first && resno_2.first)) {
	       std::cout << "WARNING: Error getting residue ends in find_deviant_geometry\n";
	    } else {

	       short int have_disulfide_residues = 0;
	       std::string altconf = "";

	       int selHnd = atom_sel.mol->NewSelection();
	       int nSelResidues;
	       mmdb::PResidue *SelResidues = NULL;

	       atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
				    chain_id.c_str(),
				    resno_1.second, "*",
				    resno_2.second, "*",
				    "*",  // residue name
				    "*",  // Residue must contain this atom name?
				    "*",  // Residue must contain this Element?
				    "*",  // altLocs
				    mmdb::SKEY_NEW // selection key
				    );
	       atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

	       // kludge in a value for icheck.
	       std::vector<std::string> kludge;
	       std::pair<int, std::vector<std::string> > icheck(1, kludge);

	       // coot::util::check_dictionary_for_residues(SelResidues, nSelResidues,
	       //                                           graphics_info_t::Geom_p());

	       if (icheck.first == 0) {
		  for (unsigned int icheck_res=0; icheck_res<icheck.second.size(); icheck_res++) {
		     std::cout << "WARNING:: Failed to find restraints for "
			       << icheck.second[icheck_res] << std::endl;
		  }
	       }

	       std::cout << "INFO:: " << nSelResidues << " residues selected for deviant object"
			 << std::endl;

	       if (nSelResidues > 0) {

		  mmdb::Manager *residues_mol =
		     create_mmdbmanager_from_res_selection(SelResidues, nSelResidues,
							   have_flanking_residue_at_start,
							   have_flanking_residue_at_end,
							   altconf,
							   chain_id,
							   0 // 0 because we are not in alt conf split
							   );
		  clipper::Xmap<float> dummy_xmap;

		  // coot::restraints_container_t
		  //    restraints(resno_1.second,
		  //               resno_2.second,
		  //               have_flanking_residue_at_start,
		  //               have_flanking_residue_at_end,
		  //               have_disulfide_residues,
		  //               altconf,
		  //               (char *) mol_chain.c_str(),
		  //               residues_mol,
		  //               fixed_atom_specs,
		  //               &dummy_xmap);

		  coot::restraints_container_t restraints(SelResidues, nSelResidues, chain_id, residues_mol, &dummy_xmap);
	       }
	    }
	 }
      }
   }
}


// ------------------------------------------------------------------
//                       sequence assignment
// ------------------------------------------------------------------

#include "high-res/sequence-assignment.hh"

// This is not the sequence assignment function that you are looking for.
// (it is the ancient high-ree function that does nothing)
//
void
molecule_class_info_t::assign_sequence(const clipper::Xmap<float> &xmap,
                                       const std::string &chain_id) {

   std::cout << "debug:: in assign_sequence() there are " << input_sequence.size() << " sequences input "
             << "for imol " << imol_no << std::endl;

   for (unsigned int i=0; i<input_sequence.size(); i++) {
      if (input_sequence[i].first == chain_id){
         coot::sequence_assignment::side_chain_score_t scs;
         std::string seq = input_sequence[i].second;
         std::cout << "calling scs.add_fasta_sequence() for chain_id \"" << chain_id << "\" seq"
                   << seq << std::endl;
         scs.add_fasta_sequence(chain_id, seq);
      }
   }
}


// Distances
std::vector<clipper::Coord_orth>
molecule_class_info_t::distances_to_point(const clipper::Coord_orth &pt,
                                          double min_dist,
                                          double max_dist) {

   std::vector<clipper::Coord_orth> v;
   if (atom_sel.n_selected_atoms > 0) {
      for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
         mmdb::Atom *at = atom_sel.atom_selection[iat];
         if (! coot::is_hydrogen_atom(at) || draw_hydrogens_flag) {
            clipper::Coord_orth atp(atom_sel.atom_selection[iat]->x,
                                    atom_sel.atom_selection[iat]->y,
                                    atom_sel.atom_selection[iat]->z);
            if (clipper::Coord_orth::length(pt, atp) <= max_dist) {
               if (clipper::Coord_orth::length(pt, atp) >= min_dist) {
                  v.push_back(atp);
               }
            }
         }
      }
   }
   return v;

}

// logical_operator_and_or_flag 0 for AND 1 for OR.
//
//
// If outlier_sigma_level is less than -50 then don't test for sigma level
// If min_dist < 0, don't test for min dist
// if max_dist < 0, don't test for max dist
// if b_factor_lim < 0, don't test for b_factor
//
std::vector <coot::atom_spec_t>
molecule_class_info_t::find_water_baddies(float b_factor_lim, const clipper::Xmap<float> &xmap_in,
                                          float map_sigma,
                                          float outlier_sigma_level,
                                          float min_dist, float max_dist,
                                          short int ignore_part_occ_contact_flag,
                                          short int ignore_zero_occ_flag,
                                          short int logical_operator_and_or_flag) {

   if (logical_operator_and_or_flag == 0)
      return find_water_baddies_AND(b_factor_lim, xmap_in, map_sigma, outlier_sigma_level,
                                    min_dist, max_dist, ignore_part_occ_contact_flag,
                                    ignore_zero_occ_flag);
   else
      return find_water_baddies_OR(b_factor_lim, xmap_in, map_sigma, outlier_sigma_level,
                                   min_dist, max_dist, ignore_part_occ_contact_flag,
                                   ignore_zero_occ_flag);

}


// This is logical opertator AND on the search criteria.
//
std::vector <coot::atom_spec_t>
molecule_class_info_t::find_water_baddies_AND(float b_factor_lim, const clipper::Xmap<float> &xmap_in,
                                              float map_sigma,
                                              float outlier_sigma_level,
                                              float min_dist, float max_dist,
                                              short int ignore_part_occ_contact_flag,
                                              short int ignore_zero_occ_flag) {

   std::vector <coot::atom_spec_t> v;
   std::vector <std::pair<int, double> > idx;
   double den;

   // put in one loop otherwise we get duplicates (or even more)

   for (int i=0; i<atom_sel.n_selected_atoms; i++) {
      if (atom_sel.atom_selection[i]->tempFactor > b_factor_lim) {
         if (! atom_sel.atom_selection[i]->isTer()) {
            std::string resname = atom_sel.atom_selection[i]->GetResName();
            if (resname == "WAT" || resname == "HOH") {
               clipper::Coord_orth a(atom_sel.atom_selection[i]->x,
                                     atom_sel.atom_selection[i]->y,
                                     atom_sel.atom_selection[i]->z);
               den = coot::util::density_at_point(xmap_in, a);

               if (den > outlier_sigma_level*map_sigma || map_sigma < 0 || outlier_sigma_level < 0) {
                  if( min_dist < 0 && max_dist < 0) {
                     idx.push_back(std::make_pair(i, den));
                  } else {
                     // need to check the distances too, only around this water!?
                     // try with mmdb
                     // first select the water
                     int SelHnd_wat;
                     SelHnd_wat = atom_sel.mol->NewSelection();
                     atom_sel.mol->Select(SelHnd_wat, mmdb::STYPE_ATOM, 0,
                                          atom_sel.atom_selection[i]->GetChainID(),
                                          atom_sel.atom_selection[i]->GetSeqNum(), "*",
                                          atom_sel.atom_selection[i]->GetSeqNum(), "*",
                                          "*",  // residue name
                                          "*",  // Residue must contain this atom name?
                                          "*",  // Residue must contain this Element?
                                          "*",  // altLocs
                                          mmdb::SKEY_NEW // selection key
                                          );

                     if (min_dist < 0)
                        min_dist = 0.0;
                     if (max_dist < 0)
                        max_dist = 10.0; // should be enough?!
                     int nSelAtoms_wat;
                     mmdb::PPAtom SelAtom_wat;
                     atom_sel.mol->GetSelIndex(SelHnd_wat, SelAtom_wat, nSelAtoms_wat);
                     atom_sel.mol->SelectNeighbours(SelHnd_wat,
                                                    mmdb::STYPE_ATOM,
                                                    SelAtom_wat,
                                                    nSelAtoms_wat,
                                                    min_dist, max_dist,
                                                    mmdb::SKEY_OR);

                     atom_sel.mol->GetSelIndex(SelHnd_wat, SelAtom_wat, nSelAtoms_wat);
                     // Do we need to remove the hydrogens?
                     if (nSelAtoms_wat == 1) {
                        // selection always contains the atom around which the neighbours are found
                        //
                        // no atoms between min and max distance found, i.e. closest contact must
                        // be outside the range.
                        idx.push_back(std::make_pair(i, den));
                     }
                  }
               }
            }
         }
      }
   }

   // now add the atoms (with all info, mabe could be done above too)
   for (unsigned int i=0; i<idx.size(); i++) {

      std::string s = "B fac: ";
      s += coot::util::float_to_string(atom_sel.atom_selection[idx[i].first]->tempFactor);
      s += "   ED: ";
      s += coot::util::float_to_string(idx[i].second);
      s += " rmsd";

      coot::atom_spec_t atom_spec(atom_sel.atom_selection[idx[i].first], s);
      atom_spec.float_user_data = atom_sel.atom_selection[idx[i].first]->occupancy;
      v.push_back(atom_spec);
   }

   return v;
}


#include "coot-utils/find-water-baddies.hh"

// This is logical opertator AND on the search criteria.
//
// If outlier_sigma_level is less than -50 then don't test for sigma level
// If min_dist < 0, don't test for min dist
// if max_dist < 0, don't test for max dist
// if b_factor_lim < 0, don't test for b_factor
//
// Hydrogens are ignored in neighbour distance search
//
std::vector <coot::atom_spec_t>
molecule_class_info_t::find_water_baddies_OR(float b_factor_lim, const clipper::Xmap<float> &xmap_in,
                                             float map_in_sigma,
                                             float outlier_sigma_level,
                                             float min_dist, float max_dist,
                                             short int ignore_part_occ_contact_flag,
                                             short int ignore_zero_occ_flag) {

   // 20231016-PE function moved into util so that it can be used by api
   //
   std::vector <coot::atom_spec_t> v = coot::find_water_baddies_OR(atom_sel, b_factor_lim, xmap_in, map_in_sigma,
                                                                   outlier_sigma_level, min_dist, max_dist,
                                                                   ignore_part_occ_contact_flag, ignore_zero_occ_flag);
   return v;
}


// shelx stuff
//
std::pair<int, std::string>
molecule_class_info_t::write_shelx_ins_file(const std::string &filename) {

   // std::cout << "DEBUG:: starting write_shelx_ins_file in molecule "<< std::endl;
   // shelxins.debug();

   std::pair<int, std::string> p(1, "");

   if (atom_sel.n_selected_atoms > 0) {
      p = shelxins.write_ins_file(atom_sel.mol, filename, is_from_shelx_ins_flag);
//       std::cout << "DEBUG:: in molecule_class_info_t::write_ins_file "
//                 << "got values " << p.first << " " << p.second
//                 << std::endl;
   } else {
      p.second = "WARNING:: No atoms to write!";
   }
   return p;
}


// Return a variable like reading a pdb file (1 on success, -1 on failure)
//
// This function doesn't get called by the normal handle_read_draw_molecule()
//
int
molecule_class_info_t::read_shelx_ins_file(const std::string &filename) {

   // returned a pair: status (0: bad), udd_afix_handle (-1 bad)

   int istat = 1;
   coot::shelx_read_file_info_t p = shelxins.read_file(filename);
   if (p.status == 0) {
      std::cout << "WARNING:: bad status in read_shelx_ins_file" << std::endl;
      istat = -1;
   } else {

      int udd_afix_handle = p.mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");
      if (false)
         std::cout << "DEBUG:: in  get_atom_selection udd_afix_handle is "
                   << udd_afix_handle << " and srf.udd_afix_handle was "
                   << udd_afix_handle << std::endl;

      if (p.udd_afix_handle == -1) {
         std::cout << "ERROR:: bad udd_afix_handle in read_shelx_ins_file"
                   << std::endl;
         istat = -1;
      }

      if (istat == 1) {
         // initialize some things.
         //
         atom_sel = make_asc(p.mol);
         // FIXME? 20070721, shelx presentation day
         //          fix_hydrogen_names(atom_sel); // including change " H0 " to " H  "
         short int is_undo_or_redo = 0;
         graphics_info_t g;

         mmdb::mat44 my_matt;
         int err = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
         if (err != mmdb::SYMOP_Ok) {
            std::cout << "!! Warning:: No symmetry available for this molecule"
                      << std::endl;
         } else {
            std::cout << "Symmetry available for this molecule" << std::endl;
         }
         is_from_shelx_ins_flag = 1;

         // std::cout << " ##### initing coord things in read_shelx_ins_file" << std::endl;
         initialize_coordinate_things_on_read_molecule_internal(filename,
                                                                is_undo_or_redo);
         // std::cout << " ##### done initing coord things in read_shelx_ins_file" << std::endl;

         set_have_unit_cell_flag_maybe(true); // but will always have symmetry

         if (molecule_is_all_c_alphas()) {
            bool force_rebonding =  true;
            ca_representation(force_rebonding);
         } else {

            short int do_rtops_flag = 0;

            // Hmmm... why is this commented?  Possibly from ncs ghost
            // toubles from many months ago?

            // 0.7 is not used (I think) if do_rtops_flag is 0
            // int nghosts = fill_ghost_info(do_rtops_flag, 0.7);
            // std::cout << "INFO:: found " << nghosts << " ghosts\n";

            // I'll reinstate it.
            int nmodels = atom_sel.mol->GetNumberOfModels();
            if (nmodels == 1) {
               // int nghosts =
               fill_ghost_info(do_rtops_flag, 0.7);
            }

            // Turn off hydrogen display if this is a protein
            // (currently the hydrogen names are different from a
            // shelx molecule, leading to a mess when refining.  So
            // let's just undisplay those hydrogens :) 20070721
            //
            // 20080126 Let reactivate the hydrogens.
            //
            //             if (p.is_protein_flag)
            //                set_draw_hydrogens_state(0);

            if (! is_undo_or_redo)
               bond_width = g.default_bond_width; // bleugh, perhaps this should
                                                  // be a passed parameter?

            // Generate bonds and save them in the graphical_bonds_container
            // which has static data members.
            //
            if (bonds_box_type == coot::UNSET_TYPE)
               bonds_box_type = coot::NORMAL_BONDS;
            make_bonds_type_checked(__FUNCTION__);
         }

         // debug();

         short int reset_rotation_centre = 1;
         //
         if (g.recentre_on_read_pdb || g.n_molecules() == 0) {
            // n_molecules is updated in c-interface.cc
            // std::cout << " ##### setting rotation centre in read_shelx_ins_file" << std::endl;
            coot::Cartesian c = ::centre_of_molecule(atom_sel);
            // std::cout << "debug:: n atoms " << atom_sel.n_selected_atoms << std::endl;
            // std::cout << "debug:: rotation centre " << c << std::endl;
            if (reset_rotation_centre)
               g.setRotationCentre(c);
            // std::cout << " ##### done setting rotation centre in read_shelx_ins_file" << std::endl;
         }

         draw_it = 1;
         if (graphics_info_t::show_symmetry == 1) {
            update_symmetry();
         }
      }

      // save state strings
      save_state_command_strings_.push_back("read-shelx-ins-file");
      save_state_command_strings_.push_back(single_quote(filename));
      save_state_command_strings_.push_back("1"); // recentre flag
   }
   // std::cout << " ##### read_shelx_ins_file returning " << istat << std::endl;
   return istat;
}

int
molecule_class_info_t::add_shelx_string_to_molecule(const std::string &str) {

   int istat = 0;
   if (is_from_shelx_ins_flag) {
      shelxins.add_pre_atoms_line(str);
      istat = 1;
   }
   return istat;
}


// -----------------------------------------------------------------
//              display list lovelies
// -----------------------------------------------------------------

bool
molecule_class_info_t::has_display_list_objects() {

   bool r = 0;
   if (draw_it) {
      if (display_list_tags.size() > 0) {
         r = 1;
      }
   }
   return r;
}


int
molecule_class_info_t::draw_display_list_objects(int GL_context) {

   int n_objects = 0;
   if (draw_it) {
      if (display_list_tags.size() > 0) {
         std::vector<coot::display_list_object_info>::const_iterator it;
         glEnable(GL_COLOR_MATERIAL);
         for (it=display_list_tags.begin(); it!=display_list_tags.end(); it++) {
            if (! it->is_closed) {
               if (it->display_it) {
                  n_objects++;
                  if (GL_context == GL_CONTEXT_MAIN) {
                     glCallList(it->tag_1);
                  }
                  if (GL_context == GL_CONTEXT_SECONDARY) {
                     glCallList(it->tag_2);
                  }
               }
            }
         }
         glDisable(GL_COLOR_MATERIAL);
      }
   }
   return n_objects;
}

// return the display list tag
int
molecule_class_info_t::make_ball_and_stick(const std::string &atom_selection_str,
                                           float bond_thickness, float sphere_size,
                                           bool do_spheres_flag, gl_context_info_t gl_info,
                                           const coot::protein_geometry *geom) {

   std::cout << "molecule make_ball_and_stick(A) called ..." << std::endl;

   coot::display_list_object_info dloi;
   // modify a copy of dloi and return it
   graphics_info_t::make_gl_context_current(graphics_info_t::GL_CONTEXT_MAIN);
   coot::display_list_object_info dloi_1 = make_ball_and_stick(atom_selection_str, bond_thickness,
                                                               sphere_size,
                                                               do_spheres_flag, 0, dloi, geom);

   if (gl_info.widget_2) {
      graphics_info_t::make_gl_context_current(graphics_info_t::GL_CONTEXT_SECONDARY);
      // modify a copy of dloi_1 and return it
      coot::display_list_object_info dloi_2 = make_ball_and_stick(atom_selection_str, bond_thickness,
                                                                  sphere_size,
                                                                  do_spheres_flag, 1, dloi_1, geom);

      if (0)
         std::cout << "pushing back dloi (2 contexts) to display_list_tags: with tag_1 = "
                   << dloi_2.tag_1 << " and tag_2 = " << dloi_2.tag_2 << std::endl;
      display_list_tags.push_back(dloi_2);
      graphics_info_t::make_gl_context_current(graphics_info_t::GL_CONTEXT_MAIN);
   } else {
      if (0)
         std::cout << "pushing back dloi (1 context)  to display_list_tags: with tag_1 = "
                   << dloi_1.tag_1 << " and tag_2 = " << dloi_1.tag_2 << std::endl;
      display_list_tags.push_back(dloi_1);
   }

   // std::cout << ".... make_ball_and_stick(A) add_reps.size() " << add_reps.size() << std::endl;
   return (display_list_tags.size() - 1);
}


// rename/delete the function

// return the display list info
coot::display_list_object_info
molecule_class_info_t::make_ball_and_stick(const std::string &atom_selection_str,
					   float bond_thickness, float sphere_size,
					   bool do_spheres_flag, bool is_second_context,
					   coot::display_list_object_info dloi,
					   const coot::protein_geometry *geom) {
#if 0

   // hack
   bool against_a_dark_background = true;

   // std::cout << "molecule make_ball_and_stick(B) called ..." << std::endl;

   // Use draw hydrogens flag that has been set already for this molecule.

   int i = -1;
   if (has_model()) {
      int SelHnd = atom_sel.mol->NewSelection();
      atom_sel.mol->Select(SelHnd, mmdb::STYPE_ATOM,
                           atom_selection_str.c_str(),
                           mmdb::SKEY_OR);
      int n_selected_atoms;
      mmdb::PPAtom atom_selection = NULL;
      atom_sel.mol->GetSelIndex(SelHnd, atom_selection, n_selected_atoms);
      atom_selection_container_t asc = atom_sel;
      asc.atom_selection = atom_selection;
      asc.n_selected_atoms = n_selected_atoms;
      asc.SelectionHandle = SelHnd;

      Bond_lines_container bonds(asc, imol_no, geom);
      graphical_bonds_container bonds_box_local = bonds.make_graphical_bonds();

      // start display list object
      GLuint bonds_tag = glGenLists(1);
      glNewList(bonds_tag, GL_COMPILE);
      if (is_second_context) {
         dloi.tag_2 = bonds_tag;
         // std::cout << "debug:: adding second context tag to dloi " << bonds_tag << std::endl;
      } else {
         dloi.tag_1 = bonds_tag;
         // std::cout << "debug:: adding first  context tag to dloi " << bonds_tag << std::endl;
      }

      GLfloat bgcolor[4] = {0.8, 0.8, 0.8, 1.0};
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glMaterialfv(GL_FRONT, GL_SPECULAR, bgcolor);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 40);

      for (int ii=0; ii<bonds_box_local.num_colours; ii++) {
         graphical_bonds_lines_list<graphics_line_t> &ll = bonds_box_local.bonds_[ii];
         set_bond_colour_by_mol_no(ii, against_a_dark_background);

         if (bond_colour_internal.size() > 2)
            for (int jj=0; jj<3; jj++)
               bgcolor[jj] = bond_colour_internal[jj];

         for (int j=0; j< bonds_box_local.bonds_[ii].num_lines; j++) {
            glPushMatrix();
            glTranslatef(ll.pair_list[j].positions.getFinish().get_x(),
                         ll.pair_list[j].positions.getFinish().get_y(),
                         ll.pair_list[j].positions.getFinish().get_z());
            double base = bond_thickness;
            double top = bond_thickness;
            if (ll.thin_lines_flag) {
               base *= 0.55;
               top  *= 0.55;
            }
            if (ll.pair_list[j].cylinder_class == graphics_line_t::DOUBLE ||
                ll.pair_list[j].cylinder_class == graphics_line_t::TRIPLE) {
               base *= 0.65; // needs optimizing
               top  *= 0.65;
            }
            coot::Cartesian bond_height =
               ll.pair_list[j].positions.getFinish() - ll.pair_list[j].positions.getStart();
            double height = bond_height.amplitude();
            int slices = 10;
            int stacks = 2;

            //             This code from ccp4mg's cprimitive.cc (but modified)
            //              -----
            double ax;
            double rx = 0;
            double ry = 0;
            double length = bond_height.length();
            double vz = bond_height.get_z();

            bool rot_x = false;
            if(fabs(vz)>1e-7){
               ax = 180.0/M_PI*acos(vz/length);
               if(vz<0.0) ax = -ax;
               rx = -bond_height.get_y()*vz;
               ry = bond_height.get_x()*vz;
            }else{
               double vx = bond_height.get_x();
               double vy = bond_height.get_y();
               ax = 180.0/M_PI*acos(vx/length);
               if(vy<0) ax = -ax;
               rot_x = true;
            }

            if (rot_x) {
               glRotated(90.0, 0.0, 1.0, 0.0);
               glRotated(ax,  -1.0, 0.0, 0.0);
            } else {
               glRotated(ax, rx, ry, 0.0);
            }
            //             --------

            GLUquadric* quad = gluNewQuadric();
            glScalef(1.0, 1.0, -1.0); // account for mg maths :-)
            gluCylinder(quad, base, top, height, slices, stacks);
            // gluQuadricNormals(quad, GL_SMOOTH);

            if (0) // debugging end-caps
               std::cout << j << " begin_end_cap: " << ll.pair_list[j].has_begin_cap
                         << " end_cap: " << ll.pair_list[j].has_end_cap
                         << std::endl;

            if (ll.pair_list[j].has_end_cap) {
               glPushMatrix();
               glScalef(1.0, 1.0, -1.0);
               gluDisk(quad, 0, base, slices, 2);
               glPopMatrix();
            }
            if (ll.pair_list[j].has_begin_cap) {
               glPushMatrix();
               glTranslated(0,0,height);
               gluDisk(quad, 0, base, slices, 2);
               glPopMatrix();
            }

            gluDeleteQuadric(quad);
            glPopMatrix();
         }
      }

      if (do_spheres_flag) {
         int slices = 20;
         int stacks = 20;
         for (int i=0; i<bonds_box_local.n_atom_centres_; i++) {
            float local_sphere_size = sphere_size;
            if (bonds_box_local.atom_centres_[i].is_hydrogen_atom)
               local_sphere_size = 0.11; // small (and cute)
            if (bonds_box_local.atom_centres_colour_[i] == HYDROGEN_GREY_BOND)
               local_sphere_size *= 0.55; // matches thin_lines_flag test above
            set_bond_colour_by_mol_no(bonds_box_local.atom_centres_colour_[i],
                                      against_a_dark_background);
            glPushMatrix();

            // GLfloat bgcolor[4]={bond_colour_internal[0],
            // bond_colour_internal[1],
            // bond_colour_internal[2],
            // 1.0};

            glMaterialfv(GL_FRONT, GL_SPECULAR, bgcolor);
            glTranslatef(bonds_box_local.atom_centres_[i].position.get_x(),
                         bonds_box_local.atom_centres_[i].position.get_y(),
                         bonds_box_local.atom_centres_[i].position.get_z());

            GLUquadric* quad = gluNewQuadric();

            gluSphere(quad, local_sphere_size, slices, stacks);
            // gluQuadricNormals(quad, GL_SMOOTH);
            gluDeleteQuadric(quad);
            glPopMatrix();
         }

      }

      if (bonds_box_local.have_rings()) {
         set_bond_colour_by_mol_no(0, against_a_dark_background);
         for (unsigned int ir=0; ir<bonds_box_local.rings.size(); ir++) {
            glPushMatrix();

             GL_matrix m(bonds_box_local.rings[ir].normal);

            glTranslatef(bonds_box_local.rings[ir].centre.x(),
                         bonds_box_local.rings[ir].centre.y(),
                         bonds_box_local.rings[ir].centre.z());
            glMultMatrixf(m.get());
	    glScalef(1.0, 1.0, 0.7);
#if 0
 	    glutSolidTorus(bonds_box_local.rings[ir].inner_radius,
			   bonds_box_local.rings[ir].outer_radius,
 			   bonds_box_local.rings[ir].n_sides,
 			   bonds_box_local.rings[ir].n_rings);
#endif
      std::cout << "Fix Solid Torus\n"	    ;
	    glPopMatrix();
	 }
      }

      glEndList();
      bonds_box_local.clear_up();
      atom_sel.mol->DeleteSelection(SelHnd);
   }
#endif
   return dloi;
}

void
molecule_class_info_t::clear_display_list_object(GLuint tag) {

   // actually, clear them all, not just those (or that one) with tag:

//    std::list<coot::display_list_object_info>::const_iterator it;
//    for (it=display_list_tags.begin(); it!=display_list_tags.end(); it++) {
//    }

   display_list_tags.clear();
}


void
molecule_class_info_t::set_occupancy_residue_range(const std::string &chain_id, int ires1, int ires2, float occ_val) {

   if (ires2 < ires1) {
      int tmp = ires1;
      ires1 = ires2;
      ires2 = tmp;
   }

   mmdb::PPAtom SelAtoms = NULL;
   int nSelAtoms;
   int SelHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(SelHnd, 0, (char *) chain_id.c_str(),
                             ires1, "*",
                             ires2, "*",
                             "*", // residue name
                             "*",
                             "*", // elements
                             "*");

   atom_sel.mol->GetSelIndex(SelHnd, SelAtoms, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "Sorry. Could not find residue range " << ires1
                << " to " << ires2 << " in this molecule: ("
                <<  imol_no << ") " << name_ << std::endl;
   } else {
      for (int i=0; i<nSelAtoms; i++) {
         mmdb::Atom *at = SelAtoms[i];
         at->occupancy = occ_val;
         if (at->WhatIsSet & mmdb::ASET_Occupancy) {
            // mmdb::ASET_Occupancy not set in mmdb yet
            at->WhatIsSet |= mmdb::ASET_Occupancy;
         }
      }
      atom_sel.mol->DeleteSelection(SelHnd);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
}

void
molecule_class_info_t::set_b_factor_residue_range(const std::string &chain_id,
                                                  int ires1, int ires2, float b_val) {

   if (ires2 < ires1) {
      int tmp = ires1;
      ires1 = ires2;
      ires2 = tmp;
   }

   mmdb::PPAtom SelAtoms = NULL;
   int nSelAtoms;
   int SelHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(SelHnd, 0, chain_id.c_str(),
                             ires1, "*",
                             ires2, "*",
                             "*", // residue name
                             "*",
                             "*", // elements
                             "*");

   atom_sel.mol->GetSelIndex(SelHnd, SelAtoms, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "Sorry. Could not find residue range " << ires1
                << " to " << ires2 << " in this molecule: ("
                <<  imol_no << ") " << name_ << std::endl;
   } else {
      for (int i=0; i<nSelAtoms; i++) {
         SelAtoms[i]->tempFactor = b_val;
         if (SelAtoms[i]->WhatIsSet & mmdb::ASET_tempFactor) {
            // mmdb::ASET_tempFactor not set in mmdb yet
            SelAtoms[i]->WhatIsSet |= mmdb::ASET_tempFactor;
         }
      }
      atom_sel.mol->DeleteSelection(SelHnd);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
}

// as for occupancies but for B-factors and use a more general
// atom_selection_container. Let's do everything else later and/or
// on scripting level.
// need to copy from moving atoms to 'real' atoms!!! or the other way round?!
void
molecule_class_info_t::set_b_factor_atom_selection(const atom_selection_container_t &asc,
                                                   float b_val, bool moving_atoms) {

  int n_atom = 0;
  int tmp_index;

  make_backup();
  for (int i=0; i<asc.n_selected_atoms; i++) {
    int idx = -1;
    mmdb::Atom *atom = asc.atom_selection[i];
    if (moving_atoms) {
      if (asc.UDDOldAtomIndexHandle >= 0) { // OK for fast atom indexing
        if (atom->GetUDData(asc.UDDOldAtomIndexHandle, tmp_index) == mmdb::UDDATA_Ok) {
          if (tmp_index >= 0) {
            if (moving_atom_matches(atom, tmp_index)) {
              idx = tmp_index;
            } else {
              idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                                 atom->residue->seqNum,
                                                 std::string(atom->GetInsCode()),
                                                 std::string(atom->name),
                                                 std::string(atom->altLoc));
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
                    <<  "), bad GetUDData for this atom " << std::endl;
        }
      } else {

        idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                           atom->residue->seqNum,
                                           std::string(atom->GetInsCode()),
                                           std::string(atom->name),
                                           std::string(atom->altLoc));
        if (idx == -1) {
          std::cout << "DEBUG:: set-bfactor idx: " << idx << "\n";
          std::cout << "ERROR:: failed to find spec for chain-id :"
                    << std::string(atom->residue->GetChainID()) <<  ": "
                    << atom->residue->seqNum << " inscode :"
                    << std::string(atom->GetInsCode()) << ": name :"
                    << std::string(atom->name) << ": altloc :"
                    << std::string(atom->altLoc) << ":" << std::endl;
        }
      }
    } else {
      idx = i;
    }
    if (idx >= 0) {
      n_atom++;
      mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
      mol_atom->SetCoordinates(atom->x,
                               atom->y,
                               atom->z,
                               atom->occupancy,
                               b_val);
    }
  }
  have_unsaved_changes_flag = 1;
  make_bonds_type_checked(__FUNCTION__);
}

// all atoms of specified
void
molecule_class_info_t::set_b_factor_residues(const std::vector<std::pair<coot::residue_spec_t, double> > &rbs) {

   make_backup();
   for (unsigned int i=0; i<rbs.size(); i++) {
      const coot::residue_spec_t &spec = rbs[i].first;
      double b = rbs[i].second;
      mmdb::Residue *residue_p = get_residue(spec);
      if (residue_p) {
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int j=0; j<n_residue_atoms; j++) {
            residue_atoms[j]->tempFactor = b;
         }
      } else {
         std::cout << "WARNING:: No residue for spec " << spec << std::endl;
      }
   }
   have_unsaved_changes_flag = 1;
   atom_sel.mol->FinishStructEdit();
   make_bonds_type_checked(__FUNCTION__);
}

void
molecule_class_info_t::set_b_factor_residue(coot::residue_spec_t spec, float bf) {

   make_backup();
    mmdb::Residue *residue_p = get_residue(spec);
   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int j=0; j<n_residue_atoms; j++) {
         residue_atoms[j]->tempFactor = bf;
      }
   }
   have_unsaved_changes_flag = 1;
   atom_sel.mol->FinishStructEdit();
   make_bonds_type_checked(__FUNCTION__);
}

void
molecule_class_info_t::change_b_factors_of_residue_by(coot::residue_spec_t spec, float bf) {

   make_backup();
    mmdb::Residue *residue_p = get_residue(spec);
   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int j=0; j<n_residue_atoms; j++) {
         residue_atoms[j]->tempFactor += bf;
         if (residue_atoms[j]->tempFactor < 2.0)
            residue_atoms[j]->tempFactor = 2.0;
      }
   }
   have_unsaved_changes_flag = 1;
   atom_sel.mol->FinishStructEdit();
   make_bonds_type_checked(__FUNCTION__);
}





// Change chain id
// return -1 on a conflict
// 1 on good.
// 0 on did nothing
//
std::pair<int, std::string>
molecule_class_info_t::change_chain_id(const std::string &from_chain_id,
                                       const std::string &to_chain_id,
                                       bool use_resno_range,
                                       int start_resno, int end_resno) {

   int istat = 0;
   std::string message("Nothing to say");

   //    std::cout << "DEBUG:: use_resno_range: " << use_resno_range << std::endl;

   if (atom_sel.n_selected_atoms > 0) {

      if (use_resno_range) {

         std::pair<int, std::string> r =
            change_chain_id_with_residue_range(from_chain_id, to_chain_id, start_resno, end_resno);
         istat = r.first;
         message = r.second;

      } else {
      // The usual case, I imagine

         bool target_chain_id_exists = false;

         int n_models = atom_sel.mol->GetNumberOfModels();
         for (int imod=1; imod<=n_models; imod++) {

            mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	    if (! model_p) continue;
            // run over chains of the existing mol
            int nchains = model_p->GetNumberOfChains();
            if (nchains <= 0) {
               std::cout << "bad nchains in molecule " << nchains
                         << std::endl;
            } else {
               for (int ichain=0; ichain<nchains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  if (chain_p == NULL) {
                     // This should not be necessary. It seem to be a
                     // result of mmdb corruption elsewhere - possibly
                     // DeleteChain in update_molecule_to().
                     std::cout << "NULL chain in change chain id" << std::endl;
                  } else {
                     std::string chain_id = chain_p->GetChainID();
                     if (to_chain_id == chain_id) {
                        target_chain_id_exists = true;
                        break;
                     }
                  }
               }
            }
         }

         if (!target_chain_id_exists) {

            n_models = atom_sel.mol->GetNumberOfModels();
            for (int imod=1; imod<=n_models; imod++) {

               mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	       if (! model_p) continue;
               // run over chains of the existing mol
               int nchains = model_p->GetNumberOfChains();
               if (nchains <= 0) {
                  std::cout << "bad nchains in molecule " << nchains
                            << std::endl;
               } else {
                  for (int ichain=0; ichain<nchains; ichain++) {
		     mmdb::Chain *chain_p = model_p->GetChain(ichain);
                     if (chain_p) {
                        std::string chain_id = chain_p->GetChainID();
                        if (from_chain_id == chain_id) {
                           make_backup();
                           chain_p->SetChainID(to_chain_id.c_str());
                           // change the links here

                           int n_changed = coot::util::change_chain_in_links(model_p, from_chain_id, to_chain_id);
                           istat = 1;
                           have_unsaved_changes_flag = 1;
                           atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                           atom_sel.mol->FinishStructEdit();
                           atom_sel = make_asc(atom_sel.mol);
                           make_bonds_type_checked(__FUNCTION__);
                        }
                     }
                  }
               }
            }

         } else {

            // OK, can we do a merge? (do we have non-overlapping residue ranges?)
            // If so, so a change_chain_id_with_residue_range()
            //
            mmdb::Chain *chain_p_from = get_chain(from_chain_id);
            mmdb::Chain *chain_p_to   = get_chain(to_chain_id);
            bool done_merge = false;
            if (chain_p_from) {
               if (chain_p_to) {
                  std::pair<bool, int> min_r_1 = coot::util::min_resno_in_chain(chain_p_from);
                  std::pair<bool, int> max_r_1 = coot::util::max_resno_in_chain(chain_p_from);

                  if (false) {
                     std::cout << "--------- here with min_r_1  " << min_r_1.first << " " << min_r_1.second << std::endl;
                     std::cout << "--------- here with max_r_1  " << max_r_1.first << " " << max_r_1.second << std::endl;
                     std::cout << "--------- here with from_chain_id " << from_chain_id << std::endl;
                     std::cout << "--------- here with to_chain_id " << to_chain_id << std::endl;
                  }

                  if (min_r_1.first) {
                     if (max_r_1.first) {
                        start_resno = min_r_1.second;
                        end_resno = max_r_1.second;
                        std::pair<int, std::string> r =
                           change_chain_id_with_residue_range(from_chain_id, to_chain_id, start_resno, end_resno);
                        istat = r.first;
                        message = r.second;
                     }
                  }
               }
            }

            if (! done_merge) {
               std::cout << "WARNING:: CONFLICT: target chain id " << to_chain_id << " already exists "
                         << "in this molecule" << std::endl;
               message = "WARNING:: CONFLICT: target chain id (";
               message += to_chain_id;
               message += ") already \nexists in this molecule!";
            }
         }
      } // residue range
   } // no atoms

   return std::pair<int, std::string> (istat, message);
}


std::pair<int, std::string>
molecule_class_info_t::change_chain_id_with_residue_range(const std::string &from_chain_id,
                                                          const std::string &to_chain_id,
                                                          int start_resno,
                                                          int end_resno) {

   if (false) {
      std::cout << "-------------------- change_chain_id_with_residue_range ---- " << std::endl;
      std::cout << "-------------------- change_chain_id_with_residue_range ---- " << start_resno << std::endl;
      std::cout << "-------------------- change_chain_id_with_residue_range ---- " <<   end_resno << std::endl;
   }

   std::string message;
   int istat = 0;

   short int target_chain_id_exists = 0;
   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) {

      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (! model_p) continue;

      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      if (nchains <= 0) {
         std::cout << "bad nchains in molecule " << nchains
                   << std::endl;
      } else {
         for (int ichain=0; ichain<nchains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (chain_p == NULL) {
               // This should not be necessary. It seem to be a
               // result of mmdb corruption elsewhere - possibly
               // DeleteChain in update_molecule_to().
               std::cout << "NULL chain in change chain id" << std::endl;
            } else {
               std::string chain_id = chain_p->GetChainID();
               if (to_chain_id == chain_id) {
                  target_chain_id_exists = 1;
                  break;
               }
            }
         }
      }
   }

   if (target_chain_id_exists == 0) {

      // So we are moving residues 12->24 of Chain A to (new) Chain
      // C.  Not very frequent, I suspect.

      // make sure start and end are a sensible way round
      if (end_resno < start_resno) {
         int tmp = end_resno;
         end_resno = start_resno;
         start_resno = tmp;
      }

      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 if (! model_p) continue;
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         if (nchains <= 0) {
            std::cout << "bad nchains in molecule " << nchains
                      << std::endl;
         } else {
            for (int ichain=0; ichain<nchains; ichain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(ichain);
               if (chain_p) {
                  std::string chain_id = chain_p->GetChainID();
                  if (from_chain_id == chain_id) {

                     // So we have the chain from which we wish to move residues
                     make_backup();
                     atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                     // Create a new chain and add it to the molecule
                     mmdb::Chain *new_chain_p = new mmdb::Chain;
                     new_chain_p->SetChainID(to_chain_id.c_str());
                     model_p->AddChain(new_chain_p);

                     int nresidues = chain_p->GetNumberOfResidues();
                     for (int ires=0; ires<nresidues; ires++) {
                        mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                        if (residue_p->GetSeqNum() >= start_resno &&
                            residue_p->GetSeqNum() <= end_resno) {
                           int iseqnum  = residue_p->GetSeqNum();
                           mmdb::pstr inscode = residue_p->GetInsCode();
                           mmdb::Residue *residue_copy = coot::util::deep_copy_this_residue(residue_p);
                           chain_p->DeleteResidue(iseqnum, inscode);
                           new_chain_p->AddResidue(residue_copy);
                        }
                     }

                     istat = 1;
                     have_unsaved_changes_flag = 1;
                     atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                     atom_sel.mol->FinishStructEdit();
                     atom_sel = make_asc(atom_sel.mol);
                     make_bonds_type_checked(__FUNCTION__);
                  }
               }
            }
         }
      }
   } else {

      // target chain alread exists.   Here is where we merge...

      // We need to check that we are not reproducing residues that
      // already exist in the chain.  If we are doing that, we stop
      // and give an error message back telling user that that
      // residues exists already.
      //
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 if (! model_p) continue;

         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         short int residue_already_exists_flag = 0;
         if (nchains <= 0) {
            std::cout << "bad nchains in molecule " << nchains
                      << std::endl;
         } else {
            for (int ichain=0; ichain<nchains; ichain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(ichain);
               if (chain_p) {
                  std::string chain_id = chain_p->GetChainID();
                  if (to_chain_id == chain_id) {

                     int nresidues = chain_p->GetNumberOfResidues();
                     int existing_residue_number = 0; // set later
                     for (int ires=0; ires<nresidues; ires++) {
                        mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                        if (residue_p->GetSeqNum() >= start_resno &&
                            residue_p->GetSeqNum() <= end_resno) {
                           residue_already_exists_flag = 1;
                           existing_residue_number = residue_p->GetSeqNum();
                           break;
                        }

                     }
                     if (residue_already_exists_flag) {
                        message += "CONFLICT!  Residue ";
                        message += coot::util::int_to_string(existing_residue_number);
                        message += " already exists in chain ";
                        message += chain_id;
                        message += "\nChange chain ID of residue range failed.\n";
                     } else {

                        // We are OK to move the residue into the existing chain
                        // (move is done by copy and delete)

                        mmdb::Chain *to_chain = NULL;
                        for (int ichain=0; ichain<nchains; ichain++) {
                           chain_p = model_p->GetChain(ichain);
                           if (chain_p) {
                              std::string chain_id = chain_p->GetChainID();
                              if (to_chain_id == chain_id) {
                                 to_chain = chain_p;
                              }
                           }
                        }


                        if (to_chain) {
                           make_backup();
                           for (int ichain=0; ichain<nchains; ichain++) {
                              chain_p = model_p->GetChain(ichain);
                              int to_chain_nresidues = chain_p->GetNumberOfResidues();
                              if (chain_p) {
                                 std::string chain_id = chain_p->GetChainID();
                                 if (from_chain_id == chain_id) {
                                    for (int ires=0; ires<to_chain_nresidues; ires++) {
                                       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                       if (residue_p->GetSeqNum() >= start_resno &&
                                           residue_p->GetSeqNum() <= end_resno) {

                                          int iseqnum  = residue_p->GetSeqNum();
                                          mmdb::pstr inscode = residue_p->GetInsCode();
                                          mmdb::Residue *residue_copy =
                                             coot::util::deep_copy_this_residue_add_chain(residue_p, "", 1, 1);
                                          // delete the residue in the "from" chain:
                                          chain_p->DeleteResidue(iseqnum, inscode);

                                          //
                                          change_chain_id_with_residue_range_helper_insert_or_add(to_chain, residue_copy);
                                          // this is done by the deep_copy
                                          // to_chain->AddResidue(residue_copy);
                                       }
                                    }
                                 }
                              }
                           }
                           istat = 1;
                           have_unsaved_changes_flag = 1;
                           atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                           atom_sel.mol->FinishStructEdit();
                           atom_sel = make_asc(atom_sel.mol);
                           make_bonds_type_checked(__FUNCTION__);
                        }
                     }
                  }
               }
               if (residue_already_exists_flag)
                  break;
            }
         }
         if (residue_already_exists_flag)
            break;
      }
   }

   return std::pair<int, std::string> (istat, message);

}


void
molecule_class_info_t::change_chain_id_with_residue_range_helper_insert_or_add(mmdb::Chain *to_chain_p, mmdb::Residue *new_residue) {

   // OK, if we can, let's try to *insert* the new_residue into the
   // right place in the to_chain_p.  If we don't manage to insert it,
   // let's fall back and simply add it.  The new residue is inserted
   // *before* the residue specified by the given serial number.

   // Let's use the serial number interface to InsResidue()

   int resno_new_residue = new_residue->GetSeqNum();
   std::string ins_code = new_residue->GetInsCode();
   int target_res_serial_number = coot::RESIDUE_NUMBER_UNSET;
   int target_res_seq_num = resno_new_residue; // simply ignore ins codes :)
   std::string target_res_ins_code = ""; // ignore this.  Fix later.
   int best_seq_num_diff = 99999999;

   int n_chain_residues;
   mmdb::PResidue *chain_residues;
   to_chain_p->GetResidueTable(chain_residues, n_chain_residues);
   for (int iserial=0; iserial<n_chain_residues; iserial++) {
      int chain_residue_seq_num = chain_residues[iserial]->GetSeqNum();
      int this_seq_num_diff = chain_residue_seq_num - resno_new_residue;
      if (this_seq_num_diff > 0) {
         if (this_seq_num_diff < best_seq_num_diff) {
            best_seq_num_diff = this_seq_num_diff;
            target_res_serial_number = iserial;
         }
      }
   }

   if (target_res_serial_number != coot::RESIDUE_NUMBER_UNSET) {
      // Good stuff
      // std::cout << "Debugging, inserting residue here...." << std::endl;
      to_chain_p->InsResidue(new_residue, target_res_serial_number);
   } else {
      // std::cout << "Debugging, adding residue here...." << std::endl;
      to_chain_p->AddResidue(new_residue);
   }
}




// nomenclature errors
// return a vector of the changed residues (used for updating the rotamer graph)
std::vector<mmdb::Residue *>
molecule_class_info_t::fix_nomenclature_errors(coot::protein_geometry *geom_p) {
                                                      // by looking for bad rotamers in
                                                      // some residue types and alter ing
                                                          // the atom names to see if they get
                                                      // more likely rotamers

   std::vector<mmdb::Residue *> r;
   if (atom_sel.n_selected_atoms > 0) {
      make_backup();
      coot::nomenclature n(atom_sel.mol);
      r = n.fix(geom_p);
      have_unsaved_changes_flag = 1;
   }
   return r;
}


// nomenclature errors
// return a vector of the changed residues (used for updating the rotamer graph)
std::vector<std::pair<std::string, coot::residue_spec_t> >
molecule_class_info_t::list_nomenclature_errors(const coot::protein_geometry *geom_p) {
                                                      // by looking for bad rotamers in
                                                      // some residue types and alter ing
                                                          // the atom names to see if they get
                                                      // more likely rotamers

   std::vector<mmdb::Residue *> r;
   if (atom_sel.n_selected_atoms > 0) {
      coot::nomenclature n(atom_sel.mol);
      r = n.list(geom_p);
   }
   std::vector<std::pair<std::string, coot::residue_spec_t> > rs;
   if (r.size()) {
      rs.resize(r.size());
      for (unsigned int i=0; i<r.size(); i++) {
         std::pair<std::string, coot::residue_spec_t> p(r[i]->GetResName(), coot::residue_spec_t(r[i]));
         rs[i] = p;
      }
   }
   return rs;
}

int
molecule_class_info_t::cis_trans_conversion(const std::string &chain_id, int resno, const std::string &inscode,
                                            mmdb::Manager *standard_residues_mol) {

   int imod = 1;

   bool found = 0;
   int r = 0; // returned value
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      if (chain_id == chain_p->GetChainID()) {
         int nres = chain_p->GetNumberOfResidues();
         mmdb::PResidue residue_p;
         mmdb::Atom *at;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            if (residue_p->GetSeqNum() == resno) {
               if (inscode == residue_p->GetInsCode()) {
                  int n_atoms = residue_p->GetNumberOfAtoms();

                  for (int iat=0; iat<n_atoms; iat++) {
                     at = residue_p->GetAtom(iat);

                     if (std::string(at->name) != " N  ") {

                        found = 1;
                        r = cis_trans_conversion(at, false, standard_residues_mol);
                        if (r) {
                           make_bonds_type_checked(__FUNCTION__);
                           have_unsaved_changes_flag = true;
                        }
                     }
                     if (found)
                        break;
                  }
                  if (found)
                     break;
               }
            }
            if (found)
               break;
         }
      }
      if (found)
         break;
   }
   return r;

}


// ---- cis <-> trans conversion
int
molecule_class_info_t::cis_trans_conversion(mmdb::Atom *at, short int is_N_flag, mmdb::Manager *standard_residues_mol) {

   // called from graphics_info_t::cis_trans_conversion()

   make_backup();
   mmdb::Manager *mol = atom_sel.mol;
   int status = coot::util::cis_trans_conversion(at, is_N_flag, mol, standard_residues_mol);
   if (status) {
      make_bonds_type_checked();
      have_unsaved_changes_flag = true; // draw bonds in caller
   }
   return status;
}




/* Reverse the direction of a the fragment of the clicked on
   atom/residue.  A fragment is a consequitive range of residues -
   where there is a gap in the numbering, that marks breaks between
   fragments in a chain.  There also needs to be a distance break - if
   the CA of the next/previous residue is more than 5A away, that also
   marks a break. Thow away all atoms in fragment other than CAs.

   Do it "in place" not return a new molecule.
*/
//
// This is a cheesy implementation that does the fragmentization based
// only on chain id and residue number, not distance.
//
short int
molecule_class_info_t::reverse_direction_of_fragment(const std::string &chain_id,
                                                     int resno) {


   // First find the fragment
   short int istat=0;

   if (atom_sel.n_selected_atoms > 0) {
      // Let's use a minimol:
      coot::minimol::molecule m_initial(atom_sel.mol);
      coot::minimol::molecule fragmented_mol(m_initial.fragmentize());
      short int found_fragment_flag = 0;

      for (unsigned int ifrag=0; ifrag<fragmented_mol.fragments.size(); ifrag++) {
         if (fragmented_mol[ifrag].fragment_id == chain_id) {
            for (int ires=fragmented_mol[ifrag].min_res_no();
                 ires<=fragmented_mol[ifrag].max_residue_number();
                 ires++) {
               if (fragmented_mol[ifrag][ires].atoms.size() > 0) {
                  if (ires == resno) {
                     // OK, we have found our fragment
                     found_fragment_flag = 1;
                     make_backup();
                     istat = 1; // we are doing something.

                     // Lets make a new fragment and replace the
                     // current fragment in fragmented_mol

                     coot::minimol::fragment f = fragmented_mol[ifrag];

                     // find the range of residues that have atoms -
                     // the is a kluge because currently (20050815)
                     // min_res_no() returns low [0,1?] for fragments
                     // that actually start (say) 1005.
                     int fragment_low_resno = fragmented_mol[ifrag].max_residue_number();
                     int fragment_high_resno = -9999;
                     for (int ires_in_frag=fragmented_mol[ifrag].min_res_no();
                          ires_in_frag<=fragmented_mol[ifrag].max_residue_number();
                          ires_in_frag++) {
                        if (fragmented_mol[ifrag][ires_in_frag].atoms.size() > 0) {
                           if (fragmented_mol[ifrag][ires_in_frag].seqnum > fragment_high_resno) {
                              fragment_high_resno = fragmented_mol[ifrag][ires_in_frag].seqnum;
                           }
                           if (fragmented_mol[ifrag][ires_in_frag].seqnum < fragment_low_resno) {
                              fragment_low_resno = fragmented_mol[ifrag][ires_in_frag].seqnum;
                           }
                        }
                     }
                     // OK between fragment_low_resno and
                     // fragment_high_resno we want to reverse the
                     // numbering
                     for (int ires_in_frag=fragment_low_resno;
                          ires_in_frag<=fragment_high_resno;
                          ires_in_frag++) {
                        // pencil and paper job:
                        int new_resno = fragment_high_resno - ires_in_frag + fragment_low_resno;
                        coot::minimol::residue r(new_resno, fragmented_mol[ifrag][ires_in_frag].name);
                        // add the CA
                        for (unsigned int iat=0;
                             iat<fragmented_mol[ifrag][ires_in_frag].atoms.size();
                             iat++) {
                           if (fragmented_mol[ifrag][ires_in_frag].atoms[iat].name == " CA ") {
                              r.addatom(fragmented_mol[ifrag][ires_in_frag].atoms[iat]);
                           }
                        }
                        try {
                           f.addresidue(r, 0);
                        }
                        catch (const std::runtime_error &rte) {
                           std::cout << "ERROR:: auto_fit_best_rotamer() " << rte.what() << std::endl;
                        }
                     }
                     // now replace the fragment in fragmented_mol
                     fragmented_mol[ifrag] = f;
                  }
               }
               if (found_fragment_flag)
                  break;
            }
         }
         if (found_fragment_flag)
            break;
      }
      if (found_fragment_flag) {

         mmdb::Manager *mol = fragmented_mol.pcmmdbmanager();
         // before we get rid of the old atom_sel lets save the cell, symm.
         mmdb::realtype a[6];
         mmdb::realtype vol;
         int orthcode;
         atom_sel.mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
         char *sg = atom_sel.mol->GetSpaceGroup();
         mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);

         if (sg)
            mol->SetSpaceGroup(sg);

         // Now, we can convert that minimol back to a mmdb::Manager
         delete atom_sel.mol;
         atom_sel = make_asc(mol);

         have_unsaved_changes_flag = 1;
         make_bonds_type_checked(__FUNCTION__);
      }
   }
   return istat;
}


// Return 1 on a successful flip.  Flip the last chi angle.
//
int
molecule_class_info_t::do_180_degree_side_chain_flip(const std::string &chain_id,
                                                     int resno,
                                                     const std::string &inscode,
                                                     const std::string &altconf,
                                                     coot::protein_geometry *geom_p) {

   auto check_for_nucleic_acid_atom_names = [] (mmdb::Residue *residue_p) {

      bool status = false;
      mmdb::PPAtom residue_atoms = nullptr;
      int nResidueAtoms = 0;
      residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
      unsigned int n_matched = 0;
      std::vector<std::string> real_names = {" C1'", " C2'", " C3'", " C4'", " O4'", " C5'", " O5'",
	 " P  ", " C2 ", " N1 ", " O1P", " O3'"};
      if (nResidueAtoms > 0) {
	 for(int iat=0; iat<nResidueAtoms; iat++) {
	    mmdb::Atom *at = residue_atoms[iat];
	    std::string atom_name = at->GetAtomName();
	    if (std::find(real_names.begin(), real_names.end(), atom_name) != real_names.end())
	       n_matched++;
	 }
	 if (n_matched > 7) status = true;
      }
      return status;
   };

   int status = 0;

   bool is_protein      = false;
   bool is_nucleic_acid = false;
   mmdb::Residue *residue_p = get_residue(chain_id, resno, inscode);
   if (residue_p) {
      std::string res_name = residue_p->GetResName();
      std::vector<std::string> protein_types =
	 { "GLY", "ALA", "CYS", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU",
	   "MET", "MSE", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP",
	   "TYR"};
      std::vector<std::string> nucleic_acid_types =
	 {"G", "A", "T", "U", "DA", "DG", "DC", "DT"};
      if (std::find(protein_types.begin(), protein_types.end(), res_name) !=
	  protein_types.end())
	 is_protein = true;
      if (std::find(nucleic_acid_types.begin(), nucleic_acid_types.end(), res_name) !=
	  nucleic_acid_types.end())
	 is_nucleic_acid = true;

      if (! is_nucleic_acid)
	 is_nucleic_acid = check_for_nucleic_acid_atom_names(residue_p);

      if (is_protein)
	 return do_180_degree_side_chain_flip_protein(chain_id, resno, inscode, altconf, geom_p);
      if (is_nucleic_acid)
	 return do_180_degree_side_chain_flip_nucleic_acid(chain_id, resno, inscode, altconf, geom_p);
      return status;
   }
   return status;
}

// Return 1 on a successful flip.  Flip the last chi angle.
//
int
molecule_class_info_t::do_180_degree_side_chain_flip_protein(const std::string &chain_id,
							     int resno,
							     const std::string &inscode,
							     const std::string &altconf,
							     coot::protein_geometry *geom_p) {

   // Notice that chi_angles has no concept of alt conf.
   //
   // chi_angles works on the atoms of a residue, with no alt conf
   // checking.
   //
   // So we need to create a synthetic residue (copy) with atoms of
   // the alt conf (and "" if needed).  We flip the residue copy and
   // then feed back the atoms with (altLoc == altconf) into the
   // original residue.

   int istatus=0;
   double diff = 180.0;

   int nth_chi = -1; // unset

   if (atom_sel.n_selected_atoms > 0) {

      int nSelResidues;
      mmdb::PResidue *SelResidues = NULL;
      int selnd = atom_sel.mol->NewSelection();
      atom_sel.mol->Select(selnd, mmdb::STYPE_RESIDUE, 0,
                           chain_id.c_str(),
                           resno, inscode.c_str(),
                           resno, inscode.c_str(),
                           "*", "*", "*", "*", mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selnd, SelResidues, nSelResidues);
      if (nSelResidues > 0 ) {
         mmdb::Residue *residue = SelResidues[0];
         std::string resname = residue->GetResName();

         if (true) {

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

            mmdb::PPAtom residue_atoms = NULL;
            int nResidueAtoms;
            residue->GetAtomTable(residue_atoms, nResidueAtoms);

            if (nth_chi != -1) {
               make_backup();
               mmdb::Residue *residue_copy =
                  coot::util::deep_copy_this_residue_add_chain(residue, altconf, 0, 0);

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

               coot::chi_angles chi_ang(residue_copy, 0);
               std::vector<std::vector<int> > contact_indices(n_atom_residue_copy);
               bool add_reverse_contacts = false;
               contact_indices = coot::util::get_contact_indices_from_restraints(residue_copy, geom_p, 1,
                                                                                 add_reverse_contacts);
               std::pair<short int, float> istat = chi_ang.change_by(nth_chi, diff, contact_indices);

               if (istat.first) { // failure
                  std::cout << "WARNING:: Failure to flip" << std::endl;
               } else {
                  istatus = 1;
                  // OK, we need transfer the coordinates of the
                  // altconfed atoms of residue_copy to residue:
                  //
                  for (int iatc=0; iatc<n_atom_residue_copy; iatc++) {
                     // std::cout << residue_atoms_copy[iat] << std::endl;
                     std::string atom_copy_altconf = residue_atoms_copy[iatc]->altLoc;
                     if (atom_copy_altconf == altconf) {
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

               have_unsaved_changes_flag = 1;
               make_bonds_type_checked(__FUNCTION__);

            } else {
               std::cout << "No flipping allowed for residue type " << resname
                         << std::endl;
            }
         }
      }
      atom_sel.mol->DeleteSelection(selnd);
   }
   return istatus;
}

#include "coot-utils/jed-flip.hh"

// Return 1 on a successful flip.  Flip the last chi angle.
//
int
molecule_class_info_t::do_180_degree_side_chain_flip_nucleic_acid(const std::string &chain_id,
								  int resno,
								  const std::string &inscode,
								  const std::string &altconf,
								  coot::protein_geometry *geom_p) {

   auto atom_names_to_atom_specs = [] (const std::vector<std::string> &atom_names, mmdb::Residue *residue_p) {
      std::vector<coot::atom_spec_t> specs;
      std::string chain_id = residue_p->GetChainID();
      int res_no = residue_p->GetSeqNum();
      std::string ins_code = residue_p->GetInsCode();
      for (const std::string &atom_name : atom_names) {
         specs.push_back(coot::atom_spec_t(chain_id, res_no, ins_code, atom_name, ""));
      }
      return specs;
   };

   int istatus = 0;
   mmdb::Residue *residue_p = get_residue(chain_id, resno, inscode);
   if (residue_p) {
      mmdb::Atom *at_C1prime = residue_p->GetAtom("C1'");
      mmdb::Atom *at_N1      = residue_p->GetAtom("N1");
      mmdb::Atom *at_N9      = residue_p->GetAtom("N9");
      mmdb::Atom *at_N = nullptr;
      if (at_C1prime) {
         coot::Cartesian C1prime_pos(at_C1prime->x, at_C1prime->y, at_C1prime->z);
         coot::Cartesian N_pos;
         bool found_the_N = false;
         bool N_is_N1 = false;
         bool N_is_N9 = false;
         if (at_N1) {
            coot::Cartesian N_posi(at_N1->x, at_N1->y, at_N1->z);
            float dd = coot::Cartesian::lengthsq(C1prime_pos, N_posi);
            if (dd < 10.0) {
               found_the_N = true;
               at_N = at_N1;
               N_pos = N_posi;
               N_is_N1 = true;
            }
         }
         if (at_N9) {
            coot::Cartesian N_posi(at_N9->x, at_N9->y, at_N9->z);
            float dd = coot::Cartesian::lengthsq(C1prime_pos, N_posi);
            if (dd < 10.0) {
               found_the_N = true;
               at_N = at_N1;
               N_pos = N_posi;
               N_is_N9 = true;
            }
         }
         if (found_the_N) {
            mmdb::Manager *mol = atom_sel.mol;
            std::vector<std::string> atom_names;
            if (N_is_N1) atom_names = { " C2'", " C1'", " N1 ", " C6 "};
            if (N_is_N9) atom_names = { " C2'", " C1'", " N9 ", " C8 "};
            std::vector<coot::atom_spec_t> torsion_atoms = atom_names_to_atom_specs(atom_names, residue_p);
            coot::torsion_general tg(residue_p, mol, torsion_atoms);
            Tree tree = tg.GetTree_0_based();
            tg.change_by(180.0, &tree);
            have_unsaved_changes_flag = 1;
            make_bonds_type_checked(__FUNCTION__);
         }
      }
   }
   return istatus;

}


// Return a vector of residues that have missing atoms by dictionary
// search.  missing_hydrogens_flag reflects if we want to count
// residues that have missing hydrogens as residues with missing
// atoms that should be part of the returned vector. Most of the
// time, we don't care about hydrogens and the flag is 0.
//
// We also return a vector of residue names for which we couldn't get
// a geometry dictionary entry.
//
coot::util::missing_atom_info
molecule_class_info_t::missing_atoms(short int missing_hydrogens_flag,
                                     coot::protein_geometry *geom_p) const {

   bool ignore_missing_OXT = true;
   bool ignore_missing_OP3 = true;

   std::vector<mmdb::Residue *> residues_with_missing_atoms;
   std::vector<std::string> residues_no_dictionary;
   std::map<mmdb::Residue *, std::vector<std::string> > residue_missing_atom_names_map;
   // and these atoms will need to be deleted when we auto-complete the residue.
   std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Atom *> > > atoms_in_coords_but_not_in_dict;

   if (atom_sel.n_selected_atoms > 0) {

      // residue_atoms is a vector of residues names (monomer comp_id)
      // together with a vector of atoms.  On each check in this list
      // for the residue type, associated with each atom name is a
      // flag which says if this is a hydrogen or not.
      std::vector<coot::util::dict_residue_atom_info_t> residue_atoms;

      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         mmdb::PResidue residue_p;
         mmdb::Atom *at;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            std::string residue_name(residue_p->GetResName());
            short int found_dict = 0;
            std::vector<coot::util::dict_atom_info_t> residue_dict_atoms;
            for (unsigned int idict=0; idict<residue_atoms.size(); idict++) {
               std::string tmp_str = residue_atoms[idict].residue_name;
               if (residue_name == tmp_str) {
                  residue_dict_atoms = residue_atoms[idict].atom_info;
                  found_dict = 1;
               }
            }
            if (!found_dict) {
               // we need to try to make it from the dict_info restraint info.
               // (dymamic_add will happen automatically).
               //
               // set found_dict if we get it and make it.
               coot::util::dict_residue_atom_info_t residue_atoms_for_a_residue(residue_name,
                                                                                geom_p);
               if (!residue_atoms_for_a_residue.is_empty_p()) {
                  residue_atoms.push_back(residue_atoms_for_a_residue);
                  residue_dict_atoms = residue_atoms_for_a_residue.atom_info;
                  found_dict = 1;
               }
            }

            if (!found_dict) {

               // No dictionary available, even after we try to make
               // it (dynamic add failed, presumably).  So add this
               // residue type to the list of residues for which we
               // can't find the restraint info:
               residues_no_dictionary.push_back(residue_name);

            } else {

               // OK, we have a dictionary. Deal with hydrogens, by
               // assigning them as present if we don't care if they
               // don't exist.
               //
               // Note: this kludges the isHydrogen? flag.  We are
               // testing if it was set (existed) or not!
               //
               std::vector<coot::util::dict_atom_info_t> dict_atom_names_pairs;
               for (unsigned int iat=0; iat<residue_dict_atoms.size(); iat++) {
                  // Let's not add hydrogens to the dictionary atoms
                  // if we don't care if the don't exist.
                  if ((missing_hydrogens_flag == 0) &&
                      (residue_dict_atoms[iat].is_Hydrogen_flag == 1)) {
                     // do nothing
                  } else {
                     // put it in the list and initially mark it as not found.
                     coot::util::dict_atom_info_t p(residue_dict_atoms[iat].name, 0);
                     // PDBv3 FIXME
                     bool really_missing = true;

                     if (ignore_missing_OXT) {
                        if (residue_dict_atoms[iat].name == " OXT")
                           really_missing = false;
                     }
                     if (ignore_missing_OP3) {
                        if (residue_dict_atoms[iat].name == " OP3")
                           really_missing = false;
                     }

                     if (really_missing) {
                        dict_atom_names_pairs.push_back(p);
                     }
                  }
               }

               // OK for every atom in the PDB residue, was it in the dictionary?
               //
               int n_atoms = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms; iat++) {
                  at = residue_p->GetAtom(iat);
                  std::string atom_name(at->name);
                  // check against each atom in the dictionary:
                  for (unsigned int idictat=0; idictat<dict_atom_names_pairs.size(); idictat++) {
                     if (atom_name == dict_atom_names_pairs[idictat].name) {
                        // kludge the is_Hydrogen_flag! Use as a marker
                        dict_atom_names_pairs[idictat].is_Hydrogen_flag = 1; // mark as found
                        break;
                     }
                  }
               }

               // OK, so we have run through all atoms in the PDB
               // residue.  How many atoms in the dictionary list for
               // this residue were not found? Counterintuitive: the
               // is_Hydrogen_flag is used as a marker of being found!
               //
               std::vector<std::string> missing_atom_names;
               for (unsigned int idictat=0; idictat<dict_atom_names_pairs.size(); idictat++) {
                  if (! dict_atom_names_pairs[idictat].is_Hydrogen_flag) {
                     missing_atom_names.push_back(dict_atom_names_pairs[idictat].name);
                  }
               }

               if (! missing_atom_names.empty()) {
                  residues_with_missing_atoms.push_back(residue_p);
                  residue_missing_atom_names_map[residue_p] = missing_atom_names;
               }
            }
         }
      }
   }

   coot::util::missing_atom_info mai(residues_no_dictionary,
                                     residues_with_missing_atoms,
                                     atoms_in_coords_but_not_in_dict);
   mai.residue_missing_atom_names_map = residue_missing_atom_names_map; // a bit kludgy.
   return mai;

}


// The function that uses missing atom info:
coot::util::missing_atom_info
molecule_class_info_t::fill_partial_residues(coot::protein_geometry *geom_p,
                                             int refinement_map_number) {

   coot::util::missing_atom_info info;
   if (atom_sel.n_selected_atoms > 0) {

      make_backup();
      int backup_state = backup_this_molecule;
      backup_this_molecule = 0;  // temporarily
      info = missing_atoms(0, geom_p);

      if (info.residues_with_missing_atoms.size() > 0) {
         std::cout << "INFO:: Residues with missing atoms:" << "\n";
         unsigned int n_per_line = 10;
         for (unsigned int i=0; i<info.residues_with_missing_atoms.size(); i+=n_per_line) {
            for (unsigned int ip=0; ip<n_per_line; ip++) {
               if ((i+ip) < info.residues_with_missing_atoms.size()) {
                  std::cout << info.residues_with_missing_atoms[i+ip]->GetResName() << " "
                            << info.residues_with_missing_atoms[i+ip]->GetSeqNum()  << " "
                            << info.residues_with_missing_atoms[i+ip]->GetChainID() << "  ";
               }
            }
            std::cout << "\n";
         }

         for (unsigned int i=0; i<info.residues_with_missing_atoms.size(); i++) {
            int resno =  info.residues_with_missing_atoms[i]->GetSeqNum();
            std::string chain_id = info.residues_with_missing_atoms[i]->GetChainID();
            std::string residue_type = info.residues_with_missing_atoms[i]->GetResName();
            std::string inscode = info.residues_with_missing_atoms[i]->GetInsCode();
            std::string altloc("");
            float lowest_probability = 0.8;
            int clash_flag = 1;

            mutate(resno, inscode, chain_id, residue_type); // fill missing atoms
            if (refinement_map_number >= 0)
               auto_fit_best_rotamer(ROTAMERSEARCHLOWRES, // backrub rotamers
                                     resno, altloc, inscode, chain_id,
                                     refinement_map_number, clash_flag,
                                     lowest_probability, *geom_p);

            // we could do refinement here, possibly, but instead,
            // return the missing atom info after this function is
            // finished so that the refinement can be done by the
            // calling function, were we can see some visual feedback.
         }
      } else {
         std::cout << " No Residues with missing atoms (that have dictionary entries)\n";
      }

      // restore backup state;
      backup_this_molecule = backup_state;
      have_unsaved_changes_flag = 1;
   }
   return info;
}


int
molecule_class_info_t::fill_partial_residue(const coot::residue_spec_t &residue_spec,
                                            const coot::protein_geometry *geom_p,
                                            int refinement_map_number) {

   int resno = residue_spec.res_no;
   std::string chain_id = residue_spec.chain_id;
   std::string inscode = residue_spec.ins_code;
   std::string altloc = "";

   mmdb::Residue *residue_p = get_residue(chain_id, resno, inscode);
   if (residue_p) {
      std::string residue_type = residue_p->GetResName();
      mutate(resno, inscode, chain_id, residue_type); // fill missing atoms
      if (refinement_map_number >= 0) {
         float lowest_probability = 0.8;
         int clash_flag = 1;
         auto_fit_best_rotamer(ROTAMERSEARCHLOWRES,
                               resno, altloc, inscode, chain_id,
                               refinement_map_number, clash_flag,
                               lowest_probability, *geom_p);
      }
   }
   return 0;
}


// ------------------------------------------------------------------------
//                       dots
// ------------------------------------------------------------------------

void
molecule_class_info_t::draw_dots(Shader *shader_p,
                                 const glm::mat4 &mvp,
                                 const glm::mat4 &view_rotation_matrix,
                                 const std::map<unsigned int, lights_info_t> &lights,
                                 const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                 const glm::vec4 &background_colour,
                                 bool do_depth_fog) {

#if 0 // 20240614-PE olden Coot
   if (! dots.empty())
      for (unsigned int i=0; i<dots.size(); i++)
         if (dots[i].is_open_p() == 1)
            dots[i].imm.draw(shader_p, mvp, view_rotation_matrix, lights, eye_position,
                             background_colour, do_depth_fog);
#endif

}


// return the status of whether or not the dots were cleared.
bool
molecule_class_info_t::clear_dots(int dots_handle) {

   bool r = 0;
   if ((dots_handle >= 0) && (dots_handle < int(dots.size()))) {
      if (dots[dots_handle].is_open_p()) {
         std::cout << "closing dots " << dots_handle << std::endl;
         dots[dots_handle].close_yourself();
         r = 1;
      }
   } else {
      std::cout << "WARNING:: bad dots_handle in clear_dots: "
                << dots_handle << " " << dots.size() << std::endl;
   }
   return r;
}

// close the first open dots object with name dots_object_name.
bool
molecule_class_info_t::clear_dots(const std::string &dots_object_name) {

   bool r = 0;
   for (unsigned int i=0; i<dots.size(); i++) {
      if (dots[i].get_name() == dots_object_name) {
         dots[i].close_yourself();
         r = 1;
         break;
      }
   }
   return r;
}



int
molecule_class_info_t::make_dots(const std::string &atom_selection_str,
                                 const std::string &dots_object_name,
                                 float dot_density, float atom_radius_scale) {


   int dots_handle = -1;

   if (has_model()) {
      int SelHnd = atom_sel.mol->NewSelection(); // yes, deleted.
      atom_sel.mol->Select(SelHnd, mmdb::STYPE_ATOM,
                           atom_selection_str.c_str(),
                           mmdb::SKEY_OR);
      int n_selected_atoms;
      mmdb::PPAtom atom_selection = NULL;
      atom_sel.mol->GetSelIndex(SelHnd, atom_selection, n_selected_atoms);

      gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));

      // dots is declared:    std::vector<coot::dots_representation_info_t> dots;
      pli::dots_representation_info_t dots_info(dots_object_name);
      dots_info.add_dots(SelHnd, atom_sel.mol, NULL, dot_density, dots_colour, dots_colour_set);

      dots.push_back(dots_info);
      dots_handle = dots.size() -1;

      atom_sel.mol->DeleteSelection(SelHnd);
   }
   return dots_handle;
}

// renumber_reidues starting at 1 and removing insertion codes
// (no backup)
void
molecule_class_info_t::simplify_numbering_internal(mmdb::Chain *chain_p) {

   if (chain_p) {
      int n_residues = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<n_residues; ires++) {
         mmdb::Residue *res_p = chain_p->GetResidue(ires);
         int residue_number = ires + 1;
         res_p->SetResID(res_p->name, residue_number, "");
      }
   }
}

int
molecule_class_info_t::renumber_waters() {

   int renumbered_waters = 0;
   if (atom_sel.n_selected_atoms > 0) {

      short int changes_made = 0;
      int n_models = atom_sel.mol->GetNumberOfModels();
      make_backup();
      unsigned int n_solvent_chains = 0;
      for (int imod=1; imod<=n_models; imod++) {

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         if (nchains <= 0) {
            std::cout << "bad nchains in molecule " << nchains
                      << std::endl;
         } else {
            for (int ichain=0; ichain<nchains; ichain++) {
               chain_p = model_p->GetChain(ichain);
               if (chain_p == NULL) {
                  // This should not be necessary. It seem to be a
                  // result of mmdb corruption elsewhere - possibly
                  // DeleteChain in update_molecule_to().
                  std::cout << "WARNING:: renumbered_waters() NULL chain " << ichain << std::endl;
               } else {

                  if (chain_p->isSolventChain()) {

                     n_solvent_chains++;
                     int resno = 1;
                     int nres = chain_p->GetNumberOfResidues();
                     mmdb::PResidue residue_p;
                     for (int ires=0; ires<nres; ires++) {
                        residue_p = chain_p->GetResidue(ires);
                        residue_p->seqNum = resno; // stuff it in
                        changes_made = 1;
                        resno++;  // for next residue
                     }
                  } else {
                     std::string chain_id(chain_p->GetChainID());
                     std::cout << "INFO:: in renumbered_waters() chain " << chain_id << " is not a SolvenChain" << std::endl;
                  }
               }
            }
         }
      }
      if (changes_made) {
         atom_sel.mol->FinishStructEdit();
         have_unsaved_changes_flag = 1;
         renumbered_waters = 1;
      }

      // maybe return this so it can be displayed in the GUI? 20210727-PE FIXME
      if (n_solvent_chains == 0)
         std::cout << "WARNING:: no SolventChains in the model " << std::endl;
   }
   return renumbered_waters;
}

std::vector<std::string>
molecule_class_info_t::get_symop_strings() const {

   std::vector<std::string> r;
   if (atom_sel.mol) {
      // coords
      mmdb::mat44 my_matt;
      int ierr = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
      if (ierr == 0) {
         // Good, we have symm info
         for (int i = 0; i < atom_sel.mol->GetNumberOfSymOps(); i++) {
            char* symop_str = atom_sel.mol->GetSymOp(i);
            r.push_back(symop_str);
         }
      }
   } else {
      // map
      int n = xmap.spacegroup().num_symops();
      for (int i=0; i<n; i++) {
         r.push_back(xmap.spacegroup().symop(i).format());
      }
   }
   return r;
}


int
molecule_class_info_t::make_map_from_cns_data(const clipper::Spacegroup &sg,
                                              const clipper::Cell &cell,
                                              std::string cns_data_filename) {


   clipper::CNS_HKLfile cns;
   cns.open_read(cns_data_filename);

   clipper::Resolution resolution(cns.resolution(cell));

   clipper::HKL_info mydata(sg, cell, resolution);
   clipper::HKL_data<clipper::datatypes::F_phi<float>   >  fphidata(mydata);

   std::cout << "importing info" << std::endl;
   cns.import_hkl_info(mydata);
   std::cout << "importing data" << std::endl;
   cns.import_hkl_data(fphidata);
   cns.close_read();

   std::string mol_name = cns_data_filename;

   initialize_map_things_on_read_molecule(mol_name, false, false, false); // not diff map

   std::cout << "initializing map...";
   xmap.init(mydata.spacegroup(),
                     mydata.cell(),
                     clipper::Grid_sampling(mydata.spacegroup(),
                                            mydata.cell(),
                                            mydata.resolution()));
   std::cout << "done."<< std::endl;
   std::cout << "doing fft..." ;
   xmap.fft_from( fphidata );                  // generate map
   std::cout << "done." << std::endl;
   update_map_in_display_control_widget();

   mean_and_variance<float> mv = map_density_distribution(xmap,0);

   std::cout << "Mean and sigma of map from CNS file: " << mv.mean
             << " and " << sqrt(mv.variance) << std::endl;

   // fill class variables
   map_mean_ = mv.mean;
   map_sigma_ = sqrt(mv.variance);

   xmap_is_diff_map = 0;
   contour_level = nearest_step(mv.mean + 1.5*sqrt(mv.variance), 0.05);

   // what is this graphics_info_t doing here.
//    graphics_info_t g;
//    // g.scroll_wheel_map = g.n_molecules(); // change the current scrollable map.
//    update_map();
//    return g.n_molecules;

   return imol_no;
}

clipper::Coord_orth
molecule_class_info_t::find_peak_along_line(const clipper::Coord_orth &p1,
                                            const clipper::Coord_orth &p2) const {

   float high_point_1 = -9999999.9;
   clipper::Coord_orth pbest;
   int istep_max = 500;

   for (int istep=0; istep<=istep_max; istep++) {
      float fr = float(istep)/float(istep_max);
      clipper::Coord_orth pc = p1 + fr*(p2-p1);
      float d = density_at_point(pc);
      // std::cout << ": " << istep << " " << fr << " " << d  << std::endl;
      if (d > high_point_1) {
         high_point_1 = d;
         pbest = pc;
      }
   }
   return pbest;
}

// Throw an exception if there is no point about the contour level of this map.
//
// This is a different algorithm to above.  Here we find the first
// blob about the contor level of the map.  The actual peak height
// does not matter.  (i.e. low level peaks at front will win over big
// peaks behind (no matter what the big peak density level is).
//
// So now we find the first peak about the contour level and find the
// highest point in that peak.  When that peak goes below the contour
// level we stop searching for the higest peak along the line.
//
clipper::Coord_orth
molecule_class_info_t::find_peak_along_line_favour_front(const clipper::Coord_orth &p1,
                                                         const clipper::Coord_orth &p2) const {

   float best_score = -9999999.9;
   clipper::Coord_orth pbest;
   int istep_max = 500;
   bool point_set = 0;

   for (int istep=0; istep<=istep_max; istep++) {
      float fr = float(istep)/float(istep_max);
      clipper::Coord_orth pc = p1 + fr*(p2-p1);
      float d = density_at_point(pc);
      if (d > contour_level) {
         // OK, so the point we want is somewhere in this peak
         for (int jstep=istep; jstep<=istep_max; jstep++) {
            fr = float(jstep)/float(istep_max);
            pc = p1 + fr*(p2-p1);
            d = density_at_point(pc);
            if (d > contour_level) {
               if (d> best_score) {
                  best_score = d;
                  pbest = pc;
                  point_set = 1;
               }
            } else {
               // the front peak is over (now below the contour level), we have the pbest.
               break;
            }
         }
         break;
      }
   }

   if (! point_set) {
      std::string mess("No peak above ");
      mess += coot::util::float_to_string(contour_level);
      mess += " found.";
      throw std::runtime_error(mess);
   }
   return pbest;
}

// replace molecule
int
molecule_class_info_t::replace_molecule(mmdb::Manager *mol) {

   int was_changed = 0;
   if (mol) {
      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
      delete atom_sel.mol;
      atom_sel = make_asc(mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
      trim_atom_label_table();
      update_symmetry();
      was_changed = 1;
   }

   // deubugging
   if (0) {
      std::cout << "DEBUG:: in replace_molecule was_changed " << was_changed << std::endl;
      std::cout << "DEBUG:: n_atoms:  " << atom_sel.n_selected_atoms << std::endl;

      std::cout << atom_sel.mol->GetNumberOfModels() << " models" << std::endl;
      for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
         std::cout <<  "   model " << imod << std::endl;
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         std::cout << "   " << nchains << " chains" << std::endl;
         for (int ichain=0; ichain<nchains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            std::cout <<  "      \"" << chain_p->GetChainID() << "\"" << std::endl;
            int nres = chain_p->GetNumberOfResidues();
            std::cout << "      " << nres << " residues" << std::endl;
            mmdb::Residue *residue_p;
            mmdb::Atom *at;
            for (int ires=0; ires<nres; ires++) {
               residue_p = chain_p->GetResidue(ires);
               std::cout <<  "         Residue: " << residue_p->GetSeqNum() << "" << std::endl;
               int n_atoms = residue_p->GetNumberOfAtoms();
               std::cout << "         " << n_atoms << " atoms" << std::endl;
               for (int iat=0; iat<n_atoms; iat++) {
                  at = residue_p->GetAtom(iat);
                  std::cout << "            " << at << std::endl;
               }
            }
         }
      }
   }
   return was_changed;
}

int
molecule_class_info_t::replace_models(std::deque<mmdb::Model *> model_list) {
   /*This function is nearly identical to replace_molecule(), but it replaces all MMDB models
     within the existing MMDB manager instead of replacing the entire MMDB manager.  This means
     that the header information associated with the current MMDB manager won't get erased.
   */
   int was_changed = 0;

   if (!model_list.empty()) {

      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);

      mmdb::Manager *mol = atom_sel.mol;
      // mol->DeleteAllModels(); // Crash 20120105
      mol->Delete(mmdb::MMDBFCM_Coord);
      while (!model_list.empty()) {
         mol->AddModel(model_list.front());
         model_list.pop_front();
      }

      mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      mol->FinishStructEdit();

      atom_sel = make_asc(mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
      trim_atom_label_table();
      update_symmetry();
      was_changed = 1;
   }

   return was_changed;
}


// EM map function
int
molecule_class_info_t::scale_cell(float fac_u, float fac_v, float fac_w) {

   int retval = 0;  // NXMAP-FIXME

   if (has_xmap()) {
      clipper::Cell cell_orig = xmap.cell();
      clipper::Cell_descr cell_d(cell_orig.a() * fac_u,
                                 cell_orig.b() * fac_v,
                                 cell_orig.c() * fac_w,
                                 cell_orig.alpha(),
                                 cell_orig.beta(),
                                 cell_orig.gamma());

      clipper::Cell new_cell(cell_d);
      clipper::Spacegroup new_spg(xmap.spacegroup());
      clipper::Xmap_base::Map_reference_index ix;
      clipper::Grid_sampling gs = xmap.grid_sampling();
      clipper::Xmap<float> new_map(new_spg, new_cell, gs);
      for (ix = xmap.first(); !ix.last(); ix.next() ) {
          new_map[ix] = xmap[ix];
      }
      xmap = new_map;
      update_map(true);
   }
   return retval;
}

// number of chains. Return -1 on failure
int
molecule_class_info_t::number_of_chains() const {

   int n = -1;
   if (has_model()) {
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         n = model_p->GetNumberOfChains();
      }
   }
   return n;
}


void
molecule_class_info_t::sort_chains() {


   if (atom_sel.mol)
      coot::sort_chains(atom_sel.mol);

}

void
molecule_class_info_t::sort_residues() {

   if (atom_sel.mol)
      // this calls clean up and FinishStructEdit()
      coot::sort_residues(atom_sel.mol);

}


std::vector<coot::residue_spec_t>
molecule_class_info_t::residues_near_residue(const coot::residue_spec_t &rspec, float radius) const {

   std::vector<coot::residue_spec_t> r =
      coot::residues_near_residue(rspec, atom_sel.mol, radius);

   return r;
}


// return the resulting torsion value.  In degrees.  The input angle tors is in degrees.
//
double
molecule_class_info_t::set_torsion(const std::string &chain_id,
                                   int res_no,
                                   const std::string &insertion_code,
                                   const std::string &alt_conf,
                                   const std::string &atom_name_1,
                                   const std::string &atom_name_2,
                                   const std::string &atom_name_3,
                                   const std::string &atom_name_4,
                                   double tors,
                                   const coot::protein_geometry &geom) {

   double r = -999.9;

   mmdb::Residue *residue_p = get_residue(chain_id, res_no, insertion_code);

   if (!residue_p) {
      std::cout << "WARNING:: failed to get residue with specs :"
                << chain_id << ": " << res_no << " :" << insertion_code
                << ":" << std::endl;
   } else {
      std::string res_name(residue_p->GetResName());
      std::pair<bool, coot::dictionary_residue_restraints_t> restraints_info =
         geom.get_monomer_restraints(res_name, imol_no);
      if (! restraints_info.first) {
         std::cout << "WARNING:: set_torsion: No restraints for " << res_name << std::endl;
      } else {
         make_backup();
         coot::atom_tree_t tree(restraints_info.second, residue_p, alt_conf);

         try {
            r = tree.set_dihedral(atom_name_1,
                                  atom_name_2,
                                  atom_name_3,
                                  atom_name_4,
                                  tors);
            atom_sel.mol->FinishStructEdit();
            make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
            have_unsaved_changes_flag = 1;
         }
         catch(const std::runtime_error &rte) {
            std::cout << "in set_torsion:: set_dihedral() error: " << rte.what() << std::endl;
         }
      }
   }
   return r;
}


// manipulate torsion angles of first residue in the molecule to
// match those of the passed (reference residue (from a different
// molecule, typically).
//
int
molecule_class_info_t::match_torsions(mmdb::Residue *res_reference,
                                      const std::vector <coot::dict_torsion_restraint_t> &tr_ref_res,
                                      const coot::protein_geometry &geom) {


   int n_torsions_moved = 0;
   make_backup();

   mmdb::Residue *res_ligand = coot::util::get_first_residue(atom_sel.mol); // this could/should be replaced
                                                                       // by something that allows
                                                                       // any residue in the molecule
                                                                       // to move to match the
                                                                       // reference residue.

   if (res_ligand) { // the local (moving) residue is xxx_ligand
      std::string res_name_ligand(res_ligand->GetResName());
      std::pair<bool, coot::dictionary_residue_restraints_t> ligand_restraints_info =
         geom.get_monomer_restraints(res_name_ligand, imol_no);
      if (ligand_restraints_info.first) {
         std::vector <coot::dict_torsion_restraint_t> tr_ligand =
            geom.get_monomer_torsions_from_geometry(res_name_ligand, imol_no, 0);
         if (tr_ligand.size()) {

            // find the matching torsion between res_ligand and res_reference and then
            // set the torsions of res_ligand to match those of res_reference.
            //
            // moving the res_ligand
            coot::match_torsions mt(res_ligand, res_reference, ligand_restraints_info.second);
            n_torsions_moved = mt.match(tr_ligand, tr_ref_res);
            atom_sel.mol->FinishStructEdit();
            make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
            have_unsaved_changes_flag = 1;
         } else {
            std::cout << "WARNING torsion restraints of ligand: size 0" << std::endl;
         }
      } else {
         std::cout << "WARNING ligand_restraints_info.first failed " << std::endl;
      }
   } else {
      std::cout << "WARNING:: null ligand residue (trying to get first) " << std::endl;
   }
   return n_torsions_moved;
}

void
molecule_class_info_t::lsq_improve(mmdb::Manager *mol_ref, const std::string &ref_selection_str,
                                   const std::string &moving_selection_str,
                                   int n_res, float dist_crit) {
   if (mol_ref) {
 
      try {
         make_backup();
         coot::lsq_improve lsq_imp(mol_ref, ref_selection_str, atom_sel.mol, moving_selection_str);
         lsq_imp.improve();
         clipper::RTop_orth rtop = lsq_imp.rtop_of_moving();
         std::cout << "rtop:\n" << rtop.format() << std::endl;
         coot::util::transform_mol(atom_sel.mol, rtop);
         have_unsaved_changes_flag = 1;
         make_bonds_type_checked(__FUNCTION__); // calls update_ghosts()
      }
      catch (const std::runtime_error &rte) {
         std::cout << "lsq_improve ERROR::" << rte.what() << std::endl;
      }
   }
}

// this function is fill_chiral_volume_outlier_markers()
void
molecule_class_info_t::fill_chiral_volume_outlier_marker_positions(int state) {

   double chiral_volume_limit_for_outlier = 6.0;
   chiral_volume_outlier_marker_positions.clear();
   if (state) {
      if (atom_sel.mol) {
         std::pair<std::vector<std::string> , std::vector<std::pair<coot::atom_spec_t, double> > > dcv =
            coot::distorted_chiral_volumes(imol_no, atom_sel.mol, graphics_info_t::Geom_p(),
                                           graphics_info_t::cif_dictionary_read_number,
                                           chiral_volume_limit_for_outlier);
         for (unsigned int i=0; i<dcv.second.size(); i++) {
            const auto &atom_spec = dcv.second[i].first;
            mmdb:: Atom *at = get_atom(atom_spec);
            if (at) {
               glm::vec3 p(at->x, at->y, at->z);
               // std::cout << "set_show_chiral_volume_outlier_markers() adding atom " << at << std::endl;
               chiral_volume_outlier_marker_positions.push_back(p);
            }
         }
      }
   }
}


void
molecule_class_info_t::set_show_non_bonded_contact_baddies_markers(int state) {

}
