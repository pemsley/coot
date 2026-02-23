/* coords/Bond_lines.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2009 by the University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
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

#include <string.h> // for strcmp
#include <stdlib.h> // for labs
#include <string>
#include <fstream>
#include <vector>
#include <algorithm> // for find
#include <set>
#include <iterator>

#include "Cartesian.hh"
#include <mmdb2/mmdb_manager.h>
#include "mmdb-extras.hh"
#include "mmdb.hh"
#include "mmdb-crystal.hh"  // should be merged with extras

#include "geometry/lbg-graph.hh"  // aromatic ring systems
#include "geometry/mol-utils.hh"
#include "utils/coot-utils.hh" // for int_to_string

#include "Bond_lines.hh"
#include "coot-utils/coot-coord-utils.hh"

#include "geometry/protein-donor-acceptors.hh"
#include "loop-path.hh"

#include "utils/logging.hh"
extern logging logger;

static std::string b_factor_bonds_scale_handle_name = "B-factor-bonds-scale";

unsigned int
coot::my_atom_colour_map_t::index_for_chain(const std::string &chain_id) {

   unsigned int isize = atom_colour_map.size();
   // counting backwards would be faster/better
   for (unsigned int i=0; i<isize; i++) {
      if (atom_colour_map[i] == chain_id) {
         return i;
      }
   }

   // add a new one then
   atom_colour_map.push_back(chain_id);

   if (isize == HYDROGEN_GREY_BOND) {
      atom_colour_map[isize] = "skip-hydrogen-grey-colour-for-chain";
      atom_colour_map.push_back(chain_id);
      isize++;
   }
   return isize;
}

void
coot::my_atom_colour_map_t::fill_chain_id_map(const atom_selection_container_t &SelAtom) {

   int imod = 1;
   mmdb::Model *model_p = SelAtom.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         if (n_res > 0) {
            std::string chain_id(chain_p->GetChainID());
            unsigned int idx = index_for_chain(chain_id);
         }
      }
   }
}


Bond_lines::Bond_lines(const graphics_line_t &line) {

   points.push_back(line);
}

// we can put other things here
void
Bond_lines_container::init() {

   // std::cout << "init() " << std::endl;
   rotamer_probability_tables_p = NULL;
   do_sticks_for_waters = false;
   use_deuteranomaly_mode = false;
   n_atoms_in_atom_selection = 0;
   geom = nullptr;
}


void
Bond_lines::update(mmdb::Atom **atom_selection, int n_atoms) {

   for (auto &bond : points) // "points" doesn't seem to be a good name
      bond.update(atom_selection, n_atoms);
}


// Constructor A
//
// We arrange things like this because the other constructor now uses
// construct_from_asc() too.
//
// This is a tiny bit confusing having so many constructors.
// Heyho - historical cruft.
//
Bond_lines_container::Bond_lines_container(const atom_selection_container_t &SelAtom,
                                           int imol,
                                           int do_disulphide_bonds_in,   // default argument, 0
                                           int do_bonds_to_hydrogens_in, // default argument, 1
                                           bool do_rama_markup,       // default argument false
                                           bool do_rota_markup,       // default argument false
                                           coot::rotamer_probability_tables *tables_p
                                           ) {
   // teehee
   // std::cout << "################################## yes this constructor A ###########################" << std::endl;

   if (false)
      std::cout << "DEBUG:: constructor A was called" << std::endl;

   init();
   do_disulfide_bonds_flag = do_disulphide_bonds_in;
   do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
   b_factor_scale = 1.0;
   geom = 0;
   have_dictionary = 0;
   for_GL_solid_model_rendering = 1;
   n_atoms_in_atom_selection = SelAtom.n_selected_atoms;
   if (tables_p)
      rotamer_probability_tables_p = tables_p;

   // 1.7 will not catch MET bonds (1.791 and 1.803) nor MSE bonds (1.95)
   // but SO4 bonds (1.46 are fine).
   // They should have special case, handle_MET_or_MSE_case
   // However, for VNP thingy, S1 has bonds to carbons of 1.67 1.77.  Baah.
   float max_dist = 1.71;
   int model_number = 0; // all models

   construct_from_asc(SelAtom, imol, 0.01, max_dist, coot::COLOUR_BY_ATOM_TYPE, 0, model_number,
                      do_rama_markup, do_rota_markup);
   verbose_reporting = 0;
   udd_has_ca_handle = -1;

}



// the constructor for bond by dictionary - should use this most of the time.
// geom_in can be null if you don't have it.
//
// if model_number is 0, display all models. If it is not 0 then
// display only the given model_number (if possible, of course).
//
// This one for intermediate atoms too.
//
// do_rama_markup, do_rota_markup, rotamer_probability_tables_p are default
// arguments, false/null.
//
Bond_lines_container::Bond_lines_container(const atom_selection_container_t &SelAtom,
                                           int imol,
                                           const std::set<int> &no_bonds_to_these_atoms_in,
                                           const coot::protein_geometry *geom_in,
                                           int do_disulphide_bonds_in,
                                           int do_bonds_to_hydrogens_in,
                                           bool draw_missing_loops_flag,
                                           int model_number,
                                           std::string dummy,
                                           bool do_rama_markup,
                                           bool do_rota_markup,
                                           bool do_sticks_for_waters_in,
                                           coot::rotamer_probability_tables *tables_p) : no_bonds_to_these_atoms(no_bonds_to_these_atoms_in) {


   if (false) {
      std::cout << "############# Bond_lines_container::Bond_lines_container() B "
                << SelAtom.mol << " " << no_bonds_to_these_atoms.size() << std::endl;
      std::cout << "   no_bonds_to_these_atoms set size: " << no_bonds_to_these_atoms.size() << std::endl;
      std::cout << "--------------------- no bonds to these atoms ------------------"  << std::endl;
      std::set<int>::const_iterator it;
      for(it=no_bonds_to_these_atoms.begin(); it!=no_bonds_to_these_atoms.end(); it++) {
         mmdb::Atom *at = SelAtom.atom_selection[*it];
         std::cout << "   Bond_lines_container constructor B: " << *it << " " << coot::atom_spec_t(at) << std::endl;
      }
   }

   if (false)
      std::cout << "DEBUG:: constructor B was passed model_number " << model_number << std::endl;

   do_disulfide_bonds_flag = do_disulphide_bonds_in;
   do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
   for_GL_solid_model_rendering = 1;
   b_factor_scale = 1.0;
   have_dictionary = 0;
   init();
   n_atoms_in_atom_selection = SelAtom.n_selected_atoms;
   if (tables_p)
      rotamer_probability_tables_p = tables_p;
   if (geom_in) {
      geom = geom_in;
      have_dictionary = 1;
   }

   do_sticks_for_waters = do_sticks_for_waters_in;

   // 1.7 will not catch MET bonds (1.791 and 1.803) nor MSE bonds (1.95)
   // but SO4 bonds (1.46 are fine).
   // They should have special case, handle_MET_or_MSE_case
   // However, for VNP thingy, S1 has bonds to carbons of 1.67 1.77.  Baah.
   float max_dist = 1.71;

   // If *every* atom is excluted, e.g. a ligand - or chain-refine
   // then do nothing
   unsigned int n_selected_atoms = SelAtom.n_selected_atoms; // signedness change
   if (n_selected_atoms == no_bonds_to_these_atoms.size()) {
      // do nothing
   } else {
      // sphere refine or active-residue marker - or something else

      construct_from_asc(SelAtom, imol, 0.01, max_dist, coot::COLOUR_BY_ATOM_TYPE, 0,
                         draw_missing_loops_flag,
                         model_number,
                         do_rama_markup, do_rota_markup);
   }
   verbose_reporting = 0;
   udd_has_ca_handle = -1;
}

Bond_lines_container::Bond_lines_container(atom_selection_container_t SelAtom,
                                           int imol,
                                           float max_dist) {

   verbose_reporting = 0;
   do_disulfide_bonds_flag = 1;
   udd_has_ca_handle = -1;
   do_bonds_to_hydrogens = 1;
   b_factor_scale = 1.0;
   have_dictionary = 0;
   for_GL_solid_model_rendering = 1;
   geom = 0;
   init();
   n_atoms_in_atom_selection = SelAtom.n_selected_atoms;
   int model_number = 0; // all models
   bool do_rama_markup = false;
   construct_from_asc(SelAtom, imol, 0.01, max_dist, coot::COLOUR_BY_ATOM_TYPE, 0, model_number,
                      do_rama_markup);
}


Bond_lines_container::Bond_lines_container(atom_selection_container_t SelAtom,
                                           int imol,
                                           float min_dist, float max_dist) {

   verbose_reporting = 0;
   do_disulfide_bonds_flag = 1;
   udd_has_ca_handle = -1;
   do_bonds_to_hydrogens = 1;
   b_factor_scale = 1.0;
   have_dictionary = 0;
   for_GL_solid_model_rendering = 1;
   init();
   // 0 is is_from_symmetry_flag
   int model_number = 0; // all models
   bool do_rama_markup = false;
   n_atoms_in_atom_selection = SelAtom.n_selected_atoms;
   construct_from_asc(SelAtom, imol, min_dist, max_dist, coot::COLOUR_BY_ATOM_TYPE, 0, model_number,
                      do_rama_markup);
}

// geom_in can be null.
//
// The constructor for ball and stick, this constructor implies that
// for_GL_solid_model_rendering is set.
//
Bond_lines_container::Bond_lines_container(atom_selection_container_t asc,
                                           int imol,
                                           const coot::protein_geometry *geom_in) {

   for_GL_solid_model_rendering = 1; // note!

   verbose_reporting = 0;
   do_disulfide_bonds_flag = 1;
   udd_has_ca_handle = -1;
   do_bonds_to_hydrogens = 1;
   b_factor_scale = 1.0;
   have_dictionary = 0;
   init();
   n_atoms_in_atom_selection = asc.n_selected_atoms;
   if (geom_in) {
      geom = geom_in;
      have_dictionary = 1;
   }
   // 0 is is_from_symmetry_flag
   int model_number = 0; // all models
   bool do_rama_markup = false;
   construct_from_asc(asc, imol, 0.01, 1.9, coot::COLOUR_BY_ATOM_TYPE, 0, model_number,
                      do_rama_markup);
}


// This is the one for occupancy and B-factor representation - and now
// all-atom user-define colouring too
//
Bond_lines_container::Bond_lines_container(const atom_selection_container_t &SelAtom,
                                           int imol,
                                           const coot::protein_geometry *protein_geom,
                                           Bond_lines_container::bond_representation_type br_type) {

   // std::cout << "*************************** Bond_lines_container() constructor with geom and type " << br_type << std::endl;

   init(); // sets geom to null pointer
   verbose_reporting = 0;
   do_disulfide_bonds_flag = 1;
   udd_has_ca_handle = -1;
   do_bonds_to_hydrogens = 1;
   b_factor_scale = 1.0;
   geom = protein_geom;
   have_dictionary = 0;
   if (geom) have_dictionary = 1;
   for_GL_solid_model_rendering = 1;
   n_atoms_in_atom_selection = SelAtom.n_selected_atoms;
   float max_dist = 1.71;
   int model_number = 0; // all models
   bool do_rama_markup = false;
   if (br_type == Bond_lines_container::COLOUR_BY_OCCUPANCY) {
      construct_from_asc(SelAtom, imol, 0.01, max_dist, coot::COLOUR_BY_OCCUPANCY, 0, model_number,
                         do_rama_markup);
   } else {
      if (br_type == Bond_lines_container::COLOUR_BY_B_FACTOR) {
         set_b_factor_colours(SelAtom.mol);
              try_set_b_factor_scale(SelAtom.mol);
              construct_from_asc(SelAtom, imol, 0.01, max_dist, coot::COLOUR_BY_B_FACTOR, 0, model_number,
                                 do_rama_markup);
      } else {
              if (br_type == Bond_lines_container::COLOUR_BY_USER_DEFINED_COLOURS) {
                 construct_from_asc(SelAtom, imol, 0.01, max_dist, coot::COLOUR_BY_USER_DEFINED_COLOURS, 0,
                                    model_number, do_rama_markup);
         }
      }
   }
}

void
Bond_lines_container::try_set_b_factor_scale(mmdb::Manager *mol) {

   int udd_b_factor_handle =  mol->GetUDDHandle(mmdb::UDR_HIERARCHY,
                                                coot::b_factor_bonds_scale_handle_name.c_str());

   if (udd_b_factor_handle > 0) {
      mmdb::realtype scale;
      if (mol->GetUDData(udd_b_factor_handle, scale) == mmdb::UDDATA_Ok) {
         b_factor_scale = scale;
      }
   }
}

void
Bond_lines_container::construct_from_atom_selection(const atom_selection_container_t &asc,
                                                    const mmdb::PPAtom atom_selection_1,
                                                    int n_selected_atoms_1,
                                                    const mmdb::PPAtom atom_selection_2,
                                                    int n_selected_atoms_2,
                                                    int imol,
                                                    float min_dist, float max_dist,
                                                    int atom_colour_type,
                                                    bool are_different_atom_selections,
                                                    bool have_udd_atoms,
                                                    int udd_done_bond_handle) {

   n_atoms_in_atom_selection = asc.n_selected_atoms;
   mmdb::Contact *contact = NULL;
   int ncontacts = 0;
   long i_contact_group = 1;

   // matrix stuff
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;

   // update my_matt;  You can't do this if you haven't set the space group.
   //
   // symm.GetTMatrix(my_matt, 0);

   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;

   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;


   // cout << "my_matt is: " << my_matt << endl; // Argh! doesnt work.
   //
//    cout << "my_matt is: " << endl
//         << my_matt[0][0] << " "  << my_matt[0][1] << " "
//         << my_matt[0][2] << " "  << my_matt[0][3] << " "  << endl
//         << my_matt[1][0] << " "  << my_matt[1][1] << " "
//         << my_matt[1][2] << " "  << my_matt[1][3] << " "  << endl
//         << my_matt[2][0] << " "  << my_matt[2][1] << " "
//         << my_matt[2][2] << " "  << my_matt[2][3] << " "  << endl
//         << my_matt[3][0] << " "  << my_matt[3][1] << " "
//         << my_matt[3][2] << " "  << my_matt[3][3] << " "  << endl;

//    std::cout << "Here are the atoms for which we will seek contacts\n";
//    for (int ii=0; ii<n_selected_atoms_1; ii++) {
//       std::cout << atom_selection_1[ii] << std::endl;
//    }

   // Note, it has happened a couple of times now, when we get a crash
   // in mmdb's MakeBricks (or something like that) from here, that's
   // because we are passing an atom that has a nan for a coordinate.

   if (false) {
      std::cout << "Seeking contact: selection 1 " << std::endl;
      for (int ii=0; ii<n_selected_atoms_1; ii++)
              std::cout << "   " << ii << " " << atom_selection_1[ii] << " :"
                             << atom_selection_1[ii]->isTer() << ":" << std::endl;
      std::cout << "Seeking contact: selection 2 " << std::endl;
      for (int ii=0; ii<n_selected_atoms_2; ii++)
              std::cout << "   " << ii << " " << atom_selection_2[ii] << " :"
                             << atom_selection_2[ii]->isTer() << ":" << std::endl;
   }

   asc.mol->SeekContacts(atom_selection_1, n_selected_atoms_1,
                         atom_selection_2, n_selected_atoms_2,
                         min_dist, max_dist, // min, max distances
                         0,        // seqDist 0 -> in same res also
                         contact, ncontacts,
                         0, &my_matt, i_contact_group);

   std::string element_1, element_2;
   int col; // atom colour

   int udd_atom_index_handle = asc.UDDAtomIndexHandle; // so that we don't draw bonds for ligands when we have
                                                       // intermediate atoms displayed.

   int udd_user_defined_atom_colour_index_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   if (ncontacts == 0) {

      // SF4 etc

   } else {

      if (contact) {

         if (false)
            for (int i=0; i< ncontacts; i++)
               std::cout << "contact " << i << " " << contact[i].id1 << " " << contact[i].id2
                         << std::endl;

         std::vector<std::pair<bool, mmdb::Residue *> > het_residues; // bond these separately.
         std::vector<std::pair<bool, mmdb::Residue *> > hoh_residues; // that have O and Hs.

         for (int i=0; i< ncontacts; i++) {

            if (are_different_atom_selections || (contact[i].id2 > contact[i].id1) ) {

               // Added 20171224. Are these correct or do I need to use UDData?
               int atom_index_1 = contact[i].id1;
               int atom_index_2 = contact[i].id2;

               mmdb::Atom *atom_p_1 = atom_selection_1[ contact[i].id1 ];
               mmdb::Atom *atom_p_2 = atom_selection_2[ contact[i].id2 ];

               // let's try UDD atom index
               int atom_index_1_save = atom_index_1;
               int atom_index_2_save = atom_index_2;
               int ierr_1 = atom_p_1->GetUDData(asc.UDDAtomIndexHandle, atom_index_1);
               int ierr_2 = atom_p_2->GetUDData(asc.UDDAtomIndexHandle, atom_index_2);

               // 20220209-PE why does GetUDData fail for intermediate atoms?

               if (true) {
                  // This happens for all intermediate atoms
                  if (ierr_1 != mmdb::UDDATA_Ok) {
                     // std::cout << "           Fail udd " << coot::atom_spec_t(atom_p_1) << std::endl;
                     atom_index_1 = atom_index_1_save;
                  }
                  if (ierr_2 != mmdb::UDDATA_Ok) {
                     // std::cout << "           Fail udd " << coot::atom_spec_t(atom_p_2) << std::endl;
                     atom_index_2 = atom_index_2_save;
                  }
               }

               if (false) {
                  // This happens for all intermediate atoms
                  if (ierr_1 != mmdb::UDDATA_Ok)
                     std::cout << "           Fail udd " << coot::atom_spec_t(atom_p_1) << std::endl;
                  if (ierr_2 != mmdb::UDDATA_Ok)
                     std::cout << "           Fail udd " << coot::atom_spec_t(atom_p_2) << std::endl;
               }

               std::string chain_id1(atom_p_1->GetChainID());
               std::string chain_id2(atom_p_2->GetChainID());

               std::string aloc_1(atom_p_1->altLoc);
               std::string aloc_2(atom_p_2->altLoc);

               element_1 = atom_p_1->element;
               element_2 = atom_p_2->element;

               coot::Cartesian atom_1_pos(atom_p_1->x, atom_p_1->y, atom_p_1->z);
               coot::Cartesian atom_2_pos(atom_p_2->x, atom_p_2->y, atom_p_2->z);

               if (chain_id1 == chain_id2) {

                  // alternate location test
                  //
                  if ( (aloc_1=="") || (aloc_2=="") || (aloc_1==aloc_2) ) {

                     int res_1 = atom_p_1->GetSeqNum();
                     int res_2 = atom_p_2->GetSeqNum();

                     bool bond_het_residue_by_dictionary =
                        add_bond_by_dictionary_maybe(imol, atom_p_1, atom_p_2, &het_residues); // add to het_residues maybe
                     if (false)
                        std::cout << atom_p_1 <<  " " << atom_p_2 << " bonded by dictionary: "
                                  << bond_het_residue_by_dictionary << std::endl;

                     if (bond_het_residue_by_dictionary) {

                        std::string res_name = atom_p_1->GetResName();
                        if (res_name == "HOH" || res_name == "DOD")
                           add_bond_by_dictionary_maybe(imol, atom_p_1, atom_p_2, &hoh_residues);

                     } else {

                        // this +/- 1 residue test, or are DUM atoms.

                        bool is_neighbour = false;
                        if (labs(res_1 - res_2) < 2)
                           is_neighbour = true;
                        if (! is_neighbour)
                           if (labs(atom_p_1->residue->index - atom_p_2->residue->index) < 2)
                              is_neighbour = true;

                        // Maybe DUM-DUM needs it's own bonding selection and drawing function

                        //                      if (! is_neighbour)
                        //                         if (strncmp(atom_p_1->name, "DUM", 3))
                        //                            if (strncmp(atom_p_2->name, "DUM", 3))
                        //                               is_neighbour = true;

                        if (is_neighbour) {

                           //                    std::cout << "Adding bond " << atom_selection_1[ contact[i].id1 ]
                           //                              << " to "
                           //                              << atom_selection_2[ contact[i].id2 ] << std::endl;

                           if (atom_selection_1[ contact[i].id1 ]->GetModel() ==
                               atom_selection_2[ contact[i].id2 ]->GetModel()) {

                              int imodel = atom_selection_1[ contact[i].id1 ]->GetModelNum();

                              bool done_bond_udd_handle = false; // set only for bonds to hydrogen
                              // we use this to set the residue in a bond so that
                              // the bond representation for Hamish contained a residue index
                              // (so that picked atoms on the client can know which other
                              // atoms are in the same residue as the picked atom - for highlighting).

                              graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                              if (element_1 != element_2) {

                                 // Bonded to different atom elements.
                                 if (! is_hydrogen(element_1) && ! is_hydrogen(element_2)) {

                                    add_half_bonds(atom_1_pos, atom_2_pos,
                                                   atom_selection_1[contact[i].id1],
                                                   atom_selection_2[contact[i].id2],
                                                   cc,
                                                   imodel,
                                                   atom_index_1, atom_index_2,
                                                   atom_colour_type,
                                                   udd_user_defined_atom_colour_index_handle,
                                                   nullptr, false, false);
                                 } else {

                                    // Bonds to hydrogens are one colour - HYDROGEN_GREY_BOND, not
                                    // half-bonds.
                                    //
                                    // Except hydrogens on waters are treated differently to other
                                    // hydrogens (if they are not then we don't get to see the
                                    // oxygen).

                                    std::string resname_1 = atom_p_1->GetResName();
                                    std::string resname_2 = atom_p_2->GetResName();
                                    if (resname_1 == "HOH" || resname_2 == "HOH" ||
                                        resname_1 == "DOD" || resname_2 == "DOD") {
                                       add_half_bonds(atom_1_pos, atom_2_pos,
                                                      atom_selection_1[contact[i].id1],
                                                      atom_selection_2[contact[i].id2],
                                                      cc,
                                                      imodel,
                                                      atom_index_1, atom_index_2,
                                                      atom_colour_type,
                                                      udd_user_defined_atom_colour_index_handle,
                                                      nullptr, false, false);
                                    } else {

                                       bool done_h_bond = false; // set when we make a half-bond between H and O.
                                       if (element_1 == " O") {
                                          int bond_udd = graphical_bonds_container::NO_BOND;
                                          atom_p_1->GetUDData(udd_done_bond_handle, bond_udd);
                                          if (bond_udd == graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND) {
                                          } else {
                                             add_half_bonds(atom_1_pos, atom_2_pos,
                                                            atom_selection_1[contact[i].id1],
                                                            atom_selection_2[contact[i].id2],
                                                            cc,
                                                            imodel,
                                                            atom_index_1, atom_index_2,
                                                            atom_colour_type,
                                                            udd_user_defined_atom_colour_index_handle,
                                                            nullptr, false, false);
                                             done_h_bond = true;
                                          }
                                       }

                                       if (element_2 == " O") {
                                          int bond_udd = graphical_bonds_container::NO_BOND;
                                          atom_p_2->GetUDData(udd_done_bond_handle, bond_udd);
                                          if (bond_udd == graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND) {
                                          } else {
                                             add_half_bonds(atom_1_pos, atom_2_pos,
                                                            atom_selection_1[contact[i].id1],
                                                            atom_selection_2[contact[i].id2],
                                                            cc,
                                                            imodel,
                                                            atom_index_1, atom_index_2,
                                                            atom_colour_type,
                                                            udd_user_defined_atom_colour_index_handle,
                                                            nullptr, false, false);
                                             done_h_bond = true;
                                          }
                                       }

                                       if (! done_h_bond) {
                                          if (atom_colour_type != coot::COLOUR_BY_USER_DEFINED_COLOURS) {
                                             graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                             int H_col = HYDROGEN_GREY_BOND;
                                             if (is_deuterium(element_2)) H_col = DEUTERIUM_PINK;
                                             addBond(H_col, atom_1_pos, atom_2_pos, cc, imodel,
                                                     atom_index_1, atom_index_2, false, false);
                                          } else {
                                             add_half_bonds(atom_1_pos, atom_2_pos,
                                                            atom_selection_1[contact[i].id1],
                                                            atom_selection_2[contact[i].id2],
                                                            cc,
                                                            imodel,
                                                            atom_index_1, atom_index_2,
                                                            atom_colour_type,
                                                            udd_user_defined_atom_colour_index_handle,
                                                            nullptr, false, false);
                                          }
                                       }
                                       done_bond_udd_handle = true;
                                       atom_p_1->PutUDData(udd_done_bond_handle, graphical_bonds_container::BONDED_WITH_BOND_TO_HYDROGEN);
                                       atom_p_2->PutUDData(udd_done_bond_handle, graphical_bonds_container::BONDED_WITH_BOND_TO_HYDROGEN);
                                    }
                                 } // not hydrogen test

                              } else {

                                 // Bonded to an atom of the same element.
                                 //

                                 // std::cout << "Bonding here with "
                                 // << coot::atom_spec_t(atom_selection_1[atom_index_1]) << " and "
                                 // << coot::atom_spec_t(atom_selection_2[atom_index_2]) << std::endl;

                                 if (is_hydrogen(element_1)) { // both are hydrogen
                                    float len2 = (atom_1_pos - atom_2_pos).amplitude_squared();
                                    if (len2 < 1.3) { // protection for weirdness, // was 1.0
                                       col = atom_colour(atom_p_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, nullptr);
                                       graphics_line_t::cylinder_class_t cc_bond = graphics_line_t::SINGLE;
                                       addBond(col, atom_1_pos, atom_2_pos, cc_bond, imodel, atom_index_1, atom_index_2, false, false);
                                    }
                                 } else {
                                    // should this test be here or further up?
                                    // Don't bond water Oxygens to each other... 20210817-PE - or anything else.
                                    bool do_it = true;
                                    if (atom_p_1->residue != atom_p_2->residue) {
                                       std::string res_name_1(atom_p_1->residue->GetResName());
                                       if (res_name_1 == "HOH")
                                          do_it = false;
                                       std::string res_name_2(atom_p_2->residue->GetResName());
                                       if (res_name_2 == "HOH")
                                          do_it = false;
                                    }

                                    if (do_it) {
                                       col = atom_colour(atom_p_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, nullptr);
                                       graphics_line_t::cylinder_class_t cc_bond = graphics_line_t::SINGLE;
                                       addBond(col, atom_1_pos, atom_2_pos, cc_bond, imodel, atom_index_1, atom_index_2, false, false);
                                    }
                                 }
                              }

                              // bools and indices
                              mark_atoms_as_bonded(atom_p_1, atom_p_2, have_udd_atoms, udd_done_bond_handle, done_bond_udd_handle);

                           }
                        }
                     }
                  }
               }
            } // contact atom is higher up the list check.
            //          else {
            //             std::cout << "debug:: ignoring contact " << i << std::endl;
            //          }

         } // i over ncontacts

         delete [] contact;

         // OK, now we can handle the het_residues: But we don't want to
         // do this every time that this function is called (X-X, X-H).
         // So do it only on X-X.
         //
         // het_residues is filled for by X-X for everything except HOHs.
         //
         if (! are_different_atom_selections) {
            add_bonds_het_residues(het_residues, asc, imol, atom_colour_type, have_udd_atoms, udd_done_bond_handle,
                                   udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
         }
         if (hoh_residues.size())
            add_bonds_het_residues(hoh_residues, asc, imol, atom_colour_type, have_udd_atoms, udd_done_bond_handle,
                                   udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
      }
   }
}


// add to het_residues maybe.
//
// to be used in conjunction with
// add_bonds_het_residues(het_residues, atom_colour_type, have_udd_atoms, udd_handle)
//
bool
Bond_lines_container::add_bond_by_dictionary_maybe(int imol,
                                                   mmdb::Atom *atom_p_1,
                                                   mmdb::Atom *atom_p_2,
                                                   std::vector<std::pair<bool, mmdb::Residue *> > *het_residues) {

   bool bond_het_residue_by_dictionary = false;
   if (have_dictionary && geom)
      if (atom_p_1->residue == atom_p_2->residue)
         if (atom_p_1->Het)
            if (atom_p_2->Het) {

               // Have we checked this residue type before and failed to find
               // a dictionary for it?  If so, add it to the vector.

               std::pair<bool, mmdb::Residue *> tp0(0, atom_p_1->residue);
               std::pair<bool, mmdb::Residue *> tp1(1, atom_p_1->residue);

               // add this residue to the vector if it is not there already)
               //
               // We have to check both pairs against
               // the cached results (where we have a
               // dictionary and where we don't).
               //
               std::vector<std::pair<bool, mmdb::Residue *> >::const_iterator it_1 =
                  std::find(het_residues->begin(), het_residues->end(), tp0);

               if (it_1 == het_residues->end()) {

                  std::vector<std::pair<bool, mmdb::Residue *> >::const_iterator it_2 =
                     std::find(het_residues->begin(), het_residues->end(), tp1);

                  if (it_2 == het_residues->end()) {

                     if (geom->have_at_least_minimal_dictionary_for_residue_type(atom_p_1->residue->GetResName(), imol)) {

                        if (geom->atoms_match_dictionary(imol, atom_p_1->residue, true, true).first) {

                           het_residues->push_back(tp1);
                           bond_het_residue_by_dictionary = true;
                        } else {
                           het_residues->push_back(tp0);
                        }
                     }  else {
                        het_residues->push_back(tp0);
                     }
                  } else {
                     // this HET group is already in the list and was maked as found in the dictionary.
                     bond_het_residue_by_dictionary = true;
                  }
               } else {
                  // this HET group is already in the list but not found
                  bond_het_residue_by_dictionary = false; // false anyway, I think.
               }
            }

   return bond_het_residue_by_dictionary;
}


void
Bond_lines_container::mark_atoms_as_bonded(mmdb::Atom *atom_p_1, mmdb::Atom *atom_p_2,
                                           bool have_udd_atoms,
                                           int udd_handle,
                                           bool done_bond_udd_handle) const {

   // mark atoms as bonded.
   //
   if (have_udd_atoms) {
      if (! done_bond_udd_handle) { // already happened for bonds to Hs.
         if (! ((!strcmp(atom_p_1->element, " S")) ||
                (!strcmp(atom_p_1->element, "SE")) ||
                (!strcmp(atom_p_1->element, "CL")) ||
                (!strcmp(atom_p_1->element, "BR")) ||
                (!strcmp(atom_p_1->element, "Cl")) ||
                (!strcmp(atom_p_1->element, "Br")) ||
                (!strcmp(atom_p_1->element, " P")))) {
            atom_p_1->PutUDData(udd_handle, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
         }

         if (! ((!strcmp(atom_p_2->element, " S")) ||
                (!strcmp(atom_p_2->element, "SE")) ||
                (!strcmp(atom_p_2->element, "CL")) ||
                (!strcmp(atom_p_2->element, "BR")) ||
                (!strcmp(atom_p_2->element, "Cl")) ||
                (!strcmp(atom_p_2->element, "Br")) ||
                (!strcmp(atom_p_2->element, " P")))) {
            atom_p_2->PutUDData(udd_handle, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
         }
      }
   }
}



void
Bond_lines_container::add_half_bonds(const coot::Cartesian &atom_1_pos,
                                     const coot::Cartesian &atom_2_pos,
                                     mmdb::Atom *at_1,
                                     mmdb::Atom *at_2,
                                     graphics_line_t::cylinder_class_t cc,
                                     int model_number,
                                     int atom_index_1,
                                     int atom_index_2,
                                     int atom_colour_type,
                                     int udd_user_defined_atom_colour_index_handle,
                                     coot::my_atom_colour_map_t *atom_colour_map_p,
                                     bool add_begin_end_cap,
                                     bool add_end_end_cap) {


   // I am pretty sure that half-bonds are single bonds by nature
   // (dashed bonds are not double by nature)
   // graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

   mmdb::Residue *residue_p_1 = at_1->residue;
   mmdb::Residue *residue_p_2 = at_2->residue;

   // is this slow? if so, pass it.
   // int udd_user_defined_atom_colour_index_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   coot::Cartesian bond_mid_point = atom_1_pos.mid_point(atom_2_pos);
   int col_1 = atom_colour(at_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
   addBond(col_1, atom_1_pos, bond_mid_point, cc, model_number, atom_index_1, atom_index_2, add_begin_end_cap, false);

   int col_2 = atom_colour(at_2, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
   addBond(col_2, bond_mid_point, atom_2_pos, cc, model_number, atom_index_1, atom_index_2, false, add_end_end_cap);

   if (false)
      std::cout << "add_half_bonds() with colour indices " << col_1 << " and " << col_2 << std::endl;

}

#include "mini-mol/atom-quads.hh"


void
Bond_lines_container::draw_bonded_quad_atoms_rings(const std::vector<bonded_quad_atoms> &ring_quads,
                                                   int imodel, int atom_colour_type,
                                                   coot::my_atom_colour_map_t *atom_colour_map_p,
                                                   int udd_atom_index_handle,
                                                   int udd_user_defined_atom_colour_index_handle) {

   for (unsigned int i=0; i<ring_quads.size(); i++) {
      const bonded_quad_atoms &bq = ring_quads[i];
      int col = atom_colour(bq.atom_2, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      coot::Cartesian p2(bq.atom_2->x, bq.atom_2->y, bq.atom_2->z);
      coot::Cartesian p3(bq.atom_3->x, bq.atom_3->y, bq.atom_3->z);
      int atom_2_index = -1;
      int atom_3_index = -1;
      bq.atom_2->GetUDData(udd_atom_index_handle, atom_2_index);
      bq.atom_3->GetUDData(udd_atom_index_handle, atom_3_index);

      mmdb::Atom *at_2 = bq.atom_2;
      mmdb::Atom *at_3 = bq.atom_3;
      std::string ele_2(at_2->element);
      std::string ele_3(at_3->element);

      // std::cout << "ring quad " << i << " " << bq << " has bond_type " << bq.bond_type << std::endl;

      if (bq.bond_type == bonded_quad_atoms::SINGLE || bq.bond_type == bonded_quad_atoms::DOUBLE) {

         // draw the single bond.

         // I put this code here because of the mix up with indexing.
         // this is cleaner - without a complete rewrite.

         if (ele_2 == ele_3) {
            int col = atom_colour(at_2, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
            if (bq.bond_type == bonded_quad_atoms::SINGLE)
               addBond(col, p2, p3, cc, imodel, atom_2_index, atom_3_index, false, false);
            if (bq.bond_type == bonded_quad_atoms::DOUBLE)
               addBond(col, p2, p3, cc, imodel, atom_2_index, atom_3_index, true, true);
         } else {
            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
            add_half_bonds(p2, p3, at_2, at_3, cc, imodel, atom_2_index, atom_3_index,
                           atom_colour_type, udd_user_defined_atom_colour_index_handle,
                           atom_colour_map_p, false, false);
         }
      }
      if (bq.bond_type == bonded_quad_atoms::DOUBLE) {
         // nice and consistent...
         coot::Cartesian pt_0(bq.atom_1->x, bq.atom_1->y, bq.atom_1->z);
         coot::Cartesian pt_1(bq.atom_2->x, bq.atom_2->y, bq.atom_2->z);
         coot::Cartesian pt_2(bq.atom_3->x, bq.atom_3->y, bq.atom_3->z);
         coot::Cartesian pt_3(bq.atom_4->x, bq.atom_4->y, bq.atom_4->z);
         coot::Cartesian mp = pt_0.mid_point(pt_3); // doesn't work for cyclopropane

         coot::Cartesian v1 = pt_1 - mp;
         coot::Cartesian v2 = pt_2 - mp;
         float displacement_scale = 0.78; // when bonds are fat this needs to be bigger.
         coot::Cartesian ip1 = mp + v1 * displacement_scale;
         coot::Cartesian ip2 = mp + v2 * displacement_scale;

         if (ele_2 == ele_3) {
            // std::cout << "phenyl ring bond between " << atom_1_index << " " << atom_2_index << std::endl;
            int col = atom_colour(at_2, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
            graphics_line_t::cylinder_class_t cc = graphics_line_t::KEK_DOUBLE_BOND_INNER_BOND;
            addBond(col, ip1, ip2, cc, imodel, atom_2_index, atom_3_index, true, true);
         } else {
            // a nucleotide base is using this function.
            //
            // I need to be able to say that these bonds should have end caps.
            graphics_line_t::cylinder_class_t cc = graphics_line_t::KEK_DOUBLE_BOND_INNER_BOND;
            add_half_bonds(ip1, ip2, at_2, at_3, cc, imodel, atom_2_index, atom_3_index,
                           atom_colour_type,
                           udd_user_defined_atom_colour_index_handle,
                           atom_colour_map_p, true, true);
         }

      }
   }

}

void
Bond_lines_container::draw_trp_rings(const std::vector<mmdb::Atom *> &ring_atoms, int imodel,
                                     int atom_colour_type,
                                     coot::my_atom_colour_map_t *atom_colour_map_p,
                                     int udd_atom_index_handle,
                                     int udd_user_defined_atom_colour_index_handle) {

   // the atoms in ring_atoms are in this order:
   // std::vector<std::string> trp_rings_atom_names = {" CG ", " CD1", " NE1", " CE2", " CD2", " CE3", " CZ3", " CH2", " CZ2"};

   if (ring_atoms.size() != 9) return;

   // single bonds
   std::vector<std::pair<int, int> > vp_single;
   vp_single.push_back(std::pair<int, int> (0, 1));
   vp_single.push_back(std::pair<int, int> (1, 2));
   vp_single.push_back(std::pair<int, int> (2, 3));
   vp_single.push_back(std::pair<int, int> (3, 4));
   vp_single.push_back(std::pair<int, int> (4, 0));

   vp_single.push_back(std::pair<int, int> (3, 8));
   vp_single.push_back(std::pair<int, int> (8, 7));
   vp_single.push_back(std::pair<int, int> (7, 6));
   vp_single.push_back(std::pair<int, int> (6, 5));
   vp_single.push_back(std::pair<int, int> (5, 4));

   // inner doubles
   std::vector<coot::atom_index_quad> inner_doubles;
   inner_doubles.push_back(coot::atom_index_quad(4,0,1,2));
   inner_doubles.push_back(coot::atom_index_quad(8,3,4,5));
   inner_doubles.push_back(coot::atom_index_quad(7,6,5,4));
   inner_doubles.push_back(coot::atom_index_quad(3,8,7,6));

   for (unsigned int i=0; i<vp_single.size(); i++) {
      int iat = vp_single[i].first;
      int jat = vp_single[i].second;
      mmdb::Atom *at_1 = ring_atoms[iat];
      mmdb::Atom *at_2 = ring_atoms[jat];
      int col = atom_colour(at_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      int atom_1_index = -1;
      int atom_2_index = -1;
      coot::Cartesian p1(ring_atoms[iat]->x, ring_atoms[iat]->y, ring_atoms[iat]->z);
      coot::Cartesian p2(ring_atoms[jat]->x, ring_atoms[jat]->y, ring_atoms[jat]->z);
      std::string ele_1(at_1->element);
      std::string ele_2(at_2->element);
      if (ele_1 == ele_2) {
         ring_atoms[iat]->GetUDData(udd_atom_index_handle, atom_1_index);
         ring_atoms[jat]->GetUDData(udd_atom_index_handle, atom_2_index);
         graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
         addBond(col, p1, p2, cc, imodel, atom_1_index, atom_2_index, false, false);
      } else {
         ring_atoms[iat]->GetUDData(udd_atom_index_handle, atom_1_index);
         ring_atoms[jat]->GetUDData(udd_atom_index_handle, atom_2_index);
         graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
         add_half_bonds(p1, p2, at_1, at_2, cc, imodel, atom_1_index, atom_2_index,
                        atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p, false, false);
      }
   }

   for (unsigned int i=0; i<inner_doubles.size(); i++) {

      int iat_1 = inner_doubles[i].index2;
      int iat_2 = inner_doubles[i].index3;
      int iat_3 = inner_doubles[i].index4;
      int iat_0 = inner_doubles[i].index1;

      // find the mid point of atom 0 and 3. the innner bond ends will be on the vector from there to
      // atoms 1 and 2.
      coot::Cartesian pt_0(ring_atoms[iat_0]->x, ring_atoms[iat_0]->y, ring_atoms[iat_0]->z);
      coot::Cartesian pt_1(ring_atoms[iat_1]->x, ring_atoms[iat_1]->y, ring_atoms[iat_1]->z);
      coot::Cartesian pt_2(ring_atoms[iat_2]->x, ring_atoms[iat_2]->y, ring_atoms[iat_2]->z);
      coot::Cartesian pt_3(ring_atoms[iat_3]->x, ring_atoms[iat_3]->y, ring_atoms[iat_3]->z);
      coot::Cartesian mp = pt_0.mid_point(pt_3);

      coot::Cartesian v1 = pt_1 - mp;
      coot::Cartesian v2 = pt_2 - mp;
      coot::Cartesian ip1 = mp + v1 * 0.8;
      coot::Cartesian ip2 = mp + v2 * 0.8;
      mmdb::Atom *at_1 = ring_atoms[iat_1];
      mmdb::Atom *at_2 = ring_atoms[iat_2];
      int col = atom_colour(at_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      int atom_1_index = -1;
      int atom_2_index = -1;
      ring_atoms[iat_1]->GetUDData(udd_atom_index_handle, atom_1_index);
      ring_atoms[iat_2]->GetUDData(udd_atom_index_handle, atom_2_index);
      std::string ele_1(at_1->element);
      std::string ele_2(at_2->element);
      graphics_line_t::cylinder_class_t cc = graphics_line_t::KEK_DOUBLE_BOND_INNER_BOND;
      if (ele_1 == ele_2) {
         bool add_end_cap = true;
         // std::cout << "in innner_doubles " << i << " " << atom_1_index << " " << atom_2_index << " " << add_end_cap<< std::endl;
         addBond(col, ip1, ip2, cc, imodel, atom_1_index, atom_2_index, add_end_cap, add_end_cap);
      } else {
         bool add_end_cap = true;
         add_half_bonds(ip1, ip2, at_1, at_2, cc, imodel, atom_1_index, atom_2_index, atom_colour_type,
                         udd_user_defined_atom_colour_index_handle, atom_colour_map_p, add_end_cap, add_end_cap);
      }
   }

}

void
Bond_lines_container::draw_GA_rings(const std::vector<mmdb::Atom *> &ring_atoms, int imodel,
                                    int atom_colour_type,
                                    coot::my_atom_colour_map_t *atom_colour_map_p,
                                    int udd_atom_index_handle,
                                    int udd_user_defined_atom_colour_index_handle) {

   // the atoms in ring_atoms are in this order:
   // std::vector<std::string> rings_atom_names = {" CG ", " CD1", " NE1", " CE2", " CD2", " CE3", " CZ3", " CH2", " CZ2"};

   if (ring_atoms.size() != 9) return;

   // I can't call draw_GA_rings_outer() here - I don't have the right parameters/arguments.

   // debug
   if (false)
      for (unsigned int i=0; i<ring_atoms.size(); i++)
         std::cout << "  " << i << " " << coot::atom_spec_t(ring_atoms[i]) << std::endl;

   std::string rt = ring_atoms[0]->residue->GetResName();

   // single bonds
   std::vector<std::pair<int, int> > vp_single;
   vp_single.push_back(std::pair<int, int> (0, 1));
   vp_single.push_back(std::pair<int, int> (1, 2));
   vp_single.push_back(std::pair<int, int> (2, 3));
   vp_single.push_back(std::pair<int, int> (3, 4));
   vp_single.push_back(std::pair<int, int> (4, 0));

   vp_single.push_back(std::pair<int, int> (3, 8));
   vp_single.push_back(std::pair<int, int> (8, 7));
   vp_single.push_back(std::pair<int, int> (7, 6));
   vp_single.push_back(std::pair<int, int> (6, 5));
   vp_single.push_back(std::pair<int, int> (5, 4));

   // inner doubles
   std::vector<coot::atom_index_quad> inner_doubles;
   inner_doubles.push_back(coot::atom_index_quad(0,1,2,3));
   inner_doubles.push_back(coot::atom_index_quad(8,3,4,5));
   inner_doubles.push_back(coot::atom_index_quad(7,6,5,4));
   // inner_doubles.push_back(coot::atom_index_quad(3,8,7,6)); // not in a G, because O6 bond is double

   if (rt == "A" || rt == "DA")
      inner_doubles.push_back(coot::atom_index_quad(3,8,7,6));

   for (unsigned int i=0; i<vp_single.size(); i++) {
      int iat = vp_single[i].first;
      int jat = vp_single[i].second;
      mmdb::Atom *at_1 = ring_atoms[iat];
      mmdb::Atom *at_2 = ring_atoms[jat];
      int col = atom_colour(at_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      int atom_1_index = -1;
      int atom_2_index = -1;
      coot::Cartesian p1(ring_atoms[iat]->x, ring_atoms[iat]->y, ring_atoms[iat]->z);
      coot::Cartesian p2(ring_atoms[jat]->x, ring_atoms[jat]->y, ring_atoms[jat]->z);
      std::string ele_1(at_1->element);
      std::string ele_2(at_2->element);
      graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
      if (ele_1 == ele_2) {
         ring_atoms[iat]->GetUDData(udd_atom_index_handle, atom_1_index);
         ring_atoms[jat]->GetUDData(udd_atom_index_handle, atom_2_index);
         addBond(col, p1, p2, cc, imodel, atom_1_index, atom_2_index);
      } else {
         bool add_end_cap = true;
         add_half_bonds(p1, p2, at_1, at_2, cc, imodel, atom_1_index, atom_2_index, atom_colour_type,
                        udd_user_defined_atom_colour_index_handle,
                        atom_colour_map_p, add_end_cap, add_end_cap);
      }
   }

   for (unsigned int i=0; i<inner_doubles.size(); i++) {
      int iat_1 = inner_doubles[i].index2;
      int iat_2 = inner_doubles[i].index3;
      int iat_3 = inner_doubles[i].index4;
      int iat_0 = inner_doubles[i].index1;

      // find the mid point of atom 0 and 3. the innner bond ends will be on the vector from there to
      // atoms 1 and 2.
      coot::Cartesian pt_0(ring_atoms[iat_0]->x, ring_atoms[iat_0]->y, ring_atoms[iat_0]->z);
      coot::Cartesian pt_1(ring_atoms[iat_1]->x, ring_atoms[iat_1]->y, ring_atoms[iat_1]->z);
      coot::Cartesian pt_2(ring_atoms[iat_2]->x, ring_atoms[iat_2]->y, ring_atoms[iat_2]->z);
      coot::Cartesian pt_3(ring_atoms[iat_3]->x, ring_atoms[iat_3]->y, ring_atoms[iat_3]->z);
      coot::Cartesian mp = pt_0.mid_point(pt_3);

      coot::Cartesian v1 = pt_1 - mp;
      coot::Cartesian v2 = pt_2 - mp;
      coot::Cartesian ip1 = mp + v1 * 0.8;
      coot::Cartesian ip2 = mp + v2 * 0.8;
      mmdb::Atom *at_1 = ring_atoms[iat_1];
      mmdb::Atom *at_2 = ring_atoms[iat_2];
      int col = atom_colour(at_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      int atom_1_index = -1;
      int atom_2_index = -1;
      ring_atoms[iat_1]->GetUDData(udd_atom_index_handle, atom_1_index);
      ring_atoms[iat_2]->GetUDData(udd_atom_index_handle, atom_2_index);
      std::string ele_1(at_1->element);
      std::string ele_2(at_2->element);
      graphics_line_t::cylinder_class_t cc = graphics_line_t::KEK_DOUBLE_BOND_INNER_BOND;
      if (ele_1 == ele_2) {
         addBond(col, ip1, ip2, cc, imodel, atom_1_index, atom_2_index, true, true);
      } else {
         add_half_bonds(ip1, ip2, at_1, at_2, cc, imodel, atom_1_index, atom_2_index,
                        atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p, true, true);
      }
   }

}

// the atoms have been added in order 0 is bonded to 1, 1 is bonded to 2, 2 is bonded to 3 etc.
// and there is a double bond between 0 and 1, 2 and 3, and 4 to 5. Or maybe we could explicitly
// add that to the the ring_atoms data.
//
void
Bond_lines_container::draw_6_membered_ring(const std::string &residue_name,
                                           const std::vector<mmdb::Atom *> &ring_atoms, int imodel,
                                           int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                                           int udd_atom_index_handle,
                                           int udd_user_defined_atom_colour_index_handle) {

   // We will make a representation with single bonds around the outside and double bonds

   if (ring_atoms.size() != 6) return;
   for (unsigned int iat=0; iat<ring_atoms.size(); iat++) {
      int jat = iat + 1;
      if (iat == 5) jat = 0;
      mmdb::Atom *at_1 = ring_atoms[iat];
      mmdb::Atom *at_2 = ring_atoms[jat];
      std::string ele_1(at_1->element);
      std::string ele_2(at_2->element);
      coot::Cartesian p1(ring_atoms[iat]->x, ring_atoms[iat]->y, ring_atoms[iat]->z);
      coot::Cartesian p2(ring_atoms[jat]->x, ring_atoms[jat]->y, ring_atoms[jat]->z);
      int atom_1_index = -1;
      int atom_2_index = -1;
      ring_atoms[iat]->GetUDData(udd_atom_index_handle, atom_1_index);
      ring_atoms[jat]->GetUDData(udd_atom_index_handle, atom_2_index);
      graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
      if (ele_1 == ele_2) {
         int col = atom_colour(at_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
         addBond(col, p1, p2, cc, imodel, atom_1_index, atom_2_index);
      } else {
         // a nucleotide base is using this function
         add_half_bonds(p1, p2, at_1, at_2, cc, imodel, atom_1_index, atom_2_index,
                        atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p, true, true);
      }
   }

   std::vector<coot::atom_index_quad> inner_doubles;

   if (residue_name == "T" || residue_name == "DT" || residue_name == "U") {
      inner_doubles.push_back(coot::atom_index_quad(3,4,5,0));
   } else {
      if (residue_name == "C" || residue_name == "DC") {
         inner_doubles.push_back(coot::atom_index_quad(1,2,3,4));
         inner_doubles.push_back(coot::atom_index_quad(3,4,5,0));
      } else {
         // inside bond between 1 and 2, use 0 and 3 to find the "inside" of the ring
         // old code for i=0, i<3; i++
         // int iat_1 = 2*i;
         // int iat_2 = iat_1+1;
         // int iat_3 = iat_1+2;
         // int iat_0 = iat_1-1;
         // if (iat_0  < 0) iat_0 = 5;
         // if (iat_3 == 6) iat_3 = 0;

         inner_doubles.push_back(coot::atom_index_quad(0,1,2,3));
         inner_doubles.push_back(coot::atom_index_quad(2,3,4,5));
         inner_doubles.push_back(coot::atom_index_quad(4,5,0,1));
      }
   }

   for (unsigned int i=0; i<inner_doubles.size(); i++) {

      int iat_0 = inner_doubles[i].index1;
      int iat_1 = inner_doubles[i].index2;
      int iat_2 = inner_doubles[i].index3;
      int iat_3 = inner_doubles[i].index4;

      // find the mid point of atom 0 and 3. the innner bond ends will be on the vector from thre to
      // atoms 1 and 2.
      coot::Cartesian pt_0(ring_atoms[iat_0]->x, ring_atoms[iat_0]->y, ring_atoms[iat_0]->z);
      coot::Cartesian pt_1(ring_atoms[iat_1]->x, ring_atoms[iat_1]->y, ring_atoms[iat_1]->z);
      coot::Cartesian pt_2(ring_atoms[iat_2]->x, ring_atoms[iat_2]->y, ring_atoms[iat_2]->z);
      coot::Cartesian pt_3(ring_atoms[iat_3]->x, ring_atoms[iat_3]->y, ring_atoms[iat_3]->z);
      coot::Cartesian mp = pt_0.mid_point(pt_3);

      coot::Cartesian v1 = pt_1 - mp;
      coot::Cartesian v2 = pt_2 - mp;
      coot::Cartesian ip1 = mp + v1 * 0.78;
      coot::Cartesian ip2 = mp + v2 * 0.78;
      mmdb::Atom *at_1 = ring_atoms[iat_1];
      mmdb::Atom *at_2 = ring_atoms[iat_2];
      std::string ele_1(at_1->element);
      std::string ele_2(at_2->element);
      int atom_1_index = -1;
      int atom_2_index = -1;
      ring_atoms[iat_1]->GetUDData(udd_atom_index_handle, atom_1_index);
      ring_atoms[iat_2]->GetUDData(udd_atom_index_handle, atom_2_index);
      graphics_line_t::cylinder_class_t cc = graphics_line_t::KEK_DOUBLE_BOND_INNER_BOND;
      if (ele_1 == ele_2) {
         // std::cout << "phenyl ring bond between " << atom_1_index << " " << atom_2_index << std::endl;
         int col = atom_colour(at_1, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
         addBond(col, ip1, ip2, cc, imodel, atom_1_index, atom_2_index, true, true);
      } else {
         // a nucleotide base is using this function
         add_half_bonds(ip1, ip2, at_1, at_2, cc, imodel, atom_1_index, atom_2_index,
                        atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p, true, true);
      }
   }

}

// is_deloc is an optional arg (default 0).
//
// The atom indices iat_1 and iat_2 are residue-atom indices, not all-molecule atom indices.
//
void
Bond_lines_container::add_double_bond(int imol, int imodel, int iat_1, int iat_2,
                                      mmdb::PPAtom residue_atoms, int n_residue_atoms, // atoms and n_atoms for the residue
                                      int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                                      int udd_atom_index_handle,
                                      int udd_user_defined_atom_colour_index_handle,
                                      const std::vector<coot::dict_bond_restraint_t> &bond_restraints,
                                      bool is_deloc) {

   std::string ele_1 = residue_atoms[iat_1]->element;
   std::string ele_2 = residue_atoms[iat_2]->element;

   graphics_line_t::cylinder_class_t cc = graphics_line_t::DOUBLE;

   int idx_1_mol_indexing = -1;
   int idx_2_mol_indexing = -1;

   residue_atoms[iat_1]->GetUDData(udd_atom_index_handle, idx_1_mol_indexing);
   residue_atoms[iat_2]->GetUDData(udd_atom_index_handle, idx_2_mol_indexing);

   try {

      // perp_n is the direction of the offset (from the atom position) of the start and
      // finish points in the plane of the double bond.
      //
      clipper::Coord_orth pos_at_1(residue_atoms[iat_1]->x, residue_atoms[iat_1]->y, residue_atoms[iat_1]->z);
      clipper::Coord_orth pos_at_2(residue_atoms[iat_2]->x, residue_atoms[iat_2]->y, residue_atoms[iat_2]->z);
      clipper::Coord_orth n_n = get_neighb_normal(imol, iat_1, iat_2, residue_atoms, n_residue_atoms);
      clipper::Coord_orth b(pos_at_1 - pos_at_2);
      clipper::Coord_orth b_n(b.unit());
      clipper::Coord_orth perp_n(clipper::Coord_orth::cross(n_n, b_n));
      // std::cout << "    perp_n " << perp_n.format() << " from " << n_n.format() << " x " << b_n.format() << std::endl;
      if (is_deloc)
         if (invert_deloc_bond_displacement_vector(perp_n, iat_1, iat_2, residue_atoms, n_residue_atoms, bond_restraints))
            perp_n = -perp_n;
      // std::cout << "now perp_n " << perp_n.format() << std::endl;
      int col = atom_colour(residue_atoms[iat_1], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      double offset = 0.13;
      clipper::Coord_orth pt_1_1 = pos_at_1 - offset * perp_n;
      clipper::Coord_orth pt_1_2 = pos_at_1 + offset * perp_n;
      clipper::Coord_orth pt_2_1 = pos_at_2 - offset * perp_n;
      clipper::Coord_orth pt_2_2 = pos_at_2 + offset * perp_n;

      if (ele_1 == ele_2) {
         // simple double bond (e.g. C=C)
         addBond(col, coot::Cartesian(pt_1_1), coot::Cartesian(pt_2_1), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing, true, true);
         if (! is_deloc)
            addBond(col, coot::Cartesian(pt_1_2), coot::Cartesian(pt_2_2), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing, true, true);
         else
            add_dashed_bond(col, coot::Cartesian(pt_1_2), coot::Cartesian(pt_2_2), NOT_HALF_BOND, cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
      } else {

         // we have to draw double half bonds, e.g. C=0
         clipper::Coord_orth bond_mid_point = 0.5 * clipper::Coord_orth(pos_at_1 + pos_at_2);
         clipper::Coord_orth mp_1 = bond_mid_point - offset * perp_n;
         clipper::Coord_orth mp_2 = bond_mid_point + offset * perp_n;
         if (! is_deloc) {

            addBond(col, coot::Cartesian(pt_1_1), coot::Cartesian(mp_1), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing, true, false);
            addBond(col, coot::Cartesian(pt_1_2), coot::Cartesian(mp_2), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing, true, false);
            col = atom_colour(residue_atoms[iat_2], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
            addBond(col, coot::Cartesian(pt_2_1), coot::Cartesian(mp_1), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing, true, false);
            addBond(col, coot::Cartesian(pt_2_2), coot::Cartesian(mp_2), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing, true, false);
         } else {
            addBond(col, coot::Cartesian(pt_1_1), coot::Cartesian(mp_1), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing, true, false);
            // add dashed bond doesn't take a residue pointer argument (yet)
            add_dashed_bond(col, coot::Cartesian(pt_1_2), coot::Cartesian(mp_2), HALF_BOND_FIRST_ATOM, cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
            col = atom_colour(residue_atoms[iat_2], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
            addBond(col, coot::Cartesian(pt_2_1), coot::Cartesian(mp_1), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing, true, false);
            add_dashed_bond(col, coot::Cartesian(pt_2_2), coot::Cartesian(mp_2), HALF_BOND_SECOND_ATOM, cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
         }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "caught exception add_double_bond(): " << rte.what() << std::endl;
   }
}


// these are residue atoms and residue n_atoms
void
Bond_lines_container::add_triple_bond(int imol, int imodel, int iat_1, int iat_2, mmdb::PPAtom atoms, int n_atoms,
                                      int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                                      int udd_atom_index_handle,
                                      int udd_user_defined_atom_colour_index_handle,
                                      const std::vector<coot::dict_bond_restraint_t> &bond_restraints) {

   graphics_line_t::cylinder_class_t cc = graphics_line_t::TRIPLE;

   //
   std::string ele_1 = atoms[iat_1]->element;
   std::string ele_2 = atoms[iat_2]->element;

   mmdb::Residue *residue_p_1 = atoms[iat_1]->residue;
   mmdb::Residue *residue_p_2 = atoms[iat_2]->residue;

   try {

      // need to look up the all-molecule atom index, not the residue atom index

      int idx_1_mol_indexing = -1;
      int idx_2_mol_indexing = -1;

      atoms[iat_1]->GetUDData(udd_atom_index_handle, idx_1_mol_indexing);
      atoms[iat_2]->GetUDData(udd_atom_index_handle, idx_2_mol_indexing);

      bool also_2nd_order = 1; // because linear nature of bonds to
                               // atoms in triple bond means we need
                               // more atoms.
      clipper::Coord_orth pos_at_1(atoms[iat_1]->x, atoms[iat_1]->y, atoms[iat_1]->z);
      clipper::Coord_orth pos_at_2(atoms[iat_2]->x, atoms[iat_2]->y, atoms[iat_2]->z);
      clipper::Coord_orth n_n = get_neighb_normal(imol, iat_1, iat_2, atoms, n_atoms, also_2nd_order);
      clipper::Coord_orth b(pos_at_1 - pos_at_2);
      clipper::Coord_orth b_n(b.unit());
      clipper::Coord_orth perp_n(clipper::Coord_orth::cross(n_n, b_n));
      int col = atom_colour(atoms[iat_1], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      double offset = 0.08;
      if (for_GL_solid_model_rendering)
         offset = 0.18;
      clipper::Coord_orth pt_1_1 = pos_at_1 - offset * perp_n;
      clipper::Coord_orth pt_1_2 = pos_at_1;
      clipper::Coord_orth pt_1_3 = pos_at_1 + offset * perp_n;
      clipper::Coord_orth pt_2_1 = pos_at_2 - offset * perp_n;
      clipper::Coord_orth pt_2_2 = pos_at_2;
      clipper::Coord_orth pt_2_3 = pos_at_2 + offset * perp_n;
      if (ele_1 == ele_2) {
         // e.g. -C#C-
         addBond(col, coot::Cartesian(pt_1_1), coot::Cartesian(pt_2_1), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
         addBond(col, coot::Cartesian(pt_1_2), coot::Cartesian(pt_2_2), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
         addBond(col, coot::Cartesian(pt_1_3), coot::Cartesian(pt_2_3), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
      } else {
         // e.g. -C#N
         clipper::Coord_orth bond_mid_point = 0.5 * clipper::Coord_orth(pos_at_1 + pos_at_2);
         clipper::Coord_orth mp_1 = bond_mid_point - offset * perp_n;
         clipper::Coord_orth mp_2 = bond_mid_point;
         clipper::Coord_orth mp_3 = bond_mid_point + offset * perp_n;
         addBond(col, coot::Cartesian(pt_1_1), coot::Cartesian(mp_1), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
         addBond(col, coot::Cartesian(pt_1_2), coot::Cartesian(mp_2), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
         addBond(col, coot::Cartesian(pt_1_3), coot::Cartesian(mp_3), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
         col = atom_colour(atoms[iat_2], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
         addBond(col, coot::Cartesian(pt_2_1), coot::Cartesian(mp_1), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
         addBond(col, coot::Cartesian(pt_2_2), coot::Cartesian(mp_2), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
         addBond(col, coot::Cartesian(pt_2_3), coot::Cartesian(mp_3), cc, imodel, idx_1_mol_indexing, idx_2_mol_indexing);
      }
   }

   catch (const std::runtime_error &rte) {
      std::cout << "caught exception add_double_bond(): " << rte.what() << std::endl;
   }
}


// also_2nd_order_neighbs_flag is a optional arg (default 0)
//
clipper::Coord_orth
Bond_lines_container::get_neighb_normal(int imol, int iat_1, int iat_2, mmdb::PPAtom atoms, int n_atoms,
                                        bool also_2nd_order_neighbs_flag) const {

   clipper::Coord_orth pt(0,0,0);
   if (have_dictionary) {
      std::string rn = atoms[iat_1]->residue->GetResName();
      std::string at_n_1 = atoms[iat_1]->name;
      std::string at_n_2 = atoms[iat_2]->name;
      std::vector<std::string> neighbours = geom->get_bonded_neighbours(rn, imol, at_n_1, at_n_2,
                                                                        also_2nd_order_neighbs_flag);

      if (0) {
         std::cout << "======== neighbours of " << at_n_1 << " and " << at_n_2 << ":" << std::endl;
         for (unsigned int i=0; i<neighbours.size(); i++)
            std::cout << "   " << neighbours[i] << std::endl;
      }

      std::string alt_conf_bond = atoms[iat_1]->altLoc; // same as iat_2 by the time we get here, I think
      if (neighbours.size() > 2) {
         std::vector<mmdb::Atom *> neighb_atoms;
         for (unsigned int i=0; i<neighbours.size(); i++) {
            for (int j=0; j<n_atoms; j++) {
               std::string atom_name = atoms[j]->name;
               if (neighbours[i] == atom_name) {
                  std::string alt_conf_atom = atoms[j]->altLoc;
                  if (alt_conf_atom == alt_conf_bond) {
                     // only add them if they are not there already (belt and braces test)
                     if (std::find(neighb_atoms.begin(), neighb_atoms.end(), atoms[j]) ==
                         neighb_atoms.end())
                        neighb_atoms.push_back(atoms[j]);
                  }
               }
            }
         }
         if (neighb_atoms.size() > 2) {

            // std::cout << " :::: found " << neighb_atoms.size() << " neibhbours" << std::endl;

            std::vector<clipper::Coord_orth> neighb_atoms_pos(neighb_atoms.size());
            for (unsigned int i=0; i<neighb_atoms.size(); i++)
               neighb_atoms_pos[i] = clipper::Coord_orth(neighb_atoms[i]->x,
                                                         neighb_atoms[i]->y,
                                                         neighb_atoms[i]->z);
            coot::lsq_plane_info_t lp(neighb_atoms_pos);
            pt = lp.normal();
         }
      } else {
         std::string m = "Not enough atoms to determine orientation of ";
         m += atoms[iat_1]->residue->GetResName();
         m += " - dictionary bonding fails";
         m += " found ";
         m += coot::util::int_to_string(neighbours.size());
         m += " neighbs: ";
         for (unsigned int i=0; i<neighbours.size(); i++) {
            m += neighbours[i];
            m += " ";
         }
         // make something up...
         pt = clipper::Coord_orth(0,0,1);
         // throw(std::runtime_error(m));
      }
   } else {
      // this should not happend
      std::string m = "No dictionary for ";
      m += atoms[iat_1]->residue->GetResName();
      m += " - dictionary bonding fails";
      throw(std::runtime_error(m));
   }
   return pt;

}

bool
Bond_lines_container::invert_deloc_bond_displacement_vector(const clipper::Coord_orth &vect,
                                                            int iat_1, int iat_2, mmdb::PPAtom residue_atoms, int n_atoms,
                                                            const std::vector<coot::dict_bond_restraint_t> &bond_restraints) const {

   bool r = false;

//    std::cout << " ==================== considering the swap of :"
//              << residue_atoms[iat_1]->name << ": to :"
//              << residue_atoms[iat_2]->name << ": =========================" << std::endl;

   std::string atom_name_iat = residue_atoms[iat_1]->name;
   std::string atom_name_jat = residue_atoms[iat_2]->name;

   std::map<std::string, int> atom_name_map;
   for (int iat=0; iat<n_atoms; iat++)
      atom_name_map[residue_atoms[iat]->name] = iat;

   for (unsigned int ib=0; ib<bond_restraints.size(); ib++) {
      if (bond_restraints[ib].atom_id_1_4c() == atom_name_iat) {
         if (bond_restraints[ib].atom_id_2_4c() != atom_name_jat) {

//             std::cout << "::::: bond 1 from :" << bond_restraints[ib].atom_id_1_4c()
//                       << ": to :" << bond_restraints[ib].atom_id_2_4c() << ": type "
//                       << bond_restraints[ib].type()
//                       << std::endl;

            if (bond_restraints[ib].type() == "deloc") {
               clipper::Coord_orth pt_1(residue_atoms[iat_1]->x,
                                        residue_atoms[iat_1]->y,
                                        residue_atoms[iat_1]->z);
               std::map<std::string, int>::const_iterator it;
               it = atom_name_map.find(bond_restraints[ib].atom_id_2_4c());
               if (it != atom_name_map.end()) {
                  clipper::Coord_orth pt_2(residue_atoms[it->second]->x,
                                           residue_atoms[it->second]->y,
                                           residue_atoms[it->second]->z);
                  clipper::Coord_orth diff = pt_2 - pt_1;
                  double d = clipper::Coord_orth::dot(vect, diff);
                  // std::cout << "    dot 1 : " << d << std::endl;
                  if (d < 0)
                     r = true;
                  break;
               }
            }
         }
      }


      // same again, restraints ordered differently.
      //
      if (bond_restraints[ib].atom_id_2_4c() == atom_name_iat) {
         if (bond_restraints[ib].atom_id_1_4c() != atom_name_jat) {

//             std::cout << "::::: bond 1 from :" << bond_restraints[ib].atom_id_1_4c()
//                       << ": to :" << bond_restraints[ib].atom_id_2_4c() << ": type "
//                       << bond_restraints[ib].type()
//                       << std::endl;

            if (bond_restraints[ib].type() == "deloc") {
               clipper::Coord_orth pt_1(residue_atoms[iat_1]->x,
                                        residue_atoms[iat_1]->y,
                                        residue_atoms[iat_1]->z);
               std::map<std::string, int>::const_iterator it;
               it = atom_name_map.find(bond_restraints[ib].atom_id_1_4c());
               if (it != atom_name_map.end()) {
                  clipper::Coord_orth pt_2(residue_atoms[it->second]->x,
                                           residue_atoms[it->second]->y,
                                           residue_atoms[it->second]->z);
                  clipper::Coord_orth diff = pt_2 - pt_1;
                  double d = clipper::Coord_orth::dot(vect, diff);
                  // std::cout << "    dot 2 : " << d << std::endl;
                  if (d < 0)
                     r = true;
                  break;
               }
            }
         }
      }
   }

   // std::cout << ":::::::::: invert_deloc_bond_displacement_vector() returns " << r << std::endl;
   return r;
}


void
Bond_lines_container::add_bonds_het_residues(const std::vector<std::pair<bool, mmdb::Residue *> > &het_residues,
                                             const atom_selection_container_t &atom_sel,
                                             int imol,
                                             int atom_colour_type,
                                             short int have_udd_handle,
                                             int udd_bond_handle,
                                             int udd_atom_index_handle,
                                             int udd_user_defined_atom_colour_index_handle) {

   coot::my_atom_colour_map_t atom_colour_map;
   atom_colour_map.fill_chain_id_map(atom_sel); // 20230206-PE should be the same "algorithm" as is
                                                // used in fill_default_colour_rules()
                                                // so that M2T ribbons have the same colours, by default
                                                // add Coot does in colour-by-chain-and-dictionary mode.

   if (! het_residues.empty()) {
      for (unsigned int ires=0; ires<het_residues.size(); ires++) {
         if (het_residues[ires].first) {
            std::string res_name = het_residues[ires].second->GetResName();
            std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
               geom->get_monomer_restraints_at_least_minimal(res_name, imol);

            if (false) // debugging
               if (res_name != "HOH")
                  std::cout << "          ============== Considering bonding HET residue: "
                            << coot::residue_spec_t(het_residues[ires].second) << " "
                            << res_name << " " << std::endl;

            if (! restraints.first) {
               std::cout << "Oooppps!  No bonding rules for residue type :" << res_name
                         << ": missing bonds! " << std::endl;
            } else {

               mmdb::PPAtom residue_atoms;
               int n_atoms;
               het_residues[ires].second->GetAtomTable(residue_atoms, n_atoms);
               int model_number = het_residues[ires].second->GetModelNum();
               for (unsigned int ib=0; ib<restraints.second.bond_restraint.size(); ib++) {
                  std::string atom_name_1 = restraints.second.bond_restraint[ib].atom_id_1_4c();
                  std::string atom_name_2 = restraints.second.bond_restraint[ib].atom_id_2_4c();
                  std::string bt = restraints.second.bond_restraint[ib].type();

                  bool added_bond = false;
                  for (int iat=0; iat<n_atoms; iat++) {
                     std::string residue_atom_name_1(residue_atoms[iat]->name);
                     if (atom_name_1 == residue_atom_name_1) {
                        for (int jat=0; jat<n_atoms; jat++) {
                           std::string residue_atom_name_2(residue_atoms[jat]->name);
                           if (atom_name_2 == residue_atom_name_2) {
                              std::string aloc_1 = residue_atoms[iat]->altLoc;
                              std::string aloc_2 = residue_atoms[jat]->altLoc;
                              if (aloc_1 == aloc_2 || aloc_1.empty() || aloc_2.empty()) {
                                 coot::Cartesian p1(residue_atoms[iat]->x,
                                                    residue_atoms[iat]->y,
                                                    residue_atoms[iat]->z);
                                 coot::Cartesian p2(residue_atoms[jat]->x,
                                                    residue_atoms[jat]->y,
                                                    residue_atoms[jat]->z);

                                 int iat_1_atom_index = -1;
                                 int iat_2_atom_index = -1;
                                 residue_atoms[iat]->GetUDData(udd_atom_index_handle, iat_1_atom_index);
                                 residue_atoms[jat]->GetUDData(udd_atom_index_handle, iat_2_atom_index);

                                 if (false)
                                    std::cout << "making bond between :" << residue_atoms[iat]->name
                                              << ": and :" << residue_atoms[jat]->name << ": "
                                              << bt << std::endl;

                                 std::string element_1 = residue_atoms[iat]->element;
                                 std::string element_2 = residue_atoms[jat]->element;

                                 if (element_1 != element_2) {

                                    // Bonded to different atom elements.

                                    // if ((element_1 != " H") && (element_2 != " H")) {
                                    if (!is_hydrogen(element_1) && !is_hydrogen(element_2)) {
                                       if (bt == "double") {
                                          add_double_bond(imol, model_number, iat, jat, residue_atoms, n_atoms, atom_colour_type,
                                                          &atom_colour_map, udd_atom_index_handle,
                                                          udd_user_defined_atom_colour_index_handle,
                                                          restraints.second.bond_restraint);
                                       } else {
                                          if (bt == "deloc") {
                                             bool is_deloc = true;
                                             add_double_bond(imol, model_number, iat, jat, residue_atoms, n_atoms,
                                                             atom_colour_type,
                                                             &atom_colour_map,
                                                             udd_atom_index_handle,
                                                             udd_user_defined_atom_colour_index_handle,
                                                             restraints.second.bond_restraint, is_deloc);
                                          } else {

                                             if (bt == "triple") {
                                                add_triple_bond(imol, model_number, iat, jat, residue_atoms, n_atoms,
                                                                atom_colour_type, &atom_colour_map,
                                                                udd_atom_index_handle,
                                                                udd_user_defined_atom_colour_index_handle,
                                                                restraints.second.bond_restraint);
                                             } else {
                                                // could be "metal"
                                                graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                                add_half_bonds(p1, p2,
                                                               residue_atoms[iat],
                                                               residue_atoms[jat],
                                                               cc,
                                                               model_number,
                                                               iat_1_atom_index, iat_2_atom_index,
                                                               atom_colour_type,
                                                               udd_user_defined_atom_colour_index_handle,
                                                               &atom_colour_map,
                                                               false, false);
                                             }
                                          }
                                       }
                                    } else {
                                       if (do_bonds_to_hydrogens) {
                                          if (res_name == "HOH" || res_name == "DOD") {
                                             graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                             add_half_bonds(p1, p2,
                                                            residue_atoms[iat],
                                                            residue_atoms[jat],
                                                            cc,
                                                            model_number,
                                                            iat_1_atom_index, iat_2_atom_index,
                                                            atom_colour_type,
                                                            udd_user_defined_atom_colour_index_handle,
                                                            &atom_colour_map, true, true);
                                          } else {
                                             graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                             addBond(HYDROGEN_GREY_BOND, p1, p2, cc, model_number, iat_1_atom_index, iat_2_atom_index); // 20171224-PE correct indices?w
                                          }
                                       }
                                    }

                                 } else {

                                    // Bonded to an atom of the same element.
                                    // add_double_bond() uses residue-atom indexing.
                                    // addBond uses all-molecule atom indexing.
                                    int col = atom_colour(residue_atoms[iat], atom_colour_type, udd_user_defined_atom_colour_index_handle,
                                                          &atom_colour_map);
                                    if (bt == "double") {
                                       add_double_bond(imol, model_number, iat, jat, residue_atoms, n_atoms,
                                                       atom_colour_type, &atom_colour_map,
                                                       udd_atom_index_handle,
                                                       udd_user_defined_atom_colour_index_handle,
                                                       restraints.second.bond_restraint);
                                    } else {
                                       if (bt == "deloc") {
                                          bool is_deloc = 1;
                                          add_double_bond(imol, model_number, iat, jat, residue_atoms, n_atoms,
                                                          atom_colour_type, &atom_colour_map,
                                                          udd_atom_index_handle,
                                                          udd_user_defined_atom_colour_index_handle,
                                                          restraints.second.bond_restraint, is_deloc);
                                       } else {
                                          if (bt == "triple") {
                                             add_triple_bond(imol, model_number, iat, jat, residue_atoms, n_atoms,
                                                             atom_colour_type, &atom_colour_map,
                                                             udd_atom_index_handle,
                                                             udd_user_defined_atom_colour_index_handle,
                                                             restraints.second.bond_restraint);
                                          } else {
                                             addBond(col, p1, p2, graphics_line_t::SINGLE, model_number, iat_1_atom_index, iat_2_atom_index);
                                          }
                                       }
                                    }
                                 }

                                 if (have_udd_handle) {
                                    residue_atoms[iat]->PutUDData(udd_bond_handle, graphical_bonds_container::BONDED_WITH_HETATM_BOND);
                                    residue_atoms[jat]->PutUDData(udd_bond_handle, graphical_bonds_container::BONDED_WITH_HETATM_BOND);
                                 }
                                 added_bond = true;
                                 // break; // nope! if this is in place: "" - "" is OK
                                 //                                      "A" - "A" is OK
                                 //                                      "A" - "B" is rejected as should be
                                 //                                      "" - "A" is OK
                                 //                             but      "" - "B" is is not considered

                              }
                           }
                        }
                     }
                     // 20110607: Nope!  We want to be able to bond
                     // het_residues with alt confs too.
                     //
                     // if (added_bond)
                     // break;
                  }
               }
            }

            // now aromatic ring systems.
            // int col = 0;
            // het_residue_aromatic_rings(het_residues[ires].second, restraints.second, udd_atom_index_handle, col);
         }
      }
   }
}


std::vector<std::pair<std::string, std::string> >
Bond_lines_container::get_aromatic_bonds(const coot::dictionary_residue_restraints_t &restraints) const {

   std::vector<std::pair<std::string, std::string> > aromatic_bonds;
   for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
      const coot::dict_bond_restraint_t &br = restraints.bond_restraint[ib];
      if (br.type() == "aromatic") { // old CCP4 dictionary style
         const std::string &atom_name_1 = br.atom_id_1_4c();
         const std::string &atom_name_2 = br.atom_id_2_4c();
         std::pair<std::string, std::string> p(atom_name_1, atom_name_2);
         aromatic_bonds.push_back(p);
      }
      if (br.aromaticity == coot::dict_bond_restraint_t::AROMATIC) { // new acedrg style
         const std::string &atom_name_1 = br.atom_id_1_4c();
         const std::string &atom_name_2 = br.atom_id_2_4c();
         std::pair<std::string, std::string> p(atom_name_1, atom_name_2);
         aromatic_bonds.push_back(p);
      }
   }
   return aromatic_bonds;
}

void
Bond_lines_container::het_residue_aromatic_rings(mmdb::Residue *res,
                                                 const coot::dictionary_residue_restraints_t &restraints,
                                                 int udd_atom_index_handle,
                                                 int col) {

   // not used?

   enum representation_style { AROMATIC, KEKULIZED };
   representation_style rs = KEKULIZED;

   std::vector<std::pair<std::string, std::string> > aromatic_bonds;
   for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
      const coot::dict_bond_restraint_t &br = restraints.bond_restraint[ib];
      if (br.type() == "aromatic") { // old CCP4 dictionary style
         std::string atom_name_1 = br.atom_id_1_4c();
         std::string atom_name_2 = br.atom_id_2_4c();
         std::pair<std::string, std::string> p(atom_name_1, atom_name_2);
         aromatic_bonds.push_back(p);
      }
      if (br.aromaticity == coot::dict_bond_restraint_t::AROMATIC) { // new acedrg style
         std::string atom_name_1 = br.atom_id_1_4c();
         std::string atom_name_2 = br.atom_id_2_4c();
         std::pair<std::string, std::string> p(atom_name_1, atom_name_2);
         aromatic_bonds.push_back(p);
      }
   }
   if (aromatic_bonds.size() > 4) {
      coot::aromatic_graph_t ag(aromatic_bonds);
      std::vector<std::vector<std::string> > rings = ag.ring_list();
      if (rs == AROMATIC) {
         for (unsigned int i=0; i<rings.size(); i++) {
            add_aromatic_ring_bond_lines(rings[i], res, udd_atom_index_handle, col);
         }
      }
   }
}

// pass a list of atom name that are part of the aromatic ring system.
void
Bond_lines_container::add_aromatic_ring_bond_lines(const std::vector<std::string> &ring_atom_names,
                                                   mmdb::Residue *residue_p,
                                                   int udd_atom_index_handle, int col) {


   // We can't have aromatic rings with more than 7 atoms (can we?)

   if (ring_atom_names.size() < 8) {

      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;

      std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs(residue_p);

      for (unsigned int i_alt_conf=0; i_alt_conf<alt_confs.size(); i_alt_conf++) {
         std::vector<mmdb::Atom *> found_atoms;

         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (unsigned int i=0; i<ring_atom_names.size(); i++) {
            for (int iat=0; iat<n_residue_atoms; iat++) {
               std::string atom_name(residue_atoms[iat]->name);
               std::string atom_alt_conf(residue_atoms[iat]->altLoc);
               if (atom_alt_conf == alt_confs[i_alt_conf]) {
                  if (atom_name == ring_atom_names[i]) {
                     found_atoms.push_back(residue_atoms[iat]);
                  }
               }
            }
         }

         if (found_atoms.size() == ring_atom_names.size()) {

            bool skip_this_ring = false;
            std::vector<clipper::Coord_orth> pts(ring_atom_names.size());
            for (unsigned int iat=0; iat<found_atoms.size(); iat++) {
               mmdb::Atom *found_atom = found_atoms[iat];
               pts[iat] = clipper::Coord_orth(found_atom->x,
                                              found_atom->y,
                                              found_atom->z);
               int idx_mol = -1;
               found_atom->GetUDData(udd_atom_index_handle, idx_mol);
               if (! skip_this_ring)
                  if (no_bonds_to_these_atoms.find(idx_mol) != no_bonds_to_these_atoms.end())
                     skip_this_ring = true;
            }
            coot::lsq_plane_info_t lp(pts);
            clipper::Coord_orth n = lp.normal();
            clipper::Coord_orth c = lp.centre();
            double radius = 0.8;
            if (ring_atom_names.size() == 5)
               radius = 0.6;

            int n_steps = 40;
            double step_frac = double(1.0/n_steps);

            // we want a point in the lsq plane that is radius A away from
            // centre.
            clipper::Coord_orth arb(0.2, 0.8, 0.1);
            clipper::Coord_orth cr(clipper::Coord_orth::cross(n, arb).unit());
            clipper::Coord_orth first_pt = c + radius * cr;

            if (! for_GL_solid_model_rendering) {

               if (! skip_this_ring) {
                  for (int istep=0; istep<n_steps; istep++) {
                     double angle_1 = step_frac * 2.0 * M_PI * istep;
                     double angle_2 = step_frac * 2.0 * M_PI * (istep + 1);
                     clipper::Coord_orth pt_1 = coot::util::rotate_around_vector(n, first_pt, c, angle_1);
                     clipper::Coord_orth pt_2 = coot::util::rotate_around_vector(n, first_pt, c, angle_2);
                     addBond(col, coot::Cartesian(pt_1), coot::Cartesian(pt_2), graphics_line_t::SINGLE, -1, -1, -1); // sort of, 20171224-PE FIXME needs more thought
                  }
               }
            } else {

               // for openGL rendering

               coot::torus_description_t ring(c, n, 0.07, radius, 14, 40);
               rings.push_back(ring);
            }

         } else {
            std::cout << "Not all ring atoms found in residue, for alt conf \""
                      << alt_confs[i_alt_conf] << "\", needed to draw aromatic ring: \n    ";
            std::cout << "    found " << found_atoms.size() << " but need " << ring_atom_names.size()
                      << std::endl;
            for (unsigned int i=0; i<ring_atom_names.size(); i++)
               std::cout << ":" << ring_atom_names[i] << ":  ";
            std::cout << std::endl;
         }
      }
   }
}


void
Bond_lines_container::construct_from_model_links(mmdb::Model *model_p,
                                                 int udd_atom_index_handle,
                                                 int udd_user_defined_atom_colour_index_handle,
                                                 int atom_colour_type) {

   if (false) { // debugging: show the model
      // udd_atom_index_handle is -1 for intermediate atoms
      std::cout << "in construct_from_model_links() udd_atom_index_handle is " << udd_atom_index_handle
                << "\n";

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
                  int atom_index = -1;
                  int udd_status = at->GetUDData(udd_atom_index_handle, atom_index);
                  std::cout << "  in construct_from_model_links() "
                            << coot::atom_spec_t(at) << " " << udd_status << " " << atom_index << "\n";
               }
            }
         }
      }
   }

   // Interestingly, when we add a LINK to a PDB file, if there are
   // less residues in a chain than is specified in a LINK line, mmdb
   // expands the residue list in a chain with NULL residues!

   // 20130812 Can we come here with NULL model_p? It seems that we can!
   // (somehow!)
   //
   if (model_p) {

      int n_links = model_p->GetNumberOfLinks();
      if (false)
         std::cout << "debug:: in construct_from_model_links() n_links: " << n_links << std::endl;
      if (n_links > 0) {
         for (int i_link=1; i_link<=n_links; i_link++) {
            mmdb::Link *link = model_p->GetLink(i_link);

            // For the moment, don't make Link dashed bonds to
            // symmetry-related molecules.
            //
            if ((link->s1 == link->s2) && (link->i1 == link->i2) &&
                (link->j1 == link->j2) && (link->k1 == link->k2)) {
               add_link_bond(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, atom_colour_type, link);
            }
         }
      }

      int n_linkrs = model_p->GetNumberOfLinkRs();
      if (n_linkrs > 0) {
         for (int i_link=1; i_link<=n_linkrs; i_link++) {
            mmdb::LinkR *link = model_p->GetLinkR(i_link);
            add_link_bond(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, atom_colour_type, link);
         }
      }
   }
}

void
Bond_lines_container::add_link_bond(mmdb::Model *model_p,
                                    int udd_atom_index_handle,
                                    int udd_user_defined_atom_colour_index_handle,
                                    int atom_colour_type,
                                    mmdb::Link *link) {

   if (false)
      std::cout << "calling add_link_bond with LINK "
                << "\"" << link->chainID1 << "\""
                << " "  << link->seqNum1   << " "
                << "\"" << link->atName1  << "\""
                << " to "
                << "\"" << link->chainID2 << "\""
                << " "  << link->seqNum2  << " "
                << "\"" << link->atName2  << "\""
                << std::endl;

   add_link_bond_templ(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, atom_colour_type, link);
}

void
Bond_lines_container::add_link_bond(mmdb::Model *model_p,
                                    int udd_atom_index_handle,
                                    int udd_user_defined_atom_colour_index_handle,
                                    int atom_colour_type,
                                    mmdb::LinkR *linkr) {

   // Missing LINKR bond is due to incorrect placement of atom names in the LINKR card
   if (false)
      std::cout << "calling add_link_bond_templ with LINKR "
                << "\"" << linkr->chainID1 << "\""
                << " "  << linkr->seqNum1   << " "
                << "\"" << linkr->seqNum1  << "\""
                << "\"" << linkr->atName1  << "\""
                << " to "
                << "\"" << linkr->chainID2 << "\""
                << " "  << linkr->seqNum2  << " "
                << "\"" << linkr->atName2  << "\""
                << std::endl;
   add_link_bond_templ(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, atom_colour_type, linkr);

}

template<class T>
void
Bond_lines_container::add_link_bond_templ(mmdb::Model *model_p, int udd_atom_index_handle,
                                          int udd_user_defined_atom_colour_index_handle,
                                          int atom_colour_type, T *link) {

   // std::cout << "----- add_link_bond_templ() " << atom_colour_type << std::endl;

   mmdb::PAtom atom_1 = NULL;
   mmdb::PAtom atom_2 = NULL;
   int model_number = model_p->GetSerNum();
   int n_chains = model_p->GetNumberOfChains();
   for (int ich=0; ich<n_chains; ich++) {
      mmdb::Chain *chain_p = model_p->GetChain(ich);
      if (chain_p) {
         if (std::string(chain_p->GetChainID()) == std::string(link->chainID1)) {
            int n_residues = model_p->GetNumberOfResidues();
            for (int i_res=0; i_res<n_residues; i_res++) {
               mmdb::Residue *res_p = chain_p->GetResidue(i_res);
               if (res_p) {
                  if (res_p->GetSeqNum() == link->seqNum1) {
                     if (std::string(res_p->GetInsCode()) == std::string(link->insCode1)) {
                        int n_atoms = res_p->GetNumberOfAtoms();
                        for (int iat=0; iat<n_atoms; iat++) {
                           mmdb::Atom *at = res_p->GetAtom(iat);
                           if (! at->isTer()) {
                              if (std::string(at->name) == std::string(link->atName1)) {
                                 if (std::string(at->altLoc) == std::string(link->aloc1)) {
                                    atom_1 = at;
                                    break;
                                 }
                              }
                           }
                           if (atom_1) break;
                        }
                     }
                  }
               } // null residue test
               if (atom_1) break;
            }
         }
      } // chain_p test
      if (atom_1) break;
   }

   if (atom_1) {
      for (int ich=0; ich<n_chains; ich++) {
         mmdb::Chain *chain_p = model_p->GetChain(ich);
         if (chain_p) {
            if (std::string(chain_p->GetChainID()) ==
                std::string(link->chainID2)) {
               int n_residues = model_p->GetNumberOfResidues();
               for (int i_res=0; i_res<n_residues; i_res++) {
                  mmdb::Residue *res_p = chain_p->GetResidue(i_res);
                  if (res_p) {
                     if (res_p->GetSeqNum() == link->seqNum2) {
                        if (std::string(res_p->GetInsCode()) ==
                            std::string(link->insCode2)) {
                           int n_atoms = res_p->GetNumberOfAtoms();
                           for (int iat=0; iat<n_atoms; iat++) {
                              mmdb::Atom *at = res_p->GetAtom(iat);
                              if (! at->isTer()) {
                                 if (std::string(at->name) ==
                                     std::string(link->atName2)) {
                                    if (std::string(at->altLoc) ==
                                        std::string(link->aloc2)) {
                                       atom_2 = at;
                                       break;
                                    }
                                 }
                              }
                              if (atom_2) break;
                           }
                        }
                     }
                  } // res_p test
                  if (atom_2) break;
               }
            }
            if (atom_2) break;
         } // chain_p test
      }
   }

   // OK, make the link bond then!
   if (atom_1 && atom_2) {
      int atom_index_1 = -1;
      int atom_index_2 = -1;

      // perhaps the handle should be passed, not extracted here?
      mmdb::Manager *mol = model_p->GetCoordHierarchy();

      // int udd_atom_index_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index"); // set in make_asc

      // std::cout << "DEBUG:: udd_atom_index_handle for atom index is "
      // << udd_atom_index_handle << std::endl;

      int udd_status_1 = atom_1->GetUDData(udd_atom_index_handle, atom_index_1);
      int udd_status_2 = atom_2->GetUDData(udd_atom_index_handle, atom_index_2);

#if 0 // lots of errors when drawing intermediate atoms
      if (udd_status_1 != mmdb::UDDATA_Ok) {
         std::cout << "ERROR:: in add_link_bond_templ() bad atom indexing 1 using udd_atom_index_handle "
                   << udd_atom_index_handle << std::endl;
      }
      if (udd_status_2 != mmdb::UDDATA_Ok) {
         std::cout << "ERROR:: in add_link_bond_templ() bad atom indexing 2 using udd_atom_index_handle "
                   << udd_atom_index_handle << std::endl;
      }
      if (no_bonds_to_these_atoms.find(atom_index_1) != no_bonds_to_these_atoms.end())
         std::cout << "Debug atom_index_1 " << atom_index_1 << " not to be excluded ";
      else
         std::cout << "Debug atom_index_1 " << atom_index_1 << " should be excluded ";
      if (no_bonds_to_these_atoms.find(atom_index_2) != no_bonds_to_these_atoms.end())
         std::cout << " atom_index_2 " << atom_index_2 << " not to be excluded\n";
      else
         std::cout << " atom_index_2 " << atom_index_2 << " should be excluded\n";
#endif



      // Even if the atom_index_1 or atom_index_2 were not correctly set, we can still draw
      // the bond - this needs to be fixed however.
      coot::Cartesian pos_1(atom_1->x, atom_1->y, atom_1->z);
      coot::Cartesian pos_2(atom_2->x, atom_2->y, atom_2->z);

      std::string ele_1 = atom_1->element;
      std::string ele_2 = atom_2->element;
      if (ele_1 == ele_2) {
         int col = atom_colour(atom_1, atom_colour_type, udd_user_defined_atom_colour_index_handle);
         add_dashed_bond(col, pos_1, pos_2, NOT_HALF_BOND, graphics_line_t::SINGLE, model_number, atom_index_1, atom_index_2);
      } else {
         coot::Cartesian bond_mid_point = pos_1.mid_point(pos_2);
         int col = atom_colour(atom_1, atom_colour_type, udd_user_defined_atom_colour_index_handle);
         // if the atom indices are -1, then the bond doesn't get drawn
         add_dashed_bond(col, pos_1, bond_mid_point, HALF_BOND_FIRST_ATOM, graphics_line_t::SINGLE, model_number, atom_index_1, atom_index_2);
         col = atom_colour(atom_2, atom_colour_type, udd_user_defined_atom_colour_index_handle);
         add_dashed_bond(col, bond_mid_point, pos_2, HALF_BOND_SECOND_ATOM, graphics_line_t::SINGLE, model_number, atom_index_1, atom_index_2);
      }
   } else {
#if 0
      // 20251219-PE because all of the links are copied, we often find that some links are not in the selection
      // for the sphere refinement. We don't want to hear that. Ideally, we'd not even get here in such cases.
      if (! atom_1)
         std::cout << "debug:: in add_link_bond_templ() failed to find atom-1 "
                   << "\"" << link->chainID1 << "\" " << link->seqNum1 << " \"" << link->atName1 << "\"" << std::endl;
      if (! atom_2)
         std::cout << "debug:: in add_link_bond_templ() failed to find atom-2 "
                   << "\"" << link->chainID2 << "\" " << link->seqNum2 << " \"" << link->atName2 << "\"" << std::endl;
#endif
   }
}



mmdb::PPAtom
coot::model_bond_atom_info_t::Hydrogen_atoms() const {

   mmdb::PPAtom H_atoms = new mmdb::PAtom[hydrogen_atoms_.size()];
   for (unsigned int i=0; i<hydrogen_atoms_.size(); i++) {
      H_atoms[i] = hydrogen_atoms_[i];
   }
   return H_atoms;
}

mmdb::PPAtom
coot::model_bond_atom_info_t::non_Hydrogen_atoms() const {

   mmdb::PPAtom non_H_atoms = new mmdb::PAtom[non_hydrogen_atoms_.size()];
   for (unsigned int i=0; i<non_hydrogen_atoms_.size(); i++) {
      non_H_atoms[i] = non_hydrogen_atoms_[i];
   }
   return non_H_atoms;
}

void
Bond_lines_container::atom_selection_missing_loops(const atom_selection_container_t &asc,
                                                   int udd_atom_index_handle, int udd_fixed_during_refinement_handle) {

   // this is not called if draw_missing_loops_flag is false

   mmdb::Manager *mol = asc.mol;
   int imodel= 1;
   mmdb::Model *model_p = mol->GetModel(imodel);
   if (model_p) {
      int n_chains= model_p->GetNumberOfChains();
      for (int ich=0; ich<n_chains; ich++) {
         mmdb::Chain *chain_p = model_p->GetChain(ich);
         if (! chain_p) continue;
         int nres = chain_p->GetNumberOfResidues();
         if (nres < 2) continue;
         for (int ires=1; ires<nres; ires++) {
            mmdb::Residue *residue_this = chain_p->GetResidue(ires);
            mmdb::Residue *residue_prev = chain_p->GetResidue(ires-1);
            if (residue_this && residue_prev) {

               int n_atoms_prev = residue_prev->GetNumberOfAtoms();
               int n_atoms_this = residue_this->GetNumberOfAtoms();

               if (n_atoms_prev == 0) continue;
               if (n_atoms_this == 0) continue;

               int res_no_1 = residue_prev->GetSeqNum();
               int res_no_2 = residue_this->GetSeqNum();
               int res_no_delta = res_no_2 - res_no_1;
               if (res_no_delta > 1) {
                     do_Ca_loop(imodel, ires, nres, chain_p, residue_prev, residue_this,
                                udd_atom_index_handle, udd_fixed_during_refinement_handle);
               }
            }
         }
      }
   }
}

// we rely on SelAtom.atom_selection being properly constructed to
// contain all atoms
void
Bond_lines_container::construct_from_asc(const atom_selection_container_t &SelAtom,
                                         int imol,
                                         float min_dist, float max_dist,
                                         int atom_colour_type,
                                         short int is_from_symmetry_flag,
                                         bool draw_missing_loops_flag,
                                         int model_number,
                                         bool do_rama_markup, // default false
                                         bool do_rota_markup) { // default false

   if (false)
      std::cout << "DEBUG:: construct_from_asc() was called with model_number " << model_number << std::endl;

   // initialize each colour in the Bond_lines_container
   //
   // There are now 13 colours in bond_colours (CPK extras)
   int n_col = 13;
   if (atom_colour_type == coot::COLOUR_BY_USER_DEFINED_COLOURS)
      n_col = 40;
   if (atom_colour_type == coot::COLOUR_BY_B_FACTOR)
      n_col = 50;

   if (bonds.size() == 0) {
      for (int i=0; i<n_col; i++) {
         Bond_lines a(i);
         bonds.push_back(a);
      }
   }
   float star_size = 0.22;

   // initialize the hydrogen bonding flag:
   //
   // do_bonds_to_hydrogens = 1;  // Fix off, is this sensible?  Valgrind
                               // // complained (rightly) about jump
                               // 80-90 lines below.
                               // (do_bonds_to_hydrogens had be unitialized).
   //
                               // 20060812: Actually, I now doubt that
                               // this is uninitialized.  The value is
                               // set from/in the constructor from the
                               // molecule_class_info_t value
                               // draw_hydrogens_flag, which is set in
                               // the molecule_class_info_t
                               // constructors.  So I don't know what
                               // is going on... I will comment out
                               // this line now and see if valgrind
                               // still complains.
                               //
                               // OK, I just checked. Valgrind no
                               // longer complains if I comment this line.
   // 20070407 Back here again! Valgrind *does* complain about
   // uninitialized do_bonds_to_hydrogens when I enable NCS ghosts.

   // and disulfides bond flag:
   // do_disulfide_bonds_flag = 1; // 20211005-PE - no, this should be set in the constructor, not here

   if (SelAtom.n_selected_atoms <= 0) return;

   // count the number of hydrogen atoms and non-hydrogen atoms:
   //
//    std::string element;
//    int n_H = 0;
//    int n_non_H = 0;
//    for (int i=0; i<SelAtom.n_selected_atoms; i++) {
//       element = SelAtom.atom_selection[i]->element;
//       if (element == " H" || element == " D") {
//          n_H++;
//       } else {
//          n_non_H++;
//       }
//    }

   // Now, let's not forget that some atoms don't have contacts, so
   // whenever we find a contact for an atom, we mark it with
   // UserDefinedData "found bond".
   //

   // 20190904-PE Note to self udd_atom_index_handle is -1 for an intermediate atoms asc
   //             All the no_bonds_to_these_atoms test will fail the atom index lookup.
   //             (that's OK because we don't want to exclude atoms from the intermediate
   //             atoms)

   int udd_atom_index_handle = SelAtom.UDDAtomIndexHandle;

   int udd_found_bond_handle = SelAtom.mol->RegisterUDInteger(mmdb::UDR_ATOM, "found bond");

   int udd_user_defined_atom_colour_index_handle = SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   // std::cout << "............................ construct_from_asc() "
   // << "udd_found_bond_handle register " << udd_found_bond_handle << std::endl;
   bool have_udd_atoms = 1;
   if (udd_found_bond_handle<0)  {
      std::cout << " atom bonding registration failed.\n";
      have_udd_atoms = 0;
   } else {
      for (int i=0; i<SelAtom.n_selected_atoms; i++) {
	 mmdb::Atom *at = SelAtom.atom_selection[i];
	 if (at)
	    at->PutUDData(udd_found_bond_handle, graphical_bonds_container::NO_BOND);
	 else
	    logger.log(log_t::ERROR, logging::function_name_t("construct_from_asc"), "Bad atom!", i);
      }
   }

   // We need to loop over each model.  Currently the SelAtom is not
   // selected on model.
   //
   // consider a vector (over all models) of these:
   //
   // class model_bond_atom_info_t {
   //       std::vector<mmdb::PAtom> hydrogen_atoms_;
   //       std::vector<mmdb::PAtom> non_hydrogen_atoms_;
   //    public:
   //    mmdb::PPAtom     Hydrogen_atoms();
   //    mmdb::PPAtom non_Hydrogen_atoms();
   //    int n_H();
   //    int n_non_H();
   // }

   int imodel = -1;
   int n_models = SelAtom.mol->GetNumberOfModels(); // models start at number 1
   std::vector<coot::model_bond_atom_info_t> atom_stuff_vec(n_models+1);

   for (int i=0; i<SelAtom.n_selected_atoms; i++) {
      mmdb::Atom *at = SelAtom.atom_selection[i];
      if (at) {
	 imodel = at->GetModelNum();
	 if ((imodel <= n_models) && (imodel > 0)) {
	    atom_stuff_vec[imodel].add_atom(SelAtom.atom_selection[i]);
	 }
      }
   }

   // default to all models:
   int imodel_start = 1;
   int imodel_end = n_models;
   // over-ride if we were passed model_number - in that case, just do
   // the following look once - for the given model number.
   if (model_number != 0) {
      if (model_number <= n_models) {
         imodel_start = model_number;
         imodel_end   = model_number;
      }
   }

   for (int imodel=imodel_start; imodel<=imodel_end; imodel++) {

      if (imodel == 1) {
         int udd_fixed_during_refinement_handle = SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM, "FixedDuringRefinement");
         if (draw_missing_loops_flag)
            atom_selection_missing_loops(SelAtom, udd_atom_index_handle, udd_fixed_during_refinement_handle);
      }

      mmdb::PPAtom     Hydrogen_atoms = atom_stuff_vec[imodel].Hydrogen_atoms();
      mmdb::PPAtom non_Hydrogen_atoms = atom_stuff_vec[imodel].non_Hydrogen_atoms();
      int n_non_H = atom_stuff_vec[imodel].n_non_H();
      int n_H = atom_stuff_vec[imodel].n_H();

      construct_from_atom_selection(SelAtom,
                                    non_Hydrogen_atoms, n_non_H,
                                    non_Hydrogen_atoms, n_non_H,
                                    imol,
                                    min_dist, max_dist, atom_colour_type,
                                    0, have_udd_atoms, udd_found_bond_handle);

      if (do_bonds_to_hydrogens && (n_H > 0)) {

         float H_min_dist = 0.7;
         float H_max_dist = 1.42;

         // H-H
         construct_from_atom_selection(SelAtom,
                                       Hydrogen_atoms, n_H,
                                       Hydrogen_atoms, n_H,
                                       imol,
                                       H_min_dist, H_max_dist, atom_colour_type,
                                       0, have_udd_atoms, udd_found_bond_handle);

         // H-X
         construct_from_atom_selection(SelAtom,
                                       non_Hydrogen_atoms, n_non_H,
                                       Hydrogen_atoms, n_H,
                                       imol,
                                       H_min_dist, H_max_dist, atom_colour_type,
                                       1, have_udd_atoms, udd_found_bond_handle);
      }

      mmdb::Model *model_p = SelAtom.mol->GetModel(imodel);

      if (model_p)
         construct_from_model_links(model_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, atom_colour_type);

      // std::cout << "DEBUG:: (post) SelAtom: mol, n_selected_atoms "
      // << SelAtom.mol << " " << SelAtom.n_selected_atoms << std::endl;

      if (do_disulfide_bonds_flag)
         do_disulphide_bonds(SelAtom, imodel);

      // now we have dealt with the bonded atoms, lets find the non-bonded
      // atoms...

      if (have_udd_atoms) {

         // for atoms with no neighbour (contacts):
         coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
         coot::Cartesian small_vec_y(0.0, star_size, 0.0);
         coot::Cartesian small_vec_z(0.0, 0.0, star_size);

         int ic; // changed by reference (UDData)
         int col;

         for (int i=0; i<n_non_H; i++) {

            if (non_Hydrogen_atoms[i]->GetUDData(udd_found_bond_handle, ic) == mmdb::UDDATA_Ok) {
               if ((ic == 0) ||
                   ((!strcmp(non_Hydrogen_atoms[i]->element, " S")) && (ic != graphical_bonds_container::BONDED_WITH_HETATM_BOND)) ||
                   ((!strcmp(non_Hydrogen_atoms[i]->element, "SE")) && (ic != graphical_bonds_container::BONDED_WITH_HETATM_BOND)) ||
                   ((!strcmp(non_Hydrogen_atoms[i]->element, "FE")) && (ic != graphical_bonds_container::BONDED_WITH_HETATM_BOND)) ||
                   ((!strcmp(non_Hydrogen_atoms[i]->element, " P")) && (ic != graphical_bonds_container::BONDED_WITH_HETATM_BOND))) {

                  // std::cout << "::::  No contact for " << non_Hydrogen_atoms[i]
                  //           << " with ic " << ic << std::endl;

                  // no contact found or was Sulphur, or Phosphor

                  // So, was this a seleno-methione?
                  //
                  mmdb::Residue *atom_residue_p = non_Hydrogen_atoms[i]->residue;
                  if (atom_residue_p) {

                     std::string resname = non_Hydrogen_atoms[i]->GetResName();
                     if ((is_from_symmetry_flag == 0) &&
                         (resname == "MSE" || resname == "MET" || resname == "MSO"
                          || resname == "CYS" )) {
                        handle_MET_or_MSE_case(non_Hydrogen_atoms[i], udd_found_bond_handle, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle,
                                               atom_colour_type);
                     } else {

                        std::vector<std::pair<bool, mmdb::Residue *> > het_residues; // bond these separately.
                        mmdb::Atom *atom_p_1 = non_Hydrogen_atoms[i];
                        bool bond_het_residue_by_dictionary = add_bond_by_dictionary_maybe(imol, atom_p_1, atom_p_1, &het_residues);
                        if (bond_het_residue_by_dictionary) {
                           add_bonds_het_residues(het_residues, SelAtom, imol, atom_colour_type, have_udd_atoms, udd_found_bond_handle,
                                                  udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
                        } else {

                           std::string ele = non_Hydrogen_atoms[i]->element;
                           if (ele == "CL" || ele == "BR" || ele == " S" ||  ele == " I"
                               || ele == "Cl" || ele == "Br" || ele == "MO" || ele == "Mo" || ele == "AL"
                               || ele == "PT" || ele == "RU" || ele == " W"
                               || ele == "AS" || ele == " P" || ele == "AU" || ele == "HG"
                               || ele == "PD" || ele == "PB" || ele == "AG") {
                              handle_long_bonded_atom(non_Hydrogen_atoms[i], udd_found_bond_handle, udd_atom_index_handle,
                                                      udd_user_defined_atom_colour_index_handle,
                                                      atom_colour_type);
                           }
                        }
                     }
                  } else {
                     std::cout << "INFO:: trapped atom without residue in non-bonded atom check: "
                               << non_Hydrogen_atoms[i] << std::endl;
                  }
               }
            } else {
               std::cout << "missing UDData for atom "
                         << non_Hydrogen_atoms[i] << "\n";
            }
         }


         // Make the stars...
         //
         for (int i=0; i<n_non_H; i++) {
            if (non_Hydrogen_atoms[i]->GetUDData(udd_found_bond_handle, ic) == mmdb::UDDATA_Ok) {
               if (ic == graphical_bonds_container::NO_BOND) {
                  // no contact found
                  mmdb::Residue *residue_p = non_Hydrogen_atoms[i]->residue;

                  std::string res_name(residue_p->GetResName());
                  if (res_name == "HOH")
                     if (! do_sticks_for_waters)
                        continue;

                  col = atom_colour(non_Hydrogen_atoms[i], atom_colour_type, udd_user_defined_atom_colour_index_handle);
                  coot::Cartesian atom_pos(non_Hydrogen_atoms[i]->x,
                                           non_Hydrogen_atoms[i]->y,
                                           non_Hydrogen_atoms[i]->z);

                  int iat_1 = -1;
                  int udd_status_1 = non_Hydrogen_atoms[i]->GetUDData(udd_atom_index_handle, iat_1);

                  graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                  addBond(col, atom_pos+small_vec_x, atom_pos-small_vec_x, cc, model_number, iat_1, iat_1, true, true);
                  addBond(col, atom_pos+small_vec_y, atom_pos-small_vec_y, cc, model_number, iat_1, iat_1, true, true);
                  addBond(col, atom_pos+small_vec_z, atom_pos-small_vec_z, cc, model_number, iat_1, iat_1, true, true);
               }
            }
         }

         if (do_bonds_to_hydrogens && (n_H > 0)) {
            for (int i=0; i<n_H; i++) {
               if (Hydrogen_atoms[i]->GetUDData(udd_found_bond_handle, ic) == mmdb::UDDATA_Ok) {
                  if (ic == graphical_bonds_container::NO_BOND) {

                     // no contact found
                     mmdb::Residue *residue_p = Hydrogen_atoms[i]->residue;
                     if (! residue_p)
                        std::cout << "ERROR:: catched condition for crashetty crash!" << std::endl;

                     std::string res_name(residue_p->GetResName());
                     if (res_name == "HOH")
                        if (! do_sticks_for_waters)
                           continue;

                     col = atom_colour(Hydrogen_atoms[i], atom_colour_type, udd_user_defined_atom_colour_index_handle);
                     coot::Cartesian atom(Hydrogen_atoms[i]->x,
                                          Hydrogen_atoms[i]->y,
                                          Hydrogen_atoms[i]->z);

                     // 20171224-PE FIXME by lookup
                     int iat_1 = -1;
                     int udd_status_1 = Hydrogen_atoms[i]->GetUDData(udd_atom_index_handle, iat_1);
                     graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                     addBond(col, atom+small_vec_x, atom-small_vec_x, cc, model_number, iat_1, iat_1);
                     addBond(col, atom+small_vec_y, atom-small_vec_y, cc, model_number, iat_1, iat_1);
                     addBond(col, atom+small_vec_z, atom-small_vec_z, cc, model_number, iat_1, iat_1);
                  }
               }
            }
         }
      }

      delete [] Hydrogen_atoms;
      delete [] non_Hydrogen_atoms;
   }

   add_zero_occ_spots(SelAtom);
   add_deuterium_spots(SelAtom);

   if (do_rama_markup) {
      // std::cout << "........... add_ramachandran_goodness_spots() " << std::endl;
      add_ramachandran_goodness_spots(SelAtom);
   }
   if (do_rota_markup) {
      // std::cout << "........... add_rotamer_goodness_spots() " << std::endl;
      add_rotamer_goodness_markup(SelAtom);
   }
   add_atom_centres(imol, SelAtom, atom_colour_type, model_number);
   add_cis_peptide_markup(SelAtom, model_number);
}

void
Bond_lines_container::handle_MET_or_MSE_case(mmdb::PAtom mse_atom,
                                             int udd_handle_bond, // udd for having bond assignment
                                             int udd_handle_atom_index,
                                             int udd_user_defined_atom_colour_index_handle,
                                             int atom_colour_type,
                                             coot::my_atom_colour_map_t *atom_colour_map_p) {

   // std::cout << "Handling MET/MSE case for atom " << mse_atom << std::endl;

   std::string atom_name(mse_atom->name);
   std::string residue_name(mse_atom->GetResName());
   int model_number = mse_atom->GetModelNum();
   if (residue_name == "MET" || residue_name == "MSE" || residue_name == "MSO") {
      if (atom_name == "SE  " || atom_name == " SD ") {
         int col = atom_colour(mse_atom, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);

         // We need to add special bonds SE -> CE and SE -> CG.
         mmdb::PPAtom residue_atoms;
         int nResidueAtoms;
         mse_atom->residue->GetAtomTable(residue_atoms, nResidueAtoms);
         for (int i=0; i<nResidueAtoms; i++) {
            std::string table_atom_name(residue_atoms[i]->name);
            if (table_atom_name == " CG " ||
                table_atom_name == " CE " ) {
               // mse_atom or met_atom now of course.
               coot::Cartesian cart_at1(mse_atom->x, mse_atom->y, mse_atom->z);
               coot::Cartesian cart_at2(residue_atoms[i]->x,
                                        residue_atoms[i]->y,
                                        residue_atoms[i]->z);

               std::string altconf1 = mse_atom->altLoc;
               std::string altconf2 = residue_atoms[i]->altLoc;
               if ( (altconf1=="") || (altconf2=="") || (altconf1==altconf2) ) {
                  coot::Cartesian bond_mid_point = cart_at1.mid_point(cart_at2);
                  int colc = atom_colour(residue_atoms[i], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
                  // int colc = atom_colour_type; // just to check

                  int iat_1 = -1;
                  int iat_2 = -1; // 20171224-PE FIXME by udd lookup
                  int udd_status_1 = mse_atom->GetUDData(udd_handle_atom_index, iat_1);
                  int udd_status_2 = residue_atoms[i]->GetUDData(udd_handle_atom_index, iat_2);

                  float bond_length = (cart_at1 - cart_at2).amplitude();
                  if (bond_length < 3.0) { // surely not longer than this?
                     graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                     addBond(col,  cart_at1, bond_mid_point, cc, model_number, iat_1, iat_2); // The SE->mid_pt bond
                     addBond(colc, bond_mid_point, cart_at2, cc, model_number, iat_1, iat_2); // The Cx->mid_pt bond
                     // mark atom as bonded.
                     residue_atoms[i]->PutUDData(udd_handle_bond, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                     mse_atom->PutUDData(udd_handle_bond, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                  }
               }
            }
         }
      }
   }
   if (residue_name == "CYS") {
      int col = atom_colour(mse_atom, atom_colour_type, udd_user_defined_atom_colour_index_handle,atom_colour_map_p);

      if (atom_name == " SG ") {
         // We need to add special bonds CB -> SG
         mmdb::PPAtom residue_atoms;
         int nResidueAtoms;
         mse_atom->residue->GetAtomTable(residue_atoms, nResidueAtoms);
         for (int i=0; i<nResidueAtoms; i++) {
            std::string table_atom_name(residue_atoms[i]->name);
            if (table_atom_name == " CB ") {
               coot::Cartesian cart_at1(mse_atom->x, mse_atom->y, mse_atom->z);
               coot::Cartesian cart_at2(residue_atoms[i]->x,
                                        residue_atoms[i]->y,
                                        residue_atoms[i]->z);

               std::string altconf1 = mse_atom->altLoc;
               std::string altconf2 = residue_atoms[i]->altLoc;
               if ( (altconf1=="") || (altconf2=="") || (altconf1==altconf2) ) {
                  float len2 = (cart_at1 - cart_at2).amplitude_squared();
                  if (len2 < 16) { // protection for weirdness
                     coot::Cartesian bond_mid_point = cart_at1.mid_point(cart_at2);
                     int colc = atom_colour(residue_atoms[i], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
                     graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

                     // 20181231-PE fix up UDData for atom indices at last
                     int iat_1 = -1;
                     int iat_2 = -1;
                     mse_atom->GetUDData(udd_handle_atom_index, iat_1);
                     residue_atoms[i]->GetUDData(udd_handle_atom_index, iat_2);
                     addBond(col,  cart_at1, bond_mid_point, cc, model_number, iat_1, iat_2);
                     addBond(colc, bond_mid_point, cart_at2, cc, model_number, iat_1, iat_2);
                     // mark atom as bonded.
                     residue_atoms[i]->PutUDData(udd_handle_bond, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                     mse_atom->PutUDData(udd_handle_bond, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                  }
               }
            }
         }
      }
   }
}

void
Bond_lines_container::handle_long_bonded_atom(mmdb::PAtom atom,
                                              int udd_handle_bond,
                                              int udd_handle_atom_index,
                                              int udd_user_defined_atom_colour_index_handle,
                                              int atom_colour_type,
                                              coot::my_atom_colour_map_t *atom_colour_map_p) {

   // std::cout << "Here in handle_long_bonded_atom() " << coot::atom_spec_t(atom) << std::endl;

   float bond_limit = 2.16; // A S-S bonds are 2.05A.  So we've added
                            // some wiggle room (2.1 was too short for
                            // some dictionary S-S).

   std::string atom_name(atom->name);
   std::string residue_name(atom->GetResName());
   std::string element(atom->element);
   mmdb::Residue *res = atom->residue;
   int model_number = atom->GetModelNum();

   // std::cout << "handling long bonds for " << atom << " ele " << element << std::endl;

   if (atom_name == "AS  ")
      bond_limit = 2.4;
   if (element == "AU")
      bond_limit = 2.4;
   if (element == "AS")
      bond_limit = 2.4;
   if (element == "HG")
      bond_limit = 2.4;
   if (element == "MO")
      bond_limit = 2.55; // was 0.23 Oxygen to Molybdenum, guess from 8M0
   if (element == "Mo")
      bond_limit = 2.55; // 4077716
   if (element == " I")
      bond_limit = 2.3; // C-I is 2.13 according to wikipedia
   float  bl2 = bond_limit * bond_limit;
   float h_bl2 = 1.8 * 1.8; // 20100208 all bonds to hydrogens are less than 1.8A ? (guess)
   short int bond_added_flag = 0;

   if (res) {
      // do the bonding by hand:
      coot::Cartesian atom_pos(atom->x, atom->y, atom->z);
      int col = atom_colour(atom, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      mmdb::PPAtom residue_atoms = 0;
      int nResidueAtoms;
      res->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int i=0; i<nResidueAtoms; i++) {
         if (residue_atoms[i] != atom) {
            coot::Cartesian res_atom_pos(residue_atoms[i]->x,
                                         residue_atoms[i]->y,
                                         residue_atoms[i]->z);

            // We compared squard bond distances (so that we don't
            // have to take the square root for everything of course).
            //
            std::string res_atom_ele = residue_atoms[i]->element;
            float len2 = (atom_pos - res_atom_pos).amplitude_squared();
            if (((len2 <   bl2) && (! is_hydrogen(res_atom_ele))) ||
                ((len2 < h_bl2) && (is_hydrogen(res_atom_ele)))) {
               std::string altconf1 = atom->altLoc;
               std::string altconf2 = residue_atoms[i]->altLoc;
               if ( (altconf1=="") || (altconf2=="") || (altconf1==altconf2) ) {
                  coot::Cartesian bond_mid_point = atom_pos.mid_point(res_atom_pos);
                  int colc = atom_colour(residue_atoms[i], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
                  graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

                  // 20171224-PE FIXME lookup iat_1, iat_2
                  int iat_1 = -1;
                  int iat_2 = -1;

                  // 20181231-PE Done.
                  int udd_status_1 = atom->GetUDData(udd_handle_atom_index, iat_1);
                  int udd_status_2 = residue_atoms[i]->GetUDData(udd_handle_atom_index, iat_2);

                  if (true) { // for debugging
                     addBond(col,  atom_pos, bond_mid_point, cc, model_number, iat_1, iat_2);
                     addBond(colc, bond_mid_point, res_atom_pos, cc, model_number, iat_2, iat_2);
                     bond_added_flag = 1;
                     // mark the atom as bonded.
                     atom->PutUDData(udd_handle_bond, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                     residue_atoms[i]->PutUDData(udd_handle_bond, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                  }
               }
            }
         }
      }
   }

   if (!bond_added_flag) {
      // bond it like a single atom then:

      float star_size = 0.3;
      coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
      coot::Cartesian small_vec_y(0.0, star_size, 0.0);
      coot::Cartesian small_vec_z(0.0, 0.0, star_size);

      int col = atom_colour(atom, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
      coot::Cartesian atom_pos(atom->x, atom->y, atom->z);
      graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
      // 20171224-PE FIXME lookup iat_1, iat_1
      int iat_1 = -1;
      // 20190921-PE Done.
      int udd_status = atom->GetUDData(udd_handle_atom_index, iat_1);
      addBond(col, atom_pos+small_vec_x, atom_pos-small_vec_x, cc, model_number, iat_1, iat_1, true, true);
      addBond(col, atom_pos+small_vec_y, atom_pos-small_vec_y, cc, model_number, iat_1, iat_1, true, true);
      addBond(col, atom_pos+small_vec_z, atom_pos-small_vec_z, cc, model_number, iat_1, iat_1, true, true);
   }
}





// This finds bonds between a residue and the protein (in SelAtom).
// It is used for the environment bonds box.
//
// What do we do about drawing hydrogens?
//
Bond_lines_container::Bond_lines_container(const atom_selection_container_t &SelAtom,
                                           mmdb::PPAtom residue_atoms,
                                           int n_residue_atoms,
                                           coot::protein_geometry *protein_geom_p, // modifiable
                                           bool residue_is_water_flag,
                                           bool draw_env_distances_to_hydrogens_flag,
                                           float min_dist,
                                           float max_dist) {

   if (0)
      std::cout << "Environment distances NO symm" << std::endl;
   do_bonds_to_hydrogens = 1;  // added 20070629

   b_factor_scale = 1.0;
   have_dictionary = 0;
   for_GL_solid_model_rendering = 1;
   n_atoms_in_atom_selection = SelAtom.n_selected_atoms;

   int ncontacts;
   mmdb::Contact *contact = NULL;
   // initialize each colour in the Bond_lines_container
   //
   if (bonds.size() == 0) {
      bonds.resize(13); // There are now 13 colours in bond_colours
      for (int i=0; i<13; i++)
         bonds[i] = Bond_lines(i);
   }


   // seqDist has to be 0 here or else failure to find contacts to new
   // waters added as pointer atoms. Why?  Don't know - a bug in mmdb,
   // I suspect.
   //
   SelAtom.mol->SeekContacts(residue_atoms, n_residue_atoms,
                             SelAtom.atom_selection, SelAtom.n_selected_atoms,
                             min_dist,
                             max_dist,
                             0,  // seqDist (in same residue allowed)
                             contact, ncontacts);

   if (0) {  // debugging seqDist
      std::cout << " DEBUG:: there are " << n_residue_atoms << " residue atoms "
                << " and " << SelAtom.n_selected_atoms << " mol atoms\n";
      for (int iat=0; iat<n_residue_atoms; iat++)
         std::cout << "residue atoms: " << iat << "/" << n_residue_atoms
                   << " " << residue_atoms[iat] << std::endl;
      for (int iat=0; iat<SelAtom.n_selected_atoms; iat++)
         std::cout << "    mol atoms: " << iat << "/" << SelAtom.n_selected_atoms
                   << " " << SelAtom.atom_selection[iat] << std::endl;
   }


   if (ncontacts > 0) {
      for (int i=0; i<ncontacts; i++) {

         int iat_1 = contact[i].id1;
         int iat_2 = contact[i].id2;

         if (draw_these_atom_contacts(residue_atoms[contact[i].id1], SelAtom.atom_selection[contact[i].id2],
                                      protein_geom_p) || residue_is_water_flag) {

            mmdb::Atom *atom_1 = residue_atoms[ contact[i].id1 ];
            mmdb::Atom *atom_2 = SelAtom.atom_selection[ contact[i].id2];

            coot::Cartesian atom_1_pos(residue_atoms[ contact[i].id1 ]->x,
                                       residue_atoms[ contact[i].id1 ]->y,
                                       residue_atoms[ contact[i].id1 ]->z);
            coot::Cartesian atom_2_pos(SelAtom.atom_selection[ contact[i].id2 ]->x,
                                       SelAtom.atom_selection[ contact[i].id2 ]->y,
                                       SelAtom.atom_selection[ contact[i].id2 ]->z);
            std::string ele1 = residue_atoms[ contact[i].id1 ]->element;
            std::string ele2 = SelAtom.atom_selection[ contact[i].id2 ]->element;
            std::string alt_conf_1 = residue_atoms[ contact[i].id1 ]->altLoc;
            std::string alt_conf_2 = SelAtom.atom_selection[ contact[i].id2 ]->altLoc;
            int model_number = residue_atoms[ contact[i].id1 ]->GetModelNum();

            // 20110119 environment distances for Hydrogens.  We don't
            // want to see a big web of contacts between hydrogens and
            // other atoms (including hydrogens and other hydrogens),
            // so if the atom is a hydrogen, only "mark" it if below
            // bonding_dist_max (and adjust bonding_dist_max to be
            // smaller if this atom is a hydrogen).
            //
            double bonding_dist_max = max_dist;
            // double shorter_bit = 0.38; // tinkered value
            double shorter_bit = 0.52; // tinkered value
            //
            if (is_hydrogen(ele1))
               bonding_dist_max -= shorter_bit;
            if (is_hydrogen(ele2))
               bonding_dist_max -= shorter_bit;

            if (0) { // debug
               std::cout << " DEBUG:: add environ dist "
                         << residue_atoms[ contact[i].id1 ] << " to "
                         << SelAtom.atom_selection[ contact[i].id2 ]
                         << std::endl;
            }

            double dist = coot::distance(residue_atoms[contact[i].id1],
                                         SelAtom.atom_selection[contact[i].id2]);

            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

            if (dist <= bonding_dist_max) {
               if ((alt_conf_1 == alt_conf_2) ||
                   (alt_conf_1 == "") ||
                   (alt_conf_2 == "")) {
                  if (draw_env_distances_to_hydrogens_flag ||
                      // ((ele1 != " H") && (ele2 != " H"))) {
                      ((! is_hydrogen(ele1)) && (! is_hydrogen(ele2)))) {
                     if (ele1 == " C")
                        addBond(0, atom_1_pos, atom_2_pos, cc, model_number, iat_1, iat_2);
                     else {
                        if (ele2 == " C") {
                           addBond(0, atom_1_pos, atom_2_pos, cc, model_number, iat_1, iat_2);
                        } else {

                           // both atoms not Carbon

                           // if (ele1 == " H" && ele2 == " H") {
                           if (is_hydrogen(ele1) && is_hydrogen(ele2)) {
                              addBond(0, atom_1_pos, atom_2_pos, cc, model_number, iat_1, iat_2); // not a charged/H-bond
                           } else {

                              // stop purple lines between (say) OD2 and mainchain O
                              coot::quick_protein_donor_acceptors pda;
                              coot::quick_protein_donor_acceptors::key k1(atom_1->GetResName(), atom_1->GetAtomName());
                              coot::quick_protein_donor_acceptors::key k2(atom_2->GetResName(), atom_2->GetAtomName());
                              // is-looked-up, is-H-bond
                              int colour_index = 1; // H-bond
                              std::pair<bool,bool> is_valid = pda.is_hydrogen_bond_by_types(k1,k2);
                              if (is_valid.first)
                                 if (! is_valid.second)
                                    colour_index = 0;
                              addBond(colour_index, atom_1_pos, atom_2_pos, cc, model_number, iat_1, iat_2); // interesting
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      delete [] contact;
   }
}

// Shall we draw environment bonds between this_residue and env_residue?
//
// No, if the residues are next to each other in sequence in the same
// chain and are a polymer.
//
// Otherwise, yes.
//
bool
Bond_lines_container::draw_these_residue_contacts(mmdb::Residue *this_residue,
                                                  mmdb::Residue *env_residue,
                                                  coot::protein_geometry *protein_geom_p // modifiable
                                                  ) {

   if (this_residue != env_residue) {
      std::string ch1(this_residue->GetChainID());
      std::string ch2(env_residue->GetChainID());
      if ((abs(this_residue->GetSeqNum() - env_residue->GetSeqNum()) > 1)
          || (ch1 != ch2)) {
         return true;
      } else {
         // are we in a polymer? if so, no draw.
         //
         std::string this_res_type = this_residue->GetResName();
         std::string env_residue_res_type = env_residue->GetResName();
         if (protein_geom_p->linkable_residue_types_p(this_res_type, env_residue_res_type)) {
            return false;
         } else {
            return true;
         }
      }
   } else {
      return false;
   }
}


// We want to filter out contact in the same residue.
//
// we want to filter out atom contacts along the main chain.  Previously we did that by
// checking that the residues were not next to each other (above) - but I want to see contacts
// between bases in DNA, so now, filter out distances based on atom names (and residue numbering)
//
bool
Bond_lines_container::draw_these_atom_contacts(mmdb::Atom *this_atom, mmdb::Atom *env_atom,
                                               coot::protein_geometry *protein_geom) {

   bool draw_flag = true;

   mmdb::Residue *this_residue = this_atom->GetResidue();
   mmdb::Residue *env_residue  =  env_atom->GetResidue();

   mmdb::Chain *ch_this = this_atom->GetChain();
   mmdb::Chain *ch_env  =  env_atom->GetChain();

   if (ch_this != ch_env) {
      return true;
   } else {
      if (this_residue == env_residue) {
         return false;
      } else {
         if (abs(this_residue->GetSeqNum() - env_residue->GetSeqNum()) > 1) {
            return true;
         } else {
            // OK, we have neighbouring residues in the same chain
            //
            std::string this_res_type = this_residue->GetResName();
            std::string env_residue_res_type = env_residue->GetResName();
            if (! protein_geom->linkable_residue_types_p(this_res_type, env_residue_res_type)) {
               return true;
            } else {
               std::string this_atom_name = this_atom->GetAtomName();
               std::string  env_atom_name =  env_atom->GetAtomName();
               // PDBv3 FIXME
               if (this_atom_name == " N  ") if (env_atom_name == " CA ") draw_flag = false;

               if ((this_atom_name == " N  ") || (this_atom_name == " CA ") ||
                   (this_atom_name == " C  ") || (this_atom_name == " O  ") ||
                   (this_atom_name == " H  "))
                  if ((env_atom_name == " N  ") || (env_atom_name == " CA ") ||
                      (env_atom_name == " C  ") || (env_atom_name == " O  ") ||
                      (env_atom_name == " H  "))
                     draw_flag = false;

               if ((this_atom_name == " O3'") || (this_atom_name == " C3'") ||
                   (this_atom_name == " P  ") || (this_atom_name == " OP1") ||
                   (this_atom_name == " OP2") || (this_atom_name == " O5'") ||
                   (this_atom_name == " C5'"))
                  if ((env_atom_name == " O3'") || (env_atom_name == " C3'") ||
                      (env_atom_name == " P  ") || (env_atom_name == " OP1") ||
                      (env_atom_name == " OP2") || (env_atom_name == " O5'") ||
                      (env_atom_name == " C5'"))
                     draw_flag = false;
            }
         }
      }
   }
   return draw_flag;
}



// This finds bonds between a residue and the protein (in SelAtom).
// It is used for the environment bonds box.
//
// It seems that this will draw bond to hydrogens, even if they are
// not being displayed.  I guess that we need to pass
// bonds_to_hydrogens flag.  And above.
//
Bond_lines_container::Bond_lines_container(const atom_selection_container_t &SelAtom,
                                           mmdb::PPAtom residue_atoms,
                                           int nResidueAtoms,
                                           float min_dist,
                                           float max_dist,
                                           bool draw_env_distances_to_hydrogens_flag,
                                           short int do_symm) {

   // std::cout << "Environment distances with symm" << std::endl;

   do_bonds_to_hydrogens = 1;  // added 20070629
   for_GL_solid_model_rendering = 1;
   have_dictionary = 0;
   init();
   n_atoms_in_atom_selection = SelAtom.n_selected_atoms;

   b_factor_scale = 1.0;
   if (bonds.size() == 0) {
      for (int i=0; i<13; i++) {  // 13 colours now in bond_colours
         Bond_lines a(i);
         bonds.push_back(a);
      }
   }

   if (min_dist> max_dist) {
      float tmp = max_dist;
      max_dist = min_dist;
      min_dist = tmp;
   }

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
   for (int iresatom=0; iresatom< nResidueAtoms; iresatom++) {
      mmdb::Atom *res_atom = residue_atoms[iresatom];
      coot::Cartesian res_atom_pos(res_atom->x, res_atom->y, res_atom->z);
      int model_number = residue_atoms[iresatom]->GetModelNum();

      molecule_extents_t mol_extents(SelAtom, max_dist);
      std::vector<std::pair<symm_trans_t, Cell_Translation> > boxes =
         mol_extents.which_boxes(res_atom_pos, SelAtom);

      for (unsigned int ibox=0; ibox<boxes.size(); ibox++) {
         mmdb::PPAtom translated = trans_sel(SelAtom, boxes[ibox]);

         for (int it=0; it<SelAtom.n_selected_atoms; it++) {

            int iat_1 = -1; // 20171224-PE FIXME. Maybe not needed (environment distances)
            int iat_2 = -1;

            mmdb::Atom *trans_atom = translated[it];
            coot::Cartesian symm_atom_pos(translated[it]->x,
                                          translated[it]->y,
                                          translated[it]->z);

            float d = coot::Cartesian::LineLength(symm_atom_pos, res_atom_pos);
             // std::cout << translated[it] << " d = " << d << std::endl;
            if (d < max_dist && d >= min_dist) {
               std::string ele1 = residue_atoms[iresatom]->element;
               std::string ele2 = translated[it]->element;

               if (draw_env_distances_to_hydrogens_flag ||
                   // ((ele1 != " H") && (ele2 != " H"))) {
                   ((!is_hydrogen(ele1)) && (! is_hydrogen(ele2)))) {

                  if (ele1 == " C")
                     addBond(0, symm_atom_pos, res_atom_pos, cc, model_number, iat_1, iat_2);
                  else {
                     if (ele2 == " C") {
                        addBond(0, symm_atom_pos, res_atom_pos, cc, model_number, iat_1, iat_2);
                     } else {
                        if (is_hydrogen(ele1) && is_hydrogen(ele2)) {
                           addBond(0, symm_atom_pos, res_atom_pos, cc, model_number, iat_1, iat_2);
                        } else {
                           coot::quick_protein_donor_acceptors pda;
                           coot::quick_protein_donor_acceptors::key k1(trans_atom->GetResName(), trans_atom->GetAtomName());
                           coot::quick_protein_donor_acceptors::key k2(  res_atom->GetResName(),   res_atom->GetAtomName());
                           // is-looked-up, is-H-bond
                           int colour_index = 1; // H-bond
                           std::pair<bool,bool> is_valid = pda.is_hydrogen_bond_by_types(k1,k2);
                           if (is_valid.first)
                              if (! is_valid.second)
                                 colour_index = 0;
                           addBond(colour_index, symm_atom_pos, res_atom_pos, cc, model_number, iat_1, iat_2);
                        }
                     }
                  }
               }
            }
         }
         // valgrind suggestion: 050803 - deleting translated[i] too.
         // c.f. pick usage/deletion of translation atoms
         for (int i=0; i<SelAtom.n_selected_atoms; i++)
            delete translated[i];
         delete [] translated;
      }
   }
   // std::cout << "returning from Bond_lines_container (symm env) constructor \n" ;
}


// What are the bonds to the symmetry-related copies of this molecule?
//
std::vector<Bond_lines_container::symmetry_atom_bond>
Bond_lines_container::find_intermolecular_symmetry(const atom_selection_container_t &SelAtom) const {

   std::vector<symmetry_atom_bond> sabv;

   int n_symm = SelAtom.mol->GetNumberOfSymOps();

   // std::cout << "in find_intermolecular_symmetry() n_symm is " << n_symm << std::endl;

   int shift_lim = 1;
   mmdb::mat44 my_matt;
   mmdb::realtype max_bond_dist = 2.25; // I guess

   for (int x_shift= -shift_lim; x_shift <= shift_lim; x_shift++) {
      for (int y_shift= -shift_lim; y_shift <= shift_lim; y_shift++) {
         for (int z_shift= -shift_lim; z_shift <= shift_lim; z_shift++) {
            for (int i_symm=0; i_symm < n_symm; i_symm++) {
               if (! (x_shift == 0 && y_shift == 0 && z_shift == 0 && i_symm == 0)) {
                  int i_status = SelAtom.mol->GetTMatrix(my_matt, i_symm, x_shift, y_shift, z_shift);

                  if (i_status == 0) {
                     // Happy

                     mmdb::Contact *contact = NULL;
                     int ncontacts = 0;
                     long i_contact_group = 1;

                     SelAtom.mol->SeekContacts(SelAtom.atom_selection, SelAtom.n_selected_atoms,
                                               SelAtom.atom_selection, SelAtom.n_selected_atoms,
                                               0.01, max_bond_dist,
                                               0,
                                               contact, ncontacts,
                                               0, &my_matt, i_contact_group);
                     if (ncontacts) {

                        symm_trans_t st(i_symm, x_shift, y_shift, z_shift);
                        for (int i=0; i< ncontacts; i++) {
                           mmdb::Atom *at_1 = SelAtom.atom_selection[contact[i].id1];
                           mmdb::Atom *at_2 = SelAtom.atom_selection[contact[i].id2];

                           std::string ele_1 = at_1->element;
                           std::string ele_2 = at_2->element;

                           if (ele_1 == " H" || ele_1 == "H") max_bond_dist -= 0.8;
                           if (ele_2 == " H" || ele_2 == "H") max_bond_dist -= 0.8;

                           st.symm_as_string = SelAtom.mol->GetSymOp(i_symm);

                           if (0) { //debug
                              int ierr = SelAtom.mol->GetTMatrix(my_matt, st.isym(),
                                                                 st.x(), st.y(), st.z());
                              if (! ierr) {
                                 mmdb::Atom t_atom2;
                                 t_atom2.Copy(at_2);
                                 t_atom2.Transform(my_matt);
                                 coot::Cartesian atom_1_pt(at_1->x, at_1->y, at_1->z);
                                 coot::Cartesian atom_2_pt(at_2->x, at_2->y, at_2->z);
                                 coot::Cartesian t_atom_2_pt(t_atom2.x, t_atom2.y, t_atom2.z);
                                 std::cout << "store at_1: " << atom_1_pt << " ";
                                 std::cout << "at_2: " << atom_2_pt << " ";
                                 std::cout << "t-at_2: " << t_atom_2_pt << " ";
                                 std::cout << st << " delta ";
                                 std::cout << (t_atom_2_pt - atom_1_pt).amplitude() << "\n";
                              }
                           }

                           Cell_Translation ct(0,0,0); // here for the ride only, ATM.
                           symmetry_atom_bond sab(at_1, at_2, st, ct);
                           sabv.push_back(sab);
                        }

                        delete [] contact;
                        contact = NULL;
                     }
                  } else {
                     std::cout << "unhappy call of GetTMatrix() " << std::endl;
                  }
               }
            }
         }
      }
   }
   std::cout << "found " << sabv.size() << " symmetry-atom-bonds" << std::endl;
   return sabv;
}

// This *looks* like a kludge - symm_trans is forced into a vector.
// But it is not a kludge. addSymmetry needs to take a vector, because
// it's needed for CA symm addition too.
//
//
graphical_bonds_container
Bond_lines_container::addSymmetry(const atom_selection_container_t &SelAtom,
                                  int imol,
                                  coot::Cartesian point,
                                  float symm_distance,
                                  const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
                                  short int symmetry_as_ca_flag,
                                  short int symmetry_whole_chain_flag) {

   graphical_bonds_container gbc;
   bool do_rama_markup = false;
   bool show_atoms_as_aniso_flag = false;
   do_sticks_for_waters = true;

   if (symmetry_as_ca_flag == 1) {
      gbc = addSymmetry_calphas(SelAtom, point, symm_distance, symm_trans);

   } else {

      if (symmetry_whole_chain_flag) {
         gbc = addSymmetry_whole_chain(SelAtom, imol, point,
                                       symm_distance, symm_trans);
      } else {

         mmdb::Contact *contact = NULL;
         int ncontacts;

         if (SelAtom.n_selected_atoms > 0) {
            mmdb::Atom *point_atom_p = new mmdb::Atom;
            point_atom_p->SetCoordinates(point.get_x(), point.get_y(),
                                         point.get_z(), 1.0, 99.9);

            for(unsigned int ii=0; ii<symm_trans.size(); ii++) { // i.e. boxes

               if (false)
                  std::cout << "-- " << ii << " " << symm_trans[ii].first << " "
                            << symm_trans[ii].second << std::endl;

               // now get an atom selection where all the atoms are moved
               // by symm_trans[ii]
               //
               // Create a selection trans_selection that is a copy of
               // SelAtom.atom_selection that has had symmetry
               // transformation symm_trans[ii] applied to it.
               //
               mmdb::PPAtom trans_selection = trans_sel(SelAtom, symm_trans[ii]); // heavyweight

               contact = NULL; // assigned next, deleted below.
               SelAtom.mol->SeekContacts(point_atom_p, trans_selection,
                                         SelAtom.n_selected_atoms,
                                         0.0, symm_distance,
                                         0,  // seqDist
                                         contact, ncontacts);

               if (ncontacts > 0 ) {

                  // the atom_selection of Contact_Sel is allocated.
                  // We give it below.
                  //
                  atom_selection_container_t Contact_Sel =
                     ContactSel(trans_selection, contact, ncontacts);
                  Contact_Sel.mol = SelAtom.mol;

                  int model_number = 0; // all models
                  construct_from_asc(Contact_Sel, imol, 0.01, 1.95,
                                     coot::COLOUR_BY_ATOM_TYPE, 1, model_number,
                                     show_atoms_as_aniso_flag, do_rama_markup);
                  gbc = make_graphical_symmetry_bonds();

                  // Now give back the atom_selection, (but not the atoms
                  // because they are created by trans_sel and are given
                  // back later.     giveback_1
                  //
                  delete [] Contact_Sel.atom_selection;

               } else {

                  // do a constructor for no symmetry bonds
                  // and reset graphical_symmetry_bonds.

                  no_symmetry_bonds();

               }

               // Tidy up.
               //
               if (trans_selection) {
                  for (int j=0; j<SelAtom.n_selected_atoms; j++)
                     // delete each of the atoms
                     delete trans_selection[j];
                  delete [] trans_selection;
               }

               delete [] contact;
            }
            delete point_atom_p;
         }
      }
   }
   return gbc;
}

// FYI: there is only one element to symm_trans, the is called from
// the addSymmetry_vector_symms wrapper
graphical_bonds_container
Bond_lines_container::addSymmetry_whole_chain(const atom_selection_container_t &SelAtom,
                                              int imol,
                                              const coot::Cartesian &point,
                                              float symm_distance,
                                              const std::vector <std::pair<symm_trans_t, Cell_Translation> > &symm_trans) {

   graphical_bonds_container gbc;
   mmdb::mat44 my_matt;
   mmdb::mat44 identity_matt;
   //   mmdb::Contact *contact = NULL;
   // int ncontacts;
   // long i_contact_group = 1;

   for (int im=0; im<4; im++) {
      for (int jm=0; jm<4; jm++) {
         identity_matt[im][jm] = 0.0;
      }
   }
   for (int im=0; im<4; im++)
      identity_matt[im][im] = 1.0;

   // There should be only one of these in the new (may 2004) system:
   for (unsigned int ii=0; ii<symm_trans.size(); ii++) {

      mmdb::mat44 mol_to_origin_matt;
      SelAtom.mol->GetTMatrix(mol_to_origin_matt, 0,
                              -symm_trans[ii].second.us,
                              -symm_trans[ii].second.vs,
                              -symm_trans[ii].second.ws);

      int ierr = SelAtom.mol->GetTMatrix(my_matt,
                                         symm_trans[ii].first.isym(),
                                         symm_trans[ii].first.x(),
                                         symm_trans[ii].first.y(),
                                         symm_trans[ii].first.z());

      if (ierr != 0) {
         std::cout << "!!!!!!!! something BAD with mmdb::CMMDBCryst.GetTMatrix"
              << std::endl;
      }

      int nmodels = SelAtom.mol->GetNumberOfModels();
      for (int imodel=1; imodel<=nmodels; imodel++) {
         // Here we want to get a selection of all atoms that have
         // been translated to this symm_trans, then calculate bonds
         // on them.
         mmdb::Atom **transsel = new mmdb::PAtom[SelAtom.n_selected_atoms];
         for (int i=0; i<SelAtom.n_selected_atoms; i++) {
            transsel[i] = new mmdb::Atom;
            transsel[i]->Copy(SelAtom.atom_selection[i]);
            transsel[i]->residue = SelAtom.atom_selection[i]->residue;
            transsel[i]->Transform(mol_to_origin_matt);
            transsel[i]->Transform(my_matt);
         }

         // So now we have transsel.  We need to calculate bonds from
         // these:

         short int atom_colour_type = coot::COLOUR_BY_ATOM_TYPE;
         construct_from_atom_selection(SelAtom,
                                       transsel, SelAtom.n_selected_atoms,
                                       transsel, SelAtom.n_selected_atoms,
                                       imol,
                                       0.1, 1.8,
                                       atom_colour_type,
                                       0, 1, SelAtom.UDDAtomIndexHandle);

         for (int i=0; i<SelAtom.n_selected_atoms; i++)
            delete transsel[i]; // delete the atoms
         delete [] transsel;  // delete the atom list
      }
   }
   if (symm_trans.size() > 0)
      gbc = make_graphical_symmetry_bonds();
   return gbc;
}


//
graphical_bonds_container
Bond_lines_container::addSymmetry_with_mmdb(const atom_selection_container_t &SelAtom,
                                            int imol,
                                            coot::Cartesian point,
                                            float symm_distance,
                                            const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
                                            short int symmetry_as_ca_flag) {

   graphical_bonds_container gbc;
   bool do_rama_markup = false;
   bool show_atoms_as_aniso_flag = false;

   if (symmetry_as_ca_flag == 1) {
      gbc = addSymmetry_calphas(SelAtom, point, symm_distance, symm_trans);

   } else {

      mmdb::Contact *contact;
      int ncontacts;

      if (SelAtom.n_selected_atoms > 0) {
         mmdb::PAtom point_atom_p = new mmdb::Atom;
         point_atom_p->SetCoordinates(point.get_x(), point.get_y(),
                                      point.get_z(), 1.0, 99.9);

         for(unsigned int ii=0; ii<symm_trans.size(); ii++) { // i.e. boxes

            //          cout << "---------------------------- "
            //               << symm_trans[ii]
            //               << " ---------------------------- "
            //               << endl;

            // now get an atom selection where all the atoms are moved
            // by symm_trans[ii]
            //
            // Create a selection trans_selection that is a copy of
            // SelAtom.atom_selection that has had symmetry
            // transformation symm_trans[ii] applied to it.
            //
            mmdb::PPAtom trans_selection = trans_sel(SelAtom, symm_trans[ii]);

            contact = NULL;
            SelAtom.mol->SeekContacts(point_atom_p, trans_selection,
                                      SelAtom.n_selected_atoms,
                                      0.0, symm_distance,
                                      0,  // seqDist
                                      contact, ncontacts);

            // cout << symm_trans[ii] << " had " << ncontacts
            // << " symmetry atom within "
            // << symm_distance << "A radius." << endl;

            if (ncontacts > 0 ) {

               // the atom_selection of Contact_Sel is allocated,
               // where do we give it up?
               //
               atom_selection_container_t Contact_Sel =
                  ContactSel(trans_selection, contact, ncontacts);
               Contact_Sel.mol = SelAtom.mol;

               int model_number = 0; // all models
               construct_from_asc(Contact_Sel, imol, 0.01, 1.95,
                                  coot::COLOUR_BY_ATOM_TYPE, 1, model_number,
                                  show_atoms_as_aniso_flag, do_rama_markup);
               gbc = make_graphical_symmetry_bonds();

               // Now give back the atom_selection, (but not the atoms
               // because they are created by trans_sel and are given
               // back later.     giveback_1
               //
               delete [] Contact_Sel.atom_selection;

            } else {

               // do a constructor for no symmetry bonds
               // and reset graphical_symmetry_bonds.

               no_symmetry_bonds();

            }

            // Tidy up.
            //
            if (trans_selection) {
               for (int j=0; j<SelAtom.n_selected_atoms; j++)
                  // delete each of the atoms
                  delete trans_selection[j];
               delete [] trans_selection;
            }

            delete [] contact;
         }
         delete point_atom_p;
      }
   }
   return gbc;
}


std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > >
Bond_lines_container::addSymmetry_vector_symms(const atom_selection_container_t &SelAtom,
                                               int imol,
                                               coot::Cartesian point,
                                               float symm_distance,
                                               const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
                                               short int symmetry_as_ca_flag,
                                               short int symmetry_whole_chain_flag,
                                               short int draw_hydrogens_flag,
                                               bool do_intermolecular_symmetry_bonds) {

   std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > > r;
   do_bonds_to_hydrogens = draw_hydrogens_flag;
   std::vector <symmetry_atom_bond> sabv;
   if (do_intermolecular_symmetry_bonds)
      // heavyweight.
      sabv = find_intermolecular_symmetry(SelAtom);

   do_sticks_for_waters = true; // does this work?
   for (unsigned int i=0; i<symm_trans.size(); i++) {
      const std::pair<symm_trans_t, Cell_Translation>  &st = symm_trans[i];
      std::vector<std::pair<symm_trans_t, Cell_Translation> > this_symm_trans_vec;
      this_symm_trans_vec.push_back(st);

      r.push_back(std::pair<graphical_bonds_container,
                  std::pair<symm_trans_t, Cell_Translation > > (addSymmetry(SelAtom, imol,
                                                                            point, symm_distance,
                                                                            this_symm_trans_vec,
                                                                            symmetry_as_ca_flag,
                                                                            symmetry_whole_chain_flag),
                                                                st));

      if (do_intermolecular_symmetry_bonds) {
         graphical_bonds_container int_sym_gbc = intermolecular_symmetry_graphical_bonds(SelAtom.mol, sabv, st);
         std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation > > p(int_sym_gbc, st);
         r.push_back(p);
      }
   }

   return r;
}

graphical_bonds_container
Bond_lines_container::intermolecular_symmetry_graphical_bonds(mmdb::Manager *mol,
                                                              const std::vector <Bond_lines_container::symmetry_atom_bond> &sabv,
                                                              const std::pair<symm_trans_t, Cell_Translation> &symm_trans) {

   graphical_bonds_container gbc;
   mmdb::mat44 my_matt;
   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

   if (0)
      std::cout << "in intermolecular_symmetry_graphical_bonds() running through " << sabv.size()
                << " symmetry atom bonds " << std::endl;

   for (unsigned int i=0; i<sabv.size(); i++) {

      coot::Cartesian atom_1_pt(sabv[i].at_1->x, sabv[i].at_1->y, sabv[i].at_1->z);

      int ierr = sabv[i].GetTMatrix(mol, &my_matt); // check that my_matt gets changed to something sensible
      int model_number = sabv[i].at_1->GetModelNum();

      if (! ierr) {

         coot::Cartesian atom_2_pt(sabv[i].at_2->x, sabv[i].at_2->y, sabv[i].at_2->z);
         mmdb::Atom t_atom2;
         t_atom2.Copy(sabv[i].at_2);
         t_atom2.Transform(my_matt);
         coot::Cartesian t_atom_2_pt(t_atom2.x, t_atom2.y, t_atom2.z);

         if (0)
            std::cout << "gbc int-symm bond: atom_1: " << atom_1_pt << " atom_2: "
                      << atom_2_pt << " t_atom2 " << t_atom_2_pt << " "
                      << symm_trans.first << " " << " delta "
                      << (t_atom_2_pt - atom_1_pt).amplitude() << "\n";

         int iat_1 = -1; // 20171224-PE FIXME by lookup. Maybe not needed.
         int iat_2 = -1;
         addBond(0, atom_1_pt, t_atom_2_pt, cc, model_number, iat_1, iat_2);

      } else {
         std::cout << "intermolecular_symmetry_graphical_bonds GetTMatrix() problem."
                   << std::endl;
      }
   }

   gbc = make_graphical_symmetry_bonds();
   return gbc;
}



graphical_bonds_container
Bond_lines_container::addSymmetry_calphas(const atom_selection_container_t &SelAtom,
                                          const coot::Cartesian &point,
                                          float symm_distance,
                                          const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans) {

   graphical_bonds_container gbc;
   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
   // std::cout << "DEBUG:: There are " << symm_trans.size() << "
   // symm_trans units" << std::endl;
   mmdb::mat44 my_matt;
   mmdb::Atom t_atom1;
   mmdb::Atom t_atom2;

   for (unsigned int ii=0; ii<symm_trans.size(); ii++) {
      int err = SelAtom.mol->GetTMatrix(my_matt, symm_trans[ii].first.isym(), symm_trans[ii].first.x(),
                                        symm_trans[ii].first.y(), symm_trans[ii].first.z());

      mmdb::mat44 mol_to_origin_matt;
      SelAtom.mol->GetTMatrix(mol_to_origin_matt, 0,
                              -symm_trans[ii].second.us,
                              -symm_trans[ii].second.vs,
                              -symm_trans[ii].second.ws);

      if (err != 0) {
         std::cout << "!!!!!!!!!!!!!! something BAD with mmdb::CMMDBCryst.GetTMatrix"
                   << std::endl;
      }

      int nmodels = SelAtom.mol->GetNumberOfModels();
      for (int imodel=1; imodel<=nmodels; imodel++) {
         int nchains = SelAtom.mol->GetNumberOfChains(imodel);
         for (int ichain=0; ichain<nchains; ichain++) {
            mmdb::PChain chain = SelAtom.mol->GetChain(imodel, ichain);
            int nres = chain->GetNumberOfResidues();
            for (int ires=0; ires<(nres-1); ires++) {
               // std::cout << "residue " << ires << " to  " << ires+1 << std::endl;
               // one residue to the next...
               mmdb::PResidue res1 = chain->GetResidue(ires);
               mmdb::PResidue res2 = chain->GetResidue(ires+1);

               // Search through the atoms of each residue looking for CAs:
               std::vector<mmdb::Atom *> ca_this;
               std::vector<mmdb::Atom *> ca_next;

               mmdb::PPAtom residue_atoms;
               int n_residue_atoms;
               res1->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  if (residue_atoms[iat]) {
                     if (std::string(residue_atoms[iat]->name) == " CA " ||
                         std::string(residue_atoms[iat]->name) == " P  ") {
                        ca_this.push_back(residue_atoms[iat]);
                     }
                  }
               }
               res2->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  if (residue_atoms[iat]) {
                     if (std::string(residue_atoms[iat]->name) == " CA " ||
                         std::string(residue_atoms[iat]->name) == " P  ") {
                        ca_next.push_back(residue_atoms[iat]);
                     }
                  }
               }
               // so now we have the 2 vectors filled with CA atoms of
               // this and the next residue.  What are the connections?
               //
               if (ca_this.size() > 0) {
                  if (ca_next.size() > 0) {
                     for (unsigned int iat=0; iat<ca_this.size(); iat++) {
                        std::string altconf1 = ca_this[iat]->altLoc;
                        for (unsigned int jat=0; jat<ca_next.size(); jat++) {
                           std::string altconf2 = ca_next[jat]->altLoc;

                           if ((altconf1 == altconf2) ||
                               (altconf1 == "") ||
                               (altconf2 == "")) {
                              coot::Cartesian ca_1(ca_this[iat]->x,
                                                   ca_this[iat]->y,
                                                   ca_this[iat]->z);
                              coot::Cartesian ca_2(ca_next[jat]->x,
                                                   ca_next[jat]->y,
                                                   ca_next[jat]->z);

                              double len = (ca_1 - ca_2).amplitude();
                              // CA-CA or P-P
                              if (((len < 4.7) && (len > 2.4)) ||
                                  ((len<8) && (len>5))) {

                                 t_atom1.Copy(ca_this[iat]);
                                 t_atom2.Copy(ca_next[jat]);
                                 t_atom1.Transform(mol_to_origin_matt);
                                 t_atom2.Transform(mol_to_origin_matt);
                                 t_atom1.Transform(my_matt);
                                 t_atom2.Transform(my_matt);

                                 coot::Cartesian atom1(t_atom1.x, t_atom1.y, t_atom1.z);
                                 coot::Cartesian atom2(t_atom2.x, t_atom2.y, t_atom2.z);

                                 int iat_1 = -1; // 20171224-PE FIXME maybe lookup
                                 int iat_2 = -1;

                                 addBond(0, atom1, atom2, cc, imodel, iat_1, iat_2); // these are C atoms.
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

   if (symm_trans.size() > 0)
      gbc = make_graphical_symmetry_bonds();

   return gbc;
}

std::vector<std::pair<graphical_bonds_container, symm_trans_t> >
Bond_lines_container::add_NCS(const atom_selection_container_t &SelAtom,
                              int imol,
                              coot::Cartesian point,
                              float symm_distance,
                              std::vector<std::pair<coot::coot_mat44, symm_trans_t> > &strict_ncs_mats,
                              short int symmetry_as_ca_flag,
                              short int symmetry_whole_chain_flag) {

   if (false) {
      std::cout << "in add_NCS() " << strict_ncs_mats.size() << std::endl;
      for(unsigned int ii=0; ii<strict_ncs_mats.size(); ii++) { // i.e. boxes
         std::cout << "-- " << ii << " "
                   << strict_ncs_mats[ii].second << std::endl;
      }
   }

   std::vector<std::pair<graphical_bonds_container, symm_trans_t> > r;
   for (unsigned int i=0; i<strict_ncs_mats.size(); i++) {
      std::pair<graphical_bonds_container, symm_trans_t>
         p(add_NCS_molecule(SelAtom, imol,
                            point, symm_distance,
                            strict_ncs_mats[i],
                            symmetry_as_ca_flag,
                            symmetry_whole_chain_flag),
           strict_ncs_mats[i].second);
      r.push_back(p);
   }
   return r;
}

graphical_bonds_container
Bond_lines_container::add_NCS_molecule(const atom_selection_container_t &SelAtom,
                                       int imol,
                                       const coot::Cartesian &point,
                                       float symm_distance,
                                       const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat,
                                       short int symmetry_as_ca_flag,
                                       short int symmetry_whole_chain_flag) {

   graphical_bonds_container gbc;
   mmdb::mat44 m;

   if (symmetry_as_ca_flag) {
      gbc = add_NCS_molecule_calphas(SelAtom, point, symm_distance, strict_ncs_mat);
   } else {
      if (symmetry_whole_chain_flag) {
         gbc = add_NCS_molecule_whole_chain(SelAtom, imol, point, symm_distance, strict_ncs_mat);
      } else {

         m[0][0] = strict_ncs_mat.first.m[0].v4[0];
         m[0][1] = strict_ncs_mat.first.m[0].v4[1];
         m[0][2] = strict_ncs_mat.first.m[0].v4[2];
         m[1][0] = strict_ncs_mat.first.m[1].v4[0];
         m[1][1] = strict_ncs_mat.first.m[1].v4[1];
         m[1][2] = strict_ncs_mat.first.m[1].v4[2];
         m[2][0] = strict_ncs_mat.first.m[2].v4[0];
         m[2][1] = strict_ncs_mat.first.m[2].v4[1];
         m[2][2] = strict_ncs_mat.first.m[2].v4[2];
         m[0][3] = strict_ncs_mat.first.m[0].v4[3];  // t
         m[1][3] = strict_ncs_mat.first.m[1].v4[3];  // t
         m[2][3] = strict_ncs_mat.first.m[2].v4[3];  // t
         m[3][0] = strict_ncs_mat.first.m[3].v4[0];
         m[3][1] = strict_ncs_mat.first.m[3].v4[1];
         m[3][2] = strict_ncs_mat.first.m[3].v4[2];
         m[3][3] = strict_ncs_mat.first.m[3].v4[3];

      }
   }
   return gbc;
}

graphical_bonds_container
Bond_lines_container::add_NCS_molecule_whole_chain(const atom_selection_container_t &SelAtom,
                                                   int imol,
                                                   const coot::Cartesian &point,
                                                   float symm_distance,
                                                   const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat) {

   graphical_bonds_container gbc;
   mmdb::mat44 my_matt;

   my_matt[0][0] = strict_ncs_mat.first.m[0].v4[0];
   my_matt[0][1] = strict_ncs_mat.first.m[0].v4[1];
   my_matt[0][2] = strict_ncs_mat.first.m[0].v4[2];
   my_matt[1][0] = strict_ncs_mat.first.m[1].v4[0];
   my_matt[1][1] = strict_ncs_mat.first.m[1].v4[1];
   my_matt[1][2] = strict_ncs_mat.first.m[1].v4[2];
   my_matt[2][0] = strict_ncs_mat.first.m[2].v4[0];
   my_matt[2][1] = strict_ncs_mat.first.m[2].v4[1];
   my_matt[2][2] = strict_ncs_mat.first.m[2].v4[2];
   my_matt[0][3] = strict_ncs_mat.first.m[0].v4[3];  // t
   my_matt[1][3] = strict_ncs_mat.first.m[1].v4[3];  // t
   my_matt[2][3] = strict_ncs_mat.first.m[2].v4[3];  // t
   my_matt[3][0] = strict_ncs_mat.first.m[3].v4[0];
   my_matt[3][1] = strict_ncs_mat.first.m[3].v4[1];
   my_matt[3][2] = strict_ncs_mat.first.m[3].v4[2];
   my_matt[3][3] = strict_ncs_mat.first.m[3].v4[3];

   int nmodels = SelAtom.mol->GetNumberOfModels();
   for (int imodel=1; imodel<=nmodels; imodel++) {
      // Here we want to get a selection of all atoms that have
      // been translated to this symm_trans, then calculate bonds
      // on them.
      mmdb::Atom **transsel = new mmdb::PAtom[SelAtom.n_selected_atoms];
      for (int i=0; i<SelAtom.n_selected_atoms; i++) {
         transsel[i] = new mmdb::Atom;
         transsel[i]->Copy(SelAtom.atom_selection[i]);
         transsel[i]->residue = SelAtom.atom_selection[i]->residue;
         transsel[i]->Transform(my_matt);
      }
      // So now we have transsel.  We need to calculate bonds from
      // these:

      short int atom_colour_type = coot::COLOUR_BY_ATOM_TYPE;
      construct_from_atom_selection(SelAtom,
                                    transsel, SelAtom.n_selected_atoms,
                                    transsel, SelAtom.n_selected_atoms,
                                    imol,
                                    0.1, 1.8,
                                    atom_colour_type,
                                    0, 1, SelAtom.UDDAtomIndexHandle);
//       for (int i=0; i<SelAtom.n_selected_atoms; i++) {
//          std::cout << "atom sel[" << i << "] " << SelAtom.atom_selection[i] << "\n";
//          std::cout << "transsel[" << i << "] " << transsel[i] << "\n";
//       }
      gbc = make_graphical_symmetry_bonds();
      for (int i=0; i<SelAtom.n_selected_atoms; i++)
         delete transsel[i];
      delete [] transsel;

   }
   return gbc;
}

graphical_bonds_container
Bond_lines_container::add_NCS_molecule_calphas(const atom_selection_container_t &SelAtom,
                                               const coot::Cartesian &point,
                                               float symm_distance,
                                               const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat) {

   graphical_bonds_container gbc;

   return gbc;

}

// We have a set of contacts in contact for the atom selection trans_sel.
//
// We want to return an atom selection that is just the contact atoms.
// We also need the number of atom (which is the the same as the
// number of atom contacts), so we fill in the 2 parts of an
// atom_selection_container_t.
//
atom_selection_container_t
Bond_lines_container::ContactSel(mmdb::PPAtom trans_sel,
                                 mmdb::Contact *contact, int ncontacts) const {

   // We need to sort the contacts.
   //
   SortContacts(contact, ncontacts, mmdb::CNSORT_2INC);
   //
   int id2, last_id2 = -1;
   int n_contact_atoms = 0;
   //
   // keep in mind that there is a contact id2->id1 as well as id1->id2.

   atom_selection_container_t TransSel;

   // You might ask: where does this get deleted?  It does get deleted: giveback_1
   //
   TransSel.atom_selection = new mmdb::PAtom[ncontacts]; //

   for (int icon=0; icon<ncontacts; icon++) {
      // use sorted contacts
      id2 = contact[icon].id2;
      TransSel.atom_selection[n_contact_atoms] = trans_sel[id2];
      n_contact_atoms++;
   }

   TransSel.n_selected_atoms = n_contact_atoms;
   if (n_contact_atoms > ncontacts){
      std::cout << "disaster n_contact_atoms (" << n_contact_atoms
                << ") > ncontacts (" << ncontacts << ")" << std::endl;
   }
   return TransSel;

}


// See above comments re const mmdb::CMMDBCryst &.
//
mmdb::PPAtom
Bond_lines_container::trans_sel(atom_selection_container_t AtomSel,
                                const std::pair<symm_trans_t, Cell_Translation> &symm_trans) const {
   mmdb::mat44 my_matt;
   //
   // modify my_matt;
   //
   //
   int err = AtomSel.mol->GetTMatrix(my_matt, symm_trans.first.isym(), symm_trans.first.x(),
                                     symm_trans.first.y(), symm_trans.first.z());

   if (err != 0) {
      std::cout << "!!!!!!!!!!!!!! something BAD with mmdb::CMMDBCryst.GetTMatrix"
                << std::endl;
   }

   // cout << "using symm_trans: " << symm_trans << endl;

   //cout << "applying my_matt in trans_sel: " << endl;
//    cout << my_matt[0][0] << " "  << my_matt[0][1] << " "
//         << my_matt[0][2] << " "  << my_matt[0][3] << " "  << endl
//         << my_matt[1][0] << " "  << my_matt[1][1] << " "
//         << my_matt[1][2] << " "  << my_matt[1][3] << " "  << endl
//         << my_matt[2][0] << " "  << my_matt[2][1] << " "
//         << my_matt[2][2] << " "  << my_matt[2][3] << " "  << endl
//         << my_matt[3][0] << " "  << my_matt[3][1] << " "
//         << my_matt[3][2] << " "  << my_matt[3][3] << " "  << endl;


   mmdb::mat44 mol_to_origin_matt;
   AtomSel.mol->GetTMatrix(mol_to_origin_matt, 0,
                           -symm_trans.second.us,
                           -symm_trans.second.vs,
                           -symm_trans.second.ws);

   // Whoah!  Big alloc!  Given back in addSymmetry
   //
   //cout << "allocating translation space for " << AtomSel.n_selected_atoms
   //        << " atoms" << endl;
   //
   //
   mmdb::PPAtom trans_selection = new mmdb::PAtom[AtomSel.n_selected_atoms];
   for (int ii=0; ii<AtomSel.n_selected_atoms; ii++) {

      trans_selection[ii] = new mmdb::Atom;

//       trans_selection[ii]->SetCoordinates(AtomSel.atom_selection[ii]->x,
//                                           AtomSel.atom_selection[ii]->y,
//                                           AtomSel.atom_selection[ii]->z,
//                                           1.0, 99.9);
//       trans_selection[ii]->SetAtomName(   AtomSel.atom_selection[ii]->name);
//       trans_selection[ii]->SetElementName(AtomSel.atom_selection[ii]->element);
//       trans_selection[ii]->SetResidue(    AtomSel.atom_selection[ii]->GetResidue());
//       trans_selection[ii]->tempFactor =   AtomSel.atom_selection[ii]->tempFactor;



      // can't set seqNum:
      //trans_selection[ii]->seqNum     =   AtomSel.atom_selection[ii]->GetSeqNum();
      // I don't know why, but this doesn't work.
      // trans_selection[ii]->SetChain(      AtomSel.atom_selection[ii]->GetChain());

      trans_selection[ii]->Copy(AtomSel.atom_selection[ii]);
      trans_selection[ii]->residue = AtomSel.atom_selection[ii]->residue;
      trans_selection[ii]->Transform(mol_to_origin_matt);
      trans_selection[ii]->Transform(my_matt);
   }
   return trans_selection;
}

unsigned int
graphical_bonds_container::n_bonds() const { // count them up

   int count = 0;
   for (int idx_col=0; idx_col<num_colours; idx_col++) {
      count += bonds_[idx_col].num_lines;
   }
   return count;

}

unsigned int
graphical_bonds_container::n_atoms() const {

   int count = n_atom_centres_;
   return count;
}

void
graphical_bonds_container::debug() const {
   std::cout << "This graphical_bonds_container has " << n_bonds() << " bonds and " << n_atoms()
             << " atoms." << std::endl;
}


// This function used by skeletonization (and nothing else)
void
graphical_bonds_container::add_colour(const std::vector<graphics_line_t> &a) {

   graphical_bonds_lines_list<graphics_line_t> *new_bonds_ =
         new graphical_bonds_lines_list<graphics_line_t>[num_colours+1];
      if ( bonds_ != NULL ) {
         for (int i = 0; i < num_colours; i++ ) new_bonds_[i] = bonds_[i];
         delete[] bonds_;
      }
      bonds_ = new_bonds_;
      // bonds_[num_colours].pair_list = new coot::CartesianPair[(a.size())];
      bonds_[num_colours].pair_list = new graphics_line_t[(a.size())];
      bonds_[num_colours].num_lines = a.size();

      // copy over
      for(unsigned int i=0; i<a.size(); i++) {
         bonds_[num_colours].pair_list[i] = a[i];
      }
      num_colours++;

      symmetry_bonds_ = NULL;
      symmetry_has_been_created = 0;
   }


void
Bond_lines_container::do_disulphide_bonds(atom_selection_container_t SelAtom, int imodel) {

   do_disulphide_bonds_by_distance(SelAtom, imodel);
}


void
Bond_lines_container::do_disulphide_bonds_by_header(atom_selection_container_t SelAtom,
                                                    int imodel) {

   // How do I see the SSBond records?
}

void
Bond_lines_container::do_disulphide_bonds_by_distance(atom_selection_container_t SelAtom,
                                                      int imodel) {

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

   // Lets make a new selection using mol.  We step back and sidestep
   // to a different atom selection.
   //
   mmdb::PPAtom Sulfur_selection;
   int n_sulfurs;
   mmdb::mat44 my_matt;
   mmdb::Contact *contact = NULL;
   int ncontacts = 0; // initially no S contacts.
   long i_contact_group = 1;
   int col;

   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;

   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   // I could pass this I suppose - that may be quicker, if this is a problem.
   int udd_atom_index_handle = SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index"); // set in make_asc

   int udd_user_defined_atom_colour_index_handle = SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   int selHnd2 = SelAtom.mol->NewSelection();

   // model 1
   // Note that now we force the resname to be CYS and atom name to be SG
   SelAtom.mol->SelectAtoms(selHnd2, imodel,"*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                            "CYS"," SG ","S","*" );

   SelAtom.mol->GetSelIndex(selHnd2, Sulfur_selection, n_sulfurs);

   float max_SS_bond_length = 2.4;

   if (n_sulfurs > 0) {
      SelAtom.mol->SeekContacts(Sulfur_selection, n_sulfurs,
                                Sulfur_selection, n_sulfurs,
                                0.01, max_SS_bond_length, // min, max dist.
                                0,         // in same res also.
                                contact, ncontacts,
                                0, &my_matt, i_contact_group);
   }

//    if (verbose_reporting)
//       cout << "Found " << ncontacts/2 << " disulfides" << endl;

   if (ncontacts > 0) {
      for (int i=0; i< ncontacts; i++) {

         if ( contact[i].id2 > contact[i].id1 ) {

            int iat_1 = -1;
            int iat_2 = -1;
            Sulfur_selection[contact[i].id1]->GetUDData(udd_atom_index_handle, iat_1);
            Sulfur_selection[contact[i].id2]->GetUDData(udd_atom_index_handle, iat_2);

            std::string aloc_1(Sulfur_selection[ contact[i].id1 ]->altLoc);
            std::string aloc_2(Sulfur_selection[ contact[i].id2 ]->altLoc);

            if ( (aloc_1=="") || (aloc_2=="") || (aloc_1==aloc_2) ) {
               coot::Cartesian atom_1(Sulfur_selection[ contact[i].id1 ]->x,
                                Sulfur_selection[ contact[i].id1 ]->y,
                                Sulfur_selection[ contact[i].id1 ]->z);

               coot::Cartesian atom_2(Sulfur_selection[ contact[i].id2 ]->x,
                                Sulfur_selection[ contact[i].id2 ]->y,
                                Sulfur_selection[ contact[i].id2 ]->z);

               col = atom_colour(Sulfur_selection[ contact[i].id1 ],
                                 coot::DISULFIDE_COLOUR, udd_user_defined_atom_colour_index_handle);

               if (! ((Sulfur_selection[ contact[i].id1 ]->GetSeqNum() ==
                       Sulfur_selection[ contact[i].id2 ]->GetSeqNum()) &&
                      (Sulfur_selection[ contact[i].id1 ]->GetChainID() ==
                       Sulfur_selection[ contact[i].id2 ]->GetChainID()))) {
                  int model_number = Sulfur_selection[ contact[i].id1 ]->GetModelNum();

                  // only add this bond if the atom is not already linked to something
                  // A Zn for example

                  bool is_linked = false;
                  mmdb::Model *model_p = SelAtom.mol->GetModel(model_number);
                  int n_links = model_p->GetNumberOfLinks();
                  if (n_links > 0) {
                     coot::atom_spec_t SS_atom_1_spec(Sulfur_selection[contact[i].id1]);
                     coot::atom_spec_t SS_atom_2_spec(Sulfur_selection[contact[i].id2]);
                     for (int i_link=1; i_link<=n_links; i_link++) {
                        mmdb::PLink link = model_p->GetLink(i_link);
                        std::pair<coot::atom_spec_t, coot::atom_spec_t> link_atom_specs = coot::link_atoms(link, model_p);
                        if (link_atom_specs.first  == SS_atom_1_spec) is_linked = true;
                        if (link_atom_specs.second == SS_atom_1_spec) is_linked = true;
                        if (link_atom_specs.first  == SS_atom_2_spec) is_linked = true;
                        if (link_atom_specs.second == SS_atom_2_spec) is_linked = true;
                        if (is_linked) break;
                     }
                  }

                  if (! is_linked)
                     addBond(col, atom_1, atom_2, cc, model_number, iat_1, iat_2);
               }
            }
         }
      }
      delete [] contact;
   }
   SelAtom.mol->DeleteSelection(selHnd2);
}

// reset bonds (delete if necessary) any graphics bonds.
//
void Bond_lines_container::no_symmetry_bonds() {

}


Bond_lines_container::Bond_lines_container(symm_keys key) {

   do_bonds_to_hydrogens = 1;  // added 20070629
   b_factor_scale = 1.0;
   for_GL_solid_model_rendering = 1;
   init();
   have_dictionary = 0;

   if (key == NO_SYMMETRY_BONDS) {

      no_symmetry_bonds();

   } else {
      std::cout << "Bond_lines_container::Bond_lines_container(symm_keys key)"
                << " no such key as " << key << std::endl;
   }
}


void Bond_lines_container::check_static() const {

   graphical_bonds_container pot;

   std::cout << "check: num_colours:"     << pot.num_colours << std::endl;
   std::cout << "check: bonds:"           << pot.bonds_ << std::endl;

}

graphical_bonds_container
Bond_lines_container::make_graphical_bonds() const {
   return make_graphical_bonds_with_thinning_flag(true); // allow simple-minded thinning
}


graphical_bonds_container
Bond_lines_container::make_graphical_bonds_no_thinning() const {
   return make_graphical_bonds_with_thinning_flag(false); // no thinning
}

graphical_bonds_container
Bond_lines_container::make_graphical_bonds_with_thinning_flag(bool do_thinning_flag) const {

   graphical_bonds_container box;

   int n_bond_colours = bonds.size();
   box.num_colours = bonds.size();
   box.bonds_ = new graphical_bonds_lines_list<graphics_line_t>[n_bond_colours];

   std::map<mmdb::Residue *, int> residue_index_map;
   mmdb::Residue *null_residue = 0; // for compiler to get the correct type in the next line
   residue_index_map[null_residue] = -1;

   // std::cout << ":::::::::::::::::::::::: n_bond_colours " << n_bond_colours << std::endl;

   for (int idx_col=0; idx_col<n_bond_colours; idx_col++) {

      box.bonds_[idx_col].num_lines = bonds[idx_col].size();
      box.bonds_[idx_col].pair_list = new graphics_line_t[bonds[idx_col].size()];
      for (int j=0; j<bonds[idx_col].size(); j++) {
         box.bonds_[idx_col].pair_list[j] = bonds[idx_col][j];
      }
      if (do_thinning_flag) {
         if (idx_col == HYDROGEN_GREY_BOND)
            box.bonds_[idx_col].thin_lines_flag = 1;
         if (idx_col == DEUTERIUM_PINK)
            box.bonds_[idx_col].thin_lines_flag = 1;
      }
   }
   box.add_zero_occ_spots(zero_occ_spots);
   box.add_deuterium_spots(deuterium_spots);
   // box.add_ramachandran_goodness_spots(ramachandran_goodness_spots); not in this function (I guess)

   box.add_rotamer_goodness_markup(dodecs);

   box.add_atom_centres(atom_centres, atom_centres_colour);
   box.rings = rings;
   box.add_bad_CA_CA_dist_spots(bad_CA_CA_dist_spots);
   box.add_cis_peptide_markup(cis_peptide_quads);
   return box;
}

graphical_bonds_container
Bond_lines_container::make_graphical_bonds(const ramachandrans_container_t &rc,
                                           bool do_rama_markup,
                                           bool do_rotamer_markup) const {


   // this one for intermediate atoms

   graphical_bonds_container box;
   bool thinning_flag = true;

   int ibs = bonds.size();
   box.num_colours = bonds.size();
   box.bonds_ = new graphical_bonds_lines_list<graphics_line_t>[ibs];

   // i is the colour index
   for (int i=0; i<ibs; i++) {

      box.bonds_[i].num_lines = bonds[i].size();
      // box.bonds_[i].pair_list = new coot::CartesianPair[bonds[i].size()];
      box.bonds_[i].pair_list = new graphics_line_t[bonds[i].size()];
      for (int j=0; j<bonds[i].size(); j++)
         box.bonds_[i].pair_list[j] = bonds[i][j];
      if (thinning_flag) {
         if (i == HYDROGEN_GREY_BOND)
            box.bonds_[i].thin_lines_flag = 1;
         if (i == DEUTERIUM_PINK) {
            box.bonds_[i].thin_lines_flag = 1;
         }
      }
   }
   box.add_zero_occ_spots(zero_occ_spots);
   box.add_deuterium_spots(deuterium_spots);
   if (do_rama_markup)
      box.add_ramachandran_goodness_spots(ramachandran_goodness_spots, rc);
   if (do_rotamer_markup)
      box.add_rotamer_goodness_markup(dodecs);
   box.add_atom_centres(atom_centres, atom_centres_colour);
   box.add_cis_peptide_markup(cis_peptide_quads);
   box.rings = rings;
   return box;
}

graphical_bonds_container
Bond_lines_container::make_graphical_symmetry_bonds() const {

   graphical_bonds_container box;
   box.num_colours = bonds.size();
   box.bonds_ = NULL;
   // box.bonds_ = new graphical_bonds_lines_list[bonds.size()]; // surely not!
   box.symmetry_has_been_created = 1;

   // This block can never happen!  box is created here and the
   // constructor sets symmetry_bonds_ to NULL.
   //
   // delete the old bonds
   //
   if (box.symmetry_bonds_ != NULL) {
      //
      //cout << "deleting symmetry bonds (bonds.size() is " << bonds.size()
      //           << ")" << endl;
      for (int i = 0; i < box.num_colours; i++) {
         if (box.symmetry_bonds_[i].num_lines > 0) {
            //cout << "box.symmetry_bonds_[" << i << "].num_lines is "
            //<< box.symmetry_bonds_[i].num_lines << endl;
            delete [] box.symmetry_bonds_[i].pair_list;
         }
      }
   }

   box.symmetry_bonds_ = new graphical_bonds_lines_list<graphics_line_t>[bonds.size()];
//    std::cout << "allocating symmetry bonds " << box.symmetry_bonds_
//              << " for " << bonds.size() << " symmetry bonds" << std::endl;

   int num_lines;
   for (int i = 0; i < box.num_colours; i++) {

      // i is the colour
      //
      num_lines =  bonds[i].size();
      box.symmetry_bonds_[i].num_lines = num_lines;

      if (num_lines > 0 ){

         // box.symmetry_bonds_[i].pair_list = new coot::CartesianPair[bonds[i].size()];
         box.symmetry_bonds_[i].pair_list = new graphics_line_t[bonds[i].size()];

         int bis = bonds[i].size();
         for (int j=0; j<bis; j++) {
            box.symmetry_bonds_[i].pair_list[j] = bonds[i][j];
         }
      }
   }
   return box;
}

void
Bond_lines_container::check() const {
   //
   std::cout << "Bond_lines_container::check() bonds.size() " << bonds.size() << std::endl;
   if (bonds.size() > 0) {
      std::cout <<  "Bond_lines_container::check() bonds[0].size(): "
                << bonds[0].size() << std::endl;
   }
   if (bonds.size() > 1) {
      std::cout <<  "Bond_lines_container::check() bonds[1].size(): "
                << bonds[1].size() << std::endl;
   }
}


// initialise
// graphical_bonds_container::bonds

//
void
Bond_lines_container::write(std::string filename) const {


   std::cout << "Write bonds to file: " << filename.c_str() << std::endl;

   std::ofstream bondsout(filename.c_str());
   if (! bondsout) {
      // error
      std::cout << "Could not open " << filename.c_str() << " for some reason\n";
   } else {

      for (unsigned int i = 0; i < bonds.size(); i++) {

         bondsout<< bonds[i].size() << " bonds of colour " << i << std::endl;

         int bis=bonds[i].size();
         for (int j = 0; j<bis ; j++) {

            // This gets it pass CC, strangely
            bondsout << bonds[i].GetStart(j);
            bondsout << " to ";
            bondsout << bonds[i].GetFinish(j) << std::endl;
         }
      }
   }
}

const graphics_line_t &
Bond_lines::operator[](unsigned int i) const {

   return points[i];
}


const coot::Cartesian &
Bond_lines::GetStart(unsigned int i) const {

   return points[i].positions.getStart();
}

const coot::Cartesian &
Bond_lines::GetFinish(unsigned int i) const {

   return points[i].positions.getFinish();
}

Bond_lines_container::Bond_lines_container(int col) {

   do_bonds_to_hydrogens = 1;  // added 20070629
   have_dictionary = 0;
   init();
   b_factor_scale = 1.0;
   for_GL_solid_model_rendering = 1;
   std::cout << "Strange Bond_lines_container(int col)" << std::endl;
   Bond_lines a;
   bonds.push_back(a);
}

// #ifdef USE_BACKWARD
// #include <utils/backward.hpp>
// #endif

// these are all-molecule atom indices (or should be) - not residue atom indices
void
Bond_lines_container::addBond(int colour_index,
                              const coot::Cartesian &start,
                              const coot::Cartesian &end,
                              graphics_line_t::cylinder_class_t cc,
                              int model_number,
                              int atom_index_1,
                              int atom_index_2,
                              bool add_begin_end_cap,   // default arg
                              bool add_end_end_cap      // default arg
                              ) {

// #ifdef USE_BACKWARD
//             backward::StackTrace st;
//             backward::Printer p;
//             st.load_here(32);
//             p.print(st);
// #endif

   // Needs further investigation

   // if (no_bonds_to_these_atoms.size() > 0)
   // std::cout << "debug in addBond() " << no_bonds_to_these_atoms.size() << " " << atom_index_1 << " " << atom_index_2 << std::endl;

   // if the atom selection has the same size as no_bonds_to_these_atoms then we don't want
   // to draw any bonds
   //
   if (static_cast<int>(no_bonds_to_these_atoms.size()) == n_atoms_in_atom_selection)
      if (n_atoms_in_atom_selection > 0)
         return;

#if 0 // 2026-02-22-PE I dont think this is right!
      // I think we want to duck out if *either* of the atoms are in the non-bonds-to-these-atoms set
   // duck out if both of these atoms are in the no-bonds to these atoms set
   //
   if (no_bonds_to_these_atoms.find(atom_index_1) != no_bonds_to_these_atoms.end()) {
      if (no_bonds_to_these_atoms.find(atom_index_2) != no_bonds_to_these_atoms.end()) {
              if (false)
                 std::cout << "debug::: addBond() ducking out with " << atom_index_1 << " " << atom_index_2
                      << " set size: " << no_bonds_to_these_atoms.size() << std::endl;
              return;
      }
   }
#endif

   // 2026-02-22 new logic
   if (no_bonds_to_these_atoms.find(atom_index_1) != no_bonds_to_these_atoms.end()) return;
   if (no_bonds_to_these_atoms.find(atom_index_2) != no_bonds_to_these_atoms.end()) return;

   if (false)
      if (no_bonds_to_these_atoms.size() > 0)
         std::cout << "debug in addBond() atom indices " << atom_index_1 << " " << atom_index_2
                   << " not found in no_bonds_to_these_atoms" << std::endl;

   coot::CartesianPair pair(start,end);
   int bonds_size = bonds.size();
   // std::cout << "in addBond()  colour_index " << colour_index << " bonds.size(): " << bonds.size() << std::endl;

   if (colour_index == -1) {
      std::cout << "ERROR:: colour_index is -1!" << std::endl;
   } else {
      if (colour_index >= bonds_size) {
         bonds.resize(colour_index+1);
         // std::cout << "in addBond() resizeing bonds to " << bonds.size() << std::endl;
         bonds[colour_index].add_bond(pair, cc, add_begin_end_cap, add_end_end_cap, model_number,
                                      atom_index_1, atom_index_2);
      } else {
         // normal path
         bonds[colour_index].add_bond(pair, cc, add_begin_end_cap, add_end_end_cap, model_number,
                                      atom_index_1, atom_index_2);
      }
   }
}


// if half_bond_type_flag is HALF_BOND_FIRST_ATOM this is called with atom-pos, mid-pod
// if half_bond_type_flag is HALF_BOND_SECOND_ATOM this is called with mid-pod, atom_pos
//
// When we have ca+ligands representation, this function can get
// called with col=20 when bonds.size() is 10.  Hmm..
//
void
Bond_lines_container::add_dashed_bond(int col,
                                      const coot::Cartesian &start_in,
                                      const coot::Cartesian &end_in,
                                      int half_bond_type_flag,
                                      graphics_line_t::cylinder_class_t cc,
                                      int imodel_number,
                                      int atom_index_1, int atom_index_2) {

   if (false) // debugging.
      std::cout << "   .... in add_dashed_bond() col is " << col
                << " and bonds.size() " << bonds.size()
                << " and half_bond_type_flag is " << half_bond_type_flag
                << std::endl;

   float dash_start = 0;
   float dash_end   = 19;

   // duck out if both of these atoms are in the no-bonds to these atoms set
   if (no_bonds_to_these_atoms.find(atom_index_1) != no_bonds_to_these_atoms.end())
      if (no_bonds_to_these_atoms.find(atom_index_2) != no_bonds_to_these_atoms.end())
         return;

   coot::Cartesian start = start_in;
   coot::Cartesian end   = end_in;

   if (half_bond_type_flag == HALF_BOND_FIRST_ATOM)
      dash_end = 9.5;

   if (half_bond_type_flag == HALF_BOND_SECOND_ATOM) {
      dash_start = 0;
      dash_end = 9.5;
      end = start_in;
      start = end_in;
   }

   float n_dash = dash_end - dash_start;
   coot::Cartesian delta = end - start;

   //                        1 1 1 1 1 1 1 1 1
   //   0 1 2 3 4 5 6 7 8 9  0 1 2 3 4 5 6 7 8
   //  X--__--__--__--__--_:_--__--__--__--__--X 19
   //

   // dash is, for example, in the range 0 to 9.5
   //
   if (col < int(bonds.size())) {
   } else {
      // 20230523-PE is this safe!? (Well, we do it elsewhere, so I guess so).
      // Anyway, why are we here with out of range colour? Aren't links drawn after metals?
      bonds.resize(col+1);
   }
   for (float dash=dash_start; dash<=dash_end; dash+=2) {

      float frac_1 = dash/n_dash;              // frac_1 in the range 0 to 1
      float frac_2 = (dash+1)/n_dash;
      coot::Cartesian f_s(start + delta * frac_1);
      coot::Cartesian f_e(start + delta * frac_2);
      coot::CartesianPair pair(f_s, f_e);
      graphics_line_t::cylinder_class_t cc = graphics_line_t::DOUBLE;
      bonds[col].add_bond(pair, cc, true, true, imodel_number, atom_index_1, atom_index_2);
   }
}

//
void
Bond_lines::add_bond(const coot::CartesianPair &p,
                     graphics_line_t::cylinder_class_t cc,
                     bool begin_end_cap,
                     bool end_end_cap,
                     int model_number,
                     int atom_index_1, int atom_index_2) {

   graphics_line_t gl(p, cc, begin_end_cap, end_end_cap, model_number, atom_index_1, atom_index_2);
   points.push_back(gl);
}


int
Bond_lines::size() const {
   return points.size();
}



// ------------------------------------------------------------------
//            Alpha Carbon Trace, C-alpha calpha
// ------------------------------------------------------------------


// The distances are the minimum and maximum distances to look check
// between for ca-ca (pseudo) bond.  Typically 3.6-3.8 (tight) or
// 0.0-5.0 (loose).
//
// This can get called for intermediate atoms before the restraints have been made
// (and the FixedDuringRefinement UDD is set), so that loops can flash on
// then off (when FixedDuringRefinement UDDs *are* set).
// Oh well, a brief flash is better than permanently on during refinement.
void
Bond_lines_container::do_Ca_bonds(atom_selection_container_t SelAtom,
                                  float min_dist, float max_dist, bool draw_missing_loops_flag) {

   if (udd_has_ca_handle == -1)
      udd_has_ca_handle = SelAtom.mol->RegisterUDInteger(mmdb::UDR_RESIDUE, "has CA");
   if (false)
      std::cout << "debug:: RegisterUDInteger for udd_has_ca_handle " << udd_has_ca_handle
                << std::endl;
   if (!udd_has_ca_handle) {
      std::cout << "ERROR getting udd_has_ca_handle\n";
   }
   coot::my_atom_colour_map_t atom_colour_map_in;
   coot::my_atom_colour_map_t atom_colour_map =
      do_Ca_or_P_bonds_internal(SelAtom, " CA ", atom_colour_map_in,  // PDBv3 FIXME
                                min_dist, max_dist, draw_missing_loops_flag, coot::COLOUR_BY_CHAIN);
   do_Ca_or_P_bonds_internal(SelAtom, " P  ",  atom_colour_map,
                             0.1,      7.5, draw_missing_loops_flag, coot::COLOUR_BY_CHAIN);
}

void
Bond_lines_container::do_Ca_loop(int imod, int ires, int nres,
                                 mmdb::Chain *chain_p,
                                 mmdb::Residue *residue_prev,
                                 mmdb::Residue *residue_this,
                                 int udd_atom_index_handle,
                                 int udd_fixed_during_refinement_handle) {

   std::string res_name_1 = residue_prev->GetResName();
   std::string res_name_2 = residue_this->GetResName();

   if (res_name_1 == "HOH") return;
   if (res_name_2 == "HOH") return;

   if (false)
      std::cout << "loop this? " << coot::residue_spec_t(residue_prev)
                << " " << coot::residue_spec_t(residue_this) << std::endl;

   // we want to represent the missing residues as a curved loop. To do so, we
   // need to find the positions of the CA of n-2 (for start of loop)
   // and CA of n+2 (for end of the loop) - can we do that?
   int res_idx_n_start_back = (ires - 1) -2;
   int res_idx_n_end_forwards = ires + 2;
   if (res_idx_n_start_back >= 0) {
      if (res_idx_n_end_forwards < nres) {
         mmdb::Residue *res_start_back   = chain_p->GetResidue(res_idx_n_start_back);
         mmdb::Residue *res_end_forwards = chain_p->GetResidue(res_idx_n_end_forwards);
         if (res_start_back) {
            if (res_end_forwards) {
               mmdb::PPAtom residue_atoms_pp_1 = 0;
               mmdb::PPAtom residue_atoms_pp_2 = 0;
               mmdb::PPAtom residue_atoms_pp_3 = 0;
               mmdb::PPAtom residue_atoms_pp_4 = 0;
               int n_atoms_pp_1 = 0;
               int n_atoms_pp_2 = 0;
               int n_atoms_pp_3 = 0;
               int n_atoms_pp_4 = 0;
               mmdb::Atom *at_pp_1 = 0;
               mmdb::Atom *at_pp_2 = 0;
               mmdb::Atom *at_pp_3 = 0;
               mmdb::Atom *at_pp_4 = 0;

               // we want to not draw loops if the atoms
               // are from a moving atoms molecule
               // (in which case udd_fixed_during_refinement_handle will
               // have been set) and the atoms are not fixed
               //
               // or
               //
               // we are not in a moving atoms molecule
               //
               bool loop_is_possible = true;
               bool these_are_moving_atoms = false;
               int  udd_is_fixed_during_refinement = 0;

               if (udd_fixed_during_refinement_handle > 0) {
                  these_are_moving_atoms = true; // a moving atoms molecule, I mean
                                                 // (that can contain fixed atoms)
               }

               res_start_back->GetAtomTable(residue_atoms_pp_1, n_atoms_pp_1);
               residue_prev->GetAtomTable(residue_atoms_pp_2, n_atoms_pp_2);
               residue_this->GetAtomTable(residue_atoms_pp_3, n_atoms_pp_3);
               res_end_forwards->GetAtomTable(residue_atoms_pp_4, n_atoms_pp_4);

               if (loop_is_possible) { // true
                  for (int ipp_1=0; ipp_1<n_atoms_pp_1; ipp_1++) {
                     mmdb::Atom *at = residue_atoms_pp_1[ipp_1];
                     std::string atom_name = at->GetAtomName();
                     if (! at->Het) {
                        if (atom_name == " CA " || atom_name == " P  ") {
                           at_pp_1 = at;
                           if (these_are_moving_atoms) {
                              at->GetUDData(udd_fixed_during_refinement_handle,
                                            udd_is_fixed_during_refinement);
                              if (udd_is_fixed_during_refinement == 1) {
                                 at_pp_1 = 0; // nullptr
                                 loop_is_possible = false;

                              }
                           } else {
                              // don't make a loop to atoms that are not draw (because they are in the moving atoms set)
                              int idx_in_mol = -1;
                              at->GetUDData(udd_atom_index_handle, idx_in_mol);
                              if (no_bonds_to_these_atoms.find(idx_in_mol) != no_bonds_to_these_atoms.end())
                                 loop_is_possible = false;
                           }
                           break;
                        }
                     }
                  }
               }

               if (udd_fixed_during_refinement_handle == -1)
                  loop_is_possible = true;

               if (loop_is_possible) {
                  for (int ipp_2=0; ipp_2<n_atoms_pp_2; ipp_2++) {
                     mmdb::Atom *at = residue_atoms_pp_2[ipp_2];
                     std::string atom_name = at->GetAtomName();
                     // start atom is preferably O3'
                     if (atom_name == " CA " || atom_name == " P  " || atom_name == " O3'") {
                        if (! at->Het) {
                           at_pp_2 = at;
                           if (these_are_moving_atoms) {
                              at->GetUDData(udd_fixed_during_refinement_handle,
                                            udd_is_fixed_during_refinement);
                              if (udd_is_fixed_during_refinement == 1) {
                                 at_pp_2 = 0; // nullptr
                                 loop_is_possible = false;
                              }
                           } else {
                              // don't make a loop to atoms that are not draw (because they are in the moving atoms set)
                              int idx_in_mol = -1;
                              at->GetUDData(udd_atom_index_handle, idx_in_mol);
                              if (no_bonds_to_these_atoms.find(idx_in_mol) != no_bonds_to_these_atoms.end())
                                 loop_is_possible = false;
                           }
                           // try to find the O3' if we can.
                           if (atom_name == " O3'")
                              break;
                        }
                     }
                  }
               }

               if (loop_is_possible) {
                  for (int ipp_3=0; ipp_3<n_atoms_pp_3; ipp_3++) {
                     mmdb::Atom *at = residue_atoms_pp_3[ipp_3];
                     std::string atom_name = at->GetAtomName();
                     if (! at->Het) {
                        if (atom_name == " CA " || atom_name == " P  ") {
                           at_pp_3 = at;
                           if (these_are_moving_atoms) {
                              at->GetUDData(udd_fixed_during_refinement_handle,
                                            udd_is_fixed_during_refinement);
                              if (udd_is_fixed_during_refinement == 1) {
                                 at_pp_3 = 0; // nullptr
                                 loop_is_possible = false;
                              }
                           } else {
                              // don't make a loop to atoms that are not draw (because they are in the moving atoms set)
                              int idx_in_mol = -1;
                              at->GetUDData(udd_atom_index_handle, idx_in_mol);
                              if (no_bonds_to_these_atoms.find(idx_in_mol) != no_bonds_to_these_atoms.end())
                                 loop_is_possible = false;
                           }
                           break;
                        }
                     }
                  }
               }

               if (loop_is_possible) {
                  for (int ipp_4=0; ipp_4<n_atoms_pp_4; ipp_4++) {
                     mmdb::Atom *at = residue_atoms_pp_4[ipp_4];
                     std::string atom_name = at->GetAtomName();
                     if (! at->Het) {
                        if (atom_name == " CA " || atom_name == " P  ") {
                           at_pp_4 = at;
                           if (these_are_moving_atoms) {
                              at->GetUDData(udd_fixed_during_refinement_handle,
                                            udd_is_fixed_during_refinement);
                              if (udd_is_fixed_during_refinement == 1) {
                                 at_pp_4 = 0; // nullptr
                                 loop_is_possible = false;
                              }
                           } else {
                              // don't make a loop to atoms that are not draw (because they are in the moving atoms set)
                              int idx_in_mol = -1;
                              at->GetUDData(udd_atom_index_handle, idx_in_mol);
                              if (no_bonds_to_these_atoms.find(idx_in_mol) != no_bonds_to_these_atoms.end())
                                 loop_is_possible = false;
                           }
                           break;
                        }
                     }
                  }
               }

               // Are C of previous and N of next close as in a peptide bond?
               // Then we don't want to draw a loop
               if (loop_is_possible) {
                  bool C_and_N_are_close = false;
                  mmdb::Atom *C_prev = residue_prev->GetAtom(" C  ");
                  if (C_prev) {
                     mmdb::Atom *N_this = residue_this->GetAtom(" N  ");
                     if (N_this) {
                        float dist_sqrd =
                           (C_prev->x - N_this->x) * (C_prev->x - N_this->x) +
                           (C_prev->y - N_this->y) * (C_prev->y - N_this->y) +
                           (C_prev->z - N_this->z) * (C_prev->z - N_this->z);
                        if (dist_sqrd < 2.5 * 2.5)
                           C_and_N_are_close = true;
                     }
                  }
                  if (C_and_N_are_close)
                     loop_is_possible = false;
               }

               // Are P of previous and O3' of prev close as in a phosphodiester?
               // Then we don't want to draw a loop
               if (loop_is_possible) {
                  bool P_and_O3prime_are_close = false;
                  mmdb::Atom *O3prime_prev = residue_prev->GetAtom(" O3'");
                  if (O3prime_prev) {
                     mmdb::Atom *P_this = residue_this->GetAtom(" P  ");
                     if (P_this) {
                        float dist_sqrd =
                           (O3prime_prev->x - P_this->x) * (O3prime_prev->x - P_this->x) +
                           (O3prime_prev->y - P_this->y) * (O3prime_prev->y - P_this->y) +
                           (O3prime_prev->z - P_this->z) * (O3prime_prev->z - P_this->z);
                        if (dist_sqrd < 2.5 * 2.5) // c.f 1.6 * 1.6
                           P_and_O3prime_are_close = true;
                     }
                  }
                  if (P_and_O3prime_are_close)
                     loop_is_possible = false;
               }

               if (loop_is_possible) {

                  if (at_pp_1 && at_pp_2 && at_pp_3 && at_pp_4) {
                     coot::Cartesian pp_2(at_pp_2->x, at_pp_2->y, at_pp_2->z);
                     coot::Cartesian pp_3(at_pp_3->x, at_pp_3->y, at_pp_3->z);
                     float a = (pp_3-pp_2).amplitude();
                     int n_line_segments = static_cast<int>(a*1.2);
                     std::pair<bool, std::vector<coot::CartesianPair> > lp =
                        coot::loop_path(at_pp_1, at_pp_2, at_pp_3, at_pp_4, n_line_segments);
                     for (unsigned int jj=0; jj<lp.second.size(); jj++) {
                        const coot::CartesianPair &cp = lp.second[jj];
                        int col = GREY_BOND; // just grey, really.
                        graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                        int iat_1 = -1111; // signifying a CA-loop
                        int iat_2 = -1111;
                        addBond(col, cp.getStart(), cp.getFinish(), cc, imod, iat_1, iat_2, true, true);
                     }

                     if (lp.first) {
                        // Add the balls of CA-CA badness
                        for (unsigned int jj=0; jj<lp.second.size(); jj++) {
                           const coot::CartesianPair &cp = lp.second[jj];
                           bad_CA_CA_dist_spots.push_back(cp.getStart());
                           bad_CA_CA_dist_spots.push_back(cp.getFinish());
                        }
                     }

                  } else {
                     if (false)
                        std::cout << "DEBUG:: do_CA_loops(): oops! "
                                  << at_pp_1 << " " << at_pp_2 << " " << at_pp_3 << " " << at_pp_4
                                  << std::endl;
                  }
               }
            }
         }
      }
   }
}


// where bond_colour_type is one of
// COLOUR_BY_CHAIN=0, COLOUR_BY_SEC_STRUCT=1
//
// Don't show HETATMs (but allow MSE in CA mode).
//
coot::my_atom_colour_map_t
Bond_lines_container::do_Ca_or_P_bonds_internal(atom_selection_container_t SelAtom,
                                                const char *backbone_atom_id,
                                                coot::my_atom_colour_map_t atom_colour_map,
                                                float min_dist, float max_dist,
                                                bool draw_missing_loops_flag,
                                                int bond_colour_type,
                                                int model_number) {


   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

   if (false)
      std::cout << "---- in do_Ca_or_P_bonds_internal() with bond_colour_type "
                << bond_colour_type << std::endl;

   int atom_colours_udd = -1; // unset/bad
   int udd_handle_for_user_defined_colours = -1;
   int udd_atom_index_handle = SelAtom.UDDAtomIndexHandle;

   int udd_fixed_during_refinement_handle =
      SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM, "FixedDuringRefinement");
   // that might fail, a return value of 0 represents failure to find that handle name
   // positive means that it was OK!

   // heuristic cut off for when user has omitted GAP cards.
   //
   float dist_max_CA_CA = 5.0;
   float dist_max_P_P   = 8.5; // seems reasonable.

   // we want to add a dotted loop where there is a break in the sequence of more than
   // 1 residue and the distance between the 2 ends is more than dist_max_CA_CA.

   if (false)
      std::cout << "------------------- debug:: do_Ca_or_P_bonds_internal() bond_colour_type "
                << bond_colour_type << " c.f. "
                << coot::COLOUR_BY_RAINBOW << " "
                << coot::COLOUR_BY_B_FACTOR << " "
                << Bond_lines_container::COLOUR_BY_B_FACTOR << std::endl;

   if (bond_colour_type == coot::COLOUR_BY_RAINBOW)
      atom_colours_udd = set_rainbow_colours(SelAtom.mol);

   // the name confusion/colision here is ridiculous
   if (bond_colour_type == Bond_lines_container::COLOUR_BY_B_FACTOR) {
      atom_colours_udd = set_b_factor_colours(SelAtom.mol);
   }

   if (bond_colour_type == coot::COLOUR_BY_USER_DEFINED_COLOURS)
      udd_handle_for_user_defined_colours =
         SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   int udd_has_bond_handle = SelAtom.mol->RegisterUDInteger(mmdb::UDR_ATOM, "found-backbone-bond");
   for (int i=0; i<SelAtom.n_selected_atoms; i++)
      SelAtom.atom_selection[i]->PutUDData(udd_has_bond_handle, 0);

   int udd_user_defined_atom_colour_index_handle = SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM,
                                                                             "user-defined-atom-colour-index");

   auto get_atom_colour_index = [udd_user_defined_atom_colour_index_handle] (mmdb:: Atom *at,
                                                                             const std::string &chain_id,
                                                                             coot::my_atom_colour_map_t &atom_colour_map) {
      int idx_col = atom_colour_map.index_for_chain(chain_id);
      int idx_col_udd;
      if (at->GetUDData(udd_user_defined_atom_colour_index_handle, idx_col_udd) == mmdb::UDDATA_Ok) {
         idx_col = idx_col_udd;
      }
      return idx_col;
   };

   auto new_style_nucleotide_backbone_chain_rep = [this, get_atom_colour_index] (int imod,  mmdb::Chain *chain_p,
                                                                                 int n_atoms_prev, mmdb::Residue *residue_prev,
                                                                                 int n_atoms_this, mmdb::Residue *residue_this,
                                                                                 const std::string &res_name_1,
                                                                                 const std::string &res_name_2,
                                                          int udd_atom_index_handle, int bond_colour_type,
                                                          int atom_colours_udd,
                                                          coot::my_atom_colour_map_t &atom_colour_map, // reference
                                                          int udd_handle_for_user_defined_colours,
                                                          int udd_has_bond_handle) {

      graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
      graphics_line_t::cylinder_class_t base_spike_cc = graphics_line_t::SINGLE;

      std::string chain_id = chain_p->GetChainID();

      // first do the P_C4_prime bonds internal to the residue_this
      //
      for (int iat=0; iat<n_atoms_this; iat++) {
         mmdb::Atom *at_1 = residue_this->GetAtom(iat);
         std::string atom_name_1(at_1->GetAtomName());
         if (atom_name_1 == " C4'") {
            for (int jat=0; jat<n_atoms_this; jat++) {
               mmdb::Atom *at_2 = residue_this->GetAtom(jat);
               std::string atom_name_2(at_2->GetAtomName());
               if (atom_name_2 == " P  ") {
                  std::string alt_conf_1 = at_1->altLoc;
                  std::string alt_conf_2 = at_2->altLoc;
                  if (alt_conf_1.empty() || alt_conf_2.empty() || alt_conf_2 == alt_conf_1) {
                     coot::Cartesian pt_1(at_1->x, at_1->y, at_1->z);
                     coot::Cartesian pt_2(at_2->x, at_2->y, at_2->z);
                     coot::Cartesian bond_mid_point = pt_1.mid_point(pt_2);
                     int iat_1 = -1;
                     int iat_2 = -1;
                     at_1->GetUDData(udd_atom_index_handle, iat_1);
                     at_2->GetUDData(udd_atom_index_handle, iat_2);
                     int col_1 = get_atom_colour_index(at_1, chain_id, atom_colour_map);
                     addBond(col_1, pt_1, bond_mid_point, cc, imod, iat_1, iat_2, false, false);
                     addBond(col_1, bond_mid_point, pt_2, cc, imod, iat_1, iat_2, false, false);
                     at_1->PutUDData(udd_has_bond_handle, 1);
                     at_2->PutUDData(udd_has_bond_handle, 1);
                     // for use with Ca+ligand mode
                     residue_this->PutUDData(udd_has_ca_handle, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                  }
               }

               bool do_base_stick_bond = false;
               if (atom_name_2 == " N1 ") {
                  if (res_name_2 == "DG" || res_name_2 == "G" || res_name_2 == "A" || res_name_2 == "DA") {
                     do_base_stick_bond = true;
                  }
               }
               if (atom_name_2 == " O4 ") {
                  if (res_name_2 == "DT" || res_name_2 == "U") {
                     do_base_stick_bond = true;
                  }
               }
               if (atom_name_2 == " N4 ") {
                  if (res_name_2 == "DC" || res_name_2 == "C") {
                     do_base_stick_bond = true;
                  }
               }
               if (do_base_stick_bond) { // at least potentially...
                  std::string alt_conf_1 = at_1->altLoc;
                  std::string alt_conf_2 = at_2->altLoc;
                  if (alt_conf_1.empty() || alt_conf_2.empty() || alt_conf_2 == alt_conf_1) {
                     coot::Cartesian pt_1(at_1->x, at_1->y, at_1->z);
                     coot::Cartesian pt_2(at_2->x, at_2->y, at_2->z);
                     int iat_1 = -1;
                     int iat_2 = -1;
                     at_1->GetUDData(udd_atom_index_handle, iat_1);
                     at_2->GetUDData(udd_atom_index_handle, iat_2);
                     int col_1 = atom_colour_map.index_for_chain(chain_id);
                     // the bools for the end caps seem to be the wrong way around
                     // but this works and the reverse does not.
                     addBond(col_1, pt_1, pt_2, cc, imod, iat_1, iat_2, false, true);
                  }
               }
            }
         }
      }

      // now inter-residue bonds

      for (int iat=0; iat<n_atoms_prev; iat++) {
         mmdb::Atom *at_1 = residue_prev->GetAtom(iat);
         std::string atom_name_1(at_1->GetAtomName());
         if (atom_name_1 == " C4'") { // PDBv3 FIXME
            for (int jat=0; jat<n_atoms_this; jat++) {
               mmdb::Atom *at_2 = residue_this->GetAtom(jat);
               std::string atom_name_2(at_2->GetAtomName());
               if (atom_name_2 == " P  ") { // PDBv3 FIXME
                  std::string alt_conf_prev = at_1->altLoc;
                  std::string alt_conf_this = at_2->altLoc;
                  if (alt_conf_prev.empty() || alt_conf_this.empty() || alt_conf_prev == alt_conf_this) {
                     coot::Cartesian pt_1(at_1->x, at_1->y, at_1->z);
                     coot::Cartesian pt_2(at_2->x, at_2->y, at_2->z);
                     int iat_1 = -1;
                     int iat_2 = -1;
                     at_1->GetUDData(udd_atom_index_handle, iat_1);
                     at_2->GetUDData(udd_atom_index_handle, iat_2);
                     int col_1 = get_atom_colour_index(at_1, chain_id, atom_colour_map);
                     addBond(col_1, pt_1, pt_2, cc, imod, iat_1, iat_2);
                     at_1->PutUDData(udd_has_bond_handle, 1);
                     at_2->PutUDData(udd_has_bond_handle, 1);
                     // for use with Ca+ligand mode
                     residue_this->PutUDData(udd_has_ca_handle, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                     residue_prev->PutUDData(udd_has_ca_handle, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                  }
               }
            }
         }
      }
   };

   auto CA_CA_or_P_P = [this, get_atom_colour_index] (int imod, mmdb::Chain *chain_p,
                               int n_atoms_prev, mmdb::Residue *residue_prev, int n_atoms_this, mmdb::Residue *residue_this,
                               const std::string &res_name_1, const std::string &res_name_2,
                               float dist_max_CA_CA, float dist_max_P_P,
                               int udd_atom_index_handle, int bond_colour_type,
                               int atom_colours_udd,
                               int udd_user_defined_atom_colour_index_handle,
                               coot::my_atom_colour_map_t &atom_colour_map, // reference
                               int udd_handle_for_user_defined_colours,
                               int udd_has_bond_handle,
                               int udd_has_ca_handle) {

      graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
      std::string chain_id = chain_p->GetChainID();

      // normal case (not a gap in the sequence)
      for (int iat=0; iat<n_atoms_prev; iat++) {
         mmdb::Atom *at_1 = residue_prev->GetAtom(iat);
         std::string atom_name_1(at_1->GetAtomName());
         if (atom_name_1 == " CA " || atom_name_1 == " P  ") { // PDBv3 FIXME
            for (int jat=0; jat<n_atoms_this; jat++) {
               mmdb::Atom *at_2 = residue_this->GetAtom(jat);
               std::string atom_name_2(at_2->GetAtomName());
               std::string alt_conf_prev = at_1->altLoc;
               std::string alt_conf_this = at_2->altLoc;
               // Allow MSE hetgroups in CA mode
               if ((!at_1->Het || res_name_1 == "MSE") && (!at_2->Het || res_name_2 == "MSE")) {
                  if (!at_1->isTer() && !at_2->isTer()) {
                     float dist_max_sqrd = dist_max_CA_CA * dist_max_CA_CA;
                     bool phosphate_pair = false;
                     bool Calpha_pair    = false;
                     if ((atom_name_1 == " P  ") && (atom_name_2 == " P  ")) {
                        phosphate_pair = 1;
                        dist_max_sqrd = dist_max_P_P * dist_max_P_P;
                     }
                     if ((atom_name_1 == " CA ") && (atom_name_2 == " CA ")) {
                        Calpha_pair = 1;
                     }

                     if (Calpha_pair || phosphate_pair) {
                        if (alt_conf_prev == alt_conf_this || alt_conf_this == "" || alt_conf_prev == "") {
                           int col = 0; // overridden.
                           coot::Cartesian ca_1(at_1->x, at_1->y, at_1->z);
                           coot::Cartesian ca_2(at_2->x, at_2->y, at_2->z);
                           int iat_1 = -1;
                           int iat_2 = -1;
                           at_1->GetUDData(udd_atom_index_handle, iat_1);
                           at_2->GetUDData(udd_atom_index_handle, iat_2);
                           if ((ca_1-ca_2).amplitude_squared() < dist_max_sqrd) {
                              if (bond_colour_type == Bond_lines_container::COLOUR_BY_B_FACTOR) {

                                 int col_1 = 0;
                                 int col_2 = 0;
                                 if (atom_colours_udd > 0) {
                                    mmdb::realtype f;
                                    if (at_1->GetUDData(atom_colours_udd, f) == mmdb::UDDATA_Ok)
                                       col_1 = atom_colour_map.index_for_b_factor(f);
                                    if (at_2->GetUDData(atom_colours_udd, f) == mmdb::UDDATA_Ok)
                                       col_2 = atom_colour_map.index_for_b_factor(f);
                                 }
                                 coot::Cartesian bond_mid_point = ca_1.mid_point(ca_2);
                                 addBond(col_1, ca_1, bond_mid_point, cc, imod, iat_1, iat_2);
                                 addBond(col_2, bond_mid_point, ca_2, cc, imod, iat_1, iat_2);
                              } else {
                                 if (bond_colour_type == coot::COLOUR_BY_SEC_STRUCT) {
                                    coot::Cartesian bond_mid_point = ca_1.mid_point(ca_2);
                                    col = atom_colour(at_1, coot::COLOUR_BY_SEC_STRUCT, udd_user_defined_atom_colour_index_handle);
                                    addBond(col, ca_1, bond_mid_point, cc, imod, iat_1, iat_2);
                                    col = atom_colour(at_2, coot::COLOUR_BY_SEC_STRUCT, udd_user_defined_atom_colour_index_handle);
                                    addBond(col, bond_mid_point, ca_2, cc, imod, iat_1, iat_2);
                                 } else {

                                    int col_1 = 0;
                                    int col_2 = 0;
                                    if (bond_colour_type == coot::COLOUR_BY_RAINBOW) {
                                       if (atom_colours_udd > 0) {
                                          mmdb::realtype f;
                                          if (at_1->GetUDData(atom_colours_udd, f) == mmdb::UDDATA_Ok) {
                                             col_1 = atom_colour_map.index_for_rainbow(f);
                                             if (at_2->GetUDData(atom_colours_udd, f) == mmdb::UDDATA_Ok) {
                                                col_2 = atom_colour_map.index_for_rainbow(f);
                                             } else {
                                                col_2 = 0;
                                             }
                                          } else {
                                             col_1 = 0;
                                          }
                                       } else {
                                          col_1 = 0;
                                       }
                                    } else {

                                       if (bond_colour_type == coot::COLOUR_BY_USER_DEFINED_COLOURS) {
                                          col_1 = this->get_user_defined_col_index(at_1, udd_handle_for_user_defined_colours);
                                          col_2 = this->get_user_defined_col_index(at_2, udd_handle_for_user_defined_colours);
                                          if (col_1 < 0) // problem
                                             col_1 = 0;
                                          if (col_2 < 0) // ditto
                                             col_2 = 0;
                                       } else {

                                          // col_1 = atom_colour_map.index_for_chain(chain_p->GetChainID());
                                          col_1 = get_atom_colour_index(at_1, chain_id, atom_colour_map);
                                          col_2 = col_1;
                                       }
                                    }
                                    this->bonds_size_colour_check(col_1);
                                    this->bonds_size_colour_check(col_2);
                                    coot::Cartesian bond_mid_point = ca_1.mid_point(ca_2);
                                    addBond(col_1, ca_1, bond_mid_point, cc, imod, iat_1, iat_2);
                                    addBond(col_2, bond_mid_point, ca_2, cc, imod, iat_1, iat_2);
                                 }
                              }
                              at_1->PutUDData(udd_has_bond_handle, 1);
                              at_2->PutUDData(udd_has_bond_handle, 1);
                              // for use with Ca+ligand mode
                              residue_this->PutUDData(udd_has_ca_handle, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                              residue_prev->PutUDData(udd_has_ca_handle, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   };

   // default to all models
   int n_models = SelAtom.mol->GetNumberOfModels();
   int imod_start = 1;
   int imod_end = n_models;
   if (model_number != 0) {
      if (model_number <= n_models) {
         imod_start = model_number;
         imod_end   = model_number;
      }
   }

   for(int imod = imod_start; imod<=imod_end; imod++) {
      mmdb::Model *model_p = SelAtom.mol->GetModel(imod);
      if (model_p) {
         int nchains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<nchains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            bool is_nucleotide_chain_flag = chain_p->isNucleotideChain();
            int nres = chain_p->GetNumberOfResidues();
            if (nres < 2) continue;
            for (int ires=1; ires<nres; ires++) {
               mmdb::Residue *residue_this = chain_p->GetResidue(ires);
               mmdb::Residue *residue_prev = chain_p->GetResidue(ires-1);
               if (residue_this && residue_prev) {

                  int n_atoms_prev = residue_prev->GetNumberOfAtoms();
                  int n_atoms_this = residue_this->GetNumberOfAtoms();
                  std::string res_name_1 = residue_prev->GetResName();
                  std::string res_name_2 = residue_this->GetResName();

                  int res_no_1 = residue_prev->GetSeqNum();
                  int res_no_2 = residue_this->GetSeqNum();
                  int res_no_delta = res_no_2 - res_no_1;
                  if (res_no_delta > 1) {

                     if (draw_missing_loops_flag)
                        do_Ca_loop(imod, ires, nres, chain_p, residue_prev, residue_this,
                                   udd_atom_index_handle, udd_fixed_during_refinement_handle);

                  } else {

                     bool old_style = false;
                     if (old_style) {
                        CA_CA_or_P_P(imod, chain_p, n_atoms_prev, residue_prev, n_atoms_this, residue_this,
                                     res_name_1, res_name_2, dist_max_CA_CA, dist_max_P_P,
                                     udd_atom_index_handle, bond_colour_type,
                                     atom_colours_udd, udd_user_defined_atom_colour_index_handle,
                                     atom_colour_map,
                                     udd_handle_for_user_defined_colours,
                                     udd_has_bond_handle,
                                     udd_has_ca_handle);
                     } else {
                        if (is_nucleotide_chain_flag) {
                           new_style_nucleotide_backbone_chain_rep(imod, chain_p, n_atoms_prev, residue_prev, n_atoms_this, residue_this,
                                                                   res_name_1, res_name_2,
                                                                   udd_atom_index_handle, bond_colour_type,
                                                                   atom_colours_udd, atom_colour_map,
                                                                   udd_handle_for_user_defined_colours,
                                                                   udd_has_bond_handle);
                        } else {
                           // std::cout << "was not nucleic acid chain .............................. " << chain_p->GetChainID() << std::endl;
                           CA_CA_or_P_P(imod, chain_p, n_atoms_prev, residue_prev, n_atoms_this, residue_this,
                                        res_name_1, res_name_2, dist_max_CA_CA, dist_max_P_P,
                                        udd_atom_index_handle, bond_colour_type,
                                        atom_colours_udd, udd_user_defined_atom_colour_index_handle,
                                        atom_colour_map,
                                        udd_handle_for_user_defined_colours,
                                        udd_has_bond_handle,
                                        udd_has_ca_handle);
                        }
                     }
                  }
               }
            }
         }
      }
   }


   // stars if needed:
   //
   float star_size = 0.2;
   // for atoms with no neighbour (contacts):
   coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
   coot::Cartesian small_vec_y(0.0, star_size, 0.0);
   coot::Cartesian small_vec_z(0.0, 0.0, star_size);

   int udd_status;
   int udd_value;
   for(int imod = imod_start; imod<=imod_end; imod++) {
      mmdb::Model *model_p = SelAtom.mol->GetModel(imod);
      if (model_p) {
         mmdb::Chain *chain_p;
         int nchains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<nchains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     std::string atom_name(at->GetAtomName());
                     if (atom_name == " CA " || atom_name == " P  ") {
                        udd_status = at->GetUDData(udd_has_bond_handle, udd_value);
                        // All atoms have OK status, some are marked as not having a bond
                        //
                        // then it has a bond, if not we need to mark it as a star.

                        if (udd_status != mmdb::UDDATA_Ok || udd_value == 0) {
                           int col = 0;
                           if (atom_colours_udd > 0) {
                              mmdb::realtype f;
                              if (at->GetUDData(atom_colours_udd, f) == mmdb::UDDATA_Ok)
                                 col = atom_colour_map.index_for_chain(chain_p->GetChainID());
                           }
                           if (bond_colour_type == Bond_lines_container::COLOUR_BY_B_FACTOR)
                              col = atom_colour(at, coot::COLOUR_BY_B_FACTOR, udd_user_defined_atom_colour_index_handle);
                           if (bond_colour_type == coot::COLOUR_BY_SEC_STRUCT)
                              col = atom_colour(at, bond_colour_type, udd_user_defined_atom_colour_index_handle);

                           int iat_1 = -1; // 20171224-PE FIXME
                           int udd_status_1 = at->GetUDData(udd_atom_index_handle, iat_1);
                           coot::Cartesian pos(at->x, at->y, at->z);
                           addBond(col, pos+small_vec_x, pos-small_vec_x, cc, imod, iat_1, iat_1, true, true);
                           addBond(col, pos+small_vec_y, pos-small_vec_y, cc, imod, iat_1, iat_1, true, true);
                           addBond(col, pos+small_vec_z, pos-small_vec_z, cc, imod, iat_1, iat_1, true, true);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return atom_colour_map;
}

// return the udd_handle of the UDReal values for "rainbow circle point"
int
Bond_lines_container::set_rainbow_colours(mmdb::Manager *mol) {

   int udd_handle = mol->RegisterUDReal(mmdb::UDR_ATOM, "rainbow circle point");
   if (udd_handle > 0) {

      int n_models = mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int nchains = model_p->GetNumberOfChains();
            for (int ich=0; ich<nchains; ich++) {
               mmdb::Chain *chain_p = model_p->GetChain(ich);
               int nres = chain_p->GetNumberOfResidues();
               int seq_no_max = mmdb::MinInt4;
               int seq_no_min = mmdb::MaxInt4;

               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  std::string res_name(residue_p->GetResName());
                  if (res_name != "HOH") {
                     if (coot::util::is_standard_residue_name(res_name)) {
                        int seq_no_this = ires;
                        if (seq_no_this < seq_no_min) {
                           seq_no_min = seq_no_this;
                        }
                        if (seq_no_this > seq_no_max) {
                           seq_no_max = seq_no_this;
                        }
                     }
                  }
               }
               if ((seq_no_max != mmdb::MinInt4) && (seq_no_min != mmdb::MaxInt4)) {

                  if (seq_no_min < seq_no_max) {
                     for (int ires=0; ires<nres; ires++) {
                        mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                        float range = seq_no_max - seq_no_min;
                        float chain_pos = float(ires)/range;
                        if (chain_pos < 0)
                           chain_pos = 0;
                        if (chain_pos > 1)
                           chain_pos = 1;
                        int n_atoms = residue_p->GetNumberOfAtoms();
                        for (int iat=0; iat<n_atoms; iat++) {
                           mmdb::Atom *atom_p = residue_p->GetAtom(iat);
                           if (! atom_p->Het) {
                              atom_p->PutUDData(udd_handle, chain_pos);
                           } else {
                              atom_p->PutUDData(udd_handle, 0.88);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return udd_handle;
}


// return the udd_handle of the UDReal values for "rainbow circle point"
int
Bond_lines_container::set_b_factor_colours(mmdb::Manager *mol) {

   const float max_b_factor = 70.0; // after scaling
   int udd_handle = mol->RegisterUDReal(mmdb::UDR_ATOM, "B-factor fraction point");
   if (udd_handle > 0) {

      int n_models = mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int nchains = model_p->GetNumberOfChains();
            for (int ich=0; ich<nchains; ich++) {
               mmdb::Chain *chain_p = model_p->GetChain(ich);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  std::string res_name(residue_p->GetResName());
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *atom_p = residue_p->GetAtom(iat);
                     if (! atom_p->Het) {
                        float b_factor = atom_p->tempFactor;
                        float bs = b_factor * b_factor_scale;
                        float f = bs/max_b_factor;
                        if (f < 0.0) f = 0.0;
                        if (f > 1.0) f = 1.0;
                        atom_p->PutUDData(udd_handle, f);
                     }
                  }
               }
            }
         }
      }
   }
   return udd_handle;
}

#include "geometry/main-chain.hh"
#include "geometry/hydrophobic.hh"


// atom_colour_map is an optional arg.  It is passed in the case of
// long_bonded atoms or MET/MSE residues. Default value 0 (NULL).
//
int
Bond_lines_container::atom_colour(mmdb::Atom *at, int bond_colour_type,
                                  int udd_user_defined_atom_colour_index_handle,
                                  coot::my_atom_colour_map_t *atom_colour_map_p) { // atom_colour_map_in is an optional arg

   if (false)
      std::cout << "in atom_colour() with at " << at
                << " bond_colour_type " << bond_colour_type << " vs (user-defined) " << coot::COLOUR_BY_USER_DEFINED_COLOURS
                << " vs (colour-by-atom-type) " << coot::COLOUR_BY_ATOM_TYPE
                << " vs (colour-by-b-factor) " << coot::COLOUR_BY_B_FACTOR
                << std::endl;

   int col = 0;

   // Does this atom have an over-riding/user-defined colour?
   // User-defined colours trump everything.
   int idx_col_udd;
   if (at->GetUDData(udd_user_defined_atom_colour_index_handle, idx_col_udd) == mmdb::UDDATA_Ok) {
      // std::cout << "in atom_colour(): for atom " << at << " using user defined colour " << idx_col_udd << std::endl;

      if (idx_col_udd == -1) { // -1 is a disaster, because bonds[col] is used downstream
         // let's ignore udd_user_defined_atom_colour indices if they are -1.
      } else {
         return idx_col_udd;
      }
   }

   if (bond_colour_type == coot::COLOUR_BY_MOLECULE) return col; // one colour fits all

   if (bond_colour_type == coot::COLOUR_BY_CHAIN) {

      if (atom_colour_map_p) {
         col = atom_colour_map_p->index_for_chain(std::string(at->GetChainID()));
         if (false)
            std::cout << " atom_colour_map->index_for_chain(\"" << at->GetChainID()
                      << "\") returns " << col << std::endl;
      } else {
         // std::cout << "no atom colour map" << std::endl;
      }
   } else {

      if (bond_colour_type == coot::COLOUR_BY_CHAIN_GOODSELL) {

         if (atom_colour_map_p) {
            std::string ch_id(std::string(at->GetChainID()));
            int col_idx = atom_colour_map_p->index_for_chain(ch_id);
            col = 2 * col_idx;
            std::string ele = at->element;
            if (ele != " C")
               col += 1;
            // std::cout << "here in atom_colour(): with goodsell colours with chain-id " << ch_id
            //           << " coL_idx " << col_idx << " col " << col << " ele " << ele << std::endl;
         }

      } else {

         if (bond_colour_type == coot::COLOUR_BY_HYDROPHOBIC_SIDE_CHAIN) {
            mmdb::Residue *r = at->residue;
            if (r) {
               std::string res_name(r->GetResName());
               if (coot::util::is_standard_amino_acid_name(res_name)) {
                  std::string atom_name(at->GetAtomName());
                  if (coot::is_main_chain_p(at)) {
                     col = 50; // or the chain indexed colour in future
                  } else {
                     if (coot::is_hydrophobic_atom(res_name, atom_name))
                        col = 1;
                     else
                        col = 2;
                  }
               }
            }
         }

         if (bond_colour_type == coot::COLOUR_BY_SEC_STRUCT) {
            int sse = at->residue->SSE;
            switch (sse)  {
            case mmdb::SSE_None:
               col = 0;
               break;
            case mmdb::SSE_Strand:
               col = 1;
               break;
            case mmdb::SSE_Bulge:
               col = 1;
               break;
            case mmdb::SSE_3Turn:
               col = 2;
               break;
            case mmdb::SSE_4Turn:
               col = 2;
               break;
            case mmdb::SSE_5Turn:
               col = 2;
               break;
            case mmdb::SSE_Helix:
               col = 2;
               break;
            default:
               col = 3;
            }
         } else {
            if (bond_colour_type == coot::COLOUR_BY_ATOM_TYPE) {
               std::string element(at->element);

               if (element == " C") {
                  return CARBON_BOND;
               } else {
                  if (element == " N") {
                     return BLUE_BOND;
                  } else {
                     if (element == " O") {
                        return RED_BOND;
                     } else {
                        if (element == " S") {
                           return YELLOW_BOND;
                        } else {
                           if (is_hydrogen(element)) {
                              if (is_deuterium(element))
                                 return DEUTERIUM_PINK;
                              else
                                 return HYDROGEN_GREY_BOND;
                           } else {
                              if (element == " P") {
                                 return ORANGE_BOND;
                              } else {
                                 if (element == " F") {
                                    return GREEN_BOND;
                                 } else {
                                    if (element == "CL" || element == "Cl") {
                                       return GREEN_BOND;
                                    } else {
                                       if (element == "BR") {
                                          return DARK_BROWN_BOND;
                                       } else {
                                          if (element == " I") {
                                             return DARK_VIOLET;
                                          } else {
                                             if (element == " B") {
                                                return BORON_PINK;
                                             } else {
                                                if (element == "MG" || element == "BE" || element == "CA" || element == "SR" || element == "BA") {
                                                   return DARK_GREEN_BOND;
                                                } else {
                                                   if (element == "FE") {
                                                      return DARK_ORANGE_BOND;
                                                   } else {
                                                      if (element == "LI" || element == "NA" || element == " K" || element == "RB" || element == "CS" || element == "FR") {
                                                         return VIOLET;
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
                     }
                  }
               }
               return GREY_BOND;
            } else {

               if (bond_colour_type == coot::COLOUR_BY_CHAIN_C_ONLY) {
                  std::string element(at->element);

                  if (element == " C") {   // PDBv3 FIXME (and below)
                     if (atom_colour_map_p) {
                        int l_col = atom_colour_map_p->index_for_chain(std::string(at->GetChainID()));
                        return l_col;
                     } else {
                        if (false)
                           std::cout << "ERROR:: Null atom_colour_map_p with COLOUR_BY_CHAIN_C_ONLY mode"
                                     << std::endl;
                        return col;
                     }
                  } else {
                     if (element == " N") {
                        return BLUE_BOND;
                     } else {
                        if (element == " O") {
                           return RED_BOND;
                        } else {
                           if (element == " S") {
                              return YELLOW_BOND;
                           } else {
                              if (element == " P") {
                                 return ORANGE_BOND;
                              } else {
                                 if (is_hydrogen(element)) {
                                    if (is_deuterium(element))
                                       return DEUTERIUM_PINK;
                                    else
                                       return HYDROGEN_GREY_BOND;
                                 } else {
                                    if (element == " F") {
                                       return GREEN_BOND;
                                    } else {
                                       if (element == "CL" || element == "Cl") {
                                          return GREEN_BOND;
                                       } else {
                                          if (element == "BR") {
                                             return DARK_BROWN_BOND;
                                          } else {
                                             if (element == " I") {
                                                return DARK_VIOLET;
                                             } else {
                                                if (element == " B") {
                                                   return BORON_PINK;
                                                } else {
                                                   if (element == "MG" || element == "BE" || element == "CA" || element == "SR" || element == "BA") {
                                                      return DARK_GREEN_BOND;
                                                   } else {
                                                      if (element == "FE") {
                                                         return DARK_ORANGE_BOND;
                                                      } else {
                                                         if (element == "LI" || element == "NA" || element == " K" || element == "RB" || element == "CS" || element == "FR") {
                                                            return VIOLET;
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
                        }
                     }
                  }
                  return GREY_BOND;

               } else {

                  if (bond_colour_type == coot::DISULFIDE_COLOUR) {
                     return YELLOW_BOND;
                  } else {
                     if (bond_colour_type == coot::COLOUR_BY_OCCUPANCY) {
                        if (at->occupancy > 0.95) {
                           return BLUE_BOND;
                        } else {
                           if (at->occupancy < 0.05) {
                              return RED_BOND;
                           } else {
                              if (at->occupancy > 0.7) {
                                 return CYAN_BOND;
                              } else {
                                 if (at->occupancy > 0.45) {
                                    return GREEN_BOND;
                                 } else {
                                    if (at->occupancy > 0.25) {
                                       return YELLOW_BOND;
                                    } else {
                                       return ORANGE_BOND;
                                    }
                                 }
                              }
                           }
                        }
                     } else {
                        if (bond_colour_type == coot::COLOUR_BY_B_FACTOR) {
                           // B-factors by atom are done this way.
                           float scaled_b = at->tempFactor * b_factor_scale;
                           float max_b = 100.0;
                           // std::cout << "here we go! scaled_b " << scaled_b << std::endl;
                           float f = scaled_b/max_b;
                           if (f > 0.999) f = 0.999;
                           if (f < 0.0)   f = 0.0;
                           col = static_cast<int>(f * 45); // 50 colours in the table.
                        } else {
                           if (bond_colour_type == coot::COLOUR_BY_RAINBOW) {
                              col = 20; // elsewhere
                           } else {
                              if (bond_colour_type == coot::COLOUR_BY_USER_DEFINED_COLOURS) {
                                 // up and down again...
                                 mmdb::Model *model_p = at->GetModel();
                                 if (model_p) {
                                    mmdb::Manager *mol = model_p->GetCoordHierarchy();
                                    if (mol) {
                                       int udd_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");
                                       int ic;
                                       if (at->GetUDData(udd_handle, ic) == mmdb::UDDATA_Ok) {
                                          col = ic;
                                       } else {
                                          if (false)
                                             std::cout << "DEBUG:: failed to get udd for udd_handle " << udd_handle
                                                       << " for user-defined-atom-colour-index" << std::endl;
                                          col = 20;
                                       }
                                    } else {
                                       col = 20; // :-)
                                    }
                                 } else {
                                    // disaster!
                                    col = 20; // (haha)
                                 }
                              } else {
                                 col = 20;
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

   // std::cout << "        atom_colour() returning col " << col << std::endl;
   return col;

}

// This gets called by ca_plus_ligands_rainbow_representation()
//
void
Bond_lines_container::do_Ca_plus_ligands_bonds(atom_selection_container_t SelAtom,
                                               int imol,
                                               coot::protein_geometry *pg,
                                               float min_dist,
                                               float max_dist,
                                               bool draw_missing_loops_flag,
                                               bool do_bonds_to_hydrogens_in,
                                               int model_number) {

   do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
   if (pg) {
      geom = pg;
      have_dictionary = true;
   }
   //do_Ca_plus_ligands_bonds(SelAtom, min_dist, max_dist, coot::COLOUR_BY_ATOM_TYPE);
   do_Ca_plus_ligands_bonds(SelAtom, imol, pg, min_dist, max_dist, draw_missing_loops_flag,
                            coot::COLOUR_BY_CHAIN, do_bonds_to_hydrogens, model_number);
}

void
Bond_lines_container::do_Ca_plus_ligands_and_sidechains_bonds(atom_selection_container_t SelAtom,
                                                              int imol,
                                                              coot::protein_geometry *pg,
                                                              float min_dist_ca, float max_dist_ca,
                                                              float min_dist, float max_dist,
                                                              bool draw_missing_loops_flag,
                                                              bool do_bonds_to_hydrogens_in,
                                                              int model_number) {

   do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
   if (pg) {
      geom = pg;
      have_dictionary = true;
   }
   //do_Ca_plus_ligands_bonds(SelAtom, min_dist, max_dist, coot::COLOUR_BY_ATOM_TYPE);
   do_Ca_plus_ligands_and_sidechains_bonds(SelAtom, imol,
                                           pg, min_dist_ca, max_dist_ca, min_dist, max_dist,
                                           draw_missing_loops_flag,
                                           coot::COLOUR_BY_CHAIN, do_bonds_to_hydrogens,
                                           model_number);
}

void
Bond_lines_container::do_Ca_plus_ligands_bonds(atom_selection_container_t SelAtom,
                                               int imol,
                                               coot::protein_geometry *pg,
                                               float min_dist,
                                               float max_dist,
                                               bool draw_missing_loops_flag,
                                               int atom_colour_type,
                                               bool do_bonds_to_hydrogens_in,
                                               int model_number) {

    if (false)
       std::cout << "---- in do_Ca_plus_ligands_bonds with atom_colour_type "
                 << atom_colour_type << std::endl;

   if (! SelAtom.mol) {
      std::cout << "ERROR:: Caught null mol in do_Ca_plus_ligands_bonds()" << std::endl;
      return;
   }

   do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
   if (pg) {
      geom = pg;
      have_dictionary = true;
   }

   int udd_user_defined_atom_colour_index_handle = SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM,
                                                                             "user-defined-atom-colour-index");

   int n_models = SelAtom.mol->GetNumberOfModels();
   int imodel_start = 1;
   int imodel_end = n_models;
   if (model_number != 0) {
      if (model_number <= n_models) {
         imodel_start = model_number;
         imodel_end   = model_number;
      }
   }

   for (int imodel=imodel_start; imodel<=imodel_end; imodel++) {
      mmdb::Model *model_p = SelAtom.mol->GetModel(imodel);
      if (model_p) {
         try_set_b_factor_scale(SelAtom.mol);
         int istat;
         // udd_has_ca_handle = SelAtom.mol->RegisterUDInteger (mmdb::UDR_RESIDUE, "has CA");
         if (udd_has_ca_handle == -1) {
            udd_has_ca_handle = SelAtom.mol->RegisterUDInteger (mmdb::UDR_RESIDUE, "has CA");
         }

         int nchains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<nchains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  istat = residue_p->PutUDData(udd_has_ca_handle, 0);
                  if (istat == mmdb::UDDATA_WrongUDRType) {
                     std::cout << "ERROR:: mmdb:UDDATA_WrongUDRType in do_Ca_plus_ligands_bonds "
                               << coot::residue_spec_t(residue_p) << " " << udd_has_ca_handle
                               << std::endl;
                  }
               }
            }
         }

         coot::my_atom_colour_map_t acm;
         acm.fill_chain_id_map(SelAtom);
         do_Ca_or_P_bonds_internal(SelAtom, " CA ", acm, min_dist, max_dist, draw_missing_loops_flag,
                                   atom_colour_type, model_number);

         // do_Ca_plus_ligands_bonds has set udd_has_ca_handle on
         // residues that have CAs.  Now let's run through the residues
         // again looking for residues that don't have this handle set.
         // We should do normal bonds for those residues (if they aren't
         // water, of course).
         //
         // Try to bond by dictionary, if not, fall back to previous
         // (add_ligand_bonds()) method. We use a cache of residue names
         // to know if this resiudue type has a dictionary or not.

         std::vector<mmdb::Atom *> ligand_atoms;
         std::vector<std::pair<bool, mmdb::Residue *> > het_residues;
         for (int ichain=0; ichain<nchains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            int ic;
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  if (residue_p->GetUDData(udd_has_ca_handle, ic) == mmdb::UDDATA_Ok) {
                     if (ic == 0) {
                        // Residue was not rendered as CA, needs normal bonds
                        std::string resname(residue_p->name);
                        if (resname != "WAT" && resname != "HOH") {

                           // can we do this residue by dictionary?
                           bool done_by_dict = false; // initially
                           if (have_dictionary) {
                              if (geom->have_at_least_minimal_dictionary_for_residue_type(resname, imol)) {
                                 het_residues.push_back(std::pair<bool, mmdb::Residue *>(true, residue_p));
                                 done_by_dict = true;
                              } else {
                                 std::cout << "Not even minimal for " << resname << std::endl;
                              }
                           }

                           if (! done_by_dict) {
                              // old method
                              int natoms = residue_p->GetNumberOfAtoms();
                              for (int iat=0; iat<natoms; iat++) {
                                 ligand_atoms.push_back(residue_p->GetAtom(iat));
                              }
                           }
                        }
                     }
                  }
               }
            }
         }

         int het_atoms_colour_type = coot::COLOUR_BY_RAINBOW;

         if (atom_colour_type == coot::COLOUR_BY_USER_DEFINED_COLOURS)
            het_atoms_colour_type = coot::COLOUR_BY_USER_DEFINED_COLOURS;

         bool have_udd_atoms = false;
         int udd_bond_handle = -1;
         int udd_atom_index_handle = SelAtom.UDDAtomIndexHandle;
         add_bonds_het_residues(het_residues, SelAtom, imol, het_atoms_colour_type,
                                have_udd_atoms, udd_bond_handle, udd_atom_index_handle,
                                udd_user_defined_atom_colour_index_handle);

         if (! ligand_atoms.empty()) {
            mmdb::PAtom *ligand_atoms_selection = new mmdb::PAtom[ligand_atoms.size()];
            for(unsigned int iat=0; iat<ligand_atoms.size(); iat++) {
               ligand_atoms_selection[iat] = ligand_atoms[iat];
            }
            add_ligand_bonds(SelAtom, imol, ligand_atoms_selection, ligand_atoms.size());
            delete [] ligand_atoms_selection;
         }
      }
   }

}

void
Bond_lines_container::do_Ca_plus_ligands_and_sidechains_bonds(atom_selection_container_t SelAtom,
                                                              int imol,
                                                              coot::protein_geometry *pg,
                                                              float min_dist_ca, float max_dist_ca,
                                                              float min_dist, float max_dist,
                                                              bool draw_missing_loops_flag,
                                                              int atom_colour_type,
                                                              bool do_bonds_to_hydrogens_in,
                                                              int model_number) {

   bool do_rama_markup = false;
   bool show_atoms_as_aniso_flag = false;

   if (! SelAtom.mol) {
      std::cout << "ERROR:: Caught null mol in do_Ca_plus_ligands_and_sidechains_bonds()"
                << std::endl;
      return;
   }

   // first do Ca plus ligand
   do_Ca_plus_ligands_bonds(SelAtom, imol, pg, min_dist_ca, max_dist_ca, draw_missing_loops_flag,
                            atom_colour_type, do_bonds_to_hydrogens_in, model_number);

   // now do normal bonds for CA+sidechain
   // mmmh are the distances correct!?
   // first select all atoms excluding MC (N,C,O)
   atom_selection_container_t asc = SelAtom;

   // Now make a new atom selection that excludes N,C,O by using mmdb::SKEY_XOR
   int newSelectionHandle = asc.mol->NewSelection();
   asc.SelectionHandle = newSelectionHandle;
   std::string solvent_res = "WAT,HOH";

   // We need to select all atoms here first, or crash when going back to all-atom view.
   // Mmmh!?
   asc.mol->SelectAtoms(asc.SelectionHandle, 0, "*",
                        mmdb::ANY_RES, "*",
                        mmdb::ANY_RES, "*",
                        "*", "*", "*", "*");

   asc.mol->Select(asc.SelectionHandle, mmdb::STYPE_ATOM, 0, "*",
                   mmdb::ANY_RES, "*",
                   mmdb::ANY_RES, "*",
                   solvent_res.c_str(), "*", "*", "*",
                   mmdb::SKEY_XOR);

   // this is better than it was - now the main-chain atoms of GLY are excluded.
   // but we still get a star on a GLY.
   asc.mol->Select(asc.SelectionHandle, mmdb::STYPE_ATOM, 0, "*",
                   mmdb::ANY_RES, "*",
                   mmdb::ANY_RES, "*",
                   "*", "[ C  ],[ N  ],[ O  ],[ H  ],[ HA ],[ HA2],[ HA3]", "*", "*",
                   mmdb::SKEY_XOR);

   asc.mol->GetSelIndex(asc.SelectionHandle, asc.atom_selection, asc.n_selected_atoms);

   // for these side chain atoms
   do_colour_by_chain_bonds(asc, true, imol, do_bonds_to_hydrogens_in,
                            draw_missing_loops_flag, 0, false, do_rama_markup, model_number);
   asc.mol->DeleteSelection(asc.SelectionHandle);

}

void
Bond_lines_container::do_normal_bonds_no_water(const atom_selection_container_t &asc_in, int imol,
                                               float min_dist,
                                               float max_dist) {

   atom_selection_container_t asc = asc_in;
   bool do_rama_markup = false;
   bool show_atoms_as_aniso_flag = false;

   // Now make a new atom selection that excludes WAT and HOH by using mmdb::SKEY_XOR
   int newSelectionHandle = asc.mol->NewSelection();
   asc.SelectionHandle = newSelectionHandle;
   std::string solvent_res = "WAT,HOH";

   // We need to select all atoms here first, or crash when going back to all-atom view.
   asc.mol->SelectAtoms(asc.SelectionHandle, 0, "*",
                        mmdb::ANY_RES, "*",
                        mmdb::ANY_RES, "*",
                        "*", "*", "*", "*");

   asc.mol->Select(asc.SelectionHandle, mmdb::STYPE_ATOM, 0, "*",
                   mmdb::ANY_RES, "*",
                   mmdb::ANY_RES, "*",
                   solvent_res.c_str(), "*", "*", "*",
                   mmdb::SKEY_XOR);

   asc.mol->GetSelIndex(asc.SelectionHandle, asc.atom_selection, asc.n_selected_atoms);
   // std::cout << "after water selection: n_selected_atoms: " << asc.n_selected_atoms << std::endl;
   int model_number = 0; // all models
   construct_from_asc(asc, imol, min_dist, max_dist, coot::COLOUR_BY_ATOM_TYPE, 0, model_number,
                      show_atoms_as_aniso_flag, do_rama_markup);
   asc.mol->DeleteSelection(asc.SelectionHandle);
}


int
Bond_lines_container::add_ligand_bonds(const atom_selection_container_t &SelAtom,
                                       int imol,
                                       mmdb::PPAtom ligand_atoms_selection,
                                       int n_ligand_atoms) {

   int ibond = 0;
   atom_selection_container_t asc = SelAtom;
   asc.atom_selection = ligand_atoms_selection;
   asc.n_selected_atoms = n_ligand_atoms;
   // std::cout << "debug:: here in add_ligand_bonds() with " << asc.n_selected_atoms
   // << " ligand atoms" << std::endl;
   int model_number = 0; // all models
   bool do_rama_markup = false;
   bool show_atoms_as_aniso_flag = false;
   construct_from_asc(asc, imol, 0.01, 1.9, coot::COLOUR_BY_ATOM_TYPE, 0, model_number,
                      show_atoms_as_aniso_flag, do_rama_markup);

   return ibond;

}

void
Bond_lines_container::do_colour_sec_struct_bonds(const atom_selection_container_t &asc,
                                                 int imol,
                                                 float min_dist, float max_dist) {

   if (asc.n_selected_atoms > 0) {
      int n_models = asc.mol->GetNumberOfModels();
      for (int imodel=0; imodel<n_models; imodel++) {
         mmdb::Model *model_p = asc.mol->GetModel(1);
         if (model_p)
            model_p->CalcSecStructure(imodel);
      }
      int model_number = 0; // all models
      bool do_rama_markup = false;
      bool show_atoms_as_aniso_flag = false;
      construct_from_asc(asc, imol, 0.01, 1.9, coot::COLOUR_BY_SEC_STRUCT, 0, model_number,
                         show_atoms_as_aniso_flag, do_rama_markup);
   }
}


void
Bond_lines_container::do_Ca_plus_ligands_colour_sec_struct_bonds(const atom_selection_container_t &asc,
                                                                 int imol,
                                                                 coot::protein_geometry *pg,
                                                                 float min_dist, float max_dist,
                                                                 bool draw_missing_loops_flag,
                                                                 bool do_bonds_to_hydrogens_in) {
   do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
   if (asc.n_selected_atoms > 0) {

      for (int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = asc.mol->GetModel(imod);
         if (model_p) {
            int aminoSelHnd = -1;
            model_p->CalcSecStructure(1, aminoSelHnd);
            do_Ca_plus_ligands_bonds(asc, imol, pg, min_dist, max_dist, draw_missing_loops_flag,
                                     coot::COLOUR_BY_SEC_STRUCT, do_bonds_to_hydrogens_in);
         }
      }
   }
}



// Hah!  I'd written this so long ago that I'd forgotten that I'd done
// it.  I have reimplemented symmetry Calphas: addSymmetry_calphas().
//
// This function no longer works.
void
Bond_lines_container::do_symmetry_Ca_bonds(atom_selection_container_t SelAtom,
                                           symm_trans_t symm_trans){

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
   Cell_Translation cell_trans;
   mmdb::PPAtom trans_ca_selection = NULL; // trans_sel(SelAtom, symm_trans);
   int n_ca;
   mmdb::mat44 my_matt;
   int ncontacts;
   long i_contact_group = 1;
   mmdb::Contact *contact = NULL;
   int col;

   // adjust my_matt to the symm_trans:
   //
   int err = SelAtom.mol->GetTMatrix(my_matt, symm_trans.isym(), symm_trans.x(),
                                     symm_trans.y(), symm_trans.z());

   if (err != 0) {
      std::cout << "!!!!!!!!!!!!!! something BAD with mmdb::CMMDBCryst.GetTMatrix"
                << std::endl;
   }

   int selHnd_ca = SelAtom.mol->NewSelection();

   SelAtom.mol->SelectAtoms(selHnd_ca, 0, "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                            "*"," CA ","C","*" );

   SelAtom.mol->GetSelIndex(selHnd_ca, trans_ca_selection, n_ca);

   SelAtom.mol->SeekContacts(trans_ca_selection, n_ca,
                             trans_ca_selection, n_ca,
                             0.01, 5.0, // min, max dist.
                             0,         // in same res also.
                             contact, ncontacts,
                             0, &my_matt, i_contact_group);

   std::cout << "INFO:: Found " << ncontacts/2 << " Ca-Ca links" << std::endl;

   if (ncontacts > 0) {
      for (int i=0; i< ncontacts; i++) {
         if ( contact[i].id2 >  contact[i].id1 ) {

            coot::Cartesian ca_1(trans_ca_selection[ contact[i].id1 ]->x,
                                 trans_ca_selection[ contact[i].id1 ]->y,
                                 trans_ca_selection[ contact[i].id1 ]->z);

            coot::Cartesian ca_2(trans_ca_selection[ contact[i].id2 ]->x,
                                 trans_ca_selection[ contact[i].id2 ]->y,
                                 trans_ca_selection[ contact[i].id2 ]->z);

            // 20171224-PE FIXME
            int iat_1 = -1;
            int iat_2 = -1;
            col = 0; // Carbon = yellow is the first color
            int model_number = trans_ca_selection[contact[i].id2]->GetModelNum();
            addBond(col, ca_1, ca_2, cc, model_number, iat_1, iat_2);
         }
      }
   }
   delete [] contact;

}

void
Bond_lines_container::draw_GA_rings_outer(mmdb::Residue *residue_p, int model_number,
                                          int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                                          int udd_atom_index_handle, int udd_user_defined_atom_colour_index_handle) {

   std::vector<std::string> G_rings_atom_names = {" N9 ", " C8 ", " N7 ", " C5 ", " C4 ", " N3 ", " C2 ", " N1 ", " C6 "};

   // we need a ring atom vector for each of the alt confs
   std::set<std::string> residue_alt_confs;

   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         std::string a(at->altLoc);
         residue_alt_confs.insert(a);
      }
   }

   std::set<std::string>::const_iterator it;

   for(it=residue_alt_confs.begin(); it!=residue_alt_confs.end(); it++) {
      const std::string &alt_loc(*it);
      std::vector<mmdb::Atom *> G_rings_atoms(9, 0);
      unsigned int n_found = 0;
      for (unsigned int i=0; i<G_rings_atom_names.size(); i++) {
         const std::string &G_rings_atom_name = G_rings_atom_names[i];
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string atom_name(at->name);
               if (atom_name == G_rings_atom_name) {
                  std::string a(at->altLoc);
                  if (a == alt_loc || a == "") {
                     G_rings_atoms[i] = at;
                     n_found++;
                  }
               }
            }
         }
      }
      if (n_found == 9) {
         draw_GA_rings(G_rings_atoms, model_number, atom_colour_type, atom_colour_map_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
      } else{
         if (n_found > 0)
            std::cout << "partial trp sidechain (sad face) " << n_found << " " << coot::residue_spec_t(residue_p) << std::endl;
      }
   }
}


void
Bond_lines_container::draw_trp_ring_outer(mmdb::Residue *residue_p, int model_number,
                                          int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                                          int udd_atom_index_handle,
                                          int udd_user_defined_atom_colour_index_handle) {

   std::vector<std::string> trp_rings_atom_names = {" CG ", " CD1", " NE1", " CE2", " CD2", " CE3", " CZ3", " CH2", " CZ2"};
   std::vector<mmdb::Atom *> trp_rings_atoms(9, 0);

   // we need a ring atom vector for each of the alt confs
   std::set<std::string> residue_alt_confs;

   { // don't shadow residue_atoms
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string a(at->altLoc);
            residue_alt_confs.insert(a);
         }
      }
   }


   std::set<std::string>::const_iterator it;
   for(it=residue_alt_confs.begin(); it!=residue_alt_confs.end(); ++it) {
      const std::string &alt_loc(*it);
      unsigned int n_found = 0;
      for (unsigned int i=0; i<trp_rings_atom_names.size(); i++) {
         const std::string &trp_rings_atom_name = trp_rings_atom_names[i];
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string atom_name(at->name);
               if (atom_name == trp_rings_atom_name) {
                  std::string a(at->altLoc);
                  if (a == alt_loc || a == "") {
                     trp_rings_atoms[i] = at;
                     n_found++;
                  }
               }
            }
         }
      }
      if (n_found == 9)
         draw_trp_rings(trp_rings_atoms, model_number, atom_colour_type, atom_colour_map_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
      else
         if (n_found > 0)
         std::cout << "partial trp sidechain (sad face) " << n_found << " " << coot::residue_spec_t(residue_p) << std::endl;
   }

}

void
Bond_lines_container::draw_CUT_ring(mmdb::Residue *residue_p, int model_number,
                                    int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                                    int udd_atom_index_handle, int udd_user_defined_atom_colour_index_handle) {


   std::vector<std::string> ring_atom_names = {" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 "};

   std::string rn = residue_p->GetResName();

   // we need a ring atom vector for each of the alt confs
   std::set<std::string> residue_alt_confs;

   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         std::string a(at->altLoc);
         residue_alt_confs.insert(a);
      }
   }


   std::set<std::string>::const_iterator it;
   for(it=residue_alt_confs.begin(); it!=residue_alt_confs.end(); ++it) {
      const std::string &alt_loc(*it);
      std::vector<mmdb::Atom *> ring_atoms(6,0);
      unsigned int n_found = 0;
      for (unsigned int i=0; i<ring_atom_names.size(); i++) {
         const std::string &ring_atom_name = ring_atom_names[i];
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string atom_name(at->name);
               if (atom_name == ring_atom_name) {
                  std::string a(at->altLoc);
                  if (a == alt_loc || a == "") {
                     ring_atoms[i] = at;
                     n_found++;
                  }
               }
            }
         }
      }
      if (n_found == 6)
         draw_6_membered_ring(rn, ring_atoms, model_number, atom_colour_type, atom_colour_map_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
      else
         if (n_found > 0)
            std::cout << "partial CUT atoms (sad face) " << n_found << " " << coot::residue_spec_t(residue_p) << std::endl;
   }
}

void
Bond_lines_container::draw_phenyl_ring_outer(mmdb::Residue *residue_p, int model_number,
                                             int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                                             int udd_atom_index_handle,
                                             int udd_user_defined_atom_colour_index_handle) {

   std::vector<std::string> residue_alt_confs = coot::util::get_residue_alt_confs(residue_p);
   std::string rn = residue_p->GetResName();
   std::vector<std::string> ring_atom_names = {" CG ", " CD1", " CE1", " CZ ", " CE2", " CD2"};
   std::vector<std::string>::const_iterator it;
   for(it=residue_alt_confs.begin(); it!=residue_alt_confs.end(); ++it) {
      const std::string &alt_loc(*it);
      std::vector<mmdb::Atom *> ring_atoms(6, 0);
      unsigned int n_found = 0;
      for (unsigned int i=0; i<ring_atom_names.size(); i++) {
         const std::string &ring_atom_name = ring_atom_names[i];
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string atom_name(at->name);
               if (atom_name == ring_atom_name) {
                  std::string a(at->altLoc);
                  if (a == alt_loc || a == "") {
                     ring_atoms[i] = at;
                     n_found++;
                  }
               }
            }
         }
      }
      if (n_found == 6)
         draw_6_membered_ring(rn, ring_atoms, model_number, atom_colour_type, atom_colour_map_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
      else
         if (n_found > 0)
            std::cout << "partial ring sidechain (sad face) " << n_found << " " << coot::residue_spec_t(residue_p) << std::endl;
   }

}

// how does this draw a benzodiazepine? Does the ring finder find 7 membered rings?
// - it works very nicely
void
Bond_lines_container::draw_het_group_rings(mmdb::Residue *residue_p,
                                           const std::vector<bonded_quad_atom_names> &bonded_quads,
                                           int model_number, int atom_colour_type,
                                           coot::my_atom_colour_map_t *atom_colour_map_p,
                                           int udd_atom_index_handle, int udd_user_defined_atom_colour_index_handle) {

   // 20200630-PE the current cif file parse make protein residues het groups. Hmm.

   std::vector<bonded_quad_atoms> bqv;

   for (unsigned int i=0; i<bonded_quads.size(); i++) {
      const bonded_quad_atom_names &names = bonded_quads[i];
      bonded_quad_atoms bq;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->name);
            if (names.atom_name(0) == atom_name) bq.atom_1 = at;
            if (names.atom_name(1) == atom_name) bq.atom_2 = at;
            if (names.atom_name(2) == atom_name) bq.atom_3 = at;
            if (names.atom_name(3) == atom_name) bq.atom_4 = at;
         }
      }
      if (bq.filled_p()) {
         if (names.bond_type == bonded_quad_atom_names::SINGLE)
            bq.bond_type = bonded_quad_atoms::SINGLE;
         if (names.bond_type == bonded_quad_atom_names::DOUBLE)
            bq.bond_type = bonded_quad_atoms::DOUBLE;
         bqv.push_back(bq);
      }
   }
   draw_bonded_quad_atoms_rings(bqv, model_number, atom_colour_type, atom_colour_map_p, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
}

void
Bond_lines_container::add_residue_monomer_bonds(const std::map<std::string, std::vector<mmdb::Residue *> > &residue_monomer_map,
                                                int imol,
                                                int model_number,
                                                int atom_colour_type,
                                                coot::my_atom_colour_map_t *atom_colour_map,
                                                int udd_atom_index_handle,
                                                int udd_bond_handle,
                                                int udd_user_defined_atom_colour_index_handle,
                                                int draw_hydrogens_flag,
                                                bool do_goodsell_colour_mode) {


   // std::cout << "in add_residue_monomer_bonds() goodsell_colour_mode: " << do_goodsell_colour_mode
   //           << " atom_colour_type " << atom_colour_type << std::endl;

   bool is_het = false;
   if (!geom) return;

   std::map<std::string, std::vector<mmdb::Residue *> >::const_iterator it;
   class atom_string_bits_t {
   public:
      std::string atom_name;
      std::string ele;
      std::string alt_loc;
      mmdb::Atom *at;
      atom_string_bits_t(const std::string &a, const std::string &e, const std::string &altloc, mmdb::Atom *at_in) : atom_name(a), ele(e), alt_loc(altloc), at(at_in) {}
   };
   // std::map<std::string, std::vector<std::vector<std::tuple<std::string, std::string, std::string> > > > atom_name_ele_map; // and alt-conf
   std::map<std::string, std::vector<std::vector<atom_string_bits_t> > > atom_name_ele_map; // and alt-conf
   for (it=residue_monomer_map.begin(); it!=residue_monomer_map.end(); ++it) {
      const std::string &monomer_name(it->first);
      const std::vector<mmdb::Residue *> &v = it->second;
      atom_name_ele_map[monomer_name].resize(v.size());
      for (unsigned int i=0; i<v.size(); i++) {
         mmdb::Residue *residue_p = v[i];
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string atom_name(at->name);
               std::string ele(at->element);
               std::string alt_loc(at->altLoc);
               atom_string_bits_t asb(atom_name, ele, alt_loc, at);
               atom_name_ele_map[monomer_name][i].push_back(asb);
               if (at->Het)
                  is_het = true;
            }
         }
      }
   }

   std::vector<std::string> phenyl_ring_atom_names = {" CG ", " CD1", " CE1", " CZ ", " CE2", " CD2"};
   std::vector<std::string>   CUT_rings_atom_names = {" N1 ", " C2 ", " N3 ", " C4 ", " C5 ", " C6 "};
   std::vector<std::string>   trp_rings_atom_names = {" CG ", " CD1", " NE1", " CE2", " CD2", " CE3", " CZ3", " CH2", " CZ2"};
   std::vector<std::string>     G_rings_atom_names = {" N9 ", " C8 ", " N7 ", " C5 ", " C4 ", " N3 ", " C2 ", " N1 ", " C6 "}; // "A" is same

   for (it=residue_monomer_map.begin(); it!=residue_monomer_map.end(); ++it) {
      const std::string &monomer_name(it->first);
      const std::string &res_name = monomer_name;
      const std::vector<mmdb::Residue *> &rv = it->second;
      bool phenyl_exception = false;
      bool trp_exception = false;
      bool GA_exception = false;
      bool CUT_exception = false;
      if (res_name == "TRP") trp_exception = true;
      if (res_name == "PHE" || res_name == "TYR") phenyl_exception = true;
      if (res_name == "G"   || res_name == "A"  ) GA_exception = true;
      if (res_name == "DG"  || res_name == "DA" ) GA_exception = true;
      if (res_name == "C"   || res_name == "U" || res_name == "DT") CUT_exception = true;

      if (geom->have_at_least_minimal_dictionary_for_residue_type(monomer_name, imol)) {

         std::pair<bool, coot::dictionary_residue_restraints_t> restraints_pair =
            geom->get_monomer_restraints_at_least_minimal(monomer_name, imol);

         if (restraints_pair.first) {
            const coot::dictionary_residue_restraints_t &restraints = restraints_pair.second;

            // Here outer, doesn't mean the outer ring, it's the function name to draw both (if needed) of the
            // ring bonds (and it has an inner function). Probably better to rename things when it works.

            if (phenyl_exception)
               for (unsigned int i=0; i<rv.size(); i++)
                  draw_phenyl_ring_outer(rv[i], model_number, atom_colour_type, atom_colour_map, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);

            if (trp_exception)
               for (unsigned int i=0; i<rv.size(); i++)
                  draw_trp_ring_outer(rv[i], model_number, atom_colour_type, atom_colour_map, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);

            if (GA_exception) {
               // std::cout << "------------------ got a GA_exception" << std::endl;
               for (unsigned int i=0; i<rv.size(); i++)
                  draw_GA_rings_outer(rv[i], model_number, atom_colour_type, atom_colour_map, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
            }

            if (CUT_exception)
               for (unsigned int i=0; i<rv.size(); i++)
                  draw_CUT_ring(rv[i], model_number, atom_colour_type, atom_colour_map, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);

            std::vector<std::vector<std::string> > rings;
            std::vector<bonded_quad_atom_names> rings_as_bonded_quads_atom_names;

            // mmdb cif parse make protein residues HET at the moment, so let's excude those too
            if (is_het && ! phenyl_exception && ! trp_exception && ! GA_exception && ! CUT_exception) {
               std::vector<std::pair<std::string, std::string> > aromatic_bonds = get_aromatic_bonds(restraints);
               coot::aromatic_graph_t ag(aromatic_bonds);
               rings = ag.ring_list();
               rings_as_bonded_quads_atom_names = ag.bonded_quad_ring_list();
               // put the bond orders into the elments of rings_as_bonded_quads_atom_names
               for (unsigned int i=0; i<rings_as_bonded_quads_atom_names.size(); i++) {
                  bonded_quad_atom_names &quad = rings_as_bonded_quads_atom_names[i];
                  std::string br = restraints.get_bond_type(quad.atom_name(1), quad.atom_name(2));
                  if (false)
                     std::cout << "look up bond type: for " << quad.atom_name(1) << " " << quad.atom_name(2)
                               << " got " << br << std::endl;
                  if (br == "double")
                     quad.bond_type = bonded_quad_atom_names::DOUBLE;
                  else
                     quad.bond_type = bonded_quad_atom_names::SINGLE;
               }
               for (unsigned int i=0; i<rv.size(); i++)
                  if (res_name != "HOH")
                     draw_het_group_rings(rv[i], rings_as_bonded_quads_atom_names, model_number,
                                          atom_colour_type, atom_colour_map, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
            }

            for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
               const coot::dict_bond_restraint_t &br = restraints.bond_restraint[ib];
               std::string bt = restraints.bond_restraint[ib].type();
               std::string dict_atom_name_1 = restraints.bond_restraint[ib].atom_id_1_4c();
               std::string dict_atom_name_2 = restraints.bond_restraint[ib].atom_id_2_4c();

               if (phenyl_exception)
                  if (std::find(phenyl_ring_atom_names.begin(), phenyl_ring_atom_names.end(), dict_atom_name_1) != phenyl_ring_atom_names.end())
                     if (std::find(phenyl_ring_atom_names.begin(), phenyl_ring_atom_names.end(), dict_atom_name_2) != phenyl_ring_atom_names.end())
                        continue;

               if (trp_exception)
                  if (std::find(trp_rings_atom_names.begin(), trp_rings_atom_names.end(), dict_atom_name_1) != trp_rings_atom_names.end())
                     if (std::find(trp_rings_atom_names.begin(), trp_rings_atom_names.end(), dict_atom_name_2) != trp_rings_atom_names.end())
                        continue;

               if (GA_exception)
                  if (std::find(G_rings_atom_names.begin(), G_rings_atom_names.end(), dict_atom_name_1) != G_rings_atom_names.end())
                     if (std::find(G_rings_atom_names.begin(), G_rings_atom_names.end(), dict_atom_name_2) != G_rings_atom_names.end())
                        continue;

               if (CUT_exception)
                  if (std::find(CUT_rings_atom_names.begin(), CUT_rings_atom_names.end(), dict_atom_name_1) != CUT_rings_atom_names.end())
                     if (std::find(CUT_rings_atom_names.begin(), CUT_rings_atom_names.end(), dict_atom_name_2) != CUT_rings_atom_names.end())
                        continue;

               bool in_het_ring = false;
               for (unsigned int j=0; j<rings.size(); j++) {
                  if (std::find(rings[j].begin(), rings[j].end(), dict_atom_name_1) != rings[j].end()) {
                     if (std::find(rings[j].begin(), rings[j].end(), dict_atom_name_2) != rings[j].end()) {
                        in_het_ring = true;
                        break;
                     }
                  }
               }
               if (in_het_ring)
                  continue;

               if (false)
                  std::cout << "If we are here then " << dict_atom_name_1 << " and " << dict_atom_name_2
                            << " are not in a ring" << std::endl;

               for (unsigned int i=0; i<rv.size(); i++) {
                  mmdb::Residue *residue_p = rv[i];

                  // we need to iterate the bonding over all alt confs
                  std::set<std::string> residue_alt_confs;
                  mmdb::Atom **residue_atoms = 0;
                  int n_residue_atoms = 0;
                  residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                  for (int iat=0; iat<n_residue_atoms; iat++) {
                     mmdb::Atom *at = residue_atoms[iat];
                     if (! at->isTer()) {
                        std::string a(at->altLoc);
                        residue_alt_confs.insert(a);
                     }
                  }

                  const std::vector<atom_string_bits_t> &atom_name_ele_vec = atom_name_ele_map[monomer_name][i];

                  std::set<std::string>::const_iterator it;
                  for(it=residue_alt_confs.begin(); it!=residue_alt_confs.end(); ++it) {

                     const std::string &alt_loc(*it);

                     // find the atom indices in the first atom name vector then try to use
                     // those indices to find the same atom names in the other vectors
                     int iat = -1;
                     int jat = -1;
                     mmdb::Atom *bond_atom_1 = 0;
                     mmdb::Atom *bond_atom_2 = 0;
                     // tuple: atom_name, ele, alt_loc
                     for (unsigned int ii=0; ii<atom_name_ele_vec.size(); ii++) {
                        if (! bond_atom_1)
                           if (atom_name_ele_vec[ii].atom_name == dict_atom_name_1)
                              if (atom_name_ele_vec[ii].alt_loc == alt_loc || atom_name_ele_vec[ii].alt_loc.empty())
                                 if (! residue_atoms[ii]->isTer()) {
                                    bond_atom_1 = atom_name_ele_vec[ii].at;
                                    iat = ii;
                                 }
                        if (! bond_atom_2)
                           if (atom_name_ele_vec[ii].atom_name == dict_atom_name_2)
                              if (atom_name_ele_vec[ii].alt_loc == alt_loc || atom_name_ele_vec[ii].alt_loc.empty())
                                 if (! residue_atoms[ii]->isTer()) {
                                    bond_atom_2 = atom_name_ele_vec[ii].at;
                                    jat = ii;
                                 }
                     }
                     if (bond_atom_1) {
                        if (bond_atom_2) {

                           residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                           int ierr = 0;
                           int atom_idx_1 = -1;  // atom index in the asc atom selection
                           int atom_idx_2 = -1;
                           ierr = residue_atoms[iat]->GetUDData(udd_atom_index_handle, atom_idx_1);
                           if (ierr != mmdb::UDDATA_Ok)
                              std::cout << "ERROR:: add_residue_monomer_bonds() UDD Index error A " << udd_atom_index_handle << " "
                                        << coot::atom_spec_t(residue_atoms[iat]) << std::endl;
                           ierr = residue_atoms[jat]->GetUDData(udd_atom_index_handle, atom_idx_2);
                           if (ierr != mmdb::UDDATA_Ok)
                              std::cout << "ERROR:: add_residue_monomer_bonds() UDD Index error B " << udd_atom_index_handle << " "
                                        << coot::atom_spec_t(residue_atoms[jat]) << std::endl;
                           if (true) {
                              if (ierr != mmdb::UDDATA_Ok) {
                              }
                           }

                           mmdb::Atom *atom_p_1 = bond_atom_1;
                           mmdb::Atom *atom_p_2 = bond_atom_2;
                           std::string element_1(atom_p_1->element);
                           std::string element_2(atom_p_2->element);
                           coot::Cartesian p1(atom_p_1->x, atom_p_1->y, atom_p_1->z);
                           coot::Cartesian p2(atom_p_2->x, atom_p_2->y, atom_p_2->z);
                           int iat_1_atom_index = -1;
                           int iat_2_atom_index = -1;

                           if (false)
                              std::cout << "debug/diag:: " << dict_atom_name_1 << " " << dict_atom_name_2 << " " << iat << " " << jat
                                        << atom_idx_1 << " " << atom_idx_2 << " "
                                        << coot::atom_spec_t(atom_p_1) << " "
                                        << coot::atom_spec_t(atom_p_2) << " "
                                        << std::endl;

                           // getting the UDD data again? A mistake I think.
                           residue_atoms[iat]->GetUDData(udd_atom_index_handle, iat_1_atom_index);
                           residue_atoms[jat]->GetUDData(udd_atom_index_handle, iat_2_atom_index);

                           if (element_1 != element_2) {
                              if (!is_hydrogen(element_1) && !is_hydrogen(element_2)) {

                                 if (bt == "double") {
                                    if (br.aromaticity == coot::dict_bond_restraint_t::AROMATIC) {

                                       graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                       add_half_bonds(p1, p2, atom_p_1, atom_p_2,
                                                      cc, model_number, iat_1_atom_index, iat_2_atom_index,
                                                      atom_colour_type, udd_user_defined_atom_colour_index_handle,
                                                      atom_colour_map, false, false);

                                    } else {
                                       add_double_bond(imol, model_number, iat, jat, residue_atoms, n_residue_atoms,
                                                       atom_colour_type, atom_colour_map,
                                                       udd_atom_index_handle,
                                                       udd_user_defined_atom_colour_index_handle,
                                                       restraints.bond_restraint);
                                    }
                                 } else {
                                    if (bt == "deloc") {
                                       bool is_deloc = true;
                                       add_double_bond(imol, model_number, iat, jat, residue_atoms, n_residue_atoms,
                                                       atom_colour_type, atom_colour_map,
                                                       udd_atom_index_handle,
                                                       udd_user_defined_atom_colour_index_handle,
                                                       restraints.bond_restraint, is_deloc);
                                    } else {

                                       if (bt == "triple") {
                                          add_triple_bond(imol, model_number, iat, jat, residue_atoms, n_residue_atoms,
                                                          atom_colour_type, atom_colour_map,
                                                          udd_atom_index_handle,
                                                          udd_user_defined_atom_colour_index_handle,
                                                          restraints.bond_restraint);
                                       } else {
                                          // could be "metal"
                                          graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                          add_half_bonds(p1, p2,
                                                         residue_atoms[iat],
                                                         residue_atoms[jat],
                                                         cc, model_number,
                                                         iat_1_atom_index, iat_2_atom_index,
                                                         atom_colour_type,
                                                         udd_user_defined_atom_colour_index_handle,
                                                         atom_colour_map, false, false);
                                       }
                                    }
                                 }
                              } else {
                                 if (do_bonds_to_hydrogens) {
                                    if (res_name == "HOH" || res_name == "DOD") {
                                       graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                       add_half_bonds(p1, p2,
                                                      residue_atoms[iat],
                                                      residue_atoms[jat],
                                                      cc, model_number,
                                                      iat_1_atom_index, iat_2_atom_index,
                                                      atom_colour_type,
                                                      udd_user_defined_atom_colour_index_handle,
                                                      atom_colour_map, false, false);
                                    } else {
                                       graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                       addBond(HYDROGEN_GREY_BOND, p1, p2, cc, model_number, iat_1_atom_index, iat_2_atom_index); // 20171224-PE correct indices?w
                                    }
                                 }
                              }

                           } else {

                              // Bonded to an atom of the same element.
                              // add_double_bond() uses residue-atom indexing.
                              // addBond uses all-molecule atom indexing.
                              int col = atom_colour(residue_atoms[iat], atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map);
                              if (bt == "double") {

                                 // if (br.aromaticity == coot::dict_bond_restraint_t::AROMATIC)
                                    // std::cout << "double and aromatic and same-ele!" << std::endl;

                                 add_double_bond(imol, model_number, iat, jat, residue_atoms, n_residue_atoms,
                                                 atom_colour_type, atom_colour_map,
                                                 udd_atom_index_handle,
                                                 udd_user_defined_atom_colour_index_handle,
                                                 restraints.bond_restraint);
                              } else {
                                 if (bt == "deloc") {
                                    bool is_deloc = 1;
                                    add_double_bond(imol, model_number, iat, jat, residue_atoms, n_residue_atoms,
                                                    atom_colour_type, atom_colour_map,
                                                    udd_atom_index_handle,
                                                    udd_user_defined_atom_colour_index_handle,
                                                    restraints.bond_restraint, is_deloc);
                                 } else {
                                    if (bt == "triple") {
                                       add_triple_bond(imol, model_number, iat, jat, residue_atoms, n_residue_atoms,
                                                       atom_colour_type, atom_colour_map,
                                                       udd_atom_index_handle,
                                                       udd_user_defined_atom_colour_index_handle,
                                                       restraints.bond_restraint);
                                    } else {
                                       addBond(col, p1, p2, graphics_line_t::SINGLE, model_number,
                                               iat_1_atom_index, iat_2_atom_index, false, false);
                                    }
                                 }
                              }
                           }

                           bool have_udd_bond_handle = false; // not sure what this does a the moment.
                           // I need to go back to a previous version to find
                           // how this is used.
                           if (have_udd_bond_handle) {
                              residue_atoms[iat]->PutUDData(udd_bond_handle, graphical_bonds_container::BONDED_WITH_HETATM_BOND);
                              residue_atoms[jat]->PutUDData(udd_bond_handle, graphical_bonds_container::BONDED_WITH_HETATM_BOND);
                           }
                        }
                     }
                  }
               }
            }
         }

      } else {
         std::cout << "DEBUG:: add_residue_monomer_bonds(): Bond this type by distance " << monomer_name << std::endl;
      }
   }
}

void
Bond_lines_container::bond_by_distance(const atom_selection_container_t &asc, int imol,
                                       std::vector<mmdb::Residue *> &residues,
                                       bool have_udd_atoms, int udd_found_bond_handle) {

   for (unsigned int i=0; i<residues.size(); i++) {
      mmdb::Residue *residue_p = residues[i];
      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      float max_dist = 2.0;
      float min_dist = 0.01;
      int atom_colour_type = coot::COLOUR_BY_CHAIN_C_ONLY;

      residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
      construct_from_atom_selection(asc,
                                    residue_atoms, nResidueAtoms,
                                    residue_atoms, nResidueAtoms,
                                    imol,
                                    min_dist, max_dist, atom_colour_type,
                                    0, have_udd_atoms, udd_found_bond_handle);
   }
}

void
Bond_lines_container::do_colour_by_dictionary_and_by_chain_bonds_carbons_only(const atom_selection_container_t &asc,
                                                                              int imol,
                                                                              int model_number,
                                                                              int draw_hydrogens_flag,
                                                                              bool draw_missing_loops_flag,
                                                                              bool do_goodsell_colour_mode,
                                                                              bool do_rotamer_markup) {

   bool make_stars = false;

   //  timer here.

   // 1) waters and metals
   // 2) I need to use some sort of UDD so that the inner short bonds are not displayed during refinement.

   // Monomer bonds.

   // std::cout << "in do_colour_by_dictionary_and_by_chain_bonds_carbons_only() " << std::endl;

   // std::cout << "in do_colour_by_dictionary_and_by_chain_bonds_carbons_only non-drawn bonds size "
   //           << no_bonds_to_these_atoms.size() << std::endl;

   coot::my_atom_colour_map_t atom_colour_map; // colour map for chain indexing
   atom_colour_map.fill_chain_id_map(asc);

   int n_chains = 0;
   std::map<std::string, std::vector<mmdb::Residue *> > residue_monomer_map;
   for (int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
      if (model_number != 0)
         if (imod != model_number)
            continue;
      mmdb::Model *model_p = asc.mol->GetModel(imod);
      if (model_p) {
         n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               std::string res_name = residue_p->GetResName();
               residue_monomer_map[res_name].push_back(residue_p);
            }
         }
      }
   }
   // resize the bonds vector now that we know the number of bonds
   int bonds_size = bonds.size(); // change type
   if (bonds_size < (50+n_chains))
      bonds.resize(50+n_chains);

   int udd_atom_index_handle = asc.UDDAtomIndexHandle;
   int atom_colour_type = coot::COLOUR_BY_CHAIN_C_ONLY;
   int udd_bond_handle = -1; // surely this is not right.

   int udd_user_defined_atom_colour_index_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   // std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@ here with udd_user_defined_atom_colour_index_handle "
   // << udd_user_defined_atom_colour_index_handle
   // << std::endl;

   if (false) {
      int imod = 1;
      mmdb::Model *model_p = asc.mol->GetModel(imod);
      if (model_p) {
         n_chains = model_p->GetNumberOfChains();
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
                        int idx;
                        if (at->GetUDData(udd_user_defined_atom_colour_index_handle, idx) == mmdb::UDDATA_Ok) {
                           std::cout << "atom " << at << " has colour index " << idx << std::endl;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   {

      // add a timer one day

      // 20240713-PE Harry-bonding:
      //
      if (geom) { // should be - right?
         std::vector<mmdb::Residue *> residues_with_no_dictionary; // fill this
         // cf Bond_lines_container::add_atom_centres()
         std::map<std::string, std::vector<mmdb::Residue *> >::const_iterator it;
         for (it=residue_monomer_map.begin(); it!=residue_monomer_map.end(); ++it) {
            const std::string &res_type = it->first;
            bool s = geom->have_at_least_minimal_dictionary_for_residue_type(res_type, imol);
            if (! s) {
               const std::vector<mmdb::Residue *> &v = it->second;
               residues_with_no_dictionary.insert(residues_with_no_dictionary.end(), v.begin(), v.end());
            }
         }
         if (! residues_with_no_dictionary.empty()) {
            int udd_found_bond_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "found bond");
            bool have_udd_atoms = true;
            bond_by_distance(asc, imol, residues_with_no_dictionary, have_udd_atoms, udd_found_bond_handle);
         }
      }
   }

   if (do_goodsell_colour_mode)
      atom_colour_type = coot::COLOUR_BY_CHAIN_GOODSELL;

   for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {

      if (imod != model_number)
         if (model_number != 0)
            continue;

      add_residue_monomer_bonds(residue_monomer_map, imol, imod,
                                atom_colour_type, &atom_colour_map, udd_atom_index_handle, udd_bond_handle,
                                udd_user_defined_atom_colour_index_handle,
                                draw_hydrogens_flag, do_goodsell_colour_mode);
   }

   std::cout << "DEBUG::do_colour_by_dictionary_and_by_chain_bonds_carbons_only() calling add_polymer_bonds()" << std::endl;
   add_polymer_bonds(asc, model_number, atom_colour_type, &atom_colour_map,  // draw_missing_loops_flag not used here
                     draw_hydrogens_flag, do_goodsell_colour_mode);

   int udd_fixed_during_refinement_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "FixedDuringRefinement");
   if (draw_missing_loops_flag)
      atom_selection_missing_loops(asc, udd_atom_index_handle, udd_fixed_during_refinement_handle);

   // -------- metals and waters

   // I am using RegisterUDInteger here (emphasis on the "Register").
   // I surely didn't mean that? Because the next loop won't find any atoms that
   // are registered.
   //
   int udd_found_bond_handle = asc.mol->RegisterUDInteger(mmdb::UDR_ATOM,"found bond");// Register! Not get.
   int ic = -1;
   for (int iat=0; iat<asc.n_selected_atoms; iat++) {
      mmdb::Atom *at = asc.atom_selection[iat];
      if (at->GetUDData(udd_found_bond_handle, ic) == mmdb::UDDATA_Ok) {
         if (ic == graphical_bonds_container::NO_BOND) {
            mmdb::Residue *residue_p = at->residue;
            std::string res_name(residue_p->GetResName());
            if (res_name == "HOH")
               if (! do_sticks_for_waters)
                  continue;

            int model_number = residue_p->GetModelNum();

            // extract this to own function
            if (make_stars) {
               //  std::cout << "making stars for " << at << " " << coot::atom_spec_t(at) << std::endl;
               float star_size = 0.3;
               coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
               coot::Cartesian small_vec_y(0.0, star_size, 0.0);
               coot::Cartesian small_vec_z(0.0, 0.0, star_size);
               int col = atom_colour(at, atom_colour_type, udd_user_defined_atom_colour_index_handle, &atom_colour_map);
               coot::Cartesian atom_pos(at->x, at->y, at->z);

               int iat_1 = -1;
               int udd_status_1 = at->GetUDData(udd_atom_index_handle, iat_1);

               if (udd_status_1 == mmdb::UDDATA_Ok) {
                  graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                  addBond(col, atom_pos+small_vec_x, atom_pos-small_vec_x, cc, model_number, iat_1, iat_1, true, true);
                  addBond(col, atom_pos+small_vec_y, atom_pos-small_vec_y, cc, model_number, iat_1, iat_1, true, true);
                  addBond(col, atom_pos+small_vec_z, atom_pos-small_vec_z, cc, model_number, iat_1, iat_1, true, true);
               }
            }
         }
      }
   }

   for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
      do_disulphide_bonds(asc, imod);
      add_cis_peptide_markup(asc, imod);
      construct_from_model_links(asc.mol->GetModel(imod), udd_atom_index_handle,
                                 udd_user_defined_atom_colour_index_handle, atom_colour_type);
   }

   add_zero_occ_spots(asc);
   add_deuterium_spots(asc);

   if (do_rotamer_markup)
      add_rotamer_goodness_markup(asc);

   if (do_goodsell_colour_mode)
      atom_colour_type = coot::COLOUR_BY_CHAIN_GOODSELL;

   add_atom_centres(imol, asc, atom_colour_type, model_number, &atom_colour_map);

}

void
Bond_lines_container::add_polymer_bonds(const atom_selection_container_t &asc,
                                        int model_number,
                                        int atom_colour_type,
                                        coot::my_atom_colour_map_t *atom_colour_map_p,
                                        int draw_hydrogens_flag,
                                        bool do_goodsell_colour_mode) {

   std::cout << "DEBUG:: add_polymer_bonds() " << model_number << std::endl;

   add_peptide_bonds(       asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, do_goodsell_colour_mode);
   add_phosphodiester_bonds(asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, do_goodsell_colour_mode);
   add_carbohydrate_bonds(  asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, do_goodsell_colour_mode);
}

void
Bond_lines_container::add_polymer_bonds_generic(const atom_selection_container_t &asc,
                                                int model_number,
                                                int atom_colour_type,
                                                coot::my_atom_colour_map_t *atom_colour_map_p,
                                                int draw_hydrogens_flag,
                                                const std::string &res_1_atom_name, // in "res1"
                                                const std::string &res_2_atom_name, // in "res2"
                                                bool allow_het_group_linking,
                                                bool do_goodsell_colour_mode) {
   if (false)
      std::cout << "DEBUG:: in add_polymer_bonds_generic(): with no_bonds_to_these_atoms size "
                << no_bonds_to_these_atoms.size() << std::endl;

   int udd_user_defined_atom_colour_index_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   // hetgroups should generally not be linked to each other, unless allow_het_group_link_bonds
   // which allows linking (of carbohydrates). In that case, we need a distance sanity check.

   for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
      if (model_number != 0)
         if (imod != model_number)
            continue;
      mmdb::Model *model_p = asc.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<(nres-1); ires++) {
               mmdb::Residue *residue_this_p = chain_p->GetResidue(ires);
               mmdb::Residue *residue_next_p = chain_p->GetResidue(ires+1);
               int n_this_res_atoms = residue_this_p->GetNumberOfAtoms();
               int n_next_res_atoms = residue_next_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_this_res_atoms; iat++) {
                  mmdb::Atom *at_1 = residue_this_p->GetAtom(iat);
                  if (! at_1->isTer()) {
                     std::string at_1_name(at_1->name);
                     if (at_1_name == res_1_atom_name) {
                        for (int jat=0; jat<n_next_res_atoms; jat++) {
                           mmdb::Atom *at_2 = residue_next_p->GetAtom(jat);
                           if (! at_2->isTer()) {
                              std::string at_2_name(at_2->name);
                              if (at_2_name == res_2_atom_name) {
                                 std::string alt_conf_1(at_1->altLoc);
                                 std::string alt_conf_2(at_2->altLoc);
                                 if (alt_conf_1 == alt_conf_2 || alt_conf_1 == "" || alt_conf_2 == "") {

                                    coot::Cartesian atom_1_pos(at_1->x, at_1->y, at_1->z);
                                    coot::Cartesian atom_2_pos(at_2->x, at_2->y, at_2->z);

                                    int res_no_delta = residue_next_p->GetSeqNum() - residue_this_p->GetSeqNum();
                                    bool do_it = true;

                                    if (at_1->Het && ! allow_het_group_linking) do_it = false;
                                    if (at_2->Het && ! allow_het_group_linking) do_it = false;

                                    if (do_it) {
                                       if (res_no_delta > 1 || (at_1->Het && allow_het_group_linking)) {
                                          float dd = coot::Cartesian::lengthsq(atom_1_pos, atom_2_pos);
                                          if (dd > 9.0)
                                             do_it = false;
                                       }
                                    }

                                    if (do_it) {
                                       int atom_index_1 = -1;
                                       int atom_index_2 = -1;
                                       int ierr_1 = at_1->GetUDData(asc.UDDAtomIndexHandle, atom_index_1);
                                       int ierr_2 = at_2->GetUDData(asc.UDDAtomIndexHandle, atom_index_2);

                                       if (ierr_1 == mmdb::UDDATA_Ok)  {
                                          if (ierr_2 == mmdb::UDDATA_Ok) {

                                             graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                                             add_half_bonds(atom_1_pos, atom_2_pos,
                                                            at_1, at_2,
                                                            cc, imod,
                                                            atom_index_1, atom_index_2,
                                                            atom_colour_type,
                                                            udd_user_defined_atom_colour_index_handle,
                                                            atom_colour_map_p, false, false);
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
         }
      }
   }

}

void
Bond_lines_container::add_peptide_bonds(const atom_selection_container_t &asc,
                                        int model_number,
                                        int atom_colour_type,
                                        coot::my_atom_colour_map_t *atom_colour_map_p,
                                        int draw_hydrogens_flag,
                                        bool do_goodsell_colour_mode) {

   bool allow_het_group_linking = true; // 20240519-PE was false. Thierry Fischmann Marc Elsliger complained.
                                        // In the same week - even though this code was written 4 years ago.
                                        // Funny.
   add_polymer_bonds_generic(asc, model_number, atom_colour_type, atom_colour_map_p,
                             draw_hydrogens_flag, " C  ", " N  ", allow_het_group_linking,
                             do_goodsell_colour_mode);

}

void
Bond_lines_container::add_phosphodiester_bonds(const atom_selection_container_t &asc,
                                               int model_number,
                                               int atom_colour_type,
                                               coot::my_atom_colour_map_t *atom_colour_map_p,
                                               int draw_hydrogens_flag,
                                               bool do_goodsell_colour_mode) {

   add_polymer_bonds_generic(asc, model_number, atom_colour_type, atom_colour_map_p,
                             draw_hydrogens_flag, " O3'", " P  ", false, do_goodsell_colour_mode);

}

void
Bond_lines_container::add_carbohydrate_bonds(const atom_selection_container_t &asc,
                                             int model_number,
                                             int atom_colour_type,
                                             coot::my_atom_colour_map_t *atom_colour_map_p,
                                             int draw_hydrogens_flag,
                                             bool do_goodsell_colour_mode) {

   bool gm = do_goodsell_colour_mode;
   add_polymer_bonds_generic(asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, " O1 ", " C1 ", true, gm);
   add_polymer_bonds_generic(asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, " O2 ", " C1 ", true, gm);
   add_polymer_bonds_generic(asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, " O3 ", " C1 ", true, gm);
   add_polymer_bonds_generic(asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, " O4 ", " C1 ", true, gm);
   add_polymer_bonds_generic(asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, " O5 ", " C1 ", true, gm);
   add_polymer_bonds_generic(asc, model_number, atom_colour_type, atom_colour_map_p, draw_hydrogens_flag, " O6 ", " C1 ", true, gm);

}

// main-chain colour as in colour-by-chain
// and either orange (hydrophobic) or blue

void
Bond_lines_container::do_colour_by_hydrophobic_side_chains(const atom_selection_container_t &asc,
                                                           int imol,
                                                           bool draw_missing_loops_flag,
                                                           int draw_hydrogens_flag) {

   // somthing here
   int atom_colour_type =  coot::COLOUR_BY_HYDROPHOBIC_SIDE_CHAIN;
}



void
Bond_lines_container::do_colour_by_dictionary_and_by_chain_bonds(const atom_selection_container_t &asc,
                                                                 int imol,
                                                                 int imodel,
                                                                 int draw_hydrogens_flag,
                                                                 bool draw_missing_loops_flag,
                                                                 short int change_c_only_flag,
                                                                 bool do_goodsell_colour_mode,
                                                                 bool do_rotamer_markup) {

   if (change_c_only_flag) {
      do_colour_by_dictionary_and_by_chain_bonds_carbons_only(asc, imol, imodel,
                                                              draw_hydrogens_flag, draw_missing_loops_flag,
                                                              do_goodsell_colour_mode, do_rotamer_markup);
   } else {
      bool use_asc_atom_selection_flag = true; // 20220226-PE I don't know
      do_colour_by_chain_bonds(asc, use_asc_atom_selection_flag, imol, draw_hydrogens_flag,
                               draw_missing_loops_flag, false, false, do_rotamer_markup);
   }
}

// I add Colour by Segment at last (28 Oct 2003)
//
// 2001130 Now we pass the change_c_only_flag at the request of Phil.
// So if change_c_only_flag is 0, then we have a single colour for
// each chain.
//
// If change_c_only_flag is 1, then we want only Carbons of the chain
// to be coloured by chain.  The other atoms should be coloured by
// atom type.
//
void
Bond_lines_container::do_colour_by_chain_bonds(const atom_selection_container_t &asc,
                                               bool use_asc_atom_selection_flag,
                                               int imol,
                                               int draw_hydrogens_flag,
                                               bool draw_missing_loops_flag,
                                               short int change_c_only_flag,
                                               bool do_goodsell_colour_mode,
                                               bool do_ramachandran_markup,
                                               int model_number) {

   if (false)
      std::cout << "............. in do_colour_by_chain_bonds() goodsell: " << do_goodsell_colour_mode
                << " change_c_only: " << change_c_only_flag << " model_number " << model_number << std::endl;

   coot::my_atom_colour_map_t atom_colour_map;
   atom_colour_map.fill_chain_id_map(asc);

   if (change_c_only_flag) {
      do_colour_by_dictionary_and_by_chain_bonds(asc,
                                                 imol,
                                                 model_number,
                                                 draw_hydrogens_flag,
                                                 draw_missing_loops_flag,
                                                 change_c_only_flag,
                                                 do_goodsell_colour_mode,
                                                 do_ramachandran_markup);

      return;
   }

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
   if (change_c_only_flag) {
      // old code path
      //int atom_colour_type = coot::COLOUR_BY_CHAIN_C_ONLY;
      //if (do_goodsell_colour_mode)
      // atom_colour_type = coot::COLOUR_BY_CHAIN_GOODSELL;
      // do_colour_by_chain_bonds_carbons_only(asc, imol, draw_missing_loops_flag, atom_colour_type, draw_hydrogens_flag);
      // return;
      if (do_goodsell_colour_mode) {
         int atom_colour_type = coot::COLOUR_BY_CHAIN_GOODSELL;
         do_colour_by_chain_bonds_carbons_only(asc, imol, draw_missing_loops_flag, atom_colour_type, draw_hydrogens_flag);
         return;
      }
   } else {
      if (do_goodsell_colour_mode) {
         int atom_colour_type = coot::COLOUR_BY_CHAIN_GOODSELL;
         do_colour_by_chain_bonds_carbons_only(asc, imol, draw_missing_loops_flag, atom_colour_type, draw_hydrogens_flag);
         return;
      }
   }

   std::cout << "DEBUG:: ****************************************************************************************** got here A" << std::endl;

   // OK, now we are on the "single-colour per chain" path

   mmdb::Contact *contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   float max_dist = 1.9;
   float min_dist = 0.01; // As in the constructor
                          // Bond_lines_container::Bond_lines_container(const
                          // atom_selection_container_t &SelAtom, int
                          // do_disulphide_bonds_in, int
                          // do_bonds_to_hydrogens_in)

   // matrix stuff
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;

   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;

   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   int col = 0; // atom (segment) colour

   int n_models = asc.mol->GetNumberOfModels();

   int udd_atom_index_handle = asc.UDDAtomIndexHandle;

   int udd_user_defined_atom_colour_index_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   // default to all models:
   // int imodel_start = 1;
   // int imodel_end = n_models;
   // override if we were passed a specific model_number
   // if (model_number != 0) {
   //    if (model_number <= n_models) {
   //       imodel_start = model_number;
   //       imodel_end   = model_number;
   //    }
   // }

   std::cout << "DEBUG:: ****************************************************************************************** got here B" << std::endl;

   for (int imodel=1; imodel<=n_models; imodel++) {

      std::cout << "DEBUG:: looping iwth imodel " << imodel << std::endl;

      if (imodel != model_number)
         if (model_number != 0)
            continue;

      std::cout << "DEBUG:: ... actually doing with model " << imodel << std::endl;

      mmdb::PPAtom atom_selection = 0;
      int n_selected_atoms = 0;
      contact = NULL;

      int SelectionHandle = asc.mol->NewSelection();

      if (use_asc_atom_selection_flag) {

         atom_selection = asc.atom_selection;
         n_selected_atoms = asc.n_selected_atoms;

      } else {

         // make a new atom selection, based on the model.
         asc.mol->SelectAtoms (SelectionHandle, imodel, "*",
                               mmdb::ANY_RES, // starting resno, an int
                               "*", // any insertion code
                               mmdb::ANY_RES, // ending resno
                               "*", // ending insertion code
                               "*", // any residue name
                               "*", // atom name
                               "*", // elements
                               "*"  // alt loc.
                               );

         asc.mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);
      }

      asc.mol->SeekContacts(atom_selection, n_selected_atoms,
                            atom_selection, n_selected_atoms,
                            min_dist, max_dist, // min, max distances
                            0,        // seqDist 0 -> in same res also
                            contact, ncontacts,
                            0, &my_matt, i_contact_group);

      int res1, res2;

      // Now, let's not forget that some atoms don't have contacts, so
      // whenever we find a contact for an atom, we mark it with
      // UserDefinedData "found bond".
      //
      int uddHnd = asc.mol->RegisterUDInteger (mmdb::UDR_ATOM,"found bond");
      if (uddHnd<0)  {
         std::cout << " atom bonding registration failed.\n";
      } else {
         for (int i=0; i<n_selected_atoms; i++)
            atom_selection[i]->PutUDData(uddHnd, graphical_bonds_container::NO_BOND);
      }

      if (ncontacts > 0) {

         mmdb::Atom *at1 = 0;
         mmdb::Atom *at2 = 0;
         std::string element1;
         std::string element2;

         for (int i=0; i< ncontacts; i++) {
            if (contact[i].id2 > contact[i].id1) {

               int iat_1 = contact[i].id1;
               int iat_2 = contact[i].id2;

               at1 = atom_selection[ contact[i].id1 ];
               at2 = atom_selection[ contact[i].id2 ];

               res1 = at1->GetSeqNum();
               res2 = at2->GetSeqNum();

               if (abs(res1 - res2) < 2) {

                  std::string chain_id_1(at1->GetChainID());
                  std::string chain_id_2(at2->GetChainID());
                  col = atom_colour_map.index_for_chain(chain_id_1);

                  if (chain_id_1 == chain_id_2) {

                     element1 = at1->element;
                     element2 = at2->element;
                     if ( (draw_hydrogens_flag == 1) ||
                          // (element1 != " H" && element1 != " D" &&
                          //  element2 != " H" && element2 != " D") ) {
                          (! is_hydrogen(element1) && ! is_hydrogen(element2))) {

                        coot::Cartesian atom_1(at1->x, at1->y, at1->z);
                        coot::Cartesian atom_2(at2->x, at2->y, at2->z);

                        // alternate location test
                        //
                        std::string aloc_1(at1->altLoc);
                        std::string aloc_2(at2->altLoc);
                        //
                        if (aloc_1 == "" || aloc_2 == "" || aloc_1 == aloc_2) {
                           bonds_size_colour_check(col);
                           addBond(col, atom_1, atom_2, cc, imodel, iat_1, iat_2);

                           if (uddHnd>=0) {
                              at1->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                              at2->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                           }
                        }
                     } else {
                        // It was a hydrogen (or bonded to Hydrogen).
                        // Mark it as bonded (we don't want to see single
                        // unbonded (stared) hydorgens.
                        if (uddHnd>=0) {
                           at1->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                           at2->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                        }
                     }
                  }
               }
            }
         }
         delete [] contact;
      }

      if (! (uddHnd>=0)) {
         std::cout << "ERROR:: do_colour_by_chain_bonds() bad uddHnd"
                   << std::endl;
      } else {

         float star_size = 0.22;
         // for atoms with no neighbour (contacts):
         coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
         coot::Cartesian small_vec_y(0.0, star_size, 0.0);
         coot::Cartesian small_vec_z(0.0, 0.0, star_size);

         int atom_colour_type = coot::COLOUR_BY_CHAIN;

         int ic; // changed by reference;
         int col;
         for (int i=0; i<n_selected_atoms; i++) {
            if ( atom_selection[i]->GetUDData(uddHnd, ic) != mmdb::UDDATA_Ok ) {
               std::cout << "ERROR:: do_colour_by_chain_bonds() failed to get ic bond info"
                         << std::endl;
            } else {
               if ((ic == graphical_bonds_container::NO_BOND) ||
                   (!strcmp(atom_selection[i]->element, " S")) ||
                   (!strcmp(atom_selection[i]->element, "SE")) ||
                   (!strcmp(atom_selection[i]->element, " P"))) {

                  std::string segid(atom_selection[i]->GetChainID());
                  col = atom_colour_map.index_for_chain(segid);

                  // no contact found or was Sulphur, or Phosphor

                  // So, was this a seleno-methione?
                  //
                  mmdb::Residue *atom_residue_p = atom_selection[i]->residue;
                  if (atom_residue_p) {
                     std::string resname = atom_selection[i]->GetResName();
                     if (resname == "MSE" || resname == "MET"
                         || resname == "MSO" || resname == "CYS") {
                        // handle_MET_or_MSE_case(atom_selection[i], uddHnd, udd_atom_index_handle, col);
                     } else {
                        std::string ele = atom_selection[i]->element;
                        if (ele == "CL" || ele == "BR" || ele == " S" ||  ele == " I"
                            || ele == "Cl" || ele == "Br" || ele == "MO"
                            || ele == "PT" || ele == "RU"
                            || ele == "AS" || ele == " P" || ele == "AU" || ele == "HG"
                            || ele == "PD" || ele == "PB" || ele == "AG") {
                           handle_long_bonded_atom(atom_selection[i], uddHnd, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, col);
                        }
                     }
                  } else {
                     std::cout << "INFO:: trapped atom without residue in non-bonded atom check: "
                               << atom_selection[i] << std::endl;
                  }
               }
            }
         }

         // Make the stars...
         //
         for (int i=0; i<n_selected_atoms; i++) {
            if (atom_selection[i]->GetUDData(uddHnd, ic) == mmdb::UDDATA_Ok) {
               if (ic == graphical_bonds_container::NO_BOND) {
                  // no contact found
                  std::string res_name(atom_selection[i]->residue->GetResName());
                  if (res_name == "HOH")
                     if (! do_sticks_for_waters)
                        continue;
                  col = atom_colour(atom_selection[i], atom_colour_type, udd_user_defined_atom_colour_index_handle);
                  std::string ele = atom_selection[i]->element;
                  if (!is_hydrogen(ele) || draw_hydrogens_flag) {
                     coot::Cartesian atom(atom_selection[i]->x,
                                          atom_selection[i]->y,
                                          atom_selection[i]->z);

                     int iat_1 = 1; // 20171224-PE FIXME real
                     addBond(col, atom+small_vec_x, atom-small_vec_x, cc, imodel, iat_1, iat_1, true, true);
                     addBond(col, atom+small_vec_y, atom-small_vec_y, cc, imodel, iat_1, iat_1, true, true);
                     addBond(col, atom+small_vec_z, atom-small_vec_z, cc, imodel, iat_1, iat_1, true, true);
                  }
               }
            }
         }
         construct_from_model_links(asc.mol->GetModel(imodel), udd_atom_index_handle,
                                    udd_user_defined_atom_colour_index_handle, atom_colour_type);
      }
      asc.mol->DeleteSelection(SelectionHandle);
      do_disulphide_bonds(asc, imodel);
   }
   add_zero_occ_spots(asc);
   add_deuterium_spots(asc);
   int atom_colour_type = coot::COLOUR_BY_CHAIN;
   if (do_goodsell_colour_mode)
      atom_colour_type = coot::COLOUR_BY_CHAIN_GOODSELL;
   add_atom_centres(imol, asc, atom_colour_type, model_number, &atom_colour_map);
   add_cis_peptide_markup(asc, model_number);
}

void
Bond_lines_container::do_colour_by_chain_bonds_carbons_only(const atom_selection_container_t &asc,
                                                            int imol,
                                                            bool draw_missing_loops_flag,
                                                            int atom_colour_type,
                                                            int draw_hydrogens_flag) {

   auto _ = [] (int atom_colour_type) {
      std::string s = std::to_string(atom_colour_type);

      if (atom_colour_type == coot::COLOUR_BY_CHAIN)          s = "COLOUR_BY_CHAIN";
      if (atom_colour_type == coot::COLOUR_BY_CHAIN_C_ONLY)   s = "COLOUR_BY_CHAIN_C_ONLY";
      if (atom_colour_type == coot::COLOUR_BY_CHAIN_GOODSELL) s = "COLOUR_BY_CHAIN_GOODSELL";
      if (atom_colour_type == coot::COLOUR_BY_ATOM_TYPE)      s = "COLOUR_BY_ATOM_TYPE";
      if (atom_colour_type == coot::COLOUR_BY_SEC_STRUCT)     s = "COLOUR_BY_SEC_STRUCT";
      if (atom_colour_type == coot::DISULFIDE_COLOUR)         s = "DISULFIDE_COLOUR";
      if (atom_colour_type == coot::COLOUR_BY_MOLECULE)       s = "COLOUR_BY_MOLECULE";
      if (atom_colour_type == coot::COLOUR_BY_RAINBOW)        s = "COLOUR_BY_RAINBOW";
      if (atom_colour_type == coot::COLOUR_BY_OCCUPANCY)      s = "COLOUR_BY_OCCUPANCY";
      if (atom_colour_type == coot::COLOUR_BY_B_FACTOR)       s = "COLOUR_BY_B_FACTOR";
      if (atom_colour_type == coot::COLOUR_BY_USER_DEFINED_COLOURS)   s = "COLOUR_BY_USER_DEFINED_COLOURS";
      if (atom_colour_type == coot::COLOUR_BY_HYDROPHOBIC_SIDE_CHAIN) s = "COLOUR_BY_HYDROPHOBIC_SIDE_CHAIN";

      return s;
   };

   // std::cout << "in do_colour_by_chain_bonds_carbons_only() atom_colour_type " << atom_colour_type << std::endl;

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

   // std::cout << "debug:: colour by chain, carbons only" << std::endl;
   float max_dist = 1.9;
   float min_dist = 0.01; // As in the constructor
                          // Bond_lines_container::Bond_lines_container(const
                          // atom_selection_container_t &SelAtom, int
                          // do_disulphide_bonds_in, int
                          // do_bonds_to_hydrogens_in)
   mmdb::Contact *contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   // matrix stuff
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;
   int col_idx = 0; // atom (segment) colour
   coot::my_atom_colour_map_t atom_colour_map;
   atom_colour_map.fill_chain_id_map(asc);

   int udd_user_defined_atom_colour_index_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   std::vector<std::pair<bool, mmdb::Residue *> > het_residues; // bond these separately.

   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;

   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   int n_models = asc.mol->GetNumberOfModels();
   int udd_atom_index_handle = asc.UDDAtomIndexHandle;
   for (int imodel=1; imodel<=n_models; imodel++) {

      mmdb::PPAtom atom_selection = 0;
      int n_selected_atoms = 0;
      contact = NULL;

      // make a new atom selection, based on the model.
      int SelectionHandle = asc.mol->NewSelection(); // d
      asc.mol->SelectAtoms (SelectionHandle, imodel, "*",
                            mmdb::ANY_RES, // starting resno, an int
                            "*", // any insertion code
                            mmdb::ANY_RES, // ending resno
                            "*", // ending insertion code
                            "*", // any residue name
                            "*", // atom name
                            "*", // elements
                            "*"  // alt loc.
                            );

      asc.mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);

      asc.mol->SeekContacts(atom_selection, n_selected_atoms,
                            atom_selection, n_selected_atoms,
                            min_dist, max_dist, // min, max distances
                            0,        // seqDist 0 -> in same res also
                            contact, ncontacts,
                            0, &my_matt, i_contact_group);

      int uddHnd = asc.mol->RegisterUDInteger (mmdb::UDR_ATOM, "found bond");
      if (uddHnd<0)  {
         std::cout << " atom bonding registration failed.\n";
      } else {
         for (int i=0; i<n_selected_atoms; i++)
            atom_selection[i]->PutUDData(uddHnd, graphical_bonds_container::NO_BOND);
      }

      if (contact && ncontacts > 0) {

         mmdb::Atom *at1 = 0;
         mmdb::Atom *at2 = 0;
         std::string element1;
         std::string element2;
         int res1, res2;
         // int atom_colour_type = coot::COLOUR_BY_CHAIN_C_ONLY;

         for (int i=0; i< ncontacts; i++) {
            if (contact[i].id2 > contact[i].id1) {

               int iat_1 = contact[i].id1; // 20171224-PE FIXME - are these correct
               int iat_2 = contact[i].id2;

               at1 = atom_selection[ contact[i].id1 ];
               at2 = atom_selection[ contact[i].id2 ];

               res1 = at1->GetSeqNum();
               res2 = at2->GetSeqNum();

               if (abs(res1 - res2) < 2) {

                  std::string segid1(at1->GetChainID());
                  std::string segid2(at2->GetChainID());
                  int chain_idx = atom_colour_map.index_for_chain(segid1);

                  if (segid1 == segid2) {

                     element1 = at1->element;
                     element2 = at2->element;
                     if ( (draw_hydrogens_flag == 1) ||

                          // (element1 != " H" && element1 != " D" &&
                          //  element2 != " H" && element2 != " D") ) {

                          (! is_hydrogen(element1) && ! is_hydrogen(element2))) {

                        coot::Cartesian atom_1(at1->x, at1->y, at1->z);
                        coot::Cartesian atom_2(at2->x, at2->y, at2->z);

                        // alternate location test
                        //
                        std::string aloc_1(at1->altLoc);
                        std::string aloc_2(at2->altLoc);
                        //
                        if (aloc_1 == "" || aloc_2 == "" || aloc_1 == aloc_2) {

                           // std::cout << "in do_colour_by_chain_bonds_carbons_only() atom_colour_type " << atom_colour_type
                           // << " " << iat_1 << " " << iat_2 << std::endl;
                           do_colour_by_chain_bonds_carbons_only_internals(imol, imodel,
                                                                           chain_idx,
                                                                           at1, at2,
                                                                           iat_1, iat_2,
                                                                           &het_residues,
                                                                           element1,
                                                                           element2,
                                                                           atom_1,
                                                                           atom_2,
                                                                           atom_colour_type,
                                                                           uddHnd,
                                                                           udd_user_defined_atom_colour_index_handle);
                        }
                     } else {

                        // It was a hydrogen (or bonded to Hydrogen).
                        // Mark it as bonded (we don't want to see single
                        // unbonded (stared) hydrogen atoms.

                        // check the distance.
                        coot::Cartesian pt_1(at1->x, at1->y, at1->z);
                        coot::Cartesian pt_2(at2->x, at2->y, at2->z);

                        double d = (pt_1-pt_2).amplitude();
                        if (d < 1.5) {
                           if (uddHnd>=0) {
                              at1->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                              at2->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                           }
                        }
                     }
                  }
               }
            }
         }
         delete [] contact;
      }


      if (uddHnd>=0) {

         float star_size = 0.22;
         // for atoms with no neighbour (contacts):
         coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
         coot::Cartesian small_vec_y(0.0, star_size, 0.0);
         coot::Cartesian small_vec_z(0.0, 0.0, star_size);

         int ic; // changed by reference;
         int col;
         int atom_colour_type = coot::COLOUR_BY_CHAIN_C_ONLY; // the atoms connected to the SE (say) are Cs.
         for (int i=0; i<n_selected_atoms; i++) {
            if ( atom_selection[i]->GetUDData(uddHnd, ic) != mmdb::UDDATA_Ok) {
               std::cout << "ERROR:: do_colour_by_chain_bonds() failed to get ic bond info"
                         << std::endl;
            } else {
               // std::cout << "debug:: got ic " << ic << " for " << atom_selection[i] << std::endl;

               // 20190921-PE no longer consider P for additional bonds. Stops bond flashing in Goodsell
               //             mode.
               if ((ic == graphical_bonds_container::NO_BOND) ||
                   (!strcmp(atom_selection[i]->element, " S")) ||
                   (!strcmp(atom_selection[i]->element, "SE"))
                   // (!strcmp(atom_selection[i]->element, " P")
                   ) {

                  std::string segid(atom_selection[i]->GetChainID());
                  col = atom_colour_map.index_for_chain(segid);

                  //                   std::cout << " No contact for " << non_Hydrogen_atoms[i]
                  //                             << std::endl;

                  // no contact found or was Sulphur, or Phosphor

                  // So, was this a seleno-methione?
                  //
                  mmdb::Residue *atom_residue_p = atom_selection[i]->residue;
                  if (atom_residue_p) {
                     std::string resname = atom_selection[i]->GetResName();
                     if (resname == "MSE" || resname == "MET" || resname == "MSO" || resname == "CYS") {
                        handle_MET_or_MSE_case(atom_selection[i], uddHnd, udd_atom_index_handle,
                                               udd_user_defined_atom_colour_index_handle,
                                               atom_colour_type, &atom_colour_map);
                     } else {
                        std::string ele = atom_selection[i]->element;
                        if (ele == "CL" || ele == "BR" || ele == " S" ||  ele == " I"
                            || ele == "Cl" || ele == "Br"  || ele == "MO"
                            || ele == "PT" || ele == "RU"
                            || ele == "AS" || ele == " P" || ele == "AU" || ele == "HG"
                            || ele == "PD" || ele == "PB" || ele == "AG") {
                           handle_long_bonded_atom(atom_selection[i], uddHnd, udd_atom_index_handle,
                                                   udd_user_defined_atom_colour_index_handle,
                                                   atom_colour_type);
                        }
                     }
                  } else {
                     std::cout << "INFO:: trapped atom without residue in non-bonded atom check: "
                               << atom_selection[i] << std::endl;
                  }
               }
            }
         }

         // Make the stars...
         //
         for (int i=0; i<n_selected_atoms; i++) {
            if (atom_selection[i]->GetUDData(uddHnd, ic) == mmdb::UDDATA_Ok) {
               if (ic == graphical_bonds_container::NO_BOND) {
                  // no contact found
                  mmdb::Residue *residue_p = atom_selection[i]->residue;
                  std::string res_name(residue_p->GetResName());
                  if (res_name == "HOH")
                     if (! do_sticks_for_waters)
                        continue;
                  col = atom_colour(atom_selection[i], atom_colour_type, udd_user_defined_atom_colour_index_handle);
                  std::string ele = atom_selection[i]->element;
                  // if (ele != " H" || draw_hydrogens_flag) {
                  if (! is_hydrogen(ele) || draw_hydrogens_flag) {
                     coot::Cartesian atom(atom_selection[i]->x,
                                          atom_selection[i]->y,
                                          atom_selection[i]->z);

                     addBond(col, atom+small_vec_x, atom-small_vec_x, cc, imodel, i, i, true, true);
                     addBond(col, atom+small_vec_y, atom-small_vec_y, cc, imodel, i, i, true, true);
                     addBond(col, atom+small_vec_z, atom-small_vec_z, cc, imodel, i, i, true, true);
                  }
               }
            }
         }
         construct_from_model_links(asc.mol->GetModel(imodel), udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, atom_colour_type);
      }
      asc.mol->DeleteSelection(SelectionHandle);
      do_disulphide_bonds(asc, imodel);
      add_cis_peptide_markup(asc, imodel);
   }


   // for ligands in colour-by-chain mode to come out with carbons the
   // same colour as the main-chain (I think that) we need to pass a
   // my_atom_colour_map_t to add_bonds_het_residues().  I don't want
   // to do that (not at the moment, anyway). Have a look at atom_colour()
   // to see what I mean (c.f. COLOUR_BY_CHAIN and COLOUR_BY_CHAIN_C_ONLY).
   //
   // int atom_colour_type = coot::COLOUR_BY_CHAIN_C_ONLY; // passed now
   short int have_udd_atoms = false;
   int udd_handle = -1;

   add_bonds_het_residues(het_residues, asc, imol, atom_colour_type, have_udd_atoms, udd_handle, udd_atom_index_handle, udd_user_defined_atom_colour_index_handle);
   add_zero_occ_spots(asc);
   add_deuterium_spots(asc);
   // atom_colour_type = coot::COLOUR_BY_CHAIN;
   // atom_colour_type == coot::COLOUR_BY_CHAIN_GOODSELL;

   // std::cout << "DEBUG:: in " << __FUNCTION__ << " with atom_colour_type " << _(atom_colour_type) << std::endl;
   int model_number = 0; // all models for colour-by-chain
   add_atom_centres(imol, asc, atom_colour_type, model_number, &atom_colour_map);

   int udd_fixed_during_refinement_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "FixedDuringRefinement");
   if (draw_missing_loops_flag)
      atom_selection_missing_loops(asc, udd_atom_index_handle, udd_fixed_during_refinement_handle);


}

void
Bond_lines_container::do_colour_by_chain_bonds_carbons_only_internals(int imol, int imodel,
                                                                      int chain_idx,
                                                                      mmdb::Atom *at1, mmdb::Atom *at2,
                                                                      int iat_1, int iat_2,
                                                                      std::vector<std::pair<bool, mmdb::Residue *> > *het_residues_p,
                                                                      const std::string &element1,
                                                                      const std::string &element2,
                                                                      const coot::Cartesian &atom_1,
                                                                      const coot::Cartesian &atom_2,
                                                                      int atom_colour_type,
                                                                      int uddHnd,
                                                                      int udd_user_defined_atom_colour_index_handle) {

   if (atom_colour_type == coot::COLOUR_BY_CHAIN_GOODSELL) {
      do_colour_by_chain_bonds_internals_goodsell_mode(imol, imodel, chain_idx,
                                                       at1, at2, iat_1, iat_2,
                                                       het_residues_p,
                                                       element1, element2,
                                                       atom_1, atom_2,
                                                       uddHnd, udd_user_defined_atom_colour_index_handle);
      return;
   }

   bool bond_het_residue_by_dictionary =
      add_bond_by_dictionary_maybe(imol, at1, at2, het_residues_p); // add to het_residues maybe

   if (! bond_het_residue_by_dictionary) {

      if (element1 != element2) {

         // Bonded to different atom elements.
         //

         bool is_H = false;
         bool draw_it = true;
         if (element1 == " H") is_H = true;
         if (element2 == " H") is_H = true;
         if (is_H) {
            double d = (atom_1-atom_2).amplitude();
            if (d>1.5)
               draw_it = false;
         }

         if (draw_it) {

            coot::Cartesian bond_mid_point = atom_1.mid_point(atom_2);

            if (element1 != " C") {  // PDBv3 FIXME

               if (element2 != " C") {

                  // half bonds, e.g. N-O, N-H, O-H

                  // add here a test for either being H. In that caes
                  // we don't want half bonds.

                  if (is_H) {
                     graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                     addBond(HYDROGEN_GREY_BOND, atom_1, atom_2, cc, imodel, iat_1, iat_2);
                  } else {
                     graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                     int non_c_col = atom_colour(at1, atom_colour_type, udd_user_defined_atom_colour_index_handle);
                     bonds_size_colour_check(non_c_col);
                     addBond(non_c_col, atom_1, bond_mid_point, cc, imodel, iat_1, iat_2);
                     non_c_col = atom_colour(at2, atom_colour_type, udd_user_defined_atom_colour_index_handle);
                     bonds_size_colour_check(non_c_col);
                     addBond(non_c_col, atom_2, bond_mid_point, cc, imodel, iat_1, iat_2);
                  }

               } else {
                  // frequent
                  graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                  int non_c_col = atom_colour(at1, atom_colour_type, udd_user_defined_atom_colour_index_handle);
                  bonds_size_colour_check(non_c_col);
                  addBond(non_c_col, atom_1, bond_mid_point, cc, imodel, iat_1, iat_2);
                  bonds_size_colour_check(chain_idx);
                  addBond(chain_idx, atom_2, bond_mid_point, cc, imodel, iat_1, iat_2);
               }

            } else {

               // element 1 *is* a C

               if (element2 != " C") {

                  // frequent

                  if (is_H) {
                     graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                     addBond(HYDROGEN_GREY_BOND, atom_1, atom_2, cc, imodel, iat_1, iat_2);
                  } else {
                     graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                     bonds_size_colour_check(chain_idx);
                     addBond(chain_idx, atom_1, bond_mid_point, cc, imodel, iat_1, iat_2);
                     int non_c_col = atom_colour(at2, atom_colour_type, udd_user_defined_atom_colour_index_handle);
                     bonds_size_colour_check(non_c_col);
                     addBond(non_c_col, atom_2, bond_mid_point, cc, imodel, iat_1, iat_2);
                  }

               } else {
                  std::cout << "impossible " << std::endl;
                  bonds_size_colour_check(chain_idx);
                  graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
                  addBond(chain_idx, atom_2, bond_mid_point, cc, imodel, iat_1, iat_2);
               }
            }
         }

      } else {

         // same element

         if (element1 == " C") {
            bonds_size_colour_check(chain_idx);
            graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
            addBond(chain_idx, atom_1, atom_2, cc, imodel, iat_1, iat_2);
         } else {

            // If we are here: same element, not a carbon, and either drawing hydrogens
            // or these are not hydrogens, so don't draw bonds between hydrogens

            // if (element1 != " H") {
            if (! is_hydrogen(element1)) {
               int col = atom_colour(at1, atom_colour_type, udd_user_defined_atom_colour_index_handle);
               bonds_size_colour_check(col);
               graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
               addBond(col, atom_1, atom_2, cc, imodel, iat_1, iat_2);
            }
         }
      }
   }

   // we drew a bond.  Mark it up.
   if (uddHnd>=0) {
      at1->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
      at2->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
   }

}

// goodsell mode means all atoms in the chain are one colour (dependent on the chain).
// the colour for the carbons is slightly more pastel than the other atoms.
// So we need to colours per chain in the colour index.
//
void
Bond_lines_container::do_colour_by_chain_bonds_internals_goodsell_mode(int imol, int imodel,
                                                                       int chain_idx,
                                                                       mmdb::Atom *at1, mmdb::Atom *at2,
                                                                       int iat_1, int iat_2,
                                                                       std::vector<std::pair<bool, mmdb::Residue *> > *het_residues_p,
                                                                       const std::string &element1,
                                                                       const std::string &element2,
                                                                       const coot::Cartesian &atom_pos_1,
                                                                       const coot::Cartesian &atom_pos_2,
                                                                       int uddHnd,
                                                                       int udd_user_defined_atom_colour_index_handle) {

   // std::cout << "in do_colour_by_chain_bonds_internals_goodsell_mode " << at1 << " " << at2 << std::endl;

   bool bond_het_residue_by_dictionary =
      add_bond_by_dictionary_maybe(imol, at1, at2, het_residues_p); // add to het_residues maybe

   if (! bond_het_residue_by_dictionary) {
      bool draw_it = true;
      bool is_H = false;
      if (element1 == " H") is_H = true;
      if (element2 == " H") is_H = true;
      if (is_H) {
         double d = (atom_pos_1-atom_pos_2).amplitude();
         if (d>1.5)
            draw_it = false;
      }

      if (draw_it) {
         coot::Cartesian bond_mid_point = atom_pos_1.mid_point(atom_pos_2);
         if (element1 != " C") {  // PDBv3 FIXME
            if (element2 != " C") {
               graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
               int non_c_col = 2 * chain_idx + 1;
               bonds_size_colour_check(non_c_col);
               addBond(non_c_col, atom_pos_1, atom_pos_2, cc, imodel, iat_1, iat_2);
            } else {
               graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
               int non_c_col = 2 * chain_idx + 1;
               bonds_size_colour_check(non_c_col);
               addBond(non_c_col, atom_pos_1, bond_mid_point, cc, imodel, iat_1, iat_2);
               int c_col = 2 * chain_idx;
               bonds_size_colour_check(c_col);
               addBond(c_col, atom_pos_2, bond_mid_point, cc, imodel, iat_1, iat_2);
            }
         } else {
            // first atom *was* carbon
            if (element2 == " C" ) {
               graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
               int c_col = 2 * chain_idx;
               bonds_size_colour_check(c_col);
               addBond(c_col, atom_pos_1, atom_pos_2, cc, imodel, iat_1, iat_2);
            } else {
               graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
               int non_c_col = 2 * chain_idx + 1;
               bonds_size_colour_check(non_c_col);
               addBond(non_c_col, atom_pos_2, bond_mid_point, cc, imodel, iat_1, iat_2);
               int c_col = 2 * chain_idx;
               bonds_size_colour_check(c_col);
               addBond(c_col, atom_pos_1, bond_mid_point, cc, imodel, iat_1, iat_2);
            }
         }
         // we drew a bond.  Mark it up.
         if (uddHnd>=0) {
            at1->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
            at2->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
         }
      }
   }
}


void
Bond_lines_container::do_colour_by_ncs_related_chain_bonds(const atom_selection_container_t &asc,
                                                           int imol,
                                                           std::vector<std::vector<mmdb::Chain *> > ncs_related_chains,
                                                           int draw_mode,
                                                           bool change_c_only_flag, bool goodsell_mode) {

   if (true) { // test draw mode
      do_colour_by_ncs_related_chains_atoms_only(asc, imol, ncs_related_chains, change_c_only_flag, goodsell_mode);
   }

}


void
Bond_lines_container::do_colour_by_ncs_related_chains_atoms_only(const atom_selection_container_t &asc,
                                                                 int imol,
                                                                 std::vector<std::vector<mmdb::Chain *> > ncs_related_chains,
                                                                 bool change_c_only_flag, bool goodsell_mode) {

   // fill these
   atom_centres.clear();
   atom_centres_colour.clear(); // vector of ints

   std::map<mmdb::Chain *, int> chain_colour_indices;

   int icol = 0;
   for (const auto &vv : ncs_related_chains) {
      for (const auto &ch : vv) {
         chain_colour_indices[ch] = icol;
      }
      icol++;
   }
   int icol_max = icol;

   auto colour_for_chain = [&chain_colour_indices] (mmdb::Chain *chain_p) {
      int icol = 0;
      std::map<mmdb::Chain *, int>::const_iterator it = chain_colour_indices.find(chain_p);
      if (it != chain_colour_indices.end())
         icol = it->second;
      return icol;
   };

   int udd_user_defined_atom_colour_index_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = asc.mol->GetModel(imod);
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
                        int icol_base = colour_for_chain(chain_p);
                        icol = icol_base * 2 + 100; // this is the way for Goodsell colours
                        bool is_C = false; // only care about this if goodsell mode
                        if (goodsell_mode)
                           is_C = strncmp(at->element, " C", 2);
                        bool is_H_flag = (is_hydrogen(std::string(at->element)));
                        coot::Cartesian pos(at->x, at->y, at->z);
                        graphical_bonds_atom_info_t gbai(pos, iat, is_H_flag);
                        gbai.atom_p = at;
                        bool make_fat_atom = false; // because atoms are rendered as BALLS_NOT_BONDS, they don't
                                                    // need fattening here.
                        gbai.set_radius_scale_for_atom(at, make_fat_atom);
                        if (std::string(at->residue->GetResName()) == "HOH") gbai.is_water = true;
                        if (goodsell_mode) {
                           if (is_C) icol += 1; // pastel versions
                        }
                        bonds_size_colour_check(icol);

                        // does UDD colour trump the NCS chain colour?
                        int idx_col_udd;
                        if (at->GetUDData(udd_user_defined_atom_colour_index_handle, idx_col_udd) == mmdb::UDDATA_Ok)
                           icol = idx_col_udd;

                        atom_centres.push_back(gbai);
                        atom_centres_colour.push_back(icol);

                     }
                  }
               }
            }
         }
      }
   }

   // bonds.resize(icol_max); // needed?
}




void
Bond_lines_container::do_colour_by_molecule_bonds(const atom_selection_container_t &asc,
                                                  int imol,
                                                  int draw_hydrogens_flag) {

   // std::cout << "---------------------------------------- here we are in do_colour_by_molecule_bonds() " << std::endl;

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

   mmdb::Contact *contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   float max_dist = 1.9;
   float min_dist = 0.01; // As in the constructor
                          // Bond_lines_container::Bond_lines_container(const
                          // atom_selection_container_t &SelAtom, int
                          // do_disulphide_bonds_in, int
                          // do_bonds_to_hydrogens_in)

   // matrix stuff
   mmdb::mat44 my_matt;
   mmdb::SymOps symm;

   for (int i=0; i<4; i++)
      for (int j=0; j<4; j++)
         my_matt[i][j] = 0.0;

   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   int col = 0; // atom (segment) colour

   int n_models = asc.mol->GetNumberOfModels();
   int udd_atom_index_handle = asc.UDDAtomIndexHandle;

   int udd_user_defined_atom_colour_index_handle = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   for (int imodel=1; imodel<=n_models; imodel++) {

      mmdb::Model *model_p = asc.mol->GetModel(imodel);

      if (model_p) {

         mmdb::PPAtom atom_selection = 0;
         int n_selected_atoms = 0;
         contact = NULL;

         // make a new atom selection, based on the model.
         int SelectionHandle = asc.mol->NewSelection(); // d
         asc.mol->SelectAtoms (SelectionHandle, imodel, "*",
                               mmdb::ANY_RES, // starting resno, an int
                               "*", // any insertion code
                               mmdb::ANY_RES, // ending resno
                               "*", // ending insertion code
                               "*", // any residue name
                               "*", // atom name
                               "*", // elements
                               "*"  // alt loc.
                               );

         asc.mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);

         asc.mol->SeekContacts(atom_selection, n_selected_atoms,
                               atom_selection, n_selected_atoms,
                               min_dist, max_dist, // min, max distances
                               0,        // seqDist 0 -> in same res also
                               contact, ncontacts,
                               0, &my_matt, i_contact_group);

         coot::my_atom_colour_map_t atom_colour_map;
         atom_colour_map.fill_chain_id_map(asc);
         int res1, res2;
         // Now, let's not forget that some atoms don't have contacts, so
         // whenever we find a contact for an atom, we mark it with
         // UserDefinedData "found bond".
         //
         int uddHnd = asc.mol->RegisterUDInteger (mmdb::UDR_ATOM, "found bond");
         if (uddHnd<0)  {
            std::cout << " atom bonding registration failed.\n";
         } else {
            for (int i=0; i<n_selected_atoms; i++)
               atom_selection[i]->PutUDData(uddHnd, graphical_bonds_container::NO_BOND);
         }

         if (contact && ncontacts > 0) {

            std::string element1;
            std::string element2;
            for (int i=0; i< ncontacts; i++) {
               if (contact[i].id2 > contact[i].id1) {

                  int iat_1 = contact[i].id1; // 20171224-PE
                  int iat_2 = contact[i].id2;

                  mmdb::Atom *at1 = atom_selection[ contact[i].id1 ];
                  mmdb::Atom *at2 = atom_selection[ contact[i].id2 ];

                  res1 = at1->GetSeqNum();
                  res2 = at2->GetSeqNum();

                  if (abs(res1 - res2) < 2) {
                     coot::Cartesian atom_1(at1->x, at1->y, at1->z);
                     coot::Cartesian atom_2(at2->x, at2->y, at2->z);

                     element1 = at1->element;
                     element2 = at2->element;
                     if ( (draw_hydrogens_flag == 1) ||
                          // (element1 != " H" && element1 != " D" &&
                          // element2 != " H" && element2 != " D") ) {
                          (! is_hydrogen(element1) && ! is_hydrogen(element2))) {

                        // alternate location test
                        //
                        std::string aloc_1(at1->altLoc);
                        std::string aloc_2(at2->altLoc);
                        //
                        if (aloc_1 == "" || aloc_2 == "" || aloc_1 == aloc_2) {

                           // OK, draw a bond! (but not between 2 Hydrogens)

                           // if (! (element1 == " H" && element2 == " H")) {
                           if (! (is_hydrogen(element1) && is_hydrogen(element2))) {
                              bonds_size_colour_check(col);
                              addBond(col, atom_1, atom_2, cc, imodel, iat_1, iat_2);
                           }

                           if (uddHnd>=0) {
                              at1->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                              at2->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                           }
                        }
                     } else {
                        // It was a hydrogen (or bonded to Hydrogen).
                        // Mark it as bonded (we don't want to see single
                        // unbonded (stared) hydorgens).
                        if (uddHnd>=0) {
                           at1->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                           at2->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
                        }
                     }
                  }
               }
            }
            delete [] contact;
            contact = NULL;
         }

         if (uddHnd>=0) {

            float star_size = 0.3;
            // for atoms with no neighbour (contacts):
            coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
            coot::Cartesian small_vec_y(0.0, star_size, 0.0);
            coot::Cartesian small_vec_z(0.0, 0.0, star_size);

            int ic; // changed by reference;
            for (int i=0; i<n_selected_atoms; i++) {
               if ( atom_selection[i]->GetUDData(uddHnd,ic) == mmdb::UDDATA_Ok ) { // uddHnd for bond state
                  if (ic == 0) {

                     std::string res_name(atom_selection[i]->residue->GetResName());
                     if (res_name == "HOH")
                        if (! do_sticks_for_waters)
                           continue;

                     std::string segid(atom_selection[i]->GetChainID());
                     int col_inner = atom_colour_map.index_for_chain(segid);
                     bonds_size_colour_check(col_inner);
                     coot::Cartesian atom(atom_selection[i]->x,
                                          atom_selection[i]->y,
                                          atom_selection[i]->z);

                     addBond(col_inner, atom+small_vec_x, atom-small_vec_x, cc, imodel, i, i, true, true);
                     addBond(col_inner, atom+small_vec_y, atom-small_vec_y, cc, imodel, i, i, true, true);
                     addBond(col_inner, atom+small_vec_z, atom-small_vec_z, cc, imodel, i, i, true, true);
                  }
               }
            }
            construct_from_model_links(asc.mol->GetModel(imodel), udd_atom_index_handle, udd_user_defined_atom_colour_index_handle, coot::COLOUR_BY_CHAIN);
         }
         asc.mol->DeleteSelection(SelectionHandle);
         add_cis_peptide_markup(asc, imodel);
      }
   }
   add_zero_occ_spots(asc);
   add_deuterium_spots(asc);
   int atom_colour_type = coot::COLOUR_BY_MOLECULE;
   int model_number = 0; // all models for colour-by-molecule
   add_atom_centres(imol, asc, atom_colour_type, model_number, nullptr);
}


void
Bond_lines_container::add_zero_occ_spots(const atom_selection_container_t &SelAtom) {

   zero_occ_spots.clear();

   for (int i=0; i<SelAtom.n_selected_atoms; i++) {
      if (SelAtom.atom_selection[i]->occupancy < 0.01 &&
          SelAtom.atom_selection[i]->occupancy > -1) { // shelx occ test
         // we don't want to see atoms with occupancy -61 from a shelx ins
         // file with zero occupancy spots.
         std::string ele(SelAtom.atom_selection[i]->element);
         if (do_bonds_to_hydrogens ||
             ((do_bonds_to_hydrogens == 0) && (! is_hydrogen(ele)))) {
            if (no_bonds_to_these_atoms.find(i) == no_bonds_to_these_atoms.end())
               zero_occ_spots.push_back(coot::Cartesian(SelAtom.atom_selection[i]->x,
                                                        SelAtom.atom_selection[i]->y,
                                                        SelAtom.atom_selection[i]->z));
         }
      }
   }
}

void
Bond_lines_container::add_deuterium_spots(const atom_selection_container_t &SelAtom) {

   deuterium_spots.clear();

   for (int i=0; i<SelAtom.n_selected_atoms; i++) {
      std::string ele(SelAtom.atom_selection[i]->element);
      if (do_bonds_to_hydrogens && ele == " D")
         deuterium_spots.push_back(coot::Cartesian(SelAtom.atom_selection[i]->x,
                                                   SelAtom.atom_selection[i]->y,
                                                   SelAtom.atom_selection[i]->z));
   }
}


bool residue_sort_function(mmdb::Residue *r1, mmdb::Residue *r2) {

   return (coot::residue_spec_t(r1) < coot::residue_spec_t(r2));

}

void
Bond_lines_container::add_ramachandran_goodness_spots(const atom_selection_container_t &SelAtom) {


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

   ramachandran_goodness_spots.clear();
   std::set<mmdb::Residue *, bool(*)(mmdb::Residue *, mmdb::Residue *)> sorted_residues_set(residue_sort_function);

   for (int i=0; i<SelAtom.n_selected_atoms; i++) {
      mmdb::Residue *this_res = SelAtom.atom_selection[i]->residue;
      if (this_res) {
         sorted_residues_set.insert(this_res);
      }
   }

   std::set<mmdb::Residue *, bool(*)(mmdb::Residue *, mmdb::Residue *)>::const_iterator it;

   // we can't do sorted_residues_set[ii] or work out what is prev() or next() for it (experimental CXX?)
   // so convert to a vector
   std::vector<mmdb::Residue *> sorted_residues_vec(sorted_residues_set.size());
   unsigned int ii = 0;
   for (it=sorted_residues_set.begin(); it!=sorted_residues_set.end(); ++it) {
      sorted_residues_vec[ii] = *it;
      ii++;
   }

   if (sorted_residues_vec.size() > 2) {
      for (ii=1; ii<(sorted_residues_vec.size()-1); ii++) {
         mmdb::Residue *prev_res = sorted_residues_vec[ii-1];
         mmdb::Residue *this_res = sorted_residues_vec[ii];
         mmdb::Residue *next_res = sorted_residues_vec[ii+1];
         if (prev_res->GetChain() == this_res->GetChain()) {
            if (next_res->GetChain() == this_res->GetChain()) {
               if ((prev_res->GetSeqNum()+1) == this_res->GetSeqNum()) {
                  if (this_res->GetSeqNum() == (next_res->GetSeqNum()-1)) {

                     try {
                        coot::util::phi_psi_t pp(prev_res, this_res, next_res); // coot-rama.hh

                        if (false)
                           std::cout << "DEBUG:: pp " << coot::residue_spec_t(this_res)
                                     << " " << pp << std::endl;

                        mmdb::Atom *at = this_res->GetAtom(" CA "); // PDBv3 FIXME
                        if (at) {
                           coot::Cartesian pos(at->x, at->y, at->z);
                           coot::Cartesian offset_in_HA_dir_uv(0,0,1);
                           auto r = get_HA_unit_vector(this_res);
                           if (r.first)
                              offset_in_HA_dir_uv = r.second;
                           else
                              std::cout << "oooppps - missing HA vector" << std::endl;
                           std::pair<coot::Cartesian, coot::util::phi_psi_t> p(pos + offset_in_HA_dir_uv * 0.5, pp);
                           ramachandran_goodness_spots.push_back(p);
                        }
                     }
                     catch (const std::runtime_error &rte) {
                        if (false)
                           std::cout << "WARNING:: in failed to get phi,psi " << rte.what()
                                     << " for "
                                     << coot::residue_spec_t(prev_res) << " "
                                     << coot::residue_spec_t(this_res) << " "
                                     << coot::residue_spec_t(next_res) << std::endl;
                     }
                  }
               }
            }
         }
      }
   }

}



void
Bond_lines_container::add_rotamer_goodness_markup(const atom_selection_container_t &SelAtom) {

   dodecs = get_rotamer_dodecs(SelAtom);
}

// you can use the atom_colour_map here - it might be null
void
Bond_lines_container::add_atom_centres(int imol,
                                       const atom_selection_container_t &SelAtom,
                                       int atom_colour_type,
                                       int model_number,
                                       coot::my_atom_colour_map_t *atom_colour_map_p) {

   // 20230224-PE we want atoms in ligands with no dictionary to be bigger
   // So let's store a list of list of residue name and if they have a dictionary:
   std::map<std::string, bool> have_at_least_minimal_dictionary;

   atom_centres.clear();
   atom_centres_colour.clear(); // vector of ints

   int udd_user_defined_atom_colour_index_handle = SelAtom.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   // coot::my_atom_colour_map_t *atom_colour_map = 0;
   bool locally_created_atom_colour_map = false;
   if (atom_colour_map_p == 0) {
      if (atom_colour_type == coot::COLOUR_BY_CHAIN ||
          atom_colour_type == coot::COLOUR_BY_CHAIN_C_ONLY ||
          atom_colour_type == coot::COLOUR_BY_CHAIN_GOODSELL) {
         atom_colour_map_p = new coot::my_atom_colour_map_t;
         atom_colour_map_p->fill_chain_id_map(SelAtom);
         locally_created_atom_colour_map = true;
      }
   }

   for (int i=0; i<SelAtom.n_selected_atoms; i++) {

      mmdb::Atom *at = SelAtom.atom_selection[i];

      // Filter by model number if specified (model_number != 0)
      if (model_number != 0) {
         int atom_model = at->GetModelNum();
         if (atom_model != model_number)
            continue;
      }
      int idx = -1;
      at->GetUDData(SelAtom.UDDAtomIndexHandle, idx);
      // std::cout << "at " << at << "  idx: " << idx << std::endl;

      if (no_bonds_to_these_atoms.find(idx) == no_bonds_to_these_atoms.end()) {

         bool is_H_flag = false;
         std::string res_type(at->GetResName());
         bool have_dict_for_this_type = false;
         std::map<std::string, bool>::const_iterator it = have_at_least_minimal_dictionary.find(res_type);
         if (geom) {
            if (it == have_at_least_minimal_dictionary.end()) { // no hit in the cache
               bool s = geom->have_at_least_minimal_dictionary_for_residue_type(res_type, imol);
               have_at_least_minimal_dictionary[res_type] = s;
               have_dict_for_this_type = s;
            } else {
               have_dict_for_this_type = it->second;
            }
         }

         if (false)
            std::cout << "   geom: " << geom << " " << coot::atom_spec_t(at)
                      << " have_dict_for_this_type: " << have_dict_for_this_type << std::endl;

         if (is_hydrogen(std::string(at->element)))
            is_H_flag = true;

         if (do_bonds_to_hydrogens || (do_bonds_to_hydrogens == 0 && (!is_H_flag))) {
            coot::Cartesian pos(at->x, at->y, at->z);
            graphical_bonds_atom_info_t gbai(pos, idx, is_H_flag);

            // Fat atoms are for atom in residues with no dictionary - except
            // colour-by-molecule mode, where the dictionary check is not a thing.
            bool make_fat_atom = false;
            if (atom_colour_type != coot::COLOUR_BY_MOLECULE)
               if (! have_dict_for_this_type)
                  if (atom_colour_type != coot::COLOUR_BY_ATOM_TYPE)
                     make_fat_atom = true;
            // std::cout << " atom_colour_type " << atom_colour_type << " c.f. " << coot::COLOUR_BY_MOLECULE
            // << " make_fat_atom: " << make_fat_atom << std::endl;
            // 20240712-PE previous: gbai.set_radius_scale_for_atom(at, make_fat_atom);
            if (atom_colour_type != coot::COLOUR_BY_MOLECULE)
               if (! have_dict_for_this_type)
                  if (atom_colour_type != coot::COLOUR_BY_ATOM_TYPE)
                     gbai.set_radius_scale_for_atom_with_no_dictionary(at);

            gbai.set_radius_scale_for_atom(at, make_fat_atom);

            // this is a bit hacky
            if (atom_colour_type == coot::COLOUR_BY_USER_DEFINED_COLOURS)
               if (is_H_flag)
                  gbai.radius_scale += 0.18; // otherwise too tiny. At 0.25 Garib said that
                                             // the spheres were too big.

            // No small atoms (H) in COLOUR_BY_B_FACTOR or COLOUR_BY_OCCUPANCY
            // because the add_bond function doesn't take a "thin" flag
            // (thinning is only currently done by bond colour)
            //
	    mmdb::Residue *r = at->residue;
	    if (r) {
	       const char *rn = r->GetResName();
	       if (rn) {
		  std::string res_name = r->GetResName();
		  if (res_name == "HOH") gbai.is_water = true;
	       }
	    }
            if (is_H_flag) gbai.is_hydrogen_atom = true;
            gbai.atom_p = at;
            if (atom_colour_type == coot::COLOUR_BY_B_FACTOR)
               gbai.is_hydrogen_atom = false;
            if (atom_colour_type == coot::COLOUR_BY_OCCUPANCY)
               gbai.is_hydrogen_atom = false;
            if (false) // debugging large atom radius
               std::cout << "add_atom_centres() pushing back: " << coot::atom_spec_t(at)
                         << " with is_H_flag " << is_H_flag << " radius_scale " << gbai.radius_scale << std::endl;
            atom_centres.push_back(gbai);
            int icol = atom_colour(at, atom_colour_type, udd_user_defined_atom_colour_index_handle, atom_colour_map_p);
            bonds_size_colour_check(icol);
            atom_centres_colour.push_back(icol);
         }
      }
   }

   if (locally_created_atom_colour_map)
      delete atom_colour_map_p;
}

#include "coot-utils/coot-coord-extras.hh"

// if model_number is 0, do all models.
//
void
Bond_lines_container::add_cis_peptide_markup(const atom_selection_container_t &SelAtom, int model_number) {

   // cis_peptide_quads is the member data
   cis_peptide_quads.clear();

   std::vector<coot::util::cis_peptide_quad_info_t> quads =
      coot::cis_peptide_quads_from_coords(SelAtom.mol, model_number, geom); // geom can be null.

   for (unsigned int i=0; i<quads.size(); i++) {
      bool keep_this = true;
      int idx1 = quads[i].index_quad.index1;
      if (idx1 >= 0)
         if (no_bonds_to_these_atoms.find(idx1) != no_bonds_to_these_atoms.end())
            keep_this = false;
      int idx4 = quads[i].index_quad.index4;
      if (idx4 >= 0)
         if (no_bonds_to_these_atoms.find(idx4) != no_bonds_to_these_atoms.end())
            keep_this = false;
      if (keep_this)
         cis_peptide_quads.push_back(quads[i]);
   }

}

void
graphical_bonds_container::add_bad_CA_CA_dist_spots(const std::vector<coot::Cartesian> &bad_CA_CA_dist_spots_in) {

   unsigned int s = bad_CA_CA_dist_spots_in.size();
   if (s > 0) {
      n_bad_CA_CA_dist_spots = s;
      bad_CA_CA_dist_spots_ptr = new coot::Cartesian[s];
      for (std::size_t i=0; i<s; i++) {
        bad_CA_CA_dist_spots_ptr[i] = bad_CA_CA_dist_spots_in[i];
      }
   }
}

void
graphical_bonds_container::add_zero_occ_spots(const std::vector<coot::Cartesian> &spots) {

   n_zero_occ_spots = spots.size();

   if (n_zero_occ_spots > 0) {
      zero_occ_spots_ptr = new coot::Cartesian[n_zero_occ_spots];
      for (int j=0; j<n_zero_occ_spots; j++) {
         zero_occ_spots_ptr[j] = spots[j];
      }
   }
}


void
graphical_bonds_container::add_deuterium_spots(const std::vector<coot::Cartesian> &spots) {

   n_deuterium_spots = spots.size();

   if (n_deuterium_spots > 0) {
      deuterium_spots_ptr = new coot::Cartesian[n_deuterium_spots];
      for (int j=0; j<n_deuterium_spots; j++)
         deuterium_spots_ptr[j] = spots[j];
   }
}

void
graphical_bonds_container::add_ramachandran_goodness_spots(const std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t > > &spots,
                                                           const ramachandrans_container_t &rc) {

   n_ramachandran_goodness_spots = spots.size();

   if (n_ramachandran_goodness_spots > 0) {
      ramachandran_goodness_spots_ptr = new std::pair<coot::Cartesian, float>[n_ramachandran_goodness_spots];
      for (unsigned int i=0; i<spots.size(); i++) {

         const clipper::Ramachandran *rama = &rc.rama;

         if (spots[i].second.residue_name() == "PRO")
            rama = &rc.rama_pro;

         if (spots[i].second.residue_name() == "GLY")
            rama = &rc.rama_gly;

#ifdef CLIPPER_HAS_TOP8000
    if (spots[i].second.residue_name() == "ILE" ||
        spots[i].second.residue_name() == "VAL" )
       rama = &rc.rama_ileval;
    if (spots[i].second.is_pre_pro())
       if (spots[i].second.residue_name() != "GLY")
          rama = &rc.rama_pre_pro;
#endif

    // phi_psi_t needs to contain the next residue type to use
    // rama.pre_pro at some stage

         float rama_score = 10;

         if (rama->allowed(clipper::Util::d2rad(spots[i].second.phi()),
                           clipper::Util::d2rad(spots[i].second.psi())))
            rama_score = 3;
         if (rama->favored(clipper::Util::d2rad(spots[i].second.phi()),
                           clipper::Util::d2rad(spots[i].second.psi())))
            rama_score = 1;

         // ----- now lets do the size by probability
         rama_score = rama->probability(clipper::Util::d2rad(spots[i].second.phi()),
                                        clipper::Util::d2rad(spots[i].second.psi()));

         std::pair<coot::Cartesian, float> p(spots[i].first, rama_score);
         ramachandran_goodness_spots_ptr[i] = p;
      }
   }
}

void
graphical_bonds_container::add_rotamer_goodness_markup(const std::vector<rotamer_markup_container_t> &ric) {

   if (ric.size() > 0) {
      n_rotamer_markups = ric.size();
      rotamer_markups = new rotamer_markup_container_t[n_rotamer_markups];
      for (unsigned int i=0; i<ric.size(); i++)
         rotamer_markups[i] = ric[i];
   }
}

void
graphical_bonds_container::add_atom_centres(const std::vector<graphical_bonds_atom_info_t> &centres,
                                            const std::vector<int> &colours) {

   if (false) {
      std::cout << "+++++++++++++++++++++ In graphical_bonds_container::add_atom_centres() adding "
                << centres.size() << " atom centres groups" << std::endl;
      for (unsigned int i=0; i<centres.size(); i++)
         std::cout << "   Adding atom centre for atom " << centres[i].atom_p << std::endl;
   }

   if (colours.size() != centres.size()) {
      std::cout << "ERROR:: !! colours.size() != centres.size() in add_atom_centres\n";
   }
   n_atom_centres_ = centres.size();
   atom_centres_ = new graphical_bonds_atom_info_t[n_atom_centres_];
   atom_centres_colour_ = new int[n_atom_centres_];
   for (int i=0; i<n_atom_centres_; i++) {
      atom_centres_[i] = centres[i];
      atom_centres_colour_[i] = colours[i];
   }

   // now consolidate those to batches of each colour
   //
   // Comment from the ancient times:
   //   doing so is orders of magnitude faster without colour changing per atom.
   //
   int col_idx_max = 1;
   for (int i=0; i<n_atom_centres_; i++) {
      // std::cout << "add_atom_centres(): checking colours i " << i << " colours[i] " << colours[i] << std::endl;
      if (colours[i] > col_idx_max) {
         col_idx_max = colours[i];
      }
   }
   col_idx_max += 1;

   if (false)
      std::cout << "debug:: in graphical_bonds_atom_info_t::add_atom_centres() col_idx_max "
                << col_idx_max << std::endl;

   std::vector<unsigned int> counts(col_idx_max, 0);

   for (int i=0; i<n_atom_centres_; i++)
      counts[colours[i]]++;

   consolidated_atom_centres = new graphical_bonds_points_list<graphical_bonds_atom_info_t>[col_idx_max];
   n_consolidated_atom_centres = col_idx_max;

   if (false)
      std::cout << "debug:: in graphical_bonds_atom_info_t::add_atom_centres() n_consolidated_atom_centres "
                << n_consolidated_atom_centres << std::endl;

   for (int i=0; i<col_idx_max; i++)
      consolidated_atom_centres[i] = graphical_bonds_points_list<graphical_bonds_atom_info_t>(counts[i]);

   for (int i=0; i<n_atom_centres_; i++) {
      consolidated_atom_centres[colours[i]].add_point(atom_centres_[i]);
      int n_points = consolidated_atom_centres[colours[i]].num_points;
      if (false) {
         for (int j=0; j<n_points; j++) {
            const auto &point = consolidated_atom_centres[colours[i]].points[j];
            std::cout << "   " << i << " " << j
                      << " " << coot::atom_spec_t(point.atom_p)
                      << " radius_scale: " << point.radius_scale << std::endl;
         }
      }
   }

   if (false)  { // debug
      for (int i=0; i<n_atom_centres_; i++)
         std::cout << "---- graphical_bonds_container::add_atom_centres() " << i << " " << atom_centres_[i].position << "\n";
   }

   if (false) {
      for (int i=0; i<col_idx_max; i++)
         std::cout << "in add_atom_centres():  col " << i << " has " << consolidated_atom_centres[i].num_points << std::endl;
   }

}

void
graphical_bonds_container::add_cis_peptide_markup(const std::vector<coot::util::cis_peptide_quad_info_t> &cis_peptide_quads) {

   if (cis_peptide_quads.size()) {
      // all cis_peptide_quads are valid, right?
      n_cis_peptide_markups = cis_peptide_quads.size();
      cis_peptide_markups = new graphical_bonds_cis_peptide_markup[n_cis_peptide_markups];

      for (unsigned int i=0; i<cis_peptide_quads.size(); i++) {
         const coot::atom_quad &q = cis_peptide_quads[i].quad;
         coot::Cartesian c_1(q.atom_1->x, q.atom_1->y, q.atom_1->z);
         coot::Cartesian c_2(q.atom_2->x, q.atom_2->y, q.atom_2->z);
         coot::Cartesian c_3(q.atom_3->x, q.atom_3->y, q.atom_3->z);
         coot::Cartesian c_4(q.atom_4->x, q.atom_4->y, q.atom_4->z);
         bool pre_pro_flag = false;
         bool twisted_trans_flag = false;
         if (cis_peptide_quads[i].type == coot::util::cis_peptide_quad_info_t::PRE_PRO_CIS)
            pre_pro_flag = true;
         if (cis_peptide_quads[i].type == coot::util::cis_peptide_quad_info_t::TWISTED_TRANS)
            twisted_trans_flag = true;
         int model_number = q.atom_1->GetModelNum();
         graphical_bonds_cis_peptide_markup m(c_1, c_2, c_3, c_4, pre_pro_flag, twisted_trans_flag, model_number);
         m.add_atom_index_quad(cis_peptide_quads[i].index_quad);
         cis_peptide_markups[i] = m;
      }
   }
}

// for user defined colours:
//
// return a colour index, and -1 on failure
//
int
Bond_lines_container::get_user_defined_col_index(mmdb::Atom *at, int udd_handle) const {

   int r = -1;

   int ic = 0;
   int ierr = at->GetUDData(udd_handle, ic);
   if (ierr == mmdb::UDDATA_Ok)
      r = ic;

   return r;

}

// lightweight update the bonds without recalculating the bonding. make_graphical_bonds() will
// still be needed for now.
void
Bond_lines_container::update(mmdb::Atom **atom_selection, int n_atoms) {

   for (auto &coloured_bonds : bonds)
      coloured_bonds.update(atom_selection, n_atoms);

   for (auto &ac : atom_centres)
      ac.update(atom_selection, n_atoms);
}

