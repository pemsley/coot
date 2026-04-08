/* src/molecule-class-info.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 by The University of York
 * Author: Paul Emsley
 * Copyright 2007 by The University of York
 * Copyright 2009 by The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street, 02110-1307
 * USA.
 */


#include "compat/coot-sysdep.h"

#ifdef _MSC_VER
#include <windows.h>
#endif

#include <stdexcept>
#include <algorithm>
#include <string.h>  // strncpy


#include <mmdb2/mmdb_manager.h>

#include "coords/Cartesian.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coot-utils/coot-map-utils.hh"
#include "molecule-class-info.h"
#include "graphics-info.h"

#include "utils/logging.hh"
extern logging logger;


// This is called by make_bonds_type_checked(), which is called by
// update_molecule_after_additions().
//
void
molecule_class_info_t::update_ghosts() {

   // std::cout << "debug:: here in update_ghosts() with show_ghosts_flag " << show_ghosts_flag << std::endl;

   if (show_ghosts_flag) {

      if (ncs_ghosts.size() > 0) {
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {

            if (ncs_ghosts[ighost].display_it_flag)
               ncs_ghosts[ighost].update_bonds(atom_sel.mol);
         }
      }
   }
}

void
molecule_class_info_t::debug_ghosts() const {

   std::cout << "debug:: in debug_ghosts() there are " << ncs_ghosts.size()
             << " ncs ghosts" << std::endl;
   std::cout << "debug:: in debug_ghosts() show_ghosts_flag " << show_ghosts_flag
             << std::endl;
   std::cout << "debug:: in debug_ghosts() ncs_ghosts_have_rtops_flag " << ncs_ghosts_have_rtops_flag
             << std::endl;
   for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
      const auto &ghost = ncs_ghosts[ighost];
      std::cout << "debug:: in debug_ghosts() Ghost idx " << ighost
                << " name: " << ghost.name
                << " chain-id: " << ghost.chain_id
                << " residue-matches-size: " << ghost.residue_matches.size()
                << " is-empty: " << ghost.is_empty()
                << " display_it_flag: " << ghost.display_it_flag
                << " SelectionHandle: " << ghost.SelectionHandle << std::endl;
      std::cout << "debug:: in debug_ghosts() rtop:\n" << ghost.rtop.format()
                << std::endl;
   }
}


void
molecule_class_info_t::delete_ghost_selections() {

   // 20060923 Bill Scott was reporting crashes in DeleteSelection
   // here.  He was not using ghosts.  But had been manipulating his
   // molecule a lot.
   //
   // So, I suspect that the selection handle for the ghosts was going
   // out of date, which means that it's bad to delete it (of course).
   //
   // So, let's not delete selection if ghosts are not being used.  I
   // use the is_empty() test (which seems to return 0 even if ghosts
   // where not turned on! - oh dear?) and so also use the displayed?
   // flag.
   //
   // Which means of course that the SelectionHandles of the NCS
   // ghosts can be bogus because we are not updating them properly
   // (and not only in this function - *any* function that changes the
   // atom selection at all is suspect (moving atom coords/bfacs/occ
   // is fine of course)).
   //
   // fill_ghost_info has a ncs_ghosts.resize(0), which is a potential
   // memory leak, but that's not as bad as a crash (which would
   // happen if we tried to DeleteSelection on SelectionHandles of out
   // of date atom selections) is it?  For the safety checks to be
   // removed all atom manipulation functions must be checked for
   // proper operation with NCS ghosts atom selection.
   //
   // Question: why is is_empty() 0 for a not-turned-on ghost?
   // (e.g. RNASA).
   //
   // Hmmm... reflection: NCS code is complex and crash-prone.

   // 20250831-PE Here we are again.
   // NCS Chain master A
   // Delete residue B72
   // We come from:
   // "molecule-class-info-other.cc", line 1032, in delete_residue [0x7fbe27f2b320]
   //    1029:                         if (ins_code == inscodestr) {
   //    1030:                            make_backup();
   //    1031:                            atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
   //   >1032:                            delete_ghost_selections();
   //    1033:                            chain->DeleteResidue(iseqno, inscode);
   //    1034:                            was_deleted = 1;
   //    1035:                            res = NULL;
   // We crash on the line:
   //                atom_sel.mol->DeleteSelection(ncs_ghosts[ighost].SelectionHandle);

   if (! ncs_ghosts.empty()) {
      if (false)
         std::cout << "debug:: in delete_ghost_selections() there are " << ncs_ghosts.size()
                   << " ncs ghosts" << std::endl;
      for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
         if (false)
            std::cout << "debug:: Ghost " << ighost << " is-empty: "
                      << ncs_ghosts[ighost].is_empty()
                      << " display_it_flag: " << ncs_ghosts[ighost].display_it_flag << std::endl;
         if (! ncs_ghosts[ighost].is_empty()) {
            if (ncs_ghosts[ighost].display_it_flag) {

               logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__),
                          "Not actually deleting selection for ghosts");

               // 20250901-PE - don't do this here and now - the selection handles may be out of date.
               // It seems tricky and time-consuming to fix. Let's live with a possible memory leak.
               // atom_sel.mol->DeleteSelection(ncs_ghosts[ighost].SelectionHandle);
            }
         }
      }
   }
}

// This is a ghost_molecule_display_t member function
void
drawn_ghost_molecule_display_t::update_bonds(mmdb::Manager *mol) {

   if (true) {
      std::cout << "ghost_molecule_display_t::update_bonds() " << std::endl;
      std::cout << "ghost_molecule_display_t::update_bonds() rtop " << std::endl;
      std::cout << rtop.format() << std::endl;
   }


   atom_selection_container_t asc;
   asc.mol = mol;

   // We should update the atom selection here: Yes, this needs to
   // happen.  Otherwise: a modification is made and when we come to
   // Bond_lines_container constructor below, which is given an atom
   // selection, it may well point to atoms that have been removed.
   // And then distaster [Bush was re-elected today].
   //
   // (I guess that we don't need to regenerate the transformation
   // matrix)
   //
   // (Note: trash the current selection first)
   //

   asc.atom_selection = NULL;
   SelectionHandle = mol->NewSelection();
   // std::cout << "ghost:: update_bonds new SelectionHandle: " << SelectionHandle << std::endl;
   int imod = 1;
   //
   mol->SelectAtoms(SelectionHandle, imod,
                    chain_id.c_str(),
                    mmdb::ANY_RES, "*",
                    mmdb::ANY_RES, "*",
                    "*", "*", "*", "*");

   asc.mol->GetSelIndex(SelectionHandle, asc.atom_selection,
                        asc.n_selected_atoms);
   asc.SelectionHandle = SelectionHandle;

   //    std::cout << "update_bonds ghost selection selected "
   // << asc.n_selected_atoms
   // << " atoms from chain " << chain_id << "\n";

   float min_dist = 0.1;
   float max_dist = 1.85;

   //    std::cout << "ghost molecule bonds molecule has " << asc.n_selected_atoms
   //              << " selected atoms" << std::endl;

   if (! graphics_info_t::use_graphics_interface_flag)
      return;

   Bond_lines_container bonds(asc, min_dist, max_dist);
   bonds_box = bonds.make_graphical_bonds();

   // now run through bonds_box and change the coordinates therein by
   // rtop.  This is not apparently not a sensible place to put this
   // functionality because it should be part of the
   // Bond_lines_container, perhaps.  But Bond_lines_container does
   // not use clipper at all and I don't want to start now.

   // debug
//    int ilines = 0;
//    for (int i=0; i<bonds_box.num_colours; i++)
//       ilines += bonds_box.bonds_[i].num_lines;
//     std::cout << "updating ghost bonds...(which has " << bonds_box.num_colours
//               << " colours and " << ilines << " lines)\n";

   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

   for (int i=0; i<bonds_box.num_colours; i++) {
      for (int j=0; j< bonds_box.bonds_[i].num_lines; j++) {

         clipper::Coord_orth a(bonds_box.bonds_[i].pair_list[j].positions.getStart().get_x(),
                               bonds_box.bonds_[i].pair_list[j].positions.getStart().get_y(),
                               bonds_box.bonds_[i].pair_list[j].positions.getStart().get_z());

         clipper::Coord_orth b(bonds_box.bonds_[i].pair_list[j].positions.getFinish().get_x(),
                               bonds_box.bonds_[i].pair_list[j].positions.getFinish().get_y(),
                               bonds_box.bonds_[i].pair_list[j].positions.getFinish().get_z());

         clipper::Coord_orth at = a.transform(rtop);
         clipper::Coord_orth bt = b.transform(rtop);

         if (false)
            std::cout << "ghost-bond-a " << a.format() << " ghost-bond-b" << b.format() << " "
                      << "ghost-bond-at " << at.format() << " ghost-bond-bt" << bt.format() << " "
                      << std::endl;

         coot::CartesianPair p(coot::Cartesian(at.x(), at.y(), at.z()),
                               coot::Cartesian(bt.x(), bt.y(), bt.z()));
         bonds_box.bonds_[i].pair_list[j] = graphics_line_t(p, cc, false, false, -1, -1, -1);
      }
   }

   auto cartesian_to_clipper = +[] (const coot::Cartesian &p) {
      return clipper::Coord_orth(p.get_x(), p.get_y(), p.get_z());
   };
   auto clipper_to_cartesian = +[] (const clipper::Coord_orth &co) {
      return coot::Cartesian(co.x(), co.y(), co.z());
   };

   for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++) {
      for (unsigned int i=0; i<bonds_box.consolidated_atom_centres[icol].num_points; i++) {
         graphical_bonds_atom_info_t &ai = bonds_box.consolidated_atom_centres[icol].points[i];
         clipper::Coord_orth pos = cartesian_to_clipper(ai.position);
         clipper::Coord_orth trans_pos = pos.transform(rtop);
         ai.position = clipper_to_cartesian(trans_pos);
      }
   }

   int bbt = coot::NORMAL_BONDS;
   std::vector<glm::vec4> colour_table;
   for (unsigned int i=0; i<15; i++) { colour_table.push_back(glm::vec4(0.4, 0.8, 0.2, 1.0)); }
   graphics_info_t::attach_buffers();

   std::cout << "ghost code needs reworking: update_bonds() for ghosts " << std::endl;

   mesh.make_graphical_bonds(bonds_box, bbt, Mesh::representation_mode_t::BALL_AND_STICK,
                             -1, false, 0.1, 0.08, 1, 8, 2, colour_table, *graphics_info_t::Geom_p());
   if (false)
      std::cout << "########### ghost mesh v and ts: " << mesh.vertices.size()
                << " " << mesh.triangles.size() << " with representation_type ball-and-stick"
                << std::endl;

}

void
molecule_class_info_t::draw_ncs_ghosts(Shader *shader_for_meshes,
                                       stereo_eye_t eye,
                                       const glm::mat4 &mvp,
                                       const glm::mat4 &model_rotation_matrix,
                                       const std::map<unsigned int, lights_info_t> &lights,
                                       const glm::vec3 &eye_position,
                                       const glm::vec4 &background_colour) {

   if (show_ghosts_flag) {
      for (auto &ghost : ncs_ghosts) {
         ghost.draw(shader_for_meshes, eye, mvp, model_rotation_matrix, lights, eye_position, background_colour);
      }
   }

}

void
drawn_ghost_molecule_display_t::draw(Shader *shader_p,
                                     stereo_eye_t eye,
                                     const glm::mat4 &mvp,
                                     const glm::mat4 &view_rotation_matrix,
                                     const std::map<unsigned int, lights_info_t> &lights,
                                     const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                     const glm::vec4 &background_colour) {

   if (false)
      std::cout << "ncs_ghosts::draw() n-verts: " << mesh.vertices.size()
                << " n-tris: " << mesh.triangles.size() << std::endl;
   glm::vec3 rc = graphics_info_t::get_rotation_centre();
   mesh.draw(shader_p, eye, mvp, view_rotation_matrix, lights, eye_position, rc, 1.0, background_colour, false, true, false);
}


// public interface
int
molecule_class_info_t::update_ncs_ghosts() {

   update_ghosts();
   return ncs_ghosts.size();
}

// Called on read pdb:
//
// fill ncs ghosts:  detect exact ncs automatically.
int
molecule_class_info_t::fill_ghost_info(short int do_rtops_flag,
                                       float homology_lev) {

   if (false)
      std::cout << "DEBUG::   --------------- in fill_ghost_info ------- with homology_lev "
                << homology_lev << std::endl;

   std::vector<std::string> chain_ids;
   std::vector<std::vector<std::pair<std::string, int> > > residue_types;
   std::vector<int> chain_atom_selection_handles;
   std::vector<short int> first_chain_of_this_type;

   bool allow_offset_flag = 0;
   if (is_from_shelx_ins_flag)
      allow_offset_flag = 1;

   // start from a blank slate:
   ncs_ghosts.clear();

   // OK, so after we have run through this function, there will be
   // rtops if we asked for them.
   if (do_rtops_flag)
      ncs_ghosts_have_rtops_flag = 1;


   if (atom_sel.n_selected_atoms > 0) {

      int n_models = atom_sel.mol->GetNumberOfModels();
      if (n_models > 0) {
         int imod = 1; // otherwise madness

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();

         // std::cout << "DEUBG:: in fill_ghost_info() nchains is " << nchains << std::endl;

         if (nchains <= 0) {
            std::cout << "bad nchains in molecule " << nchains
                      << std::endl;
         } else {
            chain_ids.resize(nchains);
            residue_types.resize(nchains);
            chain_atom_selection_handles.resize(nchains);
            first_chain_of_this_type.resize(nchains, 1);
            for (int ichain=0; ichain<nchains; ichain++) {
               chain_p = model_p->GetChain(ichain);
               if (! chain_p->isSolventChain()) { 
                  chain_ids[ichain] = chain_p->GetChainID();
                  int iselhnd = atom_sel.mol->NewSelection();
                  mmdb::PAtom *atom_selection = NULL;
                  int nSelAtoms;
                  atom_sel.mol->SelectAtoms(iselhnd, imod,
                                            chain_p->GetChainID(),
                                            mmdb::ANY_RES, "*",
                                            mmdb::ANY_RES, "*",
                                            "*", "*", "*", "*");
                  atom_sel.mol->GetSelIndex(iselhnd, atom_selection, nSelAtoms);
                  chain_atom_selection_handles[ichain] = iselhnd;
//                   std::cout << "DEBUG:: fill_ghost_info chain_atom_selection_handles[" <<
//                      ichain << "] is " << chain_atom_selection_handles[ichain] <<
//                      " for chain id :" << chain_p->GetChainID() << ":" << std::endl;

                  // debugging the atom selection
                  if (false) {
                     mmdb::PPAtom selatoms_1 = NULL;
                     int n_sel_atoms_1;
                     atom_sel.mol->GetSelIndex(iselhnd, selatoms_1, n_sel_atoms_1);
                      std::cout << "DEBUG:: fill_ghost_info: first atom of " << n_sel_atoms_1
                                << " in " << chain_p->GetChainID()
                                << "  " << iselhnd
                                << "  selection " << selatoms_1[0] << std::endl;
                  }

                  int nres = chain_p->GetNumberOfResidues();
                  residue_types[ichain].resize(nres);
//                   std::cout << "INFO:: residues_types[" << ichain << "] resized to "
//                             << residue_types[ichain].size() << std::endl;
                  mmdb::PResidue residue_p;
                  for (int ires=0; ires<nres; ires++) {
                     residue_p = chain_p->GetResidue(ires);
                     std::string resname(residue_p->name);
                     residue_types[ichain][ires] =
                        std::pair<std::string, int> (resname, residue_p->seqNum);
                  }
                  // atom_sel.mol->DeleteSelection(iselhnd); not yet!
               }
            }
         }
      }

      if (false) {
         std::cout << "DEBUG:: fill_ghost_info allow_offset_flag: " << allow_offset_flag << std::endl;
         std::cout << "DEBUG:: calling add_ncs_ghosts_no_explicit_master() with first_chain_of_this_type ";
         for (unsigned int ifc=0; ifc<first_chain_of_this_type.size(); ifc++) {
            std::cout << "   " << ifc << ": " << first_chain_of_this_type[ifc] << " ";
         }
         std::cout << std::endl;
         std::cout << "DEBUG:: calling add_ncs_ghosts_no_explicit_master() with chain_ids: ";
         for (unsigned int ich=0; ich<chain_ids.size(); ich++) {
            std::cout << chain_ids[ich] << " ";
         }
         std::cout << std::endl;
      }

      add_ncs_ghosts_no_explicit_master(chain_ids, residue_types, first_chain_of_this_type,
                                        chain_atom_selection_handles, do_rtops_flag, homology_lev,
                                        allow_offset_flag);

      if (! ncs_ghosts.empty()) {
         update_ghosts();
         // std::cout << "  INFO:: fill_ghost_info Constructed " << ncs_ghosts.size() << " ghosts\n";
         int n_ghosts = ncs_ghosts.size();
         logger.log(log_t::INFO, std::string("Constructed"), std::to_string(n_ghosts),
                    std::string("ghosts"));
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            // std::cout << "      Ghost " << ighost << " name: \"" << ncs_ghosts[ighost].name
            // << "\"" << std::endl;
            std::string name = "\"" + ncs_ghosts[ighost].name + "\"";
            logger.log(log_t::INFO, "     Ghost index:", std::to_string(ighost), "name", name);
         }
      }
   }
   return ncs_ghosts.size();
}

void
molecule_class_info_t::add_ncs_ghosts_no_explicit_master(const std::vector<std::string> &chain_ids,
                                                         const std::vector<std::vector<std::pair<std::string, int> > > &residue_types,
                                                         std::vector<short int> first_chain_of_this_type, // modified
                                                         const std::vector<int> &chain_atom_selection_handles,
                                                         short int do_rtops_flag,
                                                         float homology_level,
                                                         bool allow_offset_flag) {

   // debug input
   if (false) {
      std::cout << " in add_ncs_ghosts_no_explicit_master() ncs chain_ids are ";
      for (unsigned int i=0; i<chain_ids.size(); i++) {
         std::cout << ":" << chain_ids[i] << ": ";
      }
      std::cout << std::endl;
      std::cout << " and there are " << residue_types.size() << " elements in the residue_types vector "
                << std::endl;
   }


   // when checking for ncs, we must make sure that the number of
   // residues is greater than 0 (because we get 0 for solvent
   // chains).
   // Actually, greater than 2 would be even more sensible.
   //

   // so I have marked some code as "trickiness" -
   // Consider the case of alpha-3 beta-2, A,B,C of type alpha D,E of type beta
   //
   // We want B onto A and C onto A and E onto D - but not C onto B.
   // We when we find a match, we have to mark the matcher "second" (i.e.
   // not the reference) that this is not a unique - the matcher target
   // for this type of molecule is higher up the list, we don't want to see
   // matchers to this "second" chain.
   //

   // std::cout << "DEBUG:: Checking chains for NCS.. (no explicit master)" << std::endl;
   // So now let's check the chains against each other:
   for (unsigned int ifirst=0; ifirst<(chain_ids.size()-1); ifirst++) {
      try {
         if (first_chain_of_this_type[ifirst]) { // trickiness
            for (unsigned int isec=(ifirst+1); isec<chain_ids.size(); isec++) {
               if (false)
                  std::cout << "DEBUG:: checking chains numbers "
                            << ifirst << " (" << chain_ids[ifirst] << ") and "
                            << isec << " (" << chain_ids[isec] << ") with homology_level "
                            << homology_level << std::endl;
               if (ncs_chains_match_p(residue_types[ifirst],
                                      residue_types[isec],
                                      homology_level,
                                      allow_offset_flag)) {

                  if (false)
                     std::cout << "DEBUG:: ncs_chains match! =================  "
                               << ifirst << " (" << chain_ids[ifirst] << ") and "
                               << isec << " (" << chain_ids[isec] << ")" << std::endl;

                  first_chain_of_this_type[isec] = 0; // trickiness
                  drawn_ghost_molecule_display_t ghost;
                  // slow...
                  if (do_rtops_flag) {
                     if (false)
                        std::cout << "DEBUG:: add_ncs_ghosts_no_explicit_master: first: "
                                  << ifirst << " " << chain_atom_selection_handles[ifirst]
                                  << " and isec: " << isec << " " << chain_atom_selection_handles[isec]
                                  << std::endl;

                     coot::ncs_matrix_info_t ghost_info =
                        find_ncs_matrix(chain_atom_selection_handles[ifirst],
                                        chain_atom_selection_handles[isec]);
                     if (ghost_info.state) {
                        ghost.rtop = ghost_info.rtop;
                        ghost.display_it_flag = 1;
                        ghost.residue_matches = ghost_info.residue_matches;
                     }
                  }

                  ghost.SelectionHandle = chain_atom_selection_handles[isec];
                  ghost.target_chain_id = chain_ids[ifirst];
                  ghost.chain_id = chain_ids[isec];
                  ghost.name = "NCS found from ";
                  ghost.name += "matching Chain ";
                  ghost.name += chain_ids[isec];
                  ghost.name += " onto Chain ";
                  ghost.name += chain_ids[ifirst];
                  // ghost.bonds_box filled by update_ghosts().
                  ncs_ghosts.push_back(ghost);
               } else {
                  if (false)
                     std::cout << "DEBUG:: ncs_chains NO MATCH =================  "
                               << ifirst << " (" << chain_ids[ifirst] << ") and "
                               << isec << " (" << chain_ids[isec] << ")" << std::endl;
               }
            }
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << rte.what() << std::endl;
      }
   }

};

void
molecule_class_info_t::add_ncs_ghosts_using_ncs_master(const std::string &master_chain_id,
                                                       const std::vector<std::string> &chain_ids,
                                                       const std::vector<std::vector<std::pair<std::string, int> > > &residue_types,
                                                       const std::vector<int> &chain_atom_selection_handles,
                                                       float homology_level) {

   // float homology_level = 0.7;
   bool allow_offset_flag = 0;
   // First find imaster
   int imaster = -1;
   for (unsigned int ichain=0; ichain<chain_ids.size(); ichain++) {
      if (chain_ids[ichain] == master_chain_id) {
         imaster = ichain;
         break;
      }
   }

   std::cout << "DEBUG:: add_ncs_ghosts_using_ncs_master():   %%%%%% imaster: " << imaster << std::endl;

   if (imaster != -1) {
      std::cout << "   Checking chains for NCS matching to chain " << master_chain_id << std::endl;
      for (unsigned int ichain=0; ichain<chain_ids.size(); ichain++) {
         try {
            if (chain_ids[ichain] != chain_ids[imaster]) {
               if (ncs_chains_match_p(residue_types[imaster],
                                      residue_types[ichain],
                                      homology_level,
                                      allow_offset_flag)) {
                  coot::ghost_molecule_display_t ghost;
                  coot::ncs_matrix_info_t ghost_info = 
                     find_ncs_matrix(chain_atom_selection_handles[imaster],
                                     chain_atom_selection_handles[ichain]);
                  if (ghost_info.state) { 
                     ghost.rtop = ghost_info.rtop;
                     ghost.SelectionHandle = chain_atom_selection_handles[ichain];
                     ghost.target_chain_id = master_chain_id;
                     ghost.chain_id = chain_ids[ichain];
                     ghost.name = "NCS found from matching Chain ";
                     ghost.name += chain_ids[ichain];
                     ghost.name += " onto Chain ";
                     ghost.name += master_chain_id;
                     ghost.display_it_flag = 1;
                     std::cout << "   Adding ghost with name: " << ghost.name << std::endl;
                     ncs_ghosts.push_back(ghost);
                     ncs_ghosts_have_rtops_flag = 1;
                  }
               }
            }
         }
         catch (const std::runtime_error &rte) {
            std::cout << rte.what() << std::endl;
         } 
      }
   }
}




void
molecule_class_info_t::set_show_ghosts(short int state) {

   show_ghosts_flag = state; // bool
   // caller redraws
}


// Throw an exception when the matrix is not good.
//
// was std::pair<bool, clipper::RTop_orth>
coot::ncs_matrix_info_t
molecule_class_info_t::find_ncs_matrix(int SelHandle1, int SelHandle2) const {

   bool rtop_is_good = 0; 
   clipper::RTop_orth rtop; // random numbers
   std::vector<int> residue_matches; 
   if (graphics_info_t::ncs_matrix_flag == coot::NCS_SSM) {
#ifdef HAVE_SSMLIB

      ssm::PRECISION precision = ssm::PREC_Normal;
      ssm::CONNECTIVITY connectivity = ssm::CONNECT_Flexible;
      ssm::Align *SSMAlign = new ssm::Align();
      int rc = SSMAlign->AlignSelectedMatch(atom_sel.mol, atom_sel.mol,
                                            precision, connectivity, SelHandle2, SelHandle1);

      if (rc)  {
         std::string ws;
         switch (rc)  {

         case ssm::RC_NoHits :
            std::cout << " *** secondary structure does not match.\n";
            ws = "secondary structure does not match";
            break;
         case ssm::RC_NoSuperposition :
            std::cout << " *** structures are too remote.\n";
            ws = "structures are too remote";
            break;
         case ssm::RC_NoGraph :
            std::cout << " *** can't make graph for 1 \n";
            ws = "can't make graph for 1 ";
            ws += " structure";
            break;
         case ssm::RC_NoVertices :
            std::cout << " *** empty graph for 1\n";
            ws = "empty graph for 1";
            break;
         case ssm::RC_NoGraph2 :
            std::cout << " *** can't make graph for 2\n";
            ws = "can't make graph for 2 ";
            break;
         case ssm::RC_NoVertices2 :
            std::cout << " *** empty graph for 2\n";
            ws = "empty graph for 2";
            break;
         default :
            std::cout << " *** undocumented return code: " << rc << "\n";
         }
      } else  {

         mmdb::PPAtom selatoms_1 = NULL;
         int n_sel_atoms_1; 
         atom_sel.mol->GetSelIndex(SelHandle1, selatoms_1, n_sel_atoms_1);
         mmdb::PPAtom selatoms_2 = NULL;
         int n_sel_atoms_2; 
         atom_sel.mol->GetSelIndex(SelHandle2, selatoms_2, n_sel_atoms_2);
         // debugging the atom selection
         if (false) {
            std::cout << "First atom of " << n_sel_atoms_1 << " in first  selection "
                      << selatoms_1[0] << std::endl;
            std::cout << "First atom of " << n_sel_atoms_2 << " in second selection "
                      << selatoms_2[0] << std::endl;
         }
         rtop = coot::util::make_rtop_orth_from(SSMAlign->TMatrix);
         rtop_is_good = 1;

         // add residue_matches
         // 
         // int      selHndCa1,selHndCa2; // selection handles to used C-alphas
         // mmdb::ivector  Ca1,Ca2;      // C-alpha correspondence vectors
         // Ca1[i] corresponds to a[i], where a is
         // selection identified by selHndCa1
         mmdb::PAtom *atom_selection1 = NULL;
         int n_selected_atoms_1;
         atom_sel.mol->GetSelIndex(SSMAlign->selHndCa1,
                                   atom_selection1,
                                   n_selected_atoms_1);
         for (int iat=0; iat<n_selected_atoms_1; iat++) {
            residue_matches.push_back(SSMAlign->Ca1[iat]);
         }

         if (false) {
            // Too much noise.
            std::cout << "  Residue Matches" << std::endl;
            for (unsigned int imatch=0; imatch<residue_matches.size(); imatch++) {
               std::cout << "match " << imatch << " " << residue_matches[imatch]
                         << std::endl;
            }
         }
      }

      std::cout << "   find_ncs_matrix returns (SSM) ";
      if (! rtop_is_good) {
         std::cout << "junk";
         delete SSMAlign;
         SSMAlign = 0;
         std::string mess = "No NCS matrix defined";
         throw std::runtime_error(mess);
      }
      std::cout << "\n" << rtop.format() << std::endl;
      delete SSMAlign; // added 071101
      SSMAlign = 0;

#endif // HAVE_SSMLIB
   } else {
      // debugging the atom selection
      mmdb::PPAtom selatoms_1 = NULL;
      int n_sel_atoms_1; 
      atom_sel.mol->GetSelIndex(SelHandle1, selatoms_1, n_sel_atoms_1);
      mmdb::PPAtom selatoms_2 = NULL;
      int n_sel_atoms_2; 
      atom_sel.mol->GetSelIndex(SelHandle2, selatoms_2, n_sel_atoms_2);
      if (1) { 
         std::cout << "First atom of " << n_sel_atoms_1 << " in first  selection "
                   << selatoms_1[0] << std::endl;
         std::cout << "First atom of " << n_sel_atoms_2 << " in second selection "
                   << selatoms_2[0] << std::endl;
      }

      int match_type = 2; // CA
      char *chain_id_reference = selatoms_1[0]->GetChainID();
      char *chain_id_moving    = selatoms_2[0]->GetChainID();
      // assuming no solvents at the start
      int reference_resno_start = selatoms_1[0]->GetResidue()->GetSeqNum();
      int moving_resno_start    = selatoms_2[0]->GetResidue()->GetSeqNum();
      int reference_resno_end   = selatoms_1[n_sel_atoms_1 - 1]->GetResidue()->GetSeqNum();
      int moving_resno_end      = selatoms_2[n_sel_atoms_2 - 1]->GetResidue()->GetSeqNum();

      // find the last non-solvent residue
      for (int i=n_sel_atoms_1-1; i>0; i--) {
         if (selatoms_1[i]->GetResidue()->isAminoacid() ||
             selatoms_1[i]->GetResidue()->isNucleotide()) {
            reference_resno_end   = selatoms_1[i]->GetResidue()->GetSeqNum();
            break;
         }
      }
      for (int i=n_sel_atoms_2-1; i>0; i--) {
         if (selatoms_2[i]->GetResidue()->isAminoacid() ||
             selatoms_2[i]->GetResidue()->isNucleotide()) {
            moving_resno_end   = selatoms_2[i]->GetResidue()->GetSeqNum();
            break;
         }
      }
      // now get the lowest common denominator
      int resno_start = std::max(reference_resno_start, moving_resno_start);
      int resno_end   = std::min(reference_resno_end, moving_resno_end);
      coot::lsq_range_match_info_t matches(resno_start, resno_end,
                                         chain_id_reference,
                                         resno_start, resno_end,
                                         chain_id_moving, match_type);
//      g_print("BL DEBUG:: params for range %i %i %s %i %i %s %i\n",
//              resno_start, resno_end,
//              chain_id_reference,
//              resno_start, resno_end,
//              chain_id_moving, match_type);
      std::vector<coot::lsq_range_match_info_t> ncs_vector;
      ncs_vector.push_back(matches);
      std::pair<short int, clipper::RTop_orth> rtop_info =
         coot::util::get_lsq_matrix(atom_sel.mol,
                                    atom_sel.mol,
                                    ncs_vector,
                                    graphics_info_t::ncs_matrix_flag);

      rtop = rtop_info.second;
      rtop_is_good = rtop_info.first;

      std::cout << "   find_ncs_matrix returns (LSQ) ";
      std::cout << "\n" << rtop.format() << std::endl;
      if (! rtop_is_good) { 
         std::cout << "WARNING:: Junk NCS matrix" << std::endl;
         std::string mess = "No NCS matrix defined";
         throw std::runtime_error(mess);
      }
   }
   return coot::ncs_matrix_info_t(rtop_is_good, rtop, residue_matches);
}




void
molecule_class_info_t::set_ghost_bond_thickness(float f) {

   ghost_bond_width = f;
}


int
molecule_class_info_t::test_function() {
#if 0
   int imol;
   graphics_info_t g;

   if (ncs_ghosts.size() > 0) {
      if (ncs_ghosts_have_rtops_flag == 0) {
         float homology_lev =0.7;
         fill_ghost_info(1, homology_lev); // fill the rtops and set the flag
      }
   }
   std::cout <<  "make_dynamically_transformed_maps on " << ncs_ghosts.size()
             << " maps\n";

   std::vector<coot::ghost_molecule_display_t> local_ncs_ghosts = ncs_ghosts;
   int imol_base = graphics_info_t::n_molecules();
   for(unsigned int ighost=0; ighost<10; ighost++) {
      std::cout << "DEBUG:: pre-create molecule " << ighost << "/"
                << local_ncs_ghosts.size() << std::endl;
      std::cout << "DEBUG:: This is imol=" << imol_no << std::endl;
      imol = graphics_info_t::create_molecule();
   }

   imol = imol_base;
   std::cout << "DEBUG:: pre-second-loop: This is imol=" << imol_no << std::endl;
   for(unsigned int ighost=0; ighost<local_ncs_ghosts.size(); ighost++) {

      std::cout << "DEBUG:: This is imol=" << imol_no << std::endl;
      for (int itmp=0; itmp<=imol; itmp++)
         std::cout << "DEBUG:: molecule names: " << itmp << " :"
                   << graphics_info_t::molecules[itmp].name_ << ":" << std::endl;

      std::cout << "DEBUG:: NCS Copy to map number " << imol << std::endl;
      std::cout << "DEBUG:: pre-install of ghost map " << ighost << "/"
                << local_ncs_ghosts.size() << std::endl;
      std::cout << "DEBUG:: Post install of ghost map " << ighost << "/"
        << local_ncs_ghosts.size() << std::endl;
   }

   return imol;
#endif
   return 0;
}


std::vector<drawn_ghost_molecule_display_t>
molecule_class_info_t::NCS_ghosts() const {

   return ncs_ghosts;

}


// Not const because we may modify ncs_ghosts by adding their ncs operators.
// 
// This is a coordinates molecule:
// 
std::vector<std::pair<clipper::Xmap<float>, std::string> >
molecule_class_info_t::ncs_averaged_maps(const clipper::Xmap<float> &xmap_in,
                                         float homology_lev, std::string &imol_map_name) {

   std::vector<std::pair<clipper::Xmap<float>, std::string> > annotated_xmaps;

   // First, let's make the ncs operators if they need to be made:
   if (ncs_ghosts.size() > 0) {
      if (ncs_ghosts_have_rtops_flag == 0) {
//          std::cout << "   %%%%%%%%% calling fill_ghost_info from ncs_averaged_maps "
//                    << std::endl;
         fill_ghost_info(1, homology_lev); // fill the rtops and set the flag
      }

      // Sort the targets into those that match to molecule A,
      // those that match to molecule B ... etc
      //
      // index      0         1
      // match to   A         D
      // chainid    B C       E F
      //
      std::vector<std::string> reference_ids;

      for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
         short int ifound = 0;
         for (unsigned int ir=0; ir<reference_ids.size(); ir++) {
            if (ncs_ghosts[ighost].target_chain_id == reference_ids[ir]) {
               ifound = 1;
               break;
            }
         }
         if (ifound == 0) {
            reference_ids.push_back(ncs_ghosts[ighost].target_chain_id);
         }
      }

      // so now reference_ids contains a vector of target ids
      // e.g: A C
      int n_references = reference_ids.size();
      // So reference_matchers will contain vectors of ghost indexes
      // of molecules that match the reference id (but not the
      // reference molecule).
      //
      std::vector<std::vector<int> > reference_matchers(n_references);

      for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
         for (unsigned int ir=0; ir<reference_ids.size(); ir++) {
            if (ncs_ghosts[ighost].target_chain_id == reference_ids[ir]) {
               reference_matchers[ir].push_back(ighost);
            }
         }
      }

//       // debugging:
//       for (int irm=0; irm<reference_matchers.size(); irm++) {
//          std::cout << "DEBUG:: reference_matchers[" << irm << "] has size "
//                 << reference_matchers[irm].size() << std::endl;
//          for (int ir=0; ir<reference_matchers[irm].size(); ir++) {
//             std::cout << "DEBUG:: reference_matchers[" << irm << "]["
//                       << ir << "] element: "
//                       << reference_matchers[irm][ir] << std::endl;
//          }
//       }

      // So now reference_matchers contains a vector of ghost indexes:
      // {1} {3 4}
      // ie B to A, D to C, E to C
      //
      // Let's find the extents of reference maps: (we need to select
      // all the atoms in the reference chain).
      //
      // 
      std::vector<int> SelectionHandle(reference_ids.size(), -1);
      for (unsigned int iref=0; iref<reference_ids.size(); iref++) {

         SelectionHandle[iref] = atom_sel.mol->NewSelection();
         atom_sel.mol->SelectAtoms(SelectionHandle[iref], 0,
                                   (char *) reference_ids[iref].c_str(),
                                   mmdb::ANY_RES, "*", mmdb::ANY_RES, "*",
                                   "*", "*", "*", "*");
         std::pair<clipper::Coord_orth, clipper::Coord_orth> extents_pair =
            coot::util::extents(atom_sel.mol, SelectionHandle[iref]);
         atom_sel.mol->DeleteSelection(SelectionHandle[iref]);
         // Now lets just extend that out a bit (5A)
         clipper::Coord_orth bit(5,5,5);
         // clipper::Coord_orth bit(0,0,0); // testing
         extents_pair.first  -= bit;
         extents_pair.second += bit;
         

         // Each reference_matchers should give us a NCS new map.
         //
         // These maps are contributed to by the initial map and the
         // ncs_ghosts of the given indices.
         
         clipper::Xmap<float> start_map = xmap_in;
         clipper::Xmap<float> running   = xmap_in;

         // this was useful for debugging (as long as I remembered that
         // it was happening)...
         // clipper::Xmap_base::Map_reference_index ix;
         // for (ix = running.first(); !ix.last(); ix.next() )
         // running[ix] = 0.0;

         for (unsigned int ir=0; ir<reference_matchers[iref].size(); ir++) {
            clipper::RTop_orth rtop = ncs_ghosts[reference_matchers[iref][ir]].rtop;
            std::cout << "Reference matcher: " << iref << " and " << ir << std::endl;
            std::cout << rtop.format() << std::endl;
            std::cout << "inverse of reference matcher transformation:" << std::endl;
            std::cout << rtop.inverse().format() << std::endl;

            clipper::Coord_frac box0 = extents_pair.first.coord_frac(start_map.cell());
            clipper::Coord_frac box1 = extents_pair.second.coord_frac(start_map.cell());

            clipper::Grid_map grid(box0.coord_grid(start_map.grid_sampling()),
                                   box1.coord_grid(start_map.grid_sampling()));

            clipper::Xmap_base::Map_reference_coord ix(start_map);
            std::cout << "boxing over " << box0.format()
                      << " to " << box1.format() << std::endl;

            for (int w = grid.min().w(); w <= grid.max().w(); w++) {
               for (int v = grid.min().v(); v <= grid.max().v(); v++) {
                  ix.set_coord(clipper::Coord_grid(grid.min().u(), v, w));
                  for (int u = grid.min().u(); u <= grid.max().u(); u++) {

                     clipper::Coord_grid cg = ix.coord();
                     clipper::Coord_frac cf = cg.coord_frac(start_map.grid_sampling());
                     clipper::Coord_orth co = cf.coord_orth(start_map.cell());
                     clipper::Coord_orth cs = co.transform(rtop.inverse());
                     float f = coot::util::density_at_point(start_map, cs);
                     running[ix] += f;
                     ix.next_u();
                  }
               }
            }
         }

         // Let's create a mask.  We don't want averaged density where
         // it shouldn't be.  e.g. over the B molecule in a A, B ncs
         // system.
         //
         // So the mask is constructed from points that are in the
         // averaging target volume.  1 for average map is valid, 0
         // for average map is not.
         //
         // Let's do something like the previous bit of code, but not
         // getting density from transformed point, simply set the mask flag.
         clipper::Xmap<short int> mask;
         mask.init(start_map.spacegroup(), start_map.cell(), start_map.grid_sampling());
         clipper::Xmap_base::Map_reference_index ixm;
         for (ixm = mask.first(); !ixm.last(); ixm.next() )
            mask[ixm] = 0;

         clipper::Coord_frac box0 = extents_pair.first.coord_frac(start_map.cell());
         clipper::Coord_frac box1 = extents_pair.second.coord_frac(start_map.cell());
         clipper::Grid_map grid(box0.coord_grid(mask.grid_sampling()),
                                box1.coord_grid(mask.grid_sampling()));
         clipper::Xmap_base::Map_reference_coord ix(mask);
         for (int w = grid.min().w(); w <= grid.max().w(); w++) {
            for (int v = grid.min().v(); v <= grid.max().v(); v++) {
               ix.set_coord(clipper::Coord_grid(grid.min().u(), v, w));
               for (int u = grid.min().u(); u <= grid.max().u(); u++) {
                  mask[ix] = 1;
                  ix.next_u();
               }
            }
         }

         // so now we have a map that has been added.  Now lets divide
         // all the points by the number of maps that contribute:
         float factor = 1.0/float(1 + reference_matchers[iref].size());
         // std::cout << "INFO:: There were " << 1/factor
         //            << " maps contributing to the average\n";
         logger.log(log_t::INFO, "There were", 1/factor, "maps contributing to the average");
         // std::cout << "INFO:: rescaling by " << factor << std::endl;
         logger.log(log_t::INFO, "rescaling by", factor);
         clipper::Xmap_base::Map_reference_index irx;
         int n_masked = 0;
         int n_averaged = 0;
         int n_total = 0;
         for (irx = running.first(); !irx.last(); irx.next() ) {
            n_total++;
            if (mask[irx] == 1) { 
               running[irx] *= factor;
               n_averaged++;
            } else { 
               running[irx] = 0;
               n_masked++;
            }
         }
         // std::cout << "INFO:: " << n_averaged << " out of " << n_total
         //            << " (" << 100.0*float(n_averaged)/float(n_total)
         //            << "%) map points " << " were masked out of NCS average target volume, "
         //            << " chain " << reference_ids[iref] << std::endl;
         logger.log(log_t::INFO, std::to_string(n_averaged) + " out of " + std::to_string(n_total) +
                    " (" + std::to_string(100.0f*float(n_averaged)/float(n_total)) +
                    "%) map points were masked out of NCS average target volume, chain " +
                    reference_ids[iref]);
     //std::string name = "NCS average of Chain ";
     std::string name = "Map ";
     name += imol_map_name;
     name += " NCS average of Chain ";
         name += reference_ids[iref];
         name += " type molecules";
         annotated_xmaps.push_back(std::pair<clipper::Xmap<float>, std::string> (running, name));
      } // over each reference id
   } // if have ghosts
   return annotated_xmaps;
}

void
molecule_class_info_t::install_ghost_map(const clipper::Xmap<float> &map_in, std::string name_in,
                                         const coot::ghost_molecule_display_t &ghost_info,
                                         int is_diff_map_flag,
                                         int swap_difference_map_colours_flag,
                                         float contour_level_in) {

   // std::cout << "INFO:: installing ghost map with name :" << name_in << std::endl;
   logger.log(log_t::INFO, "installing ghost map with name:", name_in);

   is_dynamically_transformed_map_flag = 1;
   xmap = map_in;

   bool is_anomalous_flag = false;
   initialize_map_things_on_read_molecule(name_in,
                                          is_diff_map_flag, is_anomalous_flag,
                                          swap_difference_map_colours_flag);
   update_map_in_display_control_widget();
   map_ghost_info = ghost_info;
   
   // fill class variables
   mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);
   map_mean_  = mv.mean;
   map_sigma_ = sqrt(mv.variance);
   contour_level  = contour_level_in;
   update_map(true);

   std::cout << "Done install_ghost_map" << std::endl;
}



void
molecule_class_info_t::add_ncs_ghost(const std::string &chain_id,
                                     const std::string &target_chain_id,
                                     const clipper::RTop_orth &rtop) {

   std::string name = "Manual Operator for Chain ";
   name += chain_id;
   int selHnd = atom_sel.mol->NewSelection();

   coot::ghost_molecule_display_t ghost(rtop, selHnd, name);
   ghost.target_chain_id = target_chain_id;
   ghost.chain_id = chain_id;
   int imod = 1;
   atom_sel.mol->SelectAtoms(selHnd, imod,
                             (char *) chain_id.c_str(),
                             mmdb::ANY_RES, "*",
                             mmdb::ANY_RES, "*", "*", "*", "*", "*");
   ncs_ghosts_have_rtops_flag = 1;

   ghost.update_bonds(atom_sel.mol);
   ncs_ghosts.push_back(ghost);

}

void
molecule_class_info_t::clear_ncs_ghost_matrices() {

   ncs_ghosts.clear();
   ncs_ghosts_have_rtops_flag = 1; 

}


// e.g.  V 1   vs  M 0 
//       T 2       V 1
//       G 3       T 2 
//                 G 3
// In such a case, we want to find a match.
//
bool
molecule_class_info_t::ncs_chains_match_p(const std::vector<std::pair<std::string, int> > &v1,
                                          const std::vector<std::pair<std::string, int> > &v2,
                                          float exact_homology_level,
                                          bool allow_offset_flag) const {

   // std::cout << "DEBUG:: ncs_chains_match_p allow_offset_flag: " << allow_offset_flag << std::endl;
   // std::cout << "   v1 size " << v1.size() << " v2.size() " << v2.size() << std::endl;
   
   bool imatch = 0;
   if (allow_offset_flag == 1) {

      return ncs_chains_match_with_offset_p(v1, v2, exact_homology_level);

   } else { 
      
      // First set max_v1 and min_v1 to the the max and min res numbers
      // from the incoming vector.
      // 
      // Similarly for max_v2 and min_v2.
      //
      // a_offset is the min(min_v1, min_v2)
      // max_index in max(max_v1, max_v2:
      //
      // Now create 2 string vectors, a and b between max_index and
      // a_offset, either "" or "-".
      //
      // Into a, put the v1 sequence, applying offset.
      //
      //
      int min_v1 = 9999, max_v1 = -9999;
      int min_v2 = 9999, max_v2 = -9999;
      if (v1.size() > 0 && v2.size() > 0) {
         for (unsigned int i=0; i<v1.size(); i++) { 
            if (v1[i].second > max_v1)
               max_v1 = v1[i].second;
            if (v1[i].second < min_v1)
               min_v1 = v1[i].second;
         }
         for (unsigned int i=0; i<v2.size(); i++) { 
            if (v2[i].second > max_v2)
               max_v2 = v2[i].second;
            if (v2[i].second < min_v2)
               min_v2 = v2[i].second;
         }

         // 20100601 Beware of the minimum resno being set to mmdb::MinInt4 (a sign of an
         // empty chain - or something, just give up in that case) otherwise (prior to
         // now) -> crash.
         // 
         int a_offset = min_v1;
         int max_indx = max_v1;
         if (min_v2 < min_v1)
            a_offset = min_v2;
         if (max_v2 > max_v1)
            max_indx = max_v2;

         // 20100601 crash protection tests added.
         // 
         if (a_offset != mmdb::MinInt4) { 
            // std::cout << "in ncs_chains_match_p() " << max_indx << " - " << a_offset
            // << " +1 = " << max_indx-a_offset + 1 << std::endl;
            int vec_size = max_indx-a_offset + 1;
            if (vec_size > 0) { 
               std::vector<std::string> a(vec_size, "");
               std::vector<std::string> b(vec_size, "-");
               for (unsigned int i=0; i<v1.size(); i++) {
                  a[v1[i].second-a_offset] = v1[i].first;
               }
               for (unsigned int i=0; i<v2.size(); i++){ 
                  b[v2[i].second-a_offset] = v2[i].first;
                  //          std::cout << " adding residue type :" << v2[i].first
                  //                    << ": to b with index " << v2[i].second-a_offset
                  //                    << " parts " << v2[i].second << " and " << a_offset
                  //                    << " index " << i << " of " << v2.size()
                  //                    << std::endl;
               }

               int n_match = 0;
               for (unsigned int i=0; i<a.size(); i++) {
                  if (a[i] == b[i]) {
                     n_match++;
                  }
               }
               int n_count = a.size();
               if (false)
                  // std::cout << "INFO:: NCS chain comparison " << n_match << "/" << v1.size() << std::endl;
                  logger.log(log_t::INFO, "NCS chain comparison " + std::to_string(n_match) + "/" + std::to_string(v1.size()));
               if (n_count > 0) {
                  // float hit_rate = float(n_match)/float(n_count);
                  // case where protein is 1 to 123 but NAP at 500 fails.  So not n_count but v1.size():
                  if (v1.size() > 0) { 
                     float hit_rate = float(n_match)/float(v1.size());
                     if (hit_rate > exact_homology_level) {
                        imatch = 1;
                     }
                  }
               }
            }
         }
      }
   }
   // std::cout << "DEBUG:: ncs_chains_match_p returns: " << imatch << std::endl;
   return imatch;
}

// Only allow matching if the sequences of v1 and v2 match.  Ignore
// the residue numbering.  If the chains are not alligned (start at
// the same place and have the same number of residues) then we fail.
// Too bad.  Typically, this comes from shelxl molecules, so it is
// more likely the user wil have consitent chains - if not they can
// always fall back to the manual NCS method.
// 
bool
molecule_class_info_t::ncs_chains_match_with_offset_p(const std::vector<std::pair<std::string, int> > &v1,
                                                      const std::vector<std::pair<std::string, int> > &v2,
                                                      float exact_homology_level) const {


   bool match_flag = 0;
   int n_match = 0;
   if (v1.size() > 0) { 
      unsigned int lim = v1.size();
      if (v2.size() < lim)
         lim = v2.size();
      
      for (unsigned int i=0; i<lim; i++) {
         if (v1[i].first == v2[i].first) {
            n_match++;
         }
      }
      float ratio = float(n_match)/float(v1.size());
      // std::cout << "DEBUG:: ncs_chains_match_with_offset_p ratio: " << ratio << std::endl;
      if (ratio > exact_homology_level)
         match_flag = 1;
   }
   // std::cout << "DEBUG:: ncs_chains_match_with_offset_p returns: " << match_flag << std::endl;
   return match_flag;
}



// overwrite the to_chain with a copy of the from_chain, internal private function.
// 
int
molecule_class_info_t::copy_chain(mmdb::Chain *from_chain, mmdb::Chain *to_chain,
                                  clipper::RTop_orth a_to_b_transform) {
   
   int done_copy = 0;
   if (atom_sel.n_selected_atoms > 0) { 
      
      int imod = 1;
      make_backup(__FUNCTION__);
      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
      
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      // std::cout << "INFO:: moving atoms using transformation:\n"
      // << a_to_b_transform.format() << std::endl;
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         if (to_chain == chain_p) { 
            mmdb::Chain *new_chain = new mmdb::Chain;
            new_chain->Copy(from_chain);
            int nres = new_chain->GetNumberOfResidues();
            mmdb::PResidue residue_p;
            for (int ires=0; ires<nres; ires++) { 
               residue_p = new_chain->GetResidue(ires);
               int n_atoms = residue_p->GetNumberOfAtoms();
               mmdb::Atom *at;
               for (int iat=0; iat<n_atoms; iat++) {
                  at = residue_p->GetAtom(iat);
                  clipper::Coord_orth p(at->x, at->y, at->z);
                  clipper::Coord_orth tp = p.transform(a_to_b_transform);
                  at->x = tp.x();
                  at->y = tp.y();
                  at->z = tp.z();
               }
            }
            // Actually, instead of adding a new chain and deleting
            // to_chain, I'd rather replace to_chain, I wonder if I can
            // do that?
            //
            // I don't seem to be able to find the function.  Must
            // talk to Eugene.  (for residues it is InsResidue)
            // 
            new_chain->SetChainID(to_chain->GetChainID());
            model_p->DeleteChain(ichain); // to_chain
            model_p->AddChain(new_chain);
            done_copy = 1;
         }
      }
   }
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   have_unsaved_changes_flag = 1;
   make_bonds_type_checked();
   trim_atom_label_table();
   update_symmetry();
   return done_copy;
}

int
molecule_class_info_t::copy_residue_range(mmdb::Chain *from_chain, mmdb::Chain *to_chain,
                                          int residue_range_1,
                                          int residue_range_2,
                                          clipper::RTop_orth a_to_b_transform) {

   // Run through the residues in the residue range of from_chain
   // (using for loop over residues). It would be a mistake to delete
   // the "to" residue and make a copy of the "from" residue and
   // insert/append it because we then get strange sequence/selection
   // problems.
   //
   // So we look for that residue in the "to" chain.  Delete all its
   // atoms and add a transformed copy for each of the atoms in the
   // corresponding "from" residue.
   //
   // If we can't find that residue in the to_chain, then do a
   // deep_copy of the residue and add insert/append it.
   // 
   
   int done_copy = 0;
   if (atom_sel.n_selected_atoms > 0) {

      if (residue_range_1 > residue_range_2) { 
         int t = residue_range_2;
         residue_range_2 = residue_range_1;
         residue_range_1 = t;
      }
         
      make_backup(__FUNCTION__);
      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);

      std::string old_seg_id_for_chain_atoms;
      bool use_old_seg_id = 0;
      try {
         old_seg_id_for_chain_atoms = coot::chain_atoms_segid(to_chain);
         use_old_seg_id = 1;
      }
      catch (const std::runtime_error &mess) {
      }
      
      // Don't do a selection here, we have the (from) chain already.
      int n_from_chain_residues = from_chain->GetNumberOfResidues();
      for (int ifc=0; ifc<n_from_chain_residues; ifc++) {

         mmdb::Residue *from_residue = from_chain->GetResidue(ifc);
         int resno = from_residue->GetSeqNum();
         char *ins_code = from_residue->GetInsCode();

         if ((resno>=residue_range_1) && (resno<=residue_range_2)) {

            mmdb::PAtom *from_residue_atoms;
            int n_from_residue_atoms;
            
            // OK, copy this residue:
            
            // Now, does there exist a corresponding residue in the "to" chain?
            // Do a residue selection:

            int nSelResidues = 0;
            int selHnd_res = atom_sel.mol->NewSelection();
            mmdb::PPResidue SelResidues;
            atom_sel.mol->Select(selHnd_res, mmdb::STYPE_RESIDUE, 1,
                                 to_chain->GetChainID(),
                                 resno, ins_code,
                                 resno, ins_code,
                                 "*",  // residue name
                                 "*",  // Residue must contain this atom name?
                                 "*",  // Residue must contain this Element?
                                 "*",  // altLocs
                                 mmdb::SKEY_NEW // selection key
                                 );
            atom_sel.mol->GetSelIndex(selHnd_res, SelResidues, nSelResidues);

            mmdb::Residue *to_residue = 0;
            if (nSelResidues > 0) {
               to_residue = SelResidues[0];
            }
            // delete the selection before we mess with it: (maybe
            // not needed because it is a residue selection and we
            // are not deleting residues)
            atom_sel.mol->DeleteSelection(selHnd_res);

            if (nSelResidues > 0) {
               
               // OK, that residue existed.

               mmdb::PAtom *old_to_residue_atoms;
               int n_old_to_residue_atoms;
               to_residue->GetAtomTable(old_to_residue_atoms, n_old_to_residue_atoms);
               for (int iat=0; iat<n_old_to_residue_atoms; iat++) {
                  if (to_residue->GetAtom(iat)) {
                     to_residue->DeleteAtom(iat);
                  } else {
                     std::cout << "ERROR:: Trapped error, null atom in to_residue!? "
                               << iat << std::endl;
                  } 
               }
               to_residue->TrimAtomTable();

               from_residue->GetAtomTable(from_residue_atoms, n_from_residue_atoms);
               for (int iat=0; iat<n_from_residue_atoms; iat++) {
                  mmdb::Atom *at = new mmdb::Atom;
                  if (from_residue_atoms[iat]) { 
                     at->Copy(from_residue_atoms[iat]);
                     to_residue->AddAtom(at);
                     if (use_old_seg_id)
                        strcpy(at->segID, old_seg_id_for_chain_atoms.c_str());
                  } else {
                     std::cout << "ERROR:: Trapped error, null from_residue_atoms "
                               << iat << std::endl;
                  } 
               }
               to_residue->TrimAtomTable();

               // PRE bug.  If the to_residue was not the same residue
               // type as the from_residue, then the to_residue needs
               // to have its residue type updated.
               std::string to_residue_name   =   to_residue->GetResName();
               std::string from_residue_name = from_residue->GetResName();
               if (to_residue_name != from_residue_name) {
                  to_residue->SetResName(from_residue->GetResName());
               } 
               
               if ((n_old_to_residue_atoms > 0) || (n_from_residue_atoms > 0))
                  atom_sel.mol->FinishStructEdit();
               
            } else {

               // There was no such residue already existing in the
               // "to" chain.  Let's create/add one.  Perhaps I need
               // InsResidue here?.
               //
               // 20080401.  Yes.  I do.  Needed to pass NCS residue range test.
               
               to_residue = coot::util::deep_copy_this_residue(from_residue);

               std::pair<int, mmdb::Residue *> serial_number =
                  find_serial_number_for_insert(to_residue->GetSeqNum(),
                                                to_residue->GetInsCode(),
                                                to_chain->GetChainID());

               if (false)
                  std::cout << "debug:: found serial_number " << serial_number.first << " "
                            << coot::residue_spec_t(serial_number.second)
                            << " for residue " << coot::residue_spec_t(to_residue)
                            << " with index " << to_residue->index
                            << std::endl;

               if (serial_number.first != -1) {
                  to_chain->InsResidue(to_residue, serial_number.first);
                  coot::copy_segid(serial_number.second, to_residue);
               } else {
                  mmdb::Residue *res = last_residue_in_chain(to_chain);
                  to_chain->AddResidue(to_residue);
                  coot::copy_segid(res, to_residue);
               }
               atom_sel.mol->FinishStructEdit();
            }
            
            if (to_residue) {
               // Now we need to transform the atoms of to_residue:
               mmdb::PAtom *to_residue_atoms;
               int n_to_residue_atoms;
               to_residue->GetAtomTable(to_residue_atoms, n_to_residue_atoms);
               for (int iat=0; iat<n_to_residue_atoms; iat++) {
                  clipper::Coord_orth a(to_residue_atoms[iat]->x,
                                        to_residue_atoms[iat]->y,
                                        to_residue_atoms[iat]->z);
                  clipper::Coord_orth a_tr = a.transform(a_to_b_transform);
                  to_residue_atoms[iat]->x = a_tr.x();
                  to_residue_atoms[iat]->y = a_tr.y();
                  to_residue_atoms[iat]->z = a_tr.z();
               }
            } else {
               std::cout << "ERROR:: Null to_residue!!!!" << std::endl;
            }
         }
      }
      
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
      trim_atom_label_table();
      update_symmetry();
   }
   return done_copy;

}


// and the other way for CNS NCS users
//
void
molecule_class_info_t::add_strict_ncs_matrix(const std::string &chain_id,
                                             const std::string &target_chain_id,
                                             const coot::coot_mat44 &m) {

   // std::cout << "-------------------------------------------------------------- add_strict_ncs_matrix imol "
   // << imol_no << " " << chain_id << " " << target_chain_id << std::endl;


   std::string name = "Strict NCS for Chain ";
   name += chain_id;
   coot::ghost_molecule_display_t ghost;
   ghost.name = name;
   ghost.chain_id = chain_id;
   ghost.target_chain_id = target_chain_id;
   strict_ncs_info.push_back(ghost);
   strict_ncs_matrices.push_back(m);
}


void
molecule_class_info_t::update_strict_ncs_symmetry(const coot::Cartesian &centre_point,
                                                  const molecule_extents_t &extents) {

   bool debug = false;

   if (debug) {
      std::cout << "DEBUG:: Update ncs symmetry for " << strict_ncs_matrices.size()
                << " NCS matrices" << std::endl;
      for (const auto &m : strict_ncs_matrices) {
         std::cout << "NCS Matrix\n"
                   << "[ " << m.m[0].v4[0] << " " << m.m[0].v4[1] << " " << m.m[0].v4[2] << " " << m.m[0].v4[3] << " ]\n"
                   << "[ " << m.m[1].v4[0] << " " << m.m[1].v4[1] << " " << m.m[1].v4[2] << " " << m.m[1].v4[3] << " ]\n"
                   << "[ " << m.m[2].v4[0] << " " << m.m[2].v4[1] << " " << m.m[2].v4[2] << " " << m.m[2].v4[3] << " ]\n"
                   << "[ " << m.m[3].v4[0] << " " << m.m[3].v4[1] << " " << m.m[3].v4[2] << " " << m.m[3].v4[3] << " ]\n";
      }
   }

   // We need to convert from internal coot_mat44 to mmdb::mat44s, then do
   // similar things to update_symmetry()

   Cell_Translation c_t = extents.coord_to_unit_cell_translations(centre_point, atom_sel);

   if (debug) {
      std::cout << "cell translation " << c_t << std::endl;
   }

   // For real NCS, uses the unit cell
   //
   std::vector<std::pair<int, symm_trans_t> > ncs_mat_indices =
      extents.which_strict_ncs(centre_point, atom_sel, strict_ncs_matrices, c_t);

   // guaranteed to be at least one.
   if (debug)
      std::cout << "There were " << ncs_mat_indices.size() << " touching molecules\n";

   std::vector<std::pair<coot::coot_mat44, symm_trans_t> > cmats;
   for (unsigned int i=0; i<ncs_mat_indices.size(); i++)
      cmats.push_back(std::pair<coot::coot_mat44, symm_trans_t> (strict_ncs_matrices[ncs_mat_indices[i].first],
                                                                 ncs_mat_indices[i].second));

   if (debug) {
      for (unsigned int i=0; i<cmats.size(); i++) {
         const coot::coot_mat44 &m44 = cmats[i].first;
         const symm_trans_t &st      = cmats[i].second;
         std::cout << " cmat " << i <<  std::endl;
         std::cout << m44;
         std::cout << st;
      }
   }

   Bond_lines_container bonds;
   strict_ncs_bonds_box.clear();
   strict_ncs_bonds_box =
      bonds.add_NCS(atom_sel, imol_no,
                    centre_point,
                    graphics_info_t::symmetry_search_radius,
                    cmats,
                    symmetry_as_calphas,
                    symmetry_whole_chain_flag);

}



// public interface to chain copying
void
molecule_class_info_t::copy_chain(const std::string &from_chain_str,
                                  const std::string &to_chain_str) {

   if (atom_sel.n_selected_atoms > 0) {
      if (ncs_ghosts.size() > 0) {
         if (ncs_ghosts[0].is_empty() || ncs_ghosts_have_rtops_flag == 0) {
            // std::cout << "   %%%%%%%%% calling fill_ghost_info from copy_chain "
            //                       << std::endl;
            fill_ghost_info(1, 0.7); // 0.7?
         }
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            if  (ncs_ghosts[ighost].chain_id == to_chain_str) {
               if (ncs_ghosts[ighost].target_chain_id == from_chain_str) {
                  clipper::RTop_orth rtop = ncs_ghosts[ighost].rtop.inverse();
                  mmdb::PPAtom atom_selection = NULL;
                  int n_atoms = 0;
                  atom_sel.mol->GetSelIndex(ncs_ghosts[ighost].SelectionHandle,
                                            atom_selection, n_atoms);
                  if (n_atoms > 0) {
                     mmdb::Chain *to_chain = atom_selection[0]->GetChain();
                     // OK, now let's find the from_chain:

                     int imod = 1;
                     mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
                     mmdb::Chain *chain_p;
                     // run over chains of the existing mol
                     int nchains = model_p->GetNumberOfChains();
                     for (int ichain=0; ichain<nchains; ichain++) {
                        chain_p = model_p->GetChain(ichain);
                        std::string chain_id(chain_p->GetChainID());
                        if (from_chain_str == chain_id) {
                           mmdb::Chain *from_chain = chain_p;
                           copy_chain(from_chain, to_chain, rtop);
                           break;
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


std::ostream& coot::operator<<(std::ostream &s, const coot::ghost_molecule_display_t &ghost) {
   s << "[GHOST Name: " << ghost.name << ", chain-id: " << ghost.chain_id << ", target-chain-id: "
     << ghost.target_chain_id << ", displayed? " << ghost.display_it_flag << "]";
   return s;
}



int
molecule_class_info_t::copy_residue_range_from_ncs_master_to_other_using_ghost(std::string from_chain_id,
                                                                               std::string to_chain_id,
                                                                               int residue_range_1,
                                                                               int residue_range_2) {

   int done_copy_flag = 0;
   if (atom_sel.n_selected_atoms > 0) {
      if (ncs_ghosts.size() > 0) {
         if (ncs_ghosts[0].is_empty() || ncs_ghosts_have_rtops_flag == 0) {
            // std::cout << "   %%%%%%%%% calling fill_ghost_info from "
            // std::cout << "copy_residue_range_from_ncs_master_... "
            // << std::endl;
            fill_ghost_info(1, 0.7); // 0.7?
         }
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            if  (ncs_ghosts[ighost].chain_id == to_chain_id) {
               if (ncs_ghosts[ighost].target_chain_id == from_chain_id) {
                  clipper::RTop_orth rtop = ncs_ghosts[ighost].rtop.inverse();
                  mmdb::PPAtom atom_selection = NULL;
                  int n_atoms;
                  atom_sel.mol->GetSelIndex(ncs_ghosts[ighost].SelectionHandle,
                                            atom_selection, n_atoms);
                  if (n_atoms > 0) {
                     mmdb::Chain *to_chain = atom_selection[0]->GetChain();
                     // OK, now let's find the from_chain:

                     int imod = 1;
                     mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
                     mmdb::Chain *chain_p;
                     // run over chains of the existing mol
                     int nchains = model_p->GetNumberOfChains();
                     for (int ichain=0; ichain<nchains; ichain++) {
                        chain_p = model_p->GetChain(ichain);
                        std::string chain_id(chain_p->GetChainID());
                        if (from_chain_id == chain_id) {
                           mmdb::Chain *from_chain = chain_p;
                           if (false) {
                              std::cout << "Doing a copy_residue range "
                                      << from_chain->GetChainID() << " to "
                                        << to_chain->GetChainID() << " " << residue_range_1
                                        << " to "
                                        << residue_range_2 << "\n"
                                        << rtop.format() << std::endl;
                           }
                           copy_residue_range(from_chain, to_chain, residue_range_1,
                                              residue_range_2, rtop);
                           done_copy_flag = 1;
                           break;
                        }
                     }
                  }
                  if (done_copy_flag)
                     break;
               }
            }
         }
      }
   }
   return done_copy_flag;
}


int
molecule_class_info_t::copy_from_ncs_master_to_others(const std::string &master_chain_id) {

   // check in the ghosts if master_chain_id is actually a master
   // molecule, and if it is, apply the copy to all ghosts for which
   // it is a master.
   int ncopied = 0;
   if (atom_sel.n_selected_atoms > 0) {
      if (ncs_ghosts.size() > 0) {
         if (ncs_ghosts[0].is_empty() || ncs_ghosts_have_rtops_flag == 0) {
            // std::cout << "   %%%%%%%%% calling fill_ghost_info from "
            // std::cout << "copy_from_ncs_master_to_others "
            // << std::endl;
            short int do_rtops = 1;
            float homology_lev = 0.7;
            fill_ghost_info(do_rtops, homology_lev);
         }
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            std::string master = ncs_ghosts[ighost].target_chain_id;
            if (master == master_chain_id) {
               copy_chain(master, ncs_ghosts[ighost].chain_id);
            }
         }
      }
   }
   return ncopied;
}

int
molecule_class_info_t::copy_from_ncs_master_to_specific_other_chains(const std::string &master_chain_id,
                                                                     const std::vector<std::string> &other_chain_ids) {

   // check in the ghosts if master_chain_id is actually a master
   // molecule, and if it is, apply the copy to all ghosts for which
   // it is a master.
   int ncopied = 0;
   if (atom_sel.n_selected_atoms > 0) {
      if (ncs_ghosts.size() > 0) {
         if (ncs_ghosts[0].is_empty() || ncs_ghosts_have_rtops_flag == 0) {
            // std::cout << "   %%%%%%%%% calling fill_ghost_info from "
            // std::cout << "copy_from_ncs_master_to_others "
            // << std::endl;
            fill_ghost_info(1, 0.7); // 0.7?
         }
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            const std::string this_chain_id = ncs_ghosts[ighost].chain_id;
            std::string master_for_this_chain = ncs_ghosts[ighost].target_chain_id;
            if (master_for_this_chain == master_chain_id) {
               if (std::find(other_chain_ids.begin(), other_chain_ids.end(), this_chain_id) != other_chain_ids.end()) {
                  copy_chain(master_for_this_chain, this_chain_id);
               }
            }
         }
      }
   }
   return ncopied;
}

int
molecule_class_info_t::copy_residue_range_from_ncs_master_to_others(const std::string &master_chain_id,
                                                                    int residue_range_1,
                                                                    int residue_range_2) {

   // check in the ghosts if master_chain_id is actually a master
   // molecule, and if it is, apply the copy to all ghosts for which
   // it is a master.
   int ncopied = 0;
   if (atom_sel.n_selected_atoms > 0) {
      if (ncs_ghosts.size() > 0) {
         if (ncs_ghosts[0].is_empty() || ncs_ghosts_have_rtops_flag == 0) {
            // std::cout << "   %%%%%%%%% calling "
            // copy_residue_range_from_ncs_master_to_others "
            // << std::endl;
            fill_ghost_info(1, 0.7); // 0.7?
         }
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            std::string master = ncs_ghosts[ighost].target_chain_id;
            if (master == master_chain_id) {
               copy_residue_range_from_ncs_master_to_other_using_ghost(master_chain_id,
                                                                       ncs_ghosts[ighost].chain_id,
                                                                       residue_range_1,
                                                                       residue_range_2);
               ncopied++;
            }
         }
      }
      if (ncopied == 0) {
         std::cout << "WARNING:: failed to find master_chain_id \"" << master_chain_id
                   << "\" in " << ncs_ghosts.size() << " NCS ghosts" << std::endl;
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            std::cout << "    ghost: chainid \"" << ncs_ghosts[ighost].chain_id
                      << "\" has NCS master \"" << ncs_ghosts[ighost] << "\""
                      << std::endl;
         }
      }
   }
   return ncopied;
}

int
molecule_class_info_t::copy_residue_range_from_ncs_master_to_chains(const std::string &master_chain_id,
                                                                    int residue_range_1, int residue_range_2,
                                                                    const std::vector<std::string> &chain_ids) {

   // Check in the ghosts if master_chain_id is actually a master
   // molecule, and if it is, apply the copy to all ghosts for which
   // it is a master.
   int ncopied = 0;
   if (atom_sel.n_selected_atoms > 0) {
      if (ncs_ghosts.size() > 0) {
         if (ncs_ghosts[0].is_empty() || ncs_ghosts_have_rtops_flag == 0) {
            fill_ghost_info(1, 0.7); // 0.7?
         }
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            std::string master = ncs_ghosts[ighost].target_chain_id;
            
            if (std::find(chain_ids.begin(), chain_ids.end(), ncs_ghosts[ighost].chain_id) !=
                chain_ids.end()) { 
               if (master == master_chain_id) {
                  copy_residue_range_from_ncs_master_to_other_using_ghost(master_chain_id,
                                                                          ncs_ghosts[ighost].chain_id,
                                                                          residue_range_1,
                                                                          residue_range_2);
                  ncopied++;
               }
            }
         }
      }
   }
   return ncopied;
} 

int
molecule_class_info_t::copy_from_ncs_master_to_chains(const std::string &master_chain_id,
                                                                    const std::vector<std::string> &chain_ids) {

   // Check in the ghosts if master_chain_id is actually a master
   // molecule, and if it is, apply the copy to all ghosts for which
   // it is a master.
   int ncopied = 0;
   if (atom_sel.n_selected_atoms > 0) {
      if (ncs_ghosts.size() > 0) {
         if (ncs_ghosts[0].is_empty() || ncs_ghosts_have_rtops_flag == 0) {
            fill_ghost_info(1, 0.7); // 0.7?
         }
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            std::string master = ncs_ghosts[ighost].target_chain_id;
            
            if (std::find(chain_ids.begin(), chain_ids.end(), ncs_ghosts[ighost].chain_id) !=
                chain_ids.end()) { 
               if (master == master_chain_id) {
                  copy_chain(master_chain_id,
                     ncs_ghosts[ighost].chain_id);
                  ncopied++;
               }
            }
         }
      }
   }
   return ncopied;
} 




int
molecule_class_info_t::set_ncs_master_chain(const std::string &new_master_chain_id, float homology_lev) {

   int retval = 0;
   std::vector<std::string> chain_ids;
   std::vector<std::vector<std::pair<std::string, int> > > residue_types;
   std::vector<int> chain_atom_selection_handles;
   std::vector<short int> first_chain_of_this_type;

   // start from a blank slate:
   ncs_ghosts.resize(0);
   ncs_ghosts_have_rtops_flag = 0;

   if (atom_sel.n_selected_atoms > 0) {

      int n_models = atom_sel.mol->GetNumberOfModels();
      if (n_models > 0) {
         int imod = 1; // otherwise madness

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         if (nchains <= 0) {
            std::cout << "bad nchains in molecule " << nchains
                      << std::endl;
         } else {
            chain_ids.resize(nchains);
            residue_types.resize(nchains);
            chain_atom_selection_handles.resize(nchains);
            first_chain_of_this_type.resize(nchains, 1);
            for (int ichain=0; ichain<nchains; ichain++) {
               chain_p = model_p->GetChain(ichain);
               if (! chain_p->isSolventChain()) {
                  chain_ids[ichain] = chain_p->GetChainID();
                  int iselhnd = atom_sel.mol->NewSelection();
                  mmdb::PAtom *atom_selection = NULL;
                  int nSelAtoms;
                  atom_sel.mol->SelectAtoms(iselhnd, imod,
                                            chain_p->GetChainID(),
                                            mmdb::ANY_RES, "*",
                                            mmdb::ANY_RES, "*",
                                            "*", "*", "*", "*");
                  atom_sel.mol->GetSelIndex(iselhnd, atom_selection, nSelAtoms);
                  chain_atom_selection_handles[ichain] = iselhnd;

                  int nres = chain_p->GetNumberOfResidues();
                  residue_types[ichain].resize(nres);
//                   std::cout << "INFO:: residues_types[" << ichain << "] resized to "
//                             << residue_types[ichain].size() << std::endl;
                  mmdb::PResidue residue_p;
                  for (int ires=0; ires<nres; ires++) { 
                     residue_p = chain_p->GetResidue(ires);
                     std::string resname(residue_p->name);
                     residue_types[ichain][ires] =
                        std::pair<std::string, int> (resname, residue_p->seqNum);
                  }
               }
            }
         }
      }
      add_ncs_ghosts_using_ncs_master(new_master_chain_id, chain_ids, residue_types,
                                      chain_atom_selection_handles, homology_lev);

      if (! ncs_ghosts.empty()) {
         update_ghosts();
         // std::cout << "INFO:: set_ncs_master_chain Constructed " << ncs_ghosts.size() << " ghosts\n";
         logger.log(log_t::INFO, "set_ncs_master_chain Constructed", ncs_ghosts.size(), "ghosts");
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            std::cout << "   Ghost info:: " << ncs_ghosts[ighost].name << std::endl;
         }
      }
   }
   return retval;
}

void
molecule_class_info_t::set_display_ncs_ghost_chain(int ichain, int state) {

   float ncs_ghost_similarity_score = 0.7; // make a class member datum

   set_show_ghosts(state);

   if (state)
      update_ncs_ghosts();

   // The ichain that is passed is the index number of the chains in
   // the atom_selection.
   //
   // So we need to convert that to the index of ncs_ghosts which has
   // the same chain_id

//     std::cout << "%%%% starting set_display_ncs_ghost_chain" << std::endl;
//     std::cout << " %%%%%%%%% and here/now show_ghosts_flag is "
//               << show_ghosts_flag << " ichain: " << ichain
//                << " state: " << state << std::endl;

   //    std::cout << "   DEBUG:: start of set_display_ncs_ghost_chain: " << std::endl;
   //    std::cout << "        There are " << ncs_ghosts.size() << " ghosts" << std::endl;

//    for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
//       std::cout << "         ighost: " << ighost<< "\n"
//                 << "         name: \""           << ncs_ghosts[ighost].name << "\"" << "\n"
//                 << "         chainid: "         << ncs_ghosts[ighost].chain_id << "\n"
//                 << "         target chain id: " << ncs_ghosts[ighost].target_chain_id<< "\n"
//                 << "         display_it_flag "  << ncs_ghosts[ighost].display_it_flag << std::endl;
//    }

   int ghost_index = -1;
   if (atom_sel.n_selected_atoms > 0) {
      if (show_ghosts_flag) {
         if (ncs_ghosts.size() > 0) {
            if (ncs_ghosts[0].is_empty() || ncs_ghosts_have_rtops_flag == 0) {

               std::cout << "        --  set_display_ncs_ghost_chain calls fill ghost info, 1"
                         << " with ncs_ghosts_have_rtops_flag " << ncs_ghosts_have_rtops_flag
                         << std::endl;

               fill_ghost_info(1, ncs_ghost_similarity_score); // 0.7?
            }
         }
      }

      std::vector<std::string> chain_ids = coot::util::chains_in_molecule(atom_sel.mol);
      if (ichain < int(chain_ids.size())) {
         for (unsigned int ich=0; ich<ncs_ghosts.size(); ich++) {
//              std::cout << "getting ghost_index: chain_ids[" << ichain << "] is "
//                        << chain_ids[ichain] << " and ncs_ghosts[" << ich
//                        << "].chain_id is " << ncs_ghosts[ich].chain_id <<  std::endl;
            if (chain_ids[ichain] == ncs_ghosts[ich].chain_id) {
//                 std::cout << "    DEBUG chain id comparison :" << chain_ids[ichain]
//                           << ": vs :" <<  ncs_ghosts[ich].chain_id << ":" << std::endl;
               ghost_index = ich;
               break;
            }
         }

         //          std::cout << "   Here ghost_index is " << ghost_index << std::endl;
         if (ghost_index > -1 ) {
            if (int(ncs_ghosts.size()) > ghost_index) {
               ncs_ghosts[ghost_index].display_it_flag = state;
            }
         }
      } else {
         // should never happen (except for user scripting)
         std::cout << "ERROR:: out of range ichain " << ichain << std::endl;
      }
   }

   //    std::cout << "   DEBUG:: end of set_display_ncs_ghost_chain: " << std::endl;
   // std::cout << "INFO:: There are " << ncs_ghosts.size() << " ghosts" << std::endl;
   logger.log(log_t::INFO, "There are", ncs_ghosts.size(), "ghosts");
   for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
      std::cout << "         ighost: " << ighost<< "\n"
                <<  "        name: \""           << ncs_ghosts[ighost].name << "\"" << "\n"
                << "         chainid: "         << ncs_ghosts[ighost].chain_id << "\n"
                << "         target chain id: " << ncs_ghosts[ighost].target_chain_id<< "\n"
                << "         display_it_flag "  << ncs_ghosts[ighost].display_it_flag << std::endl;
   }
}


// Return a new orienation, to be used to set the view orientation/quaternion.
// 
std::pair<bool, clipper::RTop_orth>
molecule_class_info_t::apply_ncs_to_view_orientation(const clipper::Mat33<double> &current_view_mat,
                                                     const clipper::Coord_orth &current_position,
                                                     const std::string &current_chain,
                                                     const std::string &next_ncs_chain,
                                                     bool forward_flag) const {

   if (forward_flag) 
      return apply_ncs_to_view_orientation_forward(current_view_mat,
                                                   current_position,
                                                   current_chain,
                                                   next_ncs_chain);
   else
      return apply_ncs_to_view_orientation_backward(current_view_mat,
                                                    current_position,
                                                    current_chain,
                                                    next_ncs_chain);
      
}

// Return a new orienation, to be used to set the view orientation/quaternion.
// 
std::pair<bool, clipper::RTop_orth>
molecule_class_info_t::apply_ncs_to_view_orientation_forward(const clipper::Mat33<double> &current_view_mat,
                                                             const clipper::Coord_orth &current_position,
                                                             const std::string &current_chain,
                                                             const std::string &next_ncs_chain) const { 

   clipper::Mat33<double> r = current_view_mat;
   bool apply_it = 0;
   clipper::Coord_orth t(0,0,0);

   unsigned int n_ghosts = ncs_ghosts.size();

   if ((n_ghosts > 0) && ncs_ghosts_have_rtops_flag)  {

      // If current_chain is not a target_chain_id
      //    { i.e. we were not sitting on an NCS master}
      //    is there a ghost that has chain_id current_chain?
      //    If so, note its target_chain_id {NCS Master}.

      //    if next_ncs_chain (the one we want to jump to) is not
      //    the target_chain_id then {
      //       Then look for a ghost that that chain_id next_ncs_chain
      //       Note its target_chain_id.
      //       if target_chain_ids match
      //          we can find the matrix
      //    } else { next_ncs_chain *was* target_chain_id }
      //          the ncs matrix is the ghost matrix
      //    }
      //
      // else { were were sitting on a NCS master }
      //
      //     look at the next ghost
      //     if the target chain for that is current_chain
      //     then we want the inverse matrix of that next
      //     ghost

      if (! ncs_ghost_chain_is_a_target_chain_p(current_chain)) {

         //    Is there a ghost that has chain_id current_chain?
         //    If so, note its target_chain_id.
         std::string found_target_of_current_chain = "not-found";
         int i_ghost_chain_match = -1;
         for (unsigned int ighost=0; ighost<n_ghosts; ighost++) {
            if (ncs_ghosts[ighost].chain_id == current_chain) {
               found_target_of_current_chain = ncs_ghosts[ighost].target_chain_id;
               i_ghost_chain_match = ighost;
               break;
            }
         }

         if (i_ghost_chain_match != -1) {
            // should always happen

            // were we on the last ghost that has the same
            // target_chain_id as this ghost has?

            if (last_ghost_matching_target_chain_id_p(i_ghost_chain_match, ncs_ghosts)) {

               // we need to go to the target_chain
               clipper::RTop_orth ncs_mat = ncs_ghosts[i_ghost_chain_match].rtop;
                // std::cout << "DEBUG:: Last ghost from " << current_chain << " to target chain "
               // << ncs_ghosts[i_ghost_chain_match].target_chain_id << std::endl;
               r = ncs_mat.rot() * r;
               t = current_position.transform(ncs_mat);
               apply_it = 1;
            } else {
                // std::cout << "DEBUG:: Not last ghost from " << current_chain
               // << " to target chain "
               // << ncs_ghosts[i_ghost_chain_match].target_chain_id << std::endl;
               clipper::RTop_orth ncs_mat_1 = ncs_ghosts[i_ghost_chain_match  ].rtop;
               clipper::RTop_orth ncs_mat_2 = ncs_ghosts[i_ghost_chain_match+1].rtop;
               r = ncs_mat_2.rot().inverse() * (ncs_mat_1.rot() * r);
               t = current_position.transform(ncs_mat_1).transform(ncs_mat_2.inverse());
               apply_it = 1;
            }
         } else {
            std::cout << "ERROR:: An NCS reference chain finding error has occured"
                      << std::endl;
         }

      } else {
         if (ncs_ghosts_have_rtops_flag) {
            // we were sitting on an NCS master

            // find a ghost that has this current_chain_id as a target_chain_id
            clipper::Mat33<double> ncs_mat = ncs_ghosts[0].rtop.rot();
            // try to override that by the right operator:
            for (unsigned int ighost=0; ighost<n_ghosts; ighost++) {
               if (ncs_ghosts[ighost].target_chain_id == current_chain) {
//                   std::cout << "debug:: setting t to inverse of "
//                             << current_chain << " to " << ncs_ghosts[ighost].chain_id
//                             << std::endl;
                  ncs_mat = ncs_ghosts[ighost].rtop.rot();
                  t = current_position.transform(ncs_ghosts[ighost].rtop.inverse());
                  break;
               }
            }
            r = ncs_mat.inverse() * r;
            apply_it = 1;
         }
      } 
   } else {
      std::cout << "WARNING:: no ghosts (with ncs)" << std::endl;
   }
   clipper::RTop_orth rtop(r, t);
   return std::pair<bool, clipper::RTop_orth> (apply_it, rtop);
}

// Return a new orienation, to be used to set the view orientation/quaternion.
// 
// This is for going for example (B->A or C->A with A,B,C,D NCS)
// 
std::pair<bool, clipper::RTop_orth>
molecule_class_info_t::apply_ncs_to_view_orientation_backward(const clipper::Mat33<double> &current_view_mat,
                                                             const clipper::Coord_orth &current_position,
                                                             const std::string &current_chain,
                                                             const std::string &next_ncs_chain) const { 

   clipper::Mat33<double> r = current_view_mat;
   bool apply_it = 0;
   clipper::Coord_orth t(0,0,0);

   unsigned int n_ghosts = ncs_ghosts.size();
   if ((n_ghosts > 0) && (ncs_ghosts_have_rtops_flag))  {

      // If current_chain is not a target_chain_id
      //    { i.e. we were not sitting on an NCS master}
      //    is there a ghost that has chain_id current_chain?
      //    If so, note its target_chain_id {NCS Master}.

      //    if next_ncs_chain (the one we want to jump to) is not
      //    the target_chain_id then {
      //       Then look for a ghost that chain_is next_ncs_chain
      //       Note its target_chain_id.
      //       if target_chain_ids match 
      //          we can find the matrix
      //    } else { next_ncs_chain *was* target_chain_id }
      //          the ncs matrix is the ghost matrix
      //    }
      //
      // else { were were sitting on a NCS master }
      //
      //     look at the "next" ghost, starting from the back
      //     if the target chain for that is current_chain
      //     then we want the inverse matrix of that next
      //     ghost

      if (! ncs_ghost_chain_is_a_target_chain_p(current_chain)) {
         // find the ghost that has current_chain as chain_id
         int i_ghost_chain_match = -1; // unset
         for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
            if (ncs_ghosts[ighost].chain_id == current_chain) {
               i_ghost_chain_match = ighost;
               break;
            }
         }

         if (i_ghost_chain_match != -1) {
            // OK, we found (as expected) a ghost that has a chain_id
            // that is the current_chain.

            // Now, imagine we have A,B,C,D NCS ("A" is the master).
            // The B->A transformation is trivial, the A->D
            // transformation is the inverse of the ghost and D->C and
            // D->B involve 2 operations (to A then on again).
            //
            // We need to now check for the B->A case
            // 
            if (ncs_ghosts[i_ghost_chain_match].target_chain_id == next_ncs_chain) {
               clipper::RTop_orth ncs_mat = ncs_ghosts[i_ghost_chain_match].rtop;
               r = ncs_mat.rot() * r;
               t = current_position.transform(ncs_mat);
               apply_it = 1;

            } else {

               // D->C or C->B

               // OK, now what is the closest (going backward) ghost
               // that has the same target chain id as the
               // i_ghost_chain_match?
               std::string tc_id = ncs_ghosts[i_ghost_chain_match  ].target_chain_id;
               int prior_ghost_with_same_target_chain_id = -1;
               for (int ighost=(ncs_ghosts.size()-1); ighost>=0; ighost--) {
                  if (ighost<i_ghost_chain_match) {
                     if (ncs_ghosts[ighost].target_chain_id == tc_id) {
                        prior_ghost_with_same_target_chain_id = ighost;
                        break;
                     }
                  }
               }

               if (prior_ghost_with_same_target_chain_id != -1) {

                  clipper::RTop_orth ncs_mat_1 = ncs_ghosts[i_ghost_chain_match].rtop;
                  clipper::RTop_orth ncs_mat_2 = ncs_ghosts[prior_ghost_with_same_target_chain_id].rtop;

                  r = ncs_mat_2.rot().inverse() * (ncs_mat_1.rot() * r);
                  t = current_position.transform(ncs_mat_1).transform(ncs_mat_2.inverse());
                  apply_it = 1;
               }
            }
         }

      } else {

         // find a ghost that has this current_chain_id as a target_chain_id
         clipper::Mat33<double> ncs_mat = ncs_ghosts[0].rtop.rot();
         // try to override that by the right operator:
         // search from the back for the "next" NCS chain
         // that has target_chain id as current_chain
         for (int ighost=(n_ghosts-1); ighost>=0; ighost--) {
            if (ncs_ghosts[ighost].target_chain_id == current_chain) {
               ncs_mat = ncs_ghosts[ighost].rtop.rot();
               t = current_position.transform(ncs_ghosts[ighost].rtop.inverse());
               break;
            }
         }
         r = ncs_mat.inverse() * r;
         apply_it = 1;
      }
   }
   clipper::RTop_orth rtop(r, t);
   if (! apply_it)
      std::cout << "WARNING:: apply_ncs_to_view_orientation_backward() apply_it not set " << std::endl;
   return std::pair<bool, clipper::RTop_orth> (apply_it, rtop);
}


bool
molecule_class_info_t::ncs_ghost_chain_is_a_target_chain_p(const std::string &chain_id) const {

   bool r = 0;
   unsigned int n_ghosts = ncs_ghosts.size();
   if (n_ghosts > 0) {
      for (unsigned int ighost=0; ighost<n_ghosts; ighost++) {
         if (ncs_ghosts[ighost].target_chain_id == chain_id) {
            r = 1;
            break;
         }
      }
   }
   return r;
}

// Is this the last ghost chain that has the passed target_chain_id as it's target_chain_id?
//
// i.e. were we on the last ghost that has the same target_chain_id as this ghost has?
//
// for example, "D" -> "A" in A,B,C,D
//
bool
molecule_class_info_t::last_ghost_matching_target_chain_id_p(int i_ghost_chain_match,
                                                             const std::vector<drawn_ghost_molecule_display_t> &ncs_ghosts) const {

   bool is_last = 0;

   std::string match_target_chain_id = ncs_ghosts[i_ghost_chain_match].target_chain_id;
   int last_ghost_with_matching_target_chain_id = -1; // unset
   for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
      if (ncs_ghosts[ighost].target_chain_id == match_target_chain_id) {
         last_ghost_with_matching_target_chain_id = ighost;
      }
   }
   if (last_ghost_with_matching_target_chain_id != -1) {
      if (i_ghost_chain_match == last_ghost_with_matching_target_chain_id) {
         is_last = 1;
      }
   }
   return is_last;
}


// not const because we can update (fill) ncs_ghosts.
coot::ncs_differences_t
molecule_class_info_t::ncs_chain_differences(std::string master_chain_id,
                                             float main_chain_weight) {

   std::vector<coot::ncs_chain_difference_t> diffs;

   // Note to self: recall:
   //
   // class ncs_chain_difference_t {
   //    std::string peer_chain_id;
   //    std::vector<ncs_residue_info_t> residue_info;


   if (ncs_ghosts.size() > 0) {
      if (ncs_ghosts_have_rtops_flag == 0) {
         float homology_lev =0.7;
         fill_ghost_info(1, homology_lev); // fill the rtops and set the flag
      }
   }

   if (ncs_ghosts.size() > 0) {
      if (!ncs_ghosts_have_rtops_flag) {
         float homology_lev =0.7;
         fill_ghost_info(1, homology_lev); // fill the rtops and set the flag
      }

      for (unsigned int ighost = 0; ighost<ncs_ghosts.size(); ighost++) {
         int imod = 1;
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         mmdb::Chain *this_chain_p = 0;
         mmdb::Chain *master_chain_p = 0;
         for (int ichain=0; ichain<nchains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            if (std::string(chain_p->GetChainID()) == ncs_ghosts[ighost].chain_id) {
               this_chain_p = chain_p;
            }
            if (std::string(chain_p->GetChainID()) == ncs_ghosts[ighost].target_chain_id) {
               master_chain_p = chain_p;
            }
         }
         if (this_chain_p && master_chain_p) {

            int nres_this     = this_chain_p->GetNumberOfResidues();
            int nres_master = master_chain_p->GetNumberOfResidues();
            mmdb::PResidue this_residue_p;
            mmdb::PResidue master_residue_p;
            std::vector<coot::ncs_residue_info_t> residue_info;

            // return first == 0 if residues not found in chain
            std::pair<short int, int> mm_master =
               coot::util::min_resno_in_chain(master_chain_p);
            std::pair<short int, int> mm_this =
               coot::util::min_resno_in_chain(this_chain_p);

            if (mm_master.first && mm_this.first) {
               int resno_offset = mm_this.second - mm_master.second;

               for (int ires=0; ires<nres_this && ires<nres_master; ires++) {
                  master_residue_p = master_chain_p->GetResidue(ires);
                  this_residue_p = 0;
                  if ( ((ires-resno_offset)<nres_this) && ((ires-resno_offset)>=0)) { 
                     this_residue_p = this_chain_p->GetResidue(master_residue_p->GetSeqNum(),
                                                               master_residue_p->GetInsCode());
                  }
                  if (this_residue_p && master_residue_p) {
                     if (this_residue_p->GetSeqNum() == master_residue_p->GetSeqNum()) {
                        coot::ncs_residue_info_t ds = 
                           ncs_ghosts[ighost].get_differences(this_residue_p, master_residue_p,
                                                              main_chain_weight);
                        if (ds.filled) {
                           if (0)
                              std::cout << "     pushing back a residue_info with resno "
                                        << ds.resno << std::endl;
                           residue_info.push_back(ds);
                        }
                     }
                  }
               }
            }
            coot::ncs_chain_difference_t d(ncs_ghosts[ighost].chain_id, residue_info);
            if (residue_info.size() > 0)
               diffs.push_back(d);
         }
      }
   }
   return coot::ncs_differences_t(master_chain_id, diffs);
}


// Return a ncs_residue_info_t, note we cannot use that
// ncs_residue_info_t if filled is 0.
//
coot::ncs_residue_info_t
coot::ghost_molecule_display_t::get_differences(mmdb::Residue *this_residue_p,
                                                mmdb::Residue *master_residue_p,
                                                float main_chain_weight) const {
   // has access to rtop
   coot::ncs_residue_info_t r;

   if (std::string(this_residue_p->GetResName()) == std::string(master_residue_p->GetResName())) {
      std::vector<std::pair<int, int> > index_pairs =
         coot::util::pair_residue_atoms(this_residue_p, master_residue_p);
      mmdb::PPAtom residue_atoms_1 = NULL;
      mmdb::PPAtom residue_atoms_2 = NULL;
      int n_residue_atoms_1, n_residue_atoms_2;
        this_residue_p->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
      master_residue_p->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
      float n_weighted_atoms = 0.0;
      double sum_dist = 0.0;
      for (unsigned int i=0; i<index_pairs.size(); i++) {
         float atom_weight = 1.0;
         mmdb::Atom *at1 = residue_atoms_1[index_pairs[i].first];
         mmdb::Atom *at2 = residue_atoms_2[index_pairs[i].second];
         std::string resname_1 = at1->GetResName(); // PDB puts waters in protein chains, bleugh.
         std::string resname_2 = at2->GetResName();
         if ((resname_1 != "HOH") && (resname_2 != "HOH")) {
            if (!at1->isTer() && !at2->isTer()) {
               if (coot::is_main_chain_p(at1))
                  atom_weight = main_chain_weight;
               clipper::Coord_orth pt1(at1->x, at1->y, at1->z);
               clipper::Coord_orth pt2(at2->x, at2->y, at2->z);
               double len = clipper::Coord_orth::length(pt1.transform(rtop), pt2);
               sum_dist +=len;
               n_weighted_atoms += atom_weight;
            }
         }
      }
      if (n_weighted_atoms > 0.0) {
         r = coot::ncs_residue_info_t(this_residue_p->GetSeqNum(),
                                      this_residue_p->GetInsCode(),
                                      this_residue_p->index,
                                      master_residue_p->GetSeqNum(),
                                      master_residue_p->GetInsCode(),
                                      master_residue_p->index); // sets r.filled
         r.mean_diff = sum_dist/float(n_weighted_atoms);
         r.n_weighted_atoms = n_weighted_atoms;
      }
   } else {
      std::cout << "different residue types " << this_residue_p->GetSeqNum()
                << " " << this_residue_p->GetInsCode() << " "
                << this_residue_p->GetResName() << "   vs  "
                << master_residue_p->GetSeqNum() << " " << master_residue_p->GetInsCode()
                << " " << master_residue_p->GetResName() << std::endl;
   }
   return r;
}

std::vector<std::vector<std::string> >
molecule_class_info_t::ncs_ghost_chains() const {

   std::vector<std::vector<std::string> > v;
   std::vector<std::string> ghost_vec;
   if (ncs_ghosts.size() > 0) {
      std::vector<std::pair<std::string, std::vector<std::string> > > grouped_ghosts;
      for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {

         //           std::cout << "DEBUG:: doing ghost number " << ighost << " "
         //                     << ncs_ghosts[ighost].chain_id << " "
         //                     << ncs_ghosts[ighost].target_chain_id
         //                     << std::endl;

         std::string master_chain = ncs_ghosts[ighost].target_chain_id;

         // does the master_chain_id of this ghost already have a
         // group ghost for it? If yes, then simply add this ghost to
         // the grouped_ghosts for this master_chain_id, if not,
         // create a new grouped_ghosts item and add this chain id to
         // that item.
         //
         bool found = 0;
         for (unsigned int igroup=0; igroup<grouped_ghosts.size(); igroup++) {
            if (grouped_ghosts[igroup].first == master_chain) {
               grouped_ghosts[igroup].second.push_back(ncs_ghosts[ighost].chain_id);
               found = 1;
               break;
            }
         }
         if (! found) {
            std::vector<std::string> v;
            v.push_back(ncs_ghosts[ighost].chain_id);
            std::pair<std::string, std::vector<std::string> > group(master_chain, v);
            grouped_ghosts.push_back(group);
         }
      }

      if (false) {
         std::cout << "DEBUG:: There are " << grouped_ghosts.size() << " grouped ghosts"
                   << std::endl;
         for (unsigned int igroup=0; igroup<grouped_ghosts.size(); igroup++) {
            std::cout << "DEBUG:: Master chain id " <<  grouped_ghosts[igroup].first
                      << " with peers ";
            for (unsigned int ipeer=0; ipeer<grouped_ghosts[igroup].second.size(); ipeer++) {
               std::cout << grouped_ghosts[igroup].second[ipeer] << " ";
            }
            std::cout << std::endl;
         }
      }


      if (grouped_ghosts.size() > 0) {
         for (unsigned int igroup=0; igroup<grouped_ghosts.size(); igroup++) {
            std::vector<std::string> linear_ghost;
            linear_ghost.push_back(grouped_ghosts[igroup].first);
            for (unsigned int ipeer=0; ipeer<grouped_ghosts[igroup].second.size(); ipeer++) {
               linear_ghost.push_back(grouped_ghosts[igroup].second[ipeer]);
            }
            v.push_back(linear_ghost);
         }
      }
   }
   return v;
}

std::pair<bool, std::string>
molecule_class_info_t::first_ncs_master_chain_id() const {  // for ncs graphs use

   bool status = 0;
   std::string master_chain_id;
   for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
      master_chain_id = ncs_ghosts[ighost].target_chain_id;
      status = 1;
   }
   return std::pair<bool, std::string> (status, master_chain_id);
}

std::vector<std::string>
molecule_class_info_t::ncs_master_chains() const {

   std::vector<std::string> s;
   for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
      std::string master_chain_id = ncs_ghosts[ighost].target_chain_id;
      if (std::find(s.begin(), s.end(), master_chain_id) ==
          s.end())
         s.push_back(master_chain_id);
   }
   return s;
}

void
molecule_class_info_t::add_strict_ncs_from_mtrix_from_file(const std::string &file_name) {

   std::vector<clipper::RTop_orth> mv = coot::mtrix_info(file_name);
   for (unsigned int i=0; i<mv.size(); i++) {
      const clipper::RTop_orth &rt = mv[i];

      // rt.rot()(0,0), rt.rot()(0,1), rt.rot()(0,2),
      // rt.rot()(1,0), rt.rot()(1,1), rt.rot()(1,2),
      // rt.rot()(2,0), rt.rot()(2,1), rt.rot()(2,2),
      // rt.trn()[0],   rt.trn()[1],   rt.trn()[2]);

      coot::coot_mat44 cm44;

      cm44.m[0].v4[0] = rt.rot()(0,0);
      cm44.m[0].v4[1] = rt.rot()(0,1);
      cm44.m[0].v4[2] = rt.rot()(0,2);
      cm44.m[1].v4[0] = rt.rot()(1,0);
      cm44.m[1].v4[1] = rt.rot()(1,1);
      cm44.m[1].v4[2] = rt.rot()(1,2);
      cm44.m[2].v4[0] = rt.rot()(2,0);
      cm44.m[2].v4[1] = rt.rot()(2,1);
      cm44.m[2].v4[2] = rt.rot()(2,2);
      // translation
      cm44.m[0].v4[3] = rt.trn()[0];
      cm44.m[1].v4[3] = rt.trn()[1];
      cm44.m[2].v4[3] = rt.trn()[2];
      // sensibles
      cm44.m[3].v4[0] = 0.0;
      cm44.m[3].v4[1] = 0.0;
      cm44.m[3].v4[2] = 0.0;
      cm44.m[3].v4[3] = 1.0;

      add_strict_ncs_matrix( "A", "A", cm44);
   }
}


void
molecule_class_info_t::add_strict_ncs_from_mtrix_from_self_file() {

   if (has_model()) {
      std::string file_name = name_;
      if (coot::file_exists(file_name)) {
         add_strict_ncs_from_mtrix_from_file(file_name);
      } else {
         std::cout << "WARNING:: in add_strict_ncs_from_mtrix_from_self_file() "
                   << "file " << file_name << " not found" << std::endl;
      }
   }
}


// if we are are the centre of a given chain_id, how big a radius
// do we need to encompass all atoms of that chain?
//
std::pair<clipper::Coord_orth, double>
molecule_class_info_t::chain_centre_and_radius(const std::string &chain_id) const {

   clipper::Coord_orth pos(0,0,0);
   double radius = -1;

   mmdb::Manager *mol = atom_sel.mol;

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         std::string chain_id_this = chain_p->GetChainID();
         if (chain_id_this == chain_id) {
            std::vector<clipper::Coord_orth> chain_atoms_positions;
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               int n_atoms = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms; iat++) {
                  mmdb::Atom *at = residue_p->GetAtom(iat);
                  chain_atoms_positions.push_back(coot::co(at));
               }
            }

            unsigned int n_atoms_in_chain = chain_atoms_positions.size();
            if (n_atoms_in_chain > 0) {
               double sum_x = 0;
               double sum_y = 0;
               double sum_z = 0;
               for (unsigned int i=0; i<chain_atoms_positions.size(); i++) {
                  sum_x += chain_atoms_positions[i].x();
                  sum_y += chain_atoms_positions[i].y();
                  sum_z += chain_atoms_positions[i].z();
               }

               pos = clipper::Coord_orth(sum_x/double(n_atoms_in_chain),
                                         sum_y/double(n_atoms_in_chain),
                                         sum_z/double(n_atoms_in_chain));

               double max_d_sqrd = 0.0;

               for (unsigned int i=0; i<chain_atoms_positions.size(); i++) {
                  double d = (chain_atoms_positions[i] - pos).lengthsq();
                  if (d > max_d_sqrd) {
                     max_d_sqrd = d;
                  }
               }
               radius = sqrt(max_d_sqrd);
            }
         }
      }
   }

   return std::pair<clipper::Coord_orth, double> (pos, radius);
}



// process REMARK 350s
void
molecule_class_info_t::add_molecular_symmetry_matrices() {

   if (atom_sel.mol) {

      std::vector<coot::coot_mat44> biomt_matrices;
      std::vector<std::string> biomt_chain_ids;

      mmdb::TitleContainer *tc_p = atom_sel.mol->GetRemarks();
      int l = tc_p->Length();
      std::map<std::pair<int, int>, quad_d_t> biomts;
      for (int i=0; i<l; i++) {
         mmdb::Remark *cr = static_cast<mmdb::Remark *> (tc_p->GetContainerClass(i));
         if (cr->remarkNum == 350) {
            std::string rm = cr->remark;
            // std::cout << l << " " << rm << std::endl;
            std::vector<std::string> parts = coot::util::split_string_no_blanks(rm);

            // which chains?
            if (parts.size() > 5) {
               if (parts[0] == "APPLY") {
                  if (parts[1] == "THE") {
                     if (parts[2] == "FOLLOWING") {
                        if (parts[3] == "TO") {
                           if (parts[4] == "CHAINS:") {
                              unsigned int n = parts.size();
                              for (unsigned int ii=5; ii<n; ii++) {
                                 unsigned int ll = parts[ii].length();
                                 if (l > 1) {
                                    std::string chain_id = parts[ii].substr(0,ll-1);
                                    biomt_chain_ids.push_back(chain_id);
                                 } else {
                                    if (ll == 1) {
                                       // no comma
                                       std::string chain_id = parts[ii];
                                       biomt_chain_ids.push_back(chain_id);
                                    }
                                 } 
                              }
                           }
                        }
                     }
                  }
               }
            }

            if (parts.size() == 6) {

               if (parts[0].substr(0,5) == "BIOMT") {
               
                  // for (unsigned int i=0; i<parts.size(); i++)
                  //   std::cout << i << " " << parts[i] << std::endl;
               
                  try {
                     int matrix_id = coot::util::string_to_int(parts[1]);
                     int ll = parts[0].length();
                     char c = parts[0][ll-1];
                     int biomt_idx = c - 48;

                     double x = coot::util::string_to_double(parts[2]);
                     double y = coot::util::string_to_double(parts[3]);
                     double z = coot::util::string_to_double(parts[4]);
                     double t = coot::util::string_to_double(parts[5]);
                     if (false)
                        std::cout << "   " << matrix_id << " " << biomt_idx << " "
                                  << x << " " << y << " " << z << " " << t 
                                  << std::endl;
                     std::pair<int, int> key(matrix_id, biomt_idx);
                     biomts[key] = quad_d_t(x,y,z,t);
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << "WARNING:: " << rte.what() << std::endl;
                  }
               }
            }
         }
      }

      // read lines. Now construct matrices.
      //
      // first determine max matrix_id (e.g. 60)
      //
      int matrix_id_max = 0;
      std::map<std::pair<int, int>, quad_d_t>::const_iterator it;
      for (it=biomts.begin(); it!=biomts.end(); ++it) {
         if (it->first.first > matrix_id_max)
            matrix_id_max = it->first.first;
      }

      if (matrix_id_max > 0) {
         for (int idx=0; idx<=matrix_id_max; idx++) {
            std::pair<int, int> key_1(idx, 1);
            std::pair<int, int> key_2(idx, 2);
            std::pair<int, int> key_3(idx, 3);

            std::map<std::pair<int, int>, quad_d_t>::const_iterator it_1 = biomts.find(key_1);
            std::map<std::pair<int, int>, quad_d_t>::const_iterator it_2 = biomts.find(key_2);
            std::map<std::pair<int, int>, quad_d_t>::const_iterator it_3 = biomts.find(key_3);

            if (it_1 != biomts.end()) {
               if (it_2 != biomts.end()) {
                  if (it_3 != biomts.end()) {

                     coot::coot_mat44 m;
                     m.m[0].v4[0] = it_1->second.x;
                     m.m[0].v4[1] = it_1->second.y;
                     m.m[0].v4[2] = it_1->second.z;
                     m.m[0].v4[3] = it_1->second.t;
                     m.m[1].v4[0] = it_2->second.x;
                     m.m[1].v4[1] = it_2->second.y;
                     m.m[1].v4[2] = it_2->second.z;
                     m.m[1].v4[3] = it_2->second.t;
                     m.m[2].v4[0] = it_3->second.x;
                     m.m[2].v4[1] = it_3->second.y;
                     m.m[2].v4[2] = it_3->second.z;
                     m.m[2].v4[3] = it_3->second.t;
                     m.m[3].v4[0] = 0;
                     m.m[3].v4[1] = 0;
                     m.m[3].v4[2] = 0;
                     m.m[3].v4[3] = 0;

                     biomt_matrices.push_back(m);
                  }
               }
            }
         }
      }

      if (false)
         std::cout << "in add_molecular_symmetry_matrices() made "
                   << biomt_matrices.size() << " biomt matrices" << std::endl;

      for (unsigned int jj=0; jj<biomt_chain_ids.size(); jj++) {
         for (unsigned int ii=0; ii<biomt_matrices.size(); ii++) {
            if (! biomt_matrices[ii].is_close_to_unit_matrix())
               add_strict_ncs_matrix(biomt_chain_ids[jj],
                                     biomt_chain_ids[jj],
                                     biomt_matrices[ii]);
         }
      }
   }
}
