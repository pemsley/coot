/* pli/flev.cc
 * 
 * Copyright 2012 by The University of Oxford
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "geometry/residue-and-atom-specs.hh"
#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <cstdint> // fix std::uint64_t problems
#include "lidia-core/svg-container.hh"
#include "lidia-core/rdkit-interface.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "flev.hh"
#include "flev-annotations.hh"
#include "dots-representation-info.hh"
#include "protein-ligand-interactions.hh"
#include "atom-ring-centre-info.hh"
#include <string>

std::string make_arrow(const lig_build::pos_t &A, const lig_build::pos_t &B,
                       const std::string &stroke_colour, bool start_arrow, bool end_arrow,
                       const lig_build::pos_t &pos,
                       const lig_build::pos_t &lig_atom_pos);

std::map<std::string, std::string>
pli::make_flat_ligand_name_map(mmdb::Residue *flat_res) {

   double bond_to_H_dist = 1.1;

   double b2Hd2 = bond_to_H_dist * bond_to_H_dist;
   std::map<std::string, std::string> map;
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   flat_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at_i = residue_atoms[iat];
      std::string ele_i = at_i->element;
      clipper::Coord_orth pt_i(at_i->x, at_i->y, at_i->z);
      if (ele_i == " H") {
         for (int jat=0; jat<n_residue_atoms; jat++) {
            if (iat != jat) {
               mmdb::Atom *at_j = residue_atoms[jat];
               std::string ele_j = at_j->element;
               if (ele_j != " H") {
                  clipper::Coord_orth pt_j(at_j->x, at_j->y, at_j->z);
                  if ((pt_i - pt_j).lengthsq() < b2Hd2) {
                     map[at_j->name] = at_i->name;
                     break;
                  }
               }
            }
         }
      }
   }

//    std::map<std::string, std::string>::const_iterator it;
//    std::cout << "=== name map: === " << std::endl;
//    for (it=map.begin(); it!=map.end(); it++) {
//       std::cout << "  :" << it->first << ": :" << it->second << ":" << std::endl;
//    }
//    std::cout << "=== === " << std::endl;

   return map;
}

std::vector<pli::fle_residues_helper_t>
pli::get_flev_residue_centres(mmdb::Residue *residue_ligand_3d,
                               mmdb::Manager *mol_containing_residue_ligand,
                               std::vector<mmdb::Residue *> residues,
                               mmdb::Manager *flat_mol) {

   std::vector<fle_residues_helper_t> centres;

   if (flat_mol) {

      // get the lsq matrix that maps the ligand in 3D onto the flat ligand
      int res_no = residue_ligand_3d->GetSeqNum();
      std::string chain_id = residue_ligand_3d->GetChainID();
      int every_nth = 1;
      std::vector<coot::lsq_range_match_info_t> matches;
      coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id, coot::lsq_t::ALL);
       matches.push_back(match);
       std::pair<short int, clipper::RTop_orth> lsq_mat =
          coot::util::get_lsq_matrix(flat_mol, mol_containing_residue_ligand, matches, every_nth, false);
      // Now make the residues

      // std::vector<coot::fle_residues_helper_t> centres(residues.size());
      centres.resize(residues.size());
      for (unsigned int ires=0; ires<residues.size(); ires++) {
         mmdb::Residue *res_copy = coot::util::deep_copy_this_residue(residues[ires]);
         std::string res_name = residues[ires]->GetResName();
         std::pair<bool, clipper::Coord_orth> absolute_centre =
            coot::util::get_residue_centre(res_copy);
         if (absolute_centre.first) {
            coot::util::transform_atoms(res_copy, lsq_mat.second);
            std::pair<bool, clipper::Coord_orth> c =
               coot::util::get_residue_centre(res_copy);
            if (c.first) {
               fle_residues_helper_t fle_centre(c.second, coot::residue_spec_t(residues[ires]), res_name);

               // Setting the interaction position to the residue
               // centre is a hack.  What we need to be is the middle
               // of the hydrogen bond or the middle of the pi-pi
               // stacking (for instance).  To do that we need to know
               // what positions of the interacting in the ligand (and
               // of the residue).  Tricky from here?
               //
               fle_centre.set_interaction_position(absolute_centre.second);
               centres[ires] = fle_centre;
            } else {
               std::cout << "WARNING:: get_flev_residue_centres() failed to get residue centre for "
                         << coot::residue_spec_t(res_copy) << std::endl;
            }
         } else {
            std::cout << "WARNING:: get_flev_residue_centres() failed to get residue centre for "
                      << coot::residue_spec_t(res_copy) << std::endl;
         }
         delete res_copy;
      }
   }
   if (false)
      for (unsigned int ic=0; ic<centres.size(); ic++) {
         std::cout << "centre " << ic << " has TR_centre "
                   << centres[ic].transformed_relative_centre.format() << std::endl;
      }
   return centres;
}

std::vector<pli::solvent_accessible_atom_t>
flev_t::convert(const std::vector<std::pair<coot::atom_spec_t, float> > &s_a_v,
                const pli::flev_attached_hydrogens_t &ah) const {

   std::vector<pli::solvent_accessible_atom_t> r(s_a_v.size());

   for (unsigned int i=0; i<s_a_v.size(); i++) {
      std::string name = s_a_v[i].first.atom_name;
      r[i].atom_name = name;
      r[i].solvent_accessibility = s_a_v[i].second;
      // std::cout << "convert(): " << i << " " << name << " " << s_a_v[i].second << std::endl;

      std::map<std::string, std::vector<coot::bash_distance_t> >::const_iterator it =
         ah.atom_bashes.find(name);
      if (it != ah.atom_bashes.end())
         r[i].bash_distances = it->second;
   }
   return r;
}

void
flev_t::render() {

   std::cout << "render" << std::endl;

}

// return 1 on solution having problems, 0 for no problems, also
// return a list of the residue circles with problems.
//
std::pair<bool, std::vector<int> >
flev_t::solution_has_problems_p() const {

   std::vector<int> v;
   bool status = 0;

   // scale factor of 1.2 is too small, 2.0 seems good for H102 in 2aei
   // double crit_dist_2 = 2.0 * SINGLE_BOND_CANVAS_LENGTH * SINGLE_BOND_CANVAS_LENGTH;
   double crit_dist = 2.5;
   double crit_dist_sqrd = crit_dist * crit_dist;
   double crit_dist_wide_2 = 4 * crit_dist_sqrd;

   // a problem can be 2 atoms < crit_dist_2 or
   //                  5 atoms < crit_dist_wide_2
   //
   // we need the wider check to find atoms that are enclosed by a
   // system of 10 or so atoms.
   //

   for (unsigned int i=0; i<residue_circles.size(); i++) {
      int n_close = 0;
      int n_close_wide = 0;
      for (unsigned int j=0; j<mol.atoms.size(); j++) {
	 double d2 = (residue_circles[i].pos-mol.atoms[j].atom_position).lengthsq();
         std::string mark = "";
         if (d2 < crit_dist_sqrd) mark = " ***";
         if (false)
            std::cout << "in solution_has_problems_p(): comparing " << std::sqrt(d2) << " " << crit_dist
                      << " residue: " << residue_circles[i].residue_label
                      << " atom: " << mol.atoms[j].atom_name << " " << mark
                      << std::endl;
	 if (d2 < crit_dist_sqrd) {
	    n_close++;
	    if (n_close > 1) {
	       v.push_back(i);
	       status = 1;
	       break;
	    }
	 }

	 if (d2 < crit_dist_wide_2) {
	    n_close_wide++;
	    if (n_close_wide > 5) {
	       v.push_back(i);
	       status = 1;
	       break;
	    }
	 }
      }
   }
   return std::pair<bool, std::vector<int> > (status, v);
}

// check that the residues with bonds to the ligand have sane distances to the
// bonding ligand atom
std::pair<bool, std::vector<int> >
flev_t::solution_has_bonded_problems_p() const {

   bool status = false;
   std::vector<int> v;
   for (unsigned int i=0; i<residue_circles.size(); i++) {
      const auto &rc = residue_circles[i];
      if (rc.is_a_primary_residue()) {
         for (unsigned int ib=0; ib<rc.bonds_to_ligand.size(); ib++) {
            const auto &btl = rc.bonds_to_ligand[ib];
            const std::string &ligand_atom_name = btl.ligand_atom_name;
            std::pair<bool, lig_build::atom_t> ligand_atom_pair = mol.get_atom_by_name(ligand_atom_name);
            if (ligand_atom_pair.first) {
               const auto &ligand_atom = ligand_atom_pair.second;
               lig_build::pos_t ligand_atom_pos = ligand_atom.atom_position;
               double d = (rc.pos - ligand_atom_pos).length();
               if (d > (btl.bond_length + 1.9)) {
                  std::cout << "target distance for atom " << ligand_atom_name << " "
                            << btl.bond_length << " vs actual: " << d << std::endl;
                  status = true;
                  v.push_back(i);
               }
            }
         }
      }
   }
   return std::pair<bool, std::vector<int> >(status, v);
}

// fiddle with the position of problematic residues that are bonded to the ligand
//
void
flev_t::reposition_bonded_problematics_and_reoptimise(const std::vector<int> &problematics,
                                                      const std::vector<int> &primary_indices) {

   std::cout << "in reposition_bonded_problematics_and_reoptimise() we have "
             << problematics.size() << " problematics " << std::endl;
   for (unsigned int ip=0; ip<problematics.size(); ip++) {
      const auto &problematic = residue_circles[problematics[ip]];
      const lig_build::pos_t &res_pos = problematic.pos;
      std::vector<std::pair<lig_build::pos_t, double> > aps =
         problematic.get_attachment_points(mol);
      if (! aps.empty()) {
         for (unsigned int i=0; i<aps.size(); i++) {
            const auto &attachment_point_pair = aps[i];

         }
      }
   }
}


// fiddle with the position of some residues in residue_circles
//
void
flev_t::reposition_problematics_and_reoptimise(const std::vector<int> &problematics,
                                               const std::vector<int> &primary_indices) {

   // is this function used!?

   std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair = mol.ligand_extents();
   ligand_grid grid(l_e_pair.first, l_e_pair.second);
   grid.fill(mol);

   for (unsigned int ip=0; ip<problematics.size(); ip++) {

      // A trapped residue is now treated as a primary (bonding)
      // residue with a distance to the "bad" solution of 3.5A (this
      // distance does not matter) - for initial positioning.  It is
      // not added to the primaries, so this bond length is not
      // optimised.

      std::vector<std::pair<lig_build::pos_t, double> > attachment_points;
      lig_build::pos_t attach_pos = residue_circles[problematics[ip]].pos;
      attach_pos += lig_build::pos_t (5.0 * double(rand())/double(RAND_MAX),
                                      5.0 * double(rand())/double(RAND_MAX));

      std::pair<lig_build::pos_t, double> p(attach_pos, 3.5);

      attachment_points.push_back(p);
      initial_primary_residue_circles_layout(grid, problematics[ip], attachment_points);
   }

   // set the initial positions for optimisation, now that we have
   // fiddled with residue_circles
   std::vector<residue_circle_t> current_circles = residue_circles;

   for (int iround=0; iround<120; iround++) {
      std::pair<int, std::vector<residue_circle_t> > new_c =
	 optimise_residue_circle_positions(residue_circles, current_circles, primary_indices);
      current_circles = new_c.second;
      if (new_c.first == GSL_ENOPROG)
	 break;
      if (new_c.first == GSL_SUCCESS)
	 break;
   }
   // now transfer results from the updating/tmp current_circles to the class data item:
   residue_circles = current_circles;
}


bool
flev_t::annotate(const std::vector<std::pair<coot::atom_spec_t, float> > &s_a_v,
                 const std::vector<pli::fle_residues_helper_t> &centres,
                 const std::vector<int> &additional_representation_handles_in,
                 const std::vector<pli::fle_ligand_bond_t> &bonds_to_ligand,
                 const std::vector<pli::solvent_exposure_difference_helper_t> &sed,
                 const pli::flev_attached_hydrogens_t &ah,
                 const pli::pi_stacking_container_t &pi_stack_info,
                 const coot::dictionary_residue_restraints_t &restraints) {

   auto input_coords_to_canvas_coords = [] (const clipper::Coord_orth &pos) {
      // apply scale and offset also
      lig_build::pos_t p(pos.x(), pos.y());
      return p;
   };

   auto map_solvent_accessibilities_to_atoms = [] (lig_build::molecule_t<svg_atom_t, svg_bond_t> *new_mol_p,
                                                   const std::vector<pli::solvent_accessible_atom_t> &solvent_accessible_atoms) {
      for (unsigned int i=0; i<new_mol_p->atoms.size(); i++) {
         svg_atom_t &atom = new_mol_p->atoms[i];
         for (unsigned int j=0; j<solvent_accessible_atoms.size(); j++) {
            if (atom.get_atom_name() == solvent_accessible_atoms[j].atom_name) {
               atom.add_solvent_accessibility(solvent_accessible_atoms[j].solvent_accessibility);
               atom.bash_distances = solvent_accessible_atoms[j].bash_distances;
               break;
            }
         }
      }
   };

   if (true)
      std::cout << "--------------------------- debug in ::annotate() with pi_stack_info.stackings.size "
                << pi_stack_info.stackings.size() << std::endl;

   if (true) {
      std::cout << "debug::flev_t::annotate(): " << bonds_to_ligand.size() << " bonds to ligand" << std::endl;
      for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) {
         std::cout << "  =============== flev_t::annotate() bond to ligand " << ib << " "
                   << bonds_to_ligand[ib].ligand_atom_spec.atom_name
                   << " by residue " << bonds_to_ligand[ib].res_spec.chain_id << " "
                   << bonds_to_ligand[ib].res_spec.res_no << " type: "
                   << bonds_to_ligand[ib].bond_type
                   << std::endl;
      }
   }

   if (false) {
      std::cout << "--------------------------------------------------------------" << std::endl;
      std::cout << "======== flev_t::annotate() here are bash distances for atoms:" << std::endl;
      std::map<std::string, std::vector<coot::bash_distance_t> >::const_iterator it;
      for (it=ah.atom_bashes.begin(); it!=ah.atom_bashes.end(); it++) {
         std::cout << it->first << " " << it->second.size() << " bashes " << std::endl;
         for (unsigned int i=0; i<it->second.size(); i++) {
            std::cout << "   " << it->second[i] << std::endl;
         }
      }
      std::cout << "--------------------------------------------------------------" << std::endl;
   }

   bool r = false;
   std::vector<pli::solvent_accessible_atom_t> solvent_accessible_atoms = convert(s_a_v, ah);
   std::cout << "annotate(): solvent_accessible_atoms size " << solvent_accessible_atoms.size() << std::endl;

   // modify the atoms of mol
   map_solvent_accessibilities_to_atoms(&mol, solvent_accessible_atoms);

   // fill class data item std::vector<residue_circle_t> residue_circles;
   //
   residue_circles.clear();
   for (unsigned int i=0; i<centres.size(); i++) {

      if (false)
         std::cout << "debug:: in flev_t::annotate() handling circle " << i << " of "
                   << centres.size() << std::endl;

      std::string label = centres[i].spec.chain_id;
      label += std::to_string(centres[i].spec.res_no);
      label += centres[i].spec.ins_code;

      clipper::Coord_orth cp(centres[i].transformed_relative_centre.x(),
                             centres[i].transformed_relative_centre.y(),
                             centres[i].transformed_relative_centre.z());

      residue_circle_t circle(cp,
                              centres[i].interaction_position,
                              centres[i].spec,
                              centres[i].residue_name,
                              label);

      // let's have a different converter
      // lig_build::pos_t pos = mol.input_coords_to_canvas_coords(cp); // 20240601-PE
      lig_build::pos_t pos = input_coords_to_canvas_coords(cp);
      circle.set_canvas_pos(pos);
      if (false)
         std::cout << "debug:: in flev_t::annotate() residue-circle "
                   << cp.x() << " " << cp.y() << " " << cp.z()
                   << " canvas coord " << pos.x << " " << pos.y << std::endl;

      if (centres[i].residue_name == "HOH") {
         for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) {
            if (bonds_to_ligand[ib].res_spec == centres[i].spec) {
               circle.set_water_dist_to_protein(bonds_to_ligand[ib].water_protein_length);
            }
         }
      }

      // now add any H-bonds to the ligand:
      //
      for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) {
         if (bonds_to_ligand[ib].res_spec == centres[i].spec) {
            double bond_l = bonds_to_ligand[ib].bond_length;
            std::string ligand_atom_name = bonds_to_ligand[ib].ligand_atom_spec.atom_name;
            bond_to_ligand_t btl(ligand_atom_name, bond_l);
            btl.bond_type = bonds_to_ligand[ib].bond_type;
            circle.add_bond_to_ligand(btl);
         }
      }

      // Now add any pi-stacking interactions to/from the ligand:
      //
      for (unsigned int istack=0; istack<pi_stack_info.stackings.size(); istack++) {
         coot::residue_spec_t spec(pi_stack_info.stackings[istack].res);
         if (spec == centres[i].spec) {
            std::cout << "------------------- pi_stacking for residue " << spec << std::endl;
            if (pi_stack_info.stackings[istack].type == pli::pi_stacking_instance_t::PI_PI_STACKING) {
               std::vector<std::string> lra = pi_stack_info.stackings[istack].ligand_ring_atom_names;
               circle.set_stacking("pi-pi", lra, "");
            }
            if (pi_stack_info.stackings[istack].type == pli::pi_stacking_instance_t::PI_CATION_STACKING) {
               std::vector<std::string> lra = pi_stack_info.stackings[istack].ligand_ring_atom_names;
               circle.set_stacking("pi-cation", lra, "");
            }
            // cation in ligand, PHE (say)
            if (pi_stack_info.stackings[istack].type == pli::pi_stacking_instance_t::CATION_PI_STACKING) {
               std::vector<std::string> lra_null;
               circle.set_stacking("cation-pi", lra_null, pi_stack_info.stackings[istack].ligand_cationic_atom_name);
            }
         }
      }

      // solvent exposure difference annotation.
      //
      // std::vector<coot::solvent_exposure_difference_helper_t> &sed,
      //
      for (unsigned int ised=0; ised<sed.size(); ised++) {
         if (sed[ised].res_spec == centres[i].spec) {
            circle.set_solvent_exposure_diff(sed[ised].exposure_fraction_holo,
                                             sed[ised].exposure_fraction_apo);
         }
      }
      residue_circles.push_back(circle);
   }

   // additional_representation_handles transfer
   additional_representation_handles = additional_representation_handles_in;


   // ligand ring centres (some of which may be aromatic and are in pi_stack_info too).
   //
   // However, for the substitution contour and the initial layout
   // ring avoidance, we blob in a circle of radius 1/(2*sin(180/n_ring_atoms)) bond lengths.

   // sets the object variable
   ring_atoms_list = restraints.get_ligand_ring_list();

   // just checking that it was passed correctly -
   //
   if (false) {
      std::cout << "debug:: in flev_t::annotate() ------------------- ring list ------------" << std::endl;
      for (unsigned int i=0; i<ring_atoms_list.size(); i++) {
         std::cout << "ring list " << i << "   ";
         for (unsigned int j=0; j<ring_atoms_list[i].size(); j++) {
         std::cout << ring_atoms_list[i][j] << "  ";
         }
         std::cout << std::endl;
      }
   }

   if (false) {
      std::cout << "debug:: in flev_t::annotate() ------------------- residue circles ------------" << std::endl;
      for (unsigned int i=0; i<residue_circles.size(); i++) {
         const auto &rc = residue_circles[i];
         std::cout << "   " << std::setw(2) << i << " : "
                 << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.x << " "
                 << std::setw(10) << std::setprecision(5) << std::right << std::fixed << rc.pos.y << std::endl;
      }
   }

   refine_residue_circle_positions();

   std::pair<bool, std::vector<int> > bonded_problem_status = solution_has_bonded_problems_p();
   if (false)
      std::cout << "::::::::: bonded problem status: " << bonded_problem_status.first << std::endl;
   for (unsigned int ip=0; ip<bonded_problem_status.second.size(); ip++) {
      std::cout << ":::::::::::: "
                << residue_circles[bonded_problem_status.second[ip]].residue_label << " "
                << residue_circles[bonded_problem_status.second[ip]].residue_type << " "
                << residue_circles[bonded_problem_status.second[ip]].pos
                << std::endl;
   }
   if (! bonded_problem_status.second.empty()) {
      // fiddle with residue_circles and reoptimise.
      std::vector<int> primary_indices = get_primary_indices();
      // ordinary residues
      reposition_bonded_problematics_and_reoptimise(bonded_problem_status.second, primary_indices);
   }

   // Now test the ordinary residues (not bonded to the ligand)
   //
   // has the current solution problems due to residues too close to the ligand?

   std::pair<bool, std::vector<int> > problem_status = solution_has_problems_p();

   if (false) {
      std::cout << "::::::::: problem status: " << problem_status.first << std::endl;
      for (unsigned int ip=0; ip<problem_status.second.size(); ip++) {
         std::cout << ":::::::::::: "
                   << residue_circles[problem_status.second[ip]].residue_label << " "
                   << residue_circles[problem_status.second[ip]].residue_type << " "
                   << residue_circles[problem_status.second[ip]].pos
                   << std::endl;
      }
   }
   if (! problem_status.second.empty()) {
      // fiddle with residue_circles and reoptimise.
      std::vector<int> primary_indices = get_primary_indices();
      // ordinary residues
      reposition_problematics_and_reoptimise(problem_status.second, primary_indices);
   }

   return r;
}


void
flev_t::write_png(const std::string &file_name) const {

   std::cout << "write png file " << file_name << std::endl;

}

void
flev_t::write_svg(const std::string &file_name) const {

   std::cout << "write svg file " << file_name << std::endl;
}

svg_container_t
pli::make_key(const lig_build::pos_t &key_pos) {

   auto add_text = [] (const lig_build::pos_t &pos, const std::string &tt) {
      std::string text_colour = "#111111";
      std::string t("   <text ");
      t += std::string("fill=\"");
      t += text_colour;
      t += std::string("\"");
      t += std::string(" x=\"");
      t += std::to_string(pos.x);
      t += std::string("\"");
      t += std::string(" y=\"");
      t += std::to_string(pos.y);
      t += std::string("\"");
      t += std::string(" text-anchor=\"start\"");
      t += std::string(" font-family=\"Helvetica, sans-serif\" font-size=\"0.05em\">");
      t += tt;
      t += std::string("</text>\n");
      return t;
   };

   auto draw_dashed_ellipse = [] (const lig_build::pos_t &pos) {

      double radius_x = 1.14;
      double radius_y = 0.8;
      double stroke_width = 0.1;
      std::string stroke_color = "#808080";
      std::string circle_string = std::string("   ") + "<ellipse cx=\"" + std::to_string(pos.x) + std::string("\" cy=\"") +
           std::to_string(pos.y) +
           std::string("\" rx=\"") +
           std::to_string(radius_x) +
           std::string("\" ry=\"") +
           std::to_string(radius_y) +
           std::string("\"");
      circle_string += " fill=\"none\" ";
      circle_string += "style=\"stroke-linecap:round; stroke-dasharray:0.1,0.2;\" ";
      circle_string +=
           std::string(" stroke=\"") + stroke_color + std::string("\"") +
           std::string(" stroke-width=\"") + std::to_string(stroke_width) + std::string("\"") +
           std::string("/>\n");
      return circle_string;
   };

   // first column 

   svg_container_t svgc;
   svgc.add_comment("Key");
   std::string green = "darkgreen";
   std::string lime = "#888820";
   std::string metal = "#990099";
   lig_build::pos_t A = key_pos;
   lig_build::pos_t top_left = key_pos;
   lig_build::pos_t B = A + lig_build::pos_t(4.4, 0);
   lig_build::pos_t x_bit(0.7, 0);
   lig_build::pos_t x_step_for_text(7.7, 0);
   double y_line_step = 1.7;
   std::string a1 = make_arrow(A, B, "blue", false, true, A, B + x_bit);
   svgc.add(a1);
   A += lig_build::pos_t(0, -y_line_step);
   B += lig_build::pos_t(0, -y_line_step);
   std::string a2 = make_arrow(A, B, "blue", true, false, A - x_bit * 1.9, B);
   svgc.add(a2);
   A += lig_build::pos_t(0, -y_line_step);
   B += lig_build::pos_t(0, -y_line_step);
   std::string a3 = make_arrow(A, B, green, false, true, A, B + x_bit);
   svgc.add(a3);
   A += lig_build::pos_t(0, -y_line_step);
   B += lig_build::pos_t(0, -y_line_step);
   std::string a4 = make_arrow(A, B, green, true, false, A - x_bit * 1.9, B);
   svgc.add(a4);
   A += lig_build::pos_t(0, -y_line_step);
   B += lig_build::pos_t(0, -y_line_step);
   std::string a5 = make_arrow(A, B, lime, false, false, A - x_bit * 1.9, B);
   svgc.add(a5);
   A += lig_build::pos_t(0, -y_line_step);
   B += lig_build::pos_t(0, -y_line_step);
   std::string a6 = make_arrow(A, B, metal, false, false, A - x_bit * 1.9, B);
   svgc.add(a6);
   double top_left_x = key_pos.x;
   double top_left_y = key_pos.y;
   lig_build::pos_t text_pos_1(top_left_x + 5.5, -top_left_y + 0.2);
   lig_build::pos_t text_pos_2(top_left_x + 5.5, -top_left_y + 0.2 + 1.0 * y_line_step);
   lig_build::pos_t text_pos_3(top_left_x + 5.5, -top_left_y + 0.2 + 2.0 * y_line_step);
   lig_build::pos_t text_pos_4(top_left_x + 5.5, -top_left_y + 0.2 + 3.0 * y_line_step);
   lig_build::pos_t text_pos_5(top_left_x + 5.5, -top_left_y + 0.2 + 4.0 * y_line_step);
   lig_build::pos_t text_pos_6(top_left_x + 5.5, -top_left_y + 0.2 + 5.0 * y_line_step);
   std::string t1 = add_text(text_pos_1, "Main-chain acceptor");
   std::string t2 = add_text(text_pos_2, "Main-chain donor");
   std::string t3 = add_text(text_pos_3, "Side-chain acceptor");
   std::string t4 = add_text(text_pos_4, "Side-chain donor");
   std::string t5 = add_text(text_pos_5, "H-bond to Water");
   std::string t6 = add_text(text_pos_6, "Metal bond");
   svgc.add(t1);
   svgc.add(t2);
   svgc.add(t3);
   svgc.add(t4);
   svgc.add(t5);
   svgc.add(t6);

   // second column

   std::string purple = "#eeccee";
   std::string red    = "#cc0000";
   std::string blue   = "#0000cc";
   std::string metalic_grey = "#d9d9d9";

   std::string grease = "#ccffbb";
   std::string polar  = purple;
   std::string acidic = red;
   std::string basic  = blue;
   lig_build::pos_t tlc2 = top_left.invert_y() + lig_build::pos_t(15.0, 0.0);
   lig_build::pos_t grease_pos = tlc2;
   std::string circle_g = make_circle(grease_pos, 0.75, 0.1, grease, "#111111");
   svgc.add(circle_g);
   lig_build::pos_t polar_pos = tlc2 + lig_build::pos_t(0, y_line_step);
   std::string circle_p = make_circle(polar_pos, 0.75, 0.1, polar, "#111111");
   svgc.add(circle_p);
   lig_build::pos_t acidic_pos = tlc2 + lig_build::pos_t(0, 2.0 * y_line_step);
   std::string circle_a = make_circle(acidic_pos, 0.75, 0.15, polar, red);
   svgc.add(circle_a);
   lig_build::pos_t basic_pos = tlc2 + lig_build::pos_t(0, 3.0 * y_line_step);
   std::string circle_b = make_circle(basic_pos, 0.75, 0.15, polar, blue);
   svgc.add(circle_b);
   lig_build::pos_t water_pos = tlc2 + lig_build::pos_t(0, 4.0 * y_line_step);
   std::string circle_w = make_circle(water_pos, 0.75, 0.1, "white", "#111111");
   svgc.add(circle_w);
   lig_build::pos_t metal_pos = tlc2 + lig_build::pos_t(0, 5.0 * y_line_step);
   std::string circle_m = make_circle(metal_pos, 0.75, 0.1, metalic_grey, "#111111");
   svgc.add(circle_m);

   lig_build::pos_t text_pos_7( top_left_x + 16.2, -top_left_y + 0.2);
   lig_build::pos_t text_pos_8( top_left_x + 16.2, -top_left_y + 0.2 + 1.0 * y_line_step);
   lig_build::pos_t text_pos_9( top_left_x + 16.2, -top_left_y + 0.2 + 2.0 * y_line_step);
   lig_build::pos_t text_pos_10(top_left_x + 16.2, -top_left_y + 0.2 + 3.0 * y_line_step);
   lig_build::pos_t text_pos_11(top_left_x + 16.2, -top_left_y + 0.2 + 4.0 * y_line_step);
   lig_build::pos_t text_pos_12(top_left_x + 16.2, -top_left_y + 0.2 + 5.0 * y_line_step);
   std::string t7  = add_text(text_pos_7,  "Grease");
   std::string t8  = add_text(text_pos_8,  "Polar");
   std::string t9  = add_text(text_pos_9,  "Acidic");
   std::string t10 = add_text(text_pos_10, "Basic");
   std::string t11 = add_text(text_pos_11, "Water");
   std::string t12 = add_text(text_pos_12, "Metal");
   svgc.add(t7);
   svgc.add(t8);
   svgc.add(t9);
   svgc.add(t10);
   svgc.add(t11);
   svgc.add(t12);

   // third column

   lig_build::pos_t sa_pos = top_left.invert_y() + lig_build::pos_t(21.0, 0.0);
   lig_build::pos_t se_pos_1 = sa_pos   + lig_build::pos_t(0, 1.4 * y_line_step);
   lig_build::pos_t se_pos_2 = se_pos_1 + lig_build::pos_t(0.2, 0.2);
   svg_container_t solvent_accessibility = draw_solvent_accessibility_of_atom(sa_pos.invert_y(), 0.15);
   svgc.add(solvent_accessibility);
   std::string circle_se_1 = make_circle(se_pos_1, 0.85, 0.00,  "#b0c0ff", "none");
   std::string circle_se_2 = make_circle(se_pos_2, 0.55, 0.06, "white",   "#111111");
   svgc.add(circle_se_1);
   svgc.add(circle_se_2);
   lig_build::pos_t sc_pos = se_pos_1 + lig_build::pos_t(0, 1.5 * y_line_step);
   std::string substitution_contour = draw_dashed_ellipse(sc_pos);
   svgc.add(substitution_contour);

   // labels for 3rd column

   lig_build::pos_t text_pos_13(top_left_x + 22.4, -top_left_y + 0.2 + 0.0 * y_line_step);
   lig_build::pos_t text_pos_14(top_left_x + 22.4, -top_left_y + 0.2 + 1.5 * y_line_step);
   lig_build::pos_t text_pos_15(top_left_x + 22.4, -top_left_y + 0.2 + 3.0 * y_line_step);
   std::string t13 = add_text(text_pos_13, "Solvent Accessibility");
   std::string t14 = add_text(text_pos_14, "Residue Protection");
   std::string t15 = add_text(text_pos_15, "Substitution Contour");
   svgc.add(t13);
   svgc.add(t14);
   svgc.add(t15);

   return svgc;
}

svg_container_t
pli::fle_view_with_rdkit_internal(mmdb::Manager *mol,
                                  int imol,
                                  coot::protein_geometry *geom_p,
                                  const std::string &chain_id, int res_no, const std::string &ins_code,
                                  float residues_near_radius,
                                  bool add_key) {

   auto write_string_to_file = [] (const std::string &s, const std::string &fn) {

      std::cout << "write string to file! " << fn << std::endl;
      std::ofstream f(fn);
      f << s;
      f << "\n";
      f.close();
   };

   std::string output_image_file_name = "something.svg";
   std::string output_format = "svg"; // was file_format;

   svg_container_t svgc_outer;

   double scale_factor = 400.0;
   bool dark_background_flag = false; // pass this

   double weight_for_3d_distances = 0.4; // for 3d distances

   bool wrap_in_refresh_html = false;

   if (mol) {
      mmdb::Residue  *res_ref = coot::util::get_residue(chain_id, res_no, ins_code, mol);
      mmdb::Manager *mol_for_res_ref = mol;
      if (res_ref) {
         std::string ligand_res_name(res_ref->GetResName());

         std::pair<bool, coot::dictionary_residue_restraints_t> p =
            geom_p->get_monomer_restraints_at_least_minimal(ligand_res_name, imol);

         if (! p.first) {

            std::string s1 = "WARNING:: fle_view_with_rdkit(): ";
            std::string s2 = "WARNING:: ";
            s1 += "Failed to get \nmonomer_restraints for ligand of type ";
            s2 += "Failed to get \nmonomer_restraints for ligand of type ";
            s1 += ligand_res_name;
            s2 += ligand_res_name;
            std::cout << s1 << std::endl;

         } else {

            std::vector<mmdb::Residue *> residues =
               coot::residues_near_residue(res_ref, mol_for_res_ref, residues_near_radius);

            std::cout << "DEBUG:: residues_near_radius " << residues_near_radius << std::endl;
            std::cout << "DEBUG:: Here are the residues near the ligand " << residues.size() << std::endl;
            for (auto res : residues) {
               std::cout << "      residues: " << coot::residue_spec_t(res) << " " << res->GetResName() << std::endl;
            }

            // residues needs to be filtered to remove waters that
            // are not connected to a protein atom.

            std::vector<mmdb::Residue *> filtered_residues =
               coot::filter_residues_by_solvent_contact(res_ref, mol_for_res_ref, residues, 3.25);

            for (auto res : filtered_residues) {
               std::cout << "      filtered_residues: " << coot::residue_spec_t(res) << " " << res->GetResName() << std::endl;
            }

            // for the atoms in the ligand only, for the moment.
            // More would require a class containing with and
            // without solvent accessibilites for each residue.
            // Maybe that can be another function.  And that would
            // need to consider neighbours of neighbours, perhaps
            // done with a larger radius.
            //
            pli::dots_representation_info_t dots;
            std::vector<std::pair<coot::atom_spec_t, float> > s_a_v =
               dots.solvent_accessibilities(res_ref, filtered_residues);

            // for the ligand environment residues:
            std::vector<pli::solvent_exposure_difference_helper_t> sed =
               dots.solvent_exposure_differences(res_ref, filtered_residues);

            try {

               // can throw a std::runtime_error
               RDKit::RWMol rdkm = coot::rdkit_mol(res_ref, imol, *geom_p);

               // assign atom names
               if (int(rdkm.getNumAtoms()) < res_ref->GetNumberOfAtoms()) {
                  std::cout << "WARNING:: failure to construct rdkit molecule " << rdkm.getNumAtoms()
                            << " vs " << res_ref->GetNumberOfAtoms() << std::endl;
               } else {
                  mmdb::PPAtom residue_atoms = 0;
                  int n_residue_atoms;
                  res_ref->GetAtomTable(residue_atoms, n_residue_atoms);

                  // polar Hs only, that is - need new function here.
                  // (can throw a std::exception)
                  coot::undelocalise(&rdkm);
                  coot::assign_formal_charges(&rdkm);
                  coot::remove_non_polar_Hs(&rdkm);

                  // we need to sanitizeMol() after remove_non_polar_Hs, and have
                  // a kekulized representation.
                  // we need to sanitize to get ring info,
                  // then we need to kekulize because we are making a 2D chemical diagram
                  //
                  // failed_op sets set by sanitizeMol() - perhaps we shouldn't ignore it.
                  //
                  unsigned int failed_op_1 = 0;
                  unsigned int failed_op_2 = 0;
                  RDKit::MolOps::sanitizeMol(rdkm, failed_op_1, RDKit::MolOps::SANITIZE_ALL);
                  RDKit::MolOps::sanitizeMol(rdkm, failed_op_2, RDKit::MolOps::SANITIZE_KEKULIZE);

                  if (true)
                     std::cout << "DEBUG:: sanitizeMol() returned with failed_op: "
                               << failed_op_1 << " " << failed_op_2
                               << " (note 'no-failure' is value 0)." << std::endl;

                  try {
                     RDKit::RingInfo *ri = rdkm.getRingInfo();
                  }
                  catch (const std::runtime_error &rte) {
                     std::vector<std::vector<int> > ring_info;
                     RDKit::MolOps::findSSSR(rdkm, ring_info);
                  }


                  int mol_2d_depict_conformer = coot::add_2d_conformer(&rdkm, weight_for_3d_distances);
                  lig_build::molfile_molecule_t m = coot::make_molfile_molecule(rdkm, mol_2d_depict_conformer);

                  flev_t flev;

                  flev.mol.import_rdkit_mol(&rdkm, mol_2d_depict_conformer);
                  // these coordinates are in molecule coordinates based around the middle of the residue
                  if (false) {
                     std::cout << "--------------- pre make_svg() for molecule" << std::endl;
                     std::cout << "there were " << flev.mol.num_atoms() << " in svg mol" << std::endl;
                     unsigned int n = flev.mol.num_atoms();
                     for (unsigned int i=0; i<n; i++) {
                        const auto &atom = flev.mol.atoms[i];
                        std::cout << "atom " << i << " " << atom.atom_id << " " << atom.atom_position << std::endl;
                     }
                  }

                  scale_factor = 1.0;
                  bool add_background = false;
                  svg_container_t svgc_mol = flev.mol.make_svg(scale_factor, dark_background_flag, add_background);

                  if (wrap_in_refresh_html) {
                     // 20250106-PE: hack the bounds for now
                     svgc_mol.set_bounds(-9, -7, 20, 20);
                     std::string s = svgc_mol.compose(true);
                     unsigned int refresh_delta = 1;
                     std::string html_top = "<!DOCTYPE html><html><head><meta http-equiv=\"refresh\" ";
                     html_top += "content=\"" + std::to_string(refresh_delta);
                     html_top += "\"></head><body>\n";
                     std::string html_foot = "</body><html>\n";
                     std::string ss = html_top + s + html_foot;
                     s = ss;
                     write_string_to_file(ss, "flev-test-1.svg.html");
                  }

                  mmdb::Residue *residue_flat = coot::make_residue(rdkm, mol_2d_depict_conformer, "XXX");
                  mmdb::Manager *mol_for_flat_residue = coot::util::create_mmdbmanager_from_residue(residue_flat); // d

                  if (true)
                     for (unsigned int iat=0; iat<m.atoms.size(); iat++)
                        std::cout << iat << "  " << m.atoms[iat] << std::endl;

                  std::string view_name = "Molecule ";
                  view_name += std::to_string(imol);
                  view_name += " ";
                  view_name += chain_id;
                  view_name += std::to_string(res_no);
                  view_name += ins_code;
                  view_name += "    code: ";
                  view_name += ligand_res_name;

                  std::pair<bool, coot::residue_spec_t> ligand_spec_pair(1, coot::residue_spec_t(res_ref));

                  bool use_graphics_flag = false;
                  bool stand_alone_flag = true;

                  // lbg_info_t *lbg_local_p = lbg(m, ligand_spec_pair,
                  //                               NULL, view_name, ligand_res_name, imol,
                  //                               *geom_p,
                  //                               use_graphics_flag, stand_alone_flag,
                  //                               coot_get_url,
                  //                               prodrg_import_function,
                  //                               sbase_import_function,
                  //                               get_drug_mdl_via_wikipedia_and_drugbank);

                   std::map<std::string, std::string> name_map = pli::make_flat_ligand_name_map(res_ref);

                  // should this be a flev function?
                  std::vector<pli::fle_ligand_bond_t> bonds_to_ligand =
                     pli::get_fle_ligand_bonds(res_ref, filtered_residues,
                                               mol_for_res_ref, name_map, *geom_p, imol,
                                               flev.fle_water_dist_max,
                                               flev.fle_h_bond_dist_max);

                  if (true) {
                     for (const auto &btl : bonds_to_ligand) {
                        std::cout << "DEBUG:: bond to ligand xxxxx " << btl << std::endl;
                     }
                  }

                  std::vector<fle_residues_helper_t> res_centres =
                     pli::get_flev_residue_centres(res_ref,
                                                   mol_for_res_ref,
                                                   filtered_residues,
                                                   mol_for_flat_residue);

                  std::vector<int> add_reps_vec;

                  if (false) {
                     for (unsigned int ic=0; ic<res_centres.size(); ic++) {
                        const auto &res_centre = res_centres[ic];
                        std::cout << "  fle_view_with_rdkit_internal(): res_centres: " << ic
                                  << " " << res_centre.spec << " " << res_centre.residue_name
                                  << " " << res_centre.transformed_relative_centre.format()
                                  << std::endl;
                     }
                  }

                  // ----------- residue infos ----------
                  //
                  pli::pi_stacking_container_t pi_stack_info(p.second, filtered_residues, res_ref, rdkm);

                  // ----------- ligand atom infos ------
                  //
                  flev_attached_hydrogens_t ah(p.second);
                  // ah.cannonballs(res_ref, mol_for_res_ref, p.second);
                  ah.distances_to_protein_using_correct_Hs(res_ref, mol_for_res_ref, *geom_p);

                  // ------------ show annotations -----------------
                  //
                  bool annotate_status = flev.annotate(s_a_v, res_centres, add_reps_vec, bonds_to_ligand,
                                                       sed, ah, pi_stack_info, p.second);

                  svgc_outer.add(svgc_mol);
                  svg_container_t svgc_2 = flev.draw_all_flev_annotations();
                  svgc_outer.prepend(svgc_2);

                  if (wrap_in_refresh_html) {
                     // 20250106-PE: hack the bounds for now
                     svgc_outer.set_bounds(-20, -20, 40, 40);
                     std::string s = svgc_outer.compose(true);
                     unsigned int refresh_delta = 1;
                     std::string html_top = "<!DOCTYPE html><html><head><meta http-equiv=\"refresh\" ";
                     html_top += "content=\"" + std::to_string(refresh_delta);
                     html_top += "\"></head><body>\n";
                     std::string html_foot = "</body><html>\n";
                     std::string ss = html_top + s + html_foot;
                     s = ss;
                     write_string_to_file(ss, "flev-test-all-parts-svg.html");
                  }

                  // 20250303-PE this function now returns an svg - let's not
                  // write out a png part-way through now this is in production
                  if (false) {
                     if (output_format == "png") flev.write_png(output_image_file_name);
                     if (output_format == "svg") flev.write_svg(output_image_file_name);
                  }

                  delete mol_for_flat_residue;

                  if (add_key) {
                     double x =  svgc_outer.min_x + 3.0;
                     double y = -svgc_outer.max_y;
                     lig_build::pos_t key_pos(x, y); // top left
                     svg_container_t svgc_key = make_key(key_pos);
                     svgc_outer.add(svgc_key);
                     svgc_outer.add_to_y_bounds(11.0);
                  }
               }
            }
            catch (const std::runtime_error &rte) {
               std::cout << "ERROR:: (runtime error) in fle_view_with_rdkit(): "
                         << rte.what() << std::endl;
            }
            catch (const std::exception &e) {
               std::cout << "ERROR (exception) in fle_view_with_rdkit(): " << e.what() << std::endl;
            }
         }
      }
   }

   // 20250106-PE: hack the bounds for now
   // svgc_outer.set_bounds(-20, -20, 40, 40);
   return svgc_outer;
}

svg_container_t
flev_t::draw_substitution_contour() {

   auto show_grid = [] (const flev_t::ligand_grid &grid) {
      for (int ix=0; ix<grid.x_size(); ix++) {
         for (int iy=0; iy<grid.y_size(); iy++) {
            double I = grid.get(ix, iy);
            std::cout << "substitution-grid " << ix << " " << iy << " " << I << std::endl;
         }
      }
   };

   svg_container_t svgc;

   bool debug = false;

   bool draw_flev_annotations_flag = true; // why wouldn't it be?

   auto debug_bash_distances = [] (const svg_molecule_t &mol) {

      for (unsigned int i=0; i<mol.atoms.size(); i++) {
            std::cout << "in draw_substitution_contour() atom " << i << " "
                      << mol.atoms[i].get_atom_name()
                      << " has "  << mol.atoms[i].bash_distances.size()
                      << " bash distances" << std::endl;
         for (unsigned int j=0; j<mol.atoms[i].bash_distances.size(); j++) {
            std::cout << "  " << mol.atoms[i].bash_distances[j];
         }
         if (! mol.atoms[i].bash_distances.empty())
            std::cout << std::endl;
      }
   };

   auto is_unlimited = [] (const svg_atom_t &atom) {
      unsigned int n_unlimited = 0;
      unsigned int n_bash_distances = atom.bash_distances.size();
      for (unsigned int i=0; i<n_bash_distances; i++)
         if (atom.bash_distances[i].unlimited())
            n_unlimited++;
      float f = static_cast<float>(n_unlimited) / static_cast<float>(n_bash_distances);
      return (f > 0.49999);
   };

   if (draw_flev_annotations_flag) {
      if (mol.atoms.size() > 0) {
         try {
            std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair = mol.ligand_extents();
            ligand_grid grid(l_e_pair.first, l_e_pair.second);
            if (debug)
               debug_bash_distances(mol);

            for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {

               const auto &atom = mol.atoms[iat];
               const auto &atom_position = mol.atoms[iat].atom_position;
               unsigned int n_bash_distances = atom.bash_distances.size();

               if (debug)
                  std::cout << "DEBUG:: n_bash_distances: " << n_bash_distances << std::endl;

               if (n_bash_distances == 0) {
                  if (mol.atoms[iat].element != "H") // checked.
                     grid.add_for_accessibility_no_bash_dist_atom(1.0, atom_position);
               } else {
                  if (is_unlimited(atom)) {
                     grid.add_for_accessibility(1.8, 0.2, atom_position);
                  } else {
                     grid.add_for_accessibility(1.2, 0.05, atom_position);
                  }
               }
            }

            // Put some values around the ring centres too.
            //
            grid.avoid_ring_centres(ring_atoms_list, mol);

            // for debugging
            if (debug)
               show_grid(grid);

            std::vector<lig_build::atom_ring_centre_info_t> unlimited_atoms;
            for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
               const auto &atom = mol.atoms[iat];
               if (is_unlimited(atom)) {
                  lig_build::atom_ring_centre_info_t ua(atom);
                  unlimited_atoms.push_back(ua);
               }
            }

            if (false) {
               std::vector<float> contour_levels = { 0.2, 0.4, 0.6, 0.8, 1.0,
                                                     1.2, 1.4, 1.6, 1.8, 2.0,
                                                     2.2, 2.4, 2.6, 2.8, 3.0,
                                                     3.2, 3.4, 3.6, 3.8, 4.0,
                                                     4.2, 4.4, 4.6, 4.8, 5.0 };

               std::string col = "#bbbbbb";
               bool is_dashed = false;
               for (const auto &level : contour_levels)
                  svgc.add(grid.show_contour(level, is_dashed, col, unlimited_atoms, ring_atoms_list));
            }

            bool is_dashed = true;
            std::string col = "#aaaaaa";
            svgc.add(grid.show_contour(0.5, is_dashed, col, unlimited_atoms, ring_atoms_list)); // was 0.5

            // debug
            // show_unlimited_atoms(unlimited_atoms);
            // show_ring_centres(ring_atoms_list, mol);

            if (false) {
               std::cout << "Here are the "<< unlimited_atoms.size()
                         << " unlimited atoms: " << std::endl;
               for (unsigned int iat=0; iat<unlimited_atoms.size(); iat++)
                  std::cout << "   " << unlimited_atoms[iat] << std::endl;
            }

         }
         catch (const std::runtime_error &rte) {
            std::cout << rte.what() << std::endl;
         }
      }
   }
   // std::cout << "subsitution-contour: " << svgc.svg << std::endl;
   return svgc;
}



// top left and bottom right corners.
//
std::pair<lig_build::pos_t, lig_build::pos_t>
flev_t::flev_residues_extents() const {

   std::pair<lig_build::pos_t, lig_build::pos_t> p; // defaults with (-1, -1) for coordinates

   if (residue_circles.size()) {
      p.first  = lig_build::pos_t( 9999, 9999);
      p.second = lig_build::pos_t(-9999, -9999);
      for (unsigned int i=0; i<residue_circles.size(); i++) {
         const lig_build::pos_t &pos = residue_circles[i].pos;
         if (pos.x < p.first.x) p.first.x = pos.x;
         if (pos.y < p.first.y) p.first.y = pos.y;
         if (pos.x > p.second.x) p.second.x = pos.x;
         if (pos.y > p.second.y) p.second.y = pos.y;
      }
   }
   return p;
}

svg_container_t
flev_t::draw_all_flev_annotations() {

   svg_container_t svgc;
   svg_container_t svgc_ra = draw_all_flev_residue_attribs();
   svg_container_t svgc_li = draw_all_flev_ligand_annotations();
   svgc.add(svgc_ra);
   svgc.prepend(svgc_li);
   return svgc;
}

svg_container_t
flev_t::draw_all_flev_residue_attribs() {

   svg_container_t svgc;
   svg_container_t svgc_rc  = draw_residue_circles(residue_circles, additional_representation_handles);
   svg_container_t svgc_btl = draw_bonds_to_ligand();
   svg_container_t svgc_si  = draw_stacking_interactions(residue_circles);
   svgc.add(svgc_rc);
   svgc.add(svgc_btl);
   svgc.add(svgc_si);
   return svgc;
}

svg_container_t
flev_t::draw_all_flev_ligand_annotations() {

   // return should be prepended

   svg_container_t svgc;
   svg_container_t svg_sc  = draw_substitution_contour();
   svg_container_t svg_saa = draw_solvent_accessibility_of_atoms();
   svgc.add(svg_sc);
   svgc.add(svg_saa);
   return svgc;
}

svg_container_t
flev_t::draw_stacking_interactions(const std::vector<residue_circle_t> &rc) {

   svg_container_t svgc;

   for (unsigned int ires=0; ires<rc.size(); ires++) {
      int st = rc[ires].get_stacking_type();
      clipper::Coord_orth click_pos = rc[ires].residue_centre_real;

      if (false)
         std::cout << ":::::::::::::::::: stacking: ires " << ires << " "
                   << rc[ires].has_ring_stacking_interaction()
                   << std::endl;

      if (rc[ires].has_ring_stacking_interaction()) {

         std::vector<std::string> ligand_ring_atom_names =
            rc[ires].get_ligand_ring_atom_names();
         if ((st == residue_circle_t::PI_PI_STACKING) ||
             (st == residue_circle_t::PI_CATION_STACKING)) {
            try {
               lig_build::pos_t lc = mol.get_ring_centre(ligand_ring_atom_names);
               svg_container_t svgc_asl = draw_annotated_stacking_line(lc, rc[ires].pos, st, click_pos);
               svgc.add(svgc_asl);
            }
            catch (const std::runtime_error &rte) {
               std::cout << rte.what() << std::endl;
            }
         }
      }

      if (st == residue_circle_t::CATION_PI_STACKING) {
         std::string at_name = rc[ires].get_ligand_cation_atom_name();
         lig_build::pos_t atom_pos = mol.get_atom_canvas_position(at_name);
         svg_container_t svgc_asl = draw_annotated_stacking_line(atom_pos, rc[ires].pos, st, click_pos);
         svgc.add(svgc_asl);
      }
   }

   // std::cout << "draw_stacking_interactions returns an svg of size " << svgc.svg.size() << std::endl;
   return svgc;
}


// click_pos is where we recentre in 3D graphics when the annotation
// (line) is clicked.
svg_container_t
flev_t::draw_annotated_stacking_line(const lig_build::pos_t &ligand_ring_centre,
                                     const lig_build::pos_t &residue_pos,
                                     int stacking_type,
                                     const clipper::Coord_orth &click_pos) {

   auto do_polygon = [] (const std::vector<std::pair<double, double> > &hex_points) {
      std::string s = "   <polygon points=\"";
      for (unsigned int i=0; i<hex_points.size(); i++)
         s += std::to_string(hex_points[i].first) + "," + std::to_string(-hex_points[i].second) + " ";
      s += "\" style=\"fill:none;stroke:#109910;stroke-width:0.1\" />\n";
      return s;
   };

   auto make_line = [] (const lig_build::pos_t &p1, const lig_build::pos_t &p2,
                        double line_width, const std::string &stroke_colour, bool dashed) {
      std::string s;
      s += "   <line x1=\"" + std::to_string(p1.x) + "\" y1=\"" + std::to_string(-p1.y) + "\" ";
      s += "x2=\"" + std::to_string(p2.x) + "\" y2=\"" + std::to_string(-p2.y) + "\" ";
      s += "style=\"stroke:" + stroke_colour + ";stroke-width:" + std::to_string(line_width) + ";";
      s += "stroke-linecap:round;";
      if (dashed)
         s += "stroke-dasharray:0.1,0.2;";
      s += "\" />\n";
      return s;
   };

   svg_container_t svgc;
   svgc.add_comment("Stacking Line");

   double ligand_target_shortening_factor = 0.1;
   if (stacking_type == residue_circle_t::CATION_PI_STACKING)
      ligand_target_shortening_factor = 0.2;
   lig_build::pos_t a_to_b_uv = (ligand_ring_centre - residue_pos).unit_vector();
   lig_build::pos_t A = residue_pos + a_to_b_uv * 0.2;
   lig_build::pos_t B = ligand_ring_centre;
   // just short of the middle of the ring
   lig_build::pos_t C = B - a_to_b_uv * ligand_target_shortening_factor;
   lig_build::pos_t mid_pt = (A + B) * 0.5;
   lig_build::pos_t close_mid_pt_1 = mid_pt - a_to_b_uv * 0.92;
   lig_build::pos_t close_mid_pt_2 = mid_pt + a_to_b_uv * 0.92;

   std::string stroke_colour = "#008000";

   std::vector<lig_build::pos_t> hex_and_ring_centre(2);
   hex_and_ring_centre[0] = mid_pt - a_to_b_uv * 0.4;
   hex_and_ring_centre[1] = mid_pt + a_to_b_uv * 0.4;

   // for cations, we draw a 6-member ring and a "+" text.
   //
   // for pi-pi, we draw 2 6-member rings
   //
   for(int ir=0; ir<2; ir++) {

      bool do_ring  = false;
      bool do_anion = false;

      if (stacking_type == residue_circle_t::PI_PI_STACKING) {
         do_ring = true;
      }
      if (stacking_type == residue_circle_t::PI_CATION_STACKING) {
         if (ir == 0) {
            do_ring = false;
            do_anion = true;
         } else {
            do_ring = true;
            do_anion = false;
         }
      }
      if (stacking_type == residue_circle_t::CATION_PI_STACKING) {
         if (ir == 0) {
            do_ring = true;
            do_anion = false;
         } else {
            do_ring = false;
            do_anion = true;
         }
      }

      if (do_ring) {

         if (false)
            std::cout << "--------------------- in the do_ring block residue_pos " << residue_pos
                      << " hex_and_ring_centre[ir] " << hex_and_ring_centre[ir] << std::endl;

         double angle_step = 60;
         double r = 0.4; // radius
         std::vector<std::pair<double, double> > hex_points(6);
         for (int ipt=1; ipt<=6; ipt++) {
            int ipt_0 = ipt - 1;
            double theta_deg_0 = 30 + angle_step * ipt_0;
            double theta_0 = theta_deg_0 * DEG_TO_RAD;
            lig_build::pos_t pt_0 = hex_and_ring_centre[ir] + lig_build::pos_t(sin(theta_0), cos(theta_0)) * r;
            hex_points[ipt_0] = std::make_pair(pt_0.x, pt_0.y);
         }
         std::string ring_text = do_polygon(hex_points);
         svgc.add(ring_text);

         lig_build::pos_t ring_centre = hex_and_ring_centre[ir];
         lig_build::pos_t circle_centre(ring_centre.x, -ring_centre.y);
         std::string circle = pli::make_circle(circle_centre, 0.2, 0.08, "none", stroke_colour);
         svgc.add(circle);

      }

      if (do_anion) {
         // the "+" symbol for the anion.
         // note that the anchor "middle" applies to x, not y
         std::string t("   <text ");
         t += std::string("fill=\"");
         t += stroke_colour;
         t += std::string("\"");
         t += std::string(" x=\"");
         t += std::to_string(hex_and_ring_centre[ir].x);
         t += std::string("\"");
         t += std::string(" y=\"");
         t += std::to_string(-hex_and_ring_centre[ir].y+0.16);
         t += std::string("\"");
         t += std::string(" text-anchor=\"middle\"");
         t += std::string(" font-family=\"Helvetica, sans-serif\" font-size=\"0.06em\">");
         t += std::string("+");
         t += std::string("</text>\n");

         svgc.add(t);
      }
   }

   // now draw the stacking interaction dotted lines
   lig_build::pos_t A_prime = A + a_to_b_uv * 1.0;
   std::string s1 = make_line(A_prime,            close_mid_pt_1, 0.15, stroke_colour, true);
   std::string s2 = make_line(ligand_ring_centre, close_mid_pt_2, 0.15, stroke_colour, true);
   svgc.add(s1);
   svgc.add(s2);

   // Now the circle blob at the centre of the aromatic ligand ring:
   if (stacking_type != residue_circle_t::CATION_PI_STACKING) {

      float radius = 0.2; // 20250118-PE was 67.0;
      float stroke_width = 0.01; // 20250118-PE was 10.0;
      auto ligand_ring_centre_inv_y = lig_build::pos_t(ligand_ring_centre.x, -ligand_ring_centre.y);
      std::string c = pli::make_circle(ligand_ring_centre_inv_y, radius, stroke_width, stroke_colour, stroke_colour);
      svgc.add(c);
   }
   return svgc;
}

std::string make_arrow(const lig_build::pos_t &A, const lig_build::pos_t &B,
                       const std::string &stroke_colour, bool start_arrow, bool end_arrow,
                       const lig_build::pos_t &pos,
                       const lig_build::pos_t &lig_atom_pos) {

   std::string arrow_string;
   lig_build::pos_t rc_to_lig_atom = lig_atom_pos - pos;
   lig_build::pos_t rc_to_lig_atom_uv = rc_to_lig_atom.unit_vector();
   std::string s;
   std::string bond_colour = stroke_colour;
   s += "<!-- Bond to Ligand -->\n";
   s += "   <line x1=\"";
   s += std::to_string(A.x);
   s += "\" y1=\"";
   s += std::to_string(-A.y);
   s += "\" x2=\"";
   s += std::to_string(B.x);
   s += "\" y2=\"";
   s += std::to_string(-B.y);
   s += "\"";
   s += " style=\"stroke:";
   // s += "#202020";
   s += bond_colour;
   s += "; stroke-width:0.15; fill:none; stroke-linecap:round; ";
   s += "\"";
   s += "/>\n";
   arrow_string += s;

   if (start_arrow) {
      // add a triangle arrow head
      lig_build::pos_t tip = pos + rc_to_lig_atom_uv * 1.15;
      double l1 = 0.85;
      double l2 = 0.3;
      lig_build::pos_t rcla_uv_90 = rc_to_lig_atom_uv.rotate(90);
      lig_build::pos_t pt_1 = tip;
      lig_build::pos_t pt_2 = tip + rc_to_lig_atom_uv * l1 + rcla_uv_90 * l2;
      lig_build::pos_t pt_3 = tip + rc_to_lig_atom_uv * l1 - rcla_uv_90 * l2;
      // std::string ts = "<polygon points="100,10 150,190 50,190" style="fill:lime;"/>";
      std::string ts = "   <polygon points =\"";
      ts += std::to_string(pt_1.x);
      ts += ",";
      ts += std::to_string(-pt_1.y);
      ts += " ";
      ts += std::to_string(pt_2.x);
      ts += ",";
      ts += std::to_string(-pt_2.y);
      ts += " ";
      ts += std::to_string(pt_3.x);
      ts += ",";
      ts += std::to_string(-pt_3.y);
      ts += "\" style=\"fill:";
      ts += bond_colour;
      ts += ";\"/>";
      ts += "\n";
      arrow_string += ts;
   }

   if (end_arrow) {
      // add a triangle arrow head
      lig_build::pos_t tip = lig_atom_pos - rc_to_lig_atom_uv * 0.35;
      double l1 = 0.85;
      double l2 = 0.3;
      lig_build::pos_t rcla_uv_90 = rc_to_lig_atom_uv.rotate(90);
      lig_build::pos_t pt_1 = tip;
      lig_build::pos_t pt_2 = tip - rc_to_lig_atom_uv * l1 + rcla_uv_90 * l2;
      lig_build::pos_t pt_3 = tip - rc_to_lig_atom_uv * l1 - rcla_uv_90 * l2;
      // std::string ts = "<polygon points="100,10 150,190 50,190" style="fill:lime;"/>";
      std::string ts = "   <polygon points =\"";
      ts += std::to_string(pt_1.x);
      ts += ",";
      ts += std::to_string(-pt_1.y);
      ts += " ";
      ts += std::to_string(pt_2.x);
      ts += ",";
      ts += std::to_string(-pt_2.y);
      ts += " ";
      ts += std::to_string(pt_3.x);
      ts += ",";
      ts += std::to_string(-pt_3.y);
      ts += "\" style=\"fill:";
      ts += bond_colour;
      ts += ";\"/>";
      ts += "\n";
      arrow_string += ts;
   }
   return arrow_string;
}


svg_container_t
flev_t::draw_bonds_to_ligand() {

   svg_container_t svgc;

   for (unsigned int ic=0; ic<residue_circles.size(); ic++) {
      if (residue_circles[ic].bonds_to_ligand.size()) {
         for (unsigned int ib=0; ib<residue_circles[ic].bonds_to_ligand.size(); ib++) {
            const bond_to_ligand_t &bond_to_ligand = residue_circles[ic].bonds_to_ligand[ib];

            try {

               lig_build::pos_t pos = residue_circles[ic].pos;
               std::string at_name = bond_to_ligand.ligand_atom_name;
               lig_build::pos_t lig_atom_pos = mol.get_atom_canvas_position(at_name);
               lig_build::pos_t rc_to_lig_atom = lig_atom_pos - pos;
               lig_build::pos_t rc_to_lig_atom_uv = rc_to_lig_atom.unit_vector();
               lig_build::pos_t B = lig_atom_pos;

               double shorten_factor = 0.6;
               // a HOH residue cirle bond-to-ligand doesn't have an arrow-head
               if (residue_circles[ic].residue_type == "HOH") shorten_factor = 0.5;
               if (bond_to_ligand.bond_type != pli::fle_ligand_bond_t::BOND_COVALENT)
                  B -= rc_to_lig_atom_uv * shorten_factor; // put the end point (say) of a hydrogen bond
                                                           // some pixels away from the atom centre - for better
                                                           // aesthetics.
               shorten_factor = 1.3;
               // but if it has a start arrow, we need to make it shorter so that the line and
               // the tip of the arrow don't overlap. Who ever will notice, I wonder...
               if (bond_to_ligand.bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_MAINCHAIN) shorten_factor = 1.4;
               if (bond_to_ligand.bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_SIDECHAIN) shorten_factor = 1.4;
               lig_build::pos_t A = pos + rc_to_lig_atom_uv * shorten_factor;

               // some colours
               std::string blue = "blue";
               std::string green = "darkgreen";
               std::string lime = "#888820";
               std::string stroke_colour = "#111111"; // unset

               // arrows (acceptor/donor) and stroke colour (depending
               // on mainchain or sidechain interaction)
               //
               bool start_arrow = false;
               bool   end_arrow = false;
               if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_DONOR_SIDECHAIN) {
                  end_arrow = 1;
                  stroke_colour = green;
                  // dash = dash_dash;
               }
               if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_DONOR_MAINCHAIN) {
                  end_arrow = 1;
                  stroke_colour = blue;
                  // dash = dash_dash;
               }
               if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_SIDECHAIN) {
                  start_arrow = 1;
                  stroke_colour = green;
                  // dash = dash_dash;
               }
               if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_MAINCHAIN) {
                  start_arrow = 1;
                  stroke_colour = blue;
                  // dash = dash_dash;
               }
               if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::METAL_CONTACT_BOND) {
                  stroke_colour = "#990099";
                  // dash = dash_dash;
               }
               if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::BOND_COVALENT) {
                  // dash = goo_canvas_line_dash_new (2, 2.7, 0.1);
                  stroke_colour = "#bb00bb";
               }

               if (residue_circles[ic].residue_type == "HOH") {
                  stroke_colour = lime;
                  start_arrow = false;
                  end_arrow = false;
                  // dash = dash_dash;
               }

               std::string arrow_string = make_arrow(A, B, stroke_colour, start_arrow, end_arrow, pos, lig_atom_pos);
               svgc.add(arrow_string);

            }
            catch (const std::runtime_error &rte) {
               std::cout << "WARNING:: " << rte.what() << std::endl;
            }
         }

      } else {
         if (false)
            std::cout << "... no bond to ligand from residue circle "
                      << residue_circles[ic].residue_label << " "
                      << residue_circles[ic].residue_type << std::endl;
      }
   }
   return svgc;
}


svg_container_t
flev_t::draw_solvent_accessibility_of_atoms() {

   svg_container_t svgc;
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      lig_build::pos_t pos = mol.atoms[iat].atom_position;
      double sa = mol.atoms[iat].get_solvent_accessibility();
      // saa of -1 is "unset"
      if (sa  > 0.0) {
         svg_container_t saa = pli::draw_solvent_accessibility_of_atom(pos, sa);
         svgc.add(saa);
      }
   }
   return svgc;
}

svg_container_t
pli::draw_solvent_accessibility_of_atom(const lig_build::pos_t &pos, double sa) {

   svg_container_t svgc;

   int n_circles = int(sa*40) + 1;    // needs fiddling?
   if (n_circles> 10) n_circles = 10; // needs fiddling?

   for (int i=0; i<n_circles; i++) {
      double rad = 0.1 * double(i+1); // needs fiddling?

      std::string comment = "Solvent Accessibilty of Atom";
      // needs to be drawn first? Use a group for sovlent accessibility circles?
      std::string c = pli::make_circle(pos.invert_y(), rad, 0.0, "#5555cc30", "#5555cc30");
      svgc.add_comment(comment);
      svgc.add(c);
   }

   return svgc;
}

void
flev_t::ligand_grid::add_for_accessibility(double bash_dist,
                                           double exp_fac,
                                           const lig_build::pos_t &atom_pos) {

   bool debug = false;
   int feature_extent = 45;

   double inv_scale_factor = 1.0;

   for (int ipos_x= -feature_extent; ipos_x<=feature_extent; ipos_x++) {
      for (int ipos_y= -feature_extent; ipos_y<=feature_extent; ipos_y++) {
         std::pair<int, int> p = mol_space_pos_to_grid_pos(atom_pos);
         int ix_grid = ipos_x + p.first;
         int iy_grid = ipos_y + p.second;
         if ((ix_grid >= 0) && (ix_grid < x_size())) {
            if ((iy_grid >= 0) && (iy_grid < y_size())) {

               double d2 = (grid_pos_to_mol_space_pos(ix_grid, iy_grid) - atom_pos).lengthsq();
               double f = d2 / (inv_scale_factor * inv_scale_factor);
               // double val = substitution_value(f, bash_dist); // expire this 20250214-PE
               double v = 0.04 * f;
               double A = bash_dist;
               double ff = exp_fac;
               double val = A * exp(- v/ff);
               if (debug)
                  if (val > 0)
                     std::cout << "add_for_accessibility(): adding " << val
                               << " to grid " << ix_grid << " " << iy_grid
                               << " from " << sqrt(d2) << " vs " << bash_dist << std::endl;
               grid_[ix_grid][iy_grid] += val;
            }
         }
      }
   }
}

void
flev_t::ligand_grid::add_for_accessibility_no_bash_dist_atom(double scale,
                                                             const lig_build::pos_t &atom_pos) {

   int feature_extent = 40;

   for (int ipos_x= -feature_extent; ipos_x<=feature_extent; ipos_x++) {
      for (int ipos_y= -feature_extent; ipos_y<=feature_extent; ipos_y++) {
         std::pair<int, int> p = mol_space_pos_to_grid_pos(atom_pos);
         int ix_grid = ipos_x + p.first;
         int iy_grid = ipos_y + p.second;
         if ((ix_grid >= 0) && (ix_grid < x_size())) {
            if ((iy_grid >= 0) && (iy_grid < y_size())) {

               double d2 = (grid_pos_to_mol_space_pos(ix_grid, iy_grid) - atom_pos).lengthsq();
               double v = d2 / ( 1.0 * 1.0);
               double A = 1.0;
               double ff = 0.95;
               double val = A * exp(- v/ff);
               grid_[ix_grid][iy_grid] += val;
            }
         }
      }
   }
}

flev_t::ligand_grid::ligand_grid(const lig_build::pos_t &low_x_and_y,
                                 const lig_build::pos_t &high_x_and_y) {

   extra_extents = 10; // grid points, added to both sides
   n_grid_per_angstrom = 5.0;
   scale_fac = 5; // same thing? (old and not used?)

   ligand_atoms_min_x =  low_x_and_y.x;
   ligand_atoms_min_y =  low_x_and_y.y;
   ligand_atoms_max_x = high_x_and_y.x;
   ligand_atoms_max_y = high_x_and_y.y;

   mol_space_grid_min_x = ligand_atoms_min_x - extra_extents / n_grid_per_angstrom;
   mol_space_grid_min_y = ligand_atoms_min_y - extra_extents / n_grid_per_angstrom;

   double mol_space_grid_max_x = ligand_atoms_max_x + extra_extents / n_grid_per_angstrom;
   double mol_space_grid_max_y = ligand_atoms_max_y + extra_extents / n_grid_per_angstrom;

   if (false) {
      std::cout << "in constructor with ligand_atoms_min_x " << ligand_atoms_min_x << std::endl;
      std::cout << "in constructor with ligand_atoms_min_y " << ligand_atoms_min_y << std::endl;

      std::cout << "in constructor with mol_space_grid_min_x " << mol_space_grid_min_x << std::endl;
      std::cout << "in constructor with mol_space_grid_min_y " << mol_space_grid_min_x << std::endl;
   }

   double delta_x = mol_space_grid_max_x - mol_space_grid_min_x;
   double delta_y = mol_space_grid_max_y - mol_space_grid_min_y;
   x_size_ = int(delta_x * n_grid_per_angstrom + 2 * extra_extents) + 1;
   y_size_ = int(delta_y * n_grid_per_angstrom + 2 * extra_extents) + 1;

   std::vector<double> tmp_y(y_size_, 0.0);
   grid_.resize(x_size_);
   for (int i=0; i<x_size_; i++)
      grid_[i] = tmp_y;

}


void
flev_t::ligand_grid::print(int primary_index) const {

   int xs = x_size();
   int ys = y_size();

   std::string file_name = "ligand-grid-" + std::to_string(primary_index) + ".table";
   std::ofstream f(file_name);
   if (f) {
      for (int ix=0; ix<xs; ix++) {
         for (int iy=0; iy<ys; iy++) {
            double v = get(ix, iy);
            f << ix << " " << iy << " " << v << "\n";
         }
      }
      f.close();
   }
}

void
flev_t::ligand_grid::avoid_ring_centres(const std::vector<std::vector<std::string> > &ring_atoms_list,
                                        const lig_build::molecule_t<svg_atom_t, svg_bond_t> &mol) {

   // For the substitution contour we blob in a circle of radius
   // 1/(2*sin(180/n_ring_atoms)) bond lengths.

   for (unsigned int iring=0; iring<ring_atoms_list.size(); iring++) {
      try {
         lig_build::pos_t centre = mol.get_ring_centre(ring_atoms_list[iring]);
         int n_atoms = ring_atoms_list[iring].size();
         // just in case we get here with n_atoms = 1 for some reason...
         if (n_atoms < 3) n_atoms = 3;

         double radius = 1/(2*sin(M_PI/double(n_atoms))) * 1.5; // in "A" or close
         // std::cout << "avoid_ring_centres() adding ring centre at " << centre
         // << " n_atoms: " << n_atoms << " radius " << radius << std::endl;
         add_for_accessibility(radius, 0.1, centre);
      }
      catch (const std::runtime_error &rte) {
         std::cout << "Opps - failed to find ring centre for ring atom name "
                   << iring << std::endl;
      }
   }
}

void
flev_t::ligand_grid::show_contour(float contour_level) {

   std::vector<lig_build::atom_ring_centre_info_t> dummy_unlimited_atoms;
   std::vector<std::vector<std::string> > dummy_ring_atom_names;
   std::string col = "#888888";
   bool is_dashed = false;
   show_contour(contour_level, is_dashed, col, dummy_unlimited_atoms, dummy_ring_atom_names);
}

svg_container_t
flev_t::ligand_grid::show_contour(float contour_level,
                                  bool is_dashed,
                                  const std::string &col,
                                  const std::vector<lig_build::atom_ring_centre_info_t> &unlimited_atoms,
                                  const std::vector<std::vector<std::string> > &ring_atoms_list) {

   auto contour_lines_to_svg = [] (const std::vector<std::vector<lig_build::pos_t> > &contour_lines,
                                   const std::string &col,
                                   bool is_dashed) {
      svg_container_t svgc;
      for(const auto &pos_vec : contour_lines) {
         unsigned int n_pos = pos_vec.size();
         if (n_pos > 1) {
            svgc.add_comment("Substitution Contour");
            for(unsigned int i=0; i<(n_pos-1); i++) {
               const auto &p_1 = pos_vec[i];
               const auto &p_2 = pos_vec[i+1];
               svgc.add_line(p_1, p_2, 0.1, col, is_dashed);
            }
         }
      }
      return svgc;
   };

   bool debug = false;

   // fill the ring centre vector, if the unlimited atom have ring centres
   std::vector<std::pair<bool, lig_build::pos_t> > ring_centres(unlimited_atoms.size());
   for (unsigned int i=0; i<unlimited_atoms.size(); i++) {
      try {
         lig_build::pos_t p;
      }
      catch (const std::runtime_error &rte) {
         std::cout << rte.what() << std::endl;
      }
   }

   std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > line_fragments;

   for (int ix=0; ix<x_size()-1; ix+=1) {
      for (int iy=0; iy<y_size()-1; iy+=1) {
         int ms_type = square_type(ix, iy, contour_level);

         if ((ms_type != MS_NO_CROSSING) && (ms_type != MS_NO_SQUARE)) {

            double v00 = get(ix,   iy);   // the values of the grid at these positions
            double v01 = get(ix,   iy+1);
            double v10 = get(ix+1, iy);
            double v11 = get(ix+1, iy+1);

            contour_fragment cf(ms_type, contour_level, v00, v01, v10, v11);

            if (cf.coords_size() == 1) {

               std::pair<double, double> xy_1 = cf.get_coords(ix, iy, 0);
               std::pair<double, double> xy_2 = cf.get_coords(ix, iy, 1);

               if (debug)
                  std::cout << "plot_contour A "
                            << xy_1.first  << " "
                            << xy_1.second << " to "
                            << xy_2.first  << " "
                            << xy_2.second << "  " << ms_type << std::endl;

               lig_build::pos_t pos_1 = grid_pos_as_double_to_mol_space_pos(xy_1.first, xy_1.second);
               lig_build::pos_t pos_2 = grid_pos_as_double_to_mol_space_pos(xy_2.first, xy_2.second);

               std::pair<lig_build::pos_t, lig_build::pos_t> fragment_pair(pos_1, pos_2);
               if (debug)
                  std::cout << "plot_contour B "
                            << pos_1.x << " " << pos_1.y << " to "
                            << pos_2.x << " " << pos_2.y << std::endl;

               // Now filter out this fragment pair if it is too close
               // to an unlimited_atom_positions
               bool plot_it = true;
               double dist_crit = 3.3; // should be/was 4.0

               for (unsigned int i=0; i<unlimited_atoms.size(); i++) {

                  lig_build::pos_t p = unlimited_atoms[i].atom.atom_position;

                  // if this atom has a ring centre, use the ring
                  // centre to atom vector to unplot vectors only in a
                  // particular direction.
                  //
                  if ((p - pos_1).lengthsq() < (dist_crit * dist_crit)) {
                     if (1) { // for debugging
                        if (unlimited_atoms[i].has_ring_centre_flag) {
                           // std::cout << " atom " << i << " has ring_centre ";
                           lig_build::pos_t d_1 =
                              unlimited_atoms[i].ring_centre - unlimited_atoms[i].atom.atom_position;
                           lig_build::pos_t d_2 = unlimited_atoms[i].atom.atom_position - pos_1;
                           double cos_theta = lig_build::pos_t::dot(d_1, d_2)/(d_1.length()*d_2.length());
                           // std::cout << " cos_theta " << cos_theta << " for unlimited atom " << i << std::endl;
                           if (cos_theta > 0.3) { // only cut in the "forwards" direction
                              plot_it = false;
                              break;
                           }
                           // std::cout << std::endl;

                        } else {
                           plot_it = false;
                           break;
                        }
                     }
                  }
               } // end unlimited atoms loop

               if (plot_it)
                  line_fragments.push_back(fragment_pair);
            }
         }
      }
   }

   std::vector<std::vector<lig_build::pos_t> > contour_lines = make_contour_lines(line_fragments);
   svg_container_t svgc = contour_lines_to_svg(contour_lines, col, is_dashed);
   return svgc;

}



double
flev_t::ligand_grid::substitution_value(double r_squared, double bash_dist) const {

   double D = bash_dist;
   double r = sqrt(r_squared);

   if (bash_dist < 1) {
      // this is not in C&L, so this is an enhancement of their
      // algorithm.  Without this there is non-smoothness as the
      // contour follows the interface between 0 and 1.

      // So we add a nice smooth slope (similar to that of
      // conventional bash distances (below).  However we add a slope
      // between +/- 0.2 of the bash distance.
      //
      double small = 0.2;
      if (r > (D + small)) {
         return 0;
      } else {
         if (r < (D - small)) {
            return 1;
         } else {
            double m = (r-(D-small))/(2*small);
            return (0.5 * (1 + cos(M_PI * m)));
         }
      }

   } else {

      if (r<1)
         return 1;
      if (r<(D-1)) {
         return 1;
      } else {
         if (r > D) {
            return 0;
         } else {
            // double v = 0.5*(1 + cos(0.5 *M_PI * (D-r))); // C&L - eh? typo/bug
            double v = 0.5 * (1.0 + cos(M_PI * (D-1-r)));
            return v;
         }
      }
   }
}


std::vector<std::vector<lig_build::pos_t> >
flev_t::ligand_grid::make_contour_lines(const std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > &line_fragments) const {

   // Look for neighboring cells so that a continuous line fragment is
   // generated. That is fine for closed contours, but this handled
   // badly contour line fragments that go over the edge, they will be
   // split into several line fragment - depending on which fragment
   // is picked first.
   //
   // Perhaps the last items in the line_frag_queue should be (any)
   // points that are on the edge of the grid - because they are the
   // best starting place for dealing with contour line fragments that
   // go over the edge.  If you do this, then the input arg
   // line_fragments should carry with it info about this line
   // fragment being on an edge (generated in show_contour())..
   //
   
   std::vector<std::vector<lig_build::pos_t> > v; // returned item.

   std::deque<std::pair<lig_build::pos_t, lig_build::pos_t> > line_frag_queue;
   for (unsigned int i=0; i<line_fragments.size(); i++)
      line_frag_queue.push_back(line_fragments[i]);

   while (!line_frag_queue.empty()) {

      // start a new line fragment
      // 
      std::vector<lig_build::pos_t> working_points;
      std::pair<lig_build::pos_t, lig_build::pos_t> start_frag = line_frag_queue.back();
      working_points.push_back(start_frag.first);
      working_points.push_back(start_frag.second);
      line_frag_queue.pop_back();
      bool line_fragment_terminated = 0;

      while (! line_fragment_terminated) { 
	 // Is there a fragment in line_frag_queue that has a start that
	 // is working_points.back()?
	 //
	 bool found = 0;
	 std::deque<std::pair<lig_build::pos_t, lig_build::pos_t> >::iterator it;
	 for (it=line_frag_queue.begin(); it!=line_frag_queue.end(); ++it) { 
	    if (it->first.near_point(working_points.back(), 0.1)) {
	       working_points.push_back(it->second);
	       line_frag_queue.erase(it);
	       found = 1;
	       break;
	    }
	 }
	 if (! found)
	    line_fragment_terminated = 1;
      }

      v.push_back(working_points);
   }
   
   return v;
} 

void
flev_t::ligand_grid::plot_contour_lines(const std::vector<std::vector<lig_build::pos_t> > &contour_lines) {

#if 0
   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color",
						// "#111111", not so dark
						"#777777",
						NULL);
   GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 1.5, 2.5);

   for (unsigned int i=0; i<contour_lines.size(); i++) {
      for (int j=0; j<int(contour_lines[i].size()-1); j++) {

	 goo_canvas_polyline_new_line(group,
				      contour_lines[i][j].x,
				      contour_lines[i][j].y,
				      contour_lines[i][j+1].x,
				      contour_lines[i][j+1].y,
				      "line_width", 1.0,
				      "line-dash", dash,
				      NULL);
      }
   }
   goo_canvas_line_dash_unref(dash);
#endif

}

flev_t::contour_fragment::contour_fragment(int ms_type,
                                           const float &contour_level,
                                           double v00, double v01, double v10, double v11) {

   int ii_next = grid_index_t::INVALID_INDEX;
   int jj_next = grid_index_t::INVALID_INDEX;

   // now done by caller
   // float v00 = grid.get(grid_index.i(),   grid_index.j());
   // float v01 = grid.get(grid_index.i(),   grid_index.j()+1);
   // float v10 = grid.get(grid_index.i()+1, grid_index.j());
   // float v11 = grid.get(grid_index.i()+1, grid_index.j()+1);

   double frac_x1 = -1;
   double frac_y1 = -1;
   double frac_x2 = -1;  // for hideous valley
   double frac_y2 = -1;

   contour_fragment::coordinates c1(0.0, X_AXIS_LOW);
   contour_fragment::coordinates c2(Y_AXIS_LOW, 0.0);
   cp_t p(c1,c2);

   switch (ms_type) {

   case ligand_grid::MS_UP_0_0:
   case ligand_grid::MS_UP_0_1_and_1_0_and_1_1:

      frac_x1 = (v00-contour_level)/(v00-v10);
      frac_y1 = (v00-contour_level)/(v00-v01);
      c1 = contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = contour_fragment::coordinates(Y_AXIS_LOW, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case ligand_grid::MS_UP_0_1:
   case ligand_grid::MS_UP_0_0_and_1_0_and_1_1:

      // these look fine
      // std::cout << " ----- case MS_UP_0,1 " << std::endl;
      frac_y1 = (v00 - contour_level)/(v00-v01);
      frac_x2 = (v01 - contour_level)/(v01-v11);
      c1 = contour_fragment::coordinates(Y_AXIS_LOW, frac_y1);
      c2 = contour_fragment::coordinates(frac_x2, X_AXIS_HIGH);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case ligand_grid::MS_UP_1_0:
   case ligand_grid::MS_UP_0_0_and_0_1_and_1_1:

      // std::cout << " ----- case MS_UP_1,0 " << std::endl;
      frac_x1 = (contour_level - v00)/(v10-v00);
      frac_y1 = (contour_level - v10)/(v11-v10);
      c1 = contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = contour_fragment::coordinates(Y_AXIS_HIGH, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;


   case ligand_grid::MS_UP_1_1:
   case ligand_grid::MS_UP_0_0_and_0_1_and_1_0:

      // std::cout << " ----- case MS_UP_1,1 " << std::endl;
      frac_x1 = (v01-contour_level)/(v01-v11);
      frac_y1 = (v10-contour_level)/(v10-v11);
      c1 = contour_fragment::coordinates(frac_x1, X_AXIS_HIGH);
      c2 = contour_fragment::coordinates(Y_AXIS_HIGH, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case ligand_grid::MS_UP_0_0_and_0_1:
   case ligand_grid::MS_UP_1_0_and_1_1:

      // std::cout << " ----- case MS_UP_0,0 and 0,1 " << std::endl;
      frac_x1 = (v00-contour_level)/(v00-v10);
      frac_x2 = (v01-contour_level)/(v01-v11);
      c1 = contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = contour_fragment::coordinates(frac_x2, X_AXIS_HIGH);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case ligand_grid::MS_UP_0_0_and_1_0:
   case ligand_grid::MS_UP_0_1_and_1_1:

      // std::cout << " ----- case MS_UP_0,0 and 1,0 " << std::endl;
      frac_y1 = (v00-contour_level)/(v00-v01);
      frac_y2 = (v10-contour_level)/(v10-v11);
      c1 = contour_fragment::coordinates(Y_AXIS_LOW,  frac_y1);
      c2 = contour_fragment::coordinates(Y_AXIS_HIGH, frac_y2);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   default:
      std::cout << "ERROR:: unhandled square type: " << ms_type << std::endl;

   }

}


std::pair<double, double>
flev_t::contour_fragment::get_coords(int ii, int jj, int coord_indx) const {

   coordinates c;

   if (coord_indx == 0)
      if (coords.size() == 0)
         std::cout << "disaster A in get_coords()" << std::endl;

   if (coord_indx == 1)
      if (coords.size() == 0)
         std::cout << "disaster B in get_coords()" << std::endl;

   if (coord_indx == 0)
      c = coords[0].first;
   if (coord_indx == 1)
      c = coords[0].second;

   // these are for hideous value (two crossing vectors)
   if (coord_indx == 2)
      c = coords[1].first;
   if (coord_indx == 3)
      c = coords[1].second;

   double iid = static_cast<double>(ii);
   double jjd = static_cast<double>(jj);
   // std::cout << "get_coords for coord_indx " << coord_indx << " parts " << iid << " " << jjd
   //           << " fracs: " << c.get_frac_x() << " " << c.get_frac_y() << std::endl;
   return std::pair<double, double> (iid+c.get_frac_x(), jjd+c.get_frac_y());
}


// for marching squares, ii and jj are the indices of the bottom left-hand side.
int
flev_t::ligand_grid::square_type(int ii, int jj, float contour_level) const {

   int square_type = ligand_grid::MS_NO_SQUARE;
   if ((ii+1) >= x_size_) {
      return ligand_grid::MS_NO_SQUARE;
   } else {
      if ((jj+1) >= y_size_) {
	 return ligand_grid::MS_NO_SQUARE;
      } else {
	 float v00 = get(ii, jj);
	 float v01 = get(ii, jj+1);
	 float v10 = get(ii+1, jj);
	 float v11 = get(ii+1, jj+1);
	 if (v00 > contour_level) {
	    if (v01 > contour_level) {
	       if (v10 > contour_level) {
		  if (v11 > contour_level) {
		     return ligand_grid::MS_NO_CROSSING;
		  }
	       }
	    }
	 }
	 if (v00 < contour_level) { 
	    if (v01 < contour_level) { 
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return ligand_grid::MS_NO_CROSSING;
		  }
	       }
	    }
	 }

	 // OK, so it is not either of the trivial cases (no
	 // crossing), there are 14 other variants.
	 // 
	 if (v00 < contour_level) { 
	    if (v01 < contour_level) { 
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return ligand_grid::MS_NO_CROSSING;
		  } else {
		     return ligand_grid::MS_UP_1_1;
		  }
	       } else {
		  if (v11 < contour_level) {
		     return ligand_grid::MS_UP_1_0;
		  } else {
		     return ligand_grid::MS_UP_1_0_and_1_1;
		  }
	       }
	    } else {

	       // 0,1 is up
	       
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return ligand_grid::MS_UP_0_1;
		  } else {
		     return ligand_grid::MS_UP_0_1_and_1_1;
		  }
	       } else {
		  if (v11 < contour_level) {
		     return ligand_grid::MS_UP_0_1_and_1_0;      // hideous valley
		  } else {
		     return ligand_grid::MS_UP_0_1_and_1_0_and_1_1; // (only 0,0 down)
		  }
	       }
	    }
	 } else {

	    // 0,0 is up

	    if (v01 < contour_level) {
	       if (v10 < contour_level) {
		  if (v11 < contour_level) {
		     return ligand_grid::MS_UP_0_0;
		  } else {
		     return ligand_grid::MS_UP_0_0_and_1_1; // another hideous valley
		  }
	       } else {
		  // 1,0 is up
		  if (v11 < contour_level) {
		     return ligand_grid::MS_UP_0_0_and_1_0;
		  } else {
		     return ligand_grid::MS_UP_0_0_and_1_0_and_1_1; // 0,1 is down
		  }
	       }
	    } else {

	       // 0,1 is up

	       if (v10 < contour_level) {
		  if (v11 < contour_level) {
		     return ligand_grid::MS_UP_0_0_and_0_1;
		  } else {
		     return ligand_grid::MS_UP_0_0_and_0_1_and_1_1; // 1,0 is down
		  }
	       } else {
		  // if we get here, this test must pass.
		  if (v11 < contour_level) {
		     return ligand_grid::MS_UP_0_0_and_0_1_and_1_0; // only 1,1 is down
		  }
	       }
	    }
	 }
      }
   }
   return square_type;
}
// scale peak value to 1.0

void
flev_t::ligand_grid::normalize() {

   double max_int = 0.0;

   // std::cout << "normalizing grid " << x_size() << " by " << y_size() << std::endl;
   for (int ix=0; ix<x_size(); ix++) {
      for (int iy=0; iy<y_size(); iy++) {
         double intensity = grid_[ix][iy];
         if (intensity > max_int)
            max_int = intensity;
      }
   }
   if (max_int > 0.0) {
      double sc_fac = 1.0/max_int;
      for (int ix=0; ix<x_size(); ix++) {
         for (int iy=0; iy<y_size(); iy++) {
            grid_[ix][iy] *= sc_fac;
         }
      }
   }
}

// --------------------------- coordinate transformation -----------------------------
std::pair<int, int>
flev_t::ligand_grid::mol_space_pos_to_grid_pos(const lig_build::pos_t &pos) const {

   double gminx = mol_space_grid_min_x;
   double gminy = mol_space_grid_min_y;
   lig_build::pos_t grid_offset(gminx, gminy);

   lig_build::pos_t p = pos - grid_offset;
   int nx = p.x * n_grid_per_angstrom;
   int ny = p.y * n_grid_per_angstrom;
   return std::pair<int, int> (nx, ny);
}

// --------------------------- coordinate transformation -----------------------------
lig_build::pos_t
flev_t::ligand_grid::grid_pos_to_mol_space_pos(int ix, int iy) const {

   double gminx = mol_space_grid_min_x;
   double gminy = mol_space_grid_min_y;
   lig_build::pos_t grid_offset(gminx, gminy);

   double rx = static_cast<double>(ix) / n_grid_per_angstrom;
   double ry = static_cast<double>(iy) / n_grid_per_angstrom;
   lig_build::pos_t p(rx, ry);
   lig_build::pos_t d = p + grid_offset;
   return d;
}

// --------------------------- coordinate transformation -----------------------------
lig_build::pos_t
flev_t::ligand_grid::grid_pos_as_double_to_mol_space_pos(double x, double y) const {

   double gminx = mol_space_grid_min_x;
   double gminy = mol_space_grid_min_y;
   lig_build::pos_t grid_offset(gminx, gminy);

   double rx = x / n_grid_per_angstrom;
   double ry = y / n_grid_per_angstrom;
   lig_build::pos_t p(rx, ry);
   lig_build::pos_t d = p + grid_offset;
   return d;
}


// 20241005-PE note to self get_ring_centre() caches the result, so we can't user a const mol here
void
flev_t::ligand_grid::fill(svg_molecule_t mol) {

   auto print_grid = [] (const std::vector<std::vector<double> > &grid, const std::string &label) {
      int grid_size = grid.size();
      for (int ipos_x = 0; ipos_x < grid_size; ipos_x++) {
         int grid_x_size = grid[ipos_x].size();
         for (int ipos_y= 0; ipos_y<= grid_x_size; ipos_y++) {
            double g = grid[ipos_x][ipos_y];
            std::cout << "ligand_grid::fill() " << label << " " << ipos_x << " " << ipos_y << " " << g << "\n";
         }
      }
   };

   double exp_scale = 0.0011;
   exp_scale = 1.0;
   double rk = 3000.0;

   // int grid_extent = 15; // 10, 12 is not enough
   // int grid_extent = 50 ; // untraps 2wot residues?


   if (false) { // debug
      for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
         for (int ipos_x=0; ipos_x<x_size_; ipos_x++) {
            for (int ipos_y=0; ipos_y<y_size_; ipos_y++) {
               lig_build::pos_t mol_space_pos_for_grid_point = grid_pos_to_mol_space_pos(ipos_x, ipos_y);
               std::pair<int, int> gp = mol_space_pos_to_grid_pos(mol_space_pos_for_grid_point);
               std::cout << "grid: XXX: " << ipos_x << " " << gp.first << " YYY: " << ipos_y << " " << gp.second
                         << std::endl;
            }
         }
      }
   }

   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      const auto &atom_pos = mol.atoms[iat].atom_position;
      for (int ipos_x=0; ipos_x<x_size_; ipos_x++) {
         for (int ipos_y=0; ipos_y<y_size_; ipos_y++) {
            lig_build::pos_t mol_space_pos_for_grid_point = grid_pos_to_mol_space_pos(ipos_x, ipos_y);
            lig_build::pos_t delta = mol_space_pos_for_grid_point - atom_pos;
            double d2 = delta.lengthsq();
            double val =  rk * exp(-0.5*exp_scale*d2);
            grid_[ipos_x][ipos_y] += val;
            if (iat == 1) {
               double d = sqrt(d2);
               if (false)
                  std::cout << "debug-grid: ipos_x " << ipos_x << " ipos_y " << ipos_y
                            << " atom_pos.x " << atom_pos.x << " atom_pos.y " << atom_pos.y
                            << " msgp.x " << mol_space_pos_for_grid_point.x << " "
                            << " msgp.y " << mol_space_pos_for_grid_point.y << " "
                            << " delta.x " << delta.x << " "
                            << " delta.y " << delta.y << " "
                            << " d " << d << std::endl;
            }
         }
      }
   }

   if (false)
      print_grid(grid_, "A");

   std::vector<lig_build::pos_t> mol_ring_centres = mol.get_ring_centres();
   // std::cout << "DEBUG:: found " << mol_ring_centres.size() << " ring centres " << std::endl;
   for (unsigned int ir=0; ir<mol_ring_centres.size(); ir++) {
      for (int ipos_x=0; ipos_x<x_size_; ipos_x++) {
         for (int ipos_y=0; ipos_y<y_size_; ipos_y++) {
            lig_build::pos_t mol_space_pos_for_grid_point = grid_pos_to_mol_space_pos(ipos_x, ipos_y);
            lig_build::pos_t delta = mol_space_pos_for_grid_point - mol_ring_centres[ir];
            double d2 = delta.lengthsq();
            double val =  rk * exp(-0.5*exp_scale*d2);
            grid_[ipos_x][ipos_y] += val;
         }
      }
   }

   // print_grid(grid_, "B");

   //  normalize(); // scaled peak value to 1.
}




// Return the fill colour and the stroke colour.
//
std::pair<std::string, std::string>
flev_t::get_residue_circle_colour(const std::string &residue_type) const {

   std::string fill_colour = "#cccccc";
   std::string stroke_colour = "#111111";

   std::string grease = "#ccffbb";
   std::string purple = "#eeccee";
   std::string red    = "#cc0000";
   std::string blue   = "#0000cc";
   std::string metalic_grey = "#d9d9d9";

   if (residue_type == "ALA")
      fill_colour = grease;
   if (residue_type == "TRP")
      fill_colour = grease;
   if (residue_type == "PHE")
      fill_colour = grease;
   if (residue_type == "LEU")
      fill_colour = grease;
   if (residue_type == "PRO")
      fill_colour = grease;
   if (residue_type == "ILE")
      fill_colour = grease;
   if (residue_type == "VAL")
      fill_colour = grease;
   if (residue_type == "MET")
      fill_colour = grease;
   if (residue_type == "MSE")
      fill_colour = grease;

   if (residue_type == "GLY")
      fill_colour = purple;
   if (residue_type == "ASP")
      fill_colour = purple;
   if (residue_type == "ASN")
      fill_colour = purple;
   if (residue_type == "CYS")
      fill_colour = purple;
   if (residue_type == "GLN")
      fill_colour = purple;
   if (residue_type == "GLU")
      fill_colour = purple;
   if (residue_type == "HIS")
      fill_colour = purple;
   if (residue_type == "LYS")
      fill_colour = purple;
   if (residue_type == "LYS")
      fill_colour = purple;
   if (residue_type == "ARG")
      fill_colour = purple;
   if (residue_type == "SER")
      fill_colour = purple;
   if (residue_type == "THR")
      fill_colour = purple;
   if (residue_type == "TYR")
      fill_colour = purple;
   if (residue_type == "HOH")
      fill_colour = "white";

   if (residue_type == "ASP")
      stroke_colour = red;
   if (residue_type == "GLU")
      stroke_colour = red;
   if (residue_type == "LYS")
      stroke_colour = blue;
   if (residue_type == "ARG")
      stroke_colour = blue;
   if (residue_type == "HIS")
      stroke_colour = blue;

   // Bases
   if (residue_type == "U")  fill_colour = purple;
   if (residue_type == "T")  fill_colour = purple;
   if (residue_type == "C")  fill_colour = purple;
   if (residue_type == "A")  fill_colour = purple;
   if (residue_type == "G")  fill_colour = purple;
   if (residue_type == "DT")  fill_colour = purple;
   if (residue_type == "DC")  fill_colour = purple;
   if (residue_type == "DA")  fill_colour = purple;
   if (residue_type == "DG")  fill_colour = purple;

   // metals
   if (residue_type == "ZN")
      fill_colour = metalic_grey;
   if (residue_type == "MG")
      fill_colour = metalic_grey;
   if (residue_type == "NA")
      fill_colour = metalic_grey;
   if (residue_type == "CA")
      fill_colour = metalic_grey;
   if (residue_type == "K")
      fill_colour = metalic_grey;
   if (residue_type == "PO4")
      fill_colour = "#ffcc90";
   if (residue_type == "SO4")
      fill_colour = "#ffff80";

   return std::pair<std::string, std::string> (fill_colour, stroke_colour);
}



// if you don't have add_rep_handles, then pass a vector or size 0.
//
svg_container_t
flev_t::draw_residue_circles(const std::vector<residue_circle_t> &l_residue_circles,
                             const std::vector<int> &add_rep_handles) {

   svg_container_t svgc;

   if (false)
      std::cout << "debug:: Here we are in draw_residue_circles() "
                << l_residue_circles.size() << " " << add_rep_handles.size() << std::endl;

   double max_dist_water_to_ligand_atom  = 3.3; // don't draw waters that are far from ligand
   double max_dist_water_to_protein_atom = 3.3; // don't draw waters that are not somehow
                                                // attached to the protein.

   bool draw_flev_annotations_flag = true; // is this a class member now?

   if (draw_flev_annotations_flag) {
      bool draw_solvent_exposures = true;
      try {
	 lig_build::pos_t ligand_centre = mol.get_ligand_centre();

	 if (draw_solvent_exposures)
	    for (unsigned int i=0; i<l_residue_circles.size(); i++) {
	       svg_container_t svgc_ec = draw_solvent_exposure_circle(l_residue_circles[i], ligand_centre);
               svgc.add(svgc_ec);
            }

	 for (unsigned int i=0; i<l_residue_circles.size(); i++) {
            const auto &residue_circle = l_residue_circles[i];
            if (false)
               std::cout << "handling residue circle " << i << " " << residue_circle.residue_label
                         << std::endl;
	    lig_build::pos_t pos = residue_circle.pos;
            // get rid of add_rep_handles in the following function call
	    int add_rep_handle = -1; // default, no handle
	    if (add_rep_handles.size() == l_residue_circles.size())
	       add_rep_handle = add_rep_handles[i];

            svg_container_t svgc_s = draw_residue_circle_top_layer(l_residue_circles[i],
                                                                   ligand_centre, add_rep_handle);
            svgc.add(svgc_s);
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: draw_residue_circles: " << rte.what() << std::endl;
      }
   }
   return svgc;
}

#include "utils/coot-utils.hh"

// if you don't have an add_rep_handle, then pass -1 (something negative)
//
svg_container_t
flev_t::draw_residue_circle_top_layer(const residue_circle_t &residue_circle,
                                      const lig_build::pos_t &ligand_centre,
                                      int add_rep_handle) {

   svg_container_t svgc;

   if (false)
      std::cout << "   draw_residue_circle_top_layer() " << residue_circle.residue_type
                << " at init pos " << residue_circle.pos << " and canvas_drag_offset "
                << std::endl;

   lig_build::pos_t circle_pos(residue_circle.pos.x, -residue_circle.pos.y);
   lig_build::pos_t pos = circle_pos;

   // GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   //GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", "#111111", NULL);

   // don't draw waters that don't have bonds to the ligand (or the
   // bonds to the ligand atoms are too long, or the water is too far
   // from any protein atom).
   //
   if (residue_circle.residue_type == "HOH") {
      if (residue_circle.bonds_to_ligand.size() == 0) {
         return svgc;
      }
   }

   // Capitalise the residue type (takes less space than upper case).
   std::string rt = residue_circle.residue_type.substr(0,1);
   rt += coot::util::downcase(residue_circle.residue_type.substr(1));
   if (residue_circle.residue_type == "HOH") rt = "WATER";
   if (residue_circle.residue_type == "PO4") rt = "PO4";
   if (residue_circle.residue_type == "SO4") rt = "SO4";

   // correct that if we are looking at dna: DA, DT, DC, DG
   if (residue_circle.residue_type == "DA" ||
       residue_circle.residue_type == "DT" ||
       residue_circle.residue_type == "DC" ||
       residue_circle.residue_type == "DG") {
      rt = "d";
      rt += residue_circle.residue_type.substr(1);
   }

   // fill colour and stroke colour of the residue circle
   std::pair<std::string, std::string> col = get_residue_circle_colour(residue_circle.residue_type);
   double stroke_width = 0.1;
   if (col.second != "#111111") // needs checking, FIXME
      stroke_width = 0.15;

   if (col.first != "") {

      std::string circle_string = std::string("   <circle ") +
         std::string("cx=\"") + std::to_string(pos.x) + std::string("\" ") +
         std::string("cy=\"") + std::to_string(pos.y) + std::string("\" ") +
         // std::string("\" r=\"22.2\"") +
         std::string("r=\"1.2\" ") +
         std::string("fill=\"")   + col.first  + std::string("\" ") +
         std::string("stroke=\"") + col.second + std::string("\" ") +
         std::string("stroke-width=\"") + std::to_string(stroke_width) + std::string("\"") +
         std::string("/>\n");
      // expand the bounding box a bit beyond the position of the residue
      float delta = 3.7; // 20250303-PE was 3.0 // about this? // was 50.0
      svgc.add("<!-- Residue Circle " + residue_circle.residue_label + std::string(" -->\n"));
      svgc.add(circle_string);
      svgc.set_bounds(pos.x - delta, pos.y - delta,
                      pos.x + delta, pos.y + delta);

      // 20241002-PE note to self, "Phe" is too far to the right
      //                           "Ile" is too tar to the left
      //             I want to be able to specify the middle of the label, not
      //             bottom left.
      // 20241005-PE Ah, now it is. Better.
      double y_offset = 0.1;
      std::string font_size = "0.06";
      if (residue_circle.residue_type == "HOH") font_size = "0.035";
      if (residue_circle.residue_type == "PO4") font_size = "0.054";
      if (residue_circle.residue_type == "SO4") font_size = "0.054";
      if (residue_circle.residue_type == "HOH") y_offset = 0.18;
      // std::cout << "debug:: testing for residue_type " << residue_circle.residue_type << " font_size " << font_size << std::endl;
      std::string text_1("   <text ");
      text_1 += std::string("fill=\"#111111\"");
      text_1 += std::string(" x=\"");
      text_1 += std::to_string(pos.x);
      text_1 += std::string("\"");
      text_1 += std::string(" y=\"");
      text_1 += std::to_string(pos.y - y_offset);
      text_1 += std::string("\"");
      text_1 += std::string(" text-anchor=\"middle\"");
      text_1 += std::string(" font-family=\"Helvetica, sans-serif\" font-size=\"" + font_size + "em\">");
      text_1 += rt;
      text_1 += std::string("</text>\n");
      svgc.add(text_1);

      std::string text_2("   <text ");
      text_2 += std::string("fill=\"#111111\"");
      text_2 += std::string(" x=\"");
      text_2 += std::to_string(pos.x); // was -9
      text_2 += std::string("\"");
      text_2 += std::string(" y=\"");
      text_2 += std::to_string(pos.y + 0.6);
      text_2 += std::string("\"");
      text_2 += std::string(" text-anchor=\"middle\"");
      text_2 += std::string(" font-family=\"Helvetica, sans-serif\" font-size=\"0.04em\">");
      text_2 += residue_circle.residue_label;
      text_2 += std::string("</text>\n");
      svgc.add(text_2);

   } else {
      std::cout << "................ missing blank colour residue circle top" << std::endl;
   }
   return svgc;
}

// Set the fill_colour to "none" if you don't want a fill.
// Note that this function does not apply the inversion of the y axis coordinate.
std::string
pli::make_circle(const lig_build::pos_t &pos, double radius, double stroke_width,
                    const std::string &fill_color, const std::string &stroke_color) {

   std::string circle_string = std::string("   ") +
      "<circle cx=\"" + std::to_string(pos.x) + std::string("\" cy=\"") +
      std::to_string(pos.y) +
      std::string("\" r=\"") +
      std::to_string(radius) +
      std::string("\"");
      circle_string += std::string(" fill=\"")   + fill_color  + std::string("\"");
   circle_string +=
      std::string(" stroke=\"") + stroke_color + std::string("\"") +
      std::string(" stroke-width=\"") + std::to_string(stroke_width) + std::string("\"") +
      std::string("/>\n");
   float delta = 50.0; // about this?
   return circle_string;
}


// solvent exposure difference of the residue due to ligand binding
svg_container_t
flev_t::draw_solvent_exposure_circle(const residue_circle_t &residue_circle,
                                     const lig_build::pos_t &ligand_centre) {

   svg_container_t svgc;

   if (residue_circle.residue_type != "HOH") {
      if (residue_circle.se_diff_set()) {
	 std::pair<double, double> se_pair = residue_circle.solvent_exposures();
	 double radius_extra = (se_pair.second - se_pair.first) * 1.2;  // was 19, was 18, was 14, was 22.
	 if (radius_extra > 0) {
	    lig_build::pos_t to_lig_centre_uv = (ligand_centre - residue_circle.pos).unit_vector();
	    lig_build::pos_t se_circle_centre = residue_circle.pos - to_lig_centre_uv * radius_extra;

	    std::string fill_colour = get_residue_solvent_exposure_fill_colour(radius_extra);
	    double r = standard_residue_circle_radius + radius_extra;
            if (false)
               std::cout << "in draw_solvent_exposure_circle() to_lig_centre_uv " << to_lig_centre_uv << " "
                         << "se_circle_centre " << se_circle_centre << " fill_colour " << fill_colour << " "
                         << "radius_extra " << radius_extra << " "
                         << "r " << r << std::endl;
            double line_width = 0.0;
            svgc.add("<!-- Exposure Circle -->\n");
            lig_build::pos_t pos = se_circle_centre;
            pos.y = -pos.y; // to match the top layer
            pos += lig_build::pos_t(0.0002, 0.0002);
            if (false)
               std::cout << "   pos " << pos << " se_circle_centre " << se_circle_centre << std::endl;
            std::string c = pli::make_circle(pos, r, line_width, fill_colour, "black");
            svgc.add(c);
         }
      }
   }
   return svgc;
}

std::string
flev_t::get_residue_solvent_exposure_fill_colour(double r) const {

   // r is ~0.4 max
   std::string colour = "#8080ff";
   double step = 0.04;
   if (r > step)
      colour = "#e0e0ff";
   if (r > step * 2)
      colour = "#d8d8ff";
   if (r > step * 3)
      colour = "#d0d0ff";
   if (r > step * 4)
      colour = "#c0c8ff";
   if (r > step * 5)
      colour = "#b0c0ff";
   if (r > step * 6)
      colour = "#a0b8ff";
   if (r > step * 7)
      colour = "#90b0ff";
   if (r > step * 8)
      colour = "#80a8ff";
   if (r > step * 9)
      colour = "#70a0ff";

   return colour;
}


// untrap residues as needed.
void
flev_t::position_non_primaries(const ligand_grid &grid,
                               const std::vector<int> &primary_indices) {

   std::vector<grid_index_t> already_positioned;

   for (unsigned int irc=0; irc<residue_circles.size(); irc++) {
      // if not a primary...
      if (std::find(primary_indices.begin(), primary_indices.end(), irc) ==
          primary_indices.end()) {

         // if this point has rejection underneath it, find the
         // position in grid that doesn't have any rejection.
         //

         std::pair<grid_index_t, lig_build::pos_t> pos =
            grid.find_nearest_zero(residue_circles[irc].pos, already_positioned);
         residue_circles[irc].pos = pos.second;
         if (pos.first.is_valid_p()) {
            std::cout << "position_non_primaries() " << irc << " " << pos.first.i() << " " << pos.first.j() << std::endl;
            already_positioned.push_back(pos.first);
         }
      }
   }
}

// can throw a std::runtime_error.
//
flev_t::grid_index_t
flev_t::ligand_grid::grid_pos_nearest(const lig_build::pos_t &pos) const {

   lig_build::pos_t p = pos - lig_build::pos_t(ligand_atoms_min_x, ligand_atoms_min_y);
   int idx_x = int(p.x/scale_fac+0.5);
   int idx_y = int(p.y/scale_fac+0.5);

   if ((idx_x < 0) || (idx_x >= x_size()) || (idx_y < 0) || (idx_y >= y_size()))
       throw std::runtime_error("out of grid index");

   return grid_index_t(idx_x, idx_y);
}


// can throw an exception
lig_build::pos_t
flev_t::ligand_grid::find_minimum_position() const {

   double best_pos_score = 1000000;
   lig_build::pos_t best_pos;
   for (int ix=0; ix<x_size(); ix++) {
      for (int iy=0; iy<y_size(); iy++) {
         if (grid_[ix][iy] < best_pos_score) {
            best_pos_score = grid_[ix][iy];
            best_pos = grid_pos_to_mol_space_pos(ix,iy);
         }
      }
   }
   if (best_pos_score > (1000000-1))
      throw std::runtime_error("failed to get minimum position from ligand grid");
   return best_pos;
}


// actually, not exactly zero but something small.
//
// Don't return a grid-point/position that matches anything in
// already_positioned.
//
std::pair<flev_t::grid_index_t, lig_build::pos_t>
flev_t::ligand_grid::find_nearest_zero(const lig_build::pos_t &pos,
                                       const std::vector<flev_t::grid_index_t> &already_positioned) const {

   lig_build::pos_t p;
   grid_index_t rgi; // initially invalid
   double shortest_dist = 43e23;
   double crit = 0.05; // less than this is effectively zero.

   try {
      grid_index_t gi=grid_pos_nearest(pos);
      if (grid_[gi.i()][gi.j()] < crit) {
         p = pos; // fine, no change
         std::cout << "fine, no change " << p << std::endl;
      } else {
         std::cout << "search for some place else, grid lims: " << x_size() << " " << y_size() << std::endl;
         // search for someplace else
         for (int ix=0; ix<x_size(); ix++) {
            for (int iy=0; iy<y_size(); iy++) {
               if (false)
                  std::cout << "grid_value " << ix << " " << iy << " is " << grid_[ix][iy] << std::endl;
               if (grid_[ix][iy] < crit) {
                  lig_build::pos_t gp = grid_pos_to_mol_space_pos(ix, iy);
                  if (false)
                     std::cout << "   ix " << ix << " iy " << iy << " gp " << gp << " c.f. pos " << pos << std::endl;
                  double d = (gp - pos).lengthsq();
                  if (d < shortest_dist) {
                     grid_index_t candidate_grid_index(ix, iy);
                     // This is OK if there is no other previous
                     // solution at the same position (if there is,
                     // keep trying, of course).
                     if (std::find(already_positioned.begin(), already_positioned.end(),
                                   candidate_grid_index) == already_positioned.end()) {
                        shortest_dist = d;
                        p = gp;
                        rgi = candidate_grid_index;
                        std::cout << "found some place else: " << p << std::endl;
                     }
                  }
               }
            }
         }
      }
   }
   catch (const std::runtime_error &rte) {
      // the pos was off the grid.  It won't be trapped inside the
      // ligand, so just return what we were given.
      p = pos;
   }
   return std::pair<grid_index_t, lig_build::pos_t> (rgi, p);
}





void
flev_t::initial_primary_residue_circles_layout(const ligand_grid &grid,
                                               int primary_index,
                                               const std::vector<std::pair<lig_build::pos_t, double> > &attachment_points) {

   if (false)
      std::cout << "DEBUG:: starting initial_primary_residue_circles_layout() primary_index " << primary_index
                << " " << residue_circles[primary_index].residue_label << " "
                << residue_circles[primary_index].residue_type
                << " has position " << residue_circles[primary_index].pos
                << std::endl;

   if (false)
      std::cout << "DEBUG:: initial_primary_residue_circles_layout() =========== adding quadratic for residue "
                << residue_circles[primary_index].residue_label
                << " ============================"
                << std::endl;

   if (false)
      grid.print(primary_index);

   ligand_grid primary_grid = grid;

   // attachment points are points on the ligand, in ligand
   // coordinates to which this primary residue circle is
   // attached (often only one attachment, but can be 2
   // sometimes).
   primary_grid.add_quadratic(attachment_points);

   // if (true)
   //    if (primary_index == 0)
   //       show_grid(grid);

   lig_build::pos_t best_pos = primary_grid.find_minimum_position();

   // OK, consider the case where there are 2 residues bonding to the
   // same atom in the ligand.  They wil both be given exactly the
   // same best_pos and they will subsequently refine together,
   // resulting in one residue sitting on top of another - bad news.
   //
   // So let's shift the residue a bit in the direction that it came
   // from so that the residues don't refine together.
   //
   lig_build::pos_t a_to_b_uv = (residue_circles[primary_index].pos - best_pos).unit_vector();

   residue_circles[primary_index].pos = best_pos + a_to_b_uv * 4;

   if (false)
      std::cout << "DEBUG::  ending initial_primary_residue_circles_layout() primary_index "
                << primary_index
                << " " << residue_circles[primary_index].residue_label << " "
                << residue_circles[primary_index].residue_type
                << " has position " << residue_circles[primary_index].pos
                << std::endl;

}

// attachment points are points on the ligand, in ligand coordinates
// to which this primary residue circle is attached (often only one
// attachment, but can be 2 sometimes).
//
void
flev_t::ligand_grid::add_quadratic(const std::vector<std::pair<lig_build::pos_t, double> > &attachment_points) {

   if (attachment_points.size()) {
      double scale_by_n_attach = 1.0/double(attachment_points.size());

      for (unsigned int iattach=0; iattach<attachment_points.size(); iattach++) {
         for (int ix=0; ix<x_size(); ix++) {
            for (int iy=0; iy<y_size(); iy++) {
               lig_build::pos_t pos = grid_pos_to_mol_space_pos(ix, iy);
               double d2 = (pos-attachment_points[iattach].first).lengthsq();
               double val = 0.2 * d2 * scale_by_n_attach; // test scaling here FIXME
               grid_[ix][iy] += val;
            }
         }
      }
   }
}

// minimise layout energy
std::pair<int, std::vector<residue_circle_t> >
flev_t::optimise_residue_circle_positions(const std::vector<residue_circle_t> &r,
                                          const std::vector<residue_circle_t> &c,
                                          const std::vector<int> &primary_indices) const {
   if (r.size() > 0) {
      if (c.size() == r.size()) {
         pli::optimise_residue_circles orc(r, c, mol, primary_indices);
         int status = orc.get_gsl_min_status();
         std::cout << "debug:: in optimise_residue_circles() get_gsl_min_status() status "
                   << status << std::endl;
         return orc.solution();
      }
   }
   std::vector<residue_circle_t> dv; // dummy
   return std::pair<int, std::vector<residue_circle_t> > (0, dv);
}

#endif // MAKE_ENHANCED_LIGAND_TOOLS
