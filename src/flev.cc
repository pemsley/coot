/* src/flev.cc
 * 
 * Copyright 2010, 2011, 2012 The University of Oxford
 * Copyright 2014 by Medical Research Council
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street,  02110-1301, USA
 */

/*  ----------------------------------------------------------------------- */
/*               Flattened Ligand Environment View  Interface               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


// For reasons I don't understand, this should come near the top of
// includes, otherwise we get RDKit dcgettext() include file problems.
//
// #include "lbg/lbg.hh"

#include "c-interface-ligands.hh"
#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"

#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "geometry/lbg-graph.hh"

#include "coot-utils/coot-h-bonds.hh"
#include "pli/protein-ligand-interactions.hh"

// 20140311: also, we want to include pi-stacking explicitly here
// because we need it to compile flev if we don't HAVE_GOOCANVAS
// (pi-stacking.hh is included from lbg.hh if HAVE_GOOCANVAS).
// 
#include "pli/pi-stacking.hh"

#include "flev.hh"

#include "cc-interface.hh" // for add_animated_ligand_interactions


int
sprout_hydrogens(int imol,
		 const char *chain_id,
		 int res_no,
		 const char *ins_code) {
   
   int r = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::pair<bool, std::string> r_add =
	 g.molecules[imol].sprout_hydrogens(chain_id, res_no, ins_code, *(g.Geom_p()));
      r = r_add.first;
      if (r)
	 graphics_draw();
      else
	 info_dialog(r_add.second.c_str());
   }
   return r;
}

//  return true if both ligand_res and env_res have hydrogens in the
//  dictionary and (at least one) hydrogen in the model
// (if ligand_res is a SO4 (no hydrogens) that should not cause a fail).
// 
bool flev_check_for_hydrogens(int imol, mmdb::Residue *ligand_res,
			      const std::vector<mmdb::Residue *> &env_residues,
			      coot::protein_geometry *geom_p) {

   bool status = false;
   int n_hydrogens_ligand = 0;
   int n_hydrogens_others = 0;

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   ligand_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) { 
      std::string ele = residue_atoms[i]->element;
      if (ele == " H") { // needs PDBv3 update
	 n_hydrogens_ligand++;
      } 
   }
   for (unsigned int ires=0; ires<env_residues.size(); ires++) { 
      residue_atoms = 0;
      n_residue_atoms = 0;
      env_residues[ires]->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) { 
	 std::string ele = residue_atoms[iat]->element;
	 if (ele == " H") { // needs PDBv3 update
	    n_hydrogens_others++;
	 }
      }
   }
   if (n_hydrogens_others > 0)
      if (n_hydrogens_ligand > 0)
	 return true;


   if (n_hydrogens_ligand == 0) {
      // perhaps the ligand was not supposed to have any ligands
      std::string residue_type = ligand_res->GetResName();
      if (geom_p->have_dictionary_for_residue_type(residue_type, imol,
						   graphics_info_t::cif_dictionary_read_number++)) {
	 int n_hydrogens_ligand_dict = geom_p->n_hydrogens(residue_type);
	 if (n_hydrogens_ligand_dict > 0)
	    status = false;
	 else
	    if (n_hydrogens_others > 0)
	       // environment has hydrogens, ligand does not (and is
	       // not supposed to have)
	       status = true;
	    else
	       // environment does not have hydrogens (and is supposed
	       // to have (we presume))
	       status = false;
      }
   }
   return status;
} 


void fle_view_with_rdkit(int imol, const char *chain_id, int res_no,
			 const char *ins_code, float residues_near_radius) {

   fle_view_with_rdkit_internal(imol, chain_id, res_no, ins_code, residues_near_radius, "", "");
}

void fle_view_with_rdkit_to_png(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *png_file_name) {

   fle_view_with_rdkit_internal(imol, chain_id, res_no, ins_code, residues_near_radius, "png", png_file_name);
}

void fle_view_with_rdkit_to_svg(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *png_file_name) {

   fle_view_with_rdkit_internal(imol, chain_id, res_no, ins_code, residues_near_radius, "svg", png_file_name);
} 


void fle_view_with_rdkit_internal(int imol, const char *chain_id, int res_no, const char *ins_code, float residues_near_radius, const char *file_format, const char *output_image_file_name) {
  // function removed
}

void fle_view_set_water_dist_max(float dist_max) { 
   graphics_info_t::fle_water_dist_max = dist_max;   // default 3.25
} 

void fle_view_set_h_bond_dist_max(float h_bond_dist_max) {
   graphics_info_t::fle_h_bond_dist_max = h_bond_dist_max;   // default 3.9
} 


std::vector<int>
coot::make_add_reps_for_near_residues(std::vector<mmdb::Residue *> filtered_residues,
				      int imol) {

   std::vector<int> v(filtered_residues.size());

   int rep_type = coot::SIMPLE_LINES;
   int bonds_box_type = coot::NORMAL_BONDS;
   int bond_width = 8;
   int draw_hydrogens_flag = 1;
   for (unsigned int i=0; i<filtered_residues.size(); i++) {
      v[i] = additional_representation_by_attributes(imol,
						     filtered_residues[i]->GetChainID(),
						     filtered_residues[i]->GetSeqNum(),
						     filtered_residues[i]->GetSeqNum(),
						     filtered_residues[i]->GetInsCode(),
						     rep_type, bonds_box_type, bond_width,
						     draw_hydrogens_flag);
      set_show_additional_representation(imol, v[i], 0); // undisplay it
   }
   return v;
}

void
coot::add_animated_ligand_interactions(int imol,
				       const std::vector<pli::fle_ligand_bond_t> &ligand_bonds) {

   for (unsigned int i=0; i<ligand_bonds.size(); i++) {
      // std::cout << "Here....  adding animated ligand interaction " << i << std::endl;
      add_animated_ligand_interaction(imol, ligand_bonds[i]);
   }
}


std::vector<pli::fle_residues_helper_t>
coot::get_flev_residue_centres(mmdb::Residue *residue_ligand_3d,
			       mmdb::Manager *mol_containing_residue_ligand, 
			       std::vector<mmdb::Residue *> residues,
			       mmdb::Manager *flat_mol) {

   std::vector<pli::fle_residues_helper_t> centres;

   if (flat_mol) { 

      // get the lsq matrix that maps the ligand in 3D onto the flat ligand
      int res_no = residue_ligand_3d->GetSeqNum();
      std::string chain_id = residue_ligand_3d->GetChainID();
      int every_nth = 1;
      std::vector<coot::lsq_range_match_info_t> matches;
      coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id,
 					 coot::lsq_t::ALL);
       matches.push_back(match);
       std::pair<short int, clipper::RTop_orth> lsq_mat =
	  coot::util::get_lsq_matrix(flat_mol, mol_containing_residue_ligand, matches, every_nth);
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
	       pli::fle_residues_helper_t fle_centre(c.second,
                                                     coot::residue_spec_t(residues[ires]),
                                                     res_name);

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
   if (0) 
      for (unsigned int ic=0; ic<centres.size(); ic++) { 
	 std::cout << "centre " << ic << " has TR_centre "
		   << centres[ic].transformed_relative_centre.format() << std::endl;
      }
   return centres;
} 


std::map<std::string, std::string>
coot::make_flat_ligand_name_map(mmdb::Residue *flat_res) {

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


 
#ifdef VERY_OLD_FLEV_FUNCTIONS		     

void
coot::write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
			const coot::pi_stacking_container_t &stack_info,
			mmdb::Residue *res_flat) {

   if (0) {
      std::cout << "DEBUG:: in write_fle_centres() " << stack_info.stackings.size()
		<< " stackings:\n";
      for (unsigned int i=0; i<stack_info.stackings.size(); i++)
	 std::cout << "    " << i << ":  " << stack_info.stackings[i] << std::endl;
      std::cout << std::endl;
   }
   
   std::string file_name = "coot-tmp-fle-view-residue-info.txt";
   std::ofstream of(file_name.c_str());
   if (!of) {
      std::cout << "failed to open output file " << file_name << std::endl;
   } else {
      
      for (unsigned int i=0; i<v.size(); i++) {
	 if (v[i].is_set) { 
	    of << "RES " << v[i].centre.x() << " " << v[i].centre.y() << " "
	       << v[i].centre.z() << " "
	       << v[i].residue_name << " " << v[i].spec.chain << v[i].spec.res_no;

	    // Add the qualifying distance to protein if this is a water atom.
	    // 
	    if (v[i].residue_name == "HOH") { 
	       for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
		  if (bonds_to_ligand[ib].res_spec == v[i].spec) {
		     of << " water->protein-length " << bonds_to_ligand[ib].water_protein_length;
		     break;
		  }
	       }
	       of << "\n";
	    } else { 
	       of << "\n";
	    }


	    // now, was there a bond with this res spec?
	    // 
	    // if so find the x y coordinates of the atom with this
	    // name in the res_flat and write them out.
	    // 
	    for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	       if (bonds_to_ligand[ib].res_spec == v[i].spec) {
		  mmdb::PPAtom residue_atoms = 0;
		  int n_residue_atoms;
		  res_flat->GetAtomTable(residue_atoms, n_residue_atoms);
		  for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
		     mmdb::Atom *at = residue_atoms[iat];
		     std::string atom_name(at->name);
		     if (atom_name == bonds_to_ligand[ib].ligand_atom_spec.atom_name) {
			// the coordinates of the matching flat ligand
			// atom (done like solvent accessibilites).
			of << "BOND " << at->name
			   << " length: " << bonds_to_ligand[ib].bond_length
			   << " type: " << bonds_to_ligand[ib].bond_type
			   << "\n";
		     }
		  }
	       }
	    }


	    // Is there a solvent accessibility associated with this residue?
	    // 
	    for (unsigned int ised=0; ised<sed.size(); ised++) { 
	       if (sed[ised].res_spec == v[i].spec) {
		  of << "SOLVENT_ACCESSIBILITY_DIFFERENCE "
		     << sed[ised].exposure_fraction_holo << " " 
		     << sed[ised].exposure_fraction_apo << "\n";
	       }
	    }

	    // Is there a stacking interaction associated with this residue?
	    //
	    // There could, in fact, be several stackings associated
	    // with this residue.  We want to pick the best one.  So
	    // we need to run through the list of stacking first,
	    // looking for the best.
	    //
	    // 
	    double best_stacking_score = -1.0;
	    coot::pi_stacking_instance_t best_stacking(NULL, "");
	    
	    for (unsigned int istack=0; istack<stack_info.stackings.size(); istack++) {
	       coot::residue_spec_t spec(stack_info.stackings[istack].res);
	       if (spec == v[i].spec) {
		  if (stack_info.stackings[istack].overlap_score > best_stacking_score) {
		     best_stacking_score = stack_info.stackings[istack].overlap_score;
		     best_stacking = stack_info.stackings[istack];
		  }
	       }
	    }

	    if (best_stacking_score > 0.0) {
	       std::string type = "    pi-pi"; // leading spaces
					       // important in format
					       // of atom names

	       // Recall that CATION_PI_STACKING is for ligand ring
	       // systems, PI_CATION_STACKING is for protein residue
	       // ring systems.
	       // 
	       if (best_stacking.type == coot::pi_stacking_instance_t::CATION_PI_STACKING)
		  type = "cation-pi";
	       if (best_stacking.type == coot::pi_stacking_instance_t::PI_CATION_STACKING)
		  type = "pi-cation";
	       if (0) // debug
		  of << "# Best stacking for RES " 
		     << v[i].residue_name << " " << v[i].spec.chain_id
		     << v[i].spec.res_no << " is " << best_stacking.type << " "
		     << best_stacking_score << std::endl;
	       of << "STACKING " << type << " ";
	       switch (best_stacking.type) {
	       case coot::pi_stacking_instance_t::NO_STACKING:
		  break;
	       case coot::pi_stacking_instance_t::PI_PI_STACKING: 
		  for (unsigned int jat=0;
		       jat<best_stacking.ligand_ring_atom_names.size();
		       jat++) {
		     of << best_stacking.ligand_ring_atom_names[jat] << "  ";
		  }
		  break;
		  
	       case coot::pi_stacking_instance_t::CATION_PI_STACKING:
		  of << best_stacking.ligand_cationic_atom_name;
		  break;
		  
	       case coot::pi_stacking_instance_t::PI_CATION_STACKING:
		  for (unsigned int jat=0;
		       jat<best_stacking.ligand_ring_atom_names.size();
		       jat++) {
		     of << best_stacking.ligand_ring_atom_names[jat] << "  ";
		  }
		  break;
	       }
	       of << "\n";
	    }
	    
	 } else {
	    std::cout << "WARNING centre " << i << " was not set " << std::endl;
	 } 
      }
   }
}

#endif //  VERY_OLD_FLEV_FUNCTIONS		     


// the reference_residue is the flat residue
void
coot::write_ligand_atom_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
					const pli::flev_attached_hydrogens_t &attached_hydrogens,
					mmdb::Residue *reference_residue) {

   // for each of the atoms in reference_residue, find the spec in sav
   // and write out the coordinates and the atom name and the
   // solvent_accessibility.

   std::string file_name = "coot-tmp-fle-view-solvent-accessibilites.txt";
   std::ofstream of(file_name.c_str());
   if (!of) {
      std::cout << "failed to open output file " << file_name << std::endl;
   } else {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      reference_residue->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 mmdb::Atom *at = residue_atoms[i];
	 std::string atom_name(at->GetAtomName());
	 for (unsigned int j=0; j<sav.size(); j++) {
	    if (atom_name == sav[j].first.atom_name) {
	       of << "ATOM:" << atom_name << " " << at->x << " " <<  at->y << " " << at->z
		  << " " << sav[j].second << "\n";
	       std::map<std::string, std::vector<coot::bash_distance_t> >::const_iterator it =
		  attached_hydrogens.atom_bashes.find(atom_name);
	       if (it != attached_hydrogens.atom_bashes.end()) {
		  for (unsigned int ibash=0; ibash<it->second.size(); ibash++) { 
		     of << "   BASH: " << it->second[ibash] << "\n";
		  }
	       } 
	       break;
	    }
	 }
      }
   }
}



bool
coot::standard_residue_name_p(const std::string &rn) {

   std::vector<std::string> s;

   s.push_back("ALA");
   s.push_back("ASP");
   s.push_back("ASN");
   s.push_back("CYS");
   s.push_back("GLN");
   s.push_back("GLY");
   s.push_back("GLU");
   s.push_back("PHE");
   s.push_back("HIS");
   s.push_back("ILE");
   s.push_back("LYS");
   s.push_back("LEU");
   s.push_back("MET");
   s.push_back("MSE");
   s.push_back("PRO");
   s.push_back("ARG");
   s.push_back("SER");
   s.push_back("THR");
   s.push_back("VAL");
   s.push_back("TRP");
   s.push_back("TYR");

   s.push_back("AR");
   s.push_back("AD");
   s.push_back("CR");
   s.push_back("CD");
   s.push_back("GR");
   s.push_back("GD");
   s.push_back("TD");

   if (std::find(s.begin(), s.end(), rn) == s.end())
      return 0;
   else
      return 1;
}

#include "pli/flev.hh"

std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > >
coot::get_cannonball_vectors(mmdb::Residue *ligand_res_3d,
			     const coot::dictionary_residue_restraints_t &monomer_restraints) {

   std::vector<std::pair<mmdb::Atom *, std::vector<clipper::Coord_orth> > > v;

   return v;
}

// put fle_view into a C++ header.
void fle_view(int imol, const char *chain_id, int res_no, const char *ins_code, float dist_max) {

   bool add_key = true;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {
         mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
         if (residue_p) {
            svg_container_t svgc = pli::fle_view_with_rdkit_internal(mol, imol, graphics_info_t::Geom_p(),
                                                                     chain_id, res_no, ins_code, dist_max, add_key);
            std::string s = svgc.compose(true);
            if (! s.empty()) {
               std::string rn = residue_p->GetResName();
               std::string title = "Coot: 2D Ligand Environment View for ";
               title += std::string(chain_id);
               title += std::string(" ");
               title += std::to_string(res_no);
               title += std::string(" ");
               title += rn;
               display_svg_from_string_in_a_dialog(s, title);
            } else {
               std::cout << "ERROR:: failed to make depiction - empty svg" << std::endl;
            }
         }
      }
   }
}
