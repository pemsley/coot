
#include "lidia-core/rdkit-interface.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "flev.hh"
#include "flev-annotations.hh"
#include "dots-representation-info.hh"
#include "protein-ligand-interactions.hh"
#include "atom-ring-centre-info.hh"


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
      coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id,
                                          COOT_LSQ_ALL);
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
   if (true)
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


   if (true) {
      for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) {
         std::cout << "  =============== lbg::annotate() bond to ligand " << ib << " "
                   << bonds_to_ligand[ib].ligand_atom_spec.atom_name
                   << " by residue " << bonds_to_ligand[ib].res_spec.chain_id << " "
                   << bonds_to_ligand[ib].res_spec.res_no << " type: "
                   << bonds_to_ligand[ib].bond_type
                   << std::endl;
      }
   }

   if (true) {
      std::cout << "--------------------------------------------------------------" << std::endl;
      std::cout << "======== lbg_info_t::annotate() here are bash distances for atoms:" << std::endl;
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

   lig_build::molecule_t<svg_atom_t, svg_bond_t> new_mol = mol;
   // modify the atoms of new_mol
   map_solvent_accessibilities_to_atoms(&new_mol, solvent_accessible_atoms);

   // fill class data item std::vector<residue_circle_t> residue_circles;
   //
   residue_circles.clear();
   for (unsigned int i=0; i<centres.size(); i++) {

      if (true)
         std::cout << "debugg:: in lbg_info_t::annotate() handling circle " << i << " of "
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
   if (true) {
      std::cout << "------------------- ring list ------------" << std::endl;
      for (unsigned int i=0; i<ring_atoms_list.size(); i++) {
         std::cout << "ring list " << i << "   ";
         for (unsigned int j=0; j<ring_atoms_list[i].size(); j++) {
         std::cout << ring_atoms_list[i][j] << "  ";
         }
         std::cout << std::endl;
      }
   }

   render();

#ifdef COMPILE_THE_HEART

   import_from_widgeted_molecule(new_mol);
   render();
   refine_residue_circle_positions();

   // has the current solution problems due to residues too close to the ligand?
   std::pair<bool, std::vector<int> > problem_status = solution_has_problems_p();
   // std::cout << "::::::::: problem status: " << problem_status.first << std::endl;
   for (unsigned int ip=0; ip<problem_status.second.size(); ip++) {
      std::cout << ":::::::::::: "
                << residue_circles[problem_status.second[ip]].residue_label << " "
                << residue_circles[problem_status.second[ip]].residue_type << " "
                << residue_circles[problem_status.second[ip]].pos
                << std::endl;
   }

   if (problem_status.second.size()) {
      // fiddle with residue_circles and reoptimise.
      std::vector<int> primary_indices = get_primary_indices();
      reposition_problematics_and_reoptimise(problem_status.second, primary_indices);
   }

   recentre_considering_residue_centres();

#endif // COMPILE_THE_HEART

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



void
pli::fle_view_with_rdkit_internal(mmdb::Manager *mol,
                                  int imol,
                                  coot::protein_geometry *geom_p,
                                  const std::string &chain_id, int res_no, const std::string &ins_code,
                                  float residues_near_radius,
                                  const std::string &file_format, const std::string &output_image_file_name) {

#ifndef MAKE_ENHANCED_LIGAND_TOOLS

   std::cout << "WARNING:: fle_view_with_rdkit_internal() not enhanced ligand " << std::endl;

# else

   double weight_for_3d_distances = 0.4; // for 3d distances
   std::string output_format = file_format;

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

            // residues needs to be filtered to remove waters that
            // are not connected to a protein atom.

            std::vector<mmdb::Residue *> filtered_residues =
               coot::filter_residues_by_solvent_contact(res_ref, mol_for_res_ref, residues, 3.25);

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


                  int mol_2d_depict_conformer =
                     coot::add_2d_conformer(&rdkm, weight_for_3d_distances);
                  lig_build::molfile_molecule_t m =
                     coot::make_molfile_molecule(rdkm, mol_2d_depict_conformer);

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

                  flev_t flev;

                  // should this be a flev function?
                   std::vector<pli::fle_ligand_bond_t> bonds_to_ligand =
                      pli::get_fle_ligand_bonds(res_ref, filtered_residues,
                                                mol_for_res_ref, name_map, *geom_p,
                                                flev.fle_water_dist_max,
                                                flev.fle_h_bond_dist_max);

                  std::vector<fle_residues_helper_t> res_centres =
                     pli::get_flev_residue_centres(res_ref,
                                                   mol_for_res_ref,
                                                   filtered_residues,
                                                   mol_for_flat_residue);

                  std::vector<int> add_reps_vec;

                  if (true) {
                     for (unsigned int ic=0; ic<res_centres.size(); ic++)
                        std::cout << "   " << ic << "  " << res_centres[ic]
                                  << std::endl;
                  }

                  // ----------- residue infos ----------
                  //
                  pli::pi_stacking_container_t pi_stack_info(p.second, filtered_residues, res_ref, rdkm);

                  // ----------- ligand atom infos ------
                  //
                   flev_attached_hydrogens_t ah(p.second);
                  // ah.cannonballs(res_ref, mol_for_res_ref, p.second);
                   ah.distances_to_protein_using_correct_Hs(res_ref, mol_for_res_ref, *geom_p);

                  // ------------ show it -----------------
                  //
                   flev.annotate(s_a_v, res_centres, add_reps_vec, bonds_to_ligand, sed, ah, pi_stack_info, p.second);

                  flev.draw_all_flev_annotations();

                  if (output_format == "png") flev.write_png(output_image_file_name);
                  if (output_format == "svg") flev.write_svg(output_image_file_name);

                  delete mol_for_flat_residue;
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
#endif // MAKE_ENHANCED_LIGAND_TOOLS
}

void
flev_t::draw_substitution_contour() {

   bool debug = true;

   bool draw_flev_annotations_flag = true; // why wouldn't it be?

   if (draw_flev_annotations_flag) {
      if (mol.atoms.size() > 0) {

         // first of all, do we have any bash distances for the atoms of this molecule?
         bool have_bash_distances = 0;
         for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
            if (mol.atoms[iat].bash_distances.size()) {
               have_bash_distances = 1;
               break;
            }
         }

         // If we don't have bash distances, then don't grid and contour anything.  If
         // we do, we do....
         //
         if (have_bash_distances) {

            // REPLACE-ME-WITH-SVG
            // GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

            try {

               std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair =
                  mol.ligand_extents();

               ligand_grid grid(l_e_pair.first, l_e_pair.second);

               if (true) { // debug
                  for (unsigned int i=0; i<mol.atoms.size(); i++) {
                     std::cout << "in draw_substitution_contour() atom " << i << " "
                               << mol.atoms[i].get_atom_name()
                               << " has "
                               << mol.atoms[i].bash_distances.size() << " bash distances"
                               << std::endl;
                     for (unsigned int j=0; j<mol.atoms[i].bash_distances.size(); j++) {
                        std::cout << "  " << mol.atoms[i].bash_distances[j];
                     }
                     if (mol.atoms[i].bash_distances.size())
                        std::cout << std::endl;
                  }
               }


               std::vector<lig_build::atom_ring_centre_info_t> unlimited_atoms;
               // std::vector<widgeted_atom_ring_centre_info_t> unlimited_atoms;

               for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
                  int n_bash_distances = 0;
                  double sum_bash = 0.0;
                  bool unlimited = 0;
                  int n_unlimited = 0;

                  if (mol.atoms[iat].bash_distances.size()) {
                     for (unsigned int j=0; j<mol.atoms[iat].bash_distances.size(); j++) {
                        if (! mol.atoms[iat].bash_distances[j].unlimited()) {
                           sum_bash += mol.atoms[iat].bash_distances[j].dist;
                           n_bash_distances++;
                        } else {
                           unlimited = 1;
                           n_unlimited++;
                        }
                     }

                     if (! unlimited) {
                        if (n_bash_distances > 0) {
                           double bash_av = sum_bash/double(n_bash_distances);
                           if (debug)
                              std::cout << "   none unlimited, using bash_av " << bash_av
                                        << " for atom " << mol.atoms[iat].get_atom_name()
                                        << std::endl;
                           grid.add_for_accessibility(bash_av, mol.atoms[iat].atom_position);
                        }
                     } else {

                        // Now, were more than half of the bash distances
                        // unlimited?  If yes, then this is unlimited.

                        if (double(n_unlimited)/double(mol.atoms[iat].bash_distances.size()) > 0.5) {

                           // just shove some value in to make the grid values vaguely
                           // correct - the unlimited atoms properly assert themselves
                           // in the drawing of the contour (that is, the selection of
                           // the contour fragments).
                           //
                           grid.add_for_accessibility(1.2, mol.atoms[iat].atom_position);
                           //                std::cout << "adding unlimited_atom_position "
                           //                          << iat << " "
                           //                          << mol.atoms[iat].atom_position

                           // not elegant because no constructor for
                           // widgeted_atom_ring_centre_info_t (because no simple
                           // constructor for lig_build::atom_t).
                           //
                           lig_build::atom_ring_centre_info_t ua(mol.atoms[iat]);
                           unlimited_atoms.push_back(ua);
                        } else {

                           // treat as a limited:
                           double bash_av =
                              (sum_bash + 4.0 * n_unlimited)/double(n_bash_distances+n_unlimited);
                           if (debug)
                              std::cout << "   few unlimited, as limited using bash_av "
                                        << bash_av << " for atom "
                                        << mol.atoms[iat].get_atom_name()
                                        << std::endl;

                           grid.add_for_accessibility(bash_av, mol.atoms[iat].atom_position);
                        }
                     }

                  } else {

                     // we don't get here currently, now that there is an
                     // outer test for having bash distances. Current way
                     // is OK, I think.

                     // there were no bash distancs - what do we do?  Leaving out
                     // atoms means gaps over the ligand - bleugh.  Shove some
                     // value in?  1.0?  if they are not hydrogens, of course.
                     //
                     if (mol.atoms[iat].element != "H") // checked.
                        grid.add_for_accessibility_no_bash_dist_atom(1.0, mol.atoms[iat].atom_position);
                  }
               }

               // Put some values around the ring centres too.
               //
               grid.avoid_ring_centres(ring_atoms_list, mol);

               // for debugging
               // show_grid(grid);

               // std::vector<widgeted_atom_ring_centre_info_t> dummy_unlimited_atoms;
               grid.show_contour(0.5, unlimited_atoms, ring_atoms_list);
               // debug
               // show_unlimited_atoms(unlimited_atoms);
               // show_ring_centres(ring_atoms_list, mol);

               if (debug) {
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
   }
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

void
flev_t::draw_all_flev_annotations() {

   draw_all_flev_residue_attribs();
   draw_all_flev_ligand_annotations();
}

void
flev_t::draw_all_flev_residue_attribs() {

   draw_residue_circles(residue_circles, additional_representation_handles);
   draw_bonds_to_ligand();
   draw_stacking_interactions(residue_circles);
}

void
flev_t::draw_all_flev_ligand_annotations() {

   draw_substitution_contour();
   draw_solvent_accessibility_of_atoms();
}

void
flev_t::draw_stacking_interactions(const std::vector<residue_circle_t> &rc) {

   for (unsigned int ires=0; ires<rc.size(); ires++) {
      int st = rc[ires].get_stacking_type();
      clipper::Coord_orth click_pos = rc[ires].residue_centre_real;

      if (rc[ires].has_ring_stacking_interaction()) {

         std::vector<std::string> ligand_ring_atom_names =
            rc[ires].get_ligand_ring_atom_names();
         if ((st == residue_circle_t::PI_PI_STACKING) ||
             (st == residue_circle_t::PI_CATION_STACKING)) {
            try {
               lig_build::pos_t lc = mol.get_ring_centre(ligand_ring_atom_names);
               draw_annotated_stacking_line(lc, rc[ires].pos, st, click_pos);
            }
            catch (const std::runtime_error &rte) {
               std::cout << rte.what() << std::endl;
            }
         }
      }
      if (st == residue_circle_t::CATION_PI_STACKING) {
         std::string at_name = rc[ires].get_ligand_cation_atom_name();
         lig_build::pos_t atom_pos = mol.get_atom_canvas_position(at_name);
         draw_annotated_stacking_line(atom_pos, rc[ires].pos, st, click_pos);
      }
   }
}


// click_pos is where we recentre in 3D graphics when the annotation
// (line) is clicked.
void
flev_t::draw_annotated_stacking_line(const lig_build::pos_t &ligand_ring_centre,
                                         const lig_build::pos_t &residue_pos,
                                         int stacking_type,
                                         const clipper::Coord_orth &click_pos) {

   double ligand_target_shortening_factor = 3;
   if (stacking_type == residue_circle_t::CATION_PI_STACKING)
      ligand_target_shortening_factor = 6;
   lig_build::pos_t a_to_b_uv = (ligand_ring_centre - residue_pos).unit_vector();
   lig_build::pos_t A = residue_pos + a_to_b_uv * 20;
   lig_build::pos_t B = ligand_ring_centre;
   // just short of the middle of the ring
   lig_build::pos_t C = B - a_to_b_uv * ligand_target_shortening_factor;
   lig_build::pos_t mid_pt = (A + B) * 0.5;
   lig_build::pos_t close_mid_pt_1 = mid_pt - a_to_b_uv * 17;
   lig_build::pos_t close_mid_pt_2 = mid_pt + a_to_b_uv * 17;

   bool start_arrow = 0;
   bool end_arrow = 1;
   std::string stroke_colour = "#008000";

   std::vector<lig_build::pos_t> hex_and_ring_centre(2);
   hex_and_ring_centre[0] = mid_pt - a_to_b_uv * 7;
   hex_and_ring_centre[1] = mid_pt + a_to_b_uv * 7;

   // for cations, we draw a 6-member ring and a "+" text.
   //
   // for pi-pi, we draw 2 6-member rings
   //
   for(int ir=0; ir<2; ir++) {

      bool do_ring = 0;
      bool do_anion = 0;

      if (stacking_type == residue_circle_t::PI_PI_STACKING) {
         do_ring = 1;
      }
      if (stacking_type == residue_circle_t::PI_CATION_STACKING) {
         if (ir == 0) {
            do_ring = 0;
            do_anion = 1;
         } else {
            do_ring = 1;
            do_anion = 0;
         }
      }
      if (stacking_type == residue_circle_t::CATION_PI_STACKING) {
         if (ir == 0) {
            do_ring = 1;
            do_anion = 0;
         } else {
            do_ring = 0;
            do_anion = 1;
         }
      }

      if (do_ring) {
         double angle_step = 60;
         double r = 8; // radius
         for (int ipt=1; ipt<=6; ipt++) {
            int ipt_0 = ipt - 1;
            double theta_deg_1 = 30 + angle_step * ipt;
            double theta_1 = theta_deg_1 * DEG_TO_RAD;
            lig_build::pos_t pt_1 =
               hex_and_ring_centre[ir] + lig_build::pos_t(sin(theta_1), cos(theta_1)) * r;

            double theta_deg_0 = 30 + angle_step * ipt_0;
            double theta_0 = theta_deg_0 * DEG_TO_RAD;
            lig_build::pos_t pt_0 =
               hex_and_ring_centre[ir] + lig_build::pos_t(sin(theta_0), cos(theta_0)) * r;
            // GooCanvasItem *line = goo_canvas_polyline_new_line(group,
            //                                                    pt_1.x, pt_1.y,
            //                                                    pt_0.x, pt_0.y,
            //                                                    "line_width", 1.8,
            //                                                    NULL);
            clipper::Coord_orth *pos_l = new clipper::Coord_orth(click_pos);
            // g_object_set_data_full(G_OBJECT(line), "position", pos_l, g_free);
         }
         // Now the circle in the annotation aromatic ring:
         // GooCanvasItem *ring = goo_canvas_ellipse_new(group,
         //                                              hex_and_ring_centre[ir].x,
         //                                              hex_and_ring_centre[ir].y,
         //                                              4.0, 4.0,
         //                                              "line_width", 1.0,
         //                                              NULL);
         clipper::Coord_orth *pos_r = new clipper::Coord_orth(click_pos);
         // g_object_set_data_full(G_OBJECT(ring), "position", pos_r, g_free);
      }

      if (do_anion) {
         // the "+" symbol for the anion
         //
         // GooCanvasItem *text_1 = goo_canvas_text_new(group,
         //                                             "+",
         //                                             hex_and_ring_centre[ir].x,
         //                                             hex_and_ring_centre[ir].y,
         //                                             -1,
         //                                             GOO_CANVAS_ANCHOR_CENTER,
         //                                             "font", "Sans 12",
         //                                             "fill_color", stroke_colour.c_str(),
         //                                             NULL);
         clipper::Coord_orth *pos_p = new clipper::Coord_orth(click_pos);
         // g_object_set_data_full(G_OBJECT(text_1), "position", pos_p, g_free);
      }
   }


   // GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 2.5, 2.5);
   // GooCanvasItem *item_1 =
   //    goo_canvas_polyline_new_line(group,
   //                                 A.x, A.y,
   //                                 close_mid_pt_1.x, close_mid_pt_1.y,
   //                                 "line-width", 2.5, // in draw_annotated_stacking_line()
   //                                 "line-dash", dash,
   //                                 "stroke-color", stroke_colour.c_str(),
   //                                 NULL);

   // GooCanvasItem *item_2 =
   //    goo_canvas_polyline_new_line(group,
   //                                 close_mid_pt_2.x, close_mid_pt_2.y,
   //                                 C.x, C.y,
   //                                 "line-width", 2.5, // in draw_annotated_stacking_line()
   //                                 "line-dash", dash,
   //                                 // "end_arrow",   end_arrow,
   //                                 "stroke-color", stroke_colour.c_str(),
   //                                 NULL);

   clipper::Coord_orth *pos_1 = new clipper::Coord_orth(click_pos);
   clipper::Coord_orth *pos_2 = new clipper::Coord_orth(click_pos);
   // g_object_set_data_full(G_OBJECT(item_1), "position", pos_1, g_free);
   // g_object_set_data_full(G_OBJECT(item_2), "position", pos_2, g_free);

   // Now the circle blob at the centre of the aromatic ligand ring:
   if (stacking_type != residue_circle_t::CATION_PI_STACKING) {
      // GooCanvasItem *item_o = goo_canvas_ellipse_new(group,
      //                                                B.x, B.y,
      //                                                3.0, 3.0,
      //                                                "line_width", 1.0,
      //                                                "fill_color", stroke_colour.c_str(),
      //                                                NULL);
      clipper::Coord_orth *pos_o = new clipper::Coord_orth(click_pos);
      // g_object_set_data_full(G_OBJECT(item_o), "position", pos_o, g_free);
   }

   clipper::Coord_orth *pos_p = new clipper::Coord_orth(click_pos);
   // g_object_set_data_full(G_OBJECT(group), "position", pos_p, g_free);


}

void
flev_t::draw_bonds_to_ligand() {

   // GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   // GooCanvasLineDash *dash_dash     = goo_canvas_line_dash_new (2, 2.1, 2.1);
   // GooCanvasLineDash *dash_unbroken = goo_canvas_line_dash_new (1, 3);

   // GooCanvasLineDash *dash = dash_dash; // re-assigned later (hopefully)

   for (unsigned int ic=0; ic<residue_circles.size(); ic++) {
      if (residue_circles[ic].bonds_to_ligand.size()) {
         for (unsigned int ib=0; ib<residue_circles[ic].bonds_to_ligand.size(); ib++) {

            try {

               lig_build::pos_t pos = residue_circles[ic].pos;
               std::string at_name = residue_circles[ic].bonds_to_ligand[ib].ligand_atom_name;
               lig_build::pos_t lig_atom_pos = mol.get_atom_canvas_position(at_name);
               lig_build::pos_t rc_to_lig_atom = lig_atom_pos - pos;
               lig_build::pos_t rc_to_lig_atom_uv = rc_to_lig_atom.unit_vector();
               lig_build::pos_t B = lig_atom_pos;

               if (residue_circles[ic].bonds_to_ligand[ib].bond_type != pli::fle_ligand_bond_t::BOND_COVALENT)
                  B -= rc_to_lig_atom_uv * 8; // put the end point (say) of a hydrogen bond
                                              // some pixels away from the atom centre - for better
                                                    // aesthetics.
               lig_build::pos_t A = pos + rc_to_lig_atom_uv * 20;

               // some colours
               std::string blue = "blue";
               std::string green = "darkgreen";
               std::string lime = "#888820";
               std::string stroke_colour = "#111111"; // unset

               // arrows (acceptor/donor) and stroke colour (depending
               // on mainchain or sidechain interaction)
               //
               bool start_arrow = 0;
               bool   end_arrow = 0;
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
                  start_arrow = 0;
                  end_arrow = 0;
                  // dash = dash_dash;
               }
               // GooCanvasItem *item = goo_canvas_polyline_new_line(root,
               //                                                    A.x, A.y,
               //                                                    B.x, B.y,
               //                                                    "line-width", 2.5, // in draw_bonds_to_ligand()
               //                                                    "line-dash", dash,
                //                                                    "start_arrow", start_arrow,
                //                                                    "end_arrow",   end_arrow,
               //                                                    "stroke-color", stroke_colour.c_str(),
               //                                                    NULL);
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
}


void
flev_t::draw_solvent_accessibility_of_atoms() {

   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      lig_build::pos_t pos = mol.atoms[iat].atom_position;
      double sa = mol.atoms[iat].get_solvent_accessibility();
      if (sa  > 0)
         draw_solvent_accessibility_of_atom(pos, sa);
   }
}

void
flev_t::draw_solvent_accessibility_of_atom(const lig_build::pos_t &pos, double sa) {

   double LIGAND_TO_CANVAS_SCALE_FACTOR = 1.0; // this needs to go somewhere

   if (true) {
      int n_circles = int(sa*40) + 1;    // needs fiddling?
      if (n_circles> 10) n_circles = 10; // needs fiddling?

      for (int i=0; i<n_circles; i++) {
         double rad =  LIGAND_TO_CANVAS_SCALE_FACTOR/23.0 * 3.0 * double(i+1); // needs fiddling?

         // CONVERT-TO-SVG
         // GooCanvasItem *circle = goo_canvas_ellipse_new(group,
         //                                                pos.x, pos.y,
         //                                                rad, rad,
         //                                                "line_width", 0.0,
         //                                                "fill-color-rgba", 0x5555cc30,
         //                                                NULL);
         // goo_canvas_item_lower(group, NULL); // to the bottom

      }
   }
}

void
flev_t::ligand_grid::add_for_accessibility(double bash_dist, const lig_build::pos_t &atom_pos) {

   bool debug = 0;
   int grid_extent = 45;

   double LIGAND_TO_CANVAS_SCALE_FACTOR = 1.0; // this needs to go somewhere
   double inv_scale_factor = 1.0/double(LIGAND_TO_CANVAS_SCALE_FACTOR);

   for (int ipos_x= -grid_extent; ipos_x<=grid_extent; ipos_x++) {
      for (int ipos_y= -grid_extent; ipos_y<=grid_extent; ipos_y++) {
         std::pair<int, int> p = canvas_pos_to_grid_pos(atom_pos);
         int ix_grid = ipos_x + p.first;
         int iy_grid = ipos_y + p.second;
         if ((ix_grid >= 0) && (ix_grid < x_size())) {
            if ((iy_grid >= 0) && (iy_grid < y_size())) {

               double d2 = (to_canvas_pos(ix_grid, iy_grid) - atom_pos).lengthsq();
               d2 *= (inv_scale_factor * inv_scale_factor);
               double val = substitution_value(d2, bash_dist);
               if (debug)
                  if (val > 0)
                     std::cout << "adding " << val << " to grid " << ix_grid << " " << iy_grid
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

   bool debug = 1;
   int grid_extent = 40;
   double LIGAND_TO_CANVAS_SCALE_FACTOR = 1.0; // this needs to go somewhere

   double inv_scale_factor = 1.0/double(LIGAND_TO_CANVAS_SCALE_FACTOR);

   for (int ipos_x= -grid_extent; ipos_x<=grid_extent; ipos_x++) {
      for (int ipos_y= -grid_extent; ipos_y<=grid_extent; ipos_y++) {
         std::pair<int, int> p = canvas_pos_to_grid_pos(atom_pos);
         int ix_grid = ipos_x + p.first;
         int iy_grid = ipos_y + p.second;
         if ((ix_grid >= 0) && (ix_grid < x_size())) {
            if ((iy_grid >= 0) && (iy_grid < y_size())) {

               double d2 = (to_canvas_pos(ix_grid, iy_grid) - atom_pos).lengthsq();
               d2 *= (inv_scale_factor * inv_scale_factor);
               // a triangle function, 1 at the atom centre, 0 at 1.5A and beyond
               //
               double dist_at_zero = 2.6;

               double val = 0.0;
               if (d2< dist_at_zero * dist_at_zero)
                  val = - 1.0/dist_at_zero * sqrt(d2) + 1.0;
               grid_[ix_grid][iy_grid] += val;
            }
         }
      }
   }
}

flev_t::ligand_grid::ligand_grid(const lig_build::pos_t &low_x_and_y,
                                 const lig_build::pos_t &high_x_and_y) {

   double extra_extents = 100;
   top_left     = low_x_and_y  - lig_build::pos_t(extra_extents, extra_extents);
   bottom_right = high_x_and_y + lig_build::pos_t(extra_extents, extra_extents);
   scale_fac = 5; // seems good
   double delta_x = bottom_right.x - top_left.x;
   double delta_y = bottom_right.y - top_left.y;
   if (0)
      std::cout << " in making grid, got delta_x and delta_y "
                << delta_x << " " << delta_y << std::endl;
   x_size_ = int(delta_x/scale_fac+1);
   y_size_ = int(delta_y/scale_fac+1);

   std::vector<double> tmp_y (y_size_, 0.0);
   grid_.resize(x_size_);
   for (int i=0; i<x_size_; i++)
      grid_[i] = tmp_y;
   if (true)
      std::cout << "--- ligand_grid constructor, grid has extents "
                << x_size_ << " " << y_size_ << " real " << grid_.size()
                << " " << grid_[0].size()
                << std::endl;
}

std::pair<int, int>
flev_t::ligand_grid::canvas_pos_to_grid_pos(const lig_build::pos_t &pos) const {

   lig_build::pos_t p = pos - top_left;
   int ix = (int)(p.x/scale_fac);
   int iy = (int)(p.y/scale_fac);
   return std::pair<int, int> (ix, iy);
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
