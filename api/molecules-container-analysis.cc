/*
 * api/molecules-container-analysis.cc
 *
 * Copyright 2020, 2021, 2-22 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */


// Validation, analysis and measurement functions of molecules_container_t

#include <filesystem>

#include <sys/types.h> // for stating
#include <sys/stat.h>

#include "molecules-container.hh"
#include "geometry/protein-geometry.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "ideal/pepflip.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/secondary-structure-headers.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/read-sm-cif.hh"
#include "coot-utils/read-amber-trajectory.hh"
#include "coot-utils/json.hpp"

#include "coords/Bond_lines.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-extras.hh"

#include "ligand/ligand.hh"

#include "mmdb2/mmdb_atom.h"
#include "utils/logging.hh"
extern logging logger;




//! rotamer validation information
coot::validation_information_t
molecules_container_t::rotamer_analysis(int imol_model) const {

   coot::validation_information_t r;
   r.name = "Rotamer analysis";
#ifdef EMSCRIPTEN
   r.type = "PROBABILITY";
#else
   r.type = coot::PROBABILITY;
#endif

   if (is_valid_model_molecule(imol_model)) {

      mmdb::Manager *mol = molecules[imol_model].atom_sel.mol;

      // fill these
      mmdb::PResidue *SelResidues = 0;
      int nSelResidues = 0;

      int selHnd = mol->NewSelection(); // yes, it's deleted.
      int imod = 1; // multiple models don't work on validation graphs

      mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                           "*", // chain_id
                           mmdb::ANY_RES, "*",
                           mmdb::ANY_RES, "*",
                           "*",  // residue name
                           "*",  // Residue must contain this atom name?
                           "*",  // Residue must contain this Element?
                           "*",  // altLocs
                           mmdb::SKEY_NEW // selection key
                           );
      mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
         coot::residue_spec_t res_spec(residue_p);
         mmdb::PAtom *residue_atoms=0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

         // double residue_density_score = coot::util::map_score(residue_atoms, n_residue_atoms, xmap, 1);

         // Use the heavy-atom count, not the total atom count, so that the set of
         // residues plotted does not depend on whether hydrogens have been added.
         // (GLY and ALA have <= 5 heavy atoms and so are correctly skipped here.)
         int n_heavy_atoms = 0;
         for (int iat=0; iat<n_residue_atoms; iat++) {
            std::string ele(residue_atoms[iat]->element);
            if (ele != " H" && ele != " D")
               n_heavy_atoms++;
         }

         if (n_heavy_atoms > 5) {

            std::string res_name = residue_p->GetResName();
            if (true) {

               coot::rotamer rot(residue_p);
               coot::rotamer_probability_info_t rpi = rot.probability_of_this_rotamer();
               double prob = rpi.probability;

               std::string l = res_spec.label();
               std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
               const std::string &chain_id = res_spec.chain_id;
               int this_resno = res_spec.res_no;
               coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
               coot::residue_validation_information_t rvi(res_spec, atom_spec, prob, l);
               r.add_residue_validation_information(rvi, chain_id);
            }
         }
      }
      mol->DeleteSelection(selHnd);
   }
   r.set_min_max();
   return r;
}

double
molecules_container_t::phi_psi_probability(const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) const {

      const clipper::Ramachandran *rama = &rc.rama;

      if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
      if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

      // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
      // if (phi_psi.is_pre_pro())
      // if (phi_psi.residue_name() != "GLY")
      // rama = &rc.rama_pre_pro;

      double rama_prob = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                           clipper::Util::d2rad(phi_psi.psi()));
      return rama_prob;
}

//! ramachandran validation information (formatted for a graph, not 3d)
coot::validation_information_t
molecules_container_t::ramachandran_analysis(int imol_model) const {

   coot::validation_information_t vi;
   vi.name = "Ramachandran plot Probability";
#ifdef EMSCRIPTEN
   vi.type = "PROBABILITY";
#else
   vi.type = coot::PROBABILITY;
#endif
   std::vector<coot::phi_psi_prob_t> rv = ramachandran_validation(imol_model);

   for (unsigned int i=0; i<rv.size(); i++) {
      std::string chain_id = rv[i].phi_psi.chain_id;
      coot::residue_spec_t residue_spec(rv[i].phi_psi.chain_id, rv[i].phi_psi.residue_number, rv[i].phi_psi.ins_code);
      double pr = rv[i].probability;
      std::string label = rv[i].phi_psi.chain_id + std::string(" ") + std::to_string(rv[i].phi_psi.residue_number);
      if (! rv[i].phi_psi.ins_code.empty())
         label += std::string(" ") + rv[i].phi_psi.ins_code;
      coot::atom_spec_t atom_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " CA ", "");
      coot::residue_validation_information_t rvi(residue_spec, atom_spec, pr, label);
      if (false)
         std::cout << "         " << residue_spec << " " << rv[i].phi_psi.phi() << " " << rv[i].phi_psi.psi()
                   << " pr " << pr << " " << std::endl;
      vi.add_residue_validation_information(rvi, chain_id);
   }
   vi.set_min_max();
   return vi;
}

//! ramachandran validation information (formatted for a graph, not 3d) for a given chain in a given molecule
//! 20230127-PE This function does not exist yet.
//!
//! @returns a `coot::validation_information_t`
coot::validation_information_t
molecules_container_t::ramachandran_analysis_for_chain(int imol_model, const std::string &user_chain_id) const {

   coot::validation_information_t vi;
   vi.name = "Ramachandran plot Probability";
#ifdef EMSCRIPTEN
   vi.type = "PROBABILITY";
#else
   vi.type = coot::PROBABILITY;
#endif
   std::vector<coot::phi_psi_prob_t> rv = ramachandran_validation(imol_model);

   for (unsigned int i=0; i<rv.size(); i++) {
      std::string chain_id = rv[i].phi_psi.chain_id;
      if (chain_id != user_chain_id) continue;
      coot::residue_spec_t residue_spec(rv[i].phi_psi.chain_id, rv[i].phi_psi.residue_number, rv[i].phi_psi.ins_code);
      double pr = rv[i].probability;
      std::string label = rv[i].phi_psi.chain_id + std::string(" ") + std::to_string(rv[i].phi_psi.residue_number);
      if (! rv[i].phi_psi.ins_code.empty())
         label += std::string(" ") + rv[i].phi_psi.ins_code;
      coot::atom_spec_t atom_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " CA ", "");
      coot::residue_validation_information_t rvi(residue_spec, atom_spec, pr, label);
      if (false)
         std::cout << "         " << residue_spec << " " << rv[i].phi_psi.phi() << " " << rv[i].phi_psi.psi()
                   << " pr " << pr << " " << std::endl;
      vi.add_residue_validation_information(rvi, chain_id);
   }
   vi.set_min_max();
   return vi;
}


//! peptide omega validation information
//! @returns a `validation_information_t`
coot::validation_information_t
molecules_container_t::peptide_omega_analysis(int imol) const {

   coot::validation_information_t vi;
   vi.name = "Peptide Omega Deviation";
#ifdef EMSCRIPTEN
   vi.type = "TORSION_ANGLE";
#else
   vi.type = coot::TORSION_ANGLE;
#endif

   if (is_valid_model_molecule(imol)) {

      bool mark_cis_peptides_as_bad_flag = false;
      bool m = mark_cis_peptides_as_bad_flag;
      std::vector<std::string> chain_ids = molecules[imol].chains_in_model();
      for (const auto &chain_id : chain_ids) {
         coot::chain_validation_information_t cvi(chain_id);
         coot::omega_distortion_info_container_t odi = molecules.at(imol).peptide_omega_analysis(geom, chain_id, m);
         for (const auto &od : odi.omega_distortions) {
            // oops - we have forgotten about the insertion code.
            coot::residue_spec_t res_spec(chain_id, od.resno, "");
            coot::atom_spec_t atom_spec(chain_id, od.resno, "", " CA ", "");
            std::string label = od.info_string;
            coot::residue_validation_information_t rvi(res_spec, atom_spec, od.distortion, label);
            cvi.add_residue_validation_information(rvi);
         }
         vi.cviv.push_back(cvi);
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return vi;
}


std::vector<coot::phi_psi_prob_t>
molecules_container_t::ramachandran_validation(int imol) const {

   std::vector<coot::phi_psi_prob_t> v;
   if (is_valid_model_molecule(imol))
      v = molecules[imol].ramachandran_validation(ramachandrans_container);
   return v;
}

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

coot::simple_mesh_t
molecules_container_t::get_ramachandran_validation_markup_mesh(int imol) const {

   // this function should be pushed into the coot::molecule_t class
   // (which means that the mesh will be copied)

   unsigned int num_subdivisions = 2;  // pass this
   float rama_ball_radius = 0.5;

   auto prob_raw_to_colour_rotation = [] (float prob) {
                                         if (prob > 0.5) prob = 0.5; // 0.4 and 2.5 f(for q) might be better (not tested)
                                         // good probabilities have q = 0
                                         // bad probabilities have q 0.66
                                         double q = (1.0 - 2.0 * prob);
                                         q = pow(q, 20.0);
                                         return q;
   };

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
   };

   auto phi_psi_probability = [] (const coot::util::phi_psi_t &phi_psi, const ramachandrans_container_t &rc) {

      const clipper::Ramachandran *rama = &rc.rama;

      if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
      if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

      // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
      // if (phi_psi.is_pre_pro())
      // if (phi_psi.residue_name() != "GLY")
      // rama = &rc.rama_pre_pro;

      double rama_prob = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                           clipper::Util::d2rad(phi_psi.psi()));
      return rama_prob;
   };

   auto test_ramachandran_probabilities = [] (const ramachandrans_container_t &rc) {

      std::vector<const clipper::Ramachandran *> ramas = { &rc.rama, &rc.rama_gly, &rc.rama_pro, &rc.rama_non_gly_pro };

      for (unsigned int ir=0; ir<ramas.size(); ir++) {
         for (unsigned int i=0; i<10; i++) {
            for (unsigned int j=0; j<10; j++) {
               double phi = static_cast<double>(i * 36.0) - 180.0;
               double psi = static_cast<double>(j * 36.0) - 180.0;
               double p = rc.rama.probability(phi, psi);
               std::cout << ir << "   "
                         << std::setw(10) << phi << " " << std::setw(10) << psi << " "
                         << std::setw(10) << p << std::endl;
            }
         }
      }
   };

   // test_ramachandran_probabilities(ramachandrans_container); // don't use rama_pre_pro without CLIPPER_HAS_TOP8000

   coot::simple_mesh_t mesh;

   // 20221126-PE Calm down the ultra-bright rama balls:
   float sober_factor = 0.75;

   if (is_valid_model_molecule(imol)) {

      std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaball = tessellate_octasphere(num_subdivisions);

      std::vector<coot::phi_psi_prob_t> ramachandran_goodness_spots = ramachandran_validation(imol);
      // now convert positions into meshes of balls
      int n_ramachandran_goodness_spots = ramachandran_goodness_spots.size();
      for (int i=0; i<n_ramachandran_goodness_spots; i++) {
         const coot::Cartesian &position = ramachandran_goodness_spots[i].position;
         // std::cout << "goodness spot " << i << " position " << position << std::endl;
         const coot::phi_psi_prob_t &phi_psi = ramachandran_goodness_spots[i];
         double prob_raw = phi_psi.probability;
         double q = prob_raw_to_colour_rotation(prob_raw);
         coot::colour_holder col = coot::colour_holder(q, 0.0, 1.0, false, std::string(""));
         glm::vec3 ball_position = cartesian_to_glm(position);
         unsigned int idx_base = mesh.vertices.size();
         unsigned int idx_tri_base = mesh.triangles.size();
         for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
            glm::vec4 col_v4(sober_factor * col.red, sober_factor * col.green, sober_factor * col.blue, 1.0f);
            const glm::vec3 &vertex_position = octaball.first[ibv];
            coot::api::vnc_vertex vertex(ball_position + rama_ball_radius * vertex_position, vertex_position, col_v4);
            mesh.vertices.push_back(vertex);
         }
         std::vector<g_triangle> octaball_triangles = octaball.second;
         mesh.triangles.insert(mesh.triangles.end(), octaball_triangles.begin(), octaball_triangles.end());

         for (unsigned int k=idx_tri_base; k<mesh.triangles.size(); k++)
            mesh.triangles[k].rebase(idx_base);
      }
   }
   return mesh;
}


std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::difference_map_peaks(int imol_map, int imol_protein, float n_rmsd) const {

   std::vector<coot::molecule_t::interesting_place_t> v;
   if (is_valid_model_molecule(imol_protein)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *m = get_mol(imol_protein);
         v = molecules[imol_map].difference_map_peaks(m, n_rmsd);
      } else {
         std::cout << "debug:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_map << std::endl;
      }
   } else {
      std::cout << "debug:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_protein << std::endl;
   }
   return v;
}



std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::unmodelled_blobs(int imol_model, int imol_map, float rmsd_cut_off) const {

   std::vector<coot::molecule_t::interesting_place_t> v;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

         coot::ligand lig;

         short int mask_waters_flag = true;
         float sigma = molecules[imol_map].get_map_rmsd_approx();
         lig.import_map_from(molecules[imol_map].xmap, sigma);
         lig.set_map_atom_mask_radius(1.9); // Angstrom
         lig.mask_map(molecules[imol_model].atom_sel.mol, mask_waters_flag);
         float sigma_cut_off = rmsd_cut_off;
         std::cout << "Unmodelled blobs using sigma cut off " << sigma_cut_off << std::endl;
         int n_cycles = 1;
         lig.water_fit(sigma_cut_off, n_cycles);
         std::vector<std::pair<clipper::Coord_orth, double> > big_blobs = lig.big_blobs();
         int n_big_blobs = lig.big_blobs().size();
         if (n_big_blobs > 0) {
            for (unsigned int i=0; i<big_blobs.size(); i++) {
               std::string l = std::string("Blob ") + std::to_string(i+1);
               clipper::Coord_orth pt = big_blobs[i].first;
               coot::molecule_t::interesting_place_t ip("Unmodelled Blob", pt, l);
               ip.set_feature_value(big_blobs[i].second);
               v.push_back(ip);
            }
         }
      }
   }
   return v;
}





// put this in a new file molecules_container_validation.cc

#include "coot-utils/pepflip-using-difference-map.hh"


std::vector<coot::molecule_t::interesting_place_t>
molecules_container_t::pepflips_using_difference_map(int imol_coords, int imol_difference_map, float n_sigma) const {

   auto mmdb_to_clipper = [] (mmdb::Atom *at) {
      return clipper::Coord_orth(at->x, at->y, at->z);
   };

   std::vector<coot::molecule_t::interesting_place_t> v;

   if (is_valid_model_molecule(imol_coords)) {
      if (is_valid_map_molecule(imol_difference_map)) {
         if (molecules[imol_difference_map].is_difference_map_p()) {
            const clipper::Xmap<float> &diff_xmap = molecules[imol_difference_map].xmap;
            mmdb::Manager *mol = get_mol(imol_coords);
            coot::pepflip_using_difference_map pf(mol, diff_xmap);
            std::vector<coot::residue_spec_t> flips = pf.get_suggested_flips(n_sigma);
            for (std::size_t i=0; i<flips.size(); i++) {
               const auto &res_spec = flips[i];
               mmdb::Residue *residue_this_p = get_residue(imol_coords, res_spec);
               if (residue_this_p) {
                  coot::residue_spec_t res_spec_next =  res_spec.next();
                  mmdb::Residue *residue_next_p = get_residue(imol_coords, res_spec);
                  if (residue_next_p) {
                     std::string feature_type = "Difference Map Suggest Pepflip";
                     std::string label = "Flip: " + res_spec.format();
                     mmdb::Atom *at_1 = residue_this_p->GetAtom(" CA ");
                     mmdb::Atom *at_2 = residue_next_p->GetAtom(" CA ");
                     if (at_1 && at_2) {
                        clipper::Coord_orth pt_1 = mmdb_to_clipper(at_1);
                        clipper::Coord_orth pt_2 = mmdb_to_clipper(at_2);
                        clipper::Coord_orth pos = 0.5 * (pt_1 + pt_2);
                        float f = static_cast<float>(i)/static_cast<float>(flips.size());
                        float badness = 20.0 + 50.0 * (1.0 - f);
                        coot::molecule_t::interesting_place_t ip(feature_type, res_spec, pos, label);
                        ip.set_badness_value(badness);
                        v.push_back(ip);
                     }
                  }
               }
            }
         }
      }
   }
   std::cout << "DEBUG:: pepflips_using_difference_map() returns " << v.size() << " flips" << std::endl;
   return v;

}

//! return@ an object that has information about residues without dictionaries and residues with missing atom
//! in the the specified molecule
coot::util::missing_atom_info
molecules_container_t::missing_atoms_info_raw(int imol) {

   coot::util::missing_atom_info mai;

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      bool do_missing_hydrogen_atoms_flag = false;
      mai = coot::util::missing_atoms(mol, do_missing_hydrogen_atoms_flag, &geom);
   }
   return mai;
}


//! @return an object that has information about residues without dictionaries and residues with missing atom
//! in the the specified molecule
std::vector<coot::residue_spec_t>
molecules_container_t::residues_with_missing_atoms(int imol) {

   std::vector<coot::residue_spec_t> v;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      bool do_missing_hydrogen_atoms_flag = false;
      coot::util::missing_atom_info mai = coot::util::missing_atoms(mol, do_missing_hydrogen_atoms_flag, &geom);
      for (unsigned int i=0; i<mai.residues_with_missing_atoms.size(); i++) {
         mmdb::Residue *r = mai.residues_with_missing_atoms[i];
         v.push_back(coot::residue_spec_t(r));
      }
   }
   return v;
}

//! @return the instanced mesh for the specified ligand
coot::instanced_mesh_t
molecules_container_t::contact_dots_for_ligand(int imol, const std::string &cid,
                                               unsigned int num_subdivisions) const {

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      im = molecules[imol].contact_dots_for_ligand(cid, geom, num_subdivisions);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return im;
}


//! @return the instanced mesh for the specified molecule
coot::instanced_mesh_t
molecules_container_t::all_molecule_contact_dots(int imol, unsigned int num_subdivisions) {

   coot::instanced_mesh_t im;
   if (is_valid_model_molecule(imol)) {
      im = molecules[imol].all_molecule_contact_dots(geom, num_subdivisions);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return im;
}

//! get the mesh for ligand validation vs dictionary, coloured by badness.
//! greater then 3 standard deviations is fully red.
//! Less than 0.5 standard deviations is fully green.
coot::simple_mesh_t
molecules_container_t::get_mesh_for_ligand_validation_vs_dictionary(int imol, const std::string &ligand_cid) {

   coot::simple_mesh_t m;
   if (is_valid_model_molecule(imol)) {
      m = molecules[imol].get_mesh_for_ligand_validation_vs_dictionary(ligand_cid, geom, thread_pool);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return m;

}

//! ligand validation - basically we do the same as the above function, but the
//! return type is validation data, not a mesh
//!
//! @return a vector of `geometry_distortion_info_container_t`
std::vector<coot::geometry_distortion_info_pod_container_t>
molecules_container_t::get_ligand_validation_vs_dictionary(int imol,
                                                           const std::string &ligand_cid,
                                                           bool with_nbcs) {

   std::vector<coot::geometry_distortion_info_pod_container_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].geometric_distortions_for_one_residue_from_mol(ligand_cid, with_nbcs, geom, thread_pool);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! General fragment distortion analysis
//!
//! @param imol is the model molecule index
//! @param selection_cid is the selection CID e.g "//A/15-23"
//! @param include_non_bonded_contacts is the flag to include non bonded contacts
//!
//! @return a vector/list of interesting geometry
std::vector<coot::geometry_distortion_info_pod_container_t>
molecules_container_t::get_validation_vs_dictionary_for_selection(int imol,
                                                                  const std::string &selection_cid,
                                                                  bool include_non_bonded_contacts) {

   std::vector<coot::geometry_distortion_info_pod_container_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].geometric_distortions_for_selection_from_mol(selection_cid,
                                                                       include_non_bonded_contacts,
                                                                       geom,
                                                                       thread_pool);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}


//! Get ligand distortion
//!
//! a more simple interface to the above
//!
//! @return a pair: the first is the status (1 for OK, 0 for fail)
//!
// should this be const?
std::pair<int, double>
molecules_container_t::get_ligand_distortion(int imol, const std::string &ligand_cid, bool with_nbcs) {

   int status = 0;
   double d = 0;
   if (is_valid_model_molecule(imol)) {
      std::pair<int, double> p =
         molecules[imol].simple_geometric_distortions_from_mol(ligand_cid, with_nbcs, geom, thread_pool);
      status = p.first;
      d = p.second;
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return std::make_pair(status, d);
}

#include "coot-utils/find-water-baddies.hh"

//! check waters, implicit OR
//! return a vector of atom specifiers
std::vector <coot::atom_spec_t>
molecules_container_t::find_water_baddies(int imol_model, int imol_map,
                                          float b_factor_lim,
                                          float outlier_sigma_level,
                                          float min_dist, float max_dist,
                                          bool ignore_part_occ_contact_flag,
                                          bool ignore_zero_occ_flag) {

   std::vector <coot::atom_spec_t> v;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

         float map_sigma = molecules[imol_map].get_map_rmsd_approx();
         v = coot::find_water_baddies_OR(molecules[imol_model].atom_sel,
                                         b_factor_lim,
                                         molecules[imol_map].xmap,
                                         map_sigma,
                                         outlier_sigma_level,
                                         min_dist, max_dist,
                                         ignore_part_occ_contact_flag,
                                         ignore_zero_occ_flag);

         std::cout << "........... find_water_baddies_OR() returned " << v.size() << " water baddies " << std::endl;
      } else {
         std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid map molecule " << imol_model << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol_map << std::endl;
   }
   return v;

}

//! return the hb_tye for the given atom. On failure return an empty string
std::string
molecules_container_t::get_hb_type(const std::string &compound_id, int imol_enc, const std::string &atom_name) const {

   coot::hb_t hbt = geom.get_h_bond_type(atom_name, compound_id, imol_enc);
   std::string hb;
   if (hbt == coot::HB_UNASSIGNED) hb = "HB_UNASSIGNED";
   if (hbt == coot::HB_NEITHER)    hb = "HB_NEITHER";
   if (hbt == coot::HB_DONOR)      hb = "HB_DONOR";
   if (hbt == coot::HB_ACCEPTOR)   hb = "HB_ACCEPTOR";
   if (hbt == coot::HB_BOTH)       hb = "HB_BOTH";
   if (hbt == coot::HB_HYDROGEN)   hb = "HB_HYDROGEN";
   return hb;
}


//! @return a vector of string pairs that were part of a gphl_chem_comp_info.
//!  return an empty vector on failure to find any such info.
std::vector<std::pair<std::string, std::string> >
molecules_container_t::get_gphl_chem_comp_info(const std::string &compound_id, int imol_enc) {

   std::vector<std::pair<std::string, std::string> > v;
   std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
      geom.get_monomer_restraints(compound_id, imol_enc);
   if (r_p.first) {
      v = r_p.second.gphl_chem_comp_info.info;
   }
   return v;
}

//! get a list of atom names and their associated atedrg atom types
//!
//! @return a list of atom names and their associated atedrg atom types, return an empty list
//! on failure (atoms types are not in the dictionary or atom failure to look up the compound id)l
std::vector<std::pair<std::string, std::string> >
molecules_container_t::get_acedrg_atom_types(const std::string &compound_id, int imol_enc) const {

   std::vector<std::pair<std::string, std::string> > v;
   std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
      geom.get_monomer_restraints(compound_id, imol_enc);
   if (r_p.first) {
      const auto &restraints = r_p.second;
      const auto &atom_info = restraints.atom_info;
      for (unsigned int iat=0; iat<atom_info.size(); iat++) {
         const auto &atom = atom_info[iat];
         const auto &atom_id = atom.atom_id;
         const auto &acedrg_atom_type = atom.acedrg_atom_type;
         if (! acedrg_atom_type.empty()) {
            auto pair = std::make_pair(atom_id, acedrg_atom_type);
            v.push_back(pair);
         }
      }
   }
   return v;

}

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/cod-atom-types.hh"
#endif

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

std::vector<std::pair<std::string, cod::atom_type_t> >
molecules_container_t::get_computed_acedrg_atom_types(const std::string &compound_id, int imol_enc) {

   std::vector<std::pair<std::string, cod::atom_type_t> > v;

   std::pair<bool, coot::dictionary_residue_restraints_t> r_p =
      geom.get_monomer_restraints(compound_id, imol_enc);
   if (r_p.first) {
      const auto &restraints = r_p.second;
      bool idealised_flag = true;
      mmdb::Manager *mol = geom.mol_from_dictionary(compound_id, imol_enc, idealised_flag);
      if (mol) {
         mmdb::Residue *residue_p = coot::util::get_first_residue(mol);
         if (residue_p) {
            try {
               RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(residue_p, imol_enc, geom);
               cod::atom_types_t t;
               std::vector<cod::atom_type_t> types = t.get_cod_atom_types(rdkm);
               if (types.size() == restraints.atom_info.size()) {
                  // Pair by atom name, not by index: types[] is in RDKit atom order,
                  // which is not guaranteed to match the dictionary atom_info order.
                  // The RDKit atoms carry a "name" property (set in rdkit_mol's mmdb
                  // conversion) that matches the dictionary atom_id.  types[idx]
                  // corresponds to the idx-th atom of rdkm (same iteration order as
                  // get_cod_atom_types()).
                  std::map<std::string, cod::atom_type_t> name_to_type;
                  unsigned int idx = 0;
                  for (RDKit::ROMol::AtomIterator ai=rdkm.beginAtoms(); ai!=rdkm.endAtoms(); ++ai, ++idx) {
                     if ((*ai)->hasProp("name")) {
                        std::string name;
                        (*ai)->getProp("name", name);
                        name_to_type[name] = types[idx];
                     }
                  }
                  // The RDKit "name" property is the (4-char padded) mmdb atom name, i.e.
                  // dict_atom::atom_id_4c, so match on that first, then the unpadded atom_id.
                  auto find_type = [&name_to_type] (const coot::dict_atom &a) {
                     std::map<std::string, cod::atom_type_t>::const_iterator it = name_to_type.find(a.atom_id_4c);
                     if (it == name_to_type.end()) it = name_to_type.find(a.atom_id);
                     return it;
                  };
                  bool all_found = true;
                  for (unsigned int iat=0; iat<restraints.atom_info.size(); iat++)
                     if (find_type(restraints.atom_info[iat]) == name_to_type.end())
                        { all_found = false; break; }
                  for (unsigned int iat=0; iat<restraints.atom_info.size(); iat++) {
                     const auto &atom_id = restraints.atom_info[iat].atom_id;
                     if (all_found)
                        v.push_back(std::make_pair(atom_id, find_type(restraints.atom_info[iat])->second));
                     else
                        v.push_back(std::make_pair(atom_id, types[iat])); // fall back to index pairing
                  }
               }
            }
            catch (const std::runtime_error &rte) {
               std::cout << "WARNING:: get_computed_acedrg_atom_types() " << compound_id
                         << " " << rte.what() << std::endl;
            }
            catch (const std::exception &e) {
               std::cout << "WARNING:: get_computed_acedrg_atom_types() " << compound_id
                         << " " << e.what() << std::endl;
            }
         }
         delete mol;
      }
   } else {
      std::cout << "DEBUG:: no restraints found for compound_id: \"" << compound_id << "\" " << imol_enc << std::endl;
   }

   return v;
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#include "coot-utils/cremer-pople.hh"

coot::cremer_pople_info_t
molecules_container_t::get_cremer_pople(int imol,
                                        const std::string &residue_cid,
                                        const std::vector<std::string> &ordered_atom_names,
                                        const std::string &up_reference_atom_name,
                                        const std::string &alt_conf) {
   coot::cremer_pople_info_t info;
   if (! is_valid_model_molecule(imol)) return info;

   mmdb::Residue *residue_p = molecules[imol].cid_to_residue(residue_cid);
   if (! residue_p) return info;

   auto strip = [] (const std::string &s) {
      std::string r = s;
      while (!r.empty() && r.front() == ' ') r.erase(0,1);
      while (!r.empty() && r.back()  == ' ') r.pop_back();
      return r;
   };

   // When alt_conf is given explicitly, match it exactly (an alt_conf that
   // doesn't exist on this residue yields "not found", not a silent fallback).
   //
   // When alt_conf is empty, do NOT just take the first atom of that name in
   // GetAtomTable() order: on a residue with genuine alternate conformers the
   // table order need not agree from atom to atom (e.g. C1's "A" record can
   // precede its "B" while C2's "B" precedes its "A"), so "first match" can
   // silently mix atoms from two different physical conformers into one ring
   // and return a plausible-looking but meaningless result. Instead: prefer
   // a blank-altLoc atom if one exists; otherwise use a non-blank altLoc only
   // if it is the sole distinct altLoc present for that atom name -- two or
   // more distinct non-blank altLocs is ambiguous, so refuse (return false)
   // rather than guess.
   auto position_of = [&] (const std::string &name, clipper::Coord_orth &pos) -> bool {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

      if (! alt_conf.empty()) {
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (at->isTer()) continue;
            if (strip(at->GetAtomName()) != name) continue;
            if (std::string(at->altLoc) != alt_conf) continue;
            pos = clipper::Coord_orth(at->x, at->y, at->z);
            return true;
         }
         return false;
      }

      mmdb::Atom *blank_match = nullptr;
      mmdb::Atom *first_non_blank_match = nullptr;
      std::string first_non_blank_altloc;
      bool ambiguous = false;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (at->isTer()) continue;
         if (strip(at->GetAtomName()) != name) continue;
         std::string this_altloc(at->altLoc);
         if (this_altloc.empty()) {
            if (! blank_match) blank_match = at;
         } else {
            if (! first_non_blank_match) {
               first_non_blank_match = at;
               first_non_blank_altloc = this_altloc;
            } else if (this_altloc != first_non_blank_altloc) {
               ambiguous = true;
            }
         }
      }
      mmdb::Atom *chosen = blank_match ? blank_match : (ambiguous ? nullptr : first_non_blank_match);
      if (! chosen) return false;
      pos = clipper::Coord_orth(chosen->x, chosen->y, chosen->z);
      return true;
   };

   std::vector<clipper::Coord_orth> ring;
   for (const auto &name : ordered_atom_names) {
      clipper::Coord_orth p;
      if (! position_of(name, p)) return info;   // missing atom -> not filled
      ring.push_back(p);
   }

   clipper::Coord_orth up;
   const clipper::Coord_orth *up_p = nullptr;
   if (! up_reference_atom_name.empty())
      if (position_of(up_reference_atom_name, up)) up_p = &up;

   coot::cremer_pople_t cp(ring, up_p);
   if (! cp.filled) return info;

   info.filled    = true;
   info.ring_size = static_cast<int>(cp.n_ring);
   info.Q  = cp.Q;
   info.q2 = cp.q2;
   info.q3 = cp.q3;
   info.phi = clipper::Util::rad2d(cp.phi);
   info.theta = (cp.n_ring == 6) ? clipper::Util::rad2d(cp.theta) : 0.0;
   return info;
}

//! get acedrg types for ligand bonds
//! @return a vector of `acedrg_types_for_residue_t`
coot::acedrg_types_for_residue_t
molecules_container_t::get_acedrg_atom_types_for_ligand(int imol, const std::string &residue_cid) const {

   coot::acedrg_types_for_residue_t types;

   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = molecules[imol].get_residue(residue_cid);
      if (residue_p) {
         int imol_enc = imol;
         types = coot::get_acedrg_types_for_residue(residue_p, imol_enc, geom);
      }
   }
   return types;
}




//! Fourier Shell Correlation (FSC) between maps
//! @return a vector or pairs of graph points (resolution, correlation)
std::vector<std::pair<double, double> >
molecules_container_t::fourier_shell_correlation(int imol_map_1, int imol_map_2) const {

   std::vector<std::pair<double, double> > v;

   if (is_valid_map_molecule(imol_map_1)) {
      if (is_valid_map_molecule(imol_map_2)) {
         const clipper::Xmap<float> &xmap_1 = molecules[imol_map_1].xmap;
         const clipper::Xmap<float> &xmap_2 = molecules[imol_map_2].xmap;
         auto fsc = coot::util::fsc(xmap_1, xmap_2);
         if (! fsc.empty()) {
            v.resize(fsc.size());
            for (unsigned int i=0; i<fsc.size(); i++) {
               v[i].first  = fsc[i].first.invresolsq_limit();
               v[i].second = fsc[i].second;
            }
         }
      }
   }
   return v;
}


//! not const because it can dynamically add dictionaries
coot::atom_overlaps_dots_container_t
molecules_container_t::get_overlap_dots(int imol) {

   coot::atom_overlaps_dots_container_t aodc;
   if (is_valid_model_molecule(imol)) {
      aodc = molecules[imol].get_overlap_dots(&geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return aodc;
}

//! not const because it can dynamically add dictionaries
coot::atom_overlaps_dots_container_t
molecules_container_t::get_overlap_dots_for_ligand(int imol, const std::string &cid_ligand) {

   coot::atom_overlaps_dots_container_t aodc;
   if (is_valid_model_molecule(imol)) {
      aodc = molecules[imol].get_overlap_dots_for_ligand(cid_ligand, &geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return aodc;
}



//! not const because it can dynamically add dictionaries
std::vector<coot::plain_atom_overlap_t>
molecules_container_t::get_atom_overlaps(int imol) {

   std::vector<coot::plain_atom_overlap_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_atom_overlaps(&geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! Get the atom overlap score
//!
//! @param imol the model molecule index
//! @return the overlap score - a negative number indicates failure
float
molecules_container_t::get_atom_overlap_score(int imol) {

   float v = -1.0;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_atom_overlap_score(&geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }

   return v;

}



//! not const because it can dynamically add dictionaries
std::vector<coot::plain_atom_overlap_t>
molecules_container_t::get_overlaps_for_ligand(int imol, const std::string &cid_ligand) {

   std::vector<coot::plain_atom_overlap_t> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_overlaps_for_ligand(cid_ligand, &geom);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! get the median temperature factor for the model
//! @return a negative number on failure.
float
molecules_container_t::get_median_temperature_factor(int imol) const {

   float b_factor = -1.1;
   if (is_valid_model_molecule(imol)) {
      b_factor = molecules[imol].get_median_temperature_factor();
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return b_factor;
}

//! Get the atom temperature factor
//!
//! @param imol is the model molecule index
//! @param atom_cid is the selection cid for the atom
//!
//! @return a negative number on failure, otherwise the temperature factor
float
molecules_container_t::get_temperature_factor_of_atom(int imol, const std::string &atom_cid) const {

   float b_factor = -1.1;
   if (is_valid_model_molecule(imol)) {
      b_factor = molecules[imol].get_temperature_factor_of_atom(atom_cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return b_factor;
}



// return the atom name match on superposing the atoms of the given dictionaries
std::map<std::string, std::string>
molecules_container_t::dictionary_atom_name_map(const std::string &comp_id_1, int imol_1, const std::string &comp_id_2, int imol_2) {

   std::map<std::string, std::string> m;

   std::pair<bool, coot::dictionary_residue_restraints_t> r_p_1 = geom.get_monomer_restraints(comp_id_1, imol_1);
   std::pair<bool, coot::dictionary_residue_restraints_t> r_p_2 = geom.get_monomer_restraints(comp_id_2, imol_2);
   if (r_p_1.first) {
      if (r_p_2.first) {
         const coot::dictionary_residue_restraints_t &dict_1 = r_p_1.second;
         const coot::dictionary_residue_restraints_t &dict_2 = r_p_2.second;
         coot::dictionary_match_info_t dm = dict_1.match_to_reference(dict_2, nullptr, comp_id_1, "dummy");
         if (false) {
            std::cout << "There are " << dm.same_names.size() << " atoms with the same name" << std::endl;
            std::cout << "There are " << dm.name_swaps.size() << " atoms with the different names" << std::endl;
         }
         for (const auto &name : dm.same_names)
            m[name] = name;
         for (const auto &name : dm.name_swaps) {
            m[name.first] = name.second;
         }
      }
   }
   return m;
}



//! Get the residue CA position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
molecules_container_t::get_residue_CA_position(int imol, const std::string &cid) const {

   std::vector<double> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_residue_CA_position(cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;

}

//! Get the average residue position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
molecules_container_t::get_residue_average_position(int imol, const std::string &cid) const {

   std::vector<double> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_residue_average_position(cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! Get the average residue side-chain position
//!
//! @return a vector. The length of the vector is 0 on failure, otherwise it is the x,y,z values
std::vector<double>
molecules_container_t::get_residue_sidechain_average_position(int imol, const std::string &cid) const {

   std::vector<double> v;
   if (is_valid_model_molecule(imol)) {
      v = molecules[imol].get_residue_sidechain_average_position(cid);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return v;
}

//! Get the torsion of the specified atom in the specified residue
//!
//! @param imol is the model molecule index
//! @param cid is the selection CID, e.g. //A/15 (residue 15 in chain A)
//! @param atom_names is a list of atom names, e.g. ["CA", "CB", "CG", "CD"]
//!
//! @return a pair, the first of which is a succes status (1 success, 0 failure), the second is the torsion in degrees
std::pair<int, double>
molecules_container_t::get_torsion(int imol, const std::string &cid, const std::vector<std::string> &atom_names) {
   std::pair<int, double> p(0,0);

   if (is_valid_model_molecule(imol)) {
      p = molecules[imol].get_torsion(cid, atom_names);
   } else {
      std::cout << "WARNING:: " << __FUNCTION__ << "(): not a valid model molecule " << imol << std::endl;
   }
   return p;
}

