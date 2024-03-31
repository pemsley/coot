/*
 * analysis/
 *
 * Copyright 2023 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include <fstream>
#include <glm/glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "kolmogorov.hh"
#include "rna-backbone.hh"

coot::rna_backbone_t::rna_backbone_t(mmdb::Manager *mol, const clipper::Xmap<float> &xmap) {

   auto find_base_residues = [] (mmdb::Manager *mol) {
      std::vector<std::string> base_names = {"C", "G", "U", "A"};
      std::vector<mmdb::Residue *> residues;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  std::string base_name = residue_p->GetResName();
                  if (std::find(base_names.begin(), base_names.end(), base_name) != base_names.end()) {
                     residues.push_back(residue_p);
                  }
               }
            }
         }
      }
      return residues;
   };

   auto sample_map_around_base_atoms = [] (const std::vector<mmdb::Residue *> &base_residues,
                                           const clipper::Xmap<float> &xmap) {
      stats::single ss;
      std::vector<std::string> base_atom_names = { " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " O6 ", " N1 ", " C2 ", " N2 ",
         " N3 ", " C4 ", " O2 ", " N4 " };
      for (unsigned int i=0; i<base_residues.size(); i++) {
         mmdb::Residue *residue_p = base_residues[i];
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         std::vector<mmdb:: Atom *> base_atoms;
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string atom_name(at->GetAtomName());
               if (std::find(base_atom_names.begin(), base_atom_names.end(), atom_name) != base_atom_names.end()) {
                  base_atoms.push_back(at);
               }
            }
         }
         if (base_atoms.size() > 4) {
            std::vector<clipper::Coord_orth> atom_positions;
            for (unsigned int iat=0; iat<base_atoms.size(); iat++) {
               mmdb:: Atom *at = base_atoms[iat];
               atom_positions.push_back(clipper::Coord_orth(at->x, at->y, at->z));
            }
            coot::lsq_plane_info_t lsqp(atom_positions);
            clipper::Coord_orth n = lsqp.normal();
            for (unsigned int iat=0; iat<base_atoms.size(); iat++) {
               mmdb:: Atom *at = base_atoms[iat];
               clipper::Coord_orth at_pos(at->x, at->y, at->z);
               clipper::Coord_orth sample_pos_1(at_pos + 0.6 * n);
               clipper::Coord_orth sample_pos_2(at_pos - 0.6 * n);
               float d1 = util::density_at_point(xmap, sample_pos_1);
               float d2 = util::density_at_point(xmap, sample_pos_2);
               ss.add(d1);
               ss.add(d2);
            }
         }
      }
      double m = ss.mean();
      double v = ss.variance();

      std::cout << "# bases mean: " << m << " std-dev " << std::sqrt(v) << " n_samples: " << ss.size() << std::endl;
      if (false) { // output samples so that I can see what they look like
         for (unsigned int i=0; i<ss.v.size(); i++)
            std::cout << i << " " << ss.v[i] << "\n";
      }
      return ss;
   };

   auto find_tandem_residue_pairs = [] (mmdb::Manager *mol) {
      // the phosphate is on the second of the two residues in the pair
      std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > tandem_residue_pairs;
      int imod = 1;
      std::vector<std::string> base_names = {"C", "G", "U", "A"};
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=1; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  std::string base_name_1 = residue_p->GetResName();
                  if (std::find(base_names.begin(), base_names.end(), base_name_1) != base_names.end()) {
                     int iprev = ires - 1;
                     mmdb::Residue *residue_prev_p = chain_p->GetResidue(iprev);
                     if (residue_prev_p) {
                        std::string base_name_2 = residue_prev_p->GetResName();
                        if (std::find(base_names.begin(), base_names.end(), base_name_2) != base_names.end()) {
                           int resno_r_1 = residue_p->GetSeqNum();
                           int resno_r_2 = residue_prev_p->GetSeqNum();
                           if (resno_r_2 == (resno_r_1 - 1)) {
                              tandem_residue_pairs.push_back(std::make_pair(residue_prev_p, residue_p));
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      std::cout << "found " << tandem_residue_pairs.size() << " tandem residue pairs" << std::endl;
      return tandem_residue_pairs;
   };

   auto find_ribose_centre = [] (mmdb::Residue *residue_p) {
      std::vector<std::string> ribose_atom_names = {" C1'", " C2'", " C3'", " C4'", " O4'"};
      clipper::Coord_orth pos(0,0,0);
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<clipper::Coord_orth> positions;
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());
            if (std::find(ribose_atom_names.begin(), ribose_atom_names.end(), atom_name) != ribose_atom_names.end()) {
               positions.push_back(clipper::Coord_orth(at->x, at->y, at->z));
            }
         }
      }
      if (positions.size() > 4) {
         float inv_size = 1.0/static_cast<float>(positions.size());
         for (unsigned int i=0; i<positions.size(); i++)
            pos += positions[i];
         pos = clipper::Coord_orth(pos.x() * inv_size, pos.y() * inv_size, pos.z() * inv_size);
      }
      return pos;
   };

   auto sample_inner = [] (const clipper::Coord_orth &p1, const clipper::Coord_orth &p2, const clipper::Xmap<float> &xmap,
                           const std::pair<clipper::Coord_orth, clipper::Coord_orth> &first_tube,
                           bool use_first_tube_for_masking) {

      auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
         return glm::vec3(c.x(), c.y(), c.z());
      };

      auto glm_to_clipper = [] (const glm::vec3 &c) {
         return clipper::Coord_orth(c.x, c.y, c.z);
      };

      bool debug_pos = false;
      stats::single ss;
      clipper::Coord_orth vector_1 = p2 - p1;
      glm::vec3 vector_1_glm = clipper_to_glm(vector_1);
      glm::vec3 normalized = glm::normalize(vector_1_glm);
      glm::mat4 ori = glm::orientation(normalized, glm::vec3(0.0, 0.0, 1.0));
      double vector_length = std::sqrt(vector_1.lengthsq());
      // std::cout << "vector_length " << vector_length << std::endl;

      // make sample points on a cylider surface
      int n_stacks = 16; // minimum is 2
      int n_slices = 10;

      double cyl_radius_1 = 2.0;
      double cyl_radius_2 = 0.9;
      // if we are on the second tube, lets' extend a bit to cover the phosphate fully
      int n_start_stack = 0;
      if (use_first_tube_for_masking) n_start_stack = -4;
      for (int i_stack=n_start_stack; i_stack<n_stacks; i_stack++) {
         float vector_fraction = static_cast<float>(i_stack)/static_cast<float>(n_stacks-1);
         for (int i_slice=0; i_slice<n_slices; i_slice++) {
            double rot_angle_deg = 360.0 * static_cast<float>(i_slice)/static_cast<double>(n_slices);
            double ang = clipper::Util::d2rad(rot_angle_deg);
            clipper::Coord_orth p_1(cyl_radius_1 * cos(ang), cyl_radius_1 * sin(ang), vector_length * vector_fraction);
            clipper::Coord_orth p_2(cyl_radius_2 * cos(ang), cyl_radius_2 * sin(ang), vector_length * vector_fraction);
            glm::vec4 p_1_glm = glm::vec4(clipper_to_glm(p_1), 1.0);
            glm::vec4 p_2_glm = glm::vec4(clipper_to_glm(p_2), 1.0);
            glm::vec4 p_1_ori = ori * p_1_glm;  // matrix * position
            glm::vec4 p_2_ori = ori * p_2_glm;
            clipper::Coord_orth sample_point_1 = p1 + glm_to_clipper(glm::vec3(p_1_ori));
            clipper::Coord_orth sample_point_2 = p1 + glm_to_clipper(glm::vec3(p_2_ori));
            bool point_is_inside_other_tube = false;
            if (use_first_tube_for_masking) {
               // set point_is_inside_other_tube if sample_point_1 is inside the first_tube
               clipper::Coord_orth sp_1_reorigin = sample_point_1 - first_tube.first;
               clipper::Coord_orth first_tube_v = first_tube.second - first_tube.first;
               clipper::Coord_orth first_tube_uv(first_tube_v.unit());
               double first_tube_length = std::sqrt((first_tube.second - first_tube.first).lengthsq());
               double dp = clipper::Coord_orth::dot(sp_1_reorigin, first_tube_uv);
               if (dp < first_tube_length) {
                  clipper::Coord_orth projected_point = dp * first_tube_uv;
                  clipper::Coord_orth point_to_projected_point = projected_point - sp_1_reorigin;
                  double point_to_projected_point_length = std::sqrt(point_to_projected_point.lengthsq());
                  if (point_to_projected_point_length < cyl_radius_1)
                     point_is_inside_other_tube = true;
               }
            }

            if (! use_first_tube_for_masking || ! point_is_inside_other_tube) {
               if (debug_pos)
                  std::cout << "pos " << sample_point_1.x() << " " << sample_point_1.y() << " " << sample_point_1.z() << std::endl;
               float d1 = util::density_at_point(xmap, sample_point_1);
               ss.add(d1);
               if (i_slice%2 == 0) {
                  if (debug_pos)
                     std::cout << "pos " << sample_point_2.x() << " " << sample_point_2.y() << " " << sample_point_2.z() << std::endl;
                  float d2 = util::density_at_point(xmap, sample_point_2);
                  ss.add(d2);
               }
            }
         }
      }
      return ss;
   };

   auto tube_sample_inner = [sample_inner] (const clipper::Coord_orth &ribose_centre_1,
                                            const clipper::Coord_orth &P_pos,
                                            const clipper::Coord_orth &ribose_centre_2,
                                            const clipper::Xmap<float> &xmap) {

      std::pair<clipper::Coord_orth, clipper::Coord_orth> first_pair(P_pos, ribose_centre_1);
      stats::single tube_1_samples = sample_inner(P_pos, ribose_centre_1, xmap, first_pair, false);
      stats::single tube_2_samples = sample_inner(P_pos, ribose_centre_2, xmap, first_pair, true);
      stats::single ss;
      ss.add(tube_1_samples);
      ss.add(tube_2_samples);
      return ss;
   };

   auto tube_sample = [find_ribose_centre, tube_sample_inner] (std::pair<mmdb::Residue *, mmdb::Residue *>tandem_residue_pair,
                                                               const clipper::Xmap<float> &xmap) {
      // We are going to sample from the surface of 2 cylinders,
      // with central line vector from the centre of ribose of residue_p to its own phosphate P
      // and from the centre of the ribose of residue_next_p to the phosphate P
      // of residue_p.

      // the phosphate is in the second residue of the pair

      stats::single ss;

      mmdb::Residue *residue_1 = tandem_residue_pair.first;
      mmdb::Residue *residue_2 = tandem_residue_pair.second;
      clipper::Coord_orth ribose_centre_1 = find_ribose_centre(residue_1);
      clipper::Coord_orth ribose_centre_2 = find_ribose_centre(residue_2);

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_2->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());
            if (atom_name == " P  ") {
               clipper::Coord_orth P_pos = co(at);
               ss = tube_sample_inner(ribose_centre_1, P_pos, ribose_centre_2, xmap);
               break;
            }
         }
      }
      return ss;
   };

   auto write_ss_file = [] (const stats::single &ss, mmdb::Residue *residue_p) {
      std::string chain_id = residue_p->GetChainID();
      int res_no = residue_p->GetSeqNum();
      std::string fn = "rna-backbone-test-" + chain_id + std::to_string(res_no) + ".tab";
      std::ofstream f(fn);
      for (unsigned int i=0; i<ss.v.size(); i++) {
         f << i << " " << ss.v[i] << "\n";
      }
      f.close();
   };

   std::vector<std::string> base_names = {"C", "G", "U", "A"};
   std::vector<mmdb::Residue *> base_residues = find_base_residues(mol);
   std::cout << "Found " << base_residues.size() << " base residues" << std::endl;
   stats::single ss_bases = sample_map_around_base_atoms(base_residues, xmap);

   std::vector<std::pair<mmdb::Residue *, mmdb::Residue *> > tandem_residue_pairs = find_tandem_residue_pairs(mol);
   stats::single accum;
   for (unsigned int i=0; i<tandem_residue_pairs.size(); i++) {
      // if (i > 0) continue; // debugging
      stats::single ss_backbone = tube_sample(tandem_residue_pairs[i], xmap);
      if (! ss_backbone.empty()) {
         // write_ss_file(ss_backbone, tandem_residue_pairs[i].first);
         accum.add(ss_backbone);

         std::vector<double> v1(ss_bases.size());
         std::vector<double> v2(ss_backbone.size());
         for (unsigned int j=0; j<ss_bases.v.size();    j++) v1[j] = ss_bases.v[j];
         for (unsigned int j=0; j<ss_backbone.v.size(); j++) v2[j] = ss_backbone.v[j];
         std::vector<double>::const_iterator it_1 = std::min_element(v1.begin(), v1.end());
         std::vector<double>::const_iterator it_2 = std::min_element(v2.begin(), v2.end());
         double min_1 = *it_1;
         double min_2 = *it_2;
         double min_min = min_1;
         if (min_2 < min_1) min_min = min_2;
         // std::cout << "minimum values: " << min_1 << " " << min_2 << std::endl;;
         if (min_1 < 0.0) for (unsigned int j=0; j<v1.size(); j++) v1[j] -= (min_min - 0.01);
         if (min_2 < 0.0) for (unsigned int j=0; j<v2.size(); j++) v2[j] -= (min_min - 0.01);
         if (false) {
            std::cout << "-------- " << v1.size() << " " << v2.size() << std::endl;
            for (unsigned int j=0; j<v1.size(); j++) std::cout << " v1: " << j << " " << v1[j] << std::endl;
            for (unsigned int j=0; j<v2.size(); j++) std::cout << " v2: " << j << " " << v2[j] << std::endl;
         }
         std::pair<double, double> kl_div = nicholls::get_KL(v1, v2);
         std::cout << "   " << coot::residue_spec_t(tandem_residue_pairs[i].first) << " " << kl_div.first << " " << kl_div.second
                   << std::endl;
      }
   }

   double m = accum.mean();
   double v = accum.variance();
   std::cout << "accum mean: " << m << " std-dev " << std::sqrt(v) << " n_samples: " << accum.size() << std::endl;

   if (false) { // output samples so that I can see what they look like
      for (unsigned int i=0; i<accum.v.size(); i++)
         std::cout << i << " " << accum.v[i] << "\n";
   }

}

void coot::rna_backbone_t::scan_all() {

};
