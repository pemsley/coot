/*
 * src/molecular-mesh-generator.cc
 *
 * Copyright 2020 by Medical Research Council
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

#include <chrono>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>

#include <memory>

// GLM
// #include <glm/ext.hpp>
#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp> // for orientation

#ifdef USE_MOLECULES_TO_TRIANGLES
#include <MoleculesToTriangles/CXXClasses/RendererGL.h>
#include <MoleculesToTriangles/CXXClasses/Light.h>
#include <MoleculesToTriangles/CXXClasses/Camera.h>
#include <MoleculesToTriangles/CXXClasses/SceneSetup.h>
#include <MoleculesToTriangles/CXXClasses/ColorScheme.h>
#include <MoleculesToTriangles/CXXClasses/MyMolecule.h>
#include <MoleculesToTriangles/CXXClasses/RepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/VertexColorNormalPrimitive.h>
#endif

#include "coords/Bond_lines.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"

#include "molecular-mesh-generator.hh"

// std::pair<std::vector<glm::mat4>, std::vector<glm::vec4> >make_some_mats_and_colours() {
void
molecular_mesh_generator_t::fill_atom_positions() {

   atom_positions.push_back(glm::vec3(45.677,  -1.080,  18.749));  // N
   atom_positions.push_back(glm::vec3(46.868,  -0.628,  19.509));  // CA
   atom_positions.push_back(glm::vec3(46.627,  -0.827,  20.970));  // CB
   atom_positions.push_back(glm::vec3(47.862,  -0.599,  21.791));  // CG
   atom_positions.push_back(glm::vec3(48.496,  -1.654,  22.429));  // CD2
   atom_positions.push_back(glm::vec3(49.643,  -1.448,  23.153));  // CE2
   atom_positions.push_back(glm::vec3(50.152,  -0.187,  23.312));  // CZ
   atom_positions.push_back(glm::vec3(51.292,   0.018,  24.128));  // OH
   atom_positions.push_back(glm::vec3(49.554,   0.891,  22.730));  // CE1
   atom_positions.push_back(glm::vec3(48.369,   0.691,  21.972));  // CD1
   atom_positions.push_back(glm::vec3(47.970,  -1.584,  19.154));  // C
   atom_positions.push_back(glm::vec3(47.728,  -2.818,  19.139));  // O


   mmdb::Manager *mol = new mmdb::Manager;
   std::string pdb_file_name = "test.pdb";
   // pdb_file_name = "pdb6mbw-pe-build-1.pdb";
   pdb_file_name = "tutorial-modern.pdb";
   pdb_file_name = "twisted-ribbon.pdb";
   pdb_file_name = "";
   if (! pdb_file_name.empty()) {
      mmdb::ERROR_CODE err = mol->ReadPDBASCII(pdb_file_name.c_str());
      if (! err) {
         atom_positions.clear();
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
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
                     if (! at->isTer())
                        atom_positions.push_back(glm::vec3(at->x, at->y, at->z));
                  }
               }
            }
         }
      }
   }
   delete mol;

   glm::vec3 sum(0,0,0);
   for (unsigned int i=0; i<atom_positions.size(); i++)
      sum += atom_positions[i];
   glm::vec3 av = 1.0f/static_cast<float>(atom_positions.size()) * sum;
   for (unsigned int i=0; i<atom_positions.size(); i++)
      atom_positions[i] -= av;
   // for (unsigned int i=0; i<atom_positions.size(); i++)
   // atom_positions[i].x -= 6.9;

   if (false)
      std::cout << "INFO:: molecular_mesh_generator_t::fill_atom_positions() "
                << "move to origin by applying vector " << glm::to_string(av) << std::endl;

}

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::get_test_twisted_trans_peptides() {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

   std::vector<glm::vec3> cis_pep_quad = { glm::vec3(46.538,   9.221,  19.792),     // CA this
                                           glm::vec3(45.726,  10.437,  19.286),     // C this
                                           glm::vec3(46.422,  11.11,   18.459),     // N next
                                           glm::vec3(46.41,   10.53,   17.17) };    // CA next

   // real twisted trans in tutorial-modern.pdb
   if (false)
      cis_pep_quad = {
                      glm::vec3(47.072,  10.589,  17.411),
                      glm::vec3(46.485,  10.875,  16.061),
                      glm::vec3(45.151,  11.014,  16.030),
                      glm::vec3(44.346,  10.357,  15.014) };

   vp = make_twisted_trans_peptide_geom(cis_pep_quad);
   return vp;

}

float
molecular_mesh_generator_t::get_torsion_angle(const std::vector<glm::vec3> &cis_pep_quad) const {

   glm::vec3 v1 = glm::normalize(cis_pep_quad[0] - cis_pep_quad[1]);
   glm::vec3 v2 = glm::normalize(cis_pep_quad[3] - cis_pep_quad[2]);
   float dp = glm::dot(v1, v2);
   if (dp >  1.0) dp =  1.0;
   if (dp < -1.0) dp = -1.0;
   return acos(dp);
}


std::vector<glm::vec3>
molecular_mesh_generator_t::generate_spline_points(std::vector<glm::vec3> &control_points, unsigned int n_steps) {

   // return a vector of size n_steps + 1

   std::vector<glm::vec3> v;
   const glm::vec3 &P1 = control_points[0];
   const glm::vec3 &P2 = control_points[1];
   const glm::vec3 &P3 = control_points[2];
   const glm::vec3 &P4 = control_points[3];
   for (unsigned int i=0; i<n_steps; i++) {
      float t = static_cast<float>(i)/static_cast<float>(n_steps);
      glm::vec3 comp_1 = (1.0f-t)*(1.0f-t)*(1.0f-t)*P1;
      glm::vec3 comp_2 = 3.0f*(1.0f-t)*(1.0f-t)*t*P2;
      glm::vec3 comp_3 = 3.0f*(1.0f-t)*t*t*P3;
      glm::vec3 comp_4 = t*t*t*P4;
      glm::vec3 p = comp_1 + comp_2 + comp_3 + comp_4;
      v.push_back(p);
   }
   v.push_back(P4);

   return v;
}


std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::make_twisted_trans_peptide_geom(const std::vector<glm::vec3> &cis_pep_quad_in) {

   // this shape looks lovely when the angle is about 90 degrees. Looks terrible to the point of
   // deletion for 145 degrees.

   // use the proper glm function - whatever it is...
   auto vec_length = [](const glm::vec3 &v) {
                        float s = v.x * v.x + v.y * v.y + v.z * v.z;
                        return sqrtf(s);
                     };

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;
    // to match the scale in instanced-cylinders.shader, remove
   // in production code.
   const float scale_down_factor = 0.3;

   // std::vector<s_generic_vertex> vv;
   // std::vector<g_triangle> vt;

   std::vector<glm::vec3> cis_pep_quad = cis_pep_quad_in;
   for (unsigned int i=0; i<cis_pep_quad.size(); i++)
      cis_pep_quad[i] -= glm::vec3(41.595322, 8.686544, 12.909259);
   for (unsigned int i=0; i<cis_pep_quad.size(); i++)
      cis_pep_quad[i] *= scale_down_factor; // to match the scale in instanced-cylinders.shader, remove
                              // in production code.

   glm::vec3 C = cis_pep_quad[1];
   glm::vec3 N = cis_pep_quad[2];
   glm::vec3 CA_midpoint  = 0.5f * (cis_pep_quad[0] + cis_pep_quad[3]);
   glm::vec3 C_N_midpoint = 0.5f * (cis_pep_quad[1] + cis_pep_quad[2]);
   glm::vec3 C_N_midpoint_to_CA_midpoint = CA_midpoint - C_N_midpoint;
   glm::vec3 C_N_midpoint_to_CA_midpoint_uv = glm::normalize(C_N_midpoint_to_CA_midpoint);
   glm::vec3 ca_this_uv = glm::normalize(cis_pep_quad[0] - CA_midpoint);
   glm::vec3 ca_next_uv = glm::normalize(cis_pep_quad[3] - CA_midpoint);
   glm::vec3 C_to_CA_midpoint = CA_midpoint - cis_pep_quad[1];
   glm::vec3 N_to_CA_midpoint = CA_midpoint - cis_pep_quad[2];
   glm::vec3 C_to_CA_midpoint_uv = glm::normalize(C_to_CA_midpoint);
   glm::vec3 N_to_CA_midpoint_uv = glm::normalize(N_to_CA_midpoint);
   glm::vec3 C_to_N_uv = glm::normalize(cis_pep_quad[2] - cis_pep_quad[1]);

   float cut_back_frac = 0.01;
   glm::vec3 C_N_midpoint_cutback = C_N_midpoint + cut_back_frac * C_N_midpoint_to_CA_midpoint_uv;
   glm::vec3 CA_this_cutback = cis_pep_quad[0] - cut_back_frac * ca_this_uv;
   glm::vec3 CA_next_cutback = cis_pep_quad[3] - cut_back_frac * ca_next_uv;
   glm::vec3 C_this_cutback  = cis_pep_quad[1] + cut_back_frac * C_to_CA_midpoint_uv;
   glm::vec3 N_next_cutback  = cis_pep_quad[2] + cut_back_frac * N_to_CA_midpoint_uv;

   glm::vec3 CA_midpoint_for_spline = C_N_midpoint_cutback + 1.2f * C_N_midpoint_to_CA_midpoint_uv; // ?

   // find the normals
   glm::vec3 v1 = cis_pep_quad[1] - cis_pep_quad[0];
   glm::vec3 v2 = cis_pep_quad[2] - cis_pep_quad[1]; // cis_pep_quad[1] was cis_pep_quad[0]

   // this normal is close to being what whe want - but not quite.
   glm::vec3 c = glm::cross(v1, v2);
   glm::vec3 normal_1 = glm::normalize(c); // the main normal, perpendicular to C-CA(1) and C-N

   // normal_1A cross product of CA(this)-C and CA-CA-midpoint
   // normal_1B cross product of CA(next)-N and CA-CA-midpoint
   glm::vec3 normal_1A = glm::normalize(glm::cross(cis_pep_quad[0]-cis_pep_quad[1], cis_pep_quad[0]-CA_midpoint_for_spline));
   glm::vec3 normal_1B = glm::normalize(glm::cross(cis_pep_quad[3]-cis_pep_quad[2], cis_pep_quad[3]-CA_midpoint_for_spline));

   v1 = cis_pep_quad[2] - cis_pep_quad[3];
   v2 = cis_pep_quad[1] - cis_pep_quad[3];
   c = glm::cross(v1, v2);
   glm::vec3 normal_2 = glm::normalize(c); // the other main normal, perpendicular to C-CA(2) and C-N

   unsigned int n_spline_points = 28;

   // we need to apply the cutback and offsets *before* the spines, because the splines do *not*
   // take the same path, with just an offset.
   //
   glm::vec3 displacement = 0.07f * normal_1;
   glm::vec3 d = displacement;


   float bl_CN = vec_length(C_this_cutback-N_next_cutback);

   // what is the angle between normal_1 and normal_2?
   float angle = 1.4; // radians
   angle = get_torsion_angle(cis_pep_quad_in); // do I want 2pi - this?

   //std::cout << "Here with torsion angle " << angle << " " << angle * 180.0/M_PI << " degrees"  << std::endl;
   // angle *= 1.7;

   // the last points around the C-N bond, for j = 0 and j = 1:
   glm::vec3 d_base = displacement;
   glm::vec3 p_2_C_end_O_A = C_this_cutback + bl_CN * C_to_N_uv + d_base;
   glm::vec3 p_2A = p_2_C_end_O_A - C;
   glm::vec3 p_2B = glm::rotate(p_2A, angle, C_to_N_uv);
   glm::vec3 p_2_C_end_0 = p_2B + C;

   glm::vec3 p_2_C_end_1_A = C_this_cutback + bl_CN * C_to_N_uv - d_base;
   p_2A = p_2_C_end_1_A - C;
   p_2B = glm::rotate(p_2A, angle, C_to_N_uv);
   glm::vec3 p_2_C_end_1 = p_2B + C;

   // Now use those 2 points to find how the displacement has rotated, so that we can use
   // it for boxing up.
   glm::vec3 displacement_rotated = 0.5f * (p_2_C_end_1 - p_2_C_end_0);

   float spline_control_point_normal_scale_factor = 0.15;
   float ds = 0.01;

   glm::vec3 cp_1_1 = 0.5f * (cis_pep_quad[0] + CA_midpoint) + normal_1A * spline_control_point_normal_scale_factor - ds * displacement;
   glm::vec3 cp_1_2 = 0.5f * (cis_pep_quad[3] + CA_midpoint) + normal_1B * spline_control_point_normal_scale_factor - ds * displacement_rotated;

   glm::vec3 cp_2_1 = 0.5f * (cis_pep_quad[0] + CA_midpoint) + normal_1A * spline_control_point_normal_scale_factor + ds * displacement;
   glm::vec3 cp_2_2 = 0.5f * (cis_pep_quad[3] + CA_midpoint) + normal_1B * spline_control_point_normal_scale_factor + ds * displacement_rotated;

   std::vector<glm::vec3> spline_quad_1 = { CA_this_cutback + displacement, cp_1_1, cp_1_2, CA_next_cutback - displacement_rotated};
   std::vector<glm::vec3> spline_quad_2 = { CA_this_cutback - displacement, cp_2_1, cp_2_2, CA_next_cutback + displacement_rotated };
   std::vector<glm::vec3> sp_1 = generate_spline_points(spline_quad_1, n_spline_points);
   std::vector<glm::vec3> sp_2 = generate_spline_points(spline_quad_2, n_spline_points);

   float dodec_scale = 0.1;
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > ca_mid_dodec = dodec_at_position(dodec_scale, CA_midpoint);
   for (unsigned int i=0; i<ca_mid_dodec.first.size(); i++) ca_mid_dodec.first[i].color = glm::vec4(0.3,0.3,0.3,1);

   for (unsigned int i=0; i<16; i++) {
      float f = static_cast<float>(i+1)/static_cast<float>(16);
      glm::vec3 dd_a = static_cast<float>(i) * 0.5f * displacement;
      glm::vec3 dd_b = static_cast<float>(i) * 0.5f * displacement_rotated;
      glm::vec3 t_cp_1_1 = 0.5f * (cis_pep_quad[0] + CA_midpoint) + normal_1A * spline_control_point_normal_scale_factor - dd_a;
      glm::vec3 t_cp_1_2 = 0.5f * (cis_pep_quad[3] + CA_midpoint) + normal_1B * spline_control_point_normal_scale_factor + dd_b;
      glm::vec3 c_cp_1_1 = 0.5f * (cis_pep_quad[0] + CA_midpoint) + normal_1A * spline_control_point_normal_scale_factor * f - displacement;
      glm::vec3 c_cp_1_2 = 0.5f * (cis_pep_quad[3] + CA_midpoint) + normal_1B * spline_control_point_normal_scale_factor * f + displacement_rotated;
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > m11 = dodec_at_position(dodec_scale, t_cp_1_1);
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > m12 = dodec_at_position(dodec_scale, t_cp_1_2);
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > c11 = dodec_at_position(dodec_scale, c_cp_1_1);
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > c12 = dodec_at_position(dodec_scale, c_cp_1_2);
      for (unsigned int i=0; i<m11.first.size(); i++) m11.first[i].color = glm::vec4(0.4,0.2,0.2,1);
      for (unsigned int i=0; i<m11.first.size(); i++) m12.first[i].color = glm::vec4(0.3,0.5,0.2,1);
      add_to_mesh(&vp, m11);
      add_to_mesh(&vp, m12);
      add_to_mesh(&vp, c11);
      add_to_mesh(&vp, c12);
   }

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > cp_dodec_1_1 = dodec_at_position(dodec_scale, cp_1_1);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > cp_dodec_1_2 = dodec_at_position(dodec_scale, cp_1_2);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > cp_dodec_2_1 = dodec_at_position(dodec_scale, cp_2_1);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > cp_dodec_2_2 = dodec_at_position(dodec_scale, cp_2_2);
   for (unsigned int i=0; i<cp_dodec_1_1.first.size(); i++) cp_dodec_1_1.first[i].color = glm::vec4(1,0,0,1);
   for (unsigned int i=0; i<cp_dodec_1_1.first.size(); i++) cp_dodec_1_2.first[i].color = glm::vec4(0.6,0.6,0,1);
   for (unsigned int i=0; i<cp_dodec_1_1.first.size(); i++) cp_dodec_2_1.first[i].color = glm::vec4(0,0,1,1);
   for (unsigned int i=0; i<cp_dodec_1_1.first.size(); i++) cp_dodec_2_2.first[i].color = glm::vec4(0.6,0.4,0.8,1);

   add_to_mesh(&vp, ca_mid_dodec);
   add_to_mesh(&vp, cp_dodec_1_1);
   add_to_mesh(&vp, cp_dodec_1_2);
   add_to_mesh(&vp, cp_dodec_2_1);
   add_to_mesh(&vp, cp_dodec_2_2);

   glm::vec3 n = normal_1;
   glm::vec4 col(0.7f, 0.7f, 0.3f, 1.0f); // ignored, because we use material


   unsigned int idx_base_spline_begin = vp.first.size();
   glm::vec3 displacement_base = displacement;

   for (unsigned int j=0; j<2; j++) {
      if (j==1) continue; // debugging
      for (unsigned int i=0; i<n_spline_points; i++) {
         float f_1 = static_cast<float>(i)  /static_cast<float>(n_spline_points);
         float f_2 = static_cast<float>(i+1)/static_cast<float>(n_spline_points);
         // rotate displacement_base around the C-N bond to generate displacements
         glm::vec3 d_base = displacement;
         if (j == 1) d_base = -displacement;

         float a_1 = f_1 * angle;
         float a_2 = f_2 * angle;
         glm::vec3 p_1 = C_this_cutback + f_1 * bl_CN * C_to_N_uv + d_base;
         glm::vec3 p_1A = p_1 - C;
         glm::vec3 p_1B = glm::rotate(p_1A, a_1, C_to_N_uv);
         glm::vec3 p_1C = p_1B + C;

         glm::vec3 p_2 = C_this_cutback + f_2 * bl_CN * C_to_N_uv + d_base;
         glm::vec3 p_2A = p_2 - C;
         glm::vec3 p_2B = glm::rotate(p_2A, a_2, C_to_N_uv);
         glm::vec3 p_2C = p_2B + C;

         unsigned int idx_base = vp.first.size();
         glm::vec3 v1 = p_1C;
         glm::vec3 v2 = p_2C;
         glm::vec3 v3 = sp_1[i];
         glm::vec3 v4 = sp_1[i+1];
         if (j == 1) {
            v3 = sp_2[i];
            v4 = sp_2[i+1];
         }
         glm::vec3 c_1 = glm::cross(v1-v3, v2-v1);
         glm::vec3 n_1 = glm::normalize(c_1);
         glm::vec3 c_2 = glm::cross(v4-v3, v4-v2); // sign/order swap?
         glm::vec3 n_2 = glm::normalize(c_2);
         if (j == 1) {
            n_1 = -n_1;
            n_2 = -n_2;
         } 
         // this quad is not flat, so we add 2 triangles, so that, after this loop, each
         // there will be 4 representations of this position, each with its own normal
         // (and hence, vertex). Each vertex is part of a different triangle.
         // In the middle of the spine these normals are very different.
         vp.first.push_back(s_generic_vertex(v1, n_1, col));
         vp.first.push_back(s_generic_vertex(v2, n_1, col));
         vp.first.push_back(s_generic_vertex(v3, n_1, col));
         vp.first.push_back(s_generic_vertex(v2, n_2, col));
         vp.first.push_back(s_generic_vertex(v3, n_2, col));
         vp.first.push_back(s_generic_vertex(v4, n_2, col));
         vp.second.push_back(g_triangle(idx_base+0, idx_base+1, idx_base+2));
         vp.second.push_back(g_triangle(idx_base+3, idx_base+4, idx_base+5));
      }
   }


   // now smooth the normals of those vertices
   unsigned int idx_base_spline_end = vp.first.size();
   smooth_vertices(&vp.first, idx_base_spline_begin, idx_base_spline_end);

   // ...

   // now box up the edges

   unsigned int idx_base = vp.first.size();


   // box up the sides

   {
      // CA-C box
      idx_base = vp.first.size();
      n = -glm::normalize(glm::cross(displacement, C_this_cutback-CA_this_cutback));
      vp.first.push_back(s_generic_vertex( C_this_cutback+displacement, n, col));
      vp.first.push_back(s_generic_vertex (C_this_cutback-displacement, n, col));
      vp.first.push_back(s_generic_vertex(CA_this_cutback+displacement, n, col));
      vp.first.push_back(s_generic_vertex(CA_this_cutback-displacement, n, col));
      vp.second.push_back(g_triangle(idx_base+0, idx_base+1, idx_base+2));
      vp.second.push_back(g_triangle(idx_base+1, idx_base+3, idx_base+2));
   }

   {
      // N-CA box

      idx_base = vp.first.size();

      // as in the boxing up the C-N bond loop:
      glm::vec3 d_base = displacement;
      float f_2 = 1.0f;
      glm::vec3 p_2 = C_this_cutback + f_2 * bl_CN * C_to_N_uv + d_base;
      glm::vec3 p_2A = p_2 - C;
      glm::vec3 p_2B = glm::rotate(p_2A, angle, C_to_N_uv);
      glm::vec3 p_2C = p_2B + C;

      glm::vec3 p_4 = C_this_cutback + f_2 * bl_CN * C_to_N_uv - d_base;
      glm::vec3 p_4A = p_4 - C;
      glm::vec3 p_4B = glm::rotate(p_4A, angle, C_to_N_uv);
      glm::vec3 p_4C = p_4B + C;

      // v1 and v2 are at the N.
      // v3 and v4 are the end control points for the splines.
      glm::vec3 v1 =  p_2C;
      glm::vec3 v2 =  p_4C;
      glm::vec3 v3 = CA_next_cutback - displacement_rotated;
      glm::vec3 v4 = CA_next_cutback + displacement_rotated;

      glm::vec3 c_nca = glm::cross(v1-v3, v2-v1);
      glm::vec3 n_nca = glm::normalize(c_nca);

      vp.first.push_back(s_generic_vertex(v1, n_nca, col));
      vp.first.push_back(s_generic_vertex(v2, n_nca, col));
      vp.first.push_back(s_generic_vertex(v3, n_nca, col));
      vp.first.push_back(s_generic_vertex(v4, n_nca, col));
      vp.second.push_back(g_triangle(idx_base+0, idx_base+1, idx_base+2));
      vp.second.push_back(g_triangle(idx_base+1, idx_base+3, idx_base+2));
   }


   // box up the points twisting around the C-N bond
   idx_base_spline_begin = vp.first.size();
   for (unsigned int i=0; i<n_spline_points; i++) {

      float f_1 = static_cast<float>(i)  /static_cast<float>(n_spline_points);
      float f_2 = static_cast<float>(i+1)/static_cast<float>(n_spline_points);
      // rotate displacement_base around the C-N bond to generate displacements
      glm::vec3 d_base = displacement;

      float a_1 = f_1 * angle;
      float a_2 = f_2 * angle;
      glm::vec3 p_1 = C_this_cutback + f_1 * bl_CN * C_to_N_uv + d_base;
      glm::vec3 p_1A = p_1 - C;
      glm::vec3 p_1B = glm::rotate(p_1A, a_1, C_to_N_uv);
      glm::vec3 p_1C = p_1B + C;

      glm::vec3 p_2 = C_this_cutback + f_2 * bl_CN * C_to_N_uv + d_base;
      glm::vec3 p_2A = p_2 - C;
      glm::vec3 p_2B = glm::rotate(p_2A, a_2, C_to_N_uv);
      glm::vec3 p_2C = p_2B + C;

      glm::vec3 p_3 = C_this_cutback + f_1 * bl_CN * C_to_N_uv - d_base;
      glm::vec3 p_3A = p_3 - C;
      glm::vec3 p_3B = glm::rotate(p_3A, a_1, C_to_N_uv);
      glm::vec3 p_3C = p_3B + C;

      glm::vec3 p_4 = C_this_cutback + f_2 * bl_CN * C_to_N_uv - d_base;
      glm::vec3 p_4A = p_4 - C;
      glm::vec3 p_4B = glm::rotate(p_4A, a_2, C_to_N_uv);
      glm::vec3 p_4C = p_4B + C;
      
      idx_base = vp.first.size();
      c = glm::cross(p_2C-p_1C, p_2C-p_3C);
      n = glm::normalize(c);
      vp.first.push_back(s_generic_vertex(p_1C, n, col));
      vp.first.push_back(s_generic_vertex(p_2C, n, col));
      vp.first.push_back(s_generic_vertex(p_3C, n, col));
      vp.first.push_back(s_generic_vertex(p_4C, n, col));
      vp.first.push_back(s_generic_vertex(CA_next_cutback - displacement_rotated, n, col));
      vp.first.push_back(s_generic_vertex(CA_next_cutback + displacement_rotated, n, col));
      vp.second.push_back(g_triangle(idx_base+0, idx_base+1, idx_base+2));
      vp.second.push_back(g_triangle(idx_base+1, idx_base+3, idx_base+2));
   }
   // now smooth the normals of those vertices
   idx_base_spline_end = vp.first.size();
   smooth_vertices(&vp.first, idx_base_spline_begin, idx_base_spline_end);

   
   // box up the splined surface
   if (true) {
      idx_base_spline_begin = vp.first.size();
      for (unsigned int i=0; i<n_spline_points; i++) {

         glm::vec3 p_1 = sp_1[i];
         glm::vec3 p_2 = sp_1[i+1];
         glm::vec3 p_3 = sp_2[i];
         glm::vec3 p_4 = sp_2[i+1];

         idx_base = vp.first.size();
         c = glm::cross(p_1-p_2, p_2-p_3);
         n = glm::normalize(c);
         vp.first.push_back(s_generic_vertex(p_1, n, col));
         vp.first.push_back(s_generic_vertex(p_2, n, col));
         vp.first.push_back(s_generic_vertex(p_3, n, col));
         vp.first.push_back(s_generic_vertex(p_4, n, col));
         vp.second.push_back(g_triangle(idx_base+0, idx_base+1, idx_base+2));
         vp.second.push_back(g_triangle(idx_base+1, idx_base+3, idx_base+2));
      }
      idx_base_spline_end = vp.first.size();
      smooth_vertices(&vp.first, idx_base_spline_begin, idx_base_spline_end);
   }

   return vp;
}

void
molecular_mesh_generator_t::smooth_vertices(std::vector<s_generic_vertex> *v_p, unsigned int idx_begin, unsigned int idx_end) {

   float min_angle = 90.0; // degrees
   float cos_min = cosf(min_angle * 2.0 * M_PI/360.0);

   std::vector<s_generic_vertex> &verts = *v_p;
   std::map<unsigned int, std::set<unsigned int> > m;
   for (unsigned int i=idx_begin; i<idx_end; i++) {
      for (unsigned int j=i; j<idx_end; j++) {
         if (i != j) {
            if (verts[i].pos == verts[j].pos) {
               m[i].insert(j);
               m[j].insert(i);
            }
         }
      }
   }

   std::map<unsigned int, glm::vec3> updated_normal_map;

   std::map<unsigned int, std::set<unsigned int> >::const_iterator it;
   for (it=m.begin(); it!=m.end(); it++) {
      glm::vec3 sum = verts[it->first].normal;
      glm::vec3 base_normal = verts[it->first].normal;
      unsigned int count = 1;
      std::set<unsigned int>::const_iterator it_s;
      for (it_s=it->second.begin(); it_s!=it->second.end(); it_s++) {
         float dp = glm::dot(base_normal, verts[*it_s].normal);
         if (dp > cos_min) {
            sum += verts[*it_s].normal;
            count += 1;
         }
      }
      float sc = 1.0f / static_cast<float>(count);
      glm::vec3 av = sc * sum;
      updated_normal_map[it->first] = av;
   }

   std::map<unsigned int, glm::vec3>::const_iterator it_udm;
   for (it_udm=updated_normal_map.begin(); it_udm!=updated_normal_map.end(); ++it_udm) {
      verts[it_udm->first].normal = glm::normalize(it_udm->second);
   }

}


std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::get_test_cis_peptides() { // maybe should be in graphical_molecule

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

#ifdef THIS_IS_HMT
#else

   std::vector<glm::vec3> cis_pep_quad = { glm::vec3(46.538,   9.221,  19.792),     // CA this
                                           glm::vec3(45.726,  10.437,  19.286),     // C this
                                           glm::vec3(46.374,  11.225,  18.462),     // N next
                                           glm::vec3(47.072,  10.589,  17.411) };   // CA next

   coot::util::cis_peptide_quad_info_t::type_t type = coot::util::cis_peptide_quad_info_t::CIS;
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp_cis_pep =
      make_cis_peptide_geom(cis_pep_quad, type);

   unsigned int idx_base = vp.first.size();
   unsigned int idx_tri_base = vp.second.size();

   vp.first.insert(vp.first.end(), vp_cis_pep.first.begin(), vp_cis_pep.first.end());
   vp.second.insert(vp.second.end(), vp_cis_pep.second.begin(), vp_cis_pep.second.end());;
   for (unsigned int i=idx_tri_base; i<vp.second.size(); i++)
      vp.second[i].rebase(idx_base);
#endif

   return vp;

}

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::get_cis_peptides(const std::string &pdb_file_name) { // maybe should be in graphical_molecule

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

   mmdb::Manager *mol = new mmdb::Manager;
   mmdb::ERROR_CODE err = mol->ReadPDBASCII(pdb_file_name.c_str());
   if (! err) {
      // meh, torsion code - I don't want to do that at the moment.
   }

   vp = get_test_cis_peptides();

   for (unsigned int i=0; i<vp.first.size(); i++) {
      vp.first[i].pos.x -= 1.2f;
      vp.first[i].pos.z -= 1.2f;
   }

   unsigned int idx_base = vp.first.size();
   unsigned int idx_tri_base = vp.second.size();

   if (false) {
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp_twisted_trans = get_test_twisted_trans_peptides();
      vp.first.insert(vp.first.end(), vp_twisted_trans.first.begin(), vp_twisted_trans.first.end());
      vp.second.insert(vp.second.end(), vp_twisted_trans.second.begin(), vp_twisted_trans.second.end());;
      for (unsigned int i=idx_tri_base; i<vp.second.size(); i++)
         vp.second[i].rebase(idx_base);
   }

   if (false) {
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp_cis_pep = get_test_cis_peptides();
      vp.first.insert(vp.first.end(), vp_cis_pep.first.begin(), vp_cis_pep.first.end());
      vp.second.insert(vp.second.end(), vp_cis_pep.second.begin(), vp_cis_pep.second.end());;
      for (unsigned int i=idx_tri_base; i<vp.second.size(); i++)
         vp.second[i].rebase(idx_base);
   }

   return vp;
}


#ifdef THIS_IS_HMT
#else
#include "coot-utils/coot-coord-utils.hh"

// return a map with the key as the model number (normally only 1 of course)
// the first of the pair is the type (cis, pre-PRO-cis and twisted-trans).
std::map<int, std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > >
molecular_mesh_generator_t::make_cis_peptide_quads_mesh(mmdb::Manager *mol) {

   // this is old code. Probably should be deleted.

   std::map<int, std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > > cis_pep_geometry_map;
   if (! mol) std::cout << "error:: in molecular_mesh_generator_t::make_cis_peptide_quads_mesh() null mol" << std::endl;
   if (! mol) return cis_pep_geometry_map;

   auto atom_to_glm = [] (mmdb::Atom *at) { return glm::vec3(at->x, at->y, at->z); };

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      std::vector<std::vector<glm::vec3> > v_cis_peps;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         bool strictly_cis = false;
         std::vector<coot::util::cis_peptide_quad_info_t> cis_peps;
         // coot::util::cis_peptide_quads_from_coords(mol, imod, strictly_cis);

         std::vector<s_generic_vertex> sum_of_vertices;
         std::vector<g_triangle> sum_of_triangles;

         for (auto cp : cis_peps) {
            const coot::atom_quad &quad = cp.quad;
            glm::vec3 p1 = atom_to_glm(quad.atom_1);
            glm::vec3 p2 = atom_to_glm(quad.atom_2);
            glm::vec3 p3 = atom_to_glm(quad.atom_3);
            glm::vec3 p4 = atom_to_glm(quad.atom_4);
            std::vector<glm::vec3> cp_atom_positions = { p1, p2, p3, p4 };
            std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
               cis_pep_geometry = make_cis_peptide_geom(cp_atom_positions, cp.type);
            add_to_mesh(&(cis_pep_geometry_map[imod]), cis_pep_geometry);
         }
      }
   }
   return  cis_pep_geometry_map;
}
#endif // THIS_IS_HMT


#ifdef THIS_IS_HMT
#else

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::make_cis_peptide_geom(const std::vector<glm::vec3> &cis_pep_quad_in,
                                                  coot::util::cis_peptide_quad_info_t::type_t type) {

   // type: 9 unset
   // type: 1 cis
   // type: 2 pre-PRO cis
   // type: 3 twisted-trans

   const   std::vector<glm::vec3> &cis_pep_quad = cis_pep_quad_in;


   // find the normal
   glm::vec3 v1 = cis_pep_quad[1] - cis_pep_quad[0];
   glm::vec3 v2 = cis_pep_quad[2] - cis_pep_quad[3];

   glm::vec3 c = glm::cross(v1, v2);
   glm::vec3 normal_1 = glm::normalize(c); // the main normal

   std::vector<s_generic_vertex> vv;
   std::vector<g_triangle> vt;

   glm::vec3 CA_midpoint  = 0.5f * (cis_pep_quad[0] + cis_pep_quad[3]);
   glm::vec3 C_N_midpoint = 0.5f * (cis_pep_quad[1] + cis_pep_quad[2]);
   glm::vec3 C_N_midpoint_to_CA_midpoint = CA_midpoint - C_N_midpoint;
   glm::vec3 C_N_midpoint_to_CA_midpoint_uv = glm::normalize(C_N_midpoint_to_CA_midpoint);
   glm::vec3 ca_this_uv = glm::normalize(cis_pep_quad[0] - CA_midpoint);
   glm::vec3 ca_next_uv = glm::normalize(cis_pep_quad[3] - CA_midpoint);
   glm::vec3 C_to_CA_midpoint = CA_midpoint - cis_pep_quad[1];
   glm::vec3 N_to_CA_midpoint = CA_midpoint - cis_pep_quad[2];
   glm::vec3 C_to_CA_midpoint_uv = glm::normalize(C_to_CA_midpoint);
   glm::vec3 N_to_CA_midpoint_uv = glm::normalize(N_to_CA_midpoint);

   float cut_back_frac = 0.2; // was 0.05
   glm::vec3 C_N_midpoint_cutback = C_N_midpoint + cut_back_frac *  C_N_midpoint_to_CA_midpoint_uv;
   glm::vec3 CA_this_cutback = cis_pep_quad[0] - cut_back_frac * ca_this_uv;
   glm::vec3 CA_next_cutback = cis_pep_quad[3] - cut_back_frac * ca_next_uv;
   glm::vec3 C_this_cutback  = cis_pep_quad[1] + cut_back_frac * C_to_CA_midpoint_uv;
   glm::vec3 N_next_cutback  = cis_pep_quad[2] + cut_back_frac * N_to_CA_midpoint_uv;

   glm::vec3 d = 3.0f * 0.02f * normal_1;

   glm::vec4 col(0.6f, 0.2f, 0.2f, 1.0f);
   if (type == coot::util::cis_peptide_quad_info_t::PRE_PRO_CIS)
      col = glm::vec4(0.1f, 0.55f, 0.1f, 1.0f);
   if (type == coot::util::cis_peptide_quad_info_t::TWISTED_TRANS)
      col = glm::vec4(0.6f, 0.6f, 0.3f, 1.0f);

   glm::vec3 normal_1A = glm::normalize(glm::cross(cis_pep_quad[0]-cis_pep_quad[1], cis_pep_quad[0]-CA_midpoint));
   glm::vec3 normal_1B = glm::normalize(glm::cross(cis_pep_quad[3]-cis_pep_quad[2], cis_pep_quad[3]-CA_midpoint));

   glm::vec3 displacement = d;
   glm::vec3 n = normal_1;

   vv.push_back(s_generic_vertex(CA_midpoint+displacement, n, col));
   vv.push_back(s_generic_vertex(C_N_midpoint_cutback+displacement, n, col));
   vv.push_back(s_generic_vertex(CA_this_cutback + displacement, -normal_1A, col));
   vv.push_back(s_generic_vertex(CA_next_cutback + displacement,  normal_1B, col));
   vv.push_back(s_generic_vertex( C_this_cutback + displacement, -normal_1A, col));
   vv.push_back(s_generic_vertex( N_next_cutback + displacement,  normal_1B, col));

   displacement = -d;
   n = -normal_1;
   vv.push_back(s_generic_vertex(CA_midpoint+displacement, n, col));
   vv.push_back(s_generic_vertex(C_N_midpoint_cutback+displacement, n, col));
   vv.push_back(s_generic_vertex(CA_this_cutback + displacement,  normal_1A, col));
   vv.push_back(s_generic_vertex(CA_next_cutback + displacement, -normal_1B, col));
   vv.push_back(s_generic_vertex(C_this_cutback + displacement,   normal_1A, col));
   vv.push_back(s_generic_vertex(N_next_cutback + displacement,  -normal_1B, col));

   vt.push_back(g_triangle(2,0,4)); // 20230215-PE reverse the winding.
   vt.push_back(g_triangle(4,0,5));
   vt.push_back(g_triangle(0,3,5));

   vt.push_back(g_triangle(8,10,6));
   vt.push_back(g_triangle(10,11,6));
   vt.push_back(g_triangle(6,11,9));

   // now the faces of the box, noting that the faces have different vertices with the same position
   // because the normals are facing in different directions.

   n = glm::normalize(CA_midpoint - C_N_midpoint);
   unsigned int idx_base = vv.size();
   vv.push_back(s_generic_vertex(CA_this_cutback + d, n, col));
   vv.push_back(s_generic_vertex(CA_this_cutback - d, n, col));
   vv.push_back(s_generic_vertex(CA_next_cutback + d, n, col));
   vv.push_back(s_generic_vertex(CA_next_cutback - d, n, col));

   vt.push_back(g_triangle(idx_base, idx_base+1, idx_base+2));
   vt.push_back(g_triangle(idx_base+1, idx_base+3, idx_base+2));

   idx_base = vv.size();
   n = glm::normalize(0.5f * (CA_this_cutback + C_this_cutback) - CA_next_cutback);
   vv.push_back(s_generic_vertex(CA_this_cutback + d, n, col));
   vv.push_back(s_generic_vertex(CA_this_cutback - d, n, col));
   vv.push_back(s_generic_vertex(C_this_cutback + d, n, col));
   vv.push_back(s_generic_vertex(C_this_cutback - d, n, col));

   vt.push_back(g_triangle(idx_base, idx_base+1, idx_base+2));
   vt.push_back(g_triangle(idx_base+1, idx_base+3, idx_base+2));

   idx_base = vv.size();
   n = glm::normalize(C_N_midpoint - CA_midpoint);
   vv.push_back(s_generic_vertex(C_this_cutback + d, n, col));
   vv.push_back(s_generic_vertex(C_this_cutback - d, n, col));
   vv.push_back(s_generic_vertex(N_next_cutback + d, n, col));
   vv.push_back(s_generic_vertex(N_next_cutback - d, n, col));
   vt.push_back(g_triangle(idx_base, idx_base+1, idx_base+2));
   vt.push_back(g_triangle(idx_base+1, idx_base+3, idx_base+2));

   idx_base = vv.size();
   n = glm::normalize(0.5f * (CA_next_cutback + N_next_cutback) - CA_this_cutback);
   vv.push_back(s_generic_vertex(N_next_cutback + d, n, col));
   vv.push_back(s_generic_vertex(N_next_cutback - d, n, col));
   vv.push_back(s_generic_vertex(CA_next_cutback + d, n, col));
   vv.push_back(s_generic_vertex(CA_next_cutback - d, n, col));
   vt.push_back(g_triangle(idx_base, idx_base+1, idx_base+2));
   vt.push_back(g_triangle(idx_base+1, idx_base+3, idx_base+2));

   // std::move here? Or is this optimized?
   return std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >(vv, vt);

}
#endif // THIS_IS_HMT

void
molecular_mesh_generator_t::update_mats_and_colours() {


   // I want to generate a set of instanced bond matrices to apply to the cylinder mesh
   // I also want to generate a set of instanced atom matrices to apply to the atoms with only one
   // bond (for now).

   bool do_GM_restraints = false;

   // what is the glm function that does this? - it's glm::distance() not glm::length()
   auto vec_length = [](const glm::vec3 &v) {
                        float s = v.x * v.x + v.y * v.y + v.z * v.z;
                        return sqrtf(s);
                     };

   const std::vector<glm::vec3> &atom_positions_local = atom_positions;

   // If they haven't been made, the make them, else update them

   auto get_bond_matrix = [atom_positions_local, vec_length] (const unsigned int &i, const unsigned int &j, float radius) {
                        glm::vec3 delta = atom_positions_local[j] - atom_positions_local[i];
                        float l_delta = vec_length(delta);
                        glm::mat4 u(1.0f);
                        glm::mat4 sc = glm::scale(u, glm::vec3(radius, radius, l_delta));
                        // orient
                        glm::vec3 normalized_bond_orientation(glm::normalize(delta));
                        glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, 1.0)); // nice
                        // translate
                        glm::mat4 t = glm::translate(u, atom_positions_local[i]);
                        glm::mat4 m = t * ori * sc;
                        return m;
                     };

   // like above, different z axis because we want the hemisphere to extend outside the cylinder - and we don't need to
   // scale to bond length
   auto get_octahemi_matrix = [atom_positions_local] (const unsigned int &i, const unsigned int &j, float radius) {
                                 glm::vec3 delta = atom_positions_local[j] - atom_positions_local[i];
                                 glm::mat4 u(1.0f);
                                 glm::mat4 sc = glm::scale(u, glm::vec3(radius, radius, radius));
                                 // orient
                                 glm::vec3 normalized_bond_orientation(glm::normalize(delta));
                                 // I can use just delta here, I think? Maybe above too? Test later
                                 glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, -1.0));
                                 // translate
                                 glm::mat4 t = glm::translate(u, atom_positions_local[i]);
                                 glm::mat4 m = t * ori * sc;
                                 return m;
                              };

   auto delta_length_to_colour = [] (const float &delta_l) {
                                    float d = delta_l;
                                    if (d >  1.0) d =  1.0;
                                    if (d < -1.0) d = -1.0;
                                    return glm::vec4(0.5 + d * 0.5, 0.5 - d * 0.5, 0.5 + d * 0.5, 1.0);
                                 };

   glm::vec4 atom_colour(0.2f, 0.6f, 0.3f, 1.0f);

   if (bond_atom_pair_indices.empty()) {
      // first (0th) frame
      bond_instance_colours.clear(); // they should be empty
      bond_instance_matrices.clear();
      std::map<unsigned int, std::vector<unsigned int> > atom_bond_count_map; // for every atom count the number of bonds attached
      for (unsigned int i=0; i<atom_positions.size(); i++) {
         for (unsigned int j=i; j<atom_positions.size(); j++) {
            if (j > i) {
               glm::vec3 delta = atom_positions[j] - atom_positions[i];
               float bond_length = vec_length(delta);
               if (bond_length > 1.7) {
                  if (bond_length < 4.0) {
                     float radius = 0.02;
                     if (do_GM_restraints) {
                        bond_instance_matrices.push_back(get_bond_matrix(i,j, radius));
                        bond_instance_colours.push_back(delta_length_to_colour(0.0f));
                        original_bond_lengths.push_back(bond_length);
                        unsigned int idx = bond_instance_matrices.size() - 1;
                        bond_atom_pair_indices[idx] = std::tuple<unsigned int, unsigned int, float, bool>(i,j,radius,false);
                     }
                  }
               } else {
                  float radius = 0.1;
                  bond_instance_matrices.push_back(get_bond_matrix(i,j, radius));
                  bond_instance_colours.push_back(atom_colour);
                  original_bond_lengths.push_back(bond_length);
                  unsigned int idx = bond_instance_matrices.size() - 1;
                  bond_atom_pair_indices[idx] = std::tuple<unsigned int, unsigned int, float, bool>(i,j,radius,true);
                  atom_bond_count_map[i].push_back(j);
                  atom_bond_count_map[j].push_back(i);
               }
            }
         }
      }

      // now the atoms, draw a hemisphere with atoms with 1 bond
      std::map<unsigned int, std::vector<unsigned int> >::const_iterator it;
      for (it=atom_bond_count_map.begin(); it!=atom_bond_count_map.end(); it++) {
         float radius = 0.1;
         unsigned int idx_this = it->first;
         unsigned int idx_other = it->second[0];
         glm::mat4 m = get_octahemi_matrix(idx_this, idx_other, radius);
         atom_instance_matrices.push_back(m);
         atom_instance_colours.push_back(atom_colour);
         atom_atom_pair_indices.push_back(std::tuple<unsigned int, unsigned int, float>(idx_this, idx_other, radius));
      }

   } else {
      // every other frame

      // bonds
      std::map<unsigned int, std::tuple<unsigned int, unsigned int, float, bool> >::const_iterator it;
      for (it=bond_atom_pair_indices.begin(); it!=bond_atom_pair_indices.end(); it++) {
         const unsigned int idx = it->first;
         const std::tuple<unsigned int, unsigned int, float, bool> &p = it->second;
         bond_instance_matrices[idx] = get_bond_matrix(std::get<0>(p), std::get<1>(p), std::get<2>(p));
         if (std::get<3>(p)) {
            bond_instance_colours[idx] = atom_colour;
         } else {
            glm::vec3 delta = atom_positions[std::get<0>(p)] - atom_positions[std::get<1>(p)];
            float l = vec_length(delta);
            float delta_l = l - original_bond_lengths[idx];
            bond_instance_colours[idx] = delta_length_to_colour(delta_l);
         }
      }

      // hemisphere atoms
      std::vector<std::tuple<unsigned int, unsigned, float> >::const_iterator it_atoms;
      for (unsigned int i=0; i<atom_atom_pair_indices.size(); i++) {
         const std::tuple<unsigned int, unsigned, float> &ap = atom_atom_pair_indices[i];
         atom_instance_matrices[i] = get_octahemi_matrix(std::get<0>(ap), std::get<1>(ap), std::get<2>(ap));
      }
   }

}

void
molecular_mesh_generator_t::move_the_atoms_and_update_the_instancing_matrices() {

   auto t_now = std::chrono::high_resolution_clock::now();
   auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t_now - t_start).count();
   auto dt_frame = std::chrono::duration_cast<std::chrono::milliseconds>(t_now - t_previous).count();
   // std::cout << "frame time " << dt_frame << std::endl; // use optimization in compilation
                                                           // for high frame rates.
   float tm = 1.6; // time multiplier
   float movement_scale = 2.3;
   for (unsigned int i=0; i<atom_positions.size(); i++) {
      float fi = static_cast<float>(i);
      float delta_x =  0.009 * movement_scale * std::sin(0.003 * tm * (dt + 500.0 * fi));
      float delta_y =  0.012 * movement_scale * std::sin(0.007 * tm * (dt + 600.0 * fi) + 1.7);
      float delta_z =  0.005 * movement_scale * std::sin(0.009 * tm * (dt + 440.0 * fi) + 1.2);
      atom_positions[i].x += delta_x;
      atom_positions[i].y += delta_y;
      atom_positions[i].z += delta_z;
   }

   update_mats_and_colours();
   t_previous = t_now; // for next time

}

std::vector<std::vector<std::pair<mmdb::Residue *, glm::vec3> > >
molecular_mesh_generator_t::make_CA_fragments(mmdb::Manager *mol) const {

   // Too simple minded - but good enough to generate points for testing

   std::vector<std::vector<std::pair<mmdb::Residue *, glm::vec3> > > fragments;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         std::vector<std::pair<mmdb::Residue *, glm::vec3> > frag;
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            int n_atoms = residue_p->GetNumberOfAtoms();
            for (int iat=0; iat<n_atoms; iat++) {
               mmdb::Atom *at = residue_p->GetAtom(iat);
               std::string atom_name(at->name);
               if (atom_name == " CA ") {
                  std::string alt_conf(at->altLoc);
                  if (alt_conf == "") {
                     std::pair<mmdb::Residue *, glm::vec3> p(at->residue, glm::vec3(at->x, at->y, at->z));
                     frag.push_back(p);
                     break;
                  }
               }
            }
         }
         if (frag.size() > 1) {

            // bring that fragment to the origin, I only want one for testing.
            glm::vec3 sum(0,0,0);
            for (unsigned int i=0; i<frag.size(); i++)
               sum += frag[i].second;
            glm::vec3 av = 1.0f/static_cast<float>(frag.size()) * sum;
            for (unsigned int i=0; i<frag.size(); i++)
               frag[i].second -= av;
            for (unsigned int i=0; i<frag.size(); i++)
               frag[i].second.x -= 6.9;

            fragments.push_back(frag);
            break;
         }
      }
   }
   return fragments;
}

#include "utils/dodec.hh"

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::dodec_at_position(float scale, const glm::vec3 &position) const {

   // These are dodec balls for the moment, not flat shaded real dodecs.
   // In a bit I will transform them to pentakis dodecs.

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

   dodec d;
   std::vector<clipper::Coord_orth> v = d.coords();
   vp.first.resize(v.size()); // should be 20

   for (unsigned int i=0; i<v.size(); i++) {
      glm::vec3 a = scale * glm::vec3(v[i].x(), v[i].y(), v[i].z());
      float f = static_cast<float>(i+1)/static_cast<float>(v.size());
      vp.first[i].pos = 0.2f * a;
      vp.first[i].pos += position;
      vp.first[i].normal = a;
      vp.first[i].color  = glm::vec4(0.3 + 0.4 * f, 0.8 - 0.8 * f, 0.1 + 0.9 * f, 1.0);
   }

   for (unsigned int i=0; i<12; i++) {
      const std::vector<unsigned int> &face = d.face(i);
      g_triangle gt_0(face[0], face[1], face[2]);
      g_triangle gt_1(face[0], face[2], face[3]);
      g_triangle gt_2(face[0], face[3], face[4]);
      vp.second.push_back(gt_0);
      vp.second.push_back(gt_1);
      vp.second.push_back(gt_2);
   }

   return vp;
}

// update vp.
//
// The day you add std::vector<s_generic_vertex> and std::vector<g_triangle> into this class
// is the day that you merge it with Mesh.
//
void
molecular_mesh_generator_t::add_to_mesh(std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > *vp,
                                        const std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > &vp_new) const {

   unsigned int idx_base = vp->first.size();
   unsigned int idx_tri_base = vp->second.size();

   vp->first.insert(vp->first.end(), vp_new.first.begin(), vp_new.first.end());
   vp->second.insert(vp->second.end(), vp_new.second.begin(), vp_new.second.end());;
   for (unsigned int i=idx_tri_base; i<vp->second.size(); i++)
      vp->second[i].rebase(idx_base);
}

void
molecular_mesh_generator_t::add_to_mesh(std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > *vp, // update vp
                                        const std::vector<s_generic_vertex> &gv,
                                        const std::vector<g_triangle> &tris) const {

   unsigned int idx_base = vp->first.size();
   unsigned int idx_tri_base = vp->second.size();

   vp->first.insert(vp->first.end(), gv.begin(), gv.end());
   vp->second.insert(vp->second.end(), tris.begin(), tris.end());
   for (unsigned int i=idx_tri_base; i<vp->second.size(); i++)
      vp->second[i].rebase(idx_base);

}

void
molecular_triangles_mesh_t::add_to_mesh(const std::vector<s_generic_vertex> &gv,
                                        const std::vector<g_triangle> &tris) {

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();

   vertices.insert(vertices.end(), gv.begin(), gv.end());
   triangles.insert(triangles.end(), tris.begin(), tris.end());
   for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);
}




std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::get_test_molecular_triangles_mesh(mmdb::Manager *mol, int secondary_structure_usage_flag) {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;

#ifdef USE_MOLECULES_TO_TRIANGLES

   int imodel = 1;
   mmdb::Model *model_p = mol->GetModel(imodel);
   if (model_p) {
      auto my_mol = std::make_shared<MyMolecule>(mol, secondary_structure_usage_flag);
      auto ss_cs = ColorScheme::colorBySecondaryScheme();
      auto ribbon_ramp_cs = ColorScheme::colorRampChainsScheme();
      auto chain_cs = ColorScheme::colorChainsScheme();
      auto this_cs = chain_cs;

      std::string atom_selection_str = "//";
      std::string style = "Ribbon";
      // style = "MolecularSurface";
      std::shared_ptr<MolecularRepresentationInstance> molrepinst_1 =
         MolecularRepresentationInstance::create(my_mol, this_cs, atom_selection_str, style);

      // now fill vp.first and vp.second
      // use fcxxcoord_to_glm() if needed.

      std::shared_ptr<Representation> r = molrepinst_1->getRepresentation();
      r->redraw();
      std::vector<std::shared_ptr<DisplayPrimitive> > vdp = r->getDisplayPrimitives();
      auto displayPrimitiveIter = vdp.begin();
      int i=0;
      std::vector<s_generic_vertex> vertices;
      std::vector<g_triangle> triangles;
      for (; displayPrimitiveIter != vdp.end(); displayPrimitiveIter++){
         DisplayPrimitive &displayPrimitive = **displayPrimitiveIter;
         // std::cout << "Display primitive type " << displayPrimitive.type() << "\n";
         if (displayPrimitive.type() == DisplayPrimitive::PrimitiveType::SurfacePrimitive ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::BoxSectionPrimitive ||
             displayPrimitive.type() == DisplayPrimitive::PrimitiveType::CylinderPrimitive ){
            displayPrimitive.generateArrays();

            VertexColorNormalPrimitive &surface = dynamic_cast<VertexColorNormalPrimitive &>(displayPrimitive);
            // std::cout << "nVertices " << surface.nVertices() << " nTriangles " << surface.nTriangles() << "\n";
            vertices.resize(surface.nVertices());

            auto vcnArray = surface.getVertexColorNormalArray();
            for (unsigned int iVertex=0; iVertex < surface.nVertices(); iVertex++){
               s_generic_vertex &gv = vertices[iVertex];
               VertexColorNormalPrimitive::VertexColorNormal &vcn = vcnArray[iVertex];
               for (int i=0; i<3; i++) {
                  gv.pos[i] = vcn.vertex[i];
                  gv.normal[i] = vcn.normal[i];
                  gv.color[i]  = 0.0037f * vcn.color[i];
               }
               gv.color[3] = 1.0;
            }

            auto indexArray = surface.getIndexArray();
            unsigned long nIndices = 3 * surface.nTriangles();
            triangles.resize(surface.nTriangles());
            for (unsigned int iTriangle=0; iTriangle<surface.nTriangles(); iTriangle++){
               g_triangle &gt = triangles[iTriangle];
               for (int i=0; i<3; i++)
                  gt[i] = indexArray[3*iTriangle+i];
            }
            add_to_mesh(&vp, vertices, triangles);
         }
      }
   }

#endif // USE_MOLECULES_TO_TRIANGLES

   return vp;
}



// Maybe this should just be part of Mesh (20211004-PE but that doesn't as yet use mmdb)
//
// maybe this should be const?
//

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecular_mesh_generator_t::get_molecular_triangles_mesh_for_active_residue(int imol, mmdb::Manager *mol,
                                                                            int model_number,
                                                                            mmdb::Residue *residue_p,
                                                                            const coot::protein_geometry *geom_in,
                                                                            int bond_width,
                                                                            molecular_mesh_generator_t::range_t range,
                                                                            int selection_mode) {

   //! selection_mode is one of 1: residue, 2: sphere, 3: big sphere, 4: chain

   auto select_atoms_in_residues = [] (int selhnd, mmdb::Manager *mol,
                                       const std::vector<mmdb::Residue *> &rs) {
      // is there a faster way than this?
      std::vector<mmdb::Residue * >::const_iterator it;
      for (it=rs.begin(); it!=rs.end(); ++it) {
         mmdb::Residue *r(*it);
         mol->SelectAtoms(selhnd, 0, r->GetChainID(), r->GetSeqNum(), "*", r->GetSeqNum(), "*",
                          "*", "*", "*", "*", mmdb::SKEY_OR);
      }
   };

   auto residues_in_range = [] (mmdb::Chain *chain_p, int resno_start, int resno_end) {
      std::vector<mmdb::Residue * > v;
      if (chain_p) {
         int n_residues = chain_p->GetNumberOfResidues();
         for (int i=0; i<n_residues; i++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(i);
            if (residue_p) {
               int res_no = residue_p->GetSeqNum();
               if (res_no >= resno_start) {
                  if (res_no <= resno_end) {
                     v.push_back(residue_p);
                  }
               }
            }
         }
      }
      return v;
   };

   auto construct_selection = [select_atoms_in_residues,
                               residues_in_range] (mmdb::Manager *mol,
                                                   range_t range,
                                                   int selection_mode,
                                                   mmdb::Residue *residue_p) {
      // c.f. coot::molecule_t::select_residues()
      int selhnd = mol->NewSelection();
      if (selection_mode == 1) { // residue
         coot::residue_spec_t res_spec(residue_p);
         res_spec.select_atoms(mol, selhnd, mmdb::SKEY_NEW);
      }
      if (selection_mode == 2) { // sphere
         float radius = 4.2;
         std::vector<mmdb::Residue * > rs = coot::residues_near_residue(residue_p, mol, radius);
         rs.push_back(residue_p);
         select_atoms_in_residues(selhnd, mol, rs);
      }
      if (selection_mode == 3) { // big sphere
         float radius = 8.0; // matches select_residues() in api
         std::vector<mmdb::Residue * > rs = coot::residues_near_residue(residue_p, mol, radius);
         rs.push_back(residue_p);
         select_atoms_in_residues(selhnd, mol, rs);
      }
      if (selection_mode == 4) { // chain
         mmdb::Model *model_p = mol->GetModel(1);
         if (model_p) {
            mmdb::Chain *chain_p = residue_p->chain;
            std::vector<mmdb::Residue * > rs = coot::util::residues_in_chain(chain_p);
            select_atoms_in_residues(selhnd, mol, rs);
         }
      }
      if (selection_mode == 5) { // residue range
         // if (graphics_info_t::in_range_define == 2) {
         if (true) {
            mmdb::Model *model_p = mol->GetModel(1);
            if (model_p) {
               mmdb::Chain *chain_p = residue_p->chain;
               if (range.is_valid) {
                  std::vector<mmdb::Residue * > rs = residues_in_range(chain_p, range.resno_start, range.resno_end);
                  if (! rs.empty()) {
                     select_atoms_in_residues(selhnd, mol, rs);
                  }
               }
            }
         }
      }
      return selhnd;
   };

   float atom_radius_scale_factor = 1.0; // this needs to be passed

   coot::residue_spec_t res_spec(residue_p); // pass the spec!
   // int selhnd = mol->NewSelection(); // d
   // res_spec.select_atoms(mol, selhnd, mmdb::SKEY_NEW);
   // replace asc innards with the atoms of the residue selection
   atom_selection_container_t asc = make_asc(mol);
   int selhnd = construct_selection(mol, range, selection_mode, residue_p); // selection is deleted below
   int n_selected_atoms = 0;
   mmdb::PAtom *atom_selection;
   mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
   asc.SelectionHandle = selhnd;
   asc.atom_selection = atom_selection;
   asc.n_selected_atoms = n_selected_atoms;
   std::set<int> no_bonds;
   int include_disulphides = false;
   int include_hydrogens = true;
   bool draw_missing_loops_flag = false;

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
      return glm::vec3(c.x(), c.y(), c.z());
   };

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

   if (n_selected_atoms == 0) {
      std::cout << "Play Nothing feedback here" << std::endl;
   }

   auto make_generic_vertices_for_atoms = [atom_radius_scale_factor, cartesian_to_glm, vnc_vertex_to_generic_vertex]
      (const graphical_bonds_container &bonds_box,
       const std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > &octaball) {
                                             // change the vertex type later
                                             std::vector<s_generic_vertex> vertices;
                                             std::vector<g_triangle> triangles;
                                             glm::mat4 unit_matrix(1.0f);

                                             float sphere_radius = 0.085 * 1.54; // how big should atoms be?
                                             float bond_width = 4.0; // pass this
                                             float radius_scale = 0.2 * bond_width; // arbs
                                             radius_scale *= atom_radius_scale_factor;
                                             glm::vec3 sphere_scaling(sphere_radius, sphere_radius, sphere_radius);

                                             for (int icol=0; icol<bonds_box.n_consolidated_atom_centres; icol++) {
                                                for (unsigned int i=0; i<bonds_box.consolidated_atom_centres[icol].num_points; i++) {
                                                   const graphical_bonds_atom_info_t &ai = bonds_box.consolidated_atom_centres[icol].points[i];
                                                   float sphere_scale = radius_scale * ai.radius_scale * 1.18;
                                                   if (ai.is_hydrogen_atom)
                                                      sphere_scale *= 0.6;
                                                   glm::vec3 atom_position = cartesian_to_glm(ai.position);
                                                   unsigned int idx_base = vertices.size();

                                                   for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
                                                      s_generic_vertex vertex(vnc_vertex_to_generic_vertex(octaball.first[ibv]));
                                                      vertex.pos = sphere_scaling * vertex.pos + atom_position;
                                                      vertices.push_back(vertex);
                                                   }

                                                   std::vector<g_triangle> octaball_triangles = octaball.second;
                                                   for (unsigned int ii=0; ii<octaball_triangles.size(); ii++)
                                                      octaball_triangles[ii].rebase(idx_base);
                                                   triangles.insert(triangles.end(), octaball_triangles.begin(), octaball_triangles.end());
                                                }
                                             }
                                             return std::make_pair(vertices, triangles);
                                          };


   auto make_generic_vertices_for_bonds = [atom_radius_scale_factor, cartesian_to_glm,
                                           vnc_vertex_vector_to_generic_vertex_vector]
      (const graphical_bonds_container &bonds_box,
       float bond_width) {
                                             std::vector<s_generic_vertex> vertices;
                                             std::vector<g_triangle> triangles;

                                             unsigned int n_slices = 16;
                                             unsigned int n_stacks = 2;
                                             glm::vec4 col(0.5, 0.5, 0.5, 1.0);

                                             float radius = 0.02f * bond_width;

                                             for (int i=0; i<bonds_box.num_colours; i++) {
                                                graphical_bonds_lines_list<graphics_line_t> &ll = bonds_box.bonds_[i];

                                                bool do_thinning = false;
                                                if (ll.thin_lines_flag) do_thinning = true;

                                                for (int j=0; j< ll.num_lines; j++) {

                                                   glm::vec3 start(cartesian_to_glm(ll.pair_list[j].positions.getStart()));
                                                   glm::vec3 finish(cartesian_to_glm(ll.pair_list[j].positions.getFinish()));
                                                   std::pair<glm::vec3, glm::vec3> pos_pair(start, finish); // correct way round?
                                                   glm::vec3 b = finish - start;
                                                   float bl = glm::distance(b, glm::vec3(0,0,0));
                                                   float radius_scale = 1.0;
                                                   if (do_thinning) radius_scale *= 0.5;
                                                   cylinder c(pos_pair, radius * radius_scale, radius * radius_scale, bl, col, n_slices, n_stacks);

                                                   unsigned int idx_base = vertices.size();
                                                   unsigned int idx_tri_base = triangles.size();
                                                   std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(c.vertices);

                                                   vertices.insert(vertices.end(), converted_vertices.begin(), converted_vertices.end());
                                                   triangles.insert(triangles.end(), c.triangles.begin(), c.triangles.end());
                                                   for (unsigned int k=idx_tri_base; k<triangles.size(); k++)
                                                      triangles[k].rebase(idx_base);
                                                }
                                             }
                                             return std::make_pair(vertices, triangles);
                                          };

   float radius = 1.0f;
   glm::vec4 col(0.5, 0.5, 0.5, 1.0);
   glm::vec3 origin(0,0,0);
   unsigned int num_subdivisions = 2; // 2 should be the default?
   std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > octaball =
      make_octasphere(num_subdivisions, origin, radius, col);

   Bond_lines_container bonds(asc, imol, no_bonds, geom_in, include_disulphides, include_hydrogens, draw_missing_loops_flag,
                              model_number, "");

   auto bonds_box = bonds.make_graphical_bonds();

   mol->DeleteSelection(selhnd);

   std::vector<glm::vec4> atom_colour_vector(bonds_box.num_colours); // colour is not used to display the outline, but *is* used
                                                                     // to generate the vertices

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > atom_bits = make_generic_vertices_for_atoms(bonds_box, octaball);

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > bond_bits = make_generic_vertices_for_bonds(bonds_box, bond_width);

   const auto &atom_vertices = atom_bits.first;
   std::vector<s_generic_vertex> vertices(atom_vertices.size());
   for (unsigned int i=0; i<atom_vertices.size(); i++)
      vertices[i] = atom_vertices[i];
   std::vector<g_triangle> triangles = atom_bits.second;

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), bond_bits.first.begin(), bond_bits.first.end());
   triangles.insert(triangles.end(), bond_bits.second.begin(), bond_bits.second.end());
   for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);

   return std::make_pair(vertices, triangles);
}

