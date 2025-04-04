
#include <iostream>

#include <clipper/core/coords.h>
#include "MoleculesToTriangles/CXXClasses/NRStuff.h"
#include "coot-utils/cylinder.hh"

#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"
#include "geometry/residue-and-atom-specs.hh"
#include "mmdb2/mmdb_selmngr.h"
#include "tubes.hh"




void secondary_structure_header_to_residue_sse(mmdb::Manager *mol);

coot::simple_mesh_t
make_mesh_for_helical_representation(mmdb::Manager *mol,
                                     const std::string &atom_selection_str,
                                     float radius_for_helices,
                                     unsigned int n_slices_for_helices) {
   coot::simple_mesh_t m;
   class helix_t{
      public:
      std::string chain_id;
      int resno_start;
      int resno_end;
      helix_t() : chain_id(""), resno_start(-9999), resno_end(-9999) {}
      helix_t(const std::string &c, int r1, int r2) : chain_id(c), resno_start(r1), resno_end(r2) {}
   };

   auto get_ref_coords = [] (mmdb::Manager *helix_mol) {

      std::vector<clipper::Coord_orth> coords;
      int imod = 1;
      mmdb::Model *model_p = helix_mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  int res_no = residue_p->GetSeqNum();
                  if (res_no > 3) break;
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     clipper::Coord_orth pt(at->x, at->y, at->z);
                     coords.push_back(pt);
                  }
               }
            }
         }
      }
      return coords;
   };

   auto get_rtop = [] (const std::string &chain_id, int res_no, mmdb::Manager *mol, std::vector<clipper::Coord_orth> &ref_coords) {
      clipper::RTop_orth rtop;

      std::vector<clipper::Coord_orth> coords;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id_this = chain_p->GetChainID();
            if (chain_id_this ==  chain_id) {
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int seq_num_this = residue_p->GetSeqNum();
                     if (seq_num_this >= res_no) {
                        if (seq_num_this <= (res_no+2)) {
                           if (seq_num_this > (res_no+2)) break;
                           int n_atoms = residue_p->GetNumberOfAtoms();
                           for (int iat=0; iat<n_atoms; iat++) {
                              mmdb::Atom *at = residue_p->GetAtom(iat);
                              if (! at->isTer()) {
                                 std::string atom_name = at->name;
                                 // PDBv3 FIXME
                                 if (atom_name == " N  " ||
                                     atom_name == " CA " ||
                                     atom_name == " C  " ||
                                     atom_name == " O  ")  {
                                    mmdb::Atom *at = residue_p->GetAtom(iat);
                                    clipper::Coord_orth pt(at->x, at->y, at->z);
                                    coords.push_back(pt);
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
      if (ref_coords.size() == coords.size()) {
         rtop = clipper::RTop_orth(ref_coords, coords); // maybe the other way around
      }
      return rtop;
   };

   // find all the helices:
   std::vector<helix_t> helices;

   int sel_hnd = mol->NewSelection(); // d
   mol->Select(sel_hnd, mmdb::STYPE_RESIDUE, atom_selection_str.c_str(), mmdb::SKEY_NEW);

   for (int imod = 1; imod <= mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p-> GetNumberOfResidues();
            helix_t helix_running;
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  if (residue_p->SSE == mmdb::SSE_Helix) {
                     if (residue_p->isInSelection(sel_hnd)) {
                        if (helix_running.resno_start == -9999) {
                           helix_running.resno_start = residue_p->GetSeqNum();
                           helix_running.chain_id = chain_p->GetChainID();
                        } else {
                           // we are in a helix
                           helix_running.resno_end = residue_p->GetSeqNum();
                        }
                     }
                  } else {
                     if (helix_running.resno_start == -9999) {
                        // nothing
                     } else {
                        helices.push_back(helix_running);
                        helix_running = helix_t();
                     }
                  }
               }
            }
            if (helix_running.resno_start != -9999) {
               helices.push_back(helix_running);
            }
         }
      }
   }
   mol->DeleteSelection(sel_hnd);

   if (true) { // debug
      for (const auto &h : helices) {
         std::cout << "helix: " << h.chain_id << " " << h.resno_start << " " << h.resno_end
                   << std::endl;
      }
   }

   if (! helices.empty()) {
      mmdb::Manager *helix_mol = new mmdb::Manager();
      helix_mol->ReadCoorFile("theor-helix-z-ori-v2.pdb");
      std::vector<clipper::Coord_orth> ref_coords = get_ref_coords(helix_mol);

      // let's pre-make the ring points
      std::vector<clipper::Coord_orth> ring_points;
      std::vector<clipper::Coord_orth> ring_point_normals;
      for (unsigned int i_pt=0; i_pt<n_slices_for_helices; i_pt++) {
         float theta = M_PI * 2.0 * static_cast<float>(i_pt) / static_cast<float>(n_slices_for_helices);
         float x = radius_for_helices * sinf(theta);
         float y = radius_for_helices * cosf(theta);
         float z = 0;
         float xn = sinf(theta);
         float yn = cosf(theta);
         float zn = 0;
         ring_points.push_back(clipper::Coord_orth(x,y,z));
         ring_point_normals.push_back(clipper::Coord_orth(xn,yn,zn));
      }
      glm::vec4 col(0.5, 0.5, 0.5, 1.0);
      for (const auto &helix : helices) {

         std::vector<coot::api::vnc_vertex> vertices;
         std::vector<g_triangle> triangles;

         unsigned int n_rings_in_helix = 0;
         for (int ires_no=helix.resno_start; ires_no<(helix.resno_end-1); ires_no++) {
            clipper::RTop_orth rtop = get_rtop(helix.chain_id, ires_no, mol, ref_coords);
            if (true) { // test here that rtop is sane
               clipper::RTop_orth rtop_for_normals = rtop;
               rtop_for_normals.trn() = clipper::Coord_orth(0,0,0);
               n_rings_in_helix++;
               std::vector<clipper::Coord_orth> transformed_points(ring_points.size());
               std::vector<clipper::Coord_orth>    rotated_normals(ring_points.size());
               for (unsigned int i_pt=0; i_pt<ring_points.size(); i_pt++) {
                  const clipper::Coord_orth &pt   = ring_points[i_pt];
                  const clipper::Coord_orth &n_pt = ring_point_normals[i_pt];
                  clipper::Coord_orth t_pt = rtop * pt;
                  clipper::Coord_orth t_n = rtop_for_normals * n_pt;
                  transformed_points[i_pt] = t_pt;
                  rotated_normals[i_pt] = t_n; 
                  glm::vec3 v(t_pt.x(), t_pt.y(), t_pt.z());
                  glm::vec3 n(t_n.x(), t_n.y(), t_n.z());
                  coot::api::vnc_vertex vnc_v(v, n, col);
                  vertices.push_back(vnc_v);

               }
            }
         }

         // now make the triangle for that helix
         //
         // try simple-minded first
         for (unsigned int i=0; i<(n_rings_in_helix-1); i++) {
            unsigned int ring_offset = i * n_slices_for_helices;
            for (unsigned int j=0; j<n_slices_for_helices-1; j++) {
               // winding is important
               g_triangle t1(ring_offset + j + 1, ring_offset + j, ring_offset + n_slices_for_helices + j);
               g_triangle t2(ring_offset + n_slices_for_helices + j, ring_offset + n_slices_for_helices + j + 1, ring_offset + j + 1);
               triangles.push_back(t1);
               // 
               if (t1.point_id[0]>vertices.size() ||
                   t1.point_id[1]>vertices.size() ||
                   t1.point_id[2]>vertices.size()) {
                  std::cout << "ERROR:: triangle point id out of range " << t1.point_id[0] << " "
                            << t1.point_id[1] << " " << t1.point_id[2] << std::endl;
               }
            }
            // now the join the end to the start:
            g_triangle t_end_1(ring_offset, ring_offset + n_slices_for_helices - 1, ring_offset + n_slices_for_helices);
            g_triangle t_end_2(ring_offset + n_slices_for_helices, ring_offset + n_slices_for_helices -1, ring_offset + 2 * n_slices_for_helices -1);
            triangles.push_back(t_end_1);
            triangles.push_back(t_end_2);
         }
         m.add_submesh(coot::simple_mesh_t(vertices, triangles));
      }
   }

   return m;
}

coot::simple_mesh_t
make_tubes_representation(mmdb::Manager *mol,
                          const std::string &atom_selection_str,
                          const std::string &colour_scheme,
                          float radius_for_coil,
                          int Cn_for_coil, int accuracy_for_coil,
                          unsigned int n_slices_for_coil,
                          int secondaryStructureUsageFlag) {

   coot::simple_mesh_t m;
   bool remove_trace_for_helices = true;
   secondary_structure_header_to_residue_sse(mol);
   float radius_for_helices = 3.2;
   unsigned int n_slices_for_helices = 16;

   coot::simple_mesh_t mesh_for_helices = make_mesh_for_helical_representation(
      mol, atom_selection_str, radius_for_helices, n_slices_for_helices);

   m.add_submesh(mesh_for_helices);
   m.add_submesh(make_coil_for_tubes_representation(mol, atom_selection_str,
                                                    radius_for_coil,
                                                    Cn_for_coil, accuracy_for_coil,
                                                    n_slices_for_coil,
                                                    remove_trace_for_helices));
   return m;

}

coot::simple_mesh_t
make_coil_for_tubes_representation(mmdb::Manager *mol,
                                   const std::string &atom_selection_str,
                                   float radius,
                                   int Cn, int accuracy,
                                   unsigned int n_slices,
                                   bool remove_trace_for_helices) {

   int iinterp = 1;

   coot::simple_mesh_t m_all;
   if (n_slices < 3) return m_all;

   auto fcxx_to_glm = [] (const FCXXCoord &c) {
      return glm::vec3(c.x(), c.y(), c.z());
   };

   auto make_cylinders = [fcxx_to_glm] (const std::vector<FCXXCoord> &spline_points) {
      std::vector<cylinder> cylinders;
      float r = 0.3f;
      unsigned int n_slices = 10;
      unsigned int n_stacks = 2;
      if (spline_points.size() > 1) {
         std::size_t spline_points_end = spline_points.size() -1;
         for (std::size_t i=0; i<spline_points_end; i++) {
            const auto &pt_bottom = spline_points[i];
            const auto &pt_top    = spline_points[i+1];
            float base_radius = r;
            float  top_radius = r;
            glm::vec3 b = fcxx_to_glm(pt_bottom);
            glm::vec3 t = fcxx_to_glm(pt_top);
            float height = glm::distance(b,t);
            cylinder c(std::make_pair(b, t), base_radius, top_radius,
                       height, n_slices, n_stacks);
            if (i == 0) c.add_flat_start_cap();
            cylinders.push_back(c);
         }
      }
      return cylinders;
   };


   // I don't know how to do the matrix multiplication needed at the heart
   // of rotate_around_vector() using FCXXCoord (and matrices?)
   //
   // So let's use Coot functions.
   // A bit inconsistent - I know.
   //
   // I copy rotate_around_vector() here because I don't want to link
   // coot libraries from here. Way too much trouble for just that function.
   // Maybe I could have used glm::rotate() or some such.

   // angle in radians.

   auto rotate_around_vector = [] (const glm::vec3 &direction,
                                   const glm::vec3 &position,
                                   const glm::vec3 &origin_shift,
                                   double angle) {

      glm::vec3 unit_vec = glm::normalize(direction);

      double l = unit_vec[0];
      double m = unit_vec[1];
      double n = unit_vec[2];

      double ll = l*l;
      double mm = m*m;
      double nn = n*n;
      double cosk = cos(angle);
      double sink = sin(angle);
      double I_cosk = 1.0 - cosk;

      glm::mat3 r(ll+(mm+nn)*cosk,    l*m*I_cosk-n*sink,  n*l*I_cosk+m*sink,
                  l*m*I_cosk+n*sink,  mm+(ll+nn)*cosk,    m*n*I_cosk-l*sink,
                  n*l*I_cosk-m*sink,  m*n*I_cosk+l*sink,  nn+(ll+mm)*cosk );

      glm::vec3 p1 = position - origin_shift;
      glm::vec3 p2 = r * p1;
      glm::vec3 p3 = p2 + origin_shift;
      return p3;
   };

   auto get_ring_around_mid_point = [rotate_around_vector] (const glm::vec3 &p1,
                                                            const glm::vec3 &p2,
                                                            const glm::vec3 &p3,
                                                            unsigned int n_slices,
                                                            float radius) {

      std::vector<std::pair<glm::vec3, glm::vec3> > v;
      glm::vec3 v1 = glm::normalize(p2-p1);
      glm::vec3 v2 = glm::normalize(p3-p2);
      glm::vec3 delta = (v2 - v1) * 0.5f;
      glm::vec3 delta_uv = glm::normalize(delta);
      glm::vec3 pos_start = p2 + radius *  delta_uv;
      glm::vec3 rotation_ori = v1 + v2;
      glm::vec3 origin(0,0,0);
      // this is the vector to spin around the point on the circle
      for (unsigned int i_slice=0; i_slice<n_slices; i_slice++) {
         double angle = 2 * M_PI * static_cast<double>(i_slice) / static_cast<double>(n_slices);
         glm::vec3 pt = rotate_around_vector(rotation_ori, pos_start, p2, angle);
         glm::vec3 n  = rotate_around_vector(rotation_ori, delta_uv,  origin, angle);
         v.push_back(std::make_pair(pt, n));
      }
      return v;
   };

   auto make_continuous_spline_mesh = [get_ring_around_mid_point, fcxx_to_glm]
      (const std::vector<FCXXCoord> &spline_points,
       float radius,
       unsigned int n_slices) {

      std::vector<std::vector<std::pair<glm::vec3, glm::vec3> > > rings;
      coot::simple_mesh_t m;
      if (spline_points.size() > 2) {
         std::size_t spline_points_end = spline_points.size() -1;
         for (std::size_t i=1; i<spline_points_end; i++) {
            glm::vec3 start_point = fcxx_to_glm(spline_points[i-1]);
            glm::vec3 mid_point   = fcxx_to_glm(spline_points[i]);
            glm::vec3 next_point  = fcxx_to_glm(spline_points[i+1]);
            std::vector<std::pair<glm::vec3, glm::vec3> > r =
               get_ring_around_mid_point(start_point, mid_point, next_point, n_slices, radius);
            rings.push_back(r);
         }

         glm::vec4 col(0.5, 0.5, 0.5, 1.0);
         for (const auto &ring : rings) {
            for (const auto &point : ring) {
               coot::api::vnc_vertex v(point.first, point.second, col);
               m.vertices.push_back(v);
            }
         }

         // simple minded first:
         bool simple_minded = false;
         if (simple_minded) {

            // so now we have the vertices - let's make some triangles!
            for (unsigned int i=0; i<(rings.size()-1); i++) {
               unsigned int ring_offset = i * n_slices;
               for (unsigned int j=0; j<n_slices-1; j++) {
                  // winding is important
                  g_triangle t1(ring_offset + j + 1, ring_offset + j, ring_offset + n_slices + j);
                  g_triangle t2(ring_offset + n_slices + j, ring_offset + n_slices + j + 1, ring_offset + j + 1);
                  m.triangles.push_back(t1);
                  m.triangles.push_back(t2);
               }
               // now the join the end to the start:
               g_triangle t_end_1(ring_offset, ring_offset + n_slices - 1, ring_offset + n_slices);
               g_triangle t_end_2(ring_offset + n_slices, ring_offset + n_slices -1, ring_offset + 2 * n_slices -1);
               m.triangles.push_back(t_end_1);
               m.triangles.push_back(t_end_2);
            }

         } else {

            // "simple minded" mostly works - but there are places where the
            // "next" ring doesn't sit above the current one - there is a rotation.
            for (unsigned int i=0; i<(rings.size()-1); i++) {
               unsigned int ring_offset = i * n_slices;
               const glm::vec3 &p0 = rings[i][0].first;
               float dd_best = 9999999.9;
               unsigned int j_best = 13;
               for (unsigned int j=0; j<n_slices; j++) {
                  const glm::vec3 &testing_pos = rings[i+1][j].first;
                  glm::vec3 d = testing_pos - p0;
                  float dd_test = d.x * d.x + d.y * d.y + d.z * d.z;
                  if (dd_test < dd_best) {
                     dd_best = dd_test;
                     j_best = j;
                  }
               }
               if (j_best != 13) {
                  for (unsigned int j=0; j<(n_slices-1); j++) {
                     unsigned jo_1 = j_best + j;
                     unsigned jo_2 = j_best + j + 1;
                     if (jo_1 >= n_slices) jo_1 -= n_slices;
                     if (jo_2 >= n_slices) jo_2 -= n_slices;
                     // winding is important
                     g_triangle t1(j + 1, j, n_slices + jo_1);
                     g_triangle t2(n_slices + jo_1, n_slices + jo_2, j + 1);
                     t1.rebase(ring_offset);
                     t2.rebase(ring_offset);
                     m.triangles.push_back(t1);
                     m.triangles.push_back(t2);
                  }
                  // now the join the end to the start: - I had to get a pen and paper out for this...
                  unsigned int jbo = j_best + n_slices - 1;
                  if (j_best == 0) jbo = 2 * n_slices - 1;
                  g_triangle t_end_1(0, n_slices - 1, j_best + n_slices);
                  g_triangle t_end_2(j_best + n_slices, n_slices -1, jbo);
                  t_end_1.rebase(ring_offset);
                  t_end_2.rebase(ring_offset);
                  m.triangles.push_back(t_end_1);
                  m.triangles.push_back(t_end_2);
               }
            }
         }
      }
      return m;
   };

   // ------------------ main line ------------------------

   std::vector<std::vector<mmdb::Residue *> > runs_of_residues;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {

      int sel_hnd = mol->NewSelection(); // d
      mol->Select(sel_hnd, mmdb::STYPE_RESIDUE, atom_selection_str.c_str(), mmdb::SKEY_NEW);
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         std::vector<mmdb::Residue *> a_run_of_residues;
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p->isInSelection(sel_hnd)) {
               bool keep_residue = true;
               if (remove_trace_for_helices)
                  if (residue_p->SSE == mmdb::SSE_Helix)
                     keep_residue = false;
               if (keep_residue) {
                  if (a_run_of_residues.empty()) {
                     a_run_of_residues.push_back(residue_p);
                  } else {
                     int idx_r_1 = a_run_of_residues.back()->index;
                     int idx_r_2 = residue_p->index;
                     if (idx_r_2 == (idx_r_1 + 1)) {
                        a_run_of_residues.push_back(residue_p);
                     } else {
                        runs_of_residues.push_back(a_run_of_residues);
                        a_run_of_residues.clear();
                     }
                  }
               } else {
                  if (false)
                     std::cout << "cut out residue " << coot::residue_spec_t(residue_p) << std::endl;
                  if (! a_run_of_residues.empty()) {
                     runs_of_residues.push_back(a_run_of_residues);
                     a_run_of_residues.clear();
                  }
               }
            } else {
               if (! a_run_of_residues.empty()) {
                  runs_of_residues.push_back(a_run_of_residues);
                  a_run_of_residues.clear();
               }
            }
         }
         if (! a_run_of_residues.empty())
            runs_of_residues.push_back(a_run_of_residues);
      }
      mol->DeleteSelection(sel_hnd);
   }
   if (! runs_of_residues.empty()) {

      // std::cout << "::::::::::::::::::::::::: found " << runs_of_residues.size() << " runs of residues"
      // << std::endl;

      for (unsigned int i=0; i<runs_of_residues.size(); i++) {
         std::vector<mmdb::Residue *> &a_run_of_residues = runs_of_residues[i];
         if (a_run_of_residues.size() < 2) continue;
         std::vector<FCXXCoord> ctlPts; // From Martin
         for (unsigned int ir=0; ir<a_run_of_residues.size(); ir++) {
            mmdb::Residue *residue_p = a_run_of_residues[ir];
            mmdb::Atom *ca_at = residue_p->GetAtom(" CA ");
            if (ca_at) {
               FCXXCoord fc(ca_at->x, ca_at->y, ca_at->z);
               ctlPts.push_back(fc);
            }
         }
         CoordSpline cs;
         int nsteps =  accuracy * (ctlPts.size() - 1);

         try {

            if (true) {
               std::cout << "debug:: a run of residues " << i << " has "
                         << ctlPts.size() << " control points" << std::endl;
               for (unsigned int ii=0; ii<ctlPts.size(); ii++)
                  std::cout << "debug:: control point     " << ctlPts[ii] << std::endl;
            }

            std::vector<FCXXCoord> v = cs.SplineCurve(ctlPts, nsteps, Cn, iinterp);
            // std::cout << "spline-curve size " << v.size() << std::endl;

            if (true) {
               for (unsigned int ii=0; ii<v.size(); ii++)
                  std::cout << "debug:: spline-curve-point     " << v[ii] << std::endl;
            }

            bool simple_mode = false;
            if (simple_mode) {
               std::vector<cylinder> cylinders = make_cylinders(v);
               for (unsigned int i=0; i<cylinders.size(); i++) {
                  const auto &c = cylinders[i];
                  m_all.add_submesh(coot::simple_mesh_t(c.vertices, c.triangles));
               }
            } else {
               coot::simple_mesh_t cs = make_continuous_spline_mesh(v, radius, n_slices);
               m_all.add_submesh(cs);
            }
         }
         catch (const std::length_error &le) {
            std::cout << "ERROR:: length error " << le.what() << std::endl;
         }
      }
   }
   return m_all;

}
