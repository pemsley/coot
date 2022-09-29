
#include <iostream>
#include <vector>

#include <glm/ext.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>

#include "g_triangle.hh"
#include "vertex.hh"

#include "prideout-octasphere.hh"
#include "oct.hh"

std::vector<g_triangle>
make_octasphere_triangles(unsigned int i_row,
                          unsigned int geodesic_verts_size,
                          unsigned int verts_size) {

   std::vector<g_triangle> triangles;

   if (i_row != 0) {
      if (geodesic_verts_size > 1) {
         unsigned int idx_base = verts_size;
         unsigned int v_index_base = idx_base - geodesic_verts_size - 1;
         for (unsigned int j=0; j<geodesic_verts_size-1; j++) {
            unsigned int idx_lower = v_index_base + j + 1;
            g_triangle t(idx_base+j+1, idx_base+j, idx_lower);
            triangles.push_back(t);
         }

         for (unsigned int j=0; j<geodesic_verts_size-1; j++) {
            // upward triangles
            unsigned int idx_lower = v_index_base + j + 1;
            g_triangle t(idx_lower, idx_lower+1, idx_base+j+1);
            triangles.push_back(t);
         }
         // left hand edge:
         unsigned int start_prev = v_index_base;
         g_triangle t(start_prev, start_prev + 1, idx_base);
         triangles.push_back(t);
      } else {
         // top geodesic_vert
         unsigned int idx_base = verts_size;
         g_triangle t(idx_base-2, idx_base-1, idx_base);
         triangles.push_back(t);
      }
   }
   return triangles;
}


std::pair<std::vector<glm::vec3>, std::vector<g_triangle> >
tessellate_hemisphere_patch(unsigned int num_subdivisions) {

   std::vector<glm::vec3> verts;
   std::vector<g_triangle> triangles;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > t =
      tessellate_octasphere_patch(num_subdivisions);

   // rotate t 4 times to make a hemisphere

   verts = t.first;
   triangles = t.second;

   for (unsigned int irot=1; irot<4; irot++) {
      float angle = 0.5f * M_PI * static_cast<float>(irot);
      unsigned int idx_base = verts.size();
      unsigned int idx_triangle_base = triangles.size();
      for (unsigned int i=0; i<t.first.size(); i++) {
         glm::vec3 v = glm::rotate(t.first[i], angle, glm::vec3(0,0,1));
         verts.push_back(v);
      }
      // reindex the triangles
      triangles.insert(triangles.end(), t.second.begin(), t.second.end());
      for (unsigned int i=idx_triangle_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);
   }

   return std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > (verts,triangles);
}

// adjust triangles_p - the vertices are not changed - so if this is successful,
// then vertices that are not used will be sent to the graphics card.
// Not many though, in the scheme of things. Also, fixing the vertex indexing might be painful.
void
remove_redundant_vertices(std::vector<g_triangle> *triangles_p,
                          const std::map<unsigned int, std::set<unsigned int> > &redundant_map) {

   std::map<unsigned int, std::set<unsigned int> >::const_iterator it;
   unsigned int n_done = 0;
   for (it=redundant_map.begin(); it!=redundant_map.end(); ++it) {
      const std::set<unsigned int> &s = it->second;
      std::set<unsigned int>::const_iterator it_s;

      if (false) {
         std::cout << "reference idx " << it->first << " has sames ";
         for (it_s=s.begin(); it_s!=s.end(); ++it_s)
            std::cout << " " << *it_s;
         std::cout << std::endl;
      }

      for (it_s=s.begin(); it_s!=s.end(); ++it_s) {
         for (unsigned int i=0; i<triangles_p->size(); i++) {
            g_triangle &tri = triangles_p->at(i);

            if (false) {
               if (tri.point_id[0] == *it_s) n_done++;
               if (tri.point_id[1] == *it_s) n_done++;
               if (tri.point_id[2] == *it_s) n_done++;
            }

            if (tri.point_id[0] == *it_s) tri.point_id[0] = it->first;
            if (tri.point_id[1] == *it_s) tri.point_id[1] = it->first;
            if (tri.point_id[2] == *it_s) tri.point_id[2] = it->first;

         }
      }
   }
   // std::cout << "remove_redundant_vertices() n_done " << n_done << std::endl;
}


std::pair<std::vector<glm::vec3>, std::vector<g_triangle> >
tessellate_octasphere(unsigned int num_subdivisions, bool remove_redundant_vertices_flag) {

   std::vector<glm::vec3> verts;
   std::vector<g_triangle> triangles;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > t =
      tessellate_octasphere_patch(num_subdivisions);

   // rotate t 4 times to make a hemisphere
   // rotate that 2 times to make a sphere

   glm::vec3 x_axis(1.0f, 0.0f, 0.0f);
   glm::vec3 z_axis(0.0f, 0.0f, 1.0f);
   for (unsigned int ih=0; ih<2; ih++) {
      for (unsigned int irot=0; irot<4; irot++) {
         float angle = 0.5f * M_PI * static_cast<float>(irot);
         unsigned int idx_base = verts.size();
         unsigned int idx_triangle_base = triangles.size();
         for (unsigned int i=0; i<t.first.size(); i++) {
            glm::vec3 v = glm::rotate(t.first[i], angle, z_axis);
            if (ih==1)
               v = glm::rotate(v, static_cast<float>(M_PI), x_axis);
            verts.push_back(v);
         }
         // reindex the triangles
         triangles.insert(triangles.end(), t.second.begin(), t.second.end());
         for (unsigned int i=idx_triangle_base; i<triangles.size(); i++)
            triangles[i].rebase(idx_base);
      }
   }

   if (remove_redundant_vertices_flag) {
      std::map<unsigned int, std::set<unsigned int> > redundant_map = find_same_vertices(verts);
      // adjust triangles
      remove_redundant_vertices(&triangles, redundant_map);
   }

   return std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > (verts,triangles);
}

std::map<unsigned int, std::set<unsigned int> >
find_same_vertices(const std::vector<glm::vec3> &verts) {

   std::map<unsigned int, std::set<unsigned int> > m;

   const unsigned int semi_edge = 50;
   const unsigned int n_per_side = 2 * semi_edge + 1; // 0.00 and 1.0 are covered
   const unsigned int n_boxes = n_per_side * n_per_side * n_per_side;
   std::vector<std::set<unsigned int> > boxes(n_boxes);

   for (unsigned int i=0; i<verts.size(); i++) {
      const glm::vec3 &vert = verts[i];
      int x_coord = semi_edge + semi_edge * vert.x;
      int y_coord = semi_edge + semi_edge * vert.y;
      int z_coord = semi_edge + semi_edge * vert.z;
      int box_id = n_per_side * n_per_side * x_coord + n_per_side * y_coord + z_coord;
      if (true) {
         if (x_coord <   0) std::cout << "ERROR: x coord " << vert.x << " " << x_coord << std::endl;
         if (x_coord > 100) std::cout << "ERROR: x coord " << vert.x << " " << x_coord << std::endl;
         if (y_coord <   0) std::cout << "ERROR: y coord " << vert.y << " " << y_coord << std::endl;
         if (y_coord > 100) std::cout << "ERROR: y coord " << vert.y << " " << y_coord << std::endl;
         if (z_coord <   0) std::cout << "ERROR: z coord " << vert.z << " " << z_coord << std::endl;
         if (z_coord > 100) std::cout << "ERROR: z coord " << vert.z << " " << z_coord << std::endl;
      }
      boxes[box_id].insert(i);
   }

   float tiny = 0.00001;
   for (unsigned int i=0; i<n_boxes; i++) {
      if (boxes[i].size() > 1) {
         const std::set<unsigned int> &box = boxes[i];
         if (false) { // debug
            std::cout << " Box " << i << " : ";
            std::set<unsigned int>::const_iterator it_1;
            for (it_1=box.begin(); it_1!=box.end(); ++it_1) {
               std::cout << *it_1 << " ";
            }
            std::cout << std::endl;
         }

         // are any of those close to any of the others?
         std::set<unsigned int>::const_iterator it_1;
         std::set<unsigned int>::const_iterator it_2;
         for (it_1=box.begin(); it_1!=box.end(); ++it_1) {
            for (it_2=it_1; it_2!=box.end(); ++it_2) {
               if (it_2 != it_1) {
                  float delta_x = fabs(verts[*it_2].x - verts[*it_1].x);
                  float delta_y = fabs(verts[*it_2].y - verts[*it_1].y);
                  float delta_z = fabs(verts[*it_2].z - verts[*it_1].z);
                  if (delta_x < tiny) {
                     if (delta_y < tiny) {
                        if (delta_z < tiny) {
                           // check here that *it_1 is not in any of the sets already
                           bool add_it = true;
                           std::map<unsigned int, std::set<unsigned int> >::const_iterator it_map;
                           for (it_map=m.begin(); it_map!=m.end(); ++it_map) {
                              const std::set<unsigned int> &s = it_map->second;
                              std::set<unsigned int>::const_iterator it_s;
                              for (it_s=s.begin(); it_s!=s.end(); ++it_s) {
                                 if (*it_s == *it_1)
                                    add_it = false;
                                 break;
                              }
                              if (! add_it)
                                 break;
                           }
                           if (add_it)
                              m[*it_1].insert(*it_2);
                        }
                     }
                  }
               }
            }
         }
      }
   }

   return m;
}

std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> >
make_octasphere(unsigned int num_subdivisions, const glm::vec3 &centre,
                float radius, const glm::vec4 &colour_in, bool remove_redundant_vertices) {

   std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > r;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > vp =
      tessellate_octasphere(num_subdivisions, remove_redundant_vertices);
   r.first.resize(vp.first.size());
   r.second = vp.second;

   for (unsigned int i=0; i<vp.first.size(); i++) {
      r.first[i].pos = vp.first[i];
      r.first[i].pos *= radius;
      r.first[i].pos += centre;
      r.first[i].color = colour_in;
      r.first[i].normal = vp.first[i];
   }

   return r;
}

std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> >
make_octasphere_dish(unsigned int num_subdivisions, const glm::vec3 &centre,
                     float radius, float radiusAlongNormal,
                     const glm::vec3 &dish_normal,
                     const glm::vec4 &colour_in) {

   std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > r;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > vp =
      tessellate_octasphere(num_subdivisions);
   r.first.resize(vp.first.size());
   r.second = vp.second;

   for (unsigned int i=0; i<vp.first.size(); i++) {
      r.first[i].pos = vp.first[i];
      float dp = glm::dot(vp.first[i], dish_normal);
      r.first[i].pos *= radius + radiusAlongNormal * fabs(dp);
      r.first[i].pos += centre;
      r.first[i].color = colour_in;
      r.first[i].normal = vp.first[i];
   }

   return r;
}

ortep_t
tessellate_sphere_sans_octant() {

   ortep_t ortep;

   std::vector<glm::vec3>  &verts     = ortep.vertices;
   std::vector<glm::vec3>  &normals   = ortep.normals;
   std::vector<g_triangle> &triangles = ortep.triangles;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > t = tessellate_octasphere_patch(3);

   const auto &patch_vertices = t.first;

   // std::cout << ":::::::::::::::::::::: t " << t.first.size() << " " << t.second.size() << std::endl;

   // rotate t 4 times to make a hemisphere
   // rotate that 2 times to make a sphere

   glm::vec3 x_axis(1.0f, 0.0f, 0.0f);
   glm::vec3 z_axis(0.0f, 0.0f, 1.0f);
   for (unsigned int ih=0; ih<2; ih++) {
      for (unsigned int irot=0; irot<4; irot++) {
         if (ih==1 && irot==30) {
            // skip this octant (the last one)
         } else {
            // happy path
            float angle = 0.5f * M_PI * static_cast<float>(irot);
            unsigned int idx_base = verts.size();
            unsigned int idx_triangle_base = triangles.size();
            for (unsigned int i=0; i<t.first.size(); i++) {
               glm::vec3 v = glm::rotate(t.first[i], angle, z_axis);
               if (ih==1)
                  v = glm::rotate(v, static_cast<float>(M_PI), x_axis);
               verts.push_back(v);
               normals.push_back(v);
            }
            // add the vertices, but not the triangles for the missing quadrant.
            // (I will use the vertices for the quadrant planes)
            if (ih==1 && irot==3) {
            } else {
               triangles.insert(triangles.end(), t.second.begin(), t.second.end());
            }
            // reindex the triangles
            for (unsigned int i=idx_triangle_base; i<triangles.size(); i++)
               triangles[i].rebase(idx_base);
         }
      }
   }

   // now I want 3 fans, in XY, YZ and XZ to the origin.

   unsigned int ps = patch_vertices.size();

   // ------------------- XZ Plane wedge/quadrant --------------------------

   if (true) {
      unsigned int b = patch_vertices.size() * 4;

      unsigned int idx_base = verts.size();
      unsigned int idx_triangle_base = triangles.size();
      glm::vec3 v0 = glm::vec3(0,0,0);
      glm::vec3 v1 = verts[b  ];  // starts at the north pole
      glm::vec3 v2 = verts[b+1];
      glm::vec3 v3 = verts[b+2];
      glm::vec3 v4 = verts[b+3];
      glm::vec3 v5 = verts[b+4];
      glm::vec3 v6 = verts[b+5];
      glm::vec3 v7 = verts[b+6];
      glm::vec3 v8 = verts[b+7];
      glm::vec3 v9 = verts[b+8];

      verts.push_back(v0);
      verts.push_back(v1);
      verts.push_back(v2);
      verts.push_back(v3);
      verts.push_back(v4);
      verts.push_back(v5);
      verts.push_back(v6);
      verts.push_back(v7);
      verts.push_back(v8);
      verts.push_back(v9);

      for (unsigned int i=0; i<10; i++) normals.push_back(glm::vec3(0,1,0));

      triangles.push_back(g_triangle(0, 1, 2));
      triangles.push_back(g_triangle(0, 2, 3));
      triangles.push_back(g_triangle(0, 3, 4));
      triangles.push_back(g_triangle(0, 4, 5));
      triangles.push_back(g_triangle(0, 5, 6));
      triangles.push_back(g_triangle(0, 6, 7));
      triangles.push_back(g_triangle(0, 7, 8));
      triangles.push_back(g_triangle(0, 8, 9));
      for (unsigned int i=idx_triangle_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);
   }

   // ------------------- YZ Plane wedge/quadrant --------------------------

   {

      unsigned int b = ps * 7;

      unsigned int idx_base = verts.size();
      unsigned int idx_triangle_base = triangles.size();
      glm::vec3 v0 = glm::vec3(0,0,0);
      glm::vec3 v1 = verts[b  ];
      glm::vec3 v2 = verts[b+1];
      glm::vec3 v3 = verts[b+2];
      glm::vec3 v4 = verts[b+3];
      glm::vec3 v5 = verts[b+4];
      glm::vec3 v6 = verts[b+5];
      glm::vec3 v7 = verts[b+6];
      glm::vec3 v8 = verts[b+7];
      glm::vec3 v9 = verts[b+8];

      verts.push_back(v0);
      verts.push_back(v1);
      verts.push_back(v2);
      verts.push_back(v3);
      verts.push_back(v4);
      verts.push_back(v5);
      verts.push_back(v6);
      verts.push_back(v7);
      verts.push_back(v8);
      verts.push_back(v9);

      for (unsigned int i=0; i<10; i++) normals.push_back(glm::vec3(1,0,0));

      triangles.push_back(g_triangle(0, 1, 2));
      triangles.push_back(g_triangle(0, 2, 3));
      triangles.push_back(g_triangle(0, 3, 4));
      triangles.push_back(g_triangle(0, 4, 5));
      triangles.push_back(g_triangle(0, 5, 6));
      triangles.push_back(g_triangle(0, 6, 7));
      triangles.push_back(g_triangle(0, 7, 8));
      triangles.push_back(g_triangle(0, 8, 9));
      for (unsigned int i=idx_triangle_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);
   }

   // ------------------- XY Plane wedge/quadrant --------------------------

   {

      unsigned int idx_base = verts.size();
      unsigned int idx_triangle_base = triangles.size();
      glm::vec3 v0 = glm::vec3(0,0,0);
      glm::vec3 v1 = verts[ 8];
      glm::vec3 v2 = verts[16];
      glm::vec3 v3 = verts[23];
      glm::vec3 v4 = verts[29];
      glm::vec3 v5 = verts[34];
      glm::vec3 v6 = verts[38];
      glm::vec3 v7 = verts[41];
      glm::vec3 v8 = verts[43];
      glm::vec3 v9 = verts[44];

      verts.push_back(v0);
      verts.push_back(v1);
      verts.push_back(v2);
      verts.push_back(v3);
      verts.push_back(v4);
      verts.push_back(v5);
      verts.push_back(v6);
      verts.push_back(v7);
      verts.push_back(v8);
      verts.push_back(v9);

      for (unsigned int i=0; i<10; i++) normals.push_back(glm::vec3(0,0,-1));

      triangles.push_back(g_triangle(0, 1, 2));
      triangles.push_back(g_triangle(0, 2, 3));
      triangles.push_back(g_triangle(0, 3, 4));
      triangles.push_back(g_triangle(0, 4, 5));
      triangles.push_back(g_triangle(0, 5, 6));
      triangles.push_back(g_triangle(0, 6, 7));
      triangles.push_back(g_triangle(0, 7, 8));
      triangles.push_back(g_triangle(0, 8, 9));
      for (unsigned int i=idx_triangle_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);
   }

   ortep.make_lines(); // better (clearer) to factor out the above function into its own function

   return ortep;
}

#include <iomanip>

void
ortep_t::make_lines() {

   // I don't like the way the lines sit over the surface. Maybe a thin box ring would be better.


   const unsigned int n_steps = 64;

   indices_for_lines.reserve(2*n_steps);

   for (unsigned int i=0; i<n_steps; i++) {
      float theta_this = ((0.5 + static_cast<float>(i)) / static_cast<float>(n_steps)) * (2.0 * M_PI);
      float c_this = cosf(theta_this);
      float s_this = sinf(theta_this);
      glm::vec3 pt1(c_this, s_this, 0.0);
      vertices_for_lines.push_back(1.008f * pt1);
   }

   for (unsigned int i=0; i<n_steps; i++) {
      float theta_this = ((0.5 + static_cast<float>(i))) / static_cast<float>(n_steps) * (2.0 * M_PI);
      float c_this = cosf(theta_this);
      float s_this = sinf(theta_this);
      glm::vec3 pt1(c_this, 0.0, -s_this);
      vertices_for_lines.push_back(1.008f * pt1);
   }

   for (unsigned int i=0; i<n_steps; i++) {
      float theta_this = ((0.5 + static_cast<float>(i)) / static_cast<float>(n_steps)) * (2.0 * M_PI);
      float c_this = cosf(theta_this);
      float s_this = sinf(theta_this);
      glm::vec3 pt1(0.0f, c_this, s_this);
      vertices_for_lines.push_back(1.008f * pt1);
   }

   for (unsigned int i=0; i<n_steps; i++) {
      int i_next = i+1;
      if (i == (n_steps-1)) i_next = 0;
      indices_for_lines.push_back(std::pair<unsigned int, unsigned int>(i, i_next));
   }

   for (unsigned int i=0; i<n_steps; i++) {
      int i_next = i+1;
      if (i == (n_steps-1)) i_next = 0;
      indices_for_lines.push_back(std::pair<unsigned int, unsigned int>(i+n_steps, i_next+n_steps));
   }
   for (unsigned int i=0; i<n_steps; i++) {
      int i_next = i+1;
      if (i == (n_steps-1)) i_next = 0;
      indices_for_lines.push_back(std::pair<unsigned int, unsigned int>(i + 2 * n_steps, i_next + 2 * n_steps));
   }

}

void ortep_t::transform(const glm::mat4 &m) {

   for (unsigned int i=0; i<vertices.size(); i++) {
      glm::vec4 v_4_b(       vertices[i], 1.0f);
      glm::vec4 v_4_e(1.4f * vertices[i], 1.0f);
      glm::vec4 t_b = m * v_4_b;
      glm::vec4 t_e = m * v_4_e;
      glm::vec3 t1(t_b);
      glm::vec3 t2(t_e);
      glm::vec3 diff = t2 - t1;
      glm::vec3 n = glm::normalize(diff);
      vertices[i] = t1;
      normals[i] = n;
   }

   for (unsigned int i=0; i<vertices_for_lines.size(); i++) {
      glm::vec4 v(vertices_for_lines[i], 1.0f);
      glm::vec4 t = m * v;
      vertices_for_lines[i] = glm::vec3(t);
   }

}


#if 0
int main(int argc, char **argv) {

   int status = 0;
   unsigned int num_subdivisions = 2;
   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > vv =
      tessellate_octasphere_patch(num_subdivisions);
   for (unsigned int i=0; i<vv.first.size(); i++)
      std::cout << glm::to_string(vv.first[i]) << std::endl;

   for (unsigned int i=0; i<vv.second.size(); i++)
      std::cout  << vv.second[i][0] << " "
                 << vv.second[i][1] << " "
                 << vv.second[i][2] << " " << std::endl;
   return status;

}

#endif
