
#include <iostream>
#include <vector>

#include <glm/ext.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>

#include "g_triangle.hh"
#include "generic-vertex.hh"

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

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
make_octasphere(unsigned int num_subdivisions, const glm::vec3 &centre,
                float radius, const glm::vec4 &colour_in, bool remove_redundant_vertices) {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > r;

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

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
make_octasphere_dish(unsigned int num_subdivisions, const glm::vec3 &centre,
                     float radius, float radiusAlongNormal,
                     const glm::vec3 &dish_normal,
                     const glm::vec4 &colour_in) {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > r;

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
tessellate_sphere_sans_octant(unsigned int num_subdivisions) {

   ortep_t ortep;

   std::vector<glm::vec3> &verts      = ortep.vertices;
   std::vector<g_triangle> &triangles = ortep.triangles;

   verts.push_back(glm::vec3(0,0,0));

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > t =
      tessellate_octasphere_patch(num_subdivisions);

   // rotate t 4 times to make a hemisphere
   // rotate that 2 times to make a sphere

   glm::vec3 x_axis(1.0f, 0.0f, 0.0f);
   glm::vec3 z_axis(0.0f, 0.0f, 1.0f);
   for (unsigned int ih=0; ih<2; ih++) {
      for (unsigned int irot=0; irot<4; irot++) {
         if (ih==1 && irot==3) {
            // skip this octant
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
            }
            // reindex the triangles
            triangles.insert(triangles.end(), t.second.begin(), t.second.end());
            for (unsigned int i=idx_triangle_base; i<triangles.size(); i++)
               triangles[i].rebase(idx_base);
         }
      }
   }

   // now I want 3 fans, in XY, YZ and XZ to the origin.

   triangles.push_back(g_triangle(0,1,2));

   return ortep;
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
