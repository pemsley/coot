
#include <vector>
#include <glm/ext.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>

#include "g_triangle.hh"
#include "generic-vertex.hh"

#include "prideout-octasphere.hh"

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

std::pair<std::vector<glm::vec3>, std::vector<g_triangle> >
tessellate_octasphere(unsigned int num_subdivisions) {

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

   return std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > (verts,triangles);
}

std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
make_octasphere(unsigned int num_subdivisions, const glm::vec3 &centre,
                float radius, const glm::vec4 &colour_in) {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > r;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > vp =
      tessellate_octasphere(num_subdivisions);
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
