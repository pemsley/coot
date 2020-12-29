
// spikey mode is 0 for smooth shaded and
// 1 for spikey/flat shaded

#include <map>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#include "make-a-dodec.hh"
#include "utils/dodec.hh"

glm::vec3
clipper_to_glm(const clipper::Coord_orth &co) {
   return glm::vec3(co.x(), co.y(), co.z());
}

std::pair<std::vector<vn_vertex>, std::vector<g_triangle> >
make_dodec() {

   std::pair<std::vector<vn_vertex>, std::vector<g_triangle> > p;
   std::vector<vn_vertex> &vertices = p.first;
   std::vector<g_triangle> &triangles = p.second;

   if (true) {
      pentakis_dodec p_dodec;
      std::vector<clipper::Coord_orth> dodec_vertices_c = p_dodec.d.coords();
      std::vector<clipper::Coord_orth> pyrimid_vertices_c = p_dodec.pyrimid_vertices;

      if (false) {
         for (unsigned int i_face_num=0; i_face_num<12; i_face_num++) {
            std::vector<unsigned int> face = p_dodec.d.face(i_face_num);
            std::cout << "face: i " << i_face_num << " ";
            for (unsigned int j=0; j<face.size(); j++) {
               std::cout << " " << face[j];
            }
            std::cout << std::endl;
            if (i_face_num != 0) continue;
         }
      }

      // kludge some vertices - the vertices are not shared between faces/triangles
      for (unsigned int i_face_num=0; i_face_num<12; i_face_num++) {
         std::vector<unsigned int> face = p_dodec.d.face(i_face_num);
         glm::vec3 sum(0,0,0);
         for (unsigned int j=0; j<5; j++)
            sum += clipper_to_glm(dodec_vertices_c[face[j]]);
         glm::vec3 n = glm::normalize(sum);
         for (unsigned int j=0; j<5; j++) {
            vn_vertex v(clipper_to_glm(dodec_vertices_c[face[j]]), n);
            vertices.push_back(v);
         }
      }

      std::cout << "################# Here with vertices.size() " << vertices.size() << std::endl;
      for (unsigned int i=0; i<vertices.size(); i++)
         std::cout << i << " " << glm::to_string(vertices[i].pos) << std::endl;

      
      for (unsigned int i_face_num=0; i_face_num<12; i_face_num++) {
         unsigned int idx_base = i_face_num * 5;
         for (unsigned int j=1; j<=3; j++) {
            unsigned int index_this = j;
            unsigned int index_next = j + 1;
            g_triangle t(idx_base, idx_base + index_this, idx_base + index_next);
            triangles.push_back(t);
         }
      }
      
      std::cout << "################# Here with triangles.size() " << triangles.size() << std::endl;
      std::cout << "triangles:" << std::endl;
      for (unsigned int i=0; i<triangles.size(); i++) {
         std::cout << "triangle " << i << " :   "
                   << triangles[i].point_id[0] << " "
                   << triangles[i].point_id[1] << " "
                   << triangles[i].point_id[2] << " "
                   << std::endl;
      }
      std::cout << "face 0 vertices: " << std::endl;
      std::cout << " "
                << glm::to_string(vertices[triangles[0].point_id[0]].pos) << " "
                << glm::to_string(vertices[triangles[0].point_id[1]].pos) << " "
                << glm::to_string(vertices[triangles[0].point_id[2]].pos) << " "
                << std::endl;
      std::cout << " "
                << glm::to_string(vertices[triangles[1].point_id[0]].pos) << " "
                << glm::to_string(vertices[triangles[1].point_id[1]].pos) << " "
                << glm::to_string(vertices[triangles[1].point_id[2]].pos) << " "
                << std::endl;
      std::cout << " "
                << glm::to_string(vertices[triangles[2].point_id[0]].pos) << " "
                << glm::to_string(vertices[triangles[2].point_id[1]].pos) << " "
                << glm::to_string(vertices[triangles[2].point_id[2]].pos) << " "
                << std::endl;
   }
   return p;
}

std::pair<std::vector<vn_vertex>, std::vector<g_triangle> >
make_pentakis_dodec(int spikey_mode) {

   std::pair<std::vector<vn_vertex>, std::vector<g_triangle> > p;
   std::vector<vn_vertex> &vertices = p.first;
   std::vector<g_triangle> &triangles = p.second;

   auto normal_from_vertices = [] (const glm::vec3 &pt_0,
                                   const glm::vec3 &pt_1,
                                   const glm::vec3 &pt_2) {
                                  glm::vec3 normal_for_flat_shading = glm::normalize(glm::cross((pt_1-pt_0), (pt_0-pt_2)));
                                  return normal_for_flat_shading;
                               };

   if (spikey_mode == DODEC_SPIKEY_MODE) {
      pentakis_dodec p_dodec(1.5);
      std::vector<clipper::Coord_orth>   dodec_vertices_c = p_dodec.d.coords();
      std::vector<clipper::Coord_orth> pyrimid_vertices_c = p_dodec.pyrimid_vertices;

      if (false) { // debugging
         for (unsigned int i_face_num=0; i_face_num<12; i_face_num++) {
            std::vector<unsigned int> face = p_dodec.d.face(i_face_num);
            std::cout << "face: i " << i_face_num << " ";
            for (unsigned int j=0; j<face.size(); j++) {
               std::cout << " " << face[j];
            }
            std::cout << std::endl;
            if (i_face_num != 0) continue;
         }
      }

      for (unsigned int i_face_num=0; i_face_num<12; i_face_num++) {
         std::vector<unsigned int> face = p_dodec.d.face(i_face_num);
         glm::vec3 point_point = clipper_to_glm(pyrimid_vertices_c[i_face_num]);
         for (unsigned int j=0; j<5; j++) {
            unsigned int index_this = j;
            unsigned int index_next = j + 1;
            if (index_next == 5) index_next = 0;
            glm::vec3 pt_1 = clipper_to_glm(dodec_vertices_c[face[index_this]]);
            glm::vec3 pt_2 = clipper_to_glm(dodec_vertices_c[face[index_next]]);
            glm::vec3 n = normal_from_vertices(point_point, pt_1, pt_2);
            vertices.push_back(vn_vertex(point_point, n));
            vertices.push_back(vn_vertex(pt_1, n));
            vertices.push_back(vn_vertex(pt_2, n));
         }
      }

      for (unsigned int i_face_num=0; i_face_num<12; i_face_num++) {
         unsigned int idx_base = i_face_num * 15;
         for (unsigned int j=0; j<5; j++) {
            unsigned int idx_0 = 3 * j;
            unsigned int idx_1 = 3 * j + 1;
            unsigned int idx_2 = 3 * j + 2;
            if (idx_2 == 6) idx_2 = 1;
            g_triangle t(idx_base + idx_0, idx_base + idx_1, idx_base + idx_2);
            triangles.push_back(t);
         }
      }

   } else {

      // average the normals

      // go away squiggles.
      if (vertices.empty()) std::cout << "" << std::endl;
      if (triangles.empty()) std::cout << "" << std::endl;
   }

   return p;
}
