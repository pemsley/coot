
#ifdef SWIG

#include "Python.h"

#include "molecules_container.hh"

PyObject *
molecules_container_t::simple_mesh_to_pythonic_mesh(const coot::simple_mesh_t &mesh, int mode) {

   // Blender has colour/material per face. Let's convert our mesh colours and add the
   // colour vect to the returned object.

   auto make_colour_hash = [] (const glm::vec4 &c0, const glm::vec4 &c1, const glm::vec4 &c2) {
      double rr = c0.r + 1290.0 * c1.r + 1290.0 * 1290.0 * c2.r;
      double gg = c0.g + 1290.0 * c1.g + 1290.0 * 1290.0 * c2.g;
      double bb = c0.b + 1290.0 * c1.b + 1290.0 * 1290.0 * c2.b;
      double aa = c0.a + 1290.0 * c1.a + 1290.0 * 1290.0 * c2.a;
      double ss = rr + 35.0 * gg + 35.0 * 35.0 * bb + 35.0 * 35.0 * 35.0 * aa;
      long h = static_cast<long>(ss);
      return h;
   };

   auto get_colour_index = [make_colour_hash] (std::map<unsigned long, unsigned int> &colour_index_map,
                                               unsigned int &colour_index,
                                               const glm::vec4 &c0, const glm::vec4 &c1, const glm::vec4 &c2) {

      auto col_hash = make_colour_hash(c0, c1, c2);
      std::map<unsigned long, unsigned int>::const_iterator it = colour_index_map.find(col_hash);
      if (it != colour_index_map.end()) {
         return it->second;
      } else {
         colour_index_map[col_hash] = colour_index;
         colour_index++;
         return colour_index;
      }
   };

   auto average_colour = [] (const coot::api::vnc_vertex &v0,
                             const coot::api::vnc_vertex &v1,
                             const coot::api::vnc_vertex &v2) {
      glm::vec4 sum(0,0,0,0);
      sum += v0.color;
      sum += v1.color;
      sum += v0.color;
      return sum * 0.3333f;
   };

   glm::vec4 c0(0.1, 0.1, 0.1, 0.1);
   glm::vec4 c1(0.4, 0.5, 0.6, 1.0);
   glm::vec4 c2(0.6, 0.5, 0.5, 1.0);
   glm::vec4 c3(1.0, 1.0, 1.0, 1.0);
   long t0 = make_colour_hash(c0, c0, c0);
   long t1 = make_colour_hash(c1, c1, c1);
   long t2 = make_colour_hash(c2, c2, c2);
   long t3 = make_colour_hash(c3, c3, c3);

   std::cout << "t0 " << t0 << std::endl;
   std::cout << "t1 " << t1 << std::endl;
   std::cout << "t2 " << t1 << std::endl;
   std::cout << "t3 " << t1 << std::endl;

   PyObject *r_py = PyList_New(0);

   if (true) {
      const std::vector<g_triangle> &tris      = mesh.triangles;
      const std::vector<coot::api::vnc_vertex> &vertices = mesh.vertices;

      std::cout << "simple_mesh_to_pythonic_mesh(): mesh vertices size " << mesh.vertices.size() << std::endl;
      std::cout << "simple_mesh_to_pythonic_mesh(): tris size " << tris.size() << std::endl;

      std::map<unsigned long, unsigned int> colour_index_map;

      unsigned int colour_index_running = 0;

      r_py = PyList_New(3);

      std::map<int, glm::vec4> colour_index_to_colour_map;

      PyObject *vertices_py     = PyList_New(mesh.vertices.size());
      PyObject *tris_py         = PyList_New(tris.size());
      PyObject *face_colours_dict_py = PyDict_New();
      for (unsigned int i=0; i<mesh.vertices.size(); i++) {
         PyObject *vert_py = PyList_New(3);
         PyList_SetItem(vert_py, 0, PyFloat_FromDouble(mesh.vertices[i].pos[0]));
         PyList_SetItem(vert_py, 1, PyFloat_FromDouble(mesh.vertices[i].pos[1]));
         PyList_SetItem(vert_py, 2, PyFloat_FromDouble(mesh.vertices[i].pos[2]));
         PyList_SetItem(vertices_py, i, vert_py);
      }
      for (unsigned int i=0; i<tris.size(); i++) {
         PyObject *tri_py = PyList_New(4);

         unsigned int colour_index = 0; // SINGLE_COLOUR
         if (mode == MULTI_COLOUR) {
            colour_index = get_colour_index(colour_index_map, // changes reference
                                            colour_index_running, // changes reference
                                            vertices[tris[i][0]].color,
                                            vertices[tris[i][1]].color,
                                            vertices[tris[i][2]].color);
            if (false)
               std::cout << "tri: " << i << " colour_index " << colour_index << std::endl;
            colour_index_to_colour_map[colour_index] =
               average_colour(vertices[tris[i][0]], vertices[tris[i][1]], vertices[tris[i][2]]);
         }
         PyList_SetItem(tri_py, 0, PyLong_FromLong(tris[i][0]));
         PyList_SetItem(tri_py, 1, PyLong_FromLong(tris[i][1]));
         PyList_SetItem(tri_py, 2, PyLong_FromLong(tris[i][2]));
         PyList_SetItem(tri_py, 3, PyLong_FromLong(colour_index));
         PyList_SetItem(tris_py, i, tri_py);
      }

      if (mode == SINGLE_COLOUR) {
         PyObject *colour_py = PyList_New(4);
         for (unsigned int i=0; i<4; i++)
            PyList_SetItem(colour_py, i, PyFloat_FromDouble(vertices[0].color[i]));
         if (! mesh.vertices.empty()) PyDict_SetItem(face_colours_dict_py, PyLong_FromLong(0), colour_py);
      }

      if (mode == MULTI_COLOUR) {
         std::cout << "DEBUG:: in simple_mesh_to_pythonic_mesh: mesh colour_index_to_colour_map size is "
                   << colour_index_to_colour_map.size() << std::endl;

         for (const auto &col : colour_index_to_colour_map) {
            PyObject *key_py = PyLong_FromLong(col.first);
            PyObject *colour_py = PyList_New(4);
            for (unsigned int i=0; i<4; i++)
               PyList_SetItem(colour_py, i, PyFloat_FromDouble(col.second[i]));
            int success = PyDict_SetItem(face_colours_dict_py, key_py, colour_py); // 0 is good
            if (success != 0)
               std::cout << "ERROR:: simple_mesh_to_pythonic_mesh() colour map " << col.first << " " << success << std::endl;
         }
      }

      PyList_SetItem(r_py, 0, vertices_py);
      PyList_SetItem(r_py, 1, tris_py);
      PyList_SetItem(r_py, 2, face_colours_dict_py);
   }
   return r_py;
}

PyObject *
molecules_container_t::get_pythonic_bonds_mesh(int imol,
                                               const std::string &mode, bool against_a_dark_background,
                                               float bond_width, float atom_radius_to_bond_width_ratio,
                                               int smoothness_factor) {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      std::string mode("bonds");
      float ratio = atom_radius_to_bond_width_ratio;
      int sf = smoothness_factor;
      mesh = molecules[imol].get_bonds_mesh(mode, &geom, against_a_dark_background, bond_width, ratio, sf, true, true);
   }
   return simple_mesh_to_pythonic_mesh(mesh, MULTI_COLOUR);
}

PyObject *
molecules_container_t::get_pythonic_map_mesh(int imol, float x, float y, float z, float radius, float contour_level) {

   coot::simple_mesh_t mesh;
   clipper::Coord_orth pt(x,y,z);
   if (is_valid_map_molecule(imol)) {
      mesh = molecules[imol].get_map_contours_mesh(pt, radius, contour_level);
   }
   return simple_mesh_to_pythonic_mesh(mesh, SINGLE_COLOUR);
}

PyObject *
molecules_container_t::get_pythonic_molecular_representation_mesh(int imol, const std::string &atom_selection,
                                                        const std::string &colour_scheme,
                                                        const std::string &style) {
   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      mesh = molecules[imol].get_molecular_representation_mesh(atom_selection, colour_scheme, style);
   }
   return simple_mesh_to_pythonic_mesh(mesh, MULTI_COLOUR);
}

#endif // SWIG
