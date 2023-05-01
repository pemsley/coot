
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
      sum += v2.color;
      return sum * 0.3333f;
   };

   PyObject *r_py = PyList_New(0);

   const std::vector<g_triangle> &tris = mesh.triangles;
   const std::vector<coot::api::vnc_vertex> &vertices = mesh.vertices;

   if (tris.empty()) return r_py;
   if (vertices.empty()) return r_py;

   std::cout << "DEBUG:: simple_mesh_to_pythonic_mesh(): mesh vertices size " << mesh.vertices.size() << std::endl;
   std::cout << "DEBUG:: simple_mesh_to_pythonic_mesh(): tris size " << tris.size() << std::endl;

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
   return r_py;
}

PyObject *
molecules_container_t::get_pythonic_bonds_mesh(int imol,
                                               const std::string &mode, bool against_a_dark_background,
                                               float bond_width, float atom_radius_to_bond_width_ratio,
                                               int smoothness_factor) {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
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

PyObject *
molecules_container_t::get_pythonic_gaussian_surface_mesh(int imol, float sigma, float contour_level,
                                                          float box_radius, float grid_scale) {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      mesh = molecules[imol].get_gaussian_surface(sigma, contour_level, box_radius, grid_scale);
   }
   return simple_mesh_to_pythonic_mesh(mesh, MULTI_COLOUR);
}

PyObject *
molecules_container_t::get_pythonic_simple_molecule(int imol, const std::string &cid, bool include_hydrogen_atoms_flag) {

   PyObject *r = PyList_New(2);

   coot::simple::molecule_t sm = get_simple_molecule(imol, cid, include_hydrogen_atoms_flag);

   unsigned int n_atoms = sm.atoms.size();
   unsigned int n_bonds = sm.bonds.size();
   PyObject *atom_list = PyList_New(n_atoms);
   PyObject *bond_list = PyList_New(n_bonds);

   for (unsigned int i=0; i<sm.atoms.size(); i++) {
      const auto &atom = sm.atoms[i];
      PyObject *atom_parts_list = PyList_New(5);
      PyObject *pos_py = PyList_New(3);
      PyList_SetItem(pos_py, 0, PyFloat_FromDouble(atom.position.x()));
      PyList_SetItem(pos_py, 1, PyFloat_FromDouble(atom.position.y()));
      PyList_SetItem(pos_py, 2, PyFloat_FromDouble(atom.position.z()));
      PyList_SetItem(atom_parts_list, 0, PyUnicode_FromString(atom.name.c_str()));
      PyList_SetItem(atom_parts_list, 1, PyUnicode_FromString(atom.element.c_str()));
      PyList_SetItem(atom_parts_list, 2, pos_py);
      PyList_SetItem(atom_parts_list, 3, PyLong_FromLong(atom.formal_charge));
      PyList_SetItem(atom_parts_list, 4, PyBool_FromLong(atom.aromatic));
      PyList_SetItem(atom_list, i, atom_parts_list);
   }

   for (unsigned int i=0; i<sm.bonds.size(); i++) {
      const auto &bond = sm.bonds[i];
      PyObject *bond_py = PyList_New(3);
      PyList_SetItem(bond_py, 0, PyLong_FromLong(bond.atom_index_1));
      PyList_SetItem(bond_py, 1, PyLong_FromLong(bond.atom_index_2));
      PyList_SetItem(bond_py, 2, PyLong_FromLong(bond.bond_type));
      PyList_SetItem(bond_list, i, bond_py);
   }
   
   PyList_SetItem(r, 0, atom_list);
   PyList_SetItem(r, 1, bond_list);
   return r;
}


#endif // SWIG


