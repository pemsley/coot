
#ifdef SWIG

#include "Python.h"

#include "molecules_container.hh"

PyObject *
molecules_container_t::simple_mesh_to_pythonic_mesh(const coot::simple_mesh_t &mesh) {

   // Blender has colour/material per face. Let's convert our mesh colours and add the
   // colour vect to the returned object.

   PyObject *r_py = PyList_New(0);

   if (true) {
      const std::vector<g_triangle> &tris = mesh.triangles;

      std::cout << "simple_mesh_to_pythonic_mesh(): mesh vertices size " << mesh.vertices.size() << std::endl;
      std::cout << "simple_mesh_to_pythonic_mesh(): tris size " << tris.size() << std::endl;

      r_py = PyList_New(3);
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
         const int &colour_index = tris[i].colour_index;
         PyList_SetItem(tri_py, 0, PyLong_FromLong(tris[i][0]));
         PyList_SetItem(tri_py, 1, PyLong_FromLong(tris[i][1]));
         PyList_SetItem(tri_py, 2, PyLong_FromLong(tris[i][2]));
         PyList_SetItem(tri_py, 3, PyLong_FromLong(colour_index));
         PyList_SetItem(tris_py, i, tri_py);
      }
      for (const auto &col : mesh.colour_index_to_colour_map) {
         PyObject *key_py = PyLong_FromLong(col.first);
         PyObject *colour_py = PyList_New(4);
         for (unsigned int i=0; i<4; i++) {
            PyList_SetItem(colour_py, i, PyFloat_FromDouble(col.second[i]));
         }
         int success = PyDict_SetItem(face_colours_dict_py, key_py, colour_py); // 0 is good
         if (success != 0)
            std::cout << "ERROR:: simple_mesh_to_pythonic_mesh() colour map " << col.first << " " << success << std::endl;
      }
      PyList_SetItem(r_py, 0, vertices_py);
      PyList_SetItem(r_py, 1, tris_py);
      PyList_SetItem(r_py, 2, face_colours_dict_py);
   }
   return r_py;
}

PyObject *
molecules_container_t::get_pythonic_bonds_mesh(int imol) {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      std::string mode("bonds");
      mesh = molecules[imol].get_bonds_mesh(mode, &geom);
   }
   return simple_mesh_to_pythonic_mesh(mesh);
}

PyObject *
molecules_container_t::get_pythonic_model_mesh(int imol, unsigned int mesh_index) {

   coot::simple_mesh_t mesh;
   if (is_valid_model_molecule(imol)) {
      // get MoleculesToTriangles representation here
   }
   return simple_mesh_to_pythonic_mesh(mesh);
}

PyObject *
molecules_container_t::get_pythonic_map_mesh(int imol, float x, float y, float z, float radius, float contour_level) {

   coot::simple_mesh_t mesh;
   clipper::Coord_orth pt(x,y,z); 
   if (is_valid_map_molecule(imol)) {
      mesh = molecules[imol].get_map_contours_mesh(pt, radius, contour_level);
   }
   return simple_mesh_to_pythonic_mesh(mesh);
}

#endif // SWIG
