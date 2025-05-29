
#include <iostream>
#include <string>
#include <vector>
#include "coot-utils/cfc.hh"
#include "analysis/stats.hh"
#include "graphics-info.h"
#include "cfc-2025.hh"

#include "utils/logging.hh"
extern logging logger;

void
chemical_feature_clustering(std::vector<std::pair<int, std::string>> &mol_info_vec) {

   // This function processes the vector of pairs of molecule index and interesting-ligand-type


   auto make_generic_display_objects_for_features = []
      (const std::vector<cfc::typed_cluster_t> &cluster_infos) {

      if (! cluster_infos.empty()) {

	 graphics_info_t g;
	 g.attach_buffers();
	 for (const auto &ci : cluster_infos) {
	    std::string name = ci.make_name();
	    std::string col = "#a0a0a0";
	    std::cout << "ci.family " << ci.family << std::endl;
	    if (ci.family == "Acceptor")     col = "#ddaaaa";
	    if (ci.family == "Donor")        col = "#aaaaff";
	    if (ci.family == "Aromatic")     col = "#88ff88";
	    if (ci.family == "NegIonizable") col = "#ff88aa";
	    if (ci.family == "PosIonizable") col = "#bb88ff";
	    if (ci.family == "ZnBinder")     col = "#bbccff";
	    coot::colour_holder ch = coot::colour_holder_from_colour_name(col);
	    clipper::Coord_orth pt(ci.pos.x, ci.pos.y, ci.pos.z);
	    int object_number = g.new_generic_object_number(name);
	    g.generic_display_objects[object_number].add_pentakis_dodecahedron(ch, col, 1.0, 0.4, pt);
	    g.generic_display_objects[object_number].mesh.setup_buffers();
	    g.set_display_generic_object_simple(object_number, 1);
	 }
      }
   };

   auto get_centre_position = [] (const std::vector<cfc::water_info_t> &wiv) {

      clipper::Coord_orth s(0,0,0);
      unsigned int n = 0;
      for (const auto &wi : wiv) {
	 clipper::Coord_orth pt(wi.pos.x, wi.pos.y, wi.pos.z);
	 s += pt;
	 n++;
      }
      if (n > 0) {
	 double div = static_cast<double>(n);
	 s = clipper::Coord_orth(s.x()/div, s.y()/div, s.z()/div);
      }
      return s;
   };

   auto get_width = [] (const std::vector<cfc::water_info_t> &wiv) {
      int width = 4;
      unsigned int n = 0;
      clipper::Coord_orth s(0,0,0);
      for (const auto &wi : wiv) {
	 clipper::Coord_orth pt(wi.pos.x, wi.pos.y, wi.pos.z);
	 s += pt;
	 n++;
      }
      if (n > 0) {
	 double div = static_cast<double>(n);
	 clipper::Coord_orth middle = clipper::Coord_orth(s.x()/div, s.y()/div, s.z()/div);
	 coot::stats::single s;
	 for (const auto &wi : wiv) {
	    clipper::Coord_orth pt(wi.pos.x, wi.pos.y, wi.pos.z);
	    clipper::Coord_orth delta = pt - middle;
	    double dd = delta.lengthsq();
	    double d = std::sqrt(dd);
	    s.add(d);
	 }
	 double mean = s.mean();
	 width = static_cast<int>(15.0 * mean * 1.4);
      }
      return width;
   };

   auto make_generic_display_objects_for_waters = [get_centre_position, get_width]
      (const std::vector<std::vector<cfc::water_info_t> > &water_infos) {
      for (unsigned int i=0; i<water_infos.size(); i++) {
	 const auto &wi = water_infos[i];
	 if (! wi.empty()) {
	    graphics_info_t g;
	    g.attach_buffers();
	    coot::colour_holder col("#dd2222");
	    int width = get_width(wi);
	    clipper::Coord_orth pt = get_centre_position(wi);
	    meshed_generic_display_object::point_info_t pi(col, pt, width);
	    std::vector<meshed_generic_display_object::point_info_t> piv;
	    piv.push_back(pi); // just one
	    std::string wcn = std::string("Water Cluster ") + std::to_string(i);
	    int object_number = g.new_generic_object_number(wcn);
	    int num_subdivisions = 2;
	    g.generic_display_objects[object_number].add_points(piv, num_subdivisions);
	    g.generic_display_objects[object_number].mesh.setup_buffers();
	    g.set_display_generic_object_simple(object_number, 1);
	 }
      }
   };

   graphics_info_t g;
   if (! g.use_graphics_interface_flag)
      return;

   for (const auto &mol_info : mol_info_vec) {
      std::cout << "DEBUG:: cfc: Molecule Index: " << mol_info.first
		<< ", Ligand Type: " << mol_info.second << std::endl;
   }

   std::vector<cfc::input_info_t> input_infos;
   for (const auto &mol_info : mol_info_vec) {
      int imol = mol_info.first;
      if (g.is_valid_model_molecule(imol)) {
	 mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	 if (mol) {
	    std::string rn = mol_info.second;
	    cfc::input_info_t input_info(mol, imol, rn);
	    input_infos.push_back(input_info);
	 }
      }
   }
   if (! input_infos.empty()) {
      std::pair<std::vector<cfc::typed_cluster_t>, std::vector<std::vector<cfc::water_info_t> > >
	 results = cfc::chemical_feature_clustering(input_infos, *g.Geom_p());
      g.cfc_gui.setup(); // if needed
      g.cfc_gui.cluster_infos = results.first;
      g.cfc_gui.water_infos   = results.second;
      make_generic_display_objects_for_features(g.cfc_gui.cluster_infos);
      make_generic_display_objects_for_waters(g.cfc_gui.water_infos);
      GtkWidget *dialog = g.cfc_gui.get_dialog();
      if (dialog) {
	 g.cfc_gui.fill_ligands_grid();
	 g.cfc_gui.fill_waters_grid();
	 gtk_widget_set_visible(dialog, TRUE);
	 g.set_transient_for_main_window(dialog);
      }
      
   } else {
      logger.log(log_t::DEBUG, "chemical_feature_clustering(): no valid molecules");
   }

}

// a vector of pairs of molecule index and interesting-ligand-type
void
new_chemical_feature_clustering_py(PyObject *mol_infos) {

   if (!PyList_Check(mol_infos)) {
      PyErr_SetString(PyExc_TypeError, "Expected a list of molecule info");
      return;
   } else {
      std::vector<std::pair<int, std::string>> mol_info_vec;
      Py_ssize_t n_mols = PyList_Size(mol_infos);
      for (Py_ssize_t i = 0; i < n_mols; ++i) {
	 PyObject *mol_info = PyList_GetItem(mol_infos, i);
	 if (!PyList_Check(mol_info)) {
	    PyErr_SetString(PyExc_TypeError, "Expected a list of molecule info");
	    return;
	 }
	 // process each molecule info
	 // ...
	 // check that the length of mol_info is 2:
	 Py_ssize_t n_items = PyList_Size(mol_info);
	 if (n_items != 2) {
	    PyErr_SetString(PyExc_ValueError, "Expected a list of two items");
	    return;
	 } else {
	    // get the molecule index
	    PyObject *mol_index_obj = PyList_GetItem(mol_info, 0);
	    if (!PyLong_Check(mol_index_obj)) {
	       PyErr_SetString(PyExc_TypeError, "Expected a long for molecule index");
	       return;
	    } else {
	       long mol_index = PyLong_AsLong(mol_index_obj);
	       if (PyErr_Occurred()) {
		  PyErr_SetString(PyExc_RuntimeError, "Failed to convert molecule index to long");
	       } else {
		  // get the interesting-ligand-type
		  PyObject *ligand_type_obj = PyList_GetItem(mol_info, 1);
		  if (!PyUnicode_Check(ligand_type_obj)) {
		     PyErr_SetString(PyExc_TypeError, "Expected a string for ligand type");
		     return;
		  } else {
		     // convert the ligand type to a string
		     const char *ligand_type = PyUnicode_AsUTF8(ligand_type_obj);
		     if (ligand_type == NULL) {
			PyErr_SetString(PyExc_RuntimeError, "Failed to convert ligand type to string");
		     } else {
			std::string res_type(ligand_type);
			std::pair<int, std::string> mol_info_pair(mol_index, res_type);
			// add the pair to the vector
			mol_info_vec.push_back(mol_info_pair);
		     }
		  }
	       }
	    }
	 }
      }
      // Now we have a vector of pairs of molecule index and interesting-ligand-type
      // Call the C++ function to process this vector
      chemical_feature_clustering(mol_info_vec);
   }
}
