
#include <iostream>
#include <string>
#include <vector>
#include "coot-utils/cfc.hh"
#include "graphics-info.h"
#include "cfc-2025.hh"

#include "utils/logging.hh"
extern logging logger;

void
chemical_feature_clustering(std::vector<std::pair<int, std::string>> &mol_info_vec) {

   // This function processes the vector of pairs of molecule index and interesting-ligand-type

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
      GtkWidget *dialog = g.cfc_gui.get_dialog();
      if (dialog) {
	 g.cfc_gui.fill_ligands_grid();
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
