
// header-here

#ifdef USE_PYTHON


// needed to parse cc-interface.hh
#ifdef USE_GUILE
#include <libguile.h> 
#endif // USE_GUILE


#include "cfc.hh"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/chemical-feature-clusters.hh"
#include "cc-interface.hh"
#include "graphics-info.h"
#endif // MAKE_ENHANCED_LIGAND_TOOLS

// return a Python object that contains (with indexing)
//   water positions around ligand
//   chemical features of ligands
//
// environment_residues_py is a list of residue specs
// solvated_ligand_info_py is a list of
//   list mol_no ligand_spec
// which, with radius will be used to find the waters
// 
PyObject *chemical_feature_clusters_py(PyObject *environment_residues_py,
				       PyObject *solvated_ligand_info_py,
				       float radius) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   PyObject *r = Py_False;

   if (PyList_Check(environment_residues_py)) {
      if (PyList_Check(solvated_ligand_info_py)) {

	 std::vector<coot::residue_spec_t> neighbs =
	    py_to_residue_specs(environment_residues_py);

	 std::vector<coot::chem_feat_solvated_ligand_spec> ligands; // fill this from
         	                                           // solvated_ligand_info_py

	 int n = PyObject_Length(solvated_ligand_info_py);
	 for(int i=0; i<n; i++) {

	    PyObject *o = PyList_GetItem(solvated_ligand_info_py, i);

	    // o should be a list of molecule-idx, residue-spec
	    if (PyList_Check(o)) {
	       int no = PyObject_Length(o);
	       if (no == 2) {
		  PyObject *mol_idx_py  = PyList_GetItem(o, 0);
		  PyObject *lig_spec_py = PyList_GetItem(o, 1);

		  if (PyInt_Check(mol_idx_py)) {
		     int imol = PyInt_AsLong(mol_idx_py);
		     coot::residue_spec_t ligand_spec = residue_spec_from_py(lig_spec_py);

		     if (graphics_info_t::is_valid_model_molecule(imol)) {

			mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;

			std::vector<coot::residue_spec_t> neighbs_waters; // fill this
			std::vector<coot::residue_spec_t> neighbs_raw =
			   coot::residues_near_residue(ligand_spec, mol, radius);
			for (unsigned int i=0; i<neighbs_raw.size(); i++) { 
			   mmdb::Residue *res = coot::util::get_residue(neighbs_raw[i], mol);
			   if (res) {
			      std::string res_name = res->GetResName();
			      if (res_name == "HOH")
				 neighbs_waters.push_back(neighbs_raw[i]);
			   }
			}
			coot::chem_feat_solvated_ligand_spec lig(ligand_spec, neighbs_waters, mol);
			ligands.push_back(lig);
			
		     } else {
			std::cout << "ERROR:: invalid model molecule " << imol << std::endl;
		     }
		     
		  } else {
		     std::cout << "ERROR:: mol_idx is not a int " << std::endl;
		     break;
		  }
	       }
	    }
	 }

	 coot::chem_feat_clust cl(neighbs, ligands, graphics_info_t::Geom_p());

	 std::vector<coot::chem_feat_clust::water_attribs> water_positions =
	    cl.get_water_positions();

	 std::cout << "INFO:: Got " << water_positions.size() << " waters" << std::endl;

	 PyObject *water_attribs_py = PyList_New(water_positions.size());
	 for (unsigned int iw=0; iw<water_positions.size(); iw++) {
	    PyObject *o = PyList_New(3);
	    PyObject *pos_py = PyList_New(3);
	    PyList_SetItem(pos_py, 0, PyFloat_FromDouble(water_positions[iw].pos.x()));
	    PyList_SetItem(pos_py, 1, PyFloat_FromDouble(water_positions[iw].pos.y()));
	    PyList_SetItem(pos_py, 2, PyFloat_FromDouble(water_positions[iw].pos.z()));
	    PyList_SetItem(o, 0, PyInt_FromLong(water_positions[iw].ligand_idx));
	    PyList_SetItem(o, 1, PyInt_FromLong(water_positions[iw].water_spec_idx));
	    PyList_SetItem(o, 2, pos_py);
	    PyList_SetItem(water_attribs_py, iw, o);
	 }

	 r = water_attribs_py;

      } else {
	 std::cout << "ERROR:: chemical_feature_clusters_py() arg 2 is not a list"  << std::endl;
      }
   } else {
      std::cout << "ERROR:: chemical_feature_clusters_py() arg 1 is not a list"  << std::endl;
   } 

   
#endif // MAKE_ENHANCED_LIGAND_TOOLS   

   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;
}




#endif // USE_PYTHON
