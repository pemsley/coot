
#ifdef USE_PYTHON
#include "Python.h"

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
				       float radius);
// scipy has done some clustering
// we get the cluster info here
// 
void chemical_feature_clusters_accept_info_py(PyObject *env_residue,
					      PyObject *mol_ligand_specs,
					      PyObject *cluster_info);

#endif // USE_PYTHON

