
#include "Python.h"
#include <string>
#include <iostream>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Bond.h>

#include <geometry/protein-geometry.hh>
#include <analysis/mogul-interface.hh>

#include "restraints-private.hh"


namespace coot { 

   void mogul_out_to_mmcif_dict(const std::string &mogul_file_name,
				const std::string &comp_id,
				const std::string &compound_name,
				const std::vector<std::string> &atom_names,
				int n_atoms_all,
				int n_atoms_non_hydrogen,
				PyObject *bond_order_restraints_py,
				const std::string &mmcif_out_file_name);

   // return restraints
   PyObject *mogul_out_to_mmcif_dict_by_mol(const std::string &mogul_file_name,
					    const std::string &comp_id,
					    const std::string &compound_name,
					    PyObject *rdkit_mol,
					    PyObject *bond_order_restraints_py,
					    const std::string &mmcif_out_file_name);

   PyObject *
   mmcif_dict_from_mol(const std::string &comp_id,
		       const std::string &compound_name,
		       PyObject *rdkit_mol,
		       const std::string &mmcif_out_file_name);
   
   // which is a wrapper for:
   coot::dictionary_residue_restraints_t
   mmcif_dict_from_mol_inner(const std::string &comp_id,
			     const std::string &compound_name,
			     PyObject *rdkit_mol);
   
   void write_restraints(PyObject *restraints_py,
			 const std::string &monomer_type,
			 const std::string &file_name);

   // return a new rdkit molecule
   void 
   regularize(PyObject *rdkit_mol, PyObject *restraints_py,
	      const std::string &res_name);

   //
   void regularize_and_write_pdb(PyObject *rdkit_mol, PyObject *restraints_py,
				 const std::string &res_name,
				 const std::string &pdb_file_name);

   void write_pdb_from_mol(PyObject *rdkit_mol_py,
			   const std::string &res_name,
			   const std::string &file_name);

   
}

