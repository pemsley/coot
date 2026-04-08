/* pyrogen/restraints-boost.hh
 * 
 * Copyright 2011 by the University of Oxford
 * Copyright 2014 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef RESTRAINTS_HH
#define RESTRAINTS_HH

#ifdef LIBCOOTAPI_BUILD // we don't want to use Python for lib
// no python
#else
#include <Python.h>
#endif

#include <cstddef>
#include <string>
#include <vector>

#include "compat/coot-sysdep.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Bond.h>

#include <geometry/protein-geometry.hh>
#include <analysis/mogul-interface.hh>

#include "restraints-private.hh"

namespace coot {

   class quartet_set {
   public:
      quartet_set() {
	 idx.resize(4,0);
      }
      explicit quartet_set(const std::vector<unsigned int> &q_in) {
	 unsigned int n_max = 4;
	 idx.resize(4,0);
	 if (q_in.size() < n_max)
	    n_max = q_in.size();
	 for (unsigned int i=0; i<n_max; i++)
	    idx[i] = q_in[i];
      }
      std::vector<unsigned int> idx;
      const unsigned int & operator[](unsigned int i) const { return idx[i]; }
   };

   class indexed_name_and_rank_t {
   public:
      unsigned int atom_index;
      unsigned int cip_rank;
      std::string atom_name;
      indexed_name_and_rank_t(unsigned int i, unsigned int r, const std::string &n) : atom_name(n) {
         atom_index = i;
         cip_rank = r;
      }
      bool operator<(const indexed_name_and_rank_t &o) const {
         return (o.cip_rank < cip_rank);
      }
   };

   int get_volume_sign_from_coordinates(const RDKit::ROMol &mol,
                                        unsigned int idx_chiral_centre_atom,
                                        const std::vector<indexed_name_and_rank_t> &neighb_names_and_ranks);

#ifdef LIBCOOTAPI_BUILD
#else
void mogul_out_to_mmcif_dict(const std::string &mogul_file_name,
				const std::string &comp_id,
				const std::string &compound_name,
				const std::vector<std::string> &atom_names,
				int n_atoms_all,
				int n_atoms_non_hydrogen,
				PyObject *bond_order_restraints_py,
				const std::string &mmcif_out_file_name,
				bool quartet_planes, bool quartet_hydrogen_planes);

   // return restraints
   PyObject *mogul_out_to_mmcif_dict_by_mol(const std::string &mogul_file_name,
					    const std::string &comp_id,
					    const std::string &compound_name,
					    PyObject *rdkit_mol,
					    PyObject *bond_order_restraints_py,
					    const std::string &mmcif_out_file_name,
					    bool quartet_planes, bool quartet_hydrogen_planes,
					    bool replace_with_mmff_b_a_restraints=true);

   PyObject *
   mmcif_dict_from_mol(const std::string &comp_id,
		       const std::string &compound_name,
		       PyObject *rdkit_mol,
                       bool do_minimization,
		       const std::string &mmcif_out_file_name,
		       bool quartet_planes, bool quartet_hydrogen_planes,
		       bool replace_with_mmff_b_a_restraints=true);

   // which is a wrapper for:
   std::pair<bool, coot::dictionary_residue_restraints_t>
   mmcif_dict_from_mol_using_energy_lib(const std::string &comp_id,
					const std::string &compound_name,
					PyObject *rdkit_mol,
					bool quartet_planes, bool quartet_hydrogen_planes);

#endif // LIBCOOTAPI_BUILD

   // 2026-01-05 now generate restraint from an rdkit molecule
   std::pair<bool, coot::dictionary_residue_restraints_t>
   mmcif_dict_from_mol_using_energy_lib(const std::string &comp_id,
					const std::string &compound_name,
					const RDKit::ROMol &rdkit_mol,
					bool quartet_planes, bool quartet_hydrogen_planes);

#ifdef LIBCOOTAPI_BUILD
#else
   void write_restraints(PyObject *restraints_py,
			 const std::string &monomer_type,
			 const std::string &file_name);

   void regularize(PyObject *rdkit_mol, PyObject *restraints_py,
		   const std::string &res_name);

   //
   void regularize_and_write_pdb(PyObject *rdkit_mol, PyObject *restraints_py,
				 const std::string &res_name,
				 const std::string &pdb_file_name);

   void write_pdb_from_mol(PyObject *rdkit_mol_py,
			   const std::string &res_name,
			   const std::string &file_name);

   PyObject *types_from_mmcif_dictionary(const std::string &file_name);

   PyObject *test_tuple();
   PyObject *match_restraints_to_dictionaries(PyObject *restraints_in,
					      PyObject *template_comp_id_list,
					      PyObject *template_cif_dict_file_names);
   void write_restraints(PyObject *restraints_py, const std::string &file_name);

#endif

}



#endif // RESTRAINTS_HH

