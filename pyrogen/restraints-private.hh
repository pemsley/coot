/* pyrogen/restraints-private.hh
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

#ifndef RESTRAINTS_PRIVATE_HH
#define RESTRAINTS_PRIVATE_HH

#ifdef LIBCOOTAPI_BUILD // we don't want to use Python for lib
#else
#include <Python.h>
#endif

#include <string>
#include <mmdb2/mmdb_manager.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "geometry/protein-geometry.hh"

#include <lidia-core/rdkit-interface.hh>

namespace coot {

#ifdef LIBCOOTAPI_BUILD
#else
   // private (no SWIG interface)
   //
   // the engine for the above calls
   std::pair<mmdb::Manager *, mmdb::Residue *>
   regularize_inner(PyObject *rdkit_mol,
		    PyObject *restraints_py,
		    const std::string &res_name);

   std::pair<mmdb::Manager *, mmdb::Residue *>
   regularize_inner(RDKit::ROMol &mol,
		    PyObject *restraints_py,
		    const std::string &res_name);
#endif

   void
   regularize_and_update_mol_and_restraints(RDKit::RWMol *mol, dictionary_residue_restraints_t *restraints_p);

   bool is_const_torsion(const RDKit::ROMol &mol,
			 const RDKit::Atom *at_2,
			 const RDKit::Atom *at_3);
   // also private (no interface)
   //
   // now update the atom positions of the rdkit_molecule from residue_p
   // (perhaps this should be in rdkit-interface.hh?)
   //
   void update_coords(RDKit::RWMol *mol, int iconf, mmdb::Residue *residue_p);

   // after we have done minimization, we want to update the coordinates in the dictionary
   void update_chem_comp_atoms_from_residue(mmdb::Residue *residue_p,
                                            dictionary_residue_restraints_t *restraints);

   // alter restraints
   int assign_chirals(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);
   // alter restraints
   void add_chem_comp_atoms(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);
   // alter restraints
   bool add_chem_comp_planes(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints,
			     bool quartet_planes, bool quartet_hydrogen_planes);
   // alter restraints
   void add_chem_comp_aromatic_planes(const RDKit::ROMol &mol,
				      dictionary_residue_restraints_t *restraints,
				      bool quartet_planes, bool quartet_hydrogen_planes);
   // which calls:
   dict_plane_restraint_t add_chem_comp_aromatic_plane_all_plane(const RDKit::MatchVectType &match,
								 const RDKit::ROMol &mol,
								 int plane_id_idx,
								 bool quartet_hydrogen_planes);
   // and
   // modify restraints
   void add_quartet_hydrogen_planes(const RDKit::ROMol &mol,
				    dictionary_residue_restraints_t *restraints);

   // and
   //
   // return the number of added planes
   int add_chem_comp_aromatic_plane_quartet_planes(const RDKit::MatchVectType &match,
						   const RDKit::ROMol &mol,
						   dictionary_residue_restraints_t *restraints,
						   int plane_id_idx);

   // alter restraints
   void add_chem_comp_deloc_planes(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);
   // alter restraints
   void add_chem_comp_sp2_N_planes(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);
   void add_chem_comp_sp2_C_planes(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);


   // alter restraints
   bool fill_with_energy_lib_bonds(const RDKit::ROMol &mol,
				   const energy_lib_t &energy_lib,
				   dictionary_residue_restraints_t *restraints);
   // which calls
   bool add_torsion_to_restraints(dictionary_residue_restraints_t *restraints,
				  const RDKit::ROMol &mol,
				  const RDKit::Atom *at_1,
				  const RDKit::Atom *at_2,
				  const RDKit::Atom *at_3,
				  const RDKit::Atom *at_4,
				  const RDKit::Bond *bond, // between atoms 2 and 3
				  unsigned int *tors_no,
				  unsigned int *const_no,
				  const energy_lib_t &energy_lib);


   // alter restraints
   bool fill_with_energy_lib_angles(const RDKit::ROMol &mol,
				    const energy_lib_t &energy_lib,
				    dictionary_residue_restraints_t *restraints);

   // alter restraints
   bool fill_with_energy_lib_torsions(const RDKit::ROMol &mol,
				      const energy_lib_t &energy_lib,
				      dictionary_residue_restraints_t *restraints);

   std::string convert_to_energy_lib_bond_type(RDKit::Bond::BondType bt);


   int assign_chirals_rdkit_tags(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);
   int assign_chirals_mmcif_tags(const RDKit::ROMol &mol, dictionary_residue_restraints_t *restraints);

   void debug_cip_ranks(const RDKit::ROMol &mol);

   // for returning best graph-match data (for dictionary atom name map to reference)
   // 
   class matching_dict_t {
      mmdb::Residue *residue;
      bool filled_flag;
   public:
      dictionary_residue_restraints_t dict;
      matching_dict_t() {
	 residue = NULL;
	 filled_flag = false;
      }
      matching_dict_t(mmdb::Residue *res, const dictionary_residue_restraints_t &d) : dict(d) {
	 residue = res;
	 filled_flag = true;
      }
      bool filled() const { return filled_flag; }
   };

   

   // Use the pointer to test if the match was successful.
   // Old-style embedded list of test compounds
   matching_dict_t
   match_restraints_to_amino_acids(const dictionary_residue_restraints_t &restraints,
				   mmdb::Residue *residue_p);

   matching_dict_t
   match_restraints_to_reference_dictionaries(const coot::dictionary_residue_restraints_t &restraints,
					      mmdb::Residue *residue_p,
					      const std::vector<std::string> &test_comp_ids,
					      const std::vector<std::string> &test_mmcif_file_names);

   void test_ccp4srs_usage(const coot::dictionary_residue_restraints_t &restraints);



}


#endif // RESTRAINTS_PRIVATE_HH

