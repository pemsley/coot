/* pyrogen/restraints-boost.cc
 * 
 * Copyright 2011 by the University of Oxford
 * Copyright 2014, 2015 by Medical Research Council
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

#ifdef LIBCOOTAPI_BUILD
#else
#include <boost/python.hpp>
#endif

// #include "rdkit/RDBoost/Exceptions.h"
#include <lidia-core/use-rdkit.hh>
#include "restraints.hh"
#ifdef LIBCOOTAPI_BUILD
#else
#include "py-restraints.hh"
#endif
#include <lidia-core/rdkit-interface.hh>
#include <utils/coot-utils.hh>
#include <coot-utils/coot-coord-utils.hh>

#include "mmff-restraints.hh" // needed?

#include <ideal/simple-restraint.hh>

#ifdef LIBCOOTAPI_BUILD
#else
void
coot::mogul_out_to_mmcif_dict(const std::string &mogul_file_name,
			      const std::string &comp_id,
			      const std::string &compound_name,
			      const std::vector<std::string> &atom_names,
			      int n_atoms_all,
			      int n_atoms_non_hydrogen,
			      PyObject *bond_order_restraints_py,
			      const std::string &cif_file_name,
			      bool quartet_planes, bool quartet_hydrogen_planes) {

   coot::mogul mogul(mogul_file_name);
   coot::dictionary_residue_restraints_t bond_order_restraints =
      monomer_restraints_from_python(bond_order_restraints_py);
   coot::dictionary_residue_restraints_t restraints = mogul.make_restraints(comp_id,
									    compound_name,
									    atom_names,
									    n_atoms_all,
									    n_atoms_non_hydrogen,
									    bond_order_restraints);
   restraints.set_use_nuclear_distances(true);
   restraints.write_cif(cif_file_name);

}
#endif

#ifdef LIBCOOTAPI_BUILD
#else
void
coot::write_restraints(PyObject *restraints_py, const std::string &file_name) {

   coot::dictionary_residue_restraints_t restraints = monomer_restraints_from_python(restraints_py);

   if (false) { // debug
      std::size_t n_atoms = restraints.atom_info.size();
      for (std::size_t iat=0; iat<n_atoms; iat++) {
         std::cout << "write_restraints() " << iat << " " << restraints.atom_info[iat].atom_id_4c  << std::endl;
      }
   }

   if (restraints.is_filled()) {
      restraints.set_use_nuclear_distances(true);
      restraints.write_cif(file_name);
   } else {
      std::cout << "No restraints in write_restraints()" << std::endl;
   }
}
#endif

#ifdef LIBCOOTAPI_BUILD
#else
// replace_with_mmff_b_a_restraints is an optional arg, default true
//
PyObject *
coot::mogul_out_to_mmcif_dict_by_mol(const std::string &mogul_file_name,
				     const std::string &comp_id,
				     const std::string &compound_name,
				     PyObject *rdkit_mol_py,
				     PyObject *bond_order_restraints_py,
				     const std::string &mmcif_out_file_name,
				     bool quartet_planes, bool quartet_hydrogen_planes,
				     bool replace_with_mmff_b_a_restraints) {
   
   // Thanks Uwe H.
   RDKit::ROMol &mol = boost::python::extract<RDKit::ROMol&>(rdkit_mol_py);
   coot::dictionary_residue_restraints_t bond_order_restraints = 
      monomer_restraints_from_python(bond_order_restraints_py);

   // std::cout << "in mogul_out_to_mmcif_dict_by_mol() " << std::endl;
   // debug_rdkit_molecule(&mol);

   mogul mogul(mogul_file_name);
   std::vector<std::string> atom_names;
   unsigned int n_atoms_all = mol.getNumAtoms();
   unsigned int n_atoms_non_hydrogen = 0;

   for (unsigned int iat=0; iat<n_atoms_all; iat++) { 
      const RDKit::Atom *at_p = mol[iat];
      if (at_p->getAtomicNum() != 1)
	 n_atoms_non_hydrogen++;
      try {
	 std::string name = "";
	 at_p->getProp("name", name);
	 atom_names.push_back(name);
      }
      catch (const KeyErrorException &kee) {
	 std::cout << "caught no-name for atom exception in mogul_out_to_mmcif_dict_by_mol(): "
		   <<  kee.what() << std::endl;
      } 
   }

   dictionary_residue_restraints_t restraints; // gets updated
   
   dictionary_residue_restraints_t mogul_restraints =
      mogul.make_restraints(comp_id,
			    compound_name,
			    atom_names,
			    n_atoms_all,
			    n_atoms_non_hydrogen,
			    bond_order_restraints);

   if (replace_with_mmff_b_a_restraints) {

      RDKit::ROMol mol_for_mmff(mol);

      // bonds and angles
      dictionary_residue_restraints_t mmff_restraints = make_mmff_restraints(mol_for_mmff);

      // when status is false, we can return with a partially filled
      // restraints, this is worse than empty (silent failure), so if
      // that's the case replace with empty restraints holder.

      std::pair<bool, dictionary_residue_restraints_t> restraints_local =
	 mmcif_dict_from_mol_using_energy_lib(comp_id, compound_name, rdkit_mol_py,
					      quartet_planes, quartet_hydrogen_planes);

      if (restraints_local.first) {

	 restraints = restraints_local.second;
	 restraints.conservatively_replace_with( mmff_restraints);
	 restraints.conservatively_replace_with(mogul_restraints);

      } else {
	 std::cout << "ERROR:: faliure in mmcif_dict_from_mol_using_energy_lib() "
		   << std::endl;
      }

   } else {

      // dont use MMFF.

      std::pair<bool, dictionary_residue_restraints_t>
	 energy_lib_restraints =
	 mmcif_dict_from_mol_using_energy_lib(comp_id, compound_name, rdkit_mol_py,
					      quartet_planes, quartet_hydrogen_planes);

      if (energy_lib_restraints.first) {
	 restraints = energy_lib_restraints.second;
	 if (replace_with_mmff_b_a_restraints)
	    restraints.conservatively_replace_with(mogul_restraints);
      }
   }
   return monomer_restraints_to_python(restraints);

}
#endif

#ifdef LIBCOOTAPI_BUILD
#else
// return a 3-value tuple:
// 0: success-bool
// 1: new-restraints
// 2: atom-name transformation (from-name, to-name)
//
PyObject *
coot::match_restraints_to_dictionaries(PyObject *restraints_py,
				       PyObject *template_comp_id_list,
				       PyObject *template_cif_dict_file_names) {

   // default return value (failure)
   PyObject *o = PyTuple_New(3);
   PyTuple_SetItem(o, 0, PyBool_FromLong(0));
   PyTuple_SetItem(o, 1, PyLong_FromLong(-1));
   PyTuple_SetItem(o, 2, PyUnicode_FromString(""));

   coot::dictionary_residue_restraints_t restraints = monomer_restraints_from_python(restraints_py);
   std::vector<std::string> comp_ids;
   std::vector<std::string> dictionary_file_names;

   if (PyList_Check(template_comp_id_list)) {
      Py_ssize_t len = PyObject_Length(template_comp_id_list);
      for (Py_ssize_t i=0; i<len; i++) {
	 std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(template_comp_id_list, i)));
	 if (! s.empty())
	    comp_ids.push_back(s);
      }
   }

   if (PyList_Check(template_cif_dict_file_names)) {
      Py_ssize_t len = PyObject_Length(template_cif_dict_file_names);
      for (Py_ssize_t i=0; i<len; i++) {
	 std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(template_cif_dict_file_names, i)));
	 dictionary_file_names.push_back(s);
      }
   }

   if (false)
      std::cout << "debug:: -------- match_restraints_to_dictionaries() calls "
                << "match_restraints_to_reference_dictionaries()" << std::endl;

   mmdb::Residue *dummy_residue_p = NULL;
   matching_dict_t md = match_restraints_to_reference_dictionaries(restraints, dummy_residue_p,
								   comp_ids, dictionary_file_names);

   if (md.filled()) {
      PyObject *name_list_py = PyList_New(md.dict.atom_info.size());
      for (unsigned int i=0; i<md.dict.atom_info.size(); i++)
	 PyList_SetItem(name_list_py, i, PyUnicode_FromString(md.dict.atom_info[i].atom_id_4c.c_str()));
      PyTuple_SetItem(o, 0, PyBool_FromLong(true));
      PyTuple_SetItem(o, 1, monomer_restraints_to_python(md.dict));
      PyTuple_SetItem(o, 2, name_list_py); 
   }

   return o;
}
#endif

// Old function
// 
// (and sugars)
coot::matching_dict_t
coot::match_restraints_to_amino_acids(const coot::dictionary_residue_restraints_t &restraints,
				      mmdb::Residue *residue_p) {

   matching_dict_t dict;

   unsigned int n_comp_ids = 21;
   std::string comp_ids[] = { "CYS", "ASP", "GLU",        "HIS", "ILE", "LYS", "LEU", "MET",
			      "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
			      "G",   "C",   "GLC", "MAN"};

   std::vector<std::string> v(n_comp_ids);
   for (unsigned int i=0; i<n_comp_ids; i++) {
      v[i] = comp_ids[i];
   }

   std::vector<std::string> blank;
   return match_restraints_to_reference_dictionaries(restraints, residue_p, v, blank);

}

void
coot::test_ccp4srs_usage(const coot::dictionary_residue_restraints_t &restraints) {
   // nope
} 


coot::matching_dict_t
coot::match_restraints_to_reference_dictionaries(const coot::dictionary_residue_restraints_t &restraints,
						 mmdb::Residue *residue_p,
						 const std::vector<std::string> &test_comp_ids,
						 const std::vector<std::string> &test_mmcif_file_names) {

   matching_dict_t dict;
   protein_geometry pg;

   if (false) {
      std::cout << "debug:: match_restraints_to_reference_dictionaries() -- start -- " << std::endl;
      std::cout << "debug:: match_restraints_to_reference_dictionaries() with "
                <<  restraints.bond_restraint.size() << " bond restraints" << std::endl;
   }

   pg.set_verbose(false);
   int read_number = 0;

   for (unsigned int i=0; i<test_comp_ids.size(); i++) {
      pg.try_dynamic_add(test_comp_ids[i], i);
   }
   for (unsigned int i=0; i<test_mmcif_file_names.size(); i++) {
      pg.init_refmac_mon_lib(test_mmcif_file_names[i], read_number++);
   }

   std::string out_comp_id = restraints.residue_info.comp_id;
   dictionary_match_info_t best_match;
   int best_idx = -1;

   if (false)
      std::cout << "debug:: match_restraints_to_reference_dictionaries() pg has size "
                << pg.size() << std::endl;

   for (unsigned int idx=0; idx<pg.size(); idx++) { 
      const dictionary_residue_restraints_t &rest = pg.get_monomer_restraints(idx);
      dictionary_match_info_t dmi = 
	 restraints.match_to_reference(rest, NULL, out_comp_id, out_comp_id); // null residue
      std::cout << "    match_restraints_to_reference_dictionaries() testing " 
                << rest.residue_info.comp_id << " made " << dmi.n_matches << " atom matches " 
      	        // not yet << dmi.n_bond_matches << " bond matches"
		<< std::endl;
      if (dmi.n_matches > best_match.n_matches) {
	 best_match = dmi;
	 best_idx = idx;
      }
   }

   if (best_idx >= 0) {
      std::string best_comp_id = pg.get_monomer_restraints(best_idx).residue_info.comp_id;
      std::cout << "INFO:: Matched to reference dictionary of comp-id: " << best_comp_id << std::endl;
      mmdb::Residue *returned_res = util::deep_copy_this_residue(residue_p);
      int imol = 0; // dummy
      std::pair<bool, dictionary_residue_restraints_t> rest =
	 pg.get_monomer_restraints(best_comp_id, imol);
      // dictionary_match_info_t dmi = 
      // restraints.match_to_reference(rest.second, returned_res, out_comp_id, out_comp_id);
      restraints.change_names(returned_res, best_match.name_swaps, best_match.new_comp_id);
      dict = matching_dict_t(returned_res, best_match.dict);
   }

   // here we need to worry about the group - should we return a dict with a group that is L-peptide
   // that doesn't have the atoms for a peptide bond?  No.  In such a case we should change
   // dict.group to be "non-polymer".
   //
   // Let's do that anyway be default.
   //
   dict.dict.residue_info.group = "non-polymer";
   
   return dict;
}

#ifdef LIBCOOTAPI_BUILD
#else
PyObject *
coot::test_tuple() {

   PyObject *o = PyTuple_New(2);
   PyTuple_SetItem(o, 0, PyLong_FromLong(-19));
   PyTuple_SetItem(o, 1, PyUnicode_FromString("this-is-part-of-a-tuple"));

   return o;
}
#endif


#ifdef LIBCOOTAPI_BUILD
#else
PyObject *
coot::types_from_mmcif_dictionary(const std::string &file_name) {

   coot::protein_geometry geom;
   geom.set_verbose(false);
   int read_number = 0;
   geom.init_refmac_mon_lib(file_name, read_number);

   std::vector<std::string> types = geom.monomer_types();

   PyObject *l = PyList_New(types.size());
   for (unsigned int i=0; i<types.size(); i++) { 
      PyObject *o = PyUnicode_FromString(types[i].c_str());
      PyList_SetItem(l, i, o);
   }
   return l;
}
#endif

#ifdef LIBCOOTAPI_BUILD
#else
// This is the function that give pyrogen its name
//
// write restraints and return restraints
// 
// replace_with_mmff_b_a_restraints is an optional arg, default true
// 
PyObject *
coot::mmcif_dict_from_mol(const std::string &comp_id,
			  const std::string &compound_name,
			  PyObject *rdkit_mol_py,
                          bool do_minimization,
			  const std::string &mmcif_out_file_name,
			  bool quartet_planes, bool quartet_hydrogen_planes,
			  bool replace_with_mmff_b_a_restraints) {

   std::pair<bool, coot::dictionary_residue_restraints_t> restraints_pair =
      mmcif_dict_from_mol_using_energy_lib(comp_id, compound_name, rdkit_mol_py,
					   quartet_planes, quartet_hydrogen_planes);

   coot::dictionary_residue_restraints_t &restraints = restraints_pair.second;

   if (false)
      std::cout << "in mmcif_dict_from_mol, mmcif_dict_from_mol_using_energy_lib "
		<< "returns with status " << restraints_pair.first << std::endl;

   if (restraints_pair.first) { 
      if (replace_with_mmff_b_a_restraints) {

	 RDKit::ROMol &mol = boost::python::extract<RDKit::ROMol&>(rdkit_mol_py);
	 RDKit::ROMol mol_for_mmff(mol);
	 // bonds and angles 
	 dictionary_residue_restraints_t mmff_restraints = make_mmff_restraints(mol_for_mmff);
	 restraints.conservatively_replace_with(mmff_restraints);
      }
   } else {
      std::cout << "WARNING:: failure in calling mmcif_dict_from_mol_using_energy_lib() " << std::endl;
   }

   bool success = restraints_pair.first;
   if (success)
      if (! restraints.is_filled()) {
	 std::cout << "WARNING:: restraints are not filled: "
		   << restraints.atom_info.size() << " atoms "
		   << restraints.bond_restraint.size() << " bonds "
		   << std::endl;
	 success = false;
      }

   if (success) {

      if (do_minimization) {
	 RDKit::ROMol &mol = boost::python::extract<RDKit::ROMol&>(rdkit_mol_py);
	 RDKit::RWMol mol_for_mmff(mol);
         regularize_and_update_mol_and_restraints(&mol_for_mmff, &restraints);
      }
      restraints.set_use_nuclear_distances(true);
      restraints.write_cif(mmcif_out_file_name);  // this gets overwritten if dictionary
                                                         // matching is enabled.
      return monomer_restraints_to_python(restraints);
   } else {
      std::cout << "no success in mmcif_dict_from_mol() " << std::endl;
      PyObject *o = new PyObject;
      o = Py_None;
      Py_INCREF(o);
      return o;
   }
}
#endif

#ifdef LIBCOOTAPI_BUILD
#else
// return also success status, true is good
//
std::pair<bool, coot::dictionary_residue_restraints_t>
coot::mmcif_dict_from_mol_using_energy_lib(const std::string &comp_id,
					   const std::string &compound_name,
					   PyObject *rdkit_mol_py,
					   bool quartet_planes, bool quartet_hydrogen_planes) {

   bool status = true;
   coot::dictionary_residue_restraints_t restraints(comp_id, 1);
   RDKit::ROMol &mol = boost::python::extract<RDKit::ROMol&>(rdkit_mol_py);
   std::pair<bool, coot::dictionary_residue_restraints_t> p =
      mmcif_dict_from_mol_using_energy_lib(comp_id, compound_name, mol, quartet_planes, quartet_hydrogen_planes);
   return p;
}
#endif

// return also success status, true is good
//
std::pair<bool, coot::dictionary_residue_restraints_t>
coot::mmcif_dict_from_mol_using_energy_lib(const std::string &comp_id,
					   const std::string &compound_name,
					   const RDKit::ROMol &mol,
					   bool quartet_planes, bool quartet_hydrogen_planes) {

   bool status = true;
   coot::dictionary_residue_restraints_t restraints(comp_id, 1);

   // Was there a user over-ride?
   std::string env_as_string;
   const char *env = getenv("ENERGY_LIB_CIF");

   // To CCP4 standard place then:
   if (! env) {
      const char *env_1 = getenv("CLIBD");
      if (env_1) {
	 env_as_string = std::string(env_1) + "/monomers/ener_lib.cif";
      } else {
	 // Coot standard place then
	 env_as_string = std::string(PKGDATADIR) + "/lib/data/monomers/ener_lib.cif";
      }
   } else {
      env_as_string = env;
   }

   if (env_as_string.empty()) {
      // restraints.is_fillled() is false
      std::cout << "ERROR:: no ENERGY_LIB_CIF env var" << std::endl;
   } else {

      // number of atom first
      //
      unsigned int n_atoms_all = mol.getNumAtoms();
      unsigned int n_atoms_non_hydrogen = 0;
      for (unsigned int iat=0; iat<n_atoms_all; iat++)
	 if (mol[iat]->getAtomicNum() != 1)
	    n_atoms_non_hydrogen++;

      coot::energy_lib_t energy_lib(env_as_string);

      // fill with ener_lib values and then add mogul updates.
      //
      restraints.residue_info.comp_id = comp_id;
      restraints.residue_info.three_letter_code = comp_id;
      restraints.residue_info.name = compound_name;
      restraints.residue_info.number_atoms_all = n_atoms_all;
      restraints.residue_info.number_atoms_nh = n_atoms_non_hydrogen;
      restraints.residue_info.group = "non-polymer";
      restraints.residue_info.description_level = "."; // default is full <smiley>

      coot::add_chem_comp_atoms(mol, &restraints); // alter restraints
      bool status_b = coot::fill_with_energy_lib_bonds(mol, energy_lib, &restraints); // alter restraints
      bool status_a = coot::fill_with_energy_lib_angles(mol, energy_lib, &restraints); // alter restraints
      bool status_t = coot::fill_with_energy_lib_torsions(mol, energy_lib, &restraints); // alter restraints

      int n_chirals = coot::assign_chirals(mol, &restraints); // alter restraints
      if (n_chirals) 
	 restraints.assign_chiral_volume_targets();

      bool status_p = coot::add_chem_comp_planes(mol, &restraints, quartet_planes, quartet_hydrogen_planes);

      // std::cout << "here in mmcif_dict_from_mol_using_energy_lib statuses are "
      //           << status_a << " " << status_b << " " << status_p << std::endl;

      if (! status_b) status = false;
      if (! status_a) status = false;
   }

   // std::cout << "mmcif_dict_from_mol_using_energy_lib returns with status "
   // << status << std::endl;
   std::pair<bool, coot::dictionary_residue_restraints_t> p(status, restraints);
   return p;
}

// return success status - did we find something for all the bonds?
// (executable should fall over if this fails).
// 
bool
coot::fill_with_energy_lib_bonds(const RDKit::ROMol &mol,
				 const coot::energy_lib_t &energy_lib,
				 coot::dictionary_residue_restraints_t *restraints) {

   unsigned int n_bonds = mol.getNumBonds();
   for (unsigned int ib=0; ib<n_bonds; ib++) {
      const RDKit::Bond *bond_p = mol.getBondWithIdx(ib);
      int idx_1 = bond_p->getBeginAtomIdx();
      int idx_2 = bond_p->getEndAtomIdx();
      const RDKit::Atom *at_1 = mol[idx_1];
      const RDKit::Atom *at_2 = mol[idx_2];
      {
	 // put the lighter atom first (so that we find "Hxx ."  rather than "N .")
	 if (at_1->getAtomicNum() > at_2->getAtomicNum())
	    std::swap(at_1, at_2);
	 try {
	    std::string atom_type_1;
	    std::string atom_type_2;
	    std::string atom_name_1;
	    std::string atom_name_2;
	    at_1->getProp("type_energy", atom_type_1);
	    at_2->getProp("type_energy", atom_type_2);
	    at_1->getProp("name", atom_name_1);
	    at_2->getProp("name", atom_name_2);
	    try {
	       std::string bt = convert_to_energy_lib_bond_type(bond_p->getBondType());
	       energy_lib_bond bond =
		  energy_lib.get_bond(atom_type_1, atom_type_2, bt); // add bond type as arg
	       if (0) // or bond.needed_permissive
		  std::cout << "....... " << atom_name_1 << " " << atom_name_2 << " types \""
			    << atom_type_1 << "\" \"" << atom_type_2
			    << "\" got bond " << bond << " with permissive search " << std::endl;
	       std::string bond_type = bond.type;
	       dict_bond_restraint_t bondr(atom_name_1, atom_name_2, bond_type, bond.length, bond.esd, 0.0, 0.0, false);
	       restraints->bond_restraint.push_back(bondr);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "ERROR::   runtime_error when adding bond restraint for bond number "
			 << ib << " atom-names: " << atom_name_1 << " " << atom_name_2 << " "
			 << rte.what() << std::endl;
	    } 
	 
	 }
	 catch (const KeyErrorException &kee) {
	    std::cout << "ERROR:: caught KeyErrorException in fill_with_energy_lib_bonds() - "
		      << "atom types and names for bond number " << ib << std::endl;
	 }
      }
   }
   // std::cout << "returnging form fill_with_energy_lib_bonds() " << n_bonds << std::endl;
   return (n_bonds == restraints->bond_restraint.size());
}

// return success status (executable should fall over if this fails).
bool
coot::fill_with_energy_lib_angles(const RDKit::ROMol &mol,
				  const coot::energy_lib_t &energy_lib,
				  coot::dictionary_residue_restraints_t *restraints) {
   
   unsigned int n_atoms = mol.getNumAtoms();
   std::map<std::string, bool> done_angle;
   for (unsigned int iat_1=0; iat_1<n_atoms; iat_1++) { 
      const RDKit::Atom *at_1 = mol[iat_1];
      RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
      boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_1);
      while(nbr_idx_1 != end_nbrs_1){
	 const RDKit::Atom *at_2 = mol[*nbr_idx_1];

	 RDKit::ROMol::ADJ_ITER nbr_idx_2, end_nbrs_2;
	 boost::tie(nbr_idx_2, end_nbrs_2) = mol.getAtomNeighbors(at_2);
	 while(nbr_idx_2 != end_nbrs_2){
	    const RDKit::Atom *at_3 = mol[*nbr_idx_2];
	    if (at_3 != at_1) { 

	       try {
		  std::string atom_type_1;
		  std::string atom_type_2;
		  std::string atom_type_3;
		  std::string atom_name_1;
		  std::string atom_name_2;
		  std::string atom_name_3;
		  at_1->getProp("type_energy", atom_type_1);
		  at_2->getProp("type_energy", atom_type_2);
		  at_3->getProp("type_energy", atom_type_3);
		  at_1->getProp("name", atom_name_1);
		  at_2->getProp("name", atom_name_2);
		  at_3->getProp("name", atom_name_3);

		  try {

		     std::string dash("-");
		     std::string angle_key_name_1 = atom_name_1 + dash + atom_name_2 + dash + atom_name_3;
		     std::string angle_key_name_2 = atom_name_3 + dash + atom_name_2 + dash + atom_name_1;

		     if (done_angle.find(angle_key_name_1) == done_angle.end() &&
			 done_angle.find(angle_key_name_2) == done_angle.end()) { 
		     
			energy_lib_angle angle =
			   energy_lib.get_angle(atom_type_1, atom_type_2, atom_type_3);
			
			dict_angle_restraint_t angler(atom_name_1, atom_name_2, atom_name_3,
						      angle.angle, angle.angle_esd);

			restraints->angle_restraint.push_back(angler);
			done_angle[angle_key_name_1] = true;
			done_angle[angle_key_name_2] = true;
		     }
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "WARNING:: error in adding angle restraint for atoms "
			       << at_1->getIdx() << " "
			       << at_2->getIdx() << " "
			       << at_3->getIdx() << " "
			       << rte.what() << std::endl;
		  } 
	       }
	       catch (const KeyErrorException &kee) {
		  std::cout << "WARNING:: caught KeyErrorException in fill_with_energy_lib_angles() "
			    << std::endl;
	       }
	       
	       
	    }
	    ++nbr_idx_2;
	 }
	 
	 ++nbr_idx_1;
      }
   }
   return true; // placeholder
}


// return success status (executable should fall over if this fails).
bool
coot::fill_with_energy_lib_torsions(const RDKit::ROMol &mol,
				    const coot::energy_lib_t &energy_lib,
				    coot::dictionary_residue_restraints_t *restraints) {

   bool status = true;
   unsigned int n_atoms = mol.getNumAtoms();
   unsigned int tors_no = 1; // incremented on addition
   unsigned int const_no = 1; // incremented on addition.  When const_no is incremented, tors_no is not.
   std::map<std::string, bool> done_torsion;
   bool debug = false;
   
   for (unsigned int iat_1=0; iat_1<n_atoms; iat_1++) { 
      const RDKit::Atom *at_1 = mol[iat_1];

      RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
      boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_1);
      while(nbr_idx_1 != end_nbrs_1){
	 const RDKit::Atom *at_2 = mol[*nbr_idx_1];

	 RDKit::ROMol::ADJ_ITER nbr_idx_2, end_nbrs_2;
	 boost::tie(nbr_idx_2, end_nbrs_2) = mol.getAtomNeighbors(at_2);
	 while(nbr_idx_2 != end_nbrs_2){
	    const RDKit::Atom *at_3 = mol[*nbr_idx_2];
	    if (at_3 != at_1) {
	       
	       RDKit::ROMol::ADJ_ITER nbr_idx_3, end_nbrs_3;
	       boost::tie(nbr_idx_3, end_nbrs_3) = mol.getAtomNeighbors(at_3);

	       // Is there another neighbour of 3 that is not a
	       // hydrogen?  If so, we prefer that.  If not, we go
	       // with the hydrogen atom.

	       bool at_4_set = false;
	       const RDKit::Atom *at_4 = mol[*nbr_idx_3]; // best so far, (maybe its at_2 though)
	       if (at_4 != at_2 && at_4 != at_1)
		  at_4_set = true; // OK, it wasn't.
	       
	       while (nbr_idx_3 != end_nbrs_3) {

		  const RDKit::Atom *at_4_trial = mol[*nbr_idx_3];
		  if (at_4_trial != at_2 && at_4_trial != at_1) {
		     if (at_4_trial->getAtomicNum() != 1) {
			// anything not hydrogen is good enough.
			at_4 = at_4_trial;
			at_4_set = true;
			break;
		     }
		  }
		  ++nbr_idx_3++;
	       }

	       if (at_4_set) {
		  
		  try { 

		     std::string atom_name_2;
		     std::string atom_name_3;
		     at_2->getProp("name", atom_name_2);
		     at_3->getProp("name", atom_name_3);

		     // if we have not done this atom-2 <-> atom-3 torsion before...
		     //
		     std::string torsion_key_name_23;
		     std::string torsion_key_name_32;
		     torsion_key_name_23  = atom_name_2;
		     torsion_key_name_23 += "-";
		     torsion_key_name_23 += atom_name_3;
			
		     torsion_key_name_32  = atom_name_3;
		     torsion_key_name_32 += "-";
		     torsion_key_name_32 += atom_name_2;

		     bool done_this_already = true;
		     if (done_torsion.find(torsion_key_name_23) == done_torsion.end() &&
			 done_torsion.find(torsion_key_name_32) == done_torsion.end())
			done_this_already = false;

		     if (debug) { 
			std::cout << "considering torsion keys \"" << torsion_key_name_23
				  << "\" and \"" << torsion_key_name_32 << "\"";
			if (done_this_already)
			   std::cout << " done already" << std::endl;
			else
			   std::cout << std::endl;
		     }

		     if (! done_this_already) {

			const RDKit::Bond *bond = mol.getBondBetweenAtoms(*nbr_idx_1, *nbr_idx_2);
			bool success = add_torsion_to_restraints(restraints, mol,
								 at_1, at_2, at_3, at_4,
								 bond,
								 &tors_no, &const_no,
								 energy_lib);
			   
			if (success) {
			   done_torsion[torsion_key_name_23] = true;
			   done_torsion[torsion_key_name_32] = true;
			}
		     }
		  }
		  catch (const KeyErrorException &kee) {
		     std::cout << "WARNING:: caught KeyErrorException in fill_with_energy_lib_torsions() "
			       << std::endl;
		  }
	       }
	       
	    }
	    ++nbr_idx_2;
	 }
	 ++nbr_idx_1;
      }
   }
   return status;
}


bool
coot::add_torsion_to_restraints(coot::dictionary_residue_restraints_t *restraints,
				const RDKit::ROMol &mol,
				const RDKit::Atom *at_1,
				const RDKit::Atom *at_2,
				const RDKit::Atom *at_3,
				const RDKit::Atom *at_4,
				const RDKit::Bond *bond, // between atoms 2 and 3
				unsigned int *tors_no,
				unsigned int *const_no,
				const coot::energy_lib_t &energy_lib) {

   bool added_state = false;
   bool debug = false;
   
   try {
      std::string atom_type_1;
      std::string atom_type_2;
      std::string atom_type_3;
      std::string atom_type_4;
      std::string atom_name_1;
      std::string atom_name_2;
      std::string atom_name_3;
      std::string atom_name_4;
      at_1->getProp("type_energy", atom_type_1);
      at_2->getProp("type_energy", atom_type_2);
      at_3->getProp("type_energy", atom_type_3);
      at_4->getProp("type_energy", atom_type_4);
      at_1->getProp("name", atom_name_1);
      at_2->getProp("name", atom_name_2);
      at_3->getProp("name", atom_name_3);
      at_4->getProp("name", atom_name_4);


      if (debug) 
	 std::cout << "torsion-atoms..... ids: "
		   << at_1->getIdx() << " "
		   << at_2->getIdx() << " "
		   << at_3->getIdx() << " "
		   << at_4->getIdx() << "    names: "
		   << atom_name_1 << " " 
		   << atom_name_2 << " " 
		   << atom_name_3 << " " 
		   << atom_name_4 << " " 
		   << std::endl;

      // some of the time we may try to get a torsion that does not
      // correspond to anything in the dictionary (with the given
      // atom types).  In that case, we will catch a runtime_error.
      // What to do then is not clear to me yet.
      // 
      try {
	 if (debug)
	    std::cout << "trying to get torsion-from-lib for " << atom_type_2
		      << "---" << atom_type_3 << std::endl;
	 energy_lib_torsion tors =
	    energy_lib.get_torsion(atom_type_2, atom_type_3);
			      
	 bool is_const = is_const_torsion(mol, at_2, at_3);

	 if (debug)
	    std::cout << "    torsion between a " << atom_type_2 << " and a "
		      << atom_type_3 << " gave torsion " << tors << "  is_const: "
		      << is_const << std::endl;

	 if (tors.period != 0) { 
	    double esd = 20.0;
	    std::string tors_id;
	    if (! is_const) { 
	       tors_id = "var_";
	       tors_id += util::int_to_string(*tors_no);
	       (*tors_no)++;
	    } else {
	       tors_id = "CONST_";
	       tors_id += util::int_to_string(*const_no);
	       (*const_no)++;
	    }
	    dict_torsion_restraint_t torsionr(tors_id,
					      atom_name_1, atom_name_2,
					      atom_name_3, atom_name_4,
					      tors.angle, esd, tors.period);
	    restraints->torsion_restraint.push_back(torsionr);
	    added_state = true;
	 }
      }
      catch (const std::runtime_error &rte) {

	 if (debug)
	    std::cout << "failed  to get torsion-from-lib for " << atom_type_2
		      << " to " << atom_type_3 << std::endl;
			      
	 // default
	 double angle = 180;
	 double esd = 20;
	 int period = 1;
			      
	 bool is_const = is_const_torsion(mol, at_2, at_3);
	 RDKit::Atom::HybridizationType ht_2 = at_2->getHybridization();
	 RDKit::Atom::HybridizationType ht_3 = at_3->getHybridization();

	 if (ht_2 == RDKit::Atom::SP3 || ht_3 == RDKit::Atom::SP3) {
	    period = 3;
	    angle = 60;
	    if (is_const) esd = 2;
	 }

	 if ((ht_2 == RDKit::Atom::SP2 && ht_3 == RDKit::Atom::SP3) ||
	     (ht_2 == RDKit::Atom::SP3 && ht_3 == RDKit::Atom::SP2)) {
	    period = 2;
	    angle = 90;
	    if (is_const) esd = 2;
	 }


	 if (ht_2 == RDKit::Atom::SP2 && ht_3 == RDKit::Atom::SP2) {
	    period = 2;
	    if (is_const) esd = 2;

	    // is this a forced cis or trans bond though?
	    //
	    // Note that currently (20151011) wwPDB CCD files do not
	    // contain stereo information in the bond descriptions, so
	    // bonds will not be constructed as STEREOE or STEREOZ.
	    // (SMILES will though).

	    if (bond) {
	       RDKit::Bond::BondType bt = bond->getBondType();
	       if (bt == RDKit::Bond::DOUBLE) {
		  RDKit::Bond::BondStereo st = bond->getStereo();
		  if (st == RDKit::Bond::STEREOE) { // trans double bond;
		     period = 1;
		  }
		  if (st == RDKit::Bond::STEREOZ) { // cis double bond;
		     period = 1;
		     angle = 0;
		  }
	       }
	    }
	 } 

	 if (ht_2 == RDKit::Atom::SP || ht_3 == RDKit::Atom::SP) {
	    is_const = 1;
	    esd = 2;
	 }

	 std::string tors_id("var_");
	 if (is_const) {
	    tors_id="CONST_";
	    tors_id+=util::int_to_string(*const_no);
	    (*const_no)++;
	 } else {
	    tors_id+=util::int_to_string(*tors_no);
	    (*tors_no)++;
	 }

	 dict_torsion_restraint_t torsionr(tors_id,
					   atom_name_1, atom_name_2,
					   atom_name_3, atom_name_4,
					   angle, esd, period);
	 if (debug)
	    std::cout << "fallback pushback::" << torsionr << std::endl;
	 restraints->torsion_restraint.push_back(torsionr);
	 added_state = true;
			      
      }
   }
   catch (const KeyErrorException &kee) {
      std::cout << "WARNING:: caught KeyErrorException in fill_with_energy_lib_torsions() "
		<< std::endl;
   }
   return added_state;
}


bool
coot::is_const_torsion(const RDKit::ROMol &mol,
		       const RDKit::Atom *torsion_at_2,
		       const RDKit::Atom *torsion_at_3) {

   // is the bond between the 2 central atoms of the torsion a double bond?
   // (or triple or atromatic).  If so, then this is a const torsion

   bool status = false;
   
   unsigned int n_bonds = mol.getNumBonds();
   for (unsigned int ib=0; ib<n_bonds; ib++) {
      const RDKit::Bond *bond_p = mol.getBondWithIdx(ib);
      RDKit::Atom *bond_at_1 = bond_p->getBeginAtom();
      RDKit::Atom *bond_at_2 = bond_p->getEndAtom();

      bool found_torsion_bond = false;
      if (torsion_at_2 == bond_at_1)
	 if (torsion_at_3 == bond_at_2)
	    found_torsion_bond = true;
      if (torsion_at_2 == bond_at_2)
	 if (torsion_at_3 == bond_at_1)
	    found_torsion_bond = true;

      if (found_torsion_bond) { 
	 if (bond_p->getBondType() == RDKit::Bond::AROMATIC)    status = true;
	 if (bond_p->getBondType() == RDKit::Bond::DOUBLE)      status = true;
	 if (bond_p->getBondType() == RDKit::Bond::TRIPLE)      status = true;
	 if (bond_p->getBondType() == RDKit::Bond::QUADRUPLE)   status = true;
	 if (bond_p->getBondType() == RDKit::Bond::ONEANDAHALF) status = true;
	 if (bond_p->getBondType() == RDKit::Bond::TWOANDAHALF) status = true;
	 break;
      }
   }

   return status;

}

void
coot::update_chem_comp_atoms_from_residue(mmdb::Residue *residue_p,
                                          coot::dictionary_residue_restraints_t *restraints) {

   bool verbose = false;
   if (verbose)
      std::cout << "update_chem_comp_atoms_from_residue() ******************************" << std::endl;
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         std::string atom_name(at->GetAtomName());

         std::vector<dict_atom> &atom_info = restraints->atom_info;;
         for (unsigned int jat=0; jat<atom_info.size(); jat++) {
            dict_atom &da = atom_info[jat];
            if (da.atom_id_4c == atom_name) {
               clipper::Coord_orth c = co(at);
               const clipper::Coord_orth &old_c = da.model_Cartn.second;
               if (verbose)
                  std::cout << "debug:: updating " << atom_name << " from " << old_c.format()
                            << " to " << c.format() << std::endl;
               da.model_Cartn.second = c;
            }
         }
      }
   }
}


void
coot::add_chem_comp_atoms(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints) {

   int iconf = 0;
   unsigned int n_atoms = mol.getNumAtoms();
   for (unsigned int iat=0; iat<n_atoms; iat++) { 
      const RDKit::Atom *at_p = mol[iat];
      try {
	 std::string name;
	 std::string atom_type;
	 double charge;
	 bool have_charge = true; // can be clever with GetProp() KeyErrorException
	                          // if you like
	 // std::cout << "add_chem_comp_atoms() getting names... " << std::endl;
	 at_p->getProp("name", name);
	 // std::cout << "add_chem_comp_atoms() getting type_energy... " << std::endl;
	 at_p->getProp("type_energy", atom_type);
	 // std::cout << "add_chem_comp_atoms() getting charges... " << std::endl;
	 at_p->getProp("_GasteigerCharge", charge);
	 std::pair<bool, float> charge_pair(have_charge, charge);
	 dict_atom atom(name, name, at_p->getSymbol(), atom_type, charge_pair);
 	 RDKit::Conformer conf = mol.getConformer(iconf);
 	 RDGeom::Point3D &r_pos = conf.getAtomPos(iat);
 	 clipper::Coord_orth pos(r_pos.x, r_pos.y, r_pos.z);
	 atom.model_Cartn = std::pair<bool, clipper::Coord_orth> (true, pos);	 
	 restraints->atom_info.push_back(atom);
      }
      catch (const KeyErrorException &kee) {
	 std::cout << "WARNING:: caught property exception in add_chem_comp_atoms()"
		   << iat << std::endl;
      }
   }
}

// what fun!
// C++ smarts
bool
coot::add_chem_comp_planes(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints,
			   bool quartet_planes, bool quartet_hydrogen_planes) {

   bool status = true;
   add_chem_comp_aromatic_planes(mol, restraints, quartet_planes, quartet_hydrogen_planes);
   add_chem_comp_deloc_planes(mol, restraints);
   add_chem_comp_sp2_C_planes(mol, restraints);
   restraints->remove_redundant_plane_restraints();
   restraints->reweight_subplanes();
   add_chem_comp_sp2_N_planes(mol, restraints);
   return status;
}

// what fun!
// C++ smarts
void
coot::add_chem_comp_aromatic_planes(const RDKit::ROMol &mol,
				    coot::dictionary_residue_restraints_t *restraints,
				    bool quartet_planes, bool quartet_hydrogen_planes) {

   std::vector<std::string> patterns;

   // I am not sure that fused ring systems are a good idea.
   // 
   //    patterns.push_back("a12aaaaa1aaaa2");

   patterns.push_back("a12aaaaa1aaa2"); // 6-5 is OK, I think
   
   patterns.push_back("a1aaaaa1");
   patterns.push_back("a1aaaa1");
   patterns.push_back("[*;^2]1[*;^2][*;^2][A;^2][*;^2]1"); // non-aromatic 5-ring

   int plane_id_idx = 1; 
   for (unsigned int ipat=0; ipat<patterns.size(); ipat++) {
      RDKit::ROMol *query = RDKit::SmartsToMol(patterns[ipat]);
      std::vector<RDKit::MatchVectType>  matches;
      bool recursionPossible = true;
      bool useChirality = true;
      bool uniquify = true;
      int matched = RDKit::SubstructMatch(mol,*query,matches,uniquify,recursionPossible, useChirality);
      // int matched = false; //20210923-PE FIXME
      for (unsigned int imatch=0; imatch<matches.size(); imatch++) { 
	 if (matches[imatch].size() > 0) {

	    if (false) { // too noisy when debugging other things
	       std::cout << "INFO:: matched aromatic plane: " << std::setw(14) << std::right
			 << patterns[ipat];
	       std::cout << " (";
	       for (unsigned int iat=0; iat<matches[imatch].size(); iat++) {
		  unsigned int atom_idx = matches[imatch][iat].second;
		  try {
		     const RDKit::Atom *at_p = mol[atom_idx];
		     std::string atom_name;
		     at_p->getProp("name", atom_name);
		     std::cout << " " << atom_name;
		  }
		  catch (const KeyErrorException &kee) {
		     std::cout << " " << atom_idx;
		  }
	       }
	       std::cout << " )";
	       std::cout << std::endl;
	    }

	    if (! quartet_planes) { 
	       dict_plane_restraint_t plr =
		  add_chem_comp_aromatic_plane_all_plane(matches[imatch], mol,
							 plane_id_idx,
							 quartet_hydrogen_planes);
	       if (! plr.empty()) {
		  restraints->plane_restraint.push_back(plr);
		  plane_id_idx++;
	       }
	    } else {
	       // Don't add hydrogen quartets (that's done later)
	       int n_added =
		  add_chem_comp_aromatic_plane_quartet_planes(matches[imatch], mol, restraints, plane_id_idx);
	       plane_id_idx += n_added;
	    } 
	 }
      }
   }

   if (quartet_hydrogen_planes || quartet_planes) {
      add_quartet_hydrogen_planes(mol, restraints);
   }
}

// modify restraints
void
coot::add_quartet_hydrogen_planes(const RDKit::ROMol &mol,
				  coot::dictionary_residue_restraints_t *restraints) { 

   int h_plane_quartet_id_idx = 1; // for the first
      
   // Find hydrogens that are connected to an sp2 atom and make a
   // plane of the sp2 atom and its neighbours (including this
   // hydrogen of course).
   unsigned int n_atoms = mol.getNumAtoms();
   for (unsigned int iat_1=0; iat_1<n_atoms; iat_1++) { 
      const RDKit::Atom *at_1 = mol[iat_1];
      if (at_1->getAtomicNum() == 1) {
	 std::vector<unsigned int> quartet_indices;

	 RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
	 boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_1);
	 while(nbr_idx_1 != end_nbrs_1){
	    const RDKit::Atom *at_centre = mol[*nbr_idx_1];
	       
	    if (at_centre->getHybridization() == RDKit::Atom::SP2) {

	       quartet_indices.push_back(*nbr_idx_1); // the idx of atom to which the H is connected
	       RDKit::ROMol::ADJ_ITER nbr_idx_2, end_nbrs_2;
	       boost::tie(nbr_idx_2, end_nbrs_2) = mol.getAtomNeighbors(at_centre);
	       while(nbr_idx_2 != end_nbrs_2){
		  quartet_indices.push_back(*nbr_idx_2);
		  ++nbr_idx_2;
	       }
	    }
	    ++nbr_idx_1;
	 }

	 // OK! We found a H-quartet. Add it.

	 if (quartet_indices.size() == 4) {
	    try {
	       std::vector<std::string> quartet_names;
	       for (unsigned int jj=0; jj<quartet_indices.size(); jj++) {
		  std::string name;
		  mol[quartet_indices[jj]]->getProp("name", name);
		  if (0)
		     std::cout << "Quartet " << h_plane_quartet_id_idx << ": "
			       << quartet_indices[jj] << " " << name << std::endl;
		  quartet_names.push_back(name);
	       }
	       std::string quartet_plane_id = "H-quartet-" + util::int_to_string(h_plane_quartet_id_idx);

	       if (0) { // debug
		  std::cout << "Adding plane " << quartet_plane_id << " with atoms ";
		  for (unsigned int iat=0; iat<quartet_names.size(); iat++) { 
		     std::cout << quartet_names[iat] << " ";
		  }
		  std::cout << std::endl;
	       }

	       
	       double dist_esd = 0.02;
	       coot::dict_plane_restraint_t rest(quartet_plane_id, quartet_names, dist_esd);
	       restraints->plane_restraint.push_back(rest);
	       h_plane_quartet_id_idx++;
	    }
	    catch (const KeyErrorException &kee) {
	       std::cout << "Badness missing atom name in H-quartet" << std::endl;
	    }
	 }
      } 
   } 
}

coot::dict_plane_restraint_t
coot::add_chem_comp_aromatic_plane_all_plane(const RDKit::MatchVectType &match,
					     const RDKit::ROMol &mol,
					     int plane_id_idx,
					     bool quartet_hydrogen_planes) {

   coot::dict_plane_restraint_t plane_restraint; // returend value, empty initially

   std::string plane_id = "plane-arom-" + util::int_to_string(plane_id_idx);
   std::vector<std::string> plane_restraint_atoms; 
   try {
      for (unsigned int ii=0; ii<match.size(); ii++) {
	 const RDKit::Atom *at_p = mol[match[ii].second];

	 // only add this atom to a plane restraint if it not
	 // already in a plane restraint.  Test by failing to
	 // get the plane_id property.

	 bool add_atom_to_plane = true;

	 if ((at_p->getAtomicNum() != 1) || !quartet_hydrogen_planes) {

	    // add the atom to the plane if the plane that it is
	    // already in is not this plane.
	    // 
	    try {
	       std::string atom_plane;
	       at_p->getProp("plane_id", atom_plane);
	       if (atom_plane == plane_id)
		  add_atom_to_plane = false;
	    }
	    catch (const KeyErrorException &kee) {
	       add_atom_to_plane = true;
	    }
	    // the following exception is needed for my Ubuntu 10.04 machine, 
	    // don't know why: fixes:
	    // terminate called after throwing an instance of 'KeyErrorException'
	    //   what():  std::exception
	    // 
	    catch (const std::exception &stde) {
	       add_atom_to_plane = true;
	    }

	    if (add_atom_to_plane) {
		     
	       std::string name = "";
	       at_p->getProp("name", name);
	       // add name if it is not already in the vector
	       if (std::find(plane_restraint_atoms.begin(), plane_restraint_atoms.end(), name) ==
		   plane_restraint_atoms.end())
		  plane_restraint_atoms.push_back(name);
	       at_p->setProp("plane_id", plane_id);

	       // debug
	       if (0) { 
		  std::string plane_id_lookup_debug; 
		  at_p->getProp("plane_id", plane_id_lookup_debug);
		  std::cout << "debug:: set atom " << name << " to plane_id "
			    << plane_id_lookup_debug << std::endl;
	       }

	       // run through neighours, because neighours of
	       // aromatic system atoms are in the plane too.
	       // 
	       RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
	       boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_p);
	       std::vector<const RDKit::Atom *> attached_atoms;
	       while(nbr_idx_1 != end_nbrs_1) {
		  const RDKit::Atom *at_2 = mol[*nbr_idx_1];
		  // add if not a hydrogen or we are not doing quartet hydrogen planes
		  if (at_2->getAtomicNum() != 1 || !quartet_hydrogen_planes)
		     attached_atoms.push_back(at_2);
		  ++nbr_idx_1;
	       }
	       if (attached_atoms.size() == 3) {
			
		  // Yes, there was something
		  for (unsigned int iattached=0; iattached<attached_atoms.size(); iattached++) { 
		     try {
			std::string attached_atom_name;
			attached_atoms[iattached]->getProp("name", attached_atom_name);
			// add it if it is not already in a plane,
			// 
			if (std::find(plane_restraint_atoms.begin(),
				      plane_restraint_atoms.end(),
				      attached_atom_name) ==
			    plane_restraint_atoms.end())
			   plane_restraint_atoms.push_back(attached_atom_name);
		     }
		     catch (const KeyErrorException &kee) {
			// do nothing then (no name found)
		     }
		  }
	       }
	    }
	 }
      }

      // make a plane restraint with those atoms in then
      if (plane_restraint_atoms.size() > 3) {
	 mmdb::realtype dist_esd = 0.02;
	 coot::dict_plane_restraint_t rest(plane_id, plane_restraint_atoms, dist_esd);
	 plane_restraint = rest;
      } 
   }
   
   catch (const KeyErrorException &kee) {
      // this should not happen
      std::cout << "WARNING:: add_chem_comp_planes() failed to get atom name "
		<< std::endl;
   } 

   // std::cout << "returning plane_restraint with " << plane_restraint.n_atoms() << " atoms" << std::endl;
   return plane_restraint;
}


// Return the number of added planes.
// 
// Don't add hydrogen quartets (that's done later).
// 
int
coot::add_chem_comp_aromatic_plane_quartet_planes(const RDKit::MatchVectType &match,
						  const RDKit::ROMol &mol,
						  coot::dictionary_residue_restraints_t *restraints,
						  int plane_id_idx_in) {

   std::vector<quartet_set> quartet_sets_vec;
   
   int n_planes = 0;
   try {
      for (unsigned int ii=0; ii<match.size(); ii++) {
	 const RDKit::Atom *at_p = mol[match[ii].second];
	 if (at_p->getAtomicNum() != 1) {

	    if (0) {
	       std::string name;
	       at_p->getProp("name", name);
	       std::cout << "--------- considering core atom " << match[ii].second
			 << " " << name << std::endl;
	    }
	    
	    // What are the neighbour of this atom? Are there more
	    // than 2 of them?  If so, let's make a plane restraint.

	    std::vector<unsigned int> quartet_indices;
	    quartet_indices.push_back(match[ii].second);
	    RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
	    boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_p);
	    while(nbr_idx_1 != end_nbrs_1){
	       const RDKit::Atom *at_neighb = mol[*nbr_idx_1];
	       if (at_neighb->getAtomicNum() != 1) {
		  quartet_indices.push_back(*nbr_idx_1);
	       }
	       ++nbr_idx_1;
	    }

	    
	    if (quartet_indices.size() > 3) {

	       if (0) {  // debug
		  std::cout << "debug quartet_indices.size() 3 legs path: "
			    << quartet_indices.size() << std::endl;
		  for (unsigned int jj=0; jj<quartet_indices.size(); jj++) { 
		     std::string name;
		     mol[quartet_indices[jj]]->getProp("name", name);
		     std::cout << "   " << name;
		  }
		  std::cout << std::endl;
	       }


	       quartet_set q(quartet_indices);
	       quartet_sets_vec.push_back(q);

	    } else {

	       // We need neighbours of neighbours then:
	       //
	       std::vector<unsigned int> q_indices;
	       q_indices.push_back(match[ii].second);
	       
	       RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
	       boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_p);
	       while(nbr_idx_1 != end_nbrs_1){
		  if (mol[*nbr_idx_1]->getAtomicNum() != 1) {
		     q_indices.push_back(*nbr_idx_1);
		  }
		  ++nbr_idx_1;
	       }


	       // OK quartet_indices should be 3 now.  Root atom and
	       // its two neighbours.
	       //
	       if (false) { // debug
		  std::cout << "debug quartet_indices.size() (should be 3): "
			    << q_indices.size() << std::endl;
		  for (unsigned int jj=0; jj<q_indices.size(); jj++) { 
		     std::string name;
		     mol[quartet_indices[jj]]->getProp("name", name);
		     std::cout << "   " << name;
		  }
		  std::cout << std::endl;
	       }

	       // Now neighbours, then neighbours of neighbours
	       //
	       boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_p);
	       while(nbr_idx_1 != end_nbrs_1){
		  const RDKit::Atom *at_1 = mol[*nbr_idx_1];
		  if (at_1->getAtomicNum() != 1) {

		     RDKit::ROMol::ADJ_ITER nbr_idx_2, end_nbrs_2;
		     boost::tie(nbr_idx_2, end_nbrs_2) = mol.getAtomNeighbors(at_1);
		     while(nbr_idx_2 != end_nbrs_2){

			if (mol[*nbr_idx_2]->getAtomicNum() != 1) {
                           std::vector<unsigned int> local_quartet(quartet_indices);

			   // Add this atom if it's not already in the quartet
			   if (std::find(local_quartet.begin(),
					 local_quartet.end(),
					 *nbr_idx_2) == local_quartet.end()) {
			      local_quartet.push_back(*nbr_idx_2);
			      quartet_sets_vec.push_back(quartet_set(local_quartet));
			   }
			}
			++nbr_idx_2;
		     }
		  }
		  ++nbr_idx_1;
	       }
	    } 
	 } 
      }

      // std::cout << "got quartet_sets_vec.size(): " << quartet_sets_vec.size() << std::endl;
      n_planes = quartet_sets_vec.size();

      for (unsigned int i=0; i<quartet_sets_vec.size(); i++) { 
	 const quartet_set &q = quartet_sets_vec[i];
	 std::vector<std::string> atom_names;
	 for (unsigned int iat=0; iat<4; iat++) { 
	    const RDKit::Atom *at = mol[q[iat]];
	    std::string name;
	    at->getProp("name", name);
	    atom_names.push_back(name);
	 }
	 if (atom_names.size() > 3) {
	    double esd = 0.014;
	    std::string plane_id = "quartet-plane-" + util::int_to_string(plane_id_idx_in+i);

	    if (0) { // debug
	       std::cout << "Adding plane " << plane_id << " with atoms ";
	       for (unsigned int iat=0; iat<atom_names.size(); iat++) { 
		  std::cout << atom_names[iat] << " ";
	       }
	       std::cout << std::endl;
	    } 
	    
	    dict_plane_restraint_t pr(plane_id, atom_names, esd);
	    restraints->plane_restraint.push_back(pr);
	    n_planes++;
	 }
      }
   }
   catch (const KeyErrorException &kee) {
      // this should not happen
      std::cout << "WARNING:: add_chem_comp_aromatic_plane_quartet_planes() failed to get atom name "
		<< std::endl;
   }

   return n_planes;
} 



void
coot::add_chem_comp_deloc_planes(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints) {

   typedef std::pair<std::string, double> d_pat;

   // There is a problem here - if we find the COOH from A-C(=O)O (as
   // we do below) then the refmac atom types are set correctly (good)
   // but then we fail to look up the (energy-lib) bond C-OC single
   // and C-OC double.  Hmmm....

   std::vector<d_pat> patterns;
   patterns.push_back(d_pat("*C(=O)[O;H]",                 0.02));  // ASP carboxylate, valence model (H)
   patterns.push_back(d_pat("*[C;X3;^2](~O)~[O;X1]",       0.02));  // ASP carboxylate, no H
   patterns.push_back(d_pat("AC(=O)[N^2;H2,H1]([H])[A,H]", 0.02));  // ASN
   patterns.push_back(d_pat("*C(=N)[N^2;H2]([H])[A,H]",    0.02));  // amidine
   patterns.push_back(d_pat("CNC(=[NH])N([H])[H]",         0.02));  // guanidinium with H - testing
   patterns.push_back(d_pat("CNC(=[NH])N",                 0.02));  // guanidinium sans Hs
   patterns.push_back(d_pat("*[C;X3;^2](=O)[N;X3;^2;H1]([H])*", 0.02));  // amino

   // Martin's pattern, these should be weaker (than standard 0.02) though, I think
   patterns.push_back(d_pat("[*^2]=[*^2]-[*^2]=[*;X1;^2]", 0.04));
   // patterns.push_back(d_pat("[a^2]:[a^2]-[*^2]=[*;X1;^2]", 0.04)); // no. bad match for OAC,CBG in 0BU

   int n_planes = 1; 
   for (unsigned int ipat=0; ipat<patterns.size(); ipat++) {
      RDKit::ROMol *query = RDKit::SmartsToMol(patterns[ipat].first);
      std::vector<RDKit::MatchVectType>  matches;
      bool recursionPossible = true;
      bool useChirality = true;
      bool uniquify = true;
      int matched = RDKit::SubstructMatch(mol,*query,matches,uniquify,recursionPossible, useChirality);
      // int matched = false; // 20210923-PE FIXME
      for (unsigned int imatch=0; imatch<matches.size(); imatch++) { 
	 if (matches[imatch].size() > 0) {
	    
	    if (true) { // debug
	       std::cout << "INFO:: matched deloc plane: " << patterns[ipat].first << " ";
	       std::cout << " ("; 
	       for (unsigned int iat=0; iat<matches[imatch].size(); iat++) { 
		  unsigned int atom_idx = matches[imatch][iat].second;
		  try {
		     const RDKit::Atom *at_p = mol[atom_idx];
		     std::string atom_name;
		     at_p->getProp("name", atom_name);
		     std::cout << " " << atom_name;
		  }
		  catch (const KeyErrorException &kee) {
		     std::cout << " " << atom_idx;
		  } 
	       }
	       std::cout << " )";
	       std::cout << std::endl;
	    }

	    std::vector<std::string> atom_names;
	    std::string plane_id = "plane-deloc-";
	    char s[100];
	    snprintf(s,99,"%d", n_planes);
	    plane_id += std::string(s);
	    try {
	       std::vector<std::string> atom_names;
	       for (unsigned int ii=0; ii<matches[imatch].size(); ii++) {
		  const RDKit::Atom *at_p = mol[matches[imatch][ii].second];

		  // Unlike aromatics, the atoms of this type of plane
		  // can be in more than one plane.

		  std::string name = "";
		  at_p->getProp("name", name);
		  at_p->setProp("plane_id", plane_id);
		  // std::cout << "... marking " << name << " as in " << plane_id << std::endl;
		  atom_names.push_back(name);
	       }
	       if (atom_names.size() > 3) { 
		  mmdb::realtype dist_esd = patterns[ipat].second;
		  coot::dict_plane_restraint_t res(plane_id, atom_names, dist_esd);
		  restraints->plane_restraint.push_back(res);
	       }
	    }

	    catch (const KeyErrorException &kee) {
	       std::cout << "ERROR:: in add_chem_comp_planes_deloc failed to get atom name"
			 << std::endl;
	    } 
	    n_planes++;
	 }
      }
   }
} 


void
coot::add_chem_comp_sp2_C_planes(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints) {

   typedef std::pair<std::string, double> d_pat;
   std::vector<d_pat> patterns;

   patterns.push_back(d_pat("AC(=O)C",                 0.02));  // ketone, A is usually C, maybe S.

   int n_planes = 1; // counter for output text
   for (unsigned int ipat=0; ipat<patterns.size(); ipat++) {
      RDKit::ROMol *query = RDKit::SmartsToMol(patterns[ipat].first);
      std::vector<RDKit::MatchVectType>  matches;
      bool recursionPossible = true;
      bool useChirality = true;
      bool uniquify = true;
      int matched = RDKit::SubstructMatch(mol,*query,matches,uniquify,recursionPossible, useChirality);
      // int matched = false; // 20210923-PE FIXME
      if (false) // debug
	 std::cout << "Matched " << matched << " sp2 N planes" << std::endl;
      for (unsigned int imatch=0; imatch<matches.size(); imatch++) {
	 if (matches[imatch].size() > 0) {
	    std::cout << "matched sp2 N plane pattern: " << patterns[ipat].first << std::endl;
	    std::string plane_id = "plane-sp2-N-";
	    char s[100];
	    snprintf(s,99,"%d", n_planes);
	    plane_id += std::string(s);
	    try {
	       std::vector<std::string> atom_names;
	       for (unsigned int ii=0; ii<matches[imatch].size(); ii++) {
		  const RDKit::Atom *at_p = mol[matches[imatch][ii].second];

		  // Unlike aromatics, the atoms of this type of plane
		  // can be in more than one plane.

		  std::string name = "";
		  at_p->getProp("name", name);
		  at_p->setProp("plane_id", plane_id);
		  atom_names.push_back(name);
	       }
	       if (atom_names.size() > 3) {
		  mmdb::realtype dist_esd = patterns[ipat].second;
		  coot::dict_plane_restraint_t res(plane_id, atom_names, dist_esd);
		  restraints->plane_restraint.push_back(res);
	       }
	    }
	    catch (const KeyErrorException &kee) {
	    }
	    n_planes++;
	 }
      }
   }
}

void
coot::add_chem_comp_sp2_N_planes(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints) {

   typedef std::pair<std::string, double> d_pat;
   std::vector<d_pat> patterns;
   patterns.push_back(d_pat("[c,C][N^2;H2]([H])[H]", 0.02));  // N6 on Adenosine.
                                                              // Should the Hs be replaced by *s?
   int n_planes = 1; // counter for output text
   for (unsigned int ipat=0; ipat<patterns.size(); ipat++) {
      RDKit::ROMol *query = RDKit::SmartsToMol(patterns[ipat].first);
      std::vector<RDKit::MatchVectType>  matches;
      bool recursionPossible = true;
      bool useChirality = true;
      bool uniquify = true;
      int matched = RDKit::SubstructMatch(mol,*query,matches,uniquify,recursionPossible, useChirality);
      // int matched = false; // 20210923-PE FIXME
      if (false) // debug
	 std::cout << "Matched " << matched << " sp2 N planes" << std::endl;
      for (unsigned int imatch=0; imatch<matches.size(); imatch++) { 
	 if (matches[imatch].size() > 0) {
	    std::cout << "matched sp2 N plane pattern: " << patterns[ipat].first << std::endl;
	    std::string plane_id = "plane-sp2-N-";
	    char s[100];
	    snprintf(s,99,"%d", n_planes);
	    plane_id += std::string(s);
	    try {
	       std::vector<std::string> atom_names;
	       for (unsigned int ii=0; ii<matches[imatch].size(); ii++) {
		  const RDKit::Atom *at_p = mol[matches[imatch][ii].second];

		  // Unlike aromatics, the atoms of this type of plane
		  // can be in more than one plane.

		  std::string name = "";
		  at_p->getProp("name", name);
		  at_p->setProp("plane_id", plane_id);
		  atom_names.push_back(name);
	       }
	       if (atom_names.size() > 3) { 
		  mmdb::realtype dist_esd = patterns[ipat].second;
		  coot::dict_plane_restraint_t res(plane_id, atom_names, dist_esd);
		  restraints->plane_restraint.push_back(res);
	       }
	    }
	    catch (const KeyErrorException &kee) {
		  
	    }
	    n_planes++;
	 }
      }
   }
}

// alter restraints.
int 
coot::assign_chirals(const RDKit::ROMol &mol, coot::dictionary_residue_restraints_t *restraints) {

   int n_chirals = assign_chirals_mmcif_tags(mol, restraints);
   if (! n_chirals) {
      n_chirals = assign_chirals_rdkit_tags(mol, restraints);
   }
   return n_chirals;
}

   
// alter restraints.
int 
coot::assign_chirals_mmcif_tags(const RDKit::ROMol &mol,
				coot::dictionary_residue_restraints_t *restraints) {

   int n_chirals = 0;

   unsigned int n_atoms = mol.getNumAtoms();
   for (unsigned int iat=0; iat<n_atoms; iat++) { 
      const RDKit::Atom *at_p = mol[iat];
      try {
	 std::string ch;
	 std::string chiral_centre, n1, n2, n3;
	 at_p->getProp("name", chiral_centre);
	 at_p->getProp("mmcif_chiral_volume_sign", ch);
	 at_p->getProp("mmcif_chiral_N1", n1);
	 at_p->getProp("mmcif_chiral_N2", n2);
	 at_p->getProp("mmcif_chiral_N3", n3);
	 int cv = dict_chiral_restraint_t::CHIRAL_VOLUME_RESTRAINT_VOLUME_SIGN_UNASSIGNED;
	 if (ch == "positive")
	    cv = dict_chiral_restraint_t::CHIRAL_RESTRAINT_POSITIVE;
	 if (ch == "negative")
	    cv = dict_chiral_restraint_t::CHIRAL_RESTRAINT_NEGATIVE;
	 if (ch == "both")
	    cv = dict_chiral_restraint_t::CHIRAL_RESTRAINT_BOTH;
	 
	 std::string chiral_id = "chiral_" + util::int_to_string(n_chirals+1);
	 dict_chiral_restraint_t cr(chiral_id, chiral_centre, n1, n2, n3, cv);
	 restraints->chiral_restraint.push_back(cr);
	 n_chirals++;
      }
      catch (const KeyErrorException &kee) {
	 // it's OK for this to happen.  Most atoms are not chiral.
	 // Most input molecules are not from mmcif files.
	 // 
	 // std::cout << "oops caught exception " << kee.key() << " "
	 // << kee.what() << std::endl;
      }
      // This shouldn't be needed because the exception should only be
      // KeyErrorException, and that should be caught by the block above.
      catch (const std::runtime_error &rte) {

	 if (0)
	    std::cout << "assign_chirals_mmcif_tags() strange - caught runtime_error "
		      << rte.what() << std::endl;
      }
   }
   return n_chirals;
}


// return -1 or +1
int
coot::get_volume_sign_from_coordinates(const RDKit::ROMol &mol,
                                       unsigned int idx_chiral_centre_atom,
                                       const std::vector<indexed_name_and_rank_t> &neighb_names_and_ranks) {

   auto make_vector = [] (const RDGeom::Point3D &a, const RDGeom::Point3D &central) {
                         return clipper::Coord_orth(a.x-central.x, a.y-central.y, a.z-central.z);
                      };

   int n_conf = mol.getNumConformers();
   if (n_conf > 0) {
      int id_conf = n_conf -1;
      const RDKit::Conformer &conf = mol.getConformer(id_conf);
      const RDGeom::Point3D &pos_central = conf.getAtomPos(idx_chiral_centre_atom);
      const RDGeom::Point3D &pos_a = conf.getAtomPos(neighb_names_and_ranks[0].atom_index);
      const RDGeom::Point3D &pos_b = conf.getAtomPos(neighb_names_and_ranks[1].atom_index);
      const RDGeom::Point3D &pos_c = conf.getAtomPos(neighb_names_and_ranks[2].atom_index);
      clipper::Coord_orth a = make_vector(pos_a, pos_central);
      clipper::Coord_orth b = make_vector(pos_b, pos_central);
      clipper::Coord_orth c = make_vector(pos_c, pos_central);
      double vol = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
      // std::cout << "chiral vol " << vol << std::endl;
      if (vol > 0)
         return 1;
      else
         return -1;
   }
   return 1;
}

// alter restraints: RDKit/SMILES -> mmCIF chiral conversion
int
coot::assign_chirals_rdkit_tags(const RDKit::ROMol &mol,
				coot::dictionary_residue_restraints_t *restraints) {

   // debug_cip_ranks(mol);

   int n_chirals = 0;

   unsigned int n_atoms = mol.getNumAtoms();
   for (unsigned int iat=0; iat<n_atoms; iat++) {
      int vol_sign = coot::dict_chiral_restraint_t::CHIRAL_VOLUME_RESTRAINT_VOLUME_SIGN_UNASSIGNED;
      const RDKit::Atom *at_p = mol[iat];
      RDKit::Atom::ChiralType chiral_tag = at_p->getChiralTag();
      if (false)
	 std::cout << "DEBUG:: in assign_chirals_rdkit_tags() atom " << iat
		   << " chiral tag: " << chiral_tag << std::endl;

      // I think that these are round the wrong way.
      if (chiral_tag == RDKit::Atom::CHI_TETRAHEDRAL_CW)
	 vol_sign = dict_chiral_restraint_t::CHIRAL_RESTRAINT_NEGATIVE;
      if (chiral_tag == RDKit::Atom::CHI_TETRAHEDRAL_CCW)
	 vol_sign = dict_chiral_restraint_t::CHIRAL_RESTRAINT_POSITIVE;

      if (chiral_tag != RDKit::Atom::CHI_UNSPECIFIED) {
	 try {

	    if (false)
	       std::cout << "DEBUG:: in assign_chirals_rdkit_tags(): considering chiral "
			 << "for atom idx " << iat << std::endl;

            unsigned int idx_central = iat;
	    std::string chiral_centre;
	    at_p->getProp("name", chiral_centre);
	    std::string n1, n2, n3; // these need setting, c.f.
	                            // get_chiral_tag() in rdkit-interface.cc?

	    // What are the neighbours of at_p and what are their ranks?
	    //

	    // std::vector<std::pair<unsigned int, std::string> > neighb_names_and_ranks;
            std::vector<indexed_name_and_rank_t> neighb_names_and_ranks;
	    RDKit::ROMol::ADJ_ITER nbr_idx_1, end_nbrs_1;
	    boost::tie(nbr_idx_1, end_nbrs_1) = mol.getAtomNeighbors(at_p);
	    while(nbr_idx_1 != end_nbrs_1){
	       const RDKit::Atom *at_neighb = mol[*nbr_idx_1];
	       unsigned int cip_rank;
	       std::string neighb_name;
	       at_neighb->getProp(RDKit::common_properties::_CIPRank,cip_rank);
	       at_neighb->getProp("name", neighb_name);
	       indexed_name_and_rank_t p(*nbr_idx_1, cip_rank, neighb_name);
	       neighb_names_and_ranks.push_back(p);
	       ++nbr_idx_1;
	    }

	    if (false) { // debug
	       std::sort(neighb_names_and_ranks.begin(),
			 neighb_names_and_ranks.end());

	       std::reverse(neighb_names_and_ranks.begin(),
			    neighb_names_and_ranks.end());

	       std::cout << "DEBUG:: in assign_chirals_rdkit_tags() for atom "
			 << chiral_centre << " found "
			 << neighb_names_and_ranks.size() << "  neighbours: ";
	       for (unsigned int ii=0; ii<neighb_names_and_ranks.size(); ii++)
		  std::cout << " "
			    << coot::util::remove_whitespace(neighb_names_and_ranks[ii].atom_name)
                            << " (rank " << neighb_names_and_ranks[ii].cip_rank << ")";
	       std::cout << std::endl;
	    }

	    if (neighb_names_and_ranks.size() == 4) {

	       //std::sort(neighb_names_and_ranks.begin(),
               // neighb_names_and_ranks.end());

	       // std::reverse(neighb_names_and_ranks.begin(),
               // neighb_names_and_ranks.end());

	       std::string chiral_id = "chiral_" + util::int_to_string(n_chirals+1);
	       std::string n1 = neighb_names_and_ranks[0].atom_name;
	       std::string n2 = neighb_names_and_ranks[1].atom_name;
	       std::string n3 = neighb_names_and_ranks[2].atom_name;

               vol_sign = get_volume_sign_from_coordinates(mol, idx_central, neighb_names_and_ranks);

	       dict_chiral_restraint_t cr(chiral_id, chiral_centre, n1, n2, n3, vol_sign);
	       restraints->chiral_restraint.push_back(cr);

	       n_chirals++;
	    }
	 }

	 catch (KeyErrorException &kee) {
	    std::cout << "assign_chirals_rdkit_tags(): no prop name " << iat << std::endl;
	 }
      }
   }
   return n_chirals;
}

void
coot::debug_cip_ranks(const RDKit::ROMol &mol) {

   unsigned int n_atoms = mol.getNumAtoms();
   for (unsigned int iat=0; iat<n_atoms; iat++) {
      const RDKit::Atom *at_p = mol[iat];
      try {
	 unsigned int cip_rank;
	 at_p->getProp(RDKit::common_properties::_CIPRank,cip_rank);
	 std::cout << "DEBUG:: debug_cip_ranks() " << iat << " " << cip_rank << std::endl;
      }
      catch (const KeyErrorException &kee) {
	    std::cout << "caught no-cip_rank for atom exception in debug_cip_ranks(): "
		      <<  kee.what() << std::endl;
      }
   }

}



#ifdef LIBCOOTAPI_BUILD
#else
void
coot::write_restraints(PyObject *restraints_py,
		       const std::string &monomer_type,
		       const std::string &file_name) {

   coot::dictionary_residue_restraints_t rest = monomer_restraints_from_python(restraints_py);
   rest.set_use_nuclear_distances(true);
   rest.write_cif(file_name);
}
#endif


#ifdef LIBCOOTAPI_BUILD
#else
void
coot::write_pdb_from_mol(PyObject *rdkit_mol_py,
			 const std::string &res_name,
			 const std::string &file_name) {

   RDKit::ROMol &mol = boost::python::extract<RDKit::ROMol&>(rdkit_mol_py);
   mmdb::Residue *res = coot::make_residue(mol, 0, res_name);
   if (! res) {
      std::cout << "in write_pdb_from_mol() failed to make residue" << std::endl;
   } else {
      mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(res);
      int status = mol->WritePDBASCII(file_name.c_str());
      if (status == 0)
	 std::cout << "INFO:: wrote PDB   \"" << file_name << "\"" << std::endl;
      delete mol;
   }
}
#endif



#ifdef LIBCOOTAPI_BUILD
#else
void
coot::regularize_and_write_pdb(PyObject *rdkit_mol, PyObject *restraints_py,
			       const std::string &res_name,
			       const std::string &pdb_file_name) {

   std::pair<mmdb::Manager *, mmdb::Residue *> mol_res = regularize_inner(rdkit_mol, restraints_py, res_name);
   int status = mol_res.first->WritePDBASCII(pdb_file_name.c_str());
   if (status == 0)
      std::cout << "INFO:: wrote PDB   \"" << pdb_file_name << "\"" << std::endl;
}
#endif


#ifdef LIBCOOTAPI_BUILD
#else
// update the passed rdkit molecule
void
coot::regularize(PyObject *rdkit_mol_py, PyObject *restraints_py,
			   const std::string &res_name) {

   RDKit::ROMol &mol = boost::python::extract<RDKit::ROMol&>(rdkit_mol_py);
   
   std::pair<mmdb::Manager *, mmdb::Residue *> regular =
      regularize_inner(rdkit_mol_py, restraints_py, res_name);

   if (regular.second) {

      // now create a new molecule, because the one we are given is a ROMol.
      RDKit::RWMol *rw_mol = new RDKit::RWMol(mol);

      int iconf = 0;
      // this shouldn't move the atoms if bypass_refinement is true in regularize_inner().
      update_coords(rw_mol, iconf, regular.second);
   }
}
#endif


// update mol and restraints_p
void
coot::regularize_and_update_mol_and_restraints(RDKit::RWMol *mol,
                                               coot::dictionary_residue_restraints_t *restraints_p) {

   int n_conf = mol->getNumConformers();
   if (n_conf > 0) {
      int i_conf = n_conf -1;
      std::string res_name = restraints_p->residue_info.comp_id;
      mmdb::Residue *residue_p = coot::make_residue(*mol, i_conf, res_name);
      mmdb::Manager *mmdb_mol = coot::util::create_mmdbmanager_from_residue(residue_p);
      mmdb::Residue *first_residue_p = coot::util::get_first_residue(mmdb_mol);

      simple_refine(first_residue_p, mmdb_mol, *restraints_p);
      update_coords(mol, i_conf, first_residue_p);
      update_chem_comp_atoms_from_residue(first_residue_p, restraints_p);
      delete mmdb_mol;
      delete residue_p;
   } else {
      std::cout << "WARNING:: regularize_and_update_mol_and_restraints() no conformers means no minimization"
                << std::endl;
   }
}


#ifdef LIBCOOTAPI_BUILD
#else
std::pair<mmdb::Manager *, mmdb::Residue *>
coot::regularize_inner(PyObject *rdkit_mol_py,
		       PyObject *restraints_py,
		       const std::string &res_name) {

   RDKit::ROMol &mol = boost::python::extract<RDKit::ROMol&>(rdkit_mol_py);
   return regularize_inner(mol, restraints_py, res_name);
}
#endif


#ifdef LIBCOOTAPI_BUILD
#else
std::pair<mmdb::Manager *, mmdb::Residue *>
coot::regularize_inner(RDKit::ROMol &mol,
		       PyObject *restraints_py,
		       const std::string &res_name) {

   coot::dictionary_residue_restraints_t dict_restraints =
      monomer_restraints_from_python(restraints_py);
   mmdb::Residue *residue_p = coot::make_residue(mol, 0, res_name);
   // remove this NULL at some stage (soon)
   mmdb::Manager *cmmdbmanager = coot::util::create_mmdbmanager_from_residue(residue_p);

   simple_refine(residue_p, cmmdbmanager, dict_restraints);

   return std::pair<mmdb::Manager *, mmdb::Residue *> (cmmdbmanager, residue_p);
}
#endif



// PyObject *
// coot::convert_rdkit_mol_to_pyobject(RDKit::ROMol *mol) {

//    PyObject *r = Py_False;

//    return r;
// }
