/* pyrogen/py-restraints.cc
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
#include <Python.h>
#include <boost/python.hpp>
#endif

#include <cstring> // Fix strchr problems on using RDKit includes.

#include "compat/coot-sysdep.h"
#include <lidia-core/use-rdkit.hh>
#include "restraints.hh"
#ifdef LIBCOOTAPI_BUILD
#else
#include "py-restraints.hh"
#endif
#include <lidia-core/rdkit-interface.hh>
#include <coot-utils/coot-coord-utils.hh>

std::string
coot::convert_to_energy_lib_bond_type(RDKit::Bond::BondType bt) {

   std::string r = "unset";

   if (bt == RDKit::Bond::UNSPECIFIED) r = "unset";
   if (bt == RDKit::Bond::SINGLE)      r = "single";
   if (bt == RDKit::Bond::DOUBLE)      r = "double";
   if (bt == RDKit::Bond::TRIPLE)      r = "triple";
   if (bt == RDKit::Bond::QUADRUPLE)   r = "quadruple";
   if (bt == RDKit::Bond::QUINTUPLE)   r = "quintuple";
   if (bt == RDKit::Bond::HEXTUPLE)    r = "hextuple";
   if (bt == RDKit::Bond::ONEANDAHALF) r = "deloc";
   if (bt == RDKit::Bond::AROMATIC)    r = "aromatic";
   
//       TWOANDAHALF,
//       THREEANDAHALF,
//       FOURANDAHALF,
//       FIVEANDAHALF,
//       AROMATIC,
//       IONIC,
//       HYDROGEN,

   return r;
} 


#ifdef LIBCOOTAPI_BUILD
#else
coot::dictionary_residue_restraints_t
monomer_restraints_from_python(PyObject *restraints) {

   coot::dictionary_residue_restraints_t rest;

   std::vector<coot::dict_bond_restraint_t> bond_restraints;
   std::vector<coot::dict_angle_restraint_t> angle_restraints;
   std::vector<coot::dict_torsion_restraint_t> torsion_restraints;
   std::vector<coot::dict_chiral_restraint_t> chiral_restraints;
   std::vector<coot::dict_plane_restraint_t> plane_restraints;
   std::vector<coot::dict_atom> atoms;
   coot::dict_chem_comp_t residue_info;

   PyObject *key;
   PyObject *value;
   Py_ssize_t pos = 0;

   while (PyDict_Next(restraints, &pos, &key, &value)) {
      // std::cout << ":::::::key: " << PyString_AsString(key) << std::endl;

      std::string key_string = PyBytes_AS_STRING(PyUnicode_AsUTF8String(key));

      if (key_string == "_chem_comp") {
	 PyObject *chem_comp_list = value;
	 if (PyList_Check(chem_comp_list)) {
	    if (PyObject_Length(chem_comp_list) == 7) {
	       std::string comp_id  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 0)));
	       std::string tlc      = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 1)));
	       std::string name     = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 2)));
	       std::string group    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 3)));
	       long n_atoms_all      = PyLong_AsLong(PyList_GetItem(chem_comp_list, 4));
	       long n_atoms_nh       = PyLong_AsLong(PyList_GetItem(chem_comp_list, 5));
	       std::string desc_lev = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 6)));

	       coot::dict_chem_comp_t n(comp_id, tlc, name, group,
					n_atoms_all, n_atoms_nh, desc_lev);
	       residue_info = n;
	    }
	 }
      }

      if (key_string == "_chem_comp_atom") {
	 PyObject *chem_comp_atom_list = value;
	 if (PyList_Check(chem_comp_atom_list)) {
	    int n_atoms = PyObject_Length(chem_comp_atom_list);
	    for (int iat=0; iat<n_atoms; iat++) {
	       PyObject *chem_comp_atom = PyList_GetItem(chem_comp_atom_list, iat);
	       bool have_coords = false;
	       bool OK_length = false;
	       if ((PyObject_Length(chem_comp_atom) == 5)) {
		  OK_length = true;
		  have_coords = false;
	       }
	       if ((PyObject_Length(chem_comp_atom) == 6)) {
		  OK_length = true;
		  have_coords = true;
	       }
	       if (OK_length) {
		  std::string atom_id  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_atom, 0)));
		  std::string element  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_atom, 1)));
		  std::string energy_t = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_atom, 2)));
		  float part_chr        = PyFloat_AsDouble(PyList_GetItem(chem_comp_atom, 3));
		  bool flag = 0;
		  if (PyLong_AsLong(PyList_GetItem(chem_comp_atom, 4))) {
		     flag = 1;
		  }
		  std::pair<bool, float> part_charge_info(flag, part_chr);
		  coot::dict_atom at(atom_id, atom_id, element, energy_t, part_charge_info);
		  if (have_coords) {
		     PyObject *coords = PyList_GetItem(chem_comp_atom, 5);
		     if (PyList_Check(coords)) {
			int n_coords = PyObject_Length(coords);
			if (n_coords == 3) {
			   PyObject *x_py = PyList_GetItem(coords, 0);
			   PyObject *y_py = PyList_GetItem(coords, 1);
			   PyObject *z_py = PyList_GetItem(coords, 2);
			   float x = PyFloat_AsDouble(x_py);
			   float y = PyFloat_AsDouble(y_py);
			   float z = PyFloat_AsDouble(z_py);
			   clipper::Coord_orth co(x,y,z);
			   at.model_Cartn = std::pair<bool, clipper::Coord_orth> (true, co);
			}
		     }
		  }
		  atoms.push_back(at);
	       }
	    }
	 }
      }

      if (key_string == "_chem_comp_bond") {
	 PyObject *bond_restraint_list = value;
	 if (PyList_Check(bond_restraint_list)) {
	    int n_bonds = PyObject_Length(bond_restraint_list);
	    for (int i_bond=0; i_bond<n_bonds; i_bond++) {
	       PyObject *bond_restraint = PyList_GetItem(bond_restraint_list, i_bond);
	       if (PyObject_Length(bond_restraint) == 5) {
		  PyObject *atom_1_py = PyList_GetItem(bond_restraint, 0);
		  PyObject *atom_2_py = PyList_GetItem(bond_restraint, 1);
		  PyObject *type_py   = PyList_GetItem(bond_restraint, 2);
		  PyObject *dist_py   = PyList_GetItem(bond_restraint, 3);
		  PyObject *esd_py    = PyList_GetItem(bond_restraint, 4);

		  if (PyUnicode_Check(atom_1_py) &&
		      PyUnicode_Check(atom_2_py) &&
		      PyUnicode_Check(type_py) &&
		      PyFloat_Check(dist_py) &&
		      PyFloat_Check(esd_py)) {
		     std::string atom_1 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_1_py));
		     std::string atom_2 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_2_py));
		     std::string type   = PyBytes_AS_STRING(PyUnicode_AsUTF8String(type_py));
		     float  dist = PyFloat_AsDouble(dist_py);
		     float  esd  = PyFloat_AsDouble(esd_py);
		     coot::dict_bond_restraint_t rest(atom_1, atom_2, type, dist, esd, 0.0, 0.0, false);
		     bond_restraints.push_back(rest);
		  } else {
                     // ERROR?
                     std::cout << "WARNING:: Bad Pythonic type match in bond restraint" << std::endl;
                  }
	       }
	    }
	 }
      }


      if (key_string == "_chem_comp_angle") {
	 PyObject *angle_restraint_list = value;
	 if (PyList_Check(angle_restraint_list)) {
	    int n_angles = PyObject_Length(angle_restraint_list);
	    for (int i_angle=0; i_angle<n_angles; i_angle++) {
	       PyObject *angle_restraint = PyList_GetItem(angle_restraint_list, i_angle);
	       if (PyObject_Length(angle_restraint) == 5) {
		  PyObject *atom_1_py = PyList_GetItem(angle_restraint, 0);
		  PyObject *atom_2_py = PyList_GetItem(angle_restraint, 1);
		  PyObject *atom_3_py = PyList_GetItem(angle_restraint, 2);
		  PyObject *angle_py  = PyList_GetItem(angle_restraint, 3);
		  PyObject *esd_py    = PyList_GetItem(angle_restraint, 4);

		  if (PyUnicode_Check(atom_1_py) &&
		      PyUnicode_Check(atom_2_py) &&
		      PyUnicode_Check(atom_3_py) &&
		      PyFloat_Check(angle_py) && 
		      PyFloat_Check(esd_py)) {
		     std::string atom_1 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_1_py));
		     std::string atom_2 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_2_py));
		     std::string atom_3 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_3_py));
		     float  angle = PyFloat_AsDouble(angle_py);
		     float  esd   = PyFloat_AsDouble(esd_py);
		     coot::dict_angle_restraint_t rest(atom_1, atom_2, atom_3, angle, esd);
		     angle_restraints.push_back(rest);
		  }
	       }
	    }
	 }
      }

      if (key_string == "_chem_comp_tor") {
	 PyObject *torsion_restraint_list = value;
	 if (PyList_Check(torsion_restraint_list)) {
	    int n_torsions = PyObject_Length(torsion_restraint_list);
	    for (int i_torsion=0; i_torsion<n_torsions; i_torsion++) {
	       PyObject *torsion_restraint = PyList_GetItem(torsion_restraint_list, i_torsion);
	       if (PyObject_Length(torsion_restraint) == 7) { // info for Nigel.
		  std::cout << "torsions now have 8 elements starting with the torsion id\n"; 
	       } 
	       if (PyObject_Length(torsion_restraint) == 8) { 
		  PyObject *id_py     = PyList_GetItem(torsion_restraint, 0);
		  PyObject *atom_1_py = PyList_GetItem(torsion_restraint, 1);
		  PyObject *atom_2_py = PyList_GetItem(torsion_restraint, 2);
		  PyObject *atom_3_py = PyList_GetItem(torsion_restraint, 3);
		  PyObject *atom_4_py = PyList_GetItem(torsion_restraint, 4);
		  PyObject *torsion_py= PyList_GetItem(torsion_restraint, 5);
		  PyObject *esd_py    = PyList_GetItem(torsion_restraint, 6);
		  PyObject *period_py = PyList_GetItem(torsion_restraint, 7);

		  if (PyUnicode_Check(atom_1_py) &&
		      PyUnicode_Check(atom_2_py) &&
		      PyUnicode_Check(atom_3_py) &&
		      PyUnicode_Check(atom_4_py) &&
		      PyFloat_Check(torsion_py)  &&
		      PyFloat_Check(esd_py)      &&
		      PyLong_Check(period_py)) {
		     std::string id     = PyBytes_AS_STRING(PyUnicode_AsUTF8String(id_py));
		     std::string atom_1 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_1_py));
		     std::string atom_2 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_2_py));
		     std::string atom_3 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_3_py));
		     std::string atom_4 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_4_py));
		     float torsion = PyFloat_AsDouble(torsion_py);
		     float esd     = PyFloat_AsDouble(esd_py);
		     long  period  = PyLong_AsLong(period_py);
		     coot::dict_torsion_restraint_t rest(id, atom_1, atom_2, atom_3, atom_4,
							 torsion, esd, period);
		     torsion_restraints.push_back(rest);
		  }
	       }
	    }
	 }
      }

      if (key_string == "_chem_comp_chir") {
	 PyObject *chiral_restraint_list = value;
	 if (PyList_Check(chiral_restraint_list)) {
	    int n_chirals = PyObject_Length(chiral_restraint_list);
	    for (int i_chiral=0; i_chiral<n_chirals; i_chiral++) {
	       PyObject *chiral_restraint = PyList_GetItem(chiral_restraint_list, i_chiral);
	       if (PyObject_Length(chiral_restraint) == 7) { 
		  PyObject *chiral_id_py= PyList_GetItem(chiral_restraint, 0);
		  PyObject *atom_c_py   = PyList_GetItem(chiral_restraint, 1);
		  PyObject *atom_1_py   = PyList_GetItem(chiral_restraint, 2);
		  PyObject *atom_2_py   = PyList_GetItem(chiral_restraint, 3);
		  PyObject *atom_3_py   = PyList_GetItem(chiral_restraint, 4);
		  PyObject *vol_sign_py = PyList_GetItem(chiral_restraint, 5);
		  PyObject *esd_py      = PyList_GetItem(chiral_restraint, 6);

		  if (PyUnicode_Check(atom_1_py) &&
		      PyUnicode_Check(atom_2_py) &&
		      PyUnicode_Check(atom_3_py) &&
		      PyUnicode_Check(atom_c_py) &&
		      PyUnicode_Check(chiral_id_py) && 
		      PyFloat_Check(esd_py)    && 
		      PyLong_Check(vol_sign_py)) {
		     std::string chiral_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(chiral_id_py));
		     std::string atom_c    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_c_py));
		     std::string atom_1    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_1_py));
		     std::string atom_2    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_2_py));
		     std::string atom_3    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_3_py));
		     float  esd            = PyFloat_AsDouble(esd_py);
		     long vol_sign         = PyLong_AsLong(vol_sign_py);
		     coot::dict_chiral_restraint_t rest(chiral_id,
							atom_c, atom_1, atom_2, atom_3, 
							vol_sign);
		     chiral_restraints.push_back(rest);
		  }
	       }
	    }
	 }
      }


      if (key_string == "_chem_comp_plane_atom") {
	 PyObject *plane_restraint_list = value;
	 if (PyList_Check(plane_restraint_list)) {
	    int n_planes = PyObject_Length(plane_restraint_list);
	    for (int i_plane=0; i_plane<n_planes; i_plane++) {
	       PyObject *plane_restraint = PyList_GetItem(plane_restraint_list, i_plane);

	       // per-atom restraints: a plane-name and list of [atom-name, esd]s.
	       // 
	       if (PyObject_Length(plane_restraint) == 2) {
		  PyObject *plane_id_py =        PyList_GetItem(plane_restraint, 0);
		  PyObject *atom_esd_pair_list = PyList_GetItem(plane_restraint, 1);
		  if (PyUnicode_Check(plane_id_py)) {
		     if (PyList_Check(atom_esd_pair_list)) { 
			std::string plane_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(plane_id_py));
			int n_atoms = PyObject_Length(atom_esd_pair_list);
			std::vector<std::pair<std::string, double> > atom_name_esd_pair_vec;
			for (int iat=0; iat<n_atoms; iat++) {
			   PyObject *nep = PyList_GetItem(atom_esd_pair_list, iat);
			   if (PyList_Check(nep)) { 
			      PyObject *atom_name_py = PyList_GetItem(nep, 0);
			      PyObject *esd_py       = PyList_GetItem(nep, 1);
			      if (PyUnicode_Check(atom_name_py)) {
				 if (PyFloat_Check(esd_py)) {
				    std::string atom_name = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_name_py));
				    float esd = PyFloat_AsDouble(esd_py);
				    std::pair<std::string, double> p(atom_name, esd);
				    atom_name_esd_pair_vec.push_back(p);
				 }
			      }
			   }
			}
			if (atom_name_esd_pair_vec.size() > 3) {
			   coot::dict_plane_restraint_t rest(plane_id, atom_name_esd_pair_vec);
			   plane_restraints.push_back(rest);
			} 
		     }
		  }
	       } 

	       // old style: [plane-id, atom-list, esd] 
	       if (PyObject_Length(plane_restraint) == 3) {
		  std::vector<std::string> atoms;
		  PyObject *plane_id_py = PyList_GetItem(plane_restraint, 0);
		  PyObject *esd_py      = PyList_GetItem(plane_restraint, 2);
		  PyObject *py_atoms_py = PyList_GetItem(plane_restraint, 1);

		  bool atoms_pass = 1;
		  if (PyList_Check(py_atoms_py)) {
		     int n_atoms = PyObject_Length(py_atoms_py);
		     for (int iat=0; iat<n_atoms; iat++) {
			PyObject *at_py = PyList_GetItem(py_atoms_py, iat);
			if (PyUnicode_Check(at_py)) {
			   atoms.push_back(PyBytes_AS_STRING(PyUnicode_AsUTF8String(at_py)));
			} else {
			   atoms_pass = 0;
			}
		     }
		     if (atoms_pass) {
			if (PyUnicode_Check(plane_id_py)) {
			   if (PyFloat_Check(esd_py)) {
			      std::string plane_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(plane_id_py));
			      float esd = PyFloat_AsDouble(esd_py);
			      if (atoms.size() > 0) { 
				 coot::dict_plane_restraint_t rest(plane_id, atoms, esd);
				 plane_restraints.push_back(rest);
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 } else {
            std::cout << "plane_restraint_list was not a list" << std::endl;
         }
      }
   }
	
   coot::dictionary_residue_restraints_t monomer_restraints;
   monomer_restraints.bond_restraint    = bond_restraints;
   monomer_restraints.angle_restraint   = angle_restraints;
   monomer_restraints.torsion_restraint = torsion_restraints;
   monomer_restraints.chiral_restraint  = chiral_restraints;
   monomer_restraints.plane_restraint   = plane_restraints;
   monomer_restraints.residue_info      = residue_info;
   monomer_restraints.atom_info         = atoms;

   if (false)
      std::cout << "debug:: monomer_restraints_from_python() returning monomer_restraints "
                << "with bond_restraints size " << monomer_restraints.bond_restraint.size()
                << std::endl;

   return monomer_restraints;
} 
#endif

#ifdef LIBCOOTAPI_BUILD
#else
PyObject *coot::monomer_restraints_to_python(const dictionary_residue_restraints_t &restraints) {

   PyObject *r = Py_False;

   if (1) { 
      r = PyDict_New();
	 
      // ------------------ chem_comp -------------------------
      const coot::dict_chem_comp_t &info = restraints.residue_info;
      std::vector<dict_atom> atom_info = restraints.atom_info;
      unsigned int n_atoms = atom_info.size();
      
      PyObject *chem_comp_py = PyList_New(7);
      PyList_SetItem(chem_comp_py, 0, PyUnicode_FromString(info.comp_id.c_str()));
      PyList_SetItem(chem_comp_py, 1, PyUnicode_FromString(info.three_letter_code.c_str()));
      PyList_SetItem(chem_comp_py, 2, PyUnicode_FromString(info.name.c_str()));
      PyList_SetItem(chem_comp_py, 3, PyUnicode_FromString(info.group.c_str()));
      PyList_SetItem(chem_comp_py, 4, PyLong_FromLong(info.number_atoms_all));
      PyList_SetItem(chem_comp_py, 5, PyLong_FromLong(info.number_atoms_nh));
      PyList_SetItem(chem_comp_py, 6, PyUnicode_FromString(info.description_level.c_str()));

      // Put chem_comp_py into a dictionary?
      PyDict_SetItem(r, PyUnicode_FromString("_chem_comp"), chem_comp_py);


      // ------------------ chem_comp_atom -------------------------
      bool have_atom_coords = false;
      unsigned int n_with_coords = 0;
      for (unsigned int iat=0; iat<n_atoms; iat++) {
	 const dict_atom &atom = atom_info[iat];
	 if (atom.model_Cartn.first)
	    n_with_coords++;
      }
      if (n_with_coords == n_atoms)
	 have_atom_coords = true;

      PyObject *atom_info_list = PyList_New(n_atoms);
      for (unsigned int iat=0; iat<n_atoms; iat++) {
	 PyObject *atom_attributes_list = NULL;
	 if (! have_atom_coords) { 
	    atom_attributes_list = PyList_New(5);
	 } else {
	    atom_attributes_list = PyList_New(6); // additional coordinates list
	 }

	 const dict_atom &atom = atom_info[iat];
	 PyList_SetItem(atom_attributes_list, 0, PyUnicode_FromString(atom.atom_id_4c.c_str()));
	 PyList_SetItem(atom_attributes_list, 1, PyUnicode_FromString(atom.type_symbol.c_str()));
	 PyList_SetItem(atom_attributes_list, 2, PyUnicode_FromString(atom.type_energy.c_str()));
	 PyList_SetItem(atom_attributes_list, 3, PyFloat_FromDouble(atom.partial_charge.second));
	 PyObject *flag = Py_False;
	 if (atom_info[iat].partial_charge.first)
	    flag = Py_True;
	 Py_INCREF(flag);
	 PyList_SetItem(atom_attributes_list, 4, flag);
	 if (have_atom_coords) {
	    PyObject *coords = PyList_New(3);
	    PyList_SetItem(coords, 0, PyFloat_FromDouble(atom.model_Cartn.second.x()));
	    PyList_SetItem(coords, 1, PyFloat_FromDouble(atom.model_Cartn.second.y()));
	    PyList_SetItem(coords, 2, PyFloat_FromDouble(atom.model_Cartn.second.z()));
	    PyList_SetItem(atom_attributes_list, 5, coords);
	 }
	 PyList_SetItem(atom_info_list, iat, atom_attributes_list);
      }

      PyDict_SetItem(r, PyUnicode_FromString("_chem_comp_atom"), atom_info_list);

      // ------------------ Bonds -------------------------
      PyObject *bond_restraint_list = PyList_New(restraints.bond_restraint.size());
      for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
	 std::string a1   = restraints.bond_restraint[ibond].atom_id_1_4c();
	 std::string a2   = restraints.bond_restraint[ibond].atom_id_2_4c();
	 std::string type = restraints.bond_restraint[ibond].type();

	 PyObject *py_value_dist = Py_False;
	 PyObject *py_value_esd = Py_False;

	 try { 
	    double d   = restraints.bond_restraint[ibond].value_dist();
	    double esd = restraints.bond_restraint[ibond].value_esd();
	    py_value_dist = PyFloat_FromDouble(d);
	    py_value_esd  = PyFloat_FromDouble(esd);
	 }
	 catch (const std::runtime_error &rte) {

	    // Use default false values.
	    // So I suppose that I need to do this then:
	    if (PyBool_Check(py_value_dist))
	       Py_INCREF(py_value_dist);
	    if (PyBool_Check(py_value_esd))
	       Py_INCREF(py_value_esd);
	 }
	 
	 PyObject *bond_restraint = PyList_New(5);
	 PyList_SetItem(bond_restraint, 0, PyUnicode_FromString(a1.c_str()));
	 PyList_SetItem(bond_restraint, 1, PyUnicode_FromString(a2.c_str()));
	 PyList_SetItem(bond_restraint, 2, PyUnicode_FromString(type.c_str()));
	 PyList_SetItem(bond_restraint, 3, py_value_dist);
	 PyList_SetItem(bond_restraint, 4, py_value_esd);
	 PyList_SetItem(bond_restraint_list, ibond, bond_restraint);
      }

      PyDict_SetItem(r, PyUnicode_FromString("_chem_comp_bond"), bond_restraint_list);


      // ------------------ Angles -------------------------
      PyObject *angle_restraint_list = PyList_New(restraints.angle_restraint.size());
      for (unsigned int iangle=0; iangle<restraints.angle_restraint.size(); iangle++) {
	 std::string a1 = restraints.angle_restraint[iangle].atom_id_1_4c();
	 std::string a2 = restraints.angle_restraint[iangle].atom_id_2_4c();
	 std::string a3 = restraints.angle_restraint[iangle].atom_id_3_4c();
	 double d   = restraints.angle_restraint[iangle].angle();
	 double esd = restraints.angle_restraint[iangle].esd();
	 PyObject *angle_restraint = PyList_New(5);
	 PyList_SetItem(angle_restraint, 0, PyUnicode_FromString(a1.c_str()));
	 PyList_SetItem(angle_restraint, 1, PyUnicode_FromString(a2.c_str()));
	 PyList_SetItem(angle_restraint, 2, PyUnicode_FromString(a3.c_str()));
	 PyList_SetItem(angle_restraint, 3, PyFloat_FromDouble(d));
	 PyList_SetItem(angle_restraint, 4, PyFloat_FromDouble(esd));
	 PyList_SetItem(angle_restraint_list, iangle, angle_restraint);
      }

      PyDict_SetItem(r, PyUnicode_FromString("_chem_comp_angle"), angle_restraint_list);

      
      // ------------------ Torsions -------------------------
      PyObject *torsion_restraint_list = PyList_New(restraints.torsion_restraint.size());
      for (unsigned int itorsion=0; itorsion<restraints.torsion_restraint.size(); itorsion++) {
	 std::string id = restraints.torsion_restraint[itorsion].id();
	 std::string a1 = restraints.torsion_restraint[itorsion].atom_id_1_4c();
	 std::string a2 = restraints.torsion_restraint[itorsion].atom_id_2_4c();
	 std::string a3 = restraints.torsion_restraint[itorsion].atom_id_3_4c();
	 std::string a4 = restraints.torsion_restraint[itorsion].atom_id_4_4c();
	 double tor  = restraints.torsion_restraint[itorsion].angle();
	 double esd = restraints.torsion_restraint[itorsion].esd();
	 int period = restraints.torsion_restraint[itorsion].periodicity();
	 PyObject *torsion_restraint = PyList_New(8);
	 PyList_SetItem(torsion_restraint, 0, PyUnicode_FromString(id.c_str()));
	 PyList_SetItem(torsion_restraint, 1, PyUnicode_FromString(a1.c_str()));
	 PyList_SetItem(torsion_restraint, 2, PyUnicode_FromString(a2.c_str()));
	 PyList_SetItem(torsion_restraint, 3, PyUnicode_FromString(a3.c_str()));
	 PyList_SetItem(torsion_restraint, 4, PyUnicode_FromString(a4.c_str()));
	 PyList_SetItem(torsion_restraint, 5, PyFloat_FromDouble(tor));
	 PyList_SetItem(torsion_restraint, 6, PyFloat_FromDouble(esd));
	 PyList_SetItem(torsion_restraint, 7, PyLong_FromLong(period));
	 PyList_SetItem(torsion_restraint_list, itorsion, torsion_restraint);
      }

      PyDict_SetItem(r, PyUnicode_FromString("_chem_comp_tor"), torsion_restraint_list);

      // ------------------ Planes -------------------------
      PyObject *plane_restraints_list = PyList_New(restraints.plane_restraint.size());
      for (unsigned int iplane=0; iplane<restraints.plane_restraint.size(); iplane++) {
	 PyObject *atom_list = PyList_New(restraints.plane_restraint[iplane].n_atoms());
	 for (int iat=0; iat<restraints.plane_restraint[iplane].n_atoms(); iat++) { 
	    const std::string &at = restraints.plane_restraint[iplane][iat].first;
	    double atom_esd = restraints.plane_restraint[iplane][iat].second;
	    PyObject *atom_pair_list = PyList_New(2);
	    PyObject *n = PyUnicode_FromString(at.c_str());
	    PyObject *e = PyFloat_FromDouble(atom_esd);
	    PyList_SetItem(atom_pair_list, 0, n);
	    PyList_SetItem(atom_pair_list, 1, e);
	    // PyList_SetItem(atom_list, iat, n);
	    PyList_SetItem(atom_list, iat, atom_pair_list);
	 }

	 // HACK HACK! // FIXME
// 	 double esd = restraints.plane_restraint[iplane].dist_esd(0);
// 	 PyObject *plane_restraint = PyList_New(3);
// 	 PyList_SetItem(plane_restraint, 0, PyString_FromString(restraints.plane_restraint[iplane].plane_id.c_str()));
// 	 PyList_SetItem(plane_restraint, 1, atom_list);
// 	 PyList_SetItem(plane_restraint, 2, PyFloat_FromDouble(esd));
// 	 PyList_SetItem(plane_restraints_list, iplane, plane_restraint);
	 
	 PyObject *plane_restraint = PyList_New(2);
	 PyList_SetItem(plane_restraint, 0, PyUnicode_FromString(restraints.plane_restraint[iplane].plane_id.c_str()));
	 PyList_SetItem(plane_restraint, 1, atom_list);
	 PyList_SetItem(plane_restraints_list, iplane, plane_restraint);
      }

      PyDict_SetItem(r, PyUnicode_FromString("_chem_comp_plane_atom"), plane_restraints_list);

      // ------------------ Chirals -------------------------
      PyObject *chiral_restraint_list = PyList_New(restraints.chiral_restraint.size());
      for (unsigned int ichiral=0; ichiral<restraints.chiral_restraint.size(); ichiral++) {
	 
	 std::string a1 = restraints.chiral_restraint[ichiral].atom_id_1_4c();
	 std::string a2 = restraints.chiral_restraint[ichiral].atom_id_2_4c();
	 std::string a3 = restraints.chiral_restraint[ichiral].atom_id_3_4c();
	 std::string ac = restraints.chiral_restraint[ichiral].atom_id_c_4c();
	 std::string chiral_id = restraints.chiral_restraint[ichiral].Chiral_Id();

	 double esd = restraints.chiral_restraint[ichiral].volume_sigma();
	 int volume_sign = restraints.chiral_restraint[ichiral].volume_sign;
	 PyObject *chiral_restraint = PyList_New(7);
	 PyList_SetItem(chiral_restraint, 0, PyUnicode_FromString(chiral_id.c_str()));
	 PyList_SetItem(chiral_restraint, 1, PyUnicode_FromString(ac.c_str()));
	 PyList_SetItem(chiral_restraint, 2, PyUnicode_FromString(a1.c_str()));
	 PyList_SetItem(chiral_restraint, 3, PyUnicode_FromString(a2.c_str()));
	 PyList_SetItem(chiral_restraint, 4, PyUnicode_FromString(a3.c_str()));
	 PyList_SetItem(chiral_restraint, 5, PyLong_FromLong(volume_sign));
	 PyList_SetItem(chiral_restraint, 6, PyFloat_FromDouble(esd));
	 PyList_SetItem(chiral_restraint_list, ichiral, chiral_restraint);
      }

      PyDict_SetItem(r, PyUnicode_FromString("_chem_comp_chir"), chiral_restraint_list);
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // LIBCOOTAPI_BUILD
