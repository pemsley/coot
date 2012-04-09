
#include "Python.h"
#include <string>
#include <iostream>
#include <vector>

#include "protein-geometry.hh"

PyObject *set_monomer_restraints_py(const char *monomer_type, PyObject *restraints) {

   PyObject *retval = Py_False;

   if (!PyDict_Check(restraints)) {
      std::cout << " Failed to read restraints - not a list" << std::endl;
   } else {

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

      std::cout << "looping over restraint" << std::endl;
      while (PyDict_Next(restraints, &pos, &key, &value)) {
	 // std::cout << ":::::::key: " << PyString_AsString(key) << std::endl;

	 std::string key_string = PyString_AsString(key);
	 if (key_string == "_chem_comp") {
	    PyObject *chem_comp_list = value;
	    if (PyList_Check(chem_comp_list)) {
	       if (PyObject_Length(chem_comp_list) == 7) {
		  std::string comp_id  = PyString_AsString(PyList_GetItem(chem_comp_list, 0));
		  std::string tlc      = PyString_AsString(PyList_GetItem(chem_comp_list, 1));
		  std::string name     = PyString_AsString(PyList_GetItem(chem_comp_list, 2));
		  std::string group    = PyString_AsString(PyList_GetItem(chem_comp_list, 3));
		  int n_atoms_all      = PyInt_AsLong(PyList_GetItem(chem_comp_list, 4));
		  int n_atoms_nh       = PyInt_AsLong(PyList_GetItem(chem_comp_list, 5));
		  std::string desc_lev = PyString_AsString(PyList_GetItem(chem_comp_list, 6));

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
		  if (PyObject_Length(chem_comp_atom) == 5) {
		     std::string atom_id  = PyString_AsString(PyList_GetItem(chem_comp_atom, 0));
		     std::string element  = PyString_AsString(PyList_GetItem(chem_comp_atom, 1));
		     std::string energy_t = PyString_AsString(PyList_GetItem(chem_comp_atom, 2));
		     float part_chr        = PyFloat_AsDouble(PyList_GetItem(chem_comp_atom, 3));
		     bool flag = 0;
		     if (PyLong_AsLong(PyList_GetItem(chem_comp_atom, 4))) {
			flag = 1;
		     }
		     std::pair<bool, float> part_charge_info(flag, part_chr);
		     coot::dict_atom at(atom_id, atom_id, element, energy_t, part_charge_info);
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

		     if (PyString_Check(atom_1_py) &&
			 PyString_Check(atom_2_py) &&
			 PyString_Check(type_py) &&
			 PyFloat_Check(dist_py) && 
			 PyFloat_Check(esd_py)) {
			std::string atom_1 = PyString_AsString(atom_1_py);
			std::string atom_2 = PyString_AsString(atom_2_py);
			std::string type   = PyString_AsString(type_py);
			float  dist = PyFloat_AsDouble(dist_py);
			float  esd  = PyFloat_AsDouble(esd_py);
			coot::dict_bond_restraint_t rest(atom_1, atom_2, type, dist, esd);
			bond_restraints.push_back(rest);
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

		     if (PyString_Check(atom_1_py) &&
			 PyString_Check(atom_2_py) &&
			 PyString_Check(atom_3_py) &&
			 PyFloat_Check(angle_py) && 
			 PyFloat_Check(esd_py)) {
			std::string atom_1 = PyString_AsString(atom_1_py);
			std::string atom_2 = PyString_AsString(atom_2_py);
			std::string atom_3 = PyString_AsString(atom_3_py);
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

		     if (PyString_Check(atom_1_py) &&
			 PyString_Check(atom_2_py) &&
			 PyString_Check(atom_3_py) &&
			 PyString_Check(atom_4_py) &&
			 PyFloat_Check(torsion_py) && 
			 PyFloat_Check(esd_py)    && 
			 PyInt_Check(period_py)) { 
			std::string id     = PyString_AsString(id_py);
			std::string atom_1 = PyString_AsString(atom_1_py);
			std::string atom_2 = PyString_AsString(atom_2_py);
			std::string atom_3 = PyString_AsString(atom_3_py);
			std::string atom_4 = PyString_AsString(atom_4_py);
			float  torsion = PyFloat_AsDouble(torsion_py);
			float  esd     = PyFloat_AsDouble(esd_py);
			int  period    = PyInt_AsLong(period_py);
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

		     if (PyString_Check(atom_1_py) &&
			 PyString_Check(atom_2_py) &&
			 PyString_Check(atom_3_py) &&
			 PyString_Check(atom_c_py) &&
			 PyString_Check(chiral_id_py) && 
			 PyFloat_Check(esd_py)    && 
			 PyInt_Check(vol_sign_py)) {
			std::string chiral_id = PyString_AsString(chiral_id_py);
			std::string atom_c    = PyString_AsString(atom_c_py);
			std::string atom_1    = PyString_AsString(atom_1_py);
			std::string atom_2    = PyString_AsString(atom_2_py);
			std::string atom_3    = PyString_AsString(atom_3_py);
			float  esd            = PyFloat_AsDouble(esd_py);
			int  vol_sign         = PyInt_AsLong(vol_sign_py);
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
			   if (PyString_Check(at_py)) {
			      atoms.push_back(PyString_AsString(at_py));
			   } else {
			      atoms_pass = 0;
			   }
			}
			if (atoms_pass) {
			   if (PyString_Check(plane_id_py)) {
			      if (PyFloat_Check(esd_py)) {
				 std::string plane_id = PyString_AsString(plane_id_py);
				 float esd = PyFloat_AsDouble(esd_py);
				 if (atoms.size() > 0) { 
				    coot::dict_plane_restraint_t rest(plane_id, atoms[0], esd);
				    for (int i=1; i<atoms.size(); i++)
				       rest.push_back_atom(atoms[i]);
				    plane_restraints.push_back(rest);
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
	
      coot::dictionary_residue_restraints_t monomer_restraints(monomer_type, 1);
      monomer_restraints.bond_restraint    = bond_restraints;
      monomer_restraints.angle_restraint   = angle_restraints;
      monomer_restraints.torsion_restraint = torsion_restraints;
      monomer_restraints.chiral_restraint  = chiral_restraints;
      monomer_restraints.plane_restraint   = plane_restraints;
      monomer_restraints.residue_info      = residue_info;
      monomer_restraints.atom_info         = atoms; 
      
      bool s = false;
          // g.Geom_p()->replace_monomer_restraints(monomer_type, monomer_restraints);
      if (s)
	 retval = Py_True;
      
   }
   if (PyBool_Check(retval)) {
     Py_INCREF(retval);
   }
   return retval;
} 
