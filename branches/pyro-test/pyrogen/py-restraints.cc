

#include <boost/python.hpp>
#include "restraints.hh"
#include "py-restraints.hh"
#include <rdkit-interface.hh>
#include <coot-coord-utils.hh>

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
	
   coot::dictionary_residue_restraints_t monomer_restraints;
   monomer_restraints.bond_restraint    = bond_restraints;
   monomer_restraints.angle_restraint   = angle_restraints;
   monomer_restraints.torsion_restraint = torsion_restraints;
   monomer_restraints.chiral_restraint  = chiral_restraints;
   monomer_restraints.plane_restraint   = plane_restraints;
   monomer_restraints.residue_info      = residue_info;
   monomer_restraints.atom_info         = atoms; 
   return monomer_restraints;
} 


PyObject *coot::monomer_restraints_to_python(const dictionary_residue_restraints_t &restraints) {

   PyObject *r = Py_False;

   if (1) { 
      r = PyDict_New();
	 
      // ------------------ chem_comp -------------------------
      coot::dict_chem_comp_t info = restraints.residue_info;
      
      PyObject *chem_comp_py = PyList_New(7);
      PyList_SetItem(chem_comp_py, 0, PyString_FromString(info.comp_id.c_str()));
      PyList_SetItem(chem_comp_py, 1, PyString_FromString(info.three_letter_code.c_str()));
      PyList_SetItem(chem_comp_py, 2, PyString_FromString(info.name.c_str()));
      PyList_SetItem(chem_comp_py, 3, PyString_FromString(info.group.c_str()));
      PyList_SetItem(chem_comp_py, 4, PyInt_FromLong(info.number_atoms_all));
      PyList_SetItem(chem_comp_py, 5, PyInt_FromLong(info.number_atoms_nh));
      PyList_SetItem(chem_comp_py, 6, PyString_FromString(info.description_level.c_str()));
      
      // Put chem_comp_py into a dictionary?
      PyDict_SetItem(r, PyString_FromString("_chem_comp"), chem_comp_py);

      
      // ------------------ chem_comp_atom -------------------------
      std::vector<coot::dict_atom> atom_info = restraints.atom_info;
      int n_atoms = atom_info.size();
      PyObject *atom_info_list = PyList_New(n_atoms);
      for (int iat=0; iat<n_atoms; iat++) { 
	 PyObject *atom_attributes_list = PyList_New(5);
	 PyList_SetItem(atom_attributes_list, 0, PyString_FromString(atom_info[iat].atom_id_4c.c_str()));
	 PyList_SetItem(atom_attributes_list, 1, PyString_FromString(atom_info[iat].type_symbol.c_str()));
	 PyList_SetItem(atom_attributes_list, 2, PyString_FromString(atom_info[iat].type_energy.c_str()));
	 PyList_SetItem(atom_attributes_list, 3, PyFloat_FromDouble(atom_info[iat].partial_charge.second));
	 PyObject *flag = Py_False;
	 if (atom_info[iat].partial_charge.first)
	    flag = Py_True;
     Py_INCREF(flag);
	 PyList_SetItem(atom_attributes_list, 4, flag);
	 PyList_SetItem(atom_info_list, iat, atom_attributes_list);
      }

      PyDict_SetItem(r, PyString_FromString("_chem_comp_atom"), atom_info_list);

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
	 catch (std::runtime_error rte) {

	    // Use default false values.
	    // So I suppose that I need to do this then:
	    if (PyBool_Check(py_value_dist))
	       Py_INCREF(py_value_dist);
	    if (PyBool_Check(py_value_esd))
	       Py_INCREF(py_value_esd);
	 }
	 
	 PyObject *bond_restraint = PyList_New(5);
	 PyList_SetItem(bond_restraint, 0, PyString_FromString(a1.c_str()));
	 PyList_SetItem(bond_restraint, 1, PyString_FromString(a2.c_str()));
	 PyList_SetItem(bond_restraint, 2, PyString_FromString(type.c_str()));
	 PyList_SetItem(bond_restraint, 3, py_value_dist);
	 PyList_SetItem(bond_restraint, 4, py_value_esd);
	 PyList_SetItem(bond_restraint_list, ibond, bond_restraint);
      }

      PyDict_SetItem(r, PyString_FromString("_chem_comp_bond"), bond_restraint_list);


      // ------------------ Angles -------------------------
      PyObject *angle_restraint_list = PyList_New(restraints.angle_restraint.size());
      for (unsigned int iangle=0; iangle<restraints.angle_restraint.size(); iangle++) {
	 std::string a1 = restraints.angle_restraint[iangle].atom_id_1_4c();
	 std::string a2 = restraints.angle_restraint[iangle].atom_id_2_4c();
	 std::string a3 = restraints.angle_restraint[iangle].atom_id_3_4c();
	 double d   = restraints.angle_restraint[iangle].angle();
	 double esd = restraints.angle_restraint[iangle].esd();
	 PyObject *angle_restraint = PyList_New(5);
	 PyList_SetItem(angle_restraint, 0, PyString_FromString(a1.c_str()));
	 PyList_SetItem(angle_restraint, 1, PyString_FromString(a2.c_str()));
	 PyList_SetItem(angle_restraint, 2, PyString_FromString(a3.c_str()));
	 PyList_SetItem(angle_restraint, 3, PyFloat_FromDouble(d));
	 PyList_SetItem(angle_restraint, 4, PyFloat_FromDouble(esd));
	 PyList_SetItem(angle_restraint_list, iangle, angle_restraint);
      }

      PyDict_SetItem(r, PyString_FromString("_chem_comp_angle"), angle_restraint_list);

      
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
	 PyList_SetItem(torsion_restraint, 0, PyString_FromString(id.c_str()));
	 PyList_SetItem(torsion_restraint, 1, PyString_FromString(a1.c_str()));
	 PyList_SetItem(torsion_restraint, 2, PyString_FromString(a2.c_str()));
	 PyList_SetItem(torsion_restraint, 3, PyString_FromString(a3.c_str()));
	 PyList_SetItem(torsion_restraint, 4, PyString_FromString(a4.c_str()));
	 PyList_SetItem(torsion_restraint, 5, PyFloat_FromDouble(tor));
	 PyList_SetItem(torsion_restraint, 6, PyFloat_FromDouble(esd));
	 PyList_SetItem(torsion_restraint, 7, PyInt_FromLong(period));
	 PyList_SetItem(torsion_restraint_list, itorsion, torsion_restraint);
      }

      PyDict_SetItem(r, PyString_FromString("_chem_comp_tor"), torsion_restraint_list);

      // ------------------ Planes -------------------------
      PyObject *plane_restraints_list = PyList_New(restraints.plane_restraint.size());
      for (unsigned int iplane=0; iplane<restraints.plane_restraint.size(); iplane++) {
	 PyObject *atom_list = PyList_New(restraints.plane_restraint[iplane].n_atoms());
	 for (int iat=0; iat<restraints.plane_restraint[iplane].n_atoms(); iat++) { 
	    std::string at = restraints.plane_restraint[iplane][iat];
	    PyList_SetItem(atom_list, iat, PyString_FromString(at.c_str()));
	 }
	 double esd = restraints.plane_restraint[iplane].dist_esd();
	 PyObject *plane_restraint = PyList_New(3);
	 PyList_SetItem(plane_restraint, 0, PyString_FromString(restraints.plane_restraint[iplane].plane_id.c_str()));
	 PyList_SetItem(plane_restraint, 1, atom_list);
	 PyList_SetItem(plane_restraint, 2, PyFloat_FromDouble(esd));
	 PyList_SetItem(plane_restraints_list, iplane, plane_restraint);
      }

      PyDict_SetItem(r, PyString_FromString("_chem_comp_plane_atom"), plane_restraints_list);

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
	 PyList_SetItem(chiral_restraint, 0, PyString_FromString(chiral_id.c_str()));
	 PyList_SetItem(chiral_restraint, 1, PyString_FromString(ac.c_str()));
	 PyList_SetItem(chiral_restraint, 2, PyString_FromString(a1.c_str()));
	 PyList_SetItem(chiral_restraint, 3, PyString_FromString(a2.c_str()));
	 PyList_SetItem(chiral_restraint, 4, PyString_FromString(a3.c_str()));
	 PyList_SetItem(chiral_restraint, 5, PyInt_FromLong(volume_sign));
	 PyList_SetItem(chiral_restraint, 6, PyFloat_FromDouble(esd));
	 PyList_SetItem(chiral_restraint_list, ichiral, chiral_restraint);
      }

      PyDict_SetItem(r, PyString_FromString("_chem_comp_chir"), chiral_restraint_list);
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
