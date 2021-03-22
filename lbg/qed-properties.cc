/* lbg/qed-properties.cc
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2016 by Medical Research Council
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

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <boost/python.hpp>

#include "qed-interface.hh"
#include "lidia-core/rdkit-interface.hh"

std::vector<std::pair<double, double> >
get_qed_properties(PyObject *silicos_it_qed_properties_func,
		   PyObject *silicos_it_qed_pads,
		   const RDKit::ROMol &rdkm) {

   std::vector<std::pair<double, double> > p(8, std::pair<double, double>(0,0));
   
   if (silicos_it_qed_properties_func) {
      try {
	 PyObject *arg_list = PyTuple_New(1);
	 PyObject *rdkit_mol_py;
	    
	 RDKit::ROMol *mol_copy_p = new RDKit::ROMol(rdkm);
	 boost::shared_ptr<RDKit::ROMol> xx(mol_copy_p);
	 boost::python::object obj(xx);
	 rdkit_mol_py = obj.ptr();
	 PyTuple_SetItem(arg_list, 0, rdkit_mol_py);

	 PyObject *result_properties = PyEval_CallObject(silicos_it_qed_properties_func, arg_list);
	 if (PyList_Check(result_properties)) {
	    int l = PyObject_Length(result_properties);
	    std::vector<double> properties(l,0);
	    std::vector<std::string> prop_names(8);
	    prop_names[0] = "MolWeght ";
	    prop_names[1] = "MolLogP  ";
	    prop_names[2] = "#HBA     ";
	    prop_names[3] = "#Donors  ";
	    prop_names[4] = "TPSA     ";
	    prop_names[5] = "#RotBonds";
	    prop_names[6] = "#Arom    ";
	    prop_names[7] = "#Alerts  ";
	    for (long i=0; i<l; i++) {
	       PyObject *item_py = PyList_GetItem(result_properties, i);
	       // 0: MolWeght  1: LogP       2: #HBA   3: #Donors
	       // 4: TPSA      5: #RotBonds  6: Arom   7: #Alerts
	       if (PyFloat_Check(item_py)) {
		  double item_f = PyFloat_AsDouble(item_py);
		  properties[i] = item_f;
		  // std::cout << "setting float " << i << " " << item_f << std::endl;
	       } else {
		  if (PyLong_Check(item_py)) {
		     long item_i = PyLong_AsLong(item_py);
		     properties[i] = item_i;
		     // std::cout << "setting int " << i << " " << item_i << std::endl;
		  }
	       }
	    }

	    for (long i=0; i<8; i++) {
	       std::pair<double, double> s = get_qed_ads(properties, silicos_it_qed_pads, i);
	       p[i] = s;
	       if (false) // debug
		  std::cout << i << " " << prop_names[i] << "  " << properties[i]
			    << " " << s.first << " " << s.second << std::endl;
	    }
	 }
      }
      catch (const boost::python::error_already_set &e) {
	 std::cout << "catch error_already_set exception "
		   << std::endl;
	 PyObject *type_ptr = NULL, *value_ptr = NULL, *traceback_ptr = NULL;
	 PyErr_Fetch(&type_ptr, &value_ptr, &traceback_ptr);

	 PyObject *dest = PyUnicode_FromString("object: %s\n");
	 if (type_ptr)
	    std::cout << "error: type " << PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyUnicode_Format(dest, type_ptr))) << std::endl;
	 if (value_ptr)
	    std::cout << "error: value " << PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyUnicode_Format(dest, value_ptr))) << std::endl;
	 if (traceback_ptr)
	    std::cout << "error: traceback_ptr " << PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyUnicode_Format(dest, traceback_ptr)))
                      << std::endl;
	    
      }
      catch (...) {
	 std::cout << "catch all exception" << std::endl;
      }
   }
   return p;
}


// return (value, desirability)
std::pair<double, double>
get_qed_ads(const std::vector<double> &properties, PyObject *pads, long idx) {

   std::pair<double,double> r(0,0);
   if (! PyList_Check(pads)) {
      std::cout << "ERROR:: get_qed_ads() pads is not a list " << std::endl;
   } else {
      long l = PyObject_Length(pads);
      if (l != 8) {
	 std::cout << "ERROR:: get_qed_ads() pads length is not 8. " << l << std::endl;
      } else {
	 PyObject *list_py = PyList_GetItem(pads, idx);
	 if (! PyList_Check(list_py)) {
	    std::cout << "ERROR:: get_qed_ads() pads[i] is not a list " << std::endl;
	 } else {
	    long l2 = PyObject_Length(list_py);
	    if (l2 != 7) {
	       std::cout << "ERROR:: get_qed_ads() pads[i] length is not 7. " << l << std::endl;
	    } else {
	       std::vector<double> coeffs(7,0);
	       for (unsigned int i=0; i<7; i++) { 
		  PyObject *item_py = PyList_GetItem(list_py, i);
		  if (PyFloat_Check(item_py)) {
		     double f = PyFloat_AsDouble(item_py);
		     coeffs[i] = f;
		  }
	       }

	       const double &x=properties[idx];

	       const double &a=coeffs[0];
	       const double &b=coeffs[1];
	       const double &c=coeffs[2];
	       const double &d=coeffs[3];
	       const double &e=coeffs[4];
	       const double &f=coeffs[5];
	       const double &dmax=coeffs[6];

	       double desire =
		  ((a+(b/(1+exp(-1*(x-c+d/2)/e))*(1-1/(1+exp(-1*(x-c-d/2)/f)))))/dmax);

	       r = std::pair<double, double> (x, desire);
	    }
	 }
      }
   }
   return r;
}


#endif // MAKE_ENHANCED_LIGAND_TOOLS
