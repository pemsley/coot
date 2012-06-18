
#ifdef MAKE_ENTERPRISE_TOOLS
#ifdef USE_PYTHON

#include <boost/python.hpp>
using namespace boost::python;

#include "qed-interface.hh"
#include "rdkit-interface.hh"


// try to get the python function, or set it to null.
void
lbg_info_t::setup_silicos_it_qed_default_func() {

   // Build the name object
   PyObject *pName = PyString_FromString("silicos_it.descriptors.qed");
   // Load the module object
   PyObject *pModule = PyImport_Import(pName);
   if (pModule == NULL) {
      // This happens when biscu-it is not installed (which is often(?) the case)
      // 
      // std::cout << "Null pModule in get_qed() " << std::endl;
   } else { 
      // pDict is a borrowed reference 
      PyObject *pDict = PyModule_GetDict(pModule);
      if (! PyDict_Check(pDict)) {
	 std::cout << "pDict is not a dict"<< std::endl;
      } else {
	 // pFunc is also a borrowed reference 
	 PyObject *pFunc = PyDict_GetItemString(pDict, "default");
	 if (! pFunc) {
	    // std::cout << "pFunc is NULL" << std::endl;
	 } else { 
	    if (PyCallable_Check(pFunc)) {
	       silicos_it_qed_default_func = pFunc;
	    } else {
	       std::cout << "function is not callable"  << std::endl;
	    }
	 }
      }
   }
}



// needs a properly initialized (say, via rdkit_mol_sanitize()) molecule.
double
get_qed(PyObject *silicos_it_qed_default_func, const RDKit::ROMol &rdkm) {

   double r = -1.0;
   if (silicos_it_qed_default_func) { 
      try { 
	 PyObject *arg_list = PyTuple_New(1);
	 PyObject *rdkit_mol_py; // = PyInt_FromLong(6);
	    
	 RDKit::ROMol *mol_copy_p = new RDKit::ROMol(rdkm);
	 boost::shared_ptr<RDKit::ROMol> xx(mol_copy_p);
	 boost::python::object obj(xx);
	 rdkit_mol_py = obj.ptr();
	 PyTuple_SetItem(arg_list, 0, rdkit_mol_py);
	 PyObject *result = PyEval_CallObject(silicos_it_qed_default_func, arg_list);
	 if (! result) {
	    std::cout << "Null result" << std::endl;
	 } else { 
	    if (PyFloat_Check(result)) { 
	       r = PyFloat_AsDouble(result);
	    }
	 }

	 // OK, I'm finished with mol_copy_p and obj now.  How do I
	 // get rid of them?  (Does that happen automatically?)
	    
      }
      catch (boost::python::error_already_set e) {
	 std::cout << "catch error_already_set exception "
		   << std::endl;
	 PyObject *type_ptr = NULL, *value_ptr = NULL, *traceback_ptr = NULL;
	 PyErr_Fetch(&type_ptr, &value_ptr, &traceback_ptr);

	 PyObject *dest = PyString_FromString("object: %s\n");
	 if (type_ptr)
	    std::cout << "error: type "
		      << PyString_AsString(PyString_Format(dest, type_ptr)) << std::endl;
	 if (value_ptr)
	    std::cout << "error: value "
		      << PyString_AsString(PyString_Format(dest, value_ptr)) << std::endl;
	 if (traceback_ptr)
	    std::cout << "error: traceback_ptr "
		      << PyString_AsString(PyString_Format(dest, traceback_ptr)) << std::endl;
	    
      } 
      catch (...) {
	 std::cout << "catch all exception" << std::endl;
      }
   }
   return r;
}


#endif
#endif

