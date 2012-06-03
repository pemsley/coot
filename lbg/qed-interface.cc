
#ifdef MAKE_ENTERPRISE_TOOLS

#include <Python.h>

#include "qed-interface.hh"
#include "rdkit-interface.hh"

double
get_qed(const RDKit::ROMol &rdkm) {

   double r = -1.0; 
   PyObject *arg_list;
   PyObject *my_callback;
   
   // Build the name object
   PyObject *pName = PyString_FromString("silicos_it.descriptors");
   std::cout << "pName " << pName << std::endl;
   
   // Load the module object
   PyObject *pModule = PyImport_Import(pName);
   std::cout << "pModule " << pModule << std::endl;

   if (PyModule_Check(pModule)) {
      std::cout << "pModule is a module"<< std::endl;
   } else {
      std::cout << "pModule is not a module"<< std::endl;
   } 
   
   // pDict is a borrowed reference 
   PyObject *pDict = PyModule_GetDict(pModule);
   std::cout << "pDict " << pModule << std::endl;
   if (PyDict_Check(pDict)) {
      std::cout << "pDict is a dict"<< std::endl;
   } else {
      std::cout << "pDict is not a dict"<< std::endl;
   } 
   
   // pFunc is also a borrowed reference 
   my_callback = PyDict_GetItemString(pDict, "qed.default()");
   std::cout << "my_callback " << my_callback << std::endl;
   
   if (PyCallable_Check(my_callback)) {
      std::cout << "is callable!"  << std::endl;
      PyObject *result = PyEval_CallObject(my_callback, arg_list);
      if (PyFloat_Check(result)) { 
	 r = PyFloat_AsDouble(result);
      }
   } else {
      std::cout << "function is not callable"  << std::endl;
   } 
   return r;
}

#endif

