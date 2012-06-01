
namespace coot {

   // private (no SWIG interface)
   // 
   // the engine for the above calls
   std::pair<CMMDBManager *, CResidue *>
   regularize_inner(PyObject *rdkit_mol,
		    PyObject *restraints_py,
		    const std::string &res_name);

}

