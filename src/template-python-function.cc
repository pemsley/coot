
#include "Python.h"
#include "c-interface.h"

#ifdef USE_PYTHON
//! \brief return the summary info for ligand distortion
PyObject *get_ligand_distortion_summary_info_py(int imol, PyObject *residue_spec) {

   PyObject *r = Py_False;
  
   if (is_valid_model_molecule(imol)) {
      
   }
  if (PyBool_Check(r)) {
    Py_INCREF(r);
  }
  return r;
}
#endif

