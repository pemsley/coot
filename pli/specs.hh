
#ifndef PY_SPECS_HH
#define PY_SPECS_HH

#include <Python.h>

#include "geometry/residue-and-atom-specs.hh"

namespace coot {

   class py_atom_spec_t : public atom_spec_t {
   public:
      py_atom_spec_t(PyObject *obj);
      py_atom_spec_t(const atom_spec_t &spec_in) : atom_spec_t(spec_in) {}
      PyObject *pyobject() const;
   };

   class py_residue_spec_t : public residue_spec_t {
   public:
      py_residue_spec_t(PyObject *obj);
      py_residue_spec_t(const residue_spec_t &spec_in) : residue_spec_t(spec_in) {}
      PyObject *pyobject() const;
   };

}

#endif // PY_SPECS_HH
