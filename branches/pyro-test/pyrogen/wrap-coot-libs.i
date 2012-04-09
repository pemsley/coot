
%module coot_libs

%{
PyObject *set_monomer_restraints_py(const char *monomer_type, PyObject *restraints);
%} 

PyObject *set_monomer_restraints_py(const char *monomer_type, PyObject *restraints);

