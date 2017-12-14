
#ifndef C_INTERFACE_BONDS_HH
#define C_INTERFACE_BONDS_HH

// change the name of this - it's not just bonds

// c-interface-unexported or some such

#include <clipper/core/coords.h>

// not to be exported to the API
#ifdef USE_PYTHON

PyObject *pyobject_from_graphical_bonds_container(const graphical_bonds_container &gbc);
PyObject *go_to_ligand_py();
clipper::Coord_orth go_to_ligand_inner();

#endif // USE_PYTHON

#endif // C_INTERFACE_BONDS_HH
