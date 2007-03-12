
#ifndef REGULARIZE_MINIMOL_HH
#define REGULARIZE_MINIMOL_HH

#ifdef HAVE_GSL

#include "mini-mol.hh"
#include "simple-restraint.hh"


namespace coot {

   minimol::molecule
   regularize_minimol_molecule(const minimol::molecule &molin,
			       const protein_geometry &geom);

}

#endif // HAVE_GSL
#endif // REGULARIZE_MINIMOL_HH

