
#include <boost/python.hpp>
using namespace boost::python;

#define MAKE_ENTERPRISE_TOOLS

#include <GraphMol/GraphMol.h>

namespace coot {
   // try without confusion of also passing a PyObject *...
   RDKit::ROMol *new_regularize(RDKit::ROMol &r);
}

BOOST_PYTHON_MODULE(coot_libs) {
   def("new_regularize",  coot::new_regularize);
}
