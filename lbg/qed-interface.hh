
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#ifdef USE_PYTHON

#include "lbg.hh"

double get_qed(PyObject *silicos_it_qed_default_func,
	       const RDKit::ROMol &rdkm);

std::vector<double>
get_qed_properties(PyObject *silicos_it_qed_properties_func, 
                   PyObject *silicos_it_qed_pads, 
		   const RDKit::ROMol &rdkm);

double
get_qed_ads(const std::vector<double> &properties, PyObject *pads, long idx);


#endif
#endif

