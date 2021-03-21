
#ifndef QED_INTERFACE_HH
#define QED_INTERFACE_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include "lbg.hh"

double get_qed(PyObject *silicos_it_qed_default_func,
	       const RDKit::ROMol &rdkm);

// get the values and their desirabilites (0->1)
std::vector<std::pair<double, double> >
get_qed_properties(PyObject *silicos_it_qed_properties_func, 
                   PyObject *silicos_it_qed_pads, 
		   const RDKit::ROMol &rdkm);

std::pair<double, double>
get_qed_ads(const std::vector<double> &properties, PyObject *pads, long idx);


#endif



#endif // QED_INTERFACE_HH

