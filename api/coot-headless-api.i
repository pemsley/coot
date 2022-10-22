
%module coot_headless_api

%{
#include "molecules_container.hh"
}%

%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"

namespace std {
%template(vector_string) vector<std::string>;
%template(pairbf) pair<bool, float>;
%template(IntVector) vector<int>;
%template(vector_atom_spec)      vector<coot::atom_spec_t>;
%template(vector_mtz_type_label) vector<coot::mtz_type_label>;
}

%init %{
  // init_coot_as_python_module();
%}

%feature("autodoc", "1"); // add doc string for Intellisense (hopefully)

%include "molecules_container.hh"


