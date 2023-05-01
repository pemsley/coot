
%module pyrogen_swig

%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectors) vector<string>;
};

%{
#include <cstddef>
#include <lidia-core/use-rdkit.hh>
#include "pyrogen/restraints.hh"
// #include "pyrogen/geometry-store.hh"
// #include "pyrogen/geometry-store-interface.hh"
%}

%include "pyrogen/restraints.hh"

// %include "mmdb2.i"
// %include "pyrogen/geometry-store.hh"
// %include "pyrogen/geometry-store-interface.hh"
