
%module pyrogen_swig

%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectors) vector<string>;
};

%{
#include <stddef.h>
#include "pyrogen/restraints.hh"
%} 

%include "pyrogen/restraints.hh"
