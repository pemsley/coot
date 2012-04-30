
%module coot_libs

%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectors) vector<string>;
};

%{
#include "mogul-interface.hh"
#include "restraints.hh"
%} 

%include "mogul-interface.hh"
%include "restraints.hh"
