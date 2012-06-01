
%module coot_libs

%include "std_string.i"
%include "std_vector.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectors) vector<string>;
};

%{
#include "restraints.hh"
%} 

%include "restraints.hh"
