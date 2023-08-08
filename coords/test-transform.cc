
#include "mmdb-extras.h"
#include "mmdb.hh"
#include "mmdb-crystal.h" 


int main(int argc, char **argv) {

   mmdb::mat44 my_matt;
   mmdb::Atom *at = NULL;
   at->Transform(my_matt);
   return 0;
}
