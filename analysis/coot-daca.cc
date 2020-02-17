
#include "daca.hh"

int main(int argc, char **argv) {

   int status = 0;

   coot::daca daca;

   daca.write_tables_using_reference_structures_from_dir("test-pdbs");

   return status;
}
