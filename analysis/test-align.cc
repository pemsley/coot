
#include "utils/pir-alignment.hh"

int main(int argc, char **argv) {

   int status = 0;
   std::string fn = "ester/apc5-new2.pir";

   coot::pir_alignment_t a;
   a.read_file(fn);

   return status;

}
