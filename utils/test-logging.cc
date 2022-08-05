
#include <iostream>
#include "logging.hh"

int main(int argc, char **argv) {

   logging log;

   log.log(std::string("WARNING", "first message"));
   log.log(std::string("test plain message"));

   log.show();

   // coot_log << "Something here" << std::endl;

   return 0;

}
