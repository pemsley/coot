

#include "logging.hh"

int main(int argc, char **argv) {

   logging log;

   log.log("WARNING", "first message");
   log.log("test plain message");

   log.show();

   return 0;

}
