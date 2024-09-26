#include <string>
#include <iostream>

namespace mmcif_tests {
   int run_tests(bool last_test_only);
}

int main(int argc, char **argv) {

   int status = 0;
   bool last_test_only = false;
   if (argc > 1) {
      std::string arg(argv[1]);
      if (arg == "last-test-only")
         last_test_only = true;
   }
   int test_status =  mmcif_tests::run_tests(last_test_only);
   if (test_status == 0) status = 0; // convert to shell-type status
   return status;
}
