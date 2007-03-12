
#include "coot-utils.hh"
#include "coot-shelx.hh"
#include <iostream>

int
main(int argc, char **argv) {

   /*
   std::string s;
   s = "AVasdfasdfC";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
	     << std::endl;
   s = "   AVC";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
	     << std::endl;
   s = " AVC ";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
	     << std::endl;
   s = "C";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
	     << std::endl;
   s = "";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
	     << std::endl;
   */

    std::string a("ABCDefgh");
    std::cout << a << " downcased: " << coot::util::downcase(a) << std::endl;
    std::cout << a << "   upcased: " << coot::util::upcase(a) << std::endl;

   if (argc > 2) {
      coot::ShelxIns sh;
      coot::shelx_read_file_info_t p = sh.read_file(argv[1]);
      sh.write_ins_file(p.mol, std::string(argv[2]));
   } else {
      std::cout << "Usage: " << argv[0] << " shelx-ins-file-name out-file-name"
		<< std::endl;
   }
   return 0;
}
