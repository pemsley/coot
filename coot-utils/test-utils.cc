
#include <iostream>
#include "coot-utils.hh"

void split_test(const std::string &r) {

   std::vector<std::string> v = coot::util::split_string(r, " ");
   std::cout << ":" << r << ": -> ";
   for (unsigned int i=0; i<v.size(); i++) {
      std::cout << ":" << v[i] << ": ";
   }
   std::cout << "\n";
}

int main(int argv, char **argc) {


   std::string s("Cottage");
   std::string r("tag"); 
   std::cout  << "removing :" << r << ": from :" << s << ": gives :"
	      << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "tage";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
	     << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "e";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
	     << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
	     << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "ball";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
	     << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "Cottage";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
	     << coot::util::remove_string(s, r) <<  ":" << std::endl;

   r = "Columns";
   split_test(r);
   r = "Columns ";
   split_test(r);
   r = " Columns ";
   split_test(r);
   r = "Columns     of   letters  ";
   split_test(r);
   r = " Columns     of   letters  ";
   split_test(r);
   
   return 0; 
}

