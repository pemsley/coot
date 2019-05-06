
#include <iostream>
#include <logging.hh>

void logging::log(const std::string &type, const std::string &s) {

   log_item l(type, s);
   history.push_back(l);

}

void
logging::log(const std::string &s) {

   log_item l(s);
   history.push_back(l);
}


void
logging::show() const {

   for (std::size_t i=0; i<history.size(); i++) {
      const log_item &h = history[i];
      // maybe use localtime()
      std::cout << ctime(&h.t) << " " << h.type << " " << h.message << std::endl;
   }

}
