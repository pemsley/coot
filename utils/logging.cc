
#include <sys/time.h>

#include <iostream>
#include "logging.hh"

void
logging::log(logging::type_t type, const std::string &s) {

   // use:  add_log_item_to_history(log_time(type, s));
   log_item l(type, s);

   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   
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

