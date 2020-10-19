
#ifndef LOGGING_HH
#define LOGGING_HH

#include <iostream>
#include <ctime>
#include <string>
#include <vector>


class logging {

   class log_item {
   public:
      log_item(const std::string &type_in,
	       const std::string &message_in) : t(0), type(type_in), message(message_in) {}
      explicit log_item(const std::string &message_in) : t(0), type("plain"), message(message_in) {}
      time_t t;
      std::string type;
      std::string message;
      friend std::ostream& operator<<(std::ostream &o, const log_item &li);
   };

   std::vector<log_item> history;
   void operator<<(const std::string &s);
 public:
   logging() {}
   void log(const std::string &s);
   void log(const std::string &type, const std::string &s);
   void show() const;
   friend std::ostream& operator<<(std::ostream &o, const log_item &li);
};

class logging2 : public std::ostream {
public:
   std::vector<std::string> history;
   logging2() {}
   logging2& operator<<(const float &li);
   logging2& operator<<(const std::string &li);
};

logging2 coot_log;

#endif // LOGGING_HH
