
#include "pugixml.hpp"
#include <iostream>

#include "parse-pisa.hh"
#include "pugiconfig.hpp"

int parse_pisa(const std::string &file_name) {

   pugi::xml_document doc;
   pugi::xml_parse_result result = doc.load_file(file_name.c_str());
   if (!result)
      return -1;
        
   for (pugi::xml_node tool: doc.child("Profile").child("Tools").children("Tool")) {
      int timeout = tool.attribute("Timeout").as_int();
      if (timeout > 0)
         std::cout << "Tool " << tool.attribute("Filename").value() << " has timeout " << timeout << "\n";
   }
   return 0;
}


