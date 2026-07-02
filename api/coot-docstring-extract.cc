
#include <sstream>
#include <vector>
#include <iostream>

#include "coot-utils/pugixml.hpp"

#include "coot-docstring-extract.hh"

std::unordered_map<std::string, std::string>
coot::parse_docstrings_from_doxygen_xml(const std::string &xml_file_name) {

   std::unordered_map<std::string, std::string> docstring_cache;

   auto convert_type = [] (const std::string &s_in) {
      std::string s = s_in;
      if (s_in == "const std::string &") s = "str";
      if (s_in == "std::string")         s = "str";
      if (s_in == "void")                s = "None";
      if (s_in == "std::vector<")        s = "list";
      if (s_in == "std::vector< std::pair< double, double > >") s = "list";
      if (s_in == "std::vector< std::pair< std::string, std::string > >") s = "list";
      if (s_in.find("std::vector<") != std::string::npos) s = "list";
      if (s_in.compare(0,10, "std::pair<") == 0) s= "tuple"; // needs checking. Use starts_with() in C++20
      return s;
   };

   class arg_info_t {
   public:
      arg_info_t(const std::string &n, const std::string &t) : name(n), type(t) {}
      std::string name;
      std::string type;
      std::string description;
   };

   auto update_arg_in_args = [] (const std::string &arg_name, const std::string &descr,
                                 std::vector<arg_info_t> &args) {

      for(auto &arg : args) {
         if (arg.name == arg_name) {
            arg.description = descr;
         }
      }
   };

   pugi::xml_document doc;
   if (!doc.load_file(xml_file_name.c_str())) {
      std::cout << "WARNING:: doxygen file " << xml_file_name
                << " not found - nanobind API docummentation will not be generated"
                << std::endl;
      return docstring_cache;
   }
   auto compounddef = doc.child("doxygen").child("compounddef");
   for (auto sectiondef : compounddef.children("sectiondef")) {
      for (auto member : sectiondef.children("memberdef")) {
         auto name_elem = member.child("name");
         if (!name_elem) continue;
         std::string name = name_elem.child_value();
         std::ostringstream oss;
         std::vector<arg_info_t> args;

         // Collect all <para> from <briefdescription>
         auto brief = member.child("briefdescription");
         if (brief) {
            for (auto para : brief.children("para")) {
               std::string para_text = para.text().get();
               if (!para_text.empty())
                  oss << para_text << "\n";
            }
         }

         auto type = member.child("type");
         std::string type_string;
         if (type)
            type_string = convert_type(type.text().get());

         // can be many params
         for (auto param : member.children("param")) {
            auto p_type = param.child("type");
            auto p_declname = param.child("declname");
            if (p_type) {
               if (p_declname) {
                  std::string tt = convert_type(p_type.text().get());
                  arg_info_t ai(p_declname.text().get(), tt);
                  args.push_back(ai);
               }
            }
         }

         std::string return_type_docs;

         // Collect all <para> from <detaileddescription>
         auto detailed = member.child("detaileddescription");
         if (detailed) {
            unsigned int n_para = 0;
            for (auto para : detailed.children("para")) {
               n_para++;
               std::string para_text = para.text().get();
               if (!para_text.empty()) {
                  if (n_para > 1)
                     oss << "\n    ";
                  if (para_text[0] == '\n')
                     para_text.erase(0,1); // remove first char
                  oss << para_text << "\n";
               }
               for (auto parameterlist : para.children("parameterlist")) {
                  for (auto parameteritem : parameterlist.children("parameteritem")) {
                     std::string parameter_name_text;
                     for (auto parameternamelist : parameteritem.children("parameternamelist")) {
                        for (auto parametername : parameternamelist.children("parametername")) {
                           parameter_name_text = parametername.text().get();
                        }
                     }
                     for (auto parameterdescription : parameteritem.children("parameterdescription")) {
                        for (auto d_para : parameterdescription.children("para")) {
                           std::string t =  d_para.text().get();
                           if (! parameter_name_text.empty()) {
                              if (! t.empty()) {
                                 update_arg_in_args(parameter_name_text, t, args); // modify an arg in args
                              }
                           }
                        }
                     }
                  }
               }
               for (auto simplesect : para.children("simplesect")) {
                  if (std::string(simplesect.attribute("kind").value()) == "return") {
                     for(auto ss_para : simplesect.children("para")) {
                        std::string return_type_doc = ss_para.text().get();
                        return_type_docs += return_type_doc;
                     }
                  }
               }
            }
            if (! args.empty()) {
               oss << "\n";
               oss << "    Args:\n";
               for (const auto &arg : args) {
                  oss << "        " << arg.name << " (" << arg.type << "): " << arg.description << "\n";
               }
            }
            if (! type_string.empty()) {
               oss << "\n";
               oss << "    Returns:\n";
               if (return_type_docs.empty())
                  oss << "        " << type_string << "\n";
               else
                  oss << "        " << type_string << ": " << return_type_docs << "\n";
            }
         }
         docstring_cache[name] = oss.str();
      }
   }
   return docstring_cache;
}
