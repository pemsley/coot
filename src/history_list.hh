
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>

namespace coot {

   class history_list_t {
   public:
      std::vector<std::vector<std::string> > history_strings;
      std::vector<std::vector<std::string> > history_list() const {
	 return history_strings;
      }
      void add_history_command(const std::vector<std::string> &command) {
	 history_strings.push_back(command);
      }
   };

   class command_history_t {
   public:
      std::vector<std::string> commands;
      int index;
      std::string command_history_file_name;
      command_history_t() : index(0), command_history_file_name(std::string(".coot_python_commands")) {
         read_history();
      }

      void add_to_history(const std::string &s) {
         if (! s.empty())
            commands.push_back(s);
         index = commands.size();
      }

      std::string get_previous_command() {
         index -= 1;
         if (index < 0) {
            index = 0;
            return "";
         } else {
            return commands[index];
         }
      }

      std::string get_previous_command_starting_with(const std::string &search_string) {

         if (index <=0 ) return search_string;

         std::size_t l1 = search_string.size();
         for (int idx=index; idx>=0; idx--) {
            const auto &c = commands[idx];
            std::cout << "looking for " << search_string << " index " << index
                      << " in " << commands.size() << " " << c << std::endl;
            std::size_t l2 = c.size();
            if (l2 > l1) {
               if (c.substr(0,l1) == search_string) {
                  return c;
               }
            }
         }
         return search_string;
      }

      std::string get_next_command() {
         if (commands.size() == 0)
            return "";
         index += 1;
         int idx_max = commands.size() - 1;
         if (index > idx_max)
            index = idx_max;
         return commands[index];
      }

      std::vector<std::string> unique_commands() const {
         std::vector<std::string> v;
         for (unsigned int i=0; i<commands.size(); i++) {
            const auto &c = commands[i];
            std::vector<std::string>::const_iterator it = std::find(v.begin(), v.end(), c);
            if (it == v.end())
               v.push_back(c);
         }
         return v;
      }

      void write_history() {
         if (! commands.empty()) {
            std::vector<std::string> all_history;
            std::ifstream fin(command_history_file_name);
            std::string line;
            while (std::getline(fin, line)) {
               all_history.push_back(line);
            }
            fin.close();
            all_history.insert(all_history.end(), commands.begin(), commands.end());
            std::ofstream f(command_history_file_name);
            std::vector<std::string>::iterator it;
            for (it=commands.begin(); it!=commands.end(); ++it) {
               const std::string &hs = *it;
               std::vector<std::string>::iterator it_next = it + 1;
               if (std::find(it_next, commands.end(), hs) == commands.end()) {
                  f << hs << "\n";
               }
            }
            f.close();
         }
      }

      void read_history() {
         std::ifstream fin(command_history_file_name);
         std::string line;
         while (std::getline(fin, line)) {
            commands.push_back(line);
         }
         index = commands.size();
      }
   };

}
