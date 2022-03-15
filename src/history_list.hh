
#include <vector>
#include <string>
#include <fstream>

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
      command_history_t() {
         index = 0;
         command_history_file_name = std::string(".coot_python_commands");
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

      std::string get_next_command() {
         if (commands.size() == 0)
            return "";
         index += 1;
         int idx_max = commands.size() - 1;
         if (index > idx_max)
            index = idx_max;
         return commands[index];
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
            for (it=commands.begin(); it!=commands.end(); it++) {
               const std::string &hs = *it;
               std::vector<std::string>::iterator it_next = it + 1;
               std::find(it_next, commands.end(), hs);
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
