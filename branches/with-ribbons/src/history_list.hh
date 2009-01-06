
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

}
