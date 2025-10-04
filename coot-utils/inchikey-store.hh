
#include <optional>
#include <string>
#include <utility>
#include <map>

namespace coot {

   class inchikey_store_t {

      void parse_from_data_components_file();

   public:
      inchikey_store_t() {}
      void init() { parse_from_data_components_file(); }
      std::map<std::string, std::pair<std::string, std::string> > store;
      std::optional<std::pair<std::string, std::string> > get_match(const std::string &key) const;
   };
}
