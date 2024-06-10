
#include <pwd.h>

#include <cstdlib>
#include <string>
#include <filesystem>

namespace coot {
   class xdg_t {
      std::filesystem::path   data_home;
      std::filesystem::path  state_home;
      std::filesystem::path config_home;
      std::filesystem::path runtime_dir;

      std::filesystem::path get_home_dir() {
         std::string home;
#ifdef WINDOWS
#else
         struct passwd *pw = getpwuid(getuid());
         const char *home_str = pw->pw_dir;
         if (home_str)
            home = std::string(home_str);
#endif
         return std::filesystem::path(home);
      }
      public:
      xdg_t() {
         char *e;
         e = std::getenv("XDG_DATA_HOME");   if (e)   data_home = e;
         e = std::getenv("XDG_STATE_HOME");  if (e)  state_home = e;
         e = std::getenv("XDG_CONFIG_HOME"); if (e) config_home = e;
         e = std::getenv("XDG_RUNTIME_DIR"); if (e) runtime_dir = e;
         if (data_home.empty()) {
            std::filesystem::path d = get_home_dir();
            d.append(".local");
            d.append("share");
            d.append("Coot");
            data_home = d;
         }
         if (config_home.empty()) {
            std::filesystem::path d = get_home_dir();
            d.append(".config");
            d.append("Coot");
            config_home = d;
         }
         if (state_home.empty()) {
            std::filesystem::path d = get_home_dir();
            d.append(".local");
            d.append("state");
            d.append("Coot");
            state_home = d;
         }
      }
      std::filesystem::path get_state_home() const {
         if (!std::filesystem::is_directory(state_home))
            std::filesystem::create_directories(state_home);
         return state_home;
      }
      std::filesystem::path get_data_home() const {
         if (!std::filesystem::is_directory(data_home))
            std::filesystem::create_directories(data_home);
         return data_home;
      }
      std::filesystem::path get_config_home() const {
         if (!std::filesystem::is_directory(config_home))
            std::filesystem::create_directories(config_home);
         return config_home;
      }
   };
}
