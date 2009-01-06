
# awk -f gdk-keysym.awk /usr/include/gtk-2.0/gdk/gdkkeysyms.h 

/#define GDK_/ { 
  split($2, arr, "GDK_");
  split($3, bar, "0x");
  symbol = arr[2];
  binhex_str = bar[2];
  dec = strtonum($3)
  # print "(list \"" arr[2] "\"", dec")"
  print "   a.push_back(std::pair<std::string,int>(\"" arr[2] "\",", dec "));"
}

# don't forget these
#   a.push_back(std::pair<std::string,int>(":", 58));
#   a.push_back(std::pair<std::string,int>(";", 59));
#   a.push_back(std::pair<std::string,int>("<", 60));
#   a.push_back(std::pair<std::string,int>("=", 61));
#   a.push_back(std::pair<std::string,int>(">", 62));
#   a.push_back(std::pair<std::string,int>("?", 63));
#   a.push_back(std::pair<std::string,int>("@", 64));
#   a.push_back(std::pair<std::string,int>("!", 33));
#   a.push_back(std::pair<std::string,int>("$", 36));
#   a.push_back(std::pair<std::string,int>("%", 37));
#   a.push_back(std::pair<std::string,int>("&", 38));
#   a.push_back(std::pair<std::string,int>("*", 42));
#   a.push_back(std::pair<std::string,int>("+", 43));
#   a.push_back(std::pair<std::string,int>(",", 44));
#   a.push_back(std::pair<std::string,int>("-", 45));
#   a.push_back(std::pair<std::string,int>(".", 46));
#   a.push_back(std::pair<std::string,int>("/", 47));
#   a.push_back(std::pair<std::string,int>("(", 91));
#   a.push_back(std::pair<std::string,int>(")", 93));
#   a.push_back(std::pair<std::string,int>("_", 95));
#   a.push_back(std::pair<std::string,int>("|", 124));
#   a.push_back(std::pair<std::string,int>("~", 126));

 