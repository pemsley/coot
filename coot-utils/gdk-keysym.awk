
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

