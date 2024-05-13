
#include <fstream>
#include "pae.hh"

int main(int argc, char **argv) {
   int status = 0;
   int n_pixels = 600;
   std::string file_name = "AF-A0A192CA94-F1-predicted_aligned_error_v4.json";
   if (argc > 1) file_name = argv[1];
   pae_t pae(file_name, n_pixels);
   std::string s = pae.get_image();
   std::ofstream f("pae.png");
   f << s;
   f.close();
   return status;
}