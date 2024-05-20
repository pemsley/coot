
#include <string>
#include <vector>

class pae_t {
   std::string file_to_string(const std::string &file_name) const;
   //! return a png as a string
   std::string make_image(const std::vector<std::vector<int> > &pae_vecs) const;
   float get_max_value(const std::vector<std::vector<int> > &pae_vecs) const;
   public:
   pae_t(const std::string &file_name, int n_pixels);
   std::string image; // a png as a string
   int n_pixels;
   std::string get_image() const { return image; }
};
