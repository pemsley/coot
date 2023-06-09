#ifndef SVG_STORE_KEY_HH
#define SVG_STORE_KEY_HH

#include <string>

// 20230215-PE The key to this map *should* contain the imol, so have it's own class ideally (with a less-than operator).
// However, it won't be a problem frequently. Let's see how long it is before someone complains.
// (It's the multi-LIG problem).
//
// 2023040 Filo complained
//
class svg_store_key_t {
public:
   int imol;
   std::string comp_id;
   svg_store_key_t(int imol, const std::string &c) : imol(imol), comp_id(c) {}
   bool operator<(const svg_store_key_t &other) const {
      if (imol<other.imol)
         return true;
      else
         if (comp_id < other.comp_id)
            return true;
      return false;
   }
};
   


#endif // SVG_STORE_KEY_HH
