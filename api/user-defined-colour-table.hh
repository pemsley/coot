
#ifndef API_USER_DEFINED_COLOUR_TABLE_HH
#define API_USER_DEFINED_COLOUR_TABLE_HH


#include <vector>
#include <algorithm>
#include <iostream> // needed for testing

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/glm.hpp>

//! a colour table with colour stops and
//! linear interpolation - there is
//! no colour rotation here.
class user_defined_colour_table_t {

public:
   class colour_pair_t {
   public:
      float colour_frac;
      glm::vec3 colour;
      colour_pair_t(const float &fc, const glm::vec3 &c) : colour_frac(fc), colour(c) {}
   };

private:
   static bool sorter(const colour_pair_t &p1, const colour_pair_t &p2) {
      return p1.colour_frac < p2.colour_frac;
   }

   void sort_colour_table() {
      std::sort(colour_table.begin(), colour_table.end(), sorter);
   }

public:

   std::vector<colour_pair_t> colour_table;
   bool is_set() const { return (colour_table.size() > 1); }
   void add_stop(float f, const glm::vec3 &c) {
      colour_pair_t p(f, c);
      colour_table.push_back(p);
      sort_colour_table();

#if 0
      for (const auto &item : colour_table)
	 std::cout << "   " << item.colour_frac << " " << glm::to_string(item.colour) << std::endl;
#endif

   }
   void clear() {
      colour_table.clear();
   }

};


#endif // API_USER_DEFINED_COLOUR_TABLE_HH

