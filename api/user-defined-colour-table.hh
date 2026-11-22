
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

   //! @param f the fraction (typically 0 to 1) along the colour ramp
   //! @return the colour linearly interpolated between the stops that bracket f.
   //!         Fractions before the first stop or after the last stop take the
   //!         colour of that end stop. If the table is not set, return dark grey.
   glm::vec4 fraction_to_colour(float f) const {
      if (! is_set()) return glm::vec4(0.1f, 0.1f, 0.1f, 1.0f);
      const colour_pair_t &first = colour_table.front();
      const colour_pair_t &last  = colour_table.back();
      if (f <= first.colour_frac) return glm::vec4(first.colour, 1.0f);
      if (f >= last.colour_frac)  return glm::vec4(last.colour,  1.0f);
      for (unsigned int i=0; i<(colour_table.size()-1); i++) {
         const colour_pair_t &lower_point = colour_table[i];
         const colour_pair_t &upper_point = colour_table[i+1];
         if (f <= upper_point.colour_frac) {
            float range = upper_point.colour_frac - lower_point.colour_frac;
            float fraction = (range > 0.0f) ? (f - lower_point.colour_frac) / range : 0.0f;
            return glm::vec4(glm::mix(lower_point.colour, upper_point.colour, fraction), 1.0f);
         }
      }
      return glm::vec4(last.colour, 1.0f); // not reached
   }

};


#endif // API_USER_DEFINED_COLOUR_TABLE_HH

