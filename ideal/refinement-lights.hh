

#ifndef REFINEMENT_LIGHTS_HH
#define REFINEMENT_LIGHTS_HH

#include "rama-types.hh"

namespace coot {

   class refinement_lights_info_t {
   public:
      class the_worst_t {
      public:
	 the_worst_t(const unsigned int idx_in, const float &val) : restraints_index(idx_in), value(val) {
	    is_set = true;
	    restraints_index = -1;
	 }
	 the_worst_t() { is_set = false; value = -99999;}
	 int restraints_index; // can be -1
	 float value;
	 bool is_set;
	 // update if value_in is more than this value
	 void update_if_worse(const float &value_in, const int &idx_in) {
	    if (! is_set) {
	       restraints_index = idx_in;
	       value = value_in;
	       is_set = true;
	       // std::cout << "update_if_worse() setting worst to " << value << std::endl;
	    } else {
	       if (value_in > value) {
		  restraints_index = idx_in;
		  value = value_in;
		  // std::cout << "update_if_worse() updating worst to " << value << std::endl;
	       }
	    }
	 }
	 void update_if_worse(const the_worst_t &baddie_in) {
	    if (baddie_in.is_set) {
	       if (! is_set) {
		  restraints_index = baddie_in.restraints_index;
		  value = baddie_in.value;
		  is_set = true;
	       } else {
		  if (baddie_in.value > value) {
		     restraints_index = baddie_in.restraints_index;
		     value = baddie_in.value;
		  }
	       }
	    }
	 }
      };
      std::string name;   // e.g. "Bonds" or "Angles"
      std::string label;  // e.g. "Bonds:  6.543" 
      float value;        // e.g. 6.543
      int rama_type;
      the_worst_t worst_baddie;
      refinement_lights_info_t(const std::string &name_in, const std::string label_in, float value_in) {
	 name = name_in;
	 label = label_in;
	 value = value_in;
	 rama_type = RAMA_TYPE_LOGRAMA;
      }
   };

}

#endif // REFINEMENT_LIGHTS_HH
