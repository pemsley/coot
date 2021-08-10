
#ifndef FRESNEL_SETTINGS_T_HH
#define FRESNEL_SETTINGS_T_HH

#include <glm/glm.hpp>

class fresnel_settings_t {
public:
   bool state;
   float bias;
   float scale;
   float power;
   glm::vec4 colour;
   fresnel_settings_t(bool state_in, const float &f1, const float &f2, const float &f3) :
      state(state_in), bias(f1), scale(f2), power(f3), colour(glm::vec4(1,1,1,1)) {}
   fresnel_settings_t() : colour(glm::vec4(1,1,1,1)) {
      state = false;
      bias = 0.0;
      scale = 0.4;
      power = 3.5;
   }
   void update_settings(bool state_in, const float &f1, const float &f2, const float &f3) {
      state = state_in;
      bias = f1;
      scale = f2;
      power = f3;
   }
};


#endif // FRESNEL_SETTINGS_T_HH
