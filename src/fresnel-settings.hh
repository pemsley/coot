
#ifndef FRESNEL_SETTINGS_T_HH
#define FRESNEL_SETTINGS_T_HH

class fresnel_settings_t {
public:
   bool state;
   float bias;
   float scale;
   float power;
   fresnel_settings_t(bool state_in, const float &f1, const float &f2, const float &f3) :
      state(state_in), bias(f1), scale(f2), power(f3) {}
   fresnel_settings_t() {
      state = false;
      bias = 0.0;
      scale = 1.0;
      power = 1.5;
   }
};


#endif // FRESNEL_SETTINGS_T_HH
