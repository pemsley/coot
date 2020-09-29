#ifndef PULSE_DATA_HH
#define PULSE_DATA_HH

class pulse_data_t {
public:
   int n_pulse_steps;
   int n_pulse_steps_max;
   pulse_data_t(int n1, int n2) {
      n_pulse_steps = n1;
      n_pulse_steps_max = n2;
   }
};


#endif
