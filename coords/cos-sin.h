
#ifndef COS_SIN_H
#define COS_SIN_H

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

class cos_sin { 

   static int cos_to_sine_table_steps;
   static int is_table_filled; 

   //
   static float *cos_to_sine_table; 

 public:
   
   void fillTable(int nsteps); 
   float operator()(float) const; 

   cos_sin();
   cos_sin(int); 
   ~cos_sin();

   void check_table(void) const; 

};

#endif // COS_SIN_H   
 

