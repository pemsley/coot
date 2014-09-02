// -*-c++-*-

#ifndef ANGLE_INFO_H
#define ANGLE_INFO_H

#include "wrapped-ncube.h"

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

class AngleInfo { 

   MatrixSq2<float> angle_torsion_table; 

   std::vector <float> torsion; 

   float* theta_2_array; 

   int angle_torsion_table_size_; 

   void assign_angle_torsion(double angle_bin,
			     double torsion_bin,
			     double n_count); 

   float theta_2_step; 
   float step_ta_a, step_ta_t; //     torsion_angle_table_angle_step
                               // and torsion_angle_table_torsion_step

   float angle_min, torsion_min, angle_max, torsion_max;

   void 
      setup_angle_torsion_table(float step_torsion_in,
				float step_angle_in); 
   
   // simply theta_2 - a one dimentional thing.
   // 
   void construct_theta_2_table(); // 0 and 180-step are included. 

   void assign_theta_2(float angle, float n_count); // integers converted to float; 

   void theta_2_table_assign_step(float step); 


   void from_batched_angle_torsions_bits_185();
   void from_batched_angle_torsions_bits_190();
   void from_batched_angle_torsions_bits_195();
   void from_batched_angle_torsions_bits_200();
   void from_batched_angle_torsions_bits_205();
   void from_batched_angle_torsions_bits_210();
   void from_batched_angle_torsions_bits_215();
   void from_batched_angle_torsions_bits_220();
   void from_batched_angle_torsions_bits_225();
   void from_batched_angle_torsions_bits_230();
   void from_batched_angle_torsions_bits_235();
   void from_batched_angle_torsions_bits_240();
   void from_batched_angle_torsions_bits_245();
   void from_batched_angle_torsions_bits_250();
   void from_batched_angle_torsions_bits_255();
   void from_batched_angle_torsions_bits_260();
   void from_batched_angle_torsions_bits_265();
   void from_batched_angle_torsions_bits_270();
   void from_batched_angle_torsions_bits_275();
   void from_batched_angle_torsions_bits_280();
   void from_batched_angle_torsions_bits_285();
   void from_batched_angle_torsions_bits_290();
   void from_batched_angle_torsions_bits_295();
   void from_batched_angle_torsions_bits_300();
   void from_batched_angle_torsions_bits_305();
   void from_batched_angle_torsions_bits_310();
   void from_batched_angle_torsions_bits_315();
   void from_batched_angle_torsions_bits_320();
   void from_batched_angle_torsions_bits_325();
   void from_batched_angle_torsions_bits_330();
   void from_batched_angle_torsions_bits_335();
   void from_batched_angle_torsions_bits_340();
   void from_batched_angle_torsions_bits_345();
   void from_batched_angle_torsions_bits_350();
   void from_batched_angle_torsions_bits_355();
   void from_batched_angle_torsions_bits_360();
   void from_batched_angle_torsions_bits_365();
   void from_batched_angle_torsions_bits_370();
   void from_batched_angle_torsions_bits_375();
   void from_batched_angle_torsions_bits_380();
   void from_batched_angle_torsions_bits_385();
   void from_batched_angle_torsions_bits_390();
   void from_batched_angle_torsions_bits_395();
   void from_batched_angle_torsions_bits_400();
   void from_batched_angle_torsions_bits_405();
   void from_batched_angle_torsions_bits_410();
   void from_batched_angle_torsions_bits_415();
   void from_batched_angle_torsions_bits_420();
   void from_batched_angle_torsions_bits_425();
   void from_batched_angle_torsions_bits_430();
   void from_batched_angle_torsions_bits_435();
   void from_batched_angle_torsions_bits_440();
   void from_batched_angle_torsions_bits_445();
   void from_batched_angle_torsions_bits_450();
   void from_batched_angle_torsions_bits_455();
   void from_batched_angle_torsions_bits_460();
   void from_batched_angle_torsions_bits_465();
   void from_batched_angle_torsions_bits_470();
   void from_batched_angle_torsions_bits_475();
   void from_batched_angle_torsions_bits_480();
   void from_batched_angle_torsions_bits_485();
   void from_batched_angle_torsions_bits_490();
   void from_batched_angle_torsions_bits_495();
   void from_batched_angle_torsions_bits_500();
   void from_batched_angle_torsions_bits_505();
   void from_batched_angle_torsions_bits_510();
   void from_batched_angle_torsions_bits_515();
   void from_batched_angle_torsions_bits_520();
   void from_batched_angle_torsions_bits_525();
   void from_batched_angle_torsions_bits_530();
   void from_batched_angle_torsions_bits_535();
   void from_batched_angle_torsions_bits_540();

   int nint (float i) const { 
      return i < 0 ? int (i-0.5) : int (i+0.5); 
   }
   
 public:

   AngleInfo();

   int angle_torsion_table_size() { return angle_torsion_table_size_; }

   void normalize_angle_torsion(); // convert to probabilities.

   float check_sum() const; 

   // this should be renamed to from batched_angle_torsions()
   // 
   void from_batched_angle_torsions();  // created by the make-code.awk script

   float prob_angle_torsion(float angle, float torsion) const; 

   float prob_angle_torsion_by_bin(float angle, float torsion) const; 

   float prob_torsion_torsion(float tor1, float tor2) const; 


   // where we have the previous values and the previous previous values
   // 
   // i.e. we are positioning n, and we have values for n-1 and n-2
   // 
   float prob_angle_torsion_angle_torsion(float angle1, float tors1, 
					  float angle2, float tors2) const; 


   void print_angle_torsion_table() const; 

   float theta_2_score(float theta_2_in) const; 


};



// 
void a_init(); 

#endif // ANGLE_INFO_H
