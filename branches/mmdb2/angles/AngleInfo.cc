 
#include <iostream>
#include <math.h>

#include "AngleInfo.h"


AngleInfo::AngleInfo() { 
   
   construct_theta_2_table(); 
   
   setup_angle_torsion_table(10,5);
   from_batched_angle_torsions(); 

} 


void
AngleInfo::assign_angle_torsion(double angle_bin_in, 
				double torsion_bin_in, 
				double n_count) {

   //
   // int theta_bin = int (round (( double (angle_bin_in) - angle_min)/step_ta_a));
   int theta_bin = int (floor ((angle_bin_in - 0)/step_ta_a + 0.5));

   int torsion_bin = int (floor ((torsion_bin_in - -180)/step_ta_t + 0.5 ));

//    cout << "debug: assigning " << theta_bin << " " << torsion_bin 
//  	<< " counts " << n_count << std::endl; 

   angle_torsion_table(theta_bin,torsion_bin) = float (n_count); 
//    std::cout << " and " 
// 	     << angle_torsion_table(theta_bin,torsion_bin) << endl; 

   // cout << "done" << endl; 
}

void
AngleInfo::from_batched_angle_torsions() { 

   from_batched_angle_torsions_bits_185();
   from_batched_angle_torsions_bits_190();
   from_batched_angle_torsions_bits_195();
   from_batched_angle_torsions_bits_200();
   from_batched_angle_torsions_bits_205();
   from_batched_angle_torsions_bits_210();
   from_batched_angle_torsions_bits_215();
   from_batched_angle_torsions_bits_220();
   from_batched_angle_torsions_bits_225();
   from_batched_angle_torsions_bits_230();
   from_batched_angle_torsions_bits_235();
   from_batched_angle_torsions_bits_240();
   from_batched_angle_torsions_bits_245();
   from_batched_angle_torsions_bits_250();
   from_batched_angle_torsions_bits_255();
   from_batched_angle_torsions_bits_260();
   from_batched_angle_torsions_bits_265();
   from_batched_angle_torsions_bits_270();
   from_batched_angle_torsions_bits_275();
   from_batched_angle_torsions_bits_280();
   from_batched_angle_torsions_bits_285();
   from_batched_angle_torsions_bits_290();
   from_batched_angle_torsions_bits_295();
   from_batched_angle_torsions_bits_300();
   from_batched_angle_torsions_bits_305();
   from_batched_angle_torsions_bits_310();
   from_batched_angle_torsions_bits_315();
   from_batched_angle_torsions_bits_320();
   from_batched_angle_torsions_bits_325();
   from_batched_angle_torsions_bits_330();
   from_batched_angle_torsions_bits_335();
   from_batched_angle_torsions_bits_340();
   from_batched_angle_torsions_bits_345();
   from_batched_angle_torsions_bits_350();
   from_batched_angle_torsions_bits_355();
   from_batched_angle_torsions_bits_360();
   from_batched_angle_torsions_bits_365();
   from_batched_angle_torsions_bits_370();
   from_batched_angle_torsions_bits_375();
   from_batched_angle_torsions_bits_380();
   from_batched_angle_torsions_bits_385();
   from_batched_angle_torsions_bits_390();
   from_batched_angle_torsions_bits_395();
   from_batched_angle_torsions_bits_400();
   from_batched_angle_torsions_bits_405();
   from_batched_angle_torsions_bits_410();
   from_batched_angle_torsions_bits_415();
   from_batched_angle_torsions_bits_420();
   from_batched_angle_torsions_bits_425();
   from_batched_angle_torsions_bits_430();
   from_batched_angle_torsions_bits_435();
   from_batched_angle_torsions_bits_440();
   from_batched_angle_torsions_bits_445();
   from_batched_angle_torsions_bits_450();
   from_batched_angle_torsions_bits_455();
   from_batched_angle_torsions_bits_460();
   from_batched_angle_torsions_bits_465();
   from_batched_angle_torsions_bits_470();
   from_batched_angle_torsions_bits_475();
   from_batched_angle_torsions_bits_480();
   from_batched_angle_torsions_bits_485();
   from_batched_angle_torsions_bits_490();
   from_batched_angle_torsions_bits_495();
   from_batched_angle_torsions_bits_500();
   from_batched_angle_torsions_bits_505();
   from_batched_angle_torsions_bits_510();
   from_batched_angle_torsions_bits_515();
   from_batched_angle_torsions_bits_520();
   from_batched_angle_torsions_bits_525();
   from_batched_angle_torsions_bits_530();
   from_batched_angle_torsions_bits_535();
   from_batched_angle_torsions_bits_540();


} 

void 
AngleInfo::setup_angle_torsion_table(float step_torsion_in,
				     float step_angle_in) { 

   step_ta_a = step_angle_in; 
   step_ta_t = step_torsion_in; 

//    angle_min = angle_min_in;
//    angle_max = angle_max_in;

//    torsion_min = torsion_min_in;
//    torsion_max = torsion_max_in;

   int angle_resize_size = nint((180 - 0)/step_ta_a);
   int torsion_resize_size = nint((180 - -180)/step_ta_t);

//    std::cout << "angle   resize size: " << resize_size << endl;
//    std::cout << "torsion resize size: " << torsion_resize_size << endl;


   std::cout << "angle_resize_size:  " << angle_resize_size << std::endl; 
   std::cout << "torsion_resize_size " << torsion_resize_size << std::endl; 

   if ( ! (angle_resize_size == torsion_resize_size) ) { 

      std::cout << "ERROR angle_resize_size and torsion_resize_size " 
	   << "must be equal" << std::endl; 
      
   } else { 

      angle_torsion_table.resize(torsion_resize_size); 
      angle_torsion_table_size_ = torsion_resize_size; 
      
   }
   

}

void
AngleInfo::print_angle_torsion_table() const {

   // FIXME

   for (float angle = 0; angle <= 180; angle += 2) { 
      for (float torsion = -180; torsion <= 180; torsion += 4) { 

	 std::cout << angle << " " << torsion << " " 
	      << prob_angle_torsion(angle, torsion) << std::endl; 

      }
   }
}


float
AngleInfo::check_sum() const {

   float sum = 0.0; 

   for (float angle = 0; angle <= 180; angle += step_ta_t) { 
      for (float torsion = -180; torsion < 180; torsion += step_ta_t) { 

	 sum += prob_angle_torsion(angle, torsion); 
      } 
   }

   std::cout << "# angle_torsion_table sum is: " << sum << std::endl;
   return sum; 
}

void
AngleInfo::normalize_angle_torsion() {

   float sum = 0.0; 

   for (float angle = 0; angle <= 180; angle += step_ta_t) { 
      for (float torsion = -180; torsion < 180; torsion += step_ta_t) { 

	 sum += prob_angle_torsion(angle, torsion); 
      } 
   }

   for (float angle = 0; angle <= 180; angle += step_ta_t) { 
      for (float torsion = -180; torsion < 180; torsion += step_ta_t) { 

	 assign_angle_torsion(angle,torsion, 
			      prob_angle_torsion(angle,torsion)/sum); 
      } 
   }
   
}

// Return the probability of the nearest bin
// so that 89.5 will have the same return values as 90.0
// 
float
AngleInfo::prob_angle_torsion_by_bin(float angle, float torsion) const {
   //
   // int theta_bin = nint ((angle - angle_min)/step_ta_a);
   int theta_bin = nint ((angle)/step_ta_a);

   int torsion_bin = nint ((torsion - -180)/step_ta_t);

   std::cout << "debug: theta_bin: " << theta_bin << " torsion_bin: "
	     << torsion_bin << std::endl; 

   return angle_torsion_table.interp(theta_bin,torsion_bin); 
}

// As above, except we interpolate from the angles, not just the
// flat bin, 5 degrees wide.
// float
// AngleInfo::interp_prob_angle_torsion(float angle, float torsion) const {

// }
   

float 
AngleInfo::prob_angle_torsion(float angle, float torsion) const { 

   float angle_index = angle/step_ta_a; 
   float torsion_index = (torsion +180)/step_ta_t; 

//     cout << "evaluating: angle bin: " << angle_index 
//  	<< ", torsion bin: " << torsion_index << endl; 

   int a = int (floor(angle_index));
   int t = int (floor(torsion_index)); 
   
//    cout << "DEBUG: int angle bin " << a 
// 	<< " and int torsion bin " << t << " has value " 
// 	<< angle_torsion_table(a,t) << endl; 

//    cout << "DEBUG: 46, 62 " << angle_torsion_table(46,62) << endl; 

   return angle_torsion_table.interp(angle_index, torsion_index); 

} 

void 
AngleInfo::theta_2_table_assign_step(float step) { 

   theta_2_step = step;

   int nbins = int (floor(180/theta_2_step)) + 1; 
   
   // cout << "DEBUG: theta_2_table nbins: " << nbins << endl; 
   theta_2_array = new float[nbins]; 
   theta_2_array[nbins -1] = 0; // put the end point at zero.
                                // (it may not be otherwise assigned). 
}

void 
AngleInfo::assign_theta_2(float angle, float n_count) { 
   
   int i_bin = nint(angle/theta_2_step); 
//    cout << "DEBUG:  assigning bin "<< i_bin << " with " 
// 	<< n_count << " counts" << endl; 

   theta_2_array[i_bin] = n_count; 

}

// look up theta_2 and return n_count for the apropriate bin.
// Except that we are doing interpolation too.
// 
float 
AngleInfo::theta_2_score(float theta_2_in) const {

   float f = (floor (theta_2_in/theta_2_step));

   int bin = int (f); 

   float w1 = theta_2_in/theta_2_step - f;
   float w2 = 1.0 - w1; 

//    cout << theta_2_in << " corresponds to bin: " << bin << endl; 
//    cout << " combining: " << theta_2_array[bin] << " and "
// 	<< theta_2_array[bin+1] << " with " << w1 << " and " << w2 << endl; 

   return w2*theta_2_array[bin] + w1*theta_2_array[bin+1]; 

} 


void
a_init() { 
   AngleInfo ai;  
} 
