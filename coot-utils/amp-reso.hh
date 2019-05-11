
#ifndef AMP_RESO_HH
#define AMP_RESO_HH

namespace coot {

   class amplitude_vs_resolution_point {
   public:
      amplitude_vs_resolution_point() {
	 sum = 0.0; count = 0; resolution_recip_sum = 0.0; finished = false; }
      double sum;
      double average;
      unsigned int count;
      double resolution_recip_sum;
      double resolution_recip;
      bool finished;
      void add(const float &f_in, const float &inv_res_sq_in) {
	 sum += f_in;
	 resolution_recip_sum += inv_res_sq_in;
	 count += 1;
      }
      void finish() {
	 if (count > 0) {
	    average = sum/static_cast<float>(count);
	    resolution_recip = resolution_recip_sum/static_cast<float>(count);
	 }
	 finished = true;
      }
      double get_average_fsqrd() const {
	 if (finished) {
	    return average;
	 } else {
	    std::cout << "amplitude_vs_resolution_point() Not finihsed " << std::endl;
	    return 0.0;
	 }
      }
      double get_invresolsq() const {
	 if (finished) {
	    return resolution_recip;
	 } else {
	    std::cout << "amplitude_vs_resolution_point() Not finihsed " << std::endl;
	    return 0.0;
	 }
      }
   };
}

#endif // AMP_RESO_HH

