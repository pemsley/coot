
#ifndef LBG_SHARED_HH
#define LBG_SHARED_HH

namespace coot {


   class bash_distance_t {
   public:
      double dist;
      bool limited;
      bash_distance_t() {
	 limited = 0;
	 dist = -1;
      }
      bash_distance_t(double d) {
	 limited = 1;
	 dist = d;
      }
      bool unlimited() const {
	 return !limited;
      } 
      friend std::ostream& operator<< (std::ostream& s, const bash_distance_t &bd);
   };
   std::ostream& operator<< (std::ostream& s, const bash_distance_t &bd);
}

#endif // LBG_SHARED_HH
