

namespace coot {

   namespace stats {

      double cached_kurtosis;
      bool have_cached_kurtosis;

      // 1-d data
      class single {
      public:
	 std::vector<double> v;

	 single() {
	    have_cached_kurtosis = false;
	 }
	 single(const std::vector<double> &v_in) {
	    v = v_in;
	    have_cached_kurtosis = false;
	 }
	 unsigned int size() const { return v.size(); }
	 void add(const double &a) {
	    v.push_back(a);
	    have_cached_kurtosis = false;
	 }
	 
	 double mean() const {
	    double m = 0;
	    double sum = 0;
	    if (v.size() ) { 
	       for (unsigned int i=0; i<v.size(); i++)
		  sum += v[i];
	       m = sum/double(v.size());
	    }
	    return m;
	 }

	 double variance() const {
	    double var = 0;
	    double sum = 0;
	    double sum_sq = 0;
	    if (v.size() ) { 
	       for (unsigned int i=0; i<v.size(); i++) { 
		  sum += v[i];
		  sum_sq += v[i] * v[i];
	       }
	       double m = sum/double(v.size());
	       var = sum_sq/double(v.size()) - m*m;
	    }
	    if (var < 0) var = 0; // numerical stability
	    return var;
	 }

	 double kurtosis() const {

	    // recall kurtosis, $k$ of $N$ observations:
	    // k = \frac{\Sigma(x_i - \mu)^4} {N \sigma^4} - 3    
	    // (x_i - \mu)^4 = x_i^4 + 4x_i^3\mu + 6x_i^2\mu^2 + 4x_i\mu^3 + \mu^4

	    // Can't enable this! A compiler bug (maybe) is apparent when sorting
	    // (bonds_vec_k_sorter()).
	    // g++ (Ubuntu 4.4.3-4ubuntu5.1) 4.4.3
	    // 
	    // if (have_cached_kurtosis)
	    // return cached_kurtosis;
	    
	    double k = -999;
	    if (v.size() ) {

	       double sum_to_the_4 = 0;
	       double m = mean();
	       double var = variance();

	       if (var > 0) {
		  double s = sqrt(var);
		  for (unsigned int i=0; i<v.size(); i++) { 
		     double t = v[i] - m;
		     sum_to_the_4 += t * t * t * t;
		  }
		  k = sum_to_the_4/(double(v.size()) * var * var);
		  cached_kurtosis = k;
		  have_cached_kurtosis = true;
	       }
	    }
	    return k;
	 }
      };
   }
}


