

class map_statistics_t {

public:
   map_statistics_t(const double &mean_in,
		    const double &sd_in,
		    const double &skew_in,
		    const double &kurtosis_in) {
      mean = mean_in;
      sd = sd_in;
      skew = skew_in;
      kurtosis = kurtosis_in;
   } 
   double mean;
   double sd;
   double skew;
   double kurtosis;
};

