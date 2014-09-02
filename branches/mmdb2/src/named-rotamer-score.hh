#ifndef NAMED_ROTAMER_SCORE_HH
#define NAMED_ROTAMER_SCORE_HH

#include <string>

namespace coot {

   class named_rotamer_score {
   public:
      std::string name;
      float clash_score;
      float rotamer_probability_score;
      std::vector<std::pair<std::string, float> > density_score_for_atoms;
      float density_fit_score;
      named_rotamer_score(const std::string &name_in,
			  float rotamer_probability_score_in,
			  float clash_score_in,
			  const std::vector<std::pair<std::string, float> > &atom_density_in,
			  float density_fit_score_in) {
	 name = name_in;
	 clash_score = clash_score_in;
	 density_fit_score = density_fit_score_in;
	 rotamer_probability_score = rotamer_probability_score_in;
	 density_score_for_atoms = atom_density_in;
      } 
   };
} 

#endif // NAMED_ROTAMER_SCORE_HH
