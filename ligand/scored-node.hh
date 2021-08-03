#ifndef SCORED_NODE_HH
#define SCORED_NODE_HH

#include <utility>
#include <iostream>

namespace coot {
   class scored_node_t {
   public:
      unsigned int atom_idx;
      double spin_score;
      // what (best) score do we get if the peptide is in the other direction?
      std::pair<bool, double> reverse_spin_score;  // added in spin_score_pairs(), the bool marks if this is a sensible value (of course)
      double alpha; // angle at which spin_score is recorded
      bool reversed_flag; // what's this for?
      bool udd_flag;
      std::string name; // for debugging
      scored_node_t(unsigned int idx_in, double score_in, const double &angle_in) : reverse_spin_score(std::make_pair(false, 0)),
                                                                                    alpha(angle_in) {
	 atom_idx = idx_in;
	 spin_score = score_in;
	 reversed_flag = false;
	 udd_flag = false;
      }
      scored_node_t(unsigned int idx_in, double score_in, const double &angle_in, const std::string &name_in) :
         reverse_spin_score(std::make_pair(false, 0)), alpha(angle_in), name(name_in) {
	 atom_idx = idx_in;
	 spin_score = score_in;
	 reversed_flag = false;
	 udd_flag = false;
      }
      scored_node_t(unsigned int idx_in, double score_in, const double &angle_in,
		    bool reversed_flag_in, const std::string &name_in) : reverse_spin_score(std::make_pair(false, 0)),
                                                                         alpha(angle_in),
                                                                         name(name_in) { 
	 atom_idx = idx_in;
	 spin_score = score_in;
	 reversed_flag = reversed_flag_in;
	 udd_flag = false;
      }
      scored_node_t(const scored_node_t &sn, bool reversed_flag_in) : reverse_spin_score(std::make_pair(false, 0)) {
	 atom_idx   = sn.atom_idx;
	 spin_score = sn.spin_score;
	 alpha      = sn.alpha;
	 reversed_flag = reversed_flag_in;
	 udd_flag = false;
      }
      scored_node_t() : reverse_spin_score(std::make_pair(false, 0)) {
	 atom_idx = 999999;
	 spin_score = -9999;
	 alpha = -1;
	 reversed_flag = false;
	 udd_flag = false;
      }
      bool operator==(const scored_node_t &other) const
      { return (other.atom_idx == atom_idx); }
      static bool sort_scores(const scored_node_t &s1, const scored_node_t &s2) {
	 return (s2.spin_score < s1.spin_score);
      }
      static bool sort_pair_scores(const std::pair<unsigned int, scored_node_t> &s1, const std::pair<unsigned int, scored_node_t> &s2) {
	 return (s2.second.spin_score < s1.second.spin_score);
      }
      friend std::ostream& operator<<(std::ostream& s, scored_node_t sn);
   };
   std::ostream& operator<<(std::ostream &s, scored_node_t sn);

}


#endif // SCORED_NODE_HH
