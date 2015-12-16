

namespace coot {

   // ------------ molecule probability scoring ------------
   class rama_score_t {
   public:
      rama_score_t() {
	 score = 0.0;
	 score_non_sec_str = 0.0;
	 n_zeros = 0;
      }
      // for all residues
      std::vector<std::pair<residue_spec_t, double> >  scores;
      // for non-Secondary structure residues
      std::vector<std::pair<residue_spec_t, double> >  scores_non_sec_str;
      double score;
      double score_non_sec_str;
      int n_residues() const { return scores.size(); }
      int n_residues_non_sec_str() const { return scores_non_sec_str.size(); }
      int n_zeros;
      std::vector<std::pair<residue_spec_t, int> > region;
   };

   // ==-------------- all molecule rotamer scoring --------------
   class rotamer_score_t {
   public:
      rotamer_score_t() {
	 score = 0.0;
	 n_pass = 0;
      } 
      std::vector<std::pair<residue_spec_t, double> > scores;
      double score;
      int n_pass; // GLY, PRO, ALA
      int n_rotamer_residues() const { return scores.size(); }
      void add (const residue_spec_t &rs, double p) {
	 std::pair<residue_spec_t, double> pair(rs,p);
	 scores.push_back(pair);
      }
   };
   

}
