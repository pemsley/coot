
#include <string>
#include <vector>
#include <map>

namespace coot {

   class pir_t {
   public:
      int resno_start;
      int resno_end;
      std::string seq;
      pir_t(int resno_start_in, int resno_end_in, const std::string &seq_in) :
	 resno_start(resno_start_in), resno_end(resno_end_in), seq(seq_in) {}
   };

   class pir_alignment_t {
   public:
      class matched_residue_t {
      public:
	 matched_residue_t(char aligned_in, char target_in) : target(target_in), aligned(aligned_in) { }
	 matched_residue_t() {};
	 char target;
	 char aligned; // "aligned" is aligned to target/template
	 friend std::ostream &operator<< (std::ostream &s, const matched_residue_t &m);
      };
   private:
      std::ostream &operator<< (const matched_residue_t &m);
      void init(const std::string &s);
      bool is_pir_aa(char a, bool allow_gaps=true) const;
      // void store(const std::vector<std::pair<int, std::string> > &seqs);
      void store(const std::vector<pir_t> &seqs);

      // split the input string descr into fields based on ":" separator.
      // return the field_number field as an int.
      //
      std::pair<bool, int>
      description_split(const std::string &descr, int field_numpber) const;
      int description_to_resno_start(const std::string &descr) const;
   public:
      int resno_start; // the starting residue number of the alignement in the model (aligned)
      int resno_end;   // the end residue number of the alignement in the model (aligned)
      int resno_start_structure; // the starting residue number should change to this
                                 // it the target sequence (it's "our" structure - the one to which
                                 // we are mutating).
      int resno_end_structure;  // the last residue number in the range (combined with resno_end_structure)
                                // if is -1, then no "resno_end_structure" was specified in pir alignment file.
      pir_alignment_t(const std::string &s);
      void read_file(const std::string &file_name);
      pir_alignment_t();
      // any number of alignments can be read - most typically there will be
      // just 1, I imagine
      std::vector<std::vector<matched_residue_t> > matches;
      unsigned int size() const { return matches.size(); }
      unsigned int size(unsigned int idx) const { return matches[idx].size(); }
      const std::vector<matched_residue_t> &get_matches(const int i) const { return matches[i]; }
   };
   std::ostream &operator<< (std::ostream &s, const pir_alignment_t::matched_residue_t &m);

}
