
#include <string>
#include <vector>
#include <map>

namespace coot {

   // let's say data were collected on a human sample, with yeast structure known.
   //
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
      void store(const std::vector<std::pair<int, std::string> > &seqs);
      int description_to_resno_start(const std::string &descr) const;
   public:
      int resno_start;
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

}
