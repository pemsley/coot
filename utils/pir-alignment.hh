
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
	 matched_residue_t(char aligned_in, char target_in, int res_no_in) : res_no(res_no_in), target(target_in), aligned(aligned_in) { }
	 matched_residue_t() {};
	 int res_no;
	 std::string ins_code; // using this is perverse
	 char target;
	 char aligned; // aligned is aligned to target
	 // friend std::ostream &operator<< (std::ostream &s, const matched_residue_t &m);
      };
   private:
      std::ostream &operator<< (const matched_residue_t &m);
      void init(const std::string &s);
      bool is_pir_aa(char a, bool allow_gaps=true) const;
      void store(const std::vector<std::pair<int, std::string> > &seqs);
      int description_to_resno_start(const std::string &descr) const;
   public:
      pir_alignment_t(const std::string &s);
      void read_file(const std::string &file_name);
      pir_alignment_t();
      // any number of alignments can be read - most typically there will be
      // just 1, I imagine
      std::vector<std::map<int, matched_residue_t> > matches;
   };
   std::ostream &operator<< (std::ostream &s, const pir_alignment_t::matched_residue_t &m);

}
