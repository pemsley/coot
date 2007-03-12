
#ifndef COOT_FASTA_HH
#define COOT_FASTA_HH

#ifndef HAVE_STRING
#include<string>
#define HAVE_STRING
#endif


namespace coot {

   class fasta {
   public:
      std::string name;
      std::string sequence;
      short int is_fasta_aa(const std::string &a) const;
      fasta(const std::string &combined_string); // decomposition happens in constructor
      fasta(const std::string &name_in, const std::string &seq);
      std::string format() const {
	 std::string s = "> ";
	 s += name;
	 s += "\n";
	 s += sequence;
	 return s;
      }
   };

} 

#endif // COOT_FASTA_HH
