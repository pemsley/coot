

#include <fstream>
#include <iostream>
#include <stdexcept>
#include "coot-utils.hh"
#include "pir-alignment.hh"


std::ostream &
coot::operator<< (std::ostream &s, const coot::pir_alignment_t::matched_residue_t &m) {
   s << m.target << " " << m.aligned;
   return s;
}

coot::pir_alignment_t::pir_alignment_t() {

   resno_start = -1; // unset really
   resno_start_structure = -1; // unset
   resno_end_structure = -1;
}

coot::pir_alignment_t::pir_alignment_t(const std::string &s) {

   init(s);
}

bool
coot::pir_alignment_t::is_pir_aa(char a, bool allow_gaps) const {

   bool r = false;

   if (a == 'A' || a == 'G' ) { 
      r = true;
   } else {
      if (   a == 'C' || a == 'D' || a == 'E' || a == 'F' || a == 'H' || a == 'I'
	  || a == 'K' || a == 'L' || a == 'M' || a == 'N' || a == 'P' || a == 'Q' 
	  || a == 'R' || a == 'S' || a == 'T' ||             a == 'V' || a == 'W' 
          || a == 'Y' || a == 'Z' || a == 'X'
	  || a == 'U' ) {
	 r = true;
      }
   }

   if (!r)
      if (allow_gaps)
	 if (a == '-')
	    r = true;

   return r;
}


int
coot::pir_alignment_t::description_to_resno_start(const std::string &descr) const {
   int rn = -1;
   std::pair<bool, int> r = description_split(descr, 2);
   if (r.first)
      rn = r.second;
   return rn;
}

// split the input string descr into fields based on ":" separator.
// return the field_number field as an int.
//
std::pair<bool, int>
coot::pir_alignment_t::description_split(const std::string &descr, int field_numpber) const {

   std::pair<bool, int> r(false,-999);

   std::vector<std::string> v = util::split_string(descr, ":");
   int v_size = v.size(); // unsigned int conversion
   if (v_size > field_numpber) {
      const std::string &field = v[field_numpber];
      try {
	 r.second = util::string_to_int(field);
	 r.first = true;
      }
      catch (std::runtime_error &e) {
	 // who cares
      }
   }
   return r;
}




void
coot::pir_alignment_t::init(const std::string &s) {

   int n = s.size();

   bool found_newline   = false;
   bool found_textdescr = false; // true on finding a new line for the second time
   bool found_greater   = false;
   std::string seq;
   // std::vector<std::pair<int, std::string> > seqs;
   std::vector<pir_t> seqs;
   std::string running; // so that we can capture the > line and the description
   int resno_start = -1;
   int resno_end   = -1;

   for (int i=0; i<n; i++) {

      if (found_newline && found_greater && found_textdescr) {
	 char t = std::toupper(s[i]);
         // std::cout << "considering t " << t << std::endl;
	 if (is_pir_aa(t)) {
	    seq += t;
	 }
	 if (t == '*') { // end of sequence
	    if (seq.size()) {
	       int rs = description_to_resno_start(running);
	       if (rs > 0)
		  resno_start = rs;
	       std::pair<bool, int> re = description_split(running, 4);
	       if (re.first)
		  resno_end = re.second;
	       // std::pair<int, std::string> p(resno_start, seq);
	       pir_t p(resno_start, resno_end, seq);
	       seqs.push_back(p);
	       seq.clear();
	       running.clear();
	       resno_end = -1;
	    }
	 }
      } else {
	 if (found_greater && found_newline) {
	    running += s[i];
	 }
      }
      if (s[i] == '>') {
	 if (seq.size()) {
	    int rs = description_to_resno_start(running);
	    if (rs > 0)
	       resno_start = rs;
	    std::pair<bool, int> re = description_split(running, 4);
	    if (re.first)
	       resno_end = re.second;
	    // std::pair<int, std::string> p(resno_start, seq);
	    pir_t p(resno_start, resno_end, seq);
	    seqs.push_back(p);
	    seq.clear();
	    resno_end = -1;
	 }
	 found_greater = true;
	 found_textdescr = false;
	 found_newline = false;
      }
      if (s[i] == '\n') {
	 if (found_newline) {
	    found_textdescr = true;
	 }
	 if (found_greater) {
	    found_newline = true;
	 }
	 
      }
   }
   if (seq.size()) {
      int rs = description_to_resno_start(running);
      if (rs > 0)
	 resno_start = rs;
      // std::pair<int, std::string> p(resno_start, seq);
      std::pair<bool, int> re = description_split(running, 4);
      if (re.first)
	 resno_end = re.second;
      pir_t p(resno_start, resno_end, seq);
      seqs.push_back(p);
   }

   std::cout << "INFO:: pir_alignment_t::init() found " << seqs.size()
             << " sequences " << std::endl;
   for (std::size_t i=0; i<seqs.size(); i++)
      std::cout << " " << i << " " << seqs[i].seq << std::endl;

   store(seqs);
}

void
coot::pir_alignment_t::store(const std::vector<coot::pir_t> &seqs) {

   if (seqs.size() > 1) {
      const std::string &ref_seq = seqs[0].seq;
      resno_start_structure = seqs[0].resno_start;
      resno_end_structure   = seqs[0].resno_end;
      for (std::size_t i=1; i<seqs.size(); i++) { // others to first (ref_seq)
	 const std::string &seq = seqs[i].seq;
	 resno_start = seqs[i].resno_start;
	 resno_end   = seqs[i].resno_end;
	 if (seq.size() == ref_seq.size()) {
	    std::vector<matched_residue_t> res_vec(seq.size());
	    for (std::size_t j=0; j<seq.size(); j++) {
	       matched_residue_t m(seq[j], ref_seq[j]);
	       res_vec[j] = m;
	    }
	    matches.push_back(res_vec);
	 } else {
	    std::cout << "size mismatch " << seq.size() << " " << ref_seq.size()
		      << std::endl;
	 }
      }
   }
}


void
coot::pir_alignment_t::read_file(const std::string &file_name) {

   std::string s;
   if (file_exists(file_name)) {
      std::ifstream f(file_name.c_str());
      std::string line;
      while (std::getline(f, line)) {
	 s += line;
	 s += '\n';
      }
   }

   init(s);
}
