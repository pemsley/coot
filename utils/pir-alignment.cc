

#include <fstream>
#include <iostream>
#include <stdexcept>
#include "coot-utils.hh"
#include "pir-alignment.hh"


std::ostream &
coot::operator<< (std::ostream &s, const coot::pir_alignment_t::matched_residue_t &m) {
   s << m.target << " " << m.aligned << " " << m.res_no;
   return s;
}

coot::pir_alignment_t::pir_alignment_t() {

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
coot::pir_alignment_t::description_to_resno_start(const std::string &descr) const 
{
   int rn = 1;

   std::vector<std::string> v = util::split_string(descr, ":");
   if (v.size() > 2) {
      const std::string &field = v[2];
      try {
	 rn = util::string_to_int(field);
      }
      catch (std::runtime_error &e) {
	 // who cares
      }
   }

   return rn;
}



void
coot::pir_alignment_t::init(const std::string &s) {

   int n = s.size();

   bool found_newline   = false;
   bool found_textdescr = false; // true on finding a new line for the second time
   bool found_greater   = false;
   std::string seq;
   std::vector<std::pair<int, std::string> > seqs;
   std::string running; // so that we can capture the > line and the description
   int resno_start = 1;

   for (int i=0; i<n; i++) {

      if (found_newline && found_greater && found_textdescr) {
	 char t = std::toupper(s[i]);
	 if (is_pir_aa(t)) {
	    seq += t;
	 }
	 if (t == '*') { // end of sequence
	    if (seq.size()) {
	       int rs = description_to_resno_start(running);
	       if (rs > 0)
		  resno_start = rs;
	       std::pair<int, std::string> p(resno_start, seq);
	       seqs.push_back(p);
	       seq.clear();
	       running.clear();
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
	    std::pair<int, std::string> p(resno_start, seq);
	    seqs.push_back(p);
	    seq.clear();
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
      std::pair<int, std::string> p(resno_start, seq);
      seqs.push_back(p);
   }

   // std::cout << "found " << seqs.size() << " sequences " << std::endl;
   for (std::size_t i=0; i<seqs.size(); i++)
      std::cout << " " << i << " " << seqs[i].second << std::endl;

   store(seqs);
}

void
coot::pir_alignment_t::store(const std::vector<std::pair<int, std::string> > &seqs) {

   if (seqs.size() > 1) {
      const std::string &ref_seq = seqs[0].second;
      std::map<int, matched_residue_t> res_map;
      for (std::size_t i=1; i<seqs.size(); i++) { // others to first (ref_seq)
	 const std::string &seq = seqs[i].second;
	 if (seq.size() == ref_seq.size()) {
	    int res_count = 0;
	    for (std::size_t j=0; j<seq.size(); j++) {
	       int resno_start = seqs[i].first;
	       int res_no = res_count + resno_start;
	       if (seq[j] != '-')
		  res_count++;
	       matched_residue_t m(seq[j], ref_seq[j], res_no);
	       res_map[res_no] = m;
	       // std::cout << " " << m << std::endl;
	    }
	 }
      }
      matches.push_back(res_map);
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
