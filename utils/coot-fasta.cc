/* coot-utils/coot-fasta.cc
 * 
 * Copyright 2006, by The University of York
 * Copyright 2013 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <iostream>
#include <fstream>
#include "coot-utils.hh"
#include "coot-fasta.hh"

// http://ngfnblast.gbf.de/docs/fasta.html
// 
// For those programs that use amino acid query sequences (BLASTP and
// TBLASTN), the accepted amino acid codes are:
// 
//     A  alanine                         P  proline
//     B  aspartate or asparagine         Q  glutamine
//     C  cystine                         R  arginine
//     D  aspartate                       S  serine
//     E  glutamate                       T  threonine
//     F  phenylalanine                   U  selenocysteine
//     G  glycine                         V  valine
//     H  histidine                       W  tryptophane
//     I  isoleucine                      Y  tyrosine
//     K  lysine                          Z  glutamate or glutamine
//     L  leucine                         X  any
//     M  methionine                      *  translation stop
//     N  asparagine                      -  gap of indeterminate length

// sequence
coot::fasta::fasta(const std::string &sequence_chain_id_in, const std::string &fasta_seq_in, int type) {

   // format "> name\n <sequence>", we ignore everything that is not a
   // letter after the newline.

   // sequence is member data.  Let's fill it.

   if (type == SIMPLE_STRING) {

      name = sequence_chain_id_in;
      sequence = fasta_seq_in;

   } else {

      // FASTA_BLOCK

      std::string seq;

      int nchars = fasta_seq_in.length();
      short int found_greater = 0;
      short int found_newline = 0;
      std::string t;

      for (int i=0; i<nchars; i++) {

         //       std::cout << "checking character: " << seq_in[i] << std::endl;

         if (found_newline && found_greater) {
            t = toupper(fasta_seq_in[i]);
            if (is_fasta_aa(t)) {
               // std::cout << "adding character: " << seq_in[i] << std::endl;
               seq += t;
            }
         }
         if (fasta_seq_in[i] == '>') {
            // std::cout << "DEBUG:: " << seq_in[i] << " is > (greater than)\n";
            found_greater = 1;
         }
         if (fasta_seq_in[i] == '\n') {
            if (found_greater) {
               // std::cout << "DEBUG:: " << seq_in[i] << " is carriage return\n";
               found_newline = 1;
            }
         }
      }

      if (seq.length() > 0) {
         sequence = seq;
         name = sequence_chain_id_in;
      } else {
         sequence = "";
         name = "";
         std::cout << "WARNING:: fasta constructor: no sequence found or improper FASTA sequence format\n";
      }
   }
}

coot::fasta::fasta(const std::string &combined_string) { // decomposition happens in constructor

   // format "> name\n <sequence>", we ignore everything that is not a
   // letter after the newline.
   //
   // We also want to deal with just plain sequence text (that doesn't have a name).

   // sequence is member data.  Let's fill it.

   std::string seq;

   int nchars = combined_string.length();
   short int found_greater = 0;
   short int found_newline = 0;
   std::string t;

   // If the first non-blank character a ">" then we are dealing with
   // a fasta sequence, if not, then just a simple set of letters.
   //
   bool is_fasta = 0;
   for (int i=0; i<nchars; i++) {
      if (combined_string[i] == ' '  ||
	  combined_string[i] == '\n' ||
	  combined_string[i] == '\t') {
      } else {
	 if (combined_string[i] == '>') {
	    is_fasta = true;
	 }
	 break;
      }
   }

   // set the name (class variable) if this is a a fasta
   if (is_fasta) {
      name = "";
      short int in_name = 0;
      for (int i=0; i<nchars; i++) {
	 if (in_name)
	    name += combined_string[i];
	 if (combined_string[i] == '>') {
	    in_name = 1;
	 }
	 if (combined_string[i] == '\n') {
	    break;
	 }
      }
   }

   for (int i=0; i<nchars; i++) {

      if ((is_fasta && found_newline && found_greater) || !is_fasta) {
	 t = toupper(combined_string[i]);
	 if (is_fasta_aa(t)) {
	    seq += t;
	 }
      }
      if (combined_string[i] == '>') {
	 found_greater = 1;
      }
      if (combined_string[i] == '\n') {
	 if (found_greater) {
	    found_newline = 1;
	 }
      }
   }

   if (seq.length() > 0) {
      sequence = seq;
      name = "unknown";
   } else {
      sequence = "";
      name = "";
      std::cout << "WARNING:: no sequence found or improper fasta sequence format\n";
   }
}

bool
coot::fasta::is_fasta_aa(const std::string &a) const { 

   short int r = 0;
   
   if (a == "A" || a == "G" ) { 
      r = 1;
   } else { 
      if (a == "B" 
	  || a == "C" || a == "D" || a == "E" || a == "F" || a == "H" || a == "I"
	  || a == "K" || a == "L" || a == "M" || a == "N" || a == "P" || a == "Q" 
	  || a == "R" || a == "S" || a == "T" || a == "U" || a == "V" || a == "W" 
	  || a == "Y" || a == "Z" || a == "X" || a == "*" || a == "-") {
	 r = 1;
      }
   }
   return r;
}

void
coot::fasta_multi::read(const std::string &file_name) {

   if (file_exists(file_name)) {
      std::ifstream f(file_name.c_str());
      if (f) {
	 std::vector<std::string> lines;
	 std::string line;
	 while (std::getline(f, line)) { 
	    lines.push_back(line);
	 }
	 if (lines.size()) {
	    bool found_greater = false;
	    std::string running_seq;
	    std::string running_name;
	    fasta running_fasta;
	    for (unsigned int iline=0; iline<lines.size(); iline++) { 
	       std::string &line = lines[iline];
	       unsigned int n_chars = line.length();
	       for (unsigned int ic=0; ic<n_chars; ic++) {
		  if (line[ic] == '>') {
		     found_greater = true;
		     if (running_seq.length()) {
			fasta fas(running_name, running_seq, "direct");
			sequences.push_back(fas);
		     }
		     unsigned int l = line.length();
		     if (l>ic+1)
			running_name = line.substr(ic+1,l);
		     else
			running_name.clear();
		     running_seq.clear();
		     break;
		  }

		  if (found_greater) {
		     char t = toupper(line[ic]);
		     std::string ss;
		     ss += t;
		     if (running_fasta.is_fasta_aa(ss)) {
			running_seq += t;
		     } 
		  }
	       }
	    }
	    if (running_seq.length())
	       sequences.push_back(fasta(running_name, running_seq, "direct"));
	 }
      }
   } else {
      std::cout << "File not found: " << file_name << std::endl;
   }
}
