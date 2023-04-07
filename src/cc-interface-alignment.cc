
// 20180302
// This block is what is needed to include c-interface.h
// How to include c-interface.h
// --------------------------------------------
#ifdef USE_PYTHON
#include <Python.h> // add first for _XOPEN_SOURCE order issues
#endif

#include <cstddef> // define std::ptrdiff_t

#ifdef USE_GUILE
#include <libguile.h>
#endif

#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh" // 20230404-PE for add_status_bar_text()
// --------------------------------------------

#include "cc-interface-alignment.hh"

#include <fstream>

std::string get_sequence_as_fasta_for_chain(int imol, const std::string &chain_id) {

   std::string r;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].get_sequence_as_block(chain_id);
   }
   std::string n = graphics_info_t::molecules[imol].name_sans_extension(0);
   std::string  full_seq;
   full_seq = "> ";
   full_seq += n;
   full_seq += " ";
   full_seq += chain_id;
   full_seq += "\n";
   full_seq += "\n";
   full_seq += r;
   full_seq += "\n";
   return full_seq;

}

void write_sequence(int imol, const std::string &file_name) {

   if (is_valid_model_molecule(imol)) {
      std::vector<std::string> chain_ids = graphics_info_t::molecules[imol].get_chain_ids();
      std::string all_string;
      for (auto chain_id : chain_ids) {
         std::string seq = get_sequence_as_fasta_for_chain(imol, chain_id);
         all_string +=  seq;
      }
      std::ofstream f(file_name.c_str());
      if (f) {
         f << all_string;
         f.close();
      } else {
         std::cout << "Failed to open " << file_name << std::endl;
      }
   } else {
      std::string m = "Molecule " + std::to_string(imol);
      m += " is not a valid model molecule";
      add_status_bar_text(m.c_str());
      std::cout << m << std::endl;
   }
   
}


void associate_pir_alignment(int imol, std::string chain_id, std::string pir_alignment) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].associate_pir_alignment(chain_id, pir_alignment);
   }
}

void associate_pir_alignment_from_file(int imol, std::string chain_id, std::string pir_alignment_file_name) {

   if (is_valid_model_molecule(imol)) {
      if (coot::file_exists(pir_alignment_file_name)) {
	 std::string s;
	 std::ifstream f(pir_alignment_file_name.c_str());
	 std::string line;
	 while (std::getline(f, line)) {
	    s += line;
	    s += '\n';
	 }
	 graphics_info_t::molecules[imol].associate_pir_alignment(chain_id, s);
      }
   }
}

void apply_pir_alignment(int imol, std::string chain_id) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].apply_pir_alignment(chain_id);
   }
   graphics_draw();
}
