
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
// --------------------------------------------

#include "cc-interface-alignment.hh"

#include <fstream>

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
