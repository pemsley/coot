

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "utils/coot-utils.hh"
#include "graphics-info.h"
#include "cc-interface-scripting.hh"

#include "c-interface.h" // for is_valid_model_molecule()
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh" // for add_to_history()

#include "get-monomer.hh"




// return a new molecule number
int get_monomer_molecule_by_network_and_dict_gen(const std::string &text) {

   int imol = -1;

   std::string function_name = "get-pdbe-cif-for-comp-id";
   std::vector<coot::command_arg_t> args;
   args.push_back(coot::util::single_quote(text));
   coot::command_arg_t retval = coot::scripting_function(function_name, args);
   if (retval.type == coot::command_arg_t::STRING) {
      std::string file_name = retval.s;
      args.clear();
      args.push_back(coot::util::single_quote(text));
      args.push_back(coot::util::single_quote(file_name));
      retval = coot::scripting_function("generate-molecule-from-mmcif", args);
      if (retval.type == coot::command_arg_t::INT) {
	 imol = retval.type;
      }
   } 
   
   return imol;
} 


// Return the new molecule number, or else a negitive error code.
// 
int get_monomer(const char *three_letter_code) {

   int imol = -1;

   // first check if three_letter_code is valid, i.e. not empty
   if (std::string(three_letter_code) == "")
     return imol;
   // fast
   imol = get_monomer_from_dictionary(three_letter_code, 1); // idealized
   if (is_valid_model_molecule(imol)) { 
      return imol;
   } else { 
      imol = get_monomer_from_dictionary(three_letter_code, 0); // non-idealized
      if (is_valid_model_molecule(imol)) { 
	 return imol;
      }
   }
	 

   // OK, the slow path, using LIBCHECK.

   std::string function_name = "monomer-molecule-from-3-let-code";
   std::vector<coot::command_arg_t> args;
   args.push_back(coot::util::single_quote(three_letter_code));

   // now add in the bespoke cif library if it was given.  It is
   // ignored in the libcheck script if cif_lib_filename is "".
   //
   // However, we only want to pass the bespoke cif library if the
   // monomer to be generated is in the cif file.
   // 
   std::string cif_lib_filename = "";
   if (graphics_info_t::cif_dictionary_filename_vec->size() > 0) {
      std::string dict_name = (*graphics_info_t::cif_dictionary_filename_vec)[0];
      coot::simple_cif_reader r(dict_name);
      if (r.has_restraints_for(three_letter_code))
	 cif_lib_filename = dict_name;
   }
   args.push_back(coot::util::single_quote(cif_lib_filename));

   coot::command_arg_t retval = coot::scripting_function(function_name, args);
   if (retval.type == coot::command_arg_t::INT) {
      imol = retval.i;
   } 

   std::vector<std::string> command_strings;
   command_strings.push_back("get-monomer");
   command_strings.push_back(coot::util::single_quote(three_letter_code));
   add_to_history(command_strings);

   return imol;
}



/* Use the protein geometry dictionary to retrieve a set of
   coordinates quickly.  There are no restraints from this method
   though. */
int get_monomer_from_dictionary(const char *three_letter_code,
				int idealised_flag) {

   int istat = -1; // unfound molecule
   graphics_info_t g;
   mmdb::Manager *mol = g.Geom_p()->mol_from_dictionary(three_letter_code, idealised_flag);
   if (mol) {
      int imol = graphics_info_t::create_molecule();
      atom_selection_container_t asc = make_asc(mol);
      std::string name = three_letter_code;
      name += "_from_dict";
      graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      move_molecule_to_screen_centre_internal(imol);
      graphics_draw();
      istat = imol;
   }
   return istat;
}


