/* src/get-monomer.cc
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2015 by Medical Research Council
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
int get_monomer(const std::string &comp_id_in) {

   int imol = -1;

   std::string comp_id = comp_id_in;

   // first check if three_letter_code is valid, i.e. not empty
   if (comp_id.empty())
     return imol;
   // fast
   imol = get_monomer_from_dictionary(comp_id, 1); // idealized

   if (is_valid_model_molecule(imol)) { 
      return imol;
   } else {
      std::cout << "get_monomer(): trying non-idealized: " << comp_id_in << std::endl;
      imol = get_monomer_from_dictionary(comp_id, 0); // non-idealized
      std::cout << "   got imol " << imol << std::endl;
      if (is_valid_model_molecule(imol)) { 
	 return imol;
      }
   }

   if (coot::util::is_standard_residue_name(comp_id_in)) {

      // this is ugly. get_standard_residue_instance should be in coot-utils
      //
      molecule_class_info_t molci;
      mmdb::Residue *std_res = molci.get_standard_residue_instance(comp_id_in);

      if (std_res == NULL) {
	 std::cout << "WARNING:: Can't find standard residue for "
                   << comp_id_in << "\n";
      } else {
	 // happy path
	 graphics_info_t g;
	 std_res->seqNum = 1;
	 mmdb::Manager *std_res_mol = coot::util::create_mmdbmanager_from_residue(std_res);
	 atom_selection_container_t asc = make_asc(std_res_mol);
	 imol = g.create_molecule();
	 g.molecules[imol].install_model(imol, asc, g.Geom_p(), comp_id_in, true);
	 move_molecule_to_screen_centre_internal(imol);
	 graphics_draw();
      }

   } else {

      // OK, the slow path, using LIBCHECK.

      std::string function_name = "monomer-molecule-from-3-let-code";
      std::vector<coot::command_arg_t> args;
      args.push_back(coot::util::single_quote(comp_id));

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
	 if (r.has_restraints_for(comp_id))
	    cif_lib_filename = dict_name;
      }
      args.push_back(coot::util::single_quote(cif_lib_filename));
      coot::command_arg_t retval = coot::scripting_function(function_name, args);
      if (retval.type == coot::command_arg_t::INT) {
	 imol = retval.i;
      }
   }

   std::vector<std::string> command_strings;
   command_strings.push_back("get-monomer");
   command_strings.push_back(coot::util::single_quote(comp_id));
   add_to_history(command_strings);

   return imol;
}

//! get the monomer for the given 
int get_monomer_for_molecule(const std::string &comp_id, int imol) {

   graphics_info_t g;

   bool idealised_flag = true;
   mmdb::Manager *mol = g.Geom_p()->mol_from_dictionary(comp_id, imol, idealised_flag);
   if (mol) {
      imol = graphics_info_t::create_molecule();
      atom_selection_container_t asc = make_asc(mol);
      std::string name = comp_id;
      name += "_from_dict";
      graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      move_molecule_to_screen_centre_internal(imol);
      graphics_draw();
   }
   return imol;
}



/* Use the protein geometry dictionary to retrieve a set of
   coordinates quickly.  There are no restraints from this method
   though. */
int get_monomer_from_dictionary(const std::string &comp_id,
				int idealised_flag) {

   int istat = -1; // unfound molecule
   graphics_info_t g;

   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   mmdb::Manager *mol = g.Geom_p()->mol_from_dictionary(comp_id, imol_enc, idealised_flag);
   if (mol) {
      int imol = graphics_info_t::create_molecule();
      atom_selection_container_t asc = make_asc(mol);
      std::string name = comp_id;
      name += "_from_dict";
      std::cout << "debug:: get_monomer_from_dictionary() installing " << name << " into model " << imol << std::endl;
      graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      move_molecule_to_screen_centre_internal(imol);
      graphics_draw();
      istat = imol;
   } else {
      std::cout << "WARNING:: Null mol from mol_from_dictionary() with comp_id " << comp_id << " "
		<< idealised_flag << std::endl;
   }
   return istat;
}


int get_monomer_for_molecule_by_index(int dict_idx, int imol_enc) {

   graphics_info_t g;
   int imol = -1;

   int idealised_flag = true;
   mmdb::Manager *mol = g.Geom_p()->mol_from_dictionary(dict_idx, imol_enc, idealised_flag);
   if (mol) {
      imol = graphics_info_t::create_molecule();
      atom_selection_container_t asc = make_asc(mol);
      std::string name;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (! model_p) {
	 std::cout << "Null model" << std::endl;
      } else { 
	 mmdb::Chain *chain_p;
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    if (! chain_p) {
	       std::cout << "Null chain" << std::endl;
	    } else { 
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::Residue *residue_p;
	       mmdb::Atom *at;
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  name = residue_p->GetResName();
		  break;
	       }
	    }
	    if (! name.empty())
	       break;
	 }
      }
      
      name += "_from_dict";
      graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      move_molecule_to_screen_centre_internal(imol);
      graphics_draw();
   }

   return imol;

}
