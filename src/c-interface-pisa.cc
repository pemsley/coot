/* src/main.cc
 * 
 * Copyright 2009, 2010 The University of Oxford
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#include <string>
#include <iostream> // fixes undefined strchr, strchrr problems
#include <mmdb2/mmdb_manager.h>
#include <gtk/gtk.h>

#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-generic-objects.h"
#include "cc-interface.hh"
#include "c-interface-scm.hh"
#include "c-interface-python.hh"
#include "guile-fixups.h"

#include "pisa-interface.hh"

/*  ----------------------------------------------------------------------- */
/*                  PISA Interface                                      */
/*  ----------------------------------------------------------------------- */

// return the new model_number or -1;
//
// old function, made before pisa_interfaces.
// 
int pisa_interaction(int imol_1, int imol_2) {

   int imodel_new = -1;

   float dist = 4.0;
   if (is_valid_model_molecule(imol_1)) { 
      if (is_valid_model_molecule(imol_2)) {

	 mmdb::Manager *mol1 = graphics_info_t::molecules[imol_1].atom_sel.mol;
	 mmdb::Manager *mol2 = graphics_info_t::molecules[imol_2].atom_sel.mol;

	 coot::close_residues_from_different_molecules_t cr;
	 std::pair<std::vector<mmdb::Residue *>, std::vector<mmdb::Residue *> > res_pair = 
	    cr.close_residues(mol1, mol2, dist);
	 graphics_info_t g;
	 
	 if (res_pair.first.size() > 0) { 
	    std::pair<bool, mmdb::Manager *> nm =
	       coot::util::create_mmdbmanager_from_residue_vector(res_pair.first, mol1);
	    if (nm.second) {
	       int imol = graphics_info_t::create_molecule();
	       atom_selection_container_t asc = make_asc(nm.second);
	       std::string name = "interacting residues from ";
	       name += coot::util::int_to_string(imol_1);
	       graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
	       imodel_new = imol;
	    } else {
	       std::cout << "WARNING:: no molecule from create_mmdbmanager_from_residue_vector"
			 << std::endl;
	    } 
	 }
	 
	 if (res_pair.second.size() > 0) { 
	    std::pair<bool, mmdb::Manager *> nm =
	       coot::util::create_mmdbmanager_from_residue_vector(res_pair.second, mol1);
	    if (nm.second) {
	       int imol = graphics_info_t::create_molecule();
	       atom_selection_container_t asc = make_asc(nm.second);
	       std::string name = "interacting residues from ";
	       name += coot::util::int_to_string(imol_2);
	       graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
	    } else {
	       std::cout << "WARNING:: no molecule from create_mmdbmanager_from_residue_vector"
			 << std::endl;
	    } 
	 }

	 cr.clean_up();
	 graphics_draw();
      }
   }
   return imodel_new;
} 


#ifdef USE_GUILE
SCM handle_pisa_interfaces_scm(SCM interfaces_description_scm) { 

   // coot::pisa_molecule_t pisa_molecule_1;
   // coot::pisa_molecule_t pisa_molecule_2;

   // coot::pisa_interface_t pi(imol_1, imol_2, pisa_molecule_1, pisa_molecule_2);

//      std::cout << "interfaces_description_scm: "
//   	     << scm_to_locale_string(display_scm(interfaces_description_scm))
//   	     << std::endl;
//    std::cout << "pisa_molecule_record_type: "
// 	     << scm_to_locale_string(display_scm(pisa_molecule_record_type))
// 	     << std::endl;
//    std::cout << "pisa_residue_record_type: "
// 	     << scm_to_locale_string(display_scm(pisa_residue_record_type))
// 	     << std::endl;

   std::vector<coot::pisa_interface_t> pisa_treeview_info; // if this is not empty
					  		   // make a treeview at the
							   // end of this function.
   
   SCM interfaces_description_length_scm = scm_length(interfaces_description_scm);
   unsigned int interfaces_description_length = scm_to_int(interfaces_description_length_scm);
   std::cout << "INFO:: found " << interfaces_description_length << " interfaces"
	     << std::endl;

   int imol_1 = -1; // set by reading a molecule pair, used in
		    // construction of bonds
   int imol_2 = -1;
   
   for (unsigned int i=0; i<interfaces_description_length; i++) {
      SCM interface_scm = scm_list_ref(interfaces_description_scm, SCM_MAKINUM(i));
      SCM interface_length_scm = scm_length(interface_scm);
      int interface_length = scm_to_int(interface_length_scm);
      // interfaces contain a molecule pair pair.  a molecule pair is
      // simply a list the first of which is the new molecule number
      // and the second is the molecule record.
      if (interface_length == 6) {

	 SCM interface_area_scm    = scm_list_ref(interface_scm, SCM_MAKINUM(2));
	 SCM interface_solv_en_scm = scm_list_ref(interface_scm, SCM_MAKINUM(3));
	 SCM interface_pvalue_scm  = scm_list_ref(interface_scm, SCM_MAKINUM(4));
	 SCM interface_stab_en_scm = scm_list_ref(interface_scm, SCM_MAKINUM(5));

	 double interface_area    = -9999.9;
	 double interface_solv_en = -9999.9;
	 double interface_pvalue  = -9999.9;
	 double interface_stab_en = -9999.9;

 	 if (scm_is_true(scm_number_p(interface_area_scm)))
 	    interface_area = scm_to_double(interface_area_scm);
	 else
 	    std::cout << "Not a number: " << scm_to_locale_string(display_scm(interface_area_scm))
		      << std::endl;
	 
 	 if (scm_is_true(scm_number_p(interface_solv_en_scm)))
 	    interface_solv_en = scm_to_double(interface_solv_en_scm);
 	 else
 	    std::cout << "Not a number: " << scm_to_locale_string(display_scm(interface_solv_en_scm))
		      << std::endl;
	 
 	 if (scm_is_true(scm_number_p(interface_pvalue_scm)))
 	    interface_pvalue = scm_to_double(interface_pvalue_scm);
 	 else
 	    std::cout << "Not a number: " << scm_to_locale_string(display_scm(interface_pvalue_scm))
		      << std::endl;

 	 if (scm_is_true(scm_number_p(interface_stab_en_scm)))
 	    interface_solv_en = scm_to_double(interface_stab_en_scm);
 	 else
 	    std::cout << "Not a number: " << scm_to_locale_string(display_scm(interface_stab_en_scm))
		      << std::endl;

	 
	 SCM molecule_pair_pair = scm_list_ref(interface_scm, SCM_MAKINUM(0));
	 SCM molecule_pair_pair_length_scm = scm_length(molecule_pair_pair);
	 int molecule_pair_pair_length = scm_to_int(molecule_pair_pair_length_scm);

	 if (molecule_pair_pair_length == 2) {
	    // 2 molecules in an interface (good) :-)

	    SCM molecule_pair_1 =  scm_list_ref(molecule_pair_pair, SCM_MAKINUM(0));
	    SCM molecule_pair_2 =  scm_list_ref(molecule_pair_pair, SCM_MAKINUM(1));

	    SCM molecule_number_and_symop_1_scm = scm_list_ref(molecule_pair_1, SCM_MAKINUM(0));
	    SCM molecule_number_and_symop_2_scm = scm_list_ref(molecule_pair_2, SCM_MAKINUM(0));

	    SCM molecule_number_1_scm = scm_list_ref(molecule_number_and_symop_1_scm, SCM_MAKINUM(0));
	    SCM molecule_number_2_scm = scm_list_ref(molecule_number_and_symop_2_scm, SCM_MAKINUM(0));

	    SCM molecule_1_symop_scm = scm_list_ref(molecule_number_and_symop_1_scm, SCM_MAKINUM(1));
	    SCM molecule_2_symop_scm = scm_list_ref(molecule_number_and_symop_2_scm, SCM_MAKINUM(1));

	    SCM molecule_1_record_scm = scm_list_ref(molecule_pair_1, SCM_MAKINUM(1));
	    SCM molecule_2_record_scm = scm_list_ref(molecule_pair_2, SCM_MAKINUM(1));
	    
	    SCM mol_1_residue_records = pisa_molecule_record_residues(molecule_1_record_scm);
	    SCM mol_2_residue_records = pisa_molecule_record_residues(molecule_2_record_scm);

	    SCM mol_1_chain_id_scm = pisa_molecule_record_chain_id(molecule_1_record_scm);
	    SCM mol_2_chain_id_scm = pisa_molecule_record_chain_id(molecule_2_record_scm);

	    imol_1 = scm_to_int(molecule_number_1_scm);
	    imol_2 = scm_to_int(molecule_number_2_scm);

	    SCM bonds_info_scm = scm_list_ref(interface_scm, SCM_MAKINUM(1));
	    if (scm_is_true(scm_list_p(bonds_info_scm))) {
	       SCM bonds_info_length_scm = scm_length(bonds_info_scm);
	       int bonds_info_length = scm_to_int(bonds_info_length_scm);
	       for (int ibond=0; ibond<bonds_info_length; ibond++) {
		  SCM pisa_bond_scm = scm_list_ref(bonds_info_scm, SCM_MAKINUM(ibond));
		  add_pisa_interface_bond_scm(imol_1, imol_2, pisa_bond_scm, i);
	       }
	    }
	 
	    try { 

	       std::string chain_id_1 = scm_to_locale_string(mol_1_chain_id_scm);
	       std::string chain_id_2 = scm_to_locale_string(mol_2_chain_id_scm);

	       std::string symop_1 = "---"; // unset
	       std::string symop_2 = "---"; // unset

	       if (scm_is_true(scm_string_p(molecule_1_symop_scm)))
		  symop_1 = scm_to_locale_string(molecule_1_symop_scm);
	       if (scm_is_true(scm_string_p(molecule_2_symop_scm)))
		  symop_2 = scm_to_locale_string(molecule_2_symop_scm);

	       std::vector<coot::residue_spec_t> mol_1_residue_specs = 
		  residue_records_list_scm_to_residue_specs(mol_1_residue_records, chain_id_1);
	       std::vector<coot::residue_spec_t> mol_2_residue_specs = 
		  residue_records_list_scm_to_residue_specs(mol_2_residue_records, chain_id_2);

	       clipper::Coord_orth centre = 
	       make_complementary_dotted_surfaces(imol_1, imol_2,
						  mol_1_residue_specs,
						  mol_2_residue_specs);

	       
	       SCM bonds_scm = scm_list_ref(interface_scm, SCM_MAKINUM(1));

	       coot::pisa_interface_bond_info_t pibi =
		  coot::get_pisa_interface_bond_info_scm(bonds_scm);

	       // Where is the middle of the interface?  We should
	       // pass that to so that the view recentres there when
	       // the line is clicked in the treeview.  Before the end
	       // of this function, we should recentre the view to be
	       // on the centre of the top interface.

	       coot::pisa_interface_t pisa_interface_attribs(imol_1, imol_2,
							     chain_id_1, chain_id_2,
							     symop_2,
							     centre,
							     interface_area,
							     interface_solv_en,
							     interface_pvalue,
							     interface_stab_en,
							     pibi.n_h_bonds,
							     pibi.n_salt_bridges,
							     pibi.n_cov_bonds,
							     pibi.n_ss_bonds);
	       
	       pisa_treeview_info.push_back(pisa_interface_attribs);
	       
	    }
	    catch (const std::runtime_error &rte)  {
	       std::cout << "WARNING:: " << rte.what() << std::endl;
	    }
	 }
      }
   }
   if (pisa_treeview_info.size() > 0)
      coot::pisa_interfaces_gui(pisa_treeview_info);
      
   return SCM_MAKINUM(-1);
}
#endif /* USE_GUILE */


// interface_description is:
// [molecules, bonds, area, solv_en, pvalue, stab_en]
// molecules is list (2) of molecules ([imol, symop, molecule_dictionary]):
// [[imol, symop, mol_dic],[imol2,...]]
#ifdef USE_PYTHON
PyObject *handle_pisa_interfaces_py(PyObject *interfaces_description_py) { 

   // coot::pisa_molecule_t pisa_molecule_1;
   // coot::pisa_molecule_t pisa_molecule_2;

   // coot::pisa_interface_t pi(imol_1, imol_2, pisa_molecule_1, pisa_molecule_2);

//      std::cout << "interfaces_description_scm: "
//   	     << scm_to_locale_string(display_scm(interfaces_description_scm))
//   	     << std::endl;
//    std::cout << "pisa_molecule_record_type: "
// 	     << scm_to_locale_string(display_scm(pisa_molecule_record_type))
// 	     << std::endl;
//    std::cout << "pisa_residue_record_type: "
// 	     << scm_to_locale_string(display_scm(pisa_residue_record_type))
// 	     << std::endl;

   std::vector<coot::pisa_interface_t> pisa_treeview_info; // if this is not empty
					  		   // make a treeview at the
							   // end of this function.

   int interfaces_description_length = PyObject_Length(interfaces_description_py);
   std::cout << "INFO:: found " << interfaces_description_length << " interfaces"
	     << std::endl;

   int imol_1 = -1; // set by reading a molecule pair, used in
		    // construction of bonds
   int imol_2 = -1;
   
   for (int i=0; i<interfaces_description_length; i++) {
      PyObject *interface_py = PyList_GetItem(interfaces_description_py, i);
      int interface_length = PyObject_Length(interface_py);
      // interfaces contain a molecule pair pair.  a molecule pair is
      // simply a list the first of which is the new molecule number
      // and the second is the molecule record.
      if (interface_length == 6) { 
          
          PyObject *interface_area_py    = PyList_GetItem(interface_py, 2);
          PyObject *interface_solv_en_py = PyList_GetItem(interface_py, 3);
          PyObject *interface_pvalue_py  = PyList_GetItem(interface_py, 4);
          PyObject *interface_stab_en_py = PyList_GetItem(interface_py, 5);

          double interface_area    = -9999.9;
          double interface_solv_en = -9999.9;
          double interface_pvalue  = -9999.9;
          double interface_stab_en = -9999.9;

          PyObject *tmp;
          tmp = PyFloat_FromString(interface_area_py, NULL);
          if (tmp)
              interface_area = PyFloat_AsDouble(tmp);
          else
              std::cout << "Not a number: " << PyString_AsString(display_python(interface_area_py))
                        << std::endl;
          Py_XDECREF(tmp);
	 
          tmp = PyFloat_FromString(interface_solv_en_py, NULL);
          if (tmp)
              interface_solv_en = PyFloat_AsDouble(tmp);
          else
              std::cout << "Not a number: " << PyString_AsString(display_python(interface_solv_en_py))
                        << std::endl;
          Py_XDECREF(tmp);
	 
          tmp = PyFloat_FromString(interface_pvalue_py, NULL);
          if (tmp)
              interface_pvalue = PyFloat_AsDouble(tmp);
          else
              std::cout << "Not a number: " << PyString_AsString(display_python(interface_pvalue_py))
                        << std::endl;
          Py_XDECREF(tmp);

          tmp = PyFloat_FromString(interface_stab_en_py, NULL);
          if (tmp)
              interface_stab_en = PyFloat_AsDouble(tmp);
          else
              std::cout << "Not a number: " << PyString_AsString(display_python(interface_stab_en_py))
                        << std::endl;
          Py_XDECREF(tmp);

          PyObject *molecule_pair_pair = PyList_GetItem(interface_py, 0);
          int molecule_pair_pair_length = PyObject_Length(molecule_pair_pair);
          if (molecule_pair_pair_length == 2) {
              // 2 molecules in an interface (good) :-)

              PyObject *molecule_pair_1 =  PyList_GetItem(molecule_pair_pair, 0);
              PyObject *molecule_pair_2 =  PyList_GetItem(molecule_pair_pair, 1);

              PyObject *molecule_number_1_py = PyList_GetItem(molecule_pair_1, 0);
              PyObject *molecule_number_2_py = PyList_GetItem(molecule_pair_2, 0);

              PyObject *molecule_1_symop_py = PyList_GetItem(molecule_pair_1, 1);
              PyObject *molecule_2_symop_py = PyList_GetItem(molecule_pair_2, 1);

              PyObject *molecule_1_dictionary_py = PyList_GetItem(molecule_pair_1, 2);
              PyObject *molecule_2_dictionary_py = PyList_GetItem(molecule_pair_2, 2);
	    
              PyObject *mol_1_residue_records = PyDict_GetItemString(molecule_1_dictionary_py, "residues");
              PyObject *mol_2_residue_records = PyDict_GetItemString(molecule_2_dictionary_py, "residues");

              PyObject *mol_1_chain_id_py = PyDict_GetItemString(molecule_1_dictionary_py, "chain_id");
              PyObject *mol_2_chain_id_py = PyDict_GetItemString(molecule_2_dictionary_py, "chain_id");
	    
              imol_1 = PyInt_AsLong(molecule_number_1_py);
              imol_2 = PyInt_AsLong(molecule_number_2_py);
            
              PyObject *bonds_info_py = PyList_GetItem(interface_py, 1);
              if (PyList_Check(bonds_info_py)) {
                  int bonds_info_length = PyObject_Length(bonds_info_py);
                  for (int ibond=0; ibond<bonds_info_length; ibond++) {
                      PyObject *pisa_bond_py = PyList_GetItem(bonds_info_py, ibond);
                      add_pisa_interface_bond_py(imol_1, imol_2, pisa_bond_py, i);
                  }
              }

              try { 

                  std::string chain_id_1 = PyString_AsString(mol_1_chain_id_py);
                  std::string chain_id_2 = PyString_AsString(mol_2_chain_id_py);
	    
                  std::string symop_1 = "---"; // unset
                  std::string symop_2 = "---"; // unset

                  if (PyString_Check(molecule_1_symop_py))
                      symop_1 = PyString_AsString(molecule_1_symop_py);
                  if (PyString_AsString(molecule_2_symop_py))
                      symop_2 = PyString_AsString(molecule_2_symop_py);

                  std::vector<coot::residue_spec_t> mol_1_residue_specs = 
                      residue_records_list_py_to_residue_specs(mol_1_residue_records, chain_id_1);
                  std::vector<coot::residue_spec_t> mol_2_residue_specs = 
                      residue_records_list_py_to_residue_specs(mol_2_residue_records, chain_id_2);

                  clipper::Coord_orth centre = 
                      make_complementary_dotted_surfaces(imol_1, imol_2,
                                                         mol_1_residue_specs,
                                                         mol_2_residue_specs);
               
                  coot::pisa_interface_bond_info_t pibi =
                      coot::get_pisa_interface_bond_info_py(bonds_info_py);
                  // DECREF bonds_info?

                  // Where is the middle of the interface?  We should
                  // pass that to so that the view recentres there when
                  // the line is clicked in the treeview.  Before the end
                  // of this function, we should recentre the view to be
                  // on the centre of the top interface.

                  coot::pisa_interface_t pisa_interface_attribs(imol_1, imol_2,
                                                                chain_id_1, chain_id_2,
                                                                symop_2,
                                                                centre,
                                                                interface_area,
                                                                interface_solv_en,
                                                                interface_pvalue,
                                                                interface_stab_en,
                                                                pibi.n_h_bonds, 
                                                                pibi.n_salt_bridges,
                                                                pibi.n_cov_bonds,
                                                                pibi.n_ss_bonds);
	       
                  pisa_treeview_info.push_back(pisa_interface_attribs);
	       
              }
              catch (const std::runtime_error &rte)  {
                  std::cout << "WARNING:: " << rte.what() << std::endl;
              }
          }
      }
   }
   if (pisa_treeview_info.size() > 0)
       coot::pisa_interfaces_gui(pisa_treeview_info);
      
   return PyInt_FromLong(-1);
}
#endif /* USE_PYTHON */

// "[ZN]-:2" as a chain-id.  Bah. (Often just "A" though).
// 
std::string untangle_mmdb_chain_id_string(const std::string &mmdb_chain_id_in) {

   std::string s = mmdb_chain_id_in;
   std::string::size_type ibracket = mmdb_chain_id_in.find_last_of("]");
   if (ibracket != std::string::npos) {
      // it had one...
      s = mmdb_chain_id_in.substr(ibracket+1, 1);
   }
   if (s == "-")
      s = "";
   // std::cout << "=== converted \"" << mmdb_chain_id_in << "\" to \"" << s << "\"\n";
   return s;
}

   
#ifdef USE_GUILE   
coot::pisa_interface_bond_info_t
coot::get_pisa_interface_bond_info_scm(SCM bonds_info_scm) {
   coot::pisa_interface_bond_info_t pibi;

   SCM bonds_info_length_scm = scm_length(bonds_info_scm);
   int bonds_info_length = scm_to_int(bonds_info_length_scm);
   for (int ib=0; ib<bonds_info_length; ib++) {
      SCM bond_info_scm = scm_list_ref(bonds_info_scm, SCM_MAKINUM(ib));
      SCM bond_info_scm_length_scm = scm_length(bond_info_scm);
      int bond_info_scm_length = scm_to_int(bond_info_scm_length_scm);
      if (bond_info_scm_length == 3) {
	 SCM bond_type_scm = scm_list_ref(bond_info_scm, SCM_MAKINUM(0));
	 std::string bond_str = scm_to_locale_string(scm_symbol_to_string(bond_type_scm));
	 // 'h-bonds 'salt-bridges 'ss-bonds 'cov-bonds
	 if (bond_str == "h-bonds") { 
	    pibi.n_h_bonds++;
	 }
	 if (bond_str == "salt-bridges") { 
	    pibi.n_salt_bridges++;
	 }
	 if (bond_str == "ss-bonds") { 
	    pibi.n_ss_bonds++;
	 }
	 if (bond_str == "cov-bonds") { 
	    pibi.n_cov_bonds++;
	 }
      }
   }
   return pibi;
}
#endif   

#ifdef USE_PYTHON   
coot::pisa_interface_bond_info_t
coot::get_pisa_interface_bond_info_py(PyObject *bonds_info_py) {
   coot::pisa_interface_bond_info_t pibi;

   int bonds_info_length = PyObject_Length(bonds_info_py);
   for (int ib=0; ib<bonds_info_length; ib++) {
       PyObject *bond_info_py = PyList_GetItem(bonds_info_py, ib);
       int bond_info_py_length = PyObject_Length(bond_info_py);
       if (bond_info_py_length == 3) {
           PyObject *bond_type_py = PyList_GetItem(bond_info_py, 0);
           std::string bond_str = PyString_AsString(bond_type_py);
           // 'h-bonds 'salt-bridges 'ss-bonds 'cov-bonds
           if (bond_str == "h-bonds") { 
               pibi.n_h_bonds++;
           }
           if (bond_str == "salt-bridges") { 
               pibi.n_salt_bridges++;
           }
           if (bond_str == "ss-bonds") { 
               pibi.n_ss_bonds++;
           }
           if (bond_str == "cov-bonds") { 
               pibi.n_cov_bonds++;
           }
       }
   }
   return pibi;
}
#endif   



#ifdef USE_GUILE
std::vector<coot::residue_spec_t> 
residue_records_list_scm_to_residue_specs(SCM mol_1_residues, const std::string &mmdb_chain_id_in) {

   std::vector<coot::residue_spec_t> r;
   std::string chain_id = untangle_mmdb_chain_id_string(mmdb_chain_id_in);

   if (scm_is_true(scm_list_p(mol_1_residues))) { 
      SCM residues_length_scm = scm_length(mol_1_residues);
      int residues_length = scm_to_int(residues_length_scm);

      SCM record_type_descriptor_func = scm_variable_ref(scm_c_lookup("record-type-descriptor"));
      SCM record_accessor_func = scm_variable_ref(scm_c_lookup("record-accessor"));
      // 
      SCM symbol_seq_num_scm = scm_str2symbol("seq-num");
      // 
      SCM symbol_ins_code_scm = scm_str2symbol("ins-code");
      // 
      
      for (int ires=0; ires<residues_length; ires++) {
	 SCM residue_record_scm = scm_list_ref(mol_1_residues, SCM_MAKINUM(ires));
	 SCM rec_type_scm = scm_call_1(record_type_descriptor_func, residue_record_scm);
	 SCM get_seq_num_func = scm_call_2(record_accessor_func, rec_type_scm, symbol_seq_num_scm);
	 SCM get_ins_code_func = scm_call_2(record_accessor_func, rec_type_scm, symbol_ins_code_scm);
	 SCM seq_num_scm = scm_call_1(get_seq_num_func,   residue_record_scm);
	 SCM ins_code_scm = scm_call_1(get_ins_code_func, residue_record_scm);

	 try { 
	    std::string seq_num_string = scm_to_locale_string(seq_num_scm);
	    int seq_num = coot::util::string_to_int(seq_num_string);
	    std::string ins_code = scm_to_locale_string(ins_code_scm);
	    
	    coot::residue_spec_t rs(chain_id, seq_num, ins_code);
	    r.push_back(rs);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING bad seq-num from pisa interfaces xml "
		      << scm_to_locale_string(display_scm(seq_num_scm)) << std::endl;
	 }
      }
   }
   return r;
}
#endif

#ifdef USE_PYTHON
std::vector<coot::residue_spec_t> 
residue_records_list_py_to_residue_specs(PyObject *mol_1_residues, const std::string &mmdb_chain_id_in) {

   std::vector<coot::residue_spec_t> r;
   std::string chain_id = untangle_mmdb_chain_id_string(mmdb_chain_id_in);

   if (PyList_Check(mol_1_residues)) { 
      int residues_length = PyObject_Length(mol_1_residues);
      
      for (int ires=0; ires<residues_length; ires++) {
          PyObject *seq_num_py  = PyDict_GetItemString(PyList_GetItem(mol_1_residues, ires), "seq_num");
          PyObject *ins_code_py = PyDict_GetItemString(PyList_GetItem(mol_1_residues, ires), "ins_code");

	 try { 
	    std::string seq_num_string = PyString_AsString(seq_num_py);
	    int seq_num = coot::util::string_to_int(seq_num_string);
	    std::string ins_code = PyString_AsString(ins_code_py);
	    
	    coot::residue_spec_t rs(chain_id, seq_num, ins_code);
	    r.push_back(rs);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING bad seq-num from pisa interfaces xml "
		      << PyString_AsString(display_python(seq_num_py)) << std::endl;
	 }
      }
   }
   return r;
}
#endif


// return the data item associated with symbol in the record.
// 
#ifdef USE_GUILE
SCM symbol_value_from_record(SCM record_1, const std::string &symbol) {

   SCM record_type_descriptor_func = scm_variable_ref(scm_c_lookup("record-type-descriptor"));
   SCM rec_type_scm = scm_call_1(record_type_descriptor_func, record_1);

   SCM record_accessor_func = scm_variable_ref(scm_c_lookup("record-accessor"));
   SCM symbol_residues_scm = scm_str2symbol(symbol.c_str());
   SCM get_residues_func = scm_call_2(record_accessor_func, rec_type_scm, symbol_residues_scm);

   SCM record_residues_scm = scm_call_1(get_residues_func, record_1);

   return record_residues_scm;
}
#endif /* USE_GUILE */

#ifdef USE_GUILE
SCM pisa_molecule_record_residues(SCM molecule_record_1) {
   std::string symbol = "residues";
   return symbol_value_from_record(molecule_record_1, symbol);
}
#endif /* USE_GUILE */
//#ifdef USE_PYTHON
//PyObject *pisa_molecule_record_residues_py(PyObject *molecule_record_1) {
//   std::string symbol = "residues";
//   return symbol_value_from_record(molecule_record_1, symbol);
//}
//#endif /*USE_PYTHON */

#ifdef USE_GUILE
SCM pisa_molecule_record_chain_id(SCM molecule_record_1) {
   std::string symbol = "chain-id";
   return symbol_value_from_record(molecule_record_1, symbol);
}
#endif /* USE_GUILE */
//#ifdef USE_PYTHON
//PyObject *pisa_molecule_record_chain_id_py(PyObject *molecule_record_1) {
//   std::string symbol = "chain-id";
//   return symbol_value_from_record(molecule_record_1, symbol);
//}
//#endif /*USE_PYTHON */

#ifdef USE_GUILE
void add_pisa_interface_bond_scm(int imol_1, int imol_2, SCM pisa_bond_scm,
				 int interface_number) {

   // lookup generic object numbers, creating generic objects if
   // necessary. Each interface has its own h-bonds, salt-bridges,
   // ss-bonds and cov-bonds generic objects.

   std::string h_bonds_generic_objects_name = "H-bonds-interface-";
   h_bonds_generic_objects_name += coot::util::int_to_string(interface_number);
   int h_bonds_generic_objects_number =
      generic_object_index(h_bonds_generic_objects_name.c_str());
   if (h_bonds_generic_objects_number == -1)
      h_bonds_generic_objects_number =
	 new_generic_object_number(h_bonds_generic_objects_name.c_str());
   
   std::string salt_bridges_generic_objects_name = "salt-bridges-interface-";
   salt_bridges_generic_objects_name += coot::util::int_to_string(interface_number);
   int salt_bridges_generic_objects_number =
      generic_object_index(salt_bridges_generic_objects_name.c_str());
   if (salt_bridges_generic_objects_number == -1)
      salt_bridges_generic_objects_number =
	 new_generic_object_number(salt_bridges_generic_objects_name.c_str());

   std::string ss_bonds_generic_objects_name = "SS-bonds-interface-";
   ss_bonds_generic_objects_name += coot::util::int_to_string(interface_number);
   int ss_bonds_generic_objects_number =
      generic_object_index(ss_bonds_generic_objects_name.c_str());
   if (ss_bonds_generic_objects_number == -1)
      ss_bonds_generic_objects_number =
	 new_generic_object_number(ss_bonds_generic_objects_name.c_str());

   std::string cov_bonds_generic_objects_name = "Covalent-interface-";
   cov_bonds_generic_objects_name += coot::util::int_to_string(interface_number);
   int cov_bonds_generic_objects_number =
      generic_object_index(cov_bonds_generic_objects_name.c_str());
   if (cov_bonds_generic_objects_number == -1)
      cov_bonds_generic_objects_number =
	 new_generic_object_number(cov_bonds_generic_objects_name.c_str());

   set_display_generic_object(     h_bonds_generic_objects_number, 1);
   set_display_generic_object(salt_bridges_generic_objects_number, 1);
   set_display_generic_object(   cov_bonds_generic_objects_number, 1);
   set_display_generic_object(    ss_bonds_generic_objects_number, 1);

   // a pisa_bond_scm should be a list of:
   // bond-type: 'h-bonds, 'salt-bridges 'cov-bonds 'ss-bonds
   // an atom spec of an atom in imol_1
   // an atom spec of an atom in imol_2
   //
   if (scm_is_true(scm_list_p(pisa_bond_scm))) {
      SCM pisa_bond_length_scm = scm_length(pisa_bond_scm);
      int pisa_bond_length = scm_to_int(pisa_bond_length_scm);
      if (pisa_bond_length == 3) {
	 SCM bond_type_scm = scm_list_ref(pisa_bond_scm, SCM_MAKINUM(0));
	 SCM atom_spec_1 = scm_list_ref(pisa_bond_scm, SCM_MAKINUM(1));
	 SCM atom_spec_2 = scm_list_ref(pisa_bond_scm, SCM_MAKINUM(2));
	 int generic_object_number = -1;
	 string bond_type = "";
	 std::string colour = "grey";
	 if (scm_is_true(scm_eq_p(bond_type_scm, scm_str2symbol("h-bonds")))) { 
	    bond_type = "h-bond";
	    generic_object_number = h_bonds_generic_objects_number;
	    colour = "orange";
	 }
	 if (scm_is_true(scm_eq_p(bond_type_scm, scm_str2symbol("salt-bridges")))) {
	    bond_type = "salt-bridge";
	    generic_object_number = salt_bridges_generic_objects_number;
	    colour = "green";
	 }
	 if (scm_is_true(scm_eq_p(bond_type_scm, scm_str2symbol("cov-bonds")))) {
	    bond_type = "cov-bond";
	    generic_object_number = cov_bonds_generic_objects_number;
	    colour = "red";
	 }
	 if (scm_is_true(scm_eq_p(bond_type_scm, scm_str2symbol("ss-bonds")))) {
	    bond_type = "ss-bond";
	    generic_object_number = ss_bonds_generic_objects_number;
	    colour = "yellow";
	 }

	 if (bond_type != "") {
	    add_generic_object_bond(imol_1, imol_2,
				    atom_spec_from_scm_expression(atom_spec_1),
				    atom_spec_from_scm_expression(atom_spec_2),
				    generic_object_number, colour);
	 }
      }
   }
}
#endif

#ifdef USE_PYTHON
void add_pisa_interface_bond_py(int imol_1, int imol_2, PyObject *pisa_bond_py,
                                int interface_number) {

   // lookup generic object numbers, creating generic objects if
   // necessary. Each interface has its own h-bonds, salt-bridges,
   // ss-bonds and cov-bonds generic objects.

   std::string h_bonds_generic_objects_name = "H-bonds-interface-";
   h_bonds_generic_objects_name += coot::util::int_to_string(interface_number);
   int h_bonds_generic_objects_number =
      generic_object_index(h_bonds_generic_objects_name.c_str());
   if (h_bonds_generic_objects_number == -1)
      h_bonds_generic_objects_number =
	 new_generic_object_number(h_bonds_generic_objects_name.c_str());
   
   std::string salt_bridges_generic_objects_name = "salt-bridges-interface-";
   salt_bridges_generic_objects_name += coot::util::int_to_string(interface_number);
   int salt_bridges_generic_objects_number =
      generic_object_index(salt_bridges_generic_objects_name.c_str());
   if (salt_bridges_generic_objects_number == -1)
      salt_bridges_generic_objects_number =
	 new_generic_object_number(salt_bridges_generic_objects_name.c_str());

   std::string ss_bonds_generic_objects_name = "SS-bonds-interface-";
   ss_bonds_generic_objects_name += coot::util::int_to_string(interface_number);
   int ss_bonds_generic_objects_number =
      generic_object_index(ss_bonds_generic_objects_name.c_str());
   if (ss_bonds_generic_objects_number == -1)
      ss_bonds_generic_objects_number =
	 new_generic_object_number(ss_bonds_generic_objects_name.c_str());

   std::string cov_bonds_generic_objects_name = "Covalent-interface-";
   cov_bonds_generic_objects_name += coot::util::int_to_string(interface_number);
   int cov_bonds_generic_objects_number =
      generic_object_index(cov_bonds_generic_objects_name.c_str());
   if (cov_bonds_generic_objects_number == -1)
      cov_bonds_generic_objects_number =
	 new_generic_object_number(cov_bonds_generic_objects_name.c_str());

   set_display_generic_object(     h_bonds_generic_objects_number, 1);
   set_display_generic_object(salt_bridges_generic_objects_number, 1);
   set_display_generic_object(   cov_bonds_generic_objects_number, 1);
   set_display_generic_object(    ss_bonds_generic_objects_number, 1);

   // a pisa_bond_scm should be a list of:
   // bond-type: 'h-bonds, 'salt-bridges 'cov-bonds 'ss-bonds
   // an atom spec of an atom in imol_1
   // an atom spec of an atom in imol_2
   //
   if (PyList_Check(pisa_bond_py)) {
      int pisa_bond_length = PyObject_Length(pisa_bond_py);
      if (pisa_bond_length == 3) {
	 PyObject *bond_type_py = PyList_GetItem(pisa_bond_py, 0);
	 PyObject *atom_spec_1 = PyList_GetItem(pisa_bond_py, 1);
	 PyObject *atom_spec_2 = PyList_GetItem(pisa_bond_py, 2);
	 int generic_object_number = -1;
	 string bond_type = "";
	 std::string colour = "grey";
         std::string tmp;
	 if (strcmp(PyString_AsString(bond_type_py), "h-bonds") == 0) { 
	    bond_type = "h-bond";
	    generic_object_number = h_bonds_generic_objects_number;
	    colour = "orange";
	 }
	 if (strcmp(PyString_AsString(bond_type_py), "salt-bridges") == 0) {
	    bond_type = "salt-bridge";
	    generic_object_number = salt_bridges_generic_objects_number;
	    colour = "green";
	 }
	 if (strcmp(PyString_AsString(bond_type_py), "cov-bonds") == 0) {
	    bond_type = "cov-bond";
	    generic_object_number = cov_bonds_generic_objects_number;
	    colour = "red";
	 }
	 if (strcmp(PyString_AsString(bond_type_py), "ss-bonds") == 0) {
	    bond_type = "ss-bond";
	    generic_object_number = ss_bonds_generic_objects_number;
	    colour = "yellow";
	 }

	 if (bond_type != "") {
	    add_generic_object_bond(imol_1, imol_2,
				    atom_spec_from_python_expression(atom_spec_1),
				    atom_spec_from_python_expression(atom_spec_2),
				    generic_object_number, colour);
	 }
      }
   }
}
#endif

void
add_generic_object_bond(int imol1, int imol2,
			const coot::atom_spec_t &atom_spec_1,
			const coot::atom_spec_t &atom_spec_2,
			int generic_object_number,
			const std::string &colour) {

   if (is_valid_model_molecule(imol1)) {
      if (is_valid_model_molecule(imol2)) {
	 mmdb::Atom *at1 = graphics_info_t::molecules[imol1].get_atom(atom_spec_1);
	 mmdb::Atom *at2 = graphics_info_t::molecules[imol2].get_atom(atom_spec_2);
	 if (! at1)
	    std::cout << "WARNING:: failed to get atom from spec " << atom_spec_1
		      << " in molecule " << imol1 << "\n";
	 if (! at2)
	    std::cout << "WARNING:: failed to get atom from spec " << atom_spec_2
		      << " in molecule " << imol2 << "\n";
	 if (at1 && at2) {
	    clipper::Coord_orth pt1(at1->x, at1->y, at1->z);
	    clipper::Coord_orth pt2(at2->x, at2->y, at2->z);
	    to_generic_object_add_dashed_line(generic_object_number,
					      colour.c_str(),
					      5, 6,
					      at1->x, at1->y, at1->z,
					      at2->x, at2->y, at2->z);
	 }
      }
   }
}



void pisa_clear_interfaces() {
   

}


// If you want to do something with transparent surfaces in the
// future, then replace this function.
//
// This is called for each interface to be represented.
//
// Return the centre of the interface (so that we can centre on it
// when it's row is pressed in the treeview).
// 
// 
clipper::Coord_orth 
make_complementary_dotted_surfaces(int imol_1, int imol_2, 
				   std::vector<coot::residue_spec_t> &r1, 
				   std::vector<coot::residue_spec_t> &r2) {

   if (0) { 
      std::cout << "----------------------------- ";
      std::cout << "make_complementary_dotted_surfaces  " << imol_1 << " " << imol_2 << " ";
      std::cout << "----------------------------- ";
      std::cout << std::endl;
   }
      
   // symmetry is handled outside of here, which means that imol_2 can
   // be a symmetry copy of imol_1.
   
   // make synthetic molecules, dots where each residue contains one atom (dot).
   // then call std::pair<std::vector<mmdb::Residue *>, std::vector<mmdb::Residue *> >
   // coot::close_residues_from_different_molecules_t::close_residues(mmdb::Manager *mol1,
   //                                                                 mmdb::Manager *mol2,
   //                                                                 float dist)
   // 
   // consider making the dots generation a member function of dots_representation_info_t
   // and making it available here (if it isn't already).

   float close_dist = 4.3;
   clipper::Coord_orth centre_1(0,0,0);
   clipper::Coord_orth centre_2(0,0,0);

   if (is_valid_model_molecule(imol_1)) { 
      if (is_valid_model_molecule(imol_2)) {
	 mmdb::Manager *mol_1 = graphics_info_t::molecules[imol_1].atom_sel.mol;
	 mmdb::Manager *mol_2 = graphics_info_t::molecules[imol_2].atom_sel.mol;

	 mmdb::Manager *frag_mol_1 = coot::util::create_mmdbmanager_from_residue_specs(r1, mol_1);
	 mmdb::Manager *frag_mol_2 = coot::util::create_mmdbmanager_from_residue_specs(r2, mol_2);

	 std::pair<bool, clipper::Coord_orth> c_1 = coot::centre_of_molecule(frag_mol_1);
	 std::pair<bool, clipper::Coord_orth> c_2 = coot::centre_of_molecule(frag_mol_2);

	 if (c_1.first)
	    centre_1 = c_1.second;
	 if (c_2.first)
	    centre_2 = c_2.second;

	 coot::dots_representation_info_t d1(frag_mol_1, frag_mol_2);
	 coot::dots_representation_info_t d2(frag_mol_2, frag_mol_1);

	 graphics_info_t::molecules[imol_1].add_dots(d1);
	 graphics_info_t::molecules[imol_2].add_dots(d2);

	 graphics_info_t::molecules[imol_1].set_dots_colour(0.6, 0.6, 0.3);
	 graphics_info_t::molecules[imol_2].set_dots_colour(0.6, 0.3, 0.6);

	 // delete those temporary molecules
	 delete frag_mol_1;
	 delete frag_mol_2;
      }
   }
   graphics_draw();

   clipper::Coord_orth centre_av( (centre_1.x() + centre_2.x()) * 0.5,
				  (centre_1.y() + centre_2.y()) * 0.5,
				  (centre_1.z() + centre_2.z()) * 0.5);
   
   return centre_av;
} 

// undisplay all other coord molecules except imol_1 and imol_2
// 
void
pisa_interfaces_display_only(int imol_1, int imol_2, clipper::Coord_orth centre_pt) {

   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) {
	 if ((imol != imol_1) && (imol != imol_2)) {
	    set_mol_displayed(imol, 0);
	    set_mol_active(imol, 0);
	 }
      } 
   }
   set_mol_displayed(imol_1, 1);
   set_mol_displayed(imol_2, 1);
   set_mol_active(imol_1, 1);
   set_mol_active(imol_2, 1);
   coot::Cartesian pt(centre_pt.x(), centre_pt.y(), centre_pt.z());
   graphics_info_t g;
   g.setRotationCentre(pt);
   g.update_things_on_move_and_redraw();
}
