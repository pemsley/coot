/* src/c-interface-mmdb.cc
 * 
 * Copyright 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2007 Bernhard Lohkamp
 * Copyright 2007 University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <iostream>
#include <string>
#include <string.h>

#include "compat/coot-sysdep.h"

#include "c-interface-mmdb.hh"

#include "c-interface-python.hh"

#include "guile-fixups.h"

#include "graphics-info.h"

#ifdef USE_GUILE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"

mmdb::Manager *
mmdb_manager_from_scheme_expression(SCM molecule_expression) {

   mmdb::Manager *mol = 0;

   SCM nmodel = scm_length(molecule_expression);
   int inmodel = scm_to_int(nmodel);

   if (inmodel > 0) {
      mol = new mmdb::Manager;
      for(int imodel=0; imodel<inmodel; imodel++) {
	 mmdb::Model *model_p = new mmdb::Model;
	 SCM imodel_scm = scm_from_int(imodel);
	 SCM model_expression = scm_list_ref(molecule_expression, imodel_scm);
	 SCM model_expression_length = scm_length(model_expression);
	 int len_model_expression = scm_to_int(model_expression_length);
	 if (len_model_expression == 0) {
	    std::cout << "model length zero!\n";
	 } else {
	    // SCM chain_list = model_expression; // interesting
	    int nchains = len_model_expression;

	    for (int ichain=0; ichain<nchains; ichain++) {

	       SCM chain_expression = scm_list_ref(model_expression,
						   scm_from_int(ichain));
	       SCM chain_is_list_scm = scm_list_p(chain_expression);
	       if (scm_is_true(chain_is_list_scm)) {
		  // printf("chain_expression is a list\n");
	       } else {
		  printf("chain_expression is not a list\n");
	       }
	       // 	    display_scm(chain_expression);
	       SCM chain_expression_length = scm_length(chain_expression);
	       int len_chain_expr = scm_to_int(chain_expression_length);
	       if (len_chain_expr != 2) {
		  std::cout << "bad chain expression, length "
			    << len_chain_expr << std::endl;
	       } else {
		  // normal case
		  // std::cout << "good chain expression " << std::endl;
		  SCM chain_id_scm = scm_list_ref(chain_expression, scm_from_int(0));
		  SCM residues_list = scm_list_ref(chain_expression, scm_from_int(1));
		  if (scm_is_true(scm_list_p(residues_list))) {
		     // printf("residues_list is a list\n");
		  } else {
		     printf("residue_list is not a list\n");
		  }
		  SCM n_residues_scm = scm_length(residues_list);
		  int n_residues = scm_to_int(n_residues_scm);
		  // printf("there were %d residues\n", n_residues);
		  if (n_residues > 0) {
		     mmdb::Chain *chain_p = new mmdb::Chain;
		     std::string chain_id = scm_to_locale_string(chain_id_scm);
		     chain_p->SetChainID(chain_id.c_str());
		     for (int ires=0; ires<n_residues; ires++) {
			SCM ires_scm = scm_from_int(ires);
			SCM scm_residue = scm_list_ref(residues_list, ires_scm);
			SCM scm_len_residue_expr = scm_length(scm_residue);
			int len_residue_expr = scm_to_int(scm_len_residue_expr);
			if (len_residue_expr != 4) {
			   std::cout << "bad residue expression, length "
				     << len_residue_expr << std::endl;
			} else {
			   // normal case
			   // std::cout << "good residue expression" << std::endl;
			   SCM scm_residue_number = scm_list_ref(scm_residue, scm_from_int(0));
			   SCM scm_residue_inscode = scm_list_ref(scm_residue, scm_from_int(1));
			   SCM scm_residue_name = scm_list_ref(scm_residue, scm_from_int(2));
			   SCM atoms_list = scm_list_ref(scm_residue, scm_from_int(3));

			   if (scm_is_true(scm_list_p(atoms_list))) {
			      // printf("atoms_list is a list\n");
			   } else {
			      printf("atoms_list is not a list\n");
			   }
			   SCM n_atoms_scm = scm_length(atoms_list);
			   int n_atoms = scm_to_int(n_atoms_scm);

			   if (n_atoms > 0) {
			      mmdb::Residue *residue_p = new mmdb::Residue;
			      std::string resname = scm_to_locale_string(scm_residue_name);
			      int resno = scm_to_int(scm_residue_number);
			      // std::cout << "DEBUG:: Found resno:   " << resno << std::endl;
			      // std::cout << "DEBUG:: Found resname: " << resname << std::endl;
			      std::string inscode = scm_to_locale_string(scm_residue_inscode);
			      residue_p->SetResName(resname.c_str());
			      residue_p->seqNum = resno;
			      memcpy(residue_p->insCode, inscode.c_str(), sizeof(mmdb::InsCode));
			      for (int iat=0; iat<n_atoms; iat++) {
				 SCM iat_scm = scm_from_int(iat);
				 SCM atom_expression = scm_list_ref(atoms_list, iat_scm);
				 SCM len_atom_expr_scm = scm_length(atom_expression);
				 int len_atom_expr = scm_to_int(len_atom_expr_scm);
				 if (len_atom_expr != 3) {
				    std::cout << "bad atom expression, length "
					      << len_residue_expr << std::endl;
				    SCM dest = SCM_BOOL_F;
				    SCM mess = scm_from_locale_string("object: ~S\n");
				    SCM bad_scm = scm_simple_format(dest, mess, scm_list_1(atom_expression));
				    std::string bad_str = scm_to_locale_string(bad_scm);
				    std::cout << bad_str << std::endl;
				 } else {
				    // normal case
				    // std::cout << "good atom expression " << std::endl;
				    SCM name_alt_conf_pair = scm_list_ref(atom_expression, scm_from_int(0));
				    SCM occ_b_ele = scm_list_ref(atom_expression, scm_from_int(1));
				    SCM pos_expr =  scm_list_ref(atom_expression, scm_from_int(2));
				    SCM len_name_alt_conf_scm = scm_length(name_alt_conf_pair);
				    int len_name_alt_conf = scm_to_int(len_name_alt_conf_scm);
				    SCM len_occ_b_ele_scm = scm_length(occ_b_ele);
				    int len_occ_b_ele = scm_to_int(len_occ_b_ele_scm);
				    SCM len_pos_expr_scm = scm_length(pos_expr);
				    int len_pos_expr = scm_to_int(len_pos_expr_scm);
				    if (len_name_alt_conf == 2) {
				       // the occ_b_ele list can contain an optional segid
				       if ( (len_occ_b_ele == 3) || (len_occ_b_ele == 4)) {
					  if (len_pos_expr == 3) {
					     SCM atom_name_scm = SCM_CAR(name_alt_conf_pair);
					     std::string atom_name = scm_to_locale_string(atom_name_scm);
					     SCM alt_conf_scm = SCM_CAR(SCM_CDR(name_alt_conf_pair));
					     std::string alt_conf = scm_to_locale_string(alt_conf_scm);
					     SCM occ_scm = scm_list_ref(occ_b_ele, scm_from_int(0));
					     SCM b_scm   = scm_list_ref(occ_b_ele, scm_from_int(1));
					     SCM ele_scm = scm_list_ref(occ_b_ele, scm_from_int(2));
					     std::string segid;
					     bool have_segid = 0;
					     if (len_occ_b_ele == 4) {
						SCM segid_scm = scm_list_ref(occ_b_ele, scm_from_int(3));
						if (scm_is_string(segid_scm)) { 
						   have_segid = 1;
						   segid = scm_to_locale_string(segid_scm);
						} 
					     } 
					     float b = scm_to_double(b_scm);
					     float occ = scm_to_double(occ_scm);
					     std::string ele = scm_to_locale_string(ele_scm);
					     float x = scm_to_double(scm_list_ref(pos_expr, scm_from_int(0)));
					     float y = scm_to_double(scm_list_ref(pos_expr, scm_from_int(1)));
					     float z = scm_to_double(scm_list_ref(pos_expr, scm_from_int(2)));
					     mmdb::Atom *atom = new mmdb::Atom;
					     atom->SetCoordinates(x, y, z, occ, b);
					     if ( ! ((atom_name == "") && (ele == ""))) {
						atom->SetAtomName(atom_name.c_str());
						atom->SetElementName(ele.c_str());
					     } else { 
						atom->SetAtomName(atom_name.c_str());
						atom->SetElementName(ele.c_str());
					     }
					     strncpy(atom->altLoc, alt_conf.c_str(), 2);
					     if (have_segid)
						strncpy(atom->segID, segid.c_str(), 5);
					     residue_p->AddAtom(atom);
					     
					     // std::cout << "DEBUG:: adding mmdb atom " << atom << std::endl;
					  } else {
					     std::cout << "bad atom (position expression) "
						       << std::endl;
					     SCM bad_scm = display_scm(pos_expr);
					     std::string bad_str = scm_to_locale_string(bad_scm);
					     std::cout << bad_str << std::endl;
					  }
				       } else {
					  std::cout << "bad atom (occ b element expression) "
						    << std::endl;
					  SCM bad_scm = display_scm(occ_b_ele);
					  std::string bad_str = scm_to_locale_string(bad_scm);
					  std::cout << bad_str << std::endl;
				       }
				    } else {
				       std::cout << "bad atom (name alt-conf expression) "
						 << std::endl;
				       SCM bad_scm = display_scm(name_alt_conf_pair);
				       std::string bad_str = scm_to_locale_string(bad_scm);
				       std::cout << bad_str << std::endl;
				    }
				 }
			      }
			      chain_p->AddResidue(residue_p);
			   }
			}
		     }
		     model_p->AddChain(chain_p);
		  }
	       }
	       mol->AddModel(model_p);
	    }
	 }
      }
   }
   if (mol) { 
      mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      mol->FinishStructEdit();
   }
   return mol;
} 

#endif  // USE_GUILE
#ifdef USE_PYTHON

mmdb::Manager *
mmdb_manager_from_python_expression(PyObject *molecule_expression) {

   mmdb::Manager *mol = 0;
   std::deque<mmdb::Model *> model_list = mmdb_models_from_python_expression(molecule_expression);
   
   if (!model_list.empty()) {
      mol = new mmdb::Manager;
      while (!model_list.empty()) {
         mol->AddModel(model_list.front());
         model_list.pop_front();
      }
   }
   
   if (mol) { 
      mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      mol->FinishStructEdit();
   }
   return mol;
} 


std::deque<mmdb::Model *>
mmdb_models_from_python_expression(PyObject *molecule_expression) {

   std::deque<mmdb::Model *> model_list;

   int inmodel = PyObject_Length(molecule_expression);

   if (inmodel > 0) {
      for(int imodel=0; imodel<inmodel; imodel++) {
         PyObject *model_expression = PyList_GetItem(molecule_expression, imodel);
         int len_model_expression = PyObject_Length(model_expression);
         if (len_model_expression > 0) {
	    mmdb::Model *model_p = new mmdb::Model;
	   //PyObject *chain_list = model_expression; // interesting
            int nchains = len_model_expression;

            for (int ichain=0; ichain<nchains; ichain++) {
               
               PyObject *chain_expression = PyList_GetItem(model_expression,
                                                   ichain);
               int chain_is_list_python = PyList_Check(chain_expression);
               if (chain_is_list_python > 0) {
                  // printf("chain_expression is a list\n");
               } else {
                  printf("chain_expression is not a list\n");
               }
               //           display_scm(chain_expression);
               int len_chain_expr = PyList_Size(chain_expression);
               if (len_chain_expr != 2) {
                  std::cout << "bad chain expression, length "
                            << len_chain_expr << std::endl;
               } else {
                  // normal case
                  // std::cout << "good chain expression " << std::endl;
                  PyObject *chain_id_python = PyList_GetItem(chain_expression, 0);
                  PyObject *residues_list = PyList_GetItem(chain_expression, 1);
                  if (PyList_Size(residues_list) > 0) {
                     // printf("residues_list is a list\n");
                  } else {
                     printf("residue_list is not a list\n");
                  }
                  int n_residues = PyObject_Length(residues_list);
		  // printf("there were %d residues\n", n_residues);
                  if (n_residues > 0) {
                     mmdb::Chain *chain_p = new mmdb::Chain;
                     std::string chain_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(chain_id_python));
                     chain_p->SetChainID(chain_id.c_str());
                     for (int ires=0; ires<n_residues; ires++) {
                        PyObject *python_residue = PyList_GetItem(residues_list, ires);
                        int len_residue_expr = PyObject_Length(python_residue);
                        if (len_residue_expr != 4) {
                           std::cout << "bad residue expression, length "
                                     << len_residue_expr << std::endl;
                        } else {
                           // normal case
                           // std::cout << "good residue expression" << std::endl;
                           PyObject *python_residue_number = PyList_GetItem(python_residue, 0);
                           PyObject *python_residue_inscode = PyList_GetItem(python_residue, 1);
                           PyObject *python_residue_name = PyList_GetItem(python_residue, 2);
                           PyObject *atoms_list = PyList_GetItem(python_residue, 3);

                           if (PyList_Size(atoms_list) > 0) {
                              // printf("atoms_list is a list\n");
                           } else {
                              printf("atoms_list is not a list\n");
                           }
                           int n_atoms = PyObject_Length(atoms_list);
                           if (n_atoms > 0) {
                              mmdb::Residue *residue_p = new mmdb::Residue;
                              std::string resname = PyBytes_AS_STRING(PyUnicode_AsUTF8String(python_residue_name));
                              int resno = PyLong_AsLong(python_residue_number);
                              // std::cout << "DEBUG:: Found resno:   " << resno << std::endl;
                              // std::cout << "DEBUG:: Found resname: " << resname << std::endl;
                              std::string inscode = PyBytes_AS_STRING(PyUnicode_AsUTF8String(python_residue_inscode));
                              residue_p->SetResName(resname.c_str());
                              residue_p->seqNum = resno;
                              memcpy(residue_p->insCode, inscode.c_str(), sizeof(mmdb::InsCode));
                              for (int iat=0; iat<n_atoms; iat++) {
                                 PyObject *atom_expression = PyList_GetItem(atoms_list, iat);
                                 int len_atom_expr = PyObject_Length(atom_expression);
                                 if (len_atom_expr != 3) {
                                    std::cout << "bad atom expression, length "
                                              << len_residue_expr << std::endl;
                                    const char *mess =  "object: %S\n";
                                    PyObject *bad_python = PyUnicode_FromFormat(mess, atom_expression);
                                    std::string bad_str = PyBytes_AS_STRING(PyUnicode_AsUTF8String(bad_python));
                                    std::cout << bad_str << std::endl;
                                 } else {
                                    // normal case
                                    // std::cout << "good atom expression " << std::endl;
                                    PyObject *name_alt_conf_pair = PyList_GetItem(atom_expression, 0);
                                    PyObject *occ_b_ele = PyList_GetItem(atom_expression, 1);
                                    PyObject *pos_expr = PyList_GetItem(atom_expression, 2);
                                    int len_name_alt_conf = PyObject_Length(name_alt_conf_pair);
                                    int len_occ_b_ele = PyObject_Length(occ_b_ele);
                                    int len_pos_expr = PyObject_Length(pos_expr);
                                    if (len_name_alt_conf == 2) {
                                       // the occ_b_ele list can contain an optional segid
				      if ( (len_occ_b_ele == 3) || (len_occ_b_ele == 4)) {
                                          if (len_pos_expr == 3) {
                                             PyObject *atom_name_python = PyList_GetItem(name_alt_conf_pair, 0);
                                             std::string atom_name = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_name_python));
                                             PyObject *alt_conf_python = PyList_GetItem(name_alt_conf_pair, 1);
                                             std::string alt_conf = PyBytes_AS_STRING(PyUnicode_AsUTF8String(alt_conf_python));
                                             PyObject *occ_python = PyList_GetItem(occ_b_ele, 0);
                                             PyObject *b_python = PyList_GetItem(occ_b_ele, 1);
                                             PyObject *ele_python = PyList_GetItem(occ_b_ele, 2);
					     std::string segid;
					     bool have_segid = 0;
					     if (len_occ_b_ele == 4) {
                                               PyObject *segid_python = PyList_GetItem(occ_b_ele, 3);
                                               if (PyUnicode_Check(segid_python)) { 
						 have_segid = 1;
						 segid = PyBytes_AS_STRING(PyUnicode_AsUTF8String(segid_python));
                                               } 
					     }
					     
                                             float occ = PyFloat_AsDouble(occ_python);
                                             std::string ele = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ele_python));
                                             float x = PyFloat_AsDouble(PyList_GetItem(pos_expr, 0));
                                             float y = PyFloat_AsDouble(PyList_GetItem(pos_expr, 1));
                                             float z = PyFloat_AsDouble(PyList_GetItem(pos_expr, 2));
                                             
                                             //parse the B-factor information
                                             bool have_aniso = 0;
                                             float b=0, b_u11=0, b_u22=0, b_u33=0, b_u12=0, b_u13=0, b_u23=0;
                                             PyObject *b_iso_python, *b_u11_python, *b_u22_python, *b_u33_python, *b_u12_python, *b_u13_python, *b_u23_python;
                                             // if (PyObject_TypeCheck(b_python, &PyFloat_Type) || PyObject_TypeCheck(b_python, &PyLong_Type) || PyObject_TypeCheck(b_python, &PyInt_Type)) {
                                             if (true) { // FIXME Python3
                                                //if the atom has an isotropic B-factor only
                                                b = PyFloat_AsDouble(b_python);
                                             } else if (PyObject_TypeCheck(b_python, &PyList_Type) && PyObject_Length(b_python) == 7) {
                                                //if the atom has anisotropic B-factor information
                                                have_aniso = 1;
                                                
                                                b_iso_python = PyList_GetItem(b_python, 0);
                                                b_u11_python = PyList_GetItem(b_python, 1);
                                                b_u22_python = PyList_GetItem(b_python, 2);
                                                b_u33_python = PyList_GetItem(b_python, 3);
                                                b_u12_python = PyList_GetItem(b_python, 4);
                                                b_u13_python = PyList_GetItem(b_python, 5);
                                                b_u23_python = PyList_GetItem(b_python, 6);
                                                
                                                b = PyFloat_AsDouble(b_iso_python);
                                                b_u11 = PyFloat_AsDouble(b_u11_python);
                                                b_u22 = PyFloat_AsDouble(b_u22_python);
                                                b_u33 = PyFloat_AsDouble(b_u33_python);
                                                b_u12 = PyFloat_AsDouble(b_u12_python);
                                                b_u13 = PyFloat_AsDouble(b_u13_python);
                                                b_u23 = PyFloat_AsDouble(b_u23_python);
                                             } else {
                                                b = -1.0;
                                                std::cout << "bad b-factor"
                                                       << std::endl;
                                                PyObject *bad_python = display_python(b_python);
                                                std::string bad_str = PyBytes_AS_STRING(PyUnicode_AsUTF8String(bad_python));
                                                std::cout << bad_str << std::endl;
                                             }
                                             
                                             mmdb::Atom *atom = new mmdb::Atom;
                                             atom->SetCoordinates(x, y, z, occ, b);
                                             atom->SetAtomName(atom_name.c_str());
                                             atom->SetElementName(ele.c_str());
                                             strncpy(atom->altLoc, alt_conf.c_str(), 2);
					     if (have_segid)
						strncpy(atom->segID, segid.c_str(), 5);
                                             if (have_aniso) {
                                                atom->u11 = b_u11;
                                                atom->u22 = b_u22;
                                                atom->u33 = b_u33;
                                                atom->u12 = b_u12;
                                                atom->u13 = b_u13;
                                                atom->u23 = b_u23;
                                                atom->WhatIsSet |= mmdb::ASET_Anis_tFac;
                                             }
                                             residue_p->AddAtom(atom);
					     
                                             // std::cout << "DEBUG:: adding atom " << atom << std::endl;
                                          } else {
                                             std::cout << "bad atom (position expression) "
                                                       << std::endl;
                                             PyObject *bad_python = display_python(pos_expr);
                                             std::string bad_str = PyBytes_AS_STRING(PyUnicode_AsUTF8String(bad_python));
                                             std::cout << bad_str << std::endl;
                                          }
                                       } else {
                                          std::cout << "bad atom (occ b element expression) "
                                                    << std::endl;
                                          PyObject *bad_python = display_python(occ_b_ele);
                                          std::string bad_str = PyBytes_AS_STRING(PyUnicode_AsUTF8String(bad_python));
                                          std::cout << bad_str << std::endl;
                                       }
                                    } else {
                                       std::cout << "bad atom (name alt-conf expression) "
                                                 << std::endl;
                                       PyObject *bad_python = display_python(name_alt_conf_pair);
                                       std::string bad_str = PyBytes_AS_STRING(PyUnicode_AsUTF8String(bad_python));
                                       std::cout << bad_str << std::endl;
                                    }
                                 }
                              }
                              chain_p->AddResidue(residue_p);
                           }
                        }
                     }
                     model_p->AddChain(chain_p);
                  }
               }
               model_list.push_back(model_p);
            }
         }
      }
   }
   
   return model_list;
}

#endif // USE_PYTHON

// delete CONECT records from mmdb manager
void mmdb_manager_delete_conect(mmdb::Manager *mol) {

  graphics_info_t g;

  if (g.write_conect_records_flag != 1) {
     mol->Delete ( mmdb::MMDBFCM_SC );
  }
}
