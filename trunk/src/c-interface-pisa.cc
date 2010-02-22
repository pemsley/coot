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


#include "mmdb_manager.h"
#include <gtk/gtk.h>

#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "c-interface-scm.hh"

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

	 CMMDBManager *mol1 = graphics_info_t::molecules[imol_1].atom_sel.mol;
	 CMMDBManager *mol2 = graphics_info_t::molecules[imol_2].atom_sel.mol;

	 coot::close_residues_from_different_molecules_t cr;
	 std::pair<std::vector<CResidue *>, std::vector<CResidue *> > res_pair = 
	    cr.close_residues(mol1, mol2, dist);

	 if (res_pair.first.size() > 0) { 
	    std::pair<bool, CMMDBManager *> nm =
	       coot::util::create_mmdbmanager_from_residue_vector(res_pair.first);
	    if (nm.second) {
	       int imol = graphics_info_t::create_molecule();
	       atom_selection_container_t asc = make_asc(nm.second);
	       std::string name = "interacting residues from ";
	       name += coot::util::int_to_string(imol_1);
	       graphics_info_t::molecules[imol].install_model(imol, asc, name, 1);
	       imodel_new = imol;
	    } else {
	       std::cout << "WARNING:: no molecule from create_mmdbmanager_from_residue_vector"
			 << std::endl;
	    } 
	 }
	 
	 if (res_pair.second.size() > 0) { 
	    std::pair<bool, CMMDBManager *> nm =
	       coot::util::create_mmdbmanager_from_residue_vector(res_pair.second);
	    if (nm.second) {
	       int imol = graphics_info_t::create_molecule();
	       atom_selection_container_t asc = make_asc(nm.second);
	       std::string name = "interacting residues from ";
	       name += coot::util::int_to_string(imol_2);
	       graphics_info_t::molecules[imol].install_model(imol, asc, name, 1);
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

   SCM interfaces_description_length_scm = scm_length(interfaces_description_scm);
   int interfaces_description_length = scm_to_int(interfaces_description_length_scm);
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
      if (interface_length == 2) { 
	 SCM molecule_pair_pair = scm_list_ref(interface_scm, SCM_MAKINUM(0));
	 SCM molecule_pair_pair_length_scm = scm_length(molecule_pair_pair);
	 int molecule_pair_pair_length = scm_to_int(molecule_pair_pair_length_scm);
	 if (molecule_pair_pair_length == 2) {
	    // 2 molecules in an interface (good) :-)

	    SCM molecule_pair_1 =  scm_list_ref(molecule_pair_pair, SCM_MAKINUM(0));
	    SCM molecule_pair_2 =  scm_list_ref(molecule_pair_pair, SCM_MAKINUM(1));

	    SCM molecule_number_1_scm = scm_list_ref(molecule_pair_1, SCM_MAKINUM(0));
	    SCM molecule_number_2_scm = scm_list_ref(molecule_pair_2, SCM_MAKINUM(0));

	    SCM molecule_1_record_scm = scm_list_ref(molecule_pair_1, SCM_MAKINUM(1));
	    SCM molecule_2_record_scm = scm_list_ref(molecule_pair_2, SCM_MAKINUM(1));
	    
	    SCM mol_1_residue_records = pisa_molecule_record_residues(molecule_1_record_scm);
	    SCM mol_2_residue_records = pisa_molecule_record_residues(molecule_2_record_scm);

	    SCM mol_1_chain_id_scm = pisa_molecule_record_chain_id(molecule_1_record_scm);
	    SCM mol_2_chain_id_scm = pisa_molecule_record_chain_id(molecule_2_record_scm);
	    
	    imol_1 = scm_to_int(molecule_number_1_scm);
	    imol_2 = scm_to_int(molecule_number_2_scm);

	    try { 

	       std::string chain_id_1 = scm_to_locale_string(mol_1_chain_id_scm);
	       std::string chain_id_2 = scm_to_locale_string(mol_2_chain_id_scm);
	    
	       std::vector<coot::residue_spec_t> mol_1_residue_specs = 
		  residue_records_list_scm_to_residue_specs(mol_1_residue_records, chain_id_1);
	       std::vector<coot::residue_spec_t> mol_2_residue_specs = 
		  residue_records_list_scm_to_residue_specs(mol_2_residue_records, chain_id_2);

	       make_complementary_dotted_surfaces(imol_1, imol_2,
						  mol_1_residue_specs,
						  mol_2_residue_specs);;
	    }
	    catch (std::runtime_error rte)  {
	       std::cout << "WARNING:: " << rte.what() << std::endl;
	    }
					       
	    SCM bonds_info_scm = scm_list_ref(interface_scm, SCM_MAKINUM(1));
// 	    std::cout << " in handle_pisa_interfaces_scm() got interface bonds: "
// 		      << scm_to_locale_string(display_scm(bonds_info_scm))
// 		      << std::endl;
	    if (scm_is_true(scm_list_p(bonds_info_scm))) {
	       SCM bonds_info_length_scm = scm_length(bonds_info_scm);
	       int bonds_info_length = scm_to_int(bonds_info_length_scm);
	       for (int ibond=0; ibond<bonds_info_length; ibond++) {
		  SCM pisa_bond_scm = scm_list_ref(bonds_info_scm, SCM_MAKINUM(ibond));
		  add_pisa_interface_bond_scm(imol_1, imol_2, pisa_bond_scm, i);
	       }
	    }
	 }
      }
   }
   return SCM_MAKINUM(-1);
}
#endif /* USE_GUILE */


#ifdef USE_GUILE
std::vector<coot::residue_spec_t> 
residue_records_list_scm_to_residue_specs(SCM mol_1_residues, const std::string &chain_id) {

   std::vector<coot::residue_spec_t> r;

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
	 catch (std::runtime_error rte) {
	    std::cout << "WARNING bad seq-num from pisa interfaces xml "
		      << scm_to_locale_string(display_scm(seq_num_scm)) << std::endl;
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
#endif

#ifdef USE_GUILE
SCM pisa_molecule_record_residues(SCM molecule_record_1) {
   std::string symbol = "residues";
   return symbol_value_from_record(molecule_record_1, symbol);
}
#endif
#ifdef USE_GUILE
SCM pisa_molecule_record_chain_id(SCM molecule_record_1) {
   std::string symbol = "chain-id";
   return symbol_value_from_record(molecule_record_1, symbol);
}
#endif

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

	 if (bond_type != "") 
	    add_generic_object_bond(imol_1, imol_2,
				    atom_spec_from_scm_expression(atom_spec_1),
				    atom_spec_from_scm_expression(atom_spec_2),
				    generic_object_number, colour);
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
	 CAtom *at1 = graphics_info_t::molecules[imol1].get_atom(atom_spec_1);
	 CAtom *at2 = graphics_info_t::molecules[imol2].get_atom(atom_spec_2);
	 if (! at1)
	    std::cout << "WARNING:: failed to get atom from spec " << atom_spec_1
		      << " in molecule " << imol1 << "\n";
	 if (! at2)
	    std::cout << "WARNING:: failed to get atom from spec " << atom_spec_2
		      << " in molecule " << imol2 << "\n";
	 if (at1 && at2) {
	    clipper::Coord_orth pt1(at1->x, at1->y, at1->z);
	    clipper::Coord_orth pt2(at2->x, at2->y, at2->z);
	    // std::cout << "    " << at1 <<  " and " << at2 << std::endl;
	    to_generic_object_add_dashed_line(generic_object_number,
					      colour.c_str(),
					      5, 0.5,
					      at1->x, at1->y, at1->z,
					      at2->x-at1->x, at2->y-at1->y, at2->z-at1->z);
	 }
      }
   }
}



void pisa_clear_interfaces() {
   

}


// If you want to do something with transparent surfaces in the
// future, then replace this function.
//
// This is called for each interface to be respresented.
//
// Perhaps should be std::pair<int, int>, returning the dot index
// associated with imol_1 and imol_2 (so that we have a way of turning
// them off)
// 
std::pair<int, int>
make_complementary_dotted_surfaces(int imol_1, int imol_2, 
				   std::vector<coot::residue_spec_t> &r1, 
				   std::vector<coot::residue_spec_t> &r2) {

   // symmetry is handled outside of here, which means that imol_2 can
   // be a symmetry copy of imol_1.
   
   // make synthetic molecules, dots where each residue contains one atom (dot).
   // then call std::pair<std::vector<CResidue *>, std::vector<CResidue *> >
   // coot::close_residues_from_different_molecules_t::close_residues(CMMDBManager *mol1,
   //                                                                 CMMDBManager *mol2,
   //                                                                 float dist)
   // 
   // consider making the dots generation a member function of dots_representation_info_t
   // and making it available here (if it isn't already).

   float close_dist = 4.5;

   if (is_valid_model_molecule(imol_1)) { 
      if (is_valid_model_molecule(imol_2)) {
	 CMMDBManager *mol_1 = graphics_info_t::molecules[imol_1].atom_sel.mol;
	 CMMDBManager *mol_2 = graphics_info_t::molecules[imol_2].atom_sel.mol;
	 
	 CMMDBManager *frag_mol_1 = coot::util::create_mmdbmanager_from_residue_specs(r1, mol_1);
	 CMMDBManager *frag_mol_2 = coot::util::create_mmdbmanager_from_residue_specs(r2, mol_2);

	 // coot::dots_representation_info_t d1(frag_mol_1, frag_mol_2);
	 // coot::dots_representation_info_t d2(frag_mol_2, frag_mol_1);

	 coot::dots_representation_info_t d1(frag_mol_1, frag_mol_2);
	 coot::dots_representation_info_t d2(frag_mol_2, frag_mol_1);

	 graphics_info_t::molecules[imol_1].add_dots(d1);
	 graphics_info_t::molecules[imol_2].add_dots(d2);

	 graphics_info_t::molecules[imol_1].set_dots_colour(0.6, 0.6, 0.3);
	 graphics_info_t::molecules[imol_2].set_dots_colour(0.6, 0.3, 0.6);

	 // need to delete those temporary molecules
	 delete frag_mol_1;
	 delete frag_mol_2;
      }
   }
   graphics_draw();
   return std::pair<int, int> (-1, -1);
} 

