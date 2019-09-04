/* src/testing.cc
 * 
 * Copyright 2008, 2009, 2010 by the University of Oxford
 * Copyright 2014, 2015 by Medical Research Council
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

#include "compat/coot-sysdep.h"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/use-rdkit.hh"
#endif 


#include "testing.hh"

#ifdef BUILT_IN_TESTING

#include <iostream>
#include <sstream>

#include <GL/glut.h> // needed for GLUT_ELAPSED_TIME

#include "clipper/core/ramachandran.h"
#include "clipper/ccp4/ccp4_map_io.h"

#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-rama.hh"
#include "ligand/primitive-chi-angles.hh"

#include "ligand/wligand.hh"
#include "ideal/simple-restraint.hh"
#include "ideal/torsion-bonds.hh"
#include "ligand/ligand.hh"
#include "ligand/chi-angles.hh"

#ifdef HAVE_GSL
#else
#include "utils/coot-utils.hh" // usually include from simple-restraint.hh
#endif

#include "coot-utils/coot-rama.hh"
#include "coot-utils/coot-shelx.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"

#include "ligand/dipole.hh"

#include "graphics-info.h"

// for file name globbing
#include "cc-interface.hh"
#include "coot-fileselections.h"

#include "c-interface-ligands.hh"

#include "c-interface.h"

#include "coot-utils/coot-h-bonds.hh"

#ifdef HAVE_CCP4SRS
#include <ccp4srs/ccp4srs_defs.h>
#endif

#include "pli/pi-stacking.hh"


bool close_float_p(float f1, float f2) {

   return (fabsf(f1-f2) < 0.0001);

} 

#include "testing-data.hh"

coot::protein_geometry testing_data::geom;


std::string greg_test(const std::string &file_name) {

   std::string dd;
   const char *c = getenv("COOT_TEST_DATA_DIR");
   if (c) {
      dd = c;
      dd += "/";
      dd += file_name;
   } else {
      const char *d = getenv("HOME");
      if (d) {
	 dd = d;
	 dd += "/data/greg-data/";
	 dd += file_name;
      }
   }
   return dd;
}

std::string stringify(double x) {
   std::ostringstream o;
   if (!(o << x))
      throw std::runtime_error("stringify(double)");
   return o.str();
}


#ifdef STRINGIFY_SIZE_T
std::string stringify(size_t x) {
   std::ostringstream o;
   if (!(o << x))
      throw std::runtime_error("stringify(size_t)");
   return o.str();
}
#endif

std::string stringify(int i) {
   std::ostringstream o;
   if (!(o << i))
      throw std::runtime_error("stringify(int)");
   return o.str();
}

std::string stringify(unsigned int i) {
   std::ostringstream o;
   if (!(o << i))
      throw std::runtime_error("stringify(unsigned int)");
   return o.str();
}

bool close_double(double a, double b) {
   double small_diff = 0.001;
   return (fabs(a-b) < small_diff);
}

bool close_pair_test(std::pair<double, double> a, std::pair<double, double> b) {
   return (close_double(a.first, b.first) && close_double(a.second, b.second));
}

std::ostream& operator<< (std::ostream &s, std::pair<double, double> b) {
   s << "[" << b.first << "," << b.second << "]";
   return s;
}

void add_test(int (*)(), const std::string &test_name, std::vector<named_func> *functions) {
   //    named_func p(t, test_name);
   //    functions->push_back(named_func(p));
} 

// these are not run by greg.  But they are run by python-tests/09_internal.py
//
int test_internal() {

   int status = 1;
   std::vector<named_func> functions;
   functions.push_back(named_func(kdc_torsion_test, "kevin's torsion test"));
   functions.push_back(named_func(test_alt_conf_rotamers, "test_alt_conf_rotamers"));

   // re-instate test_wiggly_ligands when you can get it to pass. 
   // functions.push_back(named_func(test_wiggly_ligands, "test_wiggly_ligands"));
   
   // functions.push_back(named_func(test_ramachandran_probabilities, "test_ramachandran_probabilities"));
   functions.push_back(named_func(test_fragmemt_atom_selection, "test_fragmemt_atom_selection"));
   functions.push_back(named_func(test_add_atom, "test_add_atom"));

   // restore me at some stage
   // functions.push_back(named_func(restr_res_vector, "restraints_for_residue_vec"));

   // restore me at some stage.
   // functions.push_back(named_func(test_peptide_link, "test_peptide_link"));


   // dictionaries don't have partial charges at the moment
   //
   // functions.push_back(named_func(test_dictionary_partial_charges,
   // 				  "test dictionary partial charges"));

   // restore me at some stage
   // functions.push_back(named_func(test_dipole, "test_dipole"));

   functions.push_back(named_func(test_segid_exchange, "test segid exchange"));

   functions.push_back(named_func(test_ligand_fit_from_given_point, "test ligand from point"));

   functions.push_back(named_func(test_peaksearch_non_close_peaks,
				  "test peak search non-close"));

   functions.push_back(named_func(test_symop_card, "test symop card"));
   functions.push_back(named_func(test_rotate_around_vector, "test rotate round vector"));
   functions.push_back(named_func(test_ssm_sequence_formatting, "SSM sequence alignment output"));
   status = run_internal_tests(functions);
   return status;
   
}

// Greg testing runs these tests - not the above
//
int greg_internal_tests() {
   int status = 1;
   std::vector<named_func> functions;
   functions.push_back(named_func(test_OXT_in_restraints, "OXT in restraints?"));
   functions.push_back(named_func(test_relativise_file_name, "Relative file name"));
   functions.push_back(named_func(test_geometry_distortion_info_type, "geometry distortion comparision"));
   functions.push_back(named_func(test_translate_close_to_origin, "test symm trans to origin"));
   functions.push_back(named_func(test_lsq_plane, "test lsq plane"));
   functions.push_back(named_func(test_COO_mod, "test COO modification"));
   functions.push_back(named_func(test_remove_whitespace, "remove whitespace"));
   functions.push_back(named_func(test_new_comp_id, "New comp_ids are sane"));
   functions.push_back(named_func(test_trailing_slash, "Remove Trailing Slash"));

   // restore this at some stage
   // functions.push_back(named_func(test_copy_cell_symm_orig_scale_headers, "test copy cell, symm, orig, scale cards"));

   status = run_internal_tests(functions);
   return status;
}

// greg runs these tests too, these tests use data from the greg test
// data directory.
// 
int greg_tests_using_external_data() {

   int status = 1;
   std::vector<named_func> functions;
   functions.push_back(named_func(test_phi_psi_values,
				  "Residues for phi,psi are close enough to be considered linked"));
   status = run_internal_tests(functions);
   return status;
}

// uses greg data test data
int test_phi_psi_values() {

   std::string filename = greg_test("frag-2wot.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);
   int n_phi_psi = 0;
   
   int status = 0;
   if (atom_sel.read_success > 0) {
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=1; ires<(nres-1); ires++) {
	    mmdb::Residue *prev_res = chain_p->GetResidue(ires-1);
	    mmdb::Residue *this_res = chain_p->GetResidue(ires);
	    mmdb::Residue *next_res = chain_p->GetResidue(ires+1);
	    try { 
	       coot::util::phi_psi_t pp(prev_res, this_res, next_res);
	       n_phi_psi++;
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << rte.what() << std::endl;
	    } 
	 }
      }
   }
   if (n_phi_psi == 5)
      status = 1;
   else
      std::cout << "   should have found 5 phi,psis - found " << n_phi_psi << std::endl;

   return status;
}

int test_internal_single() {
   int status = 0;
   try { 
      // status = test_symop_card();
      // status = test_rotate_round_vector();
      // status = test_coot_atom_tree();
      // status = test_coot_atom_tree_2();
      // status = test_coot_atom_tree_proline();
      // status = test_ssm_sequence_formatting();
      // status = test_previous_water();
      // status = test_sbase();
      // status = test_coordinated_waters();
      // status = test_geometry_distortion_info_type();
      // status = test_translate_close_to_origin();
      // status = test_flev_aromatics();
      // status = test_phi_psi_values();
      // status = test_map_segmentation();
      // status = test_lsq_plane();
      // status = test_map_segmentation();
      // status = test_copy_cell_symm_orig_scale_headers();
      // status = test_residue_atom_renaming();
      // status = test_residue_atom_renaming();
      // status = test_mcd_and_thornton_h_bonds();
      // status = test_COO_mod();
      // status = test_remove_whitespace();
      // status = test_beam_in_residue();
      // status = test_multi_residue_torsion();
      // status = test_torsions_from_residue_selection();
      // status = test_read_prosmart_distance_restraints();
      // status = test_dreiding_torsion_energy();
      // status = test_parallel_plane_restraints();
      // status = test_map_tools();
      // status = test_minimol();
      // status = test_monomer_organic_set();
      // status = test_COO_mod();
      // status = test_output_link_distances_are_correct();
      // status = test_string_splitting();
      status = test_index_splitting();
   }
   catch (const std::runtime_error &mess) {
      std::cout << "FAIL: " << " " << mess.what() << std::endl;
   }
   return status;
}

// return 1 on successful completion of tests.
int
run_internal_tests(std::vector<named_func> functions) {

   int status = 1;
   for (unsigned int i_func=0; i_func<functions.size(); i_func++) {
      std::cout << "Entering test: " << functions[i_func].second << std::endl;
      try { 
	 status = functions[i_func].first();
	 if (status) {
	    std::cout << "PASS: " << functions[i_func].second << std::endl;
	 } else { 
	    std::cout << "FAIL: " << functions[i_func].second << std::endl;
	    break;
	 }
      }
      catch (const std::runtime_error &mess) {
	 std::cout << "FAIL: " << functions[i_func].second << " " << mess.what() << std::endl;
	 status = 0;
	 break;
      }
   } 
   return status;
}


int test_minimol() {

   int status = 0;
   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);
   bool ifound = 0;

   // OK, now let's make a minimol
   if (atom_sel.read_success > 0) {
      coot::minimol::molecule m_basic(atom_sel.mol);

      coot::minimol::molecule m(m_basic[0]);

      coot::minimol::atom at1(atom_sel.atom_selection[0]);
      coot::minimol::atom at2(atom_sel.atom_selection[10]);
      m.fragments[0][-100].addatom(at1);
      m.fragments[0][ -99].addatom(at2);
      

      // test for baddies

      bool found_bad = false;
      for(unsigned int ifrag=0; ifrag<m.fragments.size(); ifrag++) {
	 for(int ires=m[ifrag].min_res_no(); ires<=m[ifrag].max_residue_number(); ires++) {
	    if (m[ifrag][ires].seqnum < -100) {
	       std::cout << "  Baddie! res-idx " << ires << " "
			 << m[ifrag][ires].seqnum << std::endl;
	       found_bad = true;
	       break;
	    } 
	 }
      }
      // m.write_file("test-minimol.pdb", 20);

      if (! found_bad)
	 status = 1; // good
   }

   std::cout << "print test_minimol returns " << status << std::endl;
   return status;
}

int test_monomer_organic_set() {

   int status = 0;
   testing_data t;
   std::string types[4] = { "TYR", "CYS", "ATP", "OS4" };
   int read_number = 40;
   for (int i=0; i<4; i++) { 
      std::string res_type = types[i];
      std::pair<bool, coot::dictionary_residue_restraints_t> rp = 
         t.geom.get_monomer_restraints(res_type, 0);
      if (! rp.first) { 
         t.geom.try_dynamic_add(res_type,read_number++);
      }
      int imol = 0; // dummy
      if (t.geom.have_dictionary_for_residue_type(res_type, imol, read_number++)) { 
         bool f = rp.second.comprised_of_organic_set();
         if (f) 
            std::cout << "test: " << res_type << " is IN organic set" << std::endl;
         else
            std::cout << "test: " << res_type << " is NOT in organic set" << std::endl;
      } else { 
         std::cout << "test: " << res_type << " -- no dictionary " << std::endl;
      } 
   }
   return status;
}

int test_output_link_distances_are_correct() {

   int status = 1;

#ifdef MMDB_HAS_LINK_DISTANCE

   status = 0;

   std::string filename = greg_test("pdb4rdq.ent");
   if (coot::file_exists("pdb4rdq.ent")) {
      atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);
      if (atom_sel.mol) {
	 mmdb::Residue *r = test_get_residue(atom_sel.mol, "E", 502);
	 if (! r) {
	    std::cout << "test_output_link_distances_are_correct():::: No residue!!! " << std::endl;
	 } else { 
	    int n_atoms = r->GetNumberOfAtoms();
	    for (int iat=0; iat<n_atoms; iat++) {
	       mmdb::Atom *at = r->GetAtom(iat);
	       at->x += 3;
	       at->z += 3;
	       // std::cout << "moving " << at << std::endl;
	    }
	    atom_sel.mol->FinishStructEdit();
	    filename = "pdb4rqd-with-moved-CA.pdb";
	    coot::write_coords_pdb(atom_sel.mol, filename);
	    if (coot::file_exists(filename)) {
	       atom_sel = get_atom_selection(filename, true, true);
	       if (atom_sel.mol) {
		  mmdb::Model *model_p = atom_sel.mol->GetModel(1);
		  int n_links = model_p->GetNumberOfLinks();
		  if (n_links > 0) {
		     status = 1; // all OK so far
		     for (int i_link=1; i_link<=n_links; i_link++) {
			mmdb::Link *link = model_p->GetLink(i_link);
			std::pair<coot::atom_spec_t, coot::atom_spec_t> lp = coot::link_atoms(link, model_p);
			mmdb::Atom *at_1 = coot::util::get_atom(lp.first,  atom_sel.mol);
			mmdb::Atom *at_2 = coot::util::get_atom(lp.second, atom_sel.mol);
			if (at_1) {
			   if (at_2) {
			      double link_dist = link->dist;
			      double atom_dist = coot::distance(at_1, at_2);
			      // std::cout << "LINK " << i_link << " check read " << link_dist << std::endl;
			      double d = fabs(link_dist - atom_dist);
			      if (d > 0.01) {
				 status = 0; // a bad LINK distance!
				 std::cout << i_link << " LINK " << link << "  dist " << link_dist
					   << " but atom dist " << atom_dist << std::endl;
				 
				 break;
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      } 
   }
#endif
   return status;
} 


      
int test_alt_conf_rotamers() {

   int status = 1;

   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);
   bool ifound = 0;

   int imod = 1;
   if (atom_sel.read_success > 0) { 
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id = chain_p->GetChainID();
	 if (chain_id == "B") {
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::PResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int resno = residue_p->GetSeqNum();
	       if (resno == 72) {

		  ifound = 1;
		  coot::primitive_chi_angles prim_chis(residue_p);
		  std::vector<coot::alt_confed_chi_angles> chis = prim_chis.get_chi_angles();
		  if (chis.size() != 2) {
		     std::string mess = "chis.size() is ";
		     mess += stringify(int(chis.size()));
		     throw std::runtime_error(mess);
		  }

		  int n_rots_found = 0;
		  for (unsigned int i_rot=0; i_rot<2; i_rot++) {
// 		     std::cout << "DEBUG:: chi: " << chis[i_rot].alt_conf << " " 
// 			       << chis[i_rot].chi_angles[0].first  << " "
// 			       << chis[i_rot].chi_angles[0].second << std::endl;
		     
		     if (chis[i_rot].alt_conf == "A") {
			float chi = chis[i_rot].chi_angles[0].second;
			if (chi > 60.0)
			   if (chi < 61.0)
			      n_rots_found++;
		     }
		     if (chis[i_rot].alt_conf == "B") {
			float chi = chis[i_rot].chi_angles[0].second;
			if (chi > -75.0)
			   if (chi < -74.0)
			      n_rots_found++;
		     }
		  }
		  if (n_rots_found != 2) {
		     std::string mess = " found only ";
		     mess += stringify(n_rots_found);
		     mess += " rotamers ";
		     throw std::runtime_error(mess);
		  }
	       }

	       if (resno == 80) {
		  coot::primitive_chi_angles prim_chis(residue_p);
		  std::vector<coot::alt_confed_chi_angles> chis = prim_chis.get_chi_angles();
		  if (chis.size() != 1) {
		     std::string mess = "chis.size() is ";
		     mess += stringify(int(chis.size()));
		     mess += " for resno ";
		     mess += stringify(80);
		     throw std::runtime_error(mess);
		  }
		  int n_rots_found = 0;
		  if (chis[0].alt_conf == "") {
		     float chi_1 = chis[0].chi_angles[0].second;
		     float chi_2 = chis[0].chi_angles[1].second;
		     if (chi_1 > -58.0)
			if (chi_1 < -57.0)
			   if (chi_2 > 95.0)
			      if (chi_2 < 96.0)
			         n_rots_found++;
		  }
		  if (n_rots_found != 1) {
		     std::string mess = " Oops found ";
		     mess += stringify(n_rots_found);
		     mess += " rotamers ";
		     mess += " for resno ";
		     mess += stringify(80);
		     throw std::runtime_error(mess);
		  }
	       } 
	    }
	 }
      }
   }

   if (! ifound) {
      std::string mess = "file not found ";
      mess += filename;
      throw std::runtime_error(mess);
   }
   atom_sel.clear_up();

   return status; 
}

int test_wiggly_ligands () {

   // This is not a proper test.  It just exercises some functions.  I
   // had set wiggly_ligand_n_samples to 100 and checked that the
   // distribution of chi angles looked reasonable (it did).
   // Analysing the chi angles can be done another day.
   //
   // BUA looks great.  So why does 3GP's wiggly ligand look terrible?
   // Not sure, but I think it is because there are torsions going on
   // in the ribose (torsionable bonds (not const - they are filtered
   // out) that shouldn't be - and the ring is breaking).
   
   int r = 1;
   std::string cif_file_name = greg_test("libcheck_BUA.cif");
   coot::protein_geometry geom;
   coot::read_refmac_mon_lib_info_t rmit = geom.init_refmac_mon_lib(cif_file_name, 0);
   if (rmit.n_bonds == 0) {
      std::string m = "Critical cif dictionary reading failure.";
      std::cout << m << std::endl;
      throw std::runtime_error(m);
   }
   coot::wligand wlig;
   coot::minimol::molecule mmol;
   mmol.read_file(greg_test("monomer-BUA.pdb"));
   unsigned int wiggly_ligand_n_samples = 10;
   try { 
      bool optim_geom = 0;
      bool fill_vec = 1;
      int imol = 0; // dummy
      std::vector<coot::installed_wiggly_ligand_info_t> ms = 
	 wlig.install_simple_wiggly_ligands(&geom, mmol, imol, wiggly_ligand_n_samples,
					    optim_geom, fill_vec);

      // jump out if no returned molecules
      if (ms.size() != wiggly_ligand_n_samples) {
	 std::cout << "FAIL: ms.size() != wiggly_ligand_n_samples " << ms.size() << " "
		   << wiggly_ligand_n_samples << std::endl;
	 return 0;
      }

      //
      for (unsigned int imol=0; imol<ms.size(); imol++) {
	 std::string file_name = "test-wiggly-ligand-";
	 file_name += stringify(imol);
	 file_name += ".pdb";
	 ms[imol].mol.write_file(file_name, 10.0);
      }
   }
   catch (const std::runtime_error &mess) {
      std::cout << mess.what() << std::endl;
   } 
   return r;
}

// Return a new mol, and a residue selection.  Delete the residue
// selection and mol when you are done with it.
// 
residue_selection_t
testing_func_probabilities_refine_fragment(atom_selection_container_t atom_sel,
					   mmdb::PResidue *SelResidues,
					   int nSelResidues,
					   const std::string &chain_id,
					   int resno_mid,
					   coot::protein_geometry geom,
					   bool enable_rama_refinement,
					   int side_step,
					   bool use_flanking_residues,
					   bool output_numerical_gradients) {

#ifdef HAVE_GSL
   long t0 = glutGet(GLUT_ELAPSED_TIME);
   
   // now refine a bit of structure:
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   short int have_flanking_residue_at_start = 0;
   short int have_flanking_residue_at_end = 0;
   if (use_flanking_residues) {
       have_flanking_residue_at_start = 1;
       have_flanking_residue_at_end = 1;
   }
   short int have_disulfide_residues = 0;  // other residues are included in the
   std::string altconf = "";
   short int in_alt_conf_split_flag = 0;
   const char *chn = chain_id.c_str(); // mmdb thing.  Needs updating on new mmdb?
	       
   std::pair<mmdb::Manager *, int> residues_mol_pair = 
      coot::util::create_mmdbmanager_from_res_selection(atom_sel.mol,
							SelResidues, nSelResidues, 
							have_flanking_residue_at_start,
							have_flanking_residue_at_end,
							altconf,
							chain_id,
							in_alt_conf_split_flag);
   
   clipper::Xmap<float> dummy_xmap;
   coot::restraints_container_t restraints(resno_mid-side_step,
					   resno_mid+side_step,
					   have_flanking_residue_at_start,
					   have_flanking_residue_at_end,
					   have_disulfide_residues,
					   altconf,
					   chn,
					   residues_mol_pair.first,
					   fixed_atom_specs, &dummy_xmap);

   ctpl::thread_pool thread_pool(2);
   restraints.thread_pool(&thread_pool, 2);

   short int do_rama_restraints = 0;
   short int do_residue_internal_torsions = 1;
   short int do_link_torsions = 0;
   float rama_plot_restraint_weight = 1.0;

   coot::restraint_usage_Flags flags =
      coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;

   if (enable_rama_refinement) { 
      do_rama_restraints = 1;
      flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
      //      flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
      //flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES;
      //flags = coot::BONDS_ANGLES_TORSIONS_PLANES_AND_NON_BONDED;
      //flags = coot::BONDS_ANGLES_TORSIONS_PLANES_AND_CHIRALS;
      //flags = coot::BONDS_AND_NON_BONDED;
      //flags = coot::RAMA;
   } 


   coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
   bool do_trans_peptide_restraints = false;
   int imol = 0; // dummy
   int nrestraints = 
      restraints.make_restraints(imol, geom, flags,
				 do_residue_internal_torsions,
				 do_trans_peptide_restraints,
				 rama_plot_restraint_weight,
				 do_rama_restraints, false, false,
				 pseudos);

   if (output_numerical_gradients)
      restraints.set_do_numerical_gradients();
   restraints.minimize(flags);

   int post_refine_selHnd = residues_mol_pair.first->NewSelection();
   int post_refine_nSelResidues; 
   mmdb::PResidue *post_refine_SelResidues = NULL;
   residues_mol_pair.first->Select(post_refine_selHnd, mmdb::STYPE_RESIDUE, 0,
				   chn,
				   resno_mid-side_step, "",
				   resno_mid+side_step, "",
				   "*",  // residue name
				   "*",  // Residue must contain this atom name?
				   "*",  // Residue must contain this Element?
				   "",   // altLocs
				   mmdb::SKEY_NEW // selection key
				   );
   residues_mol_pair.first->GetSelIndex(post_refine_selHnd,
					post_refine_SelResidues,
					post_refine_nSelResidues);

   residue_selection_t res_sel;
   res_sel.mol = residues_mol_pair.first;
   res_sel.SelectionHandle = post_refine_selHnd;
   res_sel.nSelResidues = post_refine_nSelResidues;
   res_sel.SelResidues = post_refine_SelResidues;

   long t1 = glutGet(GLUT_ELAPSED_TIME);
   float seconds = float(t1-t0)/1000.0;
   std::cout << "refinement_took " << seconds << " seconds" << std::endl;
   return res_sel;
#endif // HAVE_GSL
} 

int test_ramachandran_probabilities() {

   int r = 0;

   std::string file_name = greg_test("crashes_on_cootaneering.pdb");
   file_name = "37-41.pdb";
   atom_selection_container_t atom_sel = get_atom_selection(file_name, true, true);

   if (! atom_sel.read_success)
      throw std::runtime_error(file_name + ": file not found.");


   std::string chain_id = "A";
   char *chn = (char *) chain_id.c_str(); // mmdb thing.  Needs updating on new mmdb?
   std::vector<int> resnos;
   resnos.push_back(39);

   coot::protein_geometry geom;
   geom.init_standard();
   int n_correct = 0; 
   for (unsigned int i=0; i<resnos.size(); i++) { 
      int selHnd = atom_sel.mol->NewSelection();
      int nSelResidues; 
      mmdb::PResidue *SelResidues = NULL;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
			   chn,
			   resnos[i]-2, "",
			   resnos[i]+2, "",
			   "*",  // residue name
			   "*",  // Residue must contain this atom name?
			   "*",  // Residue must contain this Element?
			   "",   // altLocs
			   mmdb::SKEY_NEW // selection key
			   );
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

      // get the probability
      std::string residue_type = SelResidues[1]->GetResName();

      clipper::Ramachandran::TYPE rama_type = clipper::Ramachandran::NonGlyPro5;
      if (residue_type == "GLY")
	 rama_type = clipper::Ramachandran::Gly5;
      if (residue_type == "PRO")
	 rama_type = clipper::Ramachandran::Pro5;
      clipper::Ramachandran rama(rama_type);

      double prob = 0;

//       coot::util::phi_psi_t angles = coot::util::ramachandran_angles(SelResidues, nSelResidues);
//       rama.probability(clipper::Util::d2rad(angles.phi()),
// 		       clipper::Util::d2rad(angles.psi()));
      
      double post_refine_prob = 0.0; // set later

      // --------------------------------------------------
      //     non-Ramachandran refinement:
      // --------------------------------------------------
      if (0) { 
	 int enable_rama_refinement = 0;
	 residue_selection_t refined_res_sel =
	    testing_func_probabilities_refine_fragment(atom_sel, SelResidues,
						       nSelResidues,
						       chain_id, resnos[i], geom,
						       enable_rama_refinement,
						       1, 0, 1);

	 if (0) { 
	    // Let's look at the refined structure. Write them out as pdb files ;-/
	    std::string tmp_file_name = "rama-test-";
	    tmp_file_name += coot::util::int_to_string(i);
	    tmp_file_name += ".pdb";
	    refined_res_sel.mol->WritePDBASCII((char *)tmp_file_name.c_str());
	 }

	 coot::util::phi_psi_t post_refine_angles =
	    coot::util::ramachandran_angles(refined_res_sel.SelResidues,
					    refined_res_sel.nSelResidues);
	 refined_res_sel.clear_up();
	       
	 post_refine_prob =
	    rama.probability(clipper::Util::d2rad(post_refine_angles.phi()),
			     clipper::Util::d2rad(post_refine_angles.psi()));
      }

	       
      // --------------------------------------------------
      //     now with Ramachandran refinement:
      // --------------------------------------------------


      int enable_rama_refinement = 1;
      bool flank = 1;
      residue_selection_t rama_refined_res_sel =
	 testing_func_probabilities_refine_fragment(atom_sel, SelResidues,
						    nSelResidues,
						    chain_id, resnos[i], geom,
						    enable_rama_refinement, 1, flank, 1);
      coot::util::phi_psi_t rama_refine_angles =
	 coot::util::ramachandran_angles(rama_refined_res_sel.SelResidues,
					 rama_refined_res_sel.nSelResidues);
      rama_refined_res_sel.clear_up();
	       
      double rama_refine_prob =
	 rama.probability(clipper::Util::d2rad(rama_refine_angles.phi()),
			  clipper::Util::d2rad(rama_refine_angles.psi()));
      std::cout << "--------------------------------------\n";
      std::cout << "Pre-refine         Rama probability residue " << resnos[i] << ": "
		<< prob << std::endl;
      std::cout << "Post-simple refine Rama probability residue " << resnos[i] << ": "
		<< post_refine_prob << std::endl;
      std::cout << "Post-Rama   refine Rama probability residue " << resnos[i] << ": "
		<< rama_refine_prob << std::endl;
      std::cout << "--------------------------------------\n";

   }
   return r;
} 


int kdc_torsion_test() { 
  int r = 1;

#ifdef HAVE_GSL
  clipper::Coord_orth co1, co2, co3, co4;
  std::vector<double> params(12);
  const int n1(5), n2(7);
  const double dx = 1.0e-3;
  for ( int t1 = 0; t1 < n1; t1++ )
     for ( int t2 = 0; t2 < n2; t2++ ) {

	// std::cout << " ===== t1 = " << t1 << "   t2 = " << t2 << " ======" << std::endl;
	
	// set up angles at either end of torsion
	double p1 = 6.283 * double(t1) / double(n1);
	double p2 = 6.283 * double(t2) / double(n2);
	// set up atomic coords up z axis
	params[0]  = cos(p1); params[1]  = sin(p1); params[2]  = 0.0;
	params[3]  =     0.0; params[4]  =     0.0; params[5]  = 1.0;
	params[6]  =     0.0; params[7]  =     0.0; params[8]  = 2.0;
	params[9]  = cos(p2); params[10] = sin(p2); params[11] = 3.0;

	// now pertub the coords (params[12] is unperturbed)
	std::vector<std::vector<double> > dparams(13,params);
	std::vector<double> dresult(13), ngrad(12), agrad(12);
	for ( int i = 0; i < 12; i++ )
	   dparams[i][i] += dx;

	// calculate torsion and derivs numerically
	std::vector<clipper::Coord_orth> co(4);
	for (unsigned int i = 0; i < dparams.size(); i++ ) {
	   const std::vector<double>& param = dparams[i];
	   // convert parameter list to coord orths
	   for (int j=0;j<4;j++)
	      for(int k=0;k<3;k++)
		 co[j][k]=param[3*j+k];
	   // calculate torsion using clipper
	   dresult[i] = clipper::Coord_orth::torsion( co[0], co[1], co[2], co[3] );
	}
	for ( int i = 0; i < 12; i++ ) {
// 	   std::cout << "i = " << i << " (" << dresult[i] << " - " << dresult[12]
// 		     << ")/" << dx << " = " << (dresult[i]-dresult[12])/dx << std::endl;
	   ngrad[i] = (dresult[i]-dresult[12])/dx;
	}
	// push the perturbation the other way
	for ( int i = 0; i < 12; i++ )
	   dparams[i][i] -= 2*dx; 
	for (unsigned int i = 0; i < dparams.size(); i++ ) {
	   const std::vector<double>& param = dparams[i];
	   // convert parameter list to coord orths
	   for (int j=0;j<4;j++)
	      for(int k=0;k<3;k++)
		 co[j][k]=param[3*j+k];
	   // calculate torsion using clipper
	   dresult[i] = clipper::Coord_orth::torsion( co[0], co[1], co[2], co[3] );
	}
	for ( int i = 0; i < 12; i++ ) {
	   ngrad[i] -= (dresult[i]-dresult[12])/dx;
	   ngrad[i] *= 0.5; // average + and - shift
	}

      
	// first check clipper torsion calc
	double tor1 = dresult[12];
	double tor2 = p2-p1;
	if ( cos( tor1 - tor2 ) < 0.999999 ) {
	   std::cout <<"TORSION ERROR (CLIPPER)"<<tor1<<" != "<<tor2<<std::endl;
	   r = 0;
	}

      
	// now check coot torsion calc
	for (int j=0;j<4;j++) for(int k=0;k<3;k++)
	   co[j][k]=params[3*j+k];
	coot::distortion_torsion_gradients_t dtg =
	   coot::fill_distortion_torsion_gradients( co[0], co[1], co[2], co[3] );
	double tor3 = clipper::Util::d2rad(dtg.theta);
	if ( cos( tor1 - tor3 ) < 0.999999 ) {
	   std::cerr <<"TORSION ERROR (COOT)  "<<tor1<<" != "<<tor3<<std::endl;
	   r = 0;
	}

	// ts = torsion_scale
	double ts = (1.0/(1+pow(tan(clipper::Util::d2rad(dtg.theta)),2)));

	// now fetch coot torsion gradients
	agrad[0] = ts*dtg.dD_dxP1; agrad[1] = ts*dtg.dD_dyP1; agrad[2] = ts*dtg.dD_dzP1;
	agrad[3] = ts*dtg.dD_dxP2; agrad[4] = ts*dtg.dD_dyP2; agrad[5] = ts*dtg.dD_dzP2;
	agrad[6] = ts*dtg.dD_dxP3; agrad[7] = ts*dtg.dD_dyP3; agrad[8] = ts*dtg.dD_dzP3;
	agrad[9] = ts*dtg.dD_dxP4; agrad[10]= ts*dtg.dD_dyP4; agrad[11]= ts*dtg.dD_dzP4;
	for ( int i = 0; i < 12; i++ )
	   if ( fabs( ngrad[i] - agrad[i] ) > 1.0e-4 ) {
	      char xyz[] = "xyz";
	      std::cerr << "TORSIONS " << i/3 << xyz[i%3] << " "
			<< ngrad[i] << " vs " << agrad[i] << std::endl;
	      r = 0;
	   }
     }
#endif // HAVE_GSL
  return r;
}


int
test_fragmemt_atom_selection() {

   int status = 0;

   bool fill_masking_molecule_flag = 1;
   std::string atom_selection_string = "//A,B/1-5";
   int n_atoms_in_frag = 64;
   // from wc, there are 64   atoms in that atom selection in tutorial-modern.pdb
   //          there are 1465 atoms in tutorial-modern.pdb
   
   std::string f = greg_test("tutorial-modern.pdb");
   atom_selection_container_t asc = get_atom_selection(f, true, true);
   
   std::pair<coot::minimol::molecule, coot::minimol::molecule> p = 
      coot::make_mols_from_atom_selection_string(asc.mol, atom_selection_string,
						 fill_masking_molecule_flag);

   // now test the number of atoms in first and second
   int n_initial = asc.n_selected_atoms;
   int n_1 = p.first.count_atoms();
   int n_2 = p.second.count_atoms();
   

   std::cout << "   n_initial: " << n_initial << "   n_1: " << n_1 << "   n_2: "
	     << n_2 << std::endl;
   if (n_1 == (n_initial - n_atoms_in_frag))
      if (n_2 == n_atoms_in_frag)
	 status = 1;
      
   return status;
}

int test_peptide_link() {
   
   std::string f = "1h4p.pdb";
   atom_selection_container_t asc = get_atom_selection(greg_test(f), true, true);
   if (! asc.read_success)
      return 0;

   std::vector<std::pair<bool,mmdb::Residue *> > residues;
   mmdb::Manager *mol = asc.mol;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      std::string chain_id = chain_p->GetChainID();
      if (chain_id == "B") { 
	 mmdb::PResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int resno = residue_p->GetSeqNum();
	    if ((resno == 1455) || (resno == 1456))
	       residues.push_back(std::pair<bool, mmdb::Residue *>(0, residue_p));
	 }
      }
   }

   if (residues.size() != 2) {
      return 0; 
   } 

   // ======== Now start the test for real ===============
   //
   coot::protein_geometry geom;
   geom.init_standard();

   try {
      std::string comp_id_1 = "MAN";
      std::string comp_id_2 = "MAN";
      std::string group_1 = "D-pyranose"; // CCD and acedrg dictionaries now use D-SACCHARIDE
      std::string group_2 = "D-pyranose";

      float weight = 1.0;
      std::vector<coot::atom_spec_t> fixed_atom_specs;
      std::vector<mmdb::Link> links;
      clipper::Xmap<float> dummy_xmap;

      coot::restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, &dummy_xmap);
      restraints.add_map(weight);
      std::string link_type = "";
      // restraints.find_link_type(residues[0].second,
      // 		residues[1].second,
      // 		geom);

      std::cout << "   link_type: " << link_type << ":" << std::endl;
      
      std::vector<std::pair<coot::chem_link, bool> > link_infos =
	 geom.matching_chem_link(comp_id_1, group_1, comp_id_2, group_2);
      unsigned int ilink = 0;
      std::cout << "Found link :" << link_infos[ilink].first.Id() << ":" << std::endl;
      if (link_infos[ilink].first.Id() != "BETA1-3") {
	 return 0;
      } 
   }

   catch (const std::runtime_error &mess) {
      std::cout << "in test_peptide_link() matching_chem_link fails" << std::endl;
   }

   return 1;

} 

 
// Restraints are correctly generated for vector of residues (new
// restraints constructor needed) 20081106
int
restr_res_vector() {

   std::string f = greg_test("tutorial-modern.pdb");
   //    f = "7_and_96_B-a-result.pdb";
//    f = "6_7_and_96_B.pdb";
//    f = "6_7.pdb";
   atom_selection_container_t asc = get_atom_selection(f, true, true);

   std::vector<std::pair<bool,mmdb::Residue *> > residues;
   mmdb::Manager *mol = asc.mol;
   std::cout << "restr_res_vector: mol: " << mol << std::endl;
   std::vector<coot::atom_spec_t> fixed_atom_specs;

   if (!asc.read_success)
      return 0;
   
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      std::string chain_id = chain_p->GetChainID();
      if (chain_id == "B") { 
	 mmdb::PResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int resno = residue_p->GetSeqNum();
	    if ((resno == 7) || (resno == 96))
	       residues.push_back(std::pair<bool, mmdb::Residue *>(0, residue_p));
	 }
      }
   }

   if (residues.size() != 2) {
      std::cout << "  Fail to find residues - found " << residues.size() << std::endl;
      return 0;
   } else {
      clipper::Xmap<float> xmap;
      coot::util::map_fill_from_mtz(&xmap, "rnasa-1.8-all_refmac1.mtz", "FWT", "PHWT", "WT", 0, 0);
      float weight = 1;
      coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      coot::protein_geometry geom;
      geom.init_standard();
      std::vector<mmdb::Link> links;
      coot::restraints_container_t
	 restraints(residues, links, geom, mol, fixed_atom_specs, &xmap);
      restraints.add_map(weight);
      bool do_trans_peptide_restraints = true;
      int imol = 0;
      restraints.make_restraints(imol, geom, flags, 0, do_trans_peptide_restraints, 0.0, 0, false, false, coot::NO_PSEUDO_BONDS);
      restraints.minimize(flags);
      restraints.write_new_atoms("ss-test.pdb");
   }
   return 0;
} 

 
int
test_add_atom() {

   int status = 0;

   std::string f = greg_test("tutorial-modern.pdb");
   atom_selection_container_t asc = get_atom_selection(f, true, true);

   int n_test_residues = 20;
   int pass_count = 0;
   int ires_count = 0;
   
   int imod = 1;
   mmdb::Model *model_p = asc.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::PResidue residue_p;
      mmdb::Atom *at;
      while (ires_count<n_test_residues) {
	 residue_p = chain_p->GetResidue(ires_count);
	 ires_count++;
	 std::string res_name(residue_p->GetResName());
	 if (res_name == "GLY" ||
	     res_name == "ALA" ) {
	    pass_count++; // automatic pass
	 } else {
	    // normal case
	    coot::chi_angles chi(residue_p, 0);
	    std::vector<std::pair<int, float> > chi_angles = chi.get_chi_angles();
	    if (chi_angles.size() == 0) {
	       std::cout << "   Failed to find chi angles in residue "
			 << coot::residue_spec_t(residue_p) << std::endl;
	    } else {
	       float bond_length = 1.53;
	       float angle = 113.8;
	       std::string ref_atom_name(" CG ");
	       if (res_name == "CYS") { 
		  ref_atom_name = " SG ";
		  bond_length = 1.808;
	       } 
	       if (res_name == "THR") { 
		  ref_atom_name = " OG1";
		  bond_length = 1.433;
		  angle = 109.600;
	       } 
	       if (res_name == "VAL")
		  ref_atom_name = " CG1";
	       if (res_name == "SER") { 
		  ref_atom_name = " OG ";
		  bond_length = 1.417;
	       }
	       if (res_name == "PRO") { 
		  bond_length = 1.492;
		  angle = 104.500;
	       }
	       mmdb::Atom *ref_atom = residue_p->GetAtom(ref_atom_name.c_str());
	       if (!ref_atom) {
		  std::cout << "   Failed to find reference CG in residue "
			    << coot::residue_spec_t(residue_p) << std::endl;
	       } else { 
		  double tors = chi_angles[0].second;
		  bool status = coot::util::add_atom(residue_p, " N  ", " CA ", " CB ", "",
						     bond_length,
						     angle, 
						     tors,
						     " XX ", " C", 1.0, 20.0);
		  if (status) {
		     clipper::Coord_orth ref_pt(ref_atom->x, ref_atom->y, ref_atom->z);
		     mmdb::Atom *new_atom = residue_p->GetAtom(" XX ");
		     if (! new_atom) {
			std::cout << "   Failed to find reference CG in residue "
				  << coot::residue_spec_t(residue_p) << std::endl;
		     } else {
			clipper::Coord_orth new_pos(new_atom->x, new_atom->y, new_atom->z);
			double d = clipper::Coord_orth::length(ref_pt, new_pos);
			if (d < 0.12) {
			   // std::cout << "   Pass make atom " << std::endl;
			   pass_count++; 
			} else {
			   std::cout << "   Failed closeness test, d = " << d << " "
				     << new_atom << " vs " << ref_atom << std::endl;
			}
		     }
		  } else {
		     std::cout << "   Failed to add atom to residue "
			       << coot::residue_spec_t(residue_p) << std::endl;
		  }
	       }
	    } 
	 } 
      }
   }
   if (pass_count == n_test_residues)
      status = 1;
   
   return status;
}

int
test_dipole() {

   int result = 0; // fail initially.
   testing_data t;
   
   std::string res_type = "TYR";
   
   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);

   std::pair<short int, coot::dictionary_residue_restraints_t> rp = 
      t.geom.get_monomer_restraints(res_type, 0);

   if (rp.first) { 
   
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::PResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    if (std::string(residue_p->GetResName()) == res_type) { 
	       coot::dipole d(rp.second, residue_p);
	       std::cout << "residue " << coot::residue_spec_t(residue_p)
			 << "   dipole: " << d << " at " << d.position().format()
			 << std::endl;
	       break;
	    }
	 }
      }
   }
   return result;
}


int
test_dictionary_partial_charges() {

   std::vector<std::string> v;
   v.push_back("ALA");
   v.push_back("ARG");
   v.push_back("TYR");
   v.push_back("TRP");
   v.push_back("VAL");
   v.push_back("SER");

   testing_data t;
   for (unsigned int iv=0; iv<v.size(); iv++) {
      std::pair<short int, coot::dictionary_residue_restraints_t> rp = 
	 t.geom.get_monomer_restraints(v[iv], 0);
      if (! rp.first) {
	 std::cout << " Fail - no restraints for " << v[iv] << std::endl;
	 return 0;
      } else {
	 for (unsigned int iat=0; iat<rp.second.atom_info.size(); iat++) {
	    if (! rp.second.atom_info[iat].partial_charge.first) {
	       std::cout << " Fail - no partial charge for "
			 << rp.second.atom_info[iat].atom_id << " in "
			 << v[iv] << std::endl;
	       return 0;
	    }
	 }
      } 
   }
   return 1;
}

int test_segid_exchange() {

   int status = 0;
   
   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);
   bool ifound = 0;

   std::vector<mmdb::Residue *> residues;

   int imod = 1;
   if (atom_sel.read_success > 0) { 
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id = chain_p->GetChainID();
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::PResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    residues.push_back(residue_p);
	    if (residues.size() == 3)
	       break;
	 }
	 if (residues.size() == 3)
	    break;
      }

      if (residues.size() == 3) {
	 mmdb::PPAtom residue_atoms_1;
	 mmdb::PPAtom residue_atoms_2;
	 mmdb::PPAtom residue_atoms_3;
	 int n_residue_atoms_1;
	 int n_residue_atoms_2;
	 int n_residue_atoms_3;
	 
	 std::string new_seg_id  = "N";
	 
	 residues[0]->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
	 for (int iat=0; iat<n_residue_atoms_1; iat++) {
	    mmdb::Atom *at = residue_atoms_1[iat];
	    at->SetAtomName(at->GetIndex(),
			    at->serNum,
			    at->GetAtomName(),
			    at->altLoc,
			    new_seg_id.c_str(),
			    at->GetElementName());
	 }

	 // testing this function:
	 coot::copy_segid(residues[0], residues[1]);
	 
	 residues[1]->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
	 for (int iat=0; iat<n_residue_atoms_2; iat++) {
	    mmdb::Atom *at = residue_atoms_2[iat];
	    std::string this_seg_id = at->segID;
	    if (this_seg_id != new_seg_id) {
	       std::cout << "   Failed to copy seg id.  Was :"
			 << this_seg_id << ": should be :" << new_seg_id
			 << ":\n for atom " << at << std::endl;
	       return 0; // fail
	    }
	 }


	 // Now something harder, test that the segids are unchanged
	 // when the segids of provider are not consistent.
	 //

	 std::cout << "   Test with a rogue segid " << std::endl;
	 // Add a rogue:
	 // 
	 mmdb::Atom *at = residue_atoms_1[2];
	 at->SetAtomName(at->GetIndex(),
			 at->serNum,
			 at->GetAtomName(),
			 at->altLoc,
			 "C",
			 at->GetElementName());

	 residues[2]->GetAtomTable(residue_atoms_3, n_residue_atoms_3);

	 std::vector<std::string> orig_seg_ids;
	 for (int iat=0; iat<n_residue_atoms_2; iat++) {
	    mmdb::Atom *at = residue_atoms_2[iat];
	    std::string this_seg_id = at->segID;
	    orig_seg_ids.push_back(this_seg_id);
	 }

	 coot::copy_segid(residues[0], residues[2]);
	 bool fail_this = 0;
	 for (int iat=0; iat<n_residue_atoms_2; iat++) {
	    mmdb::Atom *at = residue_atoms_2[iat];
	    std::string this_seg_id = at->segID;
	    if (this_seg_id != orig_seg_ids[iat]) {
	       std::cout << "  Failed: segid changed when it shouldn't"
			 << " have, for " << at << std::endl;
	       fail_this = 1;
	       break;
	    }
	 }

	 // if we go to here it should be ok if fail_this didn't fail
	 if (! fail_this)
	    status = 1;
      } 
   }
   return status;
}

int test_ligand_fit_from_given_point() {

   int status = 0;
   int n_conformers = 5;
   testing_data t;

   std::string cif_file_name = "libcheck_3GP-torsion-filtered.cif";
   coot::read_refmac_mon_lib_info_t rmit = t.geom.init_refmac_mon_lib(greg_test(cif_file_name), 0);
   if (rmit.n_bonds == 0) {
      std::string m = "Critical cif dictionary reading failure.";
      std::cout << m << std::endl;
      throw std::runtime_error(m);
   }
   
   std::string f = greg_test("tutorial-modern.pdb");
   atom_selection_container_t asc = get_atom_selection(f, true, true);
   if (!asc.read_success)
      return 0;
   std::string l = greg_test("monomer-3GP.pdb");
   atom_selection_container_t l_asc = get_atom_selection(l, true, true);
   if (!l_asc.read_success)
      return 0;
   
   
   clipper::Xmap<float> xmap;
   std::string mtz_file_name;
   mtz_file_name = getenv("HOME");
   mtz_file_name += "/data/greg-data/rnasa-1.8-all_refmac1.mtz";
   
   bool stat = coot::util::map_fill_from_mtz(&xmap, mtz_file_name,
					     "FWT", "PHWT", "WT", 0, 0);

   if (!stat) {
      std::cout << "   ERROR:: Bad map fill from " << mtz_file_name << "\n";
      return 0; 
   }

   
   coot::wligand wlig;
   wlig.set_verbose_reporting();
   wlig.set_debug_wiggly_ligands();
   wlig.import_map_from(xmap);
   bool optim_geom = 1;
   bool fill_vec = 0;
   coot::minimol::molecule mmol(l_asc.mol);
   int imol = 0; // dummy
   wlig.install_simple_wiggly_ligands(&t.geom, mmol, imol, n_conformers, optim_geom, fill_vec);
   short int mask_waters_flag = 1;
   wlig.mask_map(asc.mol, mask_waters_flag);
   clipper::Coord_orth pt(55.06, 10.16, 21.73); // close to 3GP peak (not in it).
   float n_sigma = 1.0; // cluster points must be more than this.
   wlig.cluster_from_point(pt, n_sigma);
   wlig.fit_ligands_to_clusters(1);
   int n_final_ligands = wlig.n_clusters_final();
   if (n_final_ligands == 0) {
      return 0;
   } else {
      unsigned int iclust = 0;
      unsigned int isol   = 0;
      coot::minimol::molecule m = wlig.get_solution(isol, iclust);
      clipper::Coord_orth centre = m.centre();
      clipper::Coord_orth ref_pt(55.5, 9.36, 20.7); // coords of centred
                                                    // correct solution
      double d = clipper::Coord_orth::length(centre, ref_pt);
      if (d < 1.0) {
	 std::cout << " found distance from reference centre "
		   << d << std::endl;
	 status = 1;
      }
   }

   return status;
}


int test_ligand_conformer_torsion_angles() { 

   int status = 0;
   testing_data t;

   std::string cif_file_name = "libcheck_3GP-torsion-filtered.cif";
   coot::read_refmac_mon_lib_info_t rmit = t.geom.init_refmac_mon_lib(greg_test(cif_file_name), 0);
   if (rmit.n_bonds == 0) {
      std::string m = "Critical cif dictionary reading failure.";
      std::cout << m << std::endl;
      throw std::runtime_error(m);
   }
   
   std::string l = greg_test("monomer-3GP.pdb");
   atom_selection_container_t l_asc = get_atom_selection(l, true, true);
   if (!l_asc.read_success)
      return 0;
   
   coot::wligand wlig;
   wlig.set_verbose_reporting();
   wlig.set_debug_wiggly_ligands();
   bool optim_geom = 0;
   bool fill_vec = 1;
   int n_conformers = 200;
   coot::minimol::molecule mmol(l_asc.mol);
   int imol = 0; // dummy
   std::vector<coot::installed_wiggly_ligand_info_t> conformer_info = 
      wlig.install_simple_wiggly_ligands(&t.geom, mmol, imol, n_conformers,
					 optim_geom, fill_vec);
   std::cout << "INFO:: there were " << conformer_info.size()
	     << " returned conformers" << std::endl;
   for (unsigned int i=0; i<conformer_info.size(); i++) {
      unsigned int itor=0;
      std::pair<float, float> p = conformer_info[i].get_set_and_real_torsions(itor);
      std::cout << "   " << i << " " << itor << "  set: " << p.first << " real: "
		<< p.second << std::endl;
   }

   return 1;
}

#include "coot-utils/peak-search.hh"
int test_peaksearch_non_close_peaks() {

   clipper::Xmap<float> xmap;
   std::string mtz_file_name;
   mtz_file_name = getenv("HOME");
   mtz_file_name += "/data/greg-data/rnasa-1.8-all_refmac1.mtz";
   
   bool stat = coot::util::map_fill_from_mtz(&xmap, mtz_file_name,
					     "FWT", "PHWT", "WT", 0, 0);
   if (!stat) {
      std::cout << "   ERROR:: Bad map fill from " << mtz_file_name << "\n";
      return 0; 
   }

   double d_crit = 2.0; // don't return peaks that have a higher
		       // (absolute) peak less than 4.0 A away.


   coot::peak_search ps(xmap);
   ps.set_max_closeness(d_crit);
   std::vector<std::pair<clipper::Coord_orth, float> > peaks = 
      ps.get_peaks(xmap, 0.5, 1, 1);

   if (peaks.size() < 20) {
      std::cout << "   Not enough peaks! " << peaks.size() << std::endl;
      return 0;
   } 

   std::vector<std::pair<clipper::Coord_orth, float> > problems;


   for (unsigned int ipeak=0; ipeak<(peaks.size()-1); ipeak++) {
      for (unsigned int jpeak=(ipeak+1); jpeak<peaks.size(); jpeak++) {
	 double d = clipper::Coord_orth::length(peaks[ipeak].first, peaks[jpeak].first);
	 if (d < d_crit) { 
	    problems.push_back(peaks[jpeak]);
	    break;
	 }
      }
   }

   std::cout << "   There are " << peaks.size() << " peaks and "
	     << problems.size() << " problem peaks" << std::endl;
   if (problems.size() == 0) {
      return 1;
   } else {
      return 0;
   } 
   
} 

int test_symop_card() {

   int r = 0;

   std::string symop_string = "X-1,Y,Z";
   coot::symm_card_composition_t sc(symop_string);
   std::cout << sc;

   float rev[9] = {1,0,0,0,1,0,0,0,1};
   float tev[3] = {-1,0,0};

   if ((close_float_p(sc.x_element[0], rev[0])) && 
       (close_float_p(sc.x_element[1], rev[1])) && 
       (close_float_p(sc.x_element[2], rev[2]))) {
      if ((close_float_p(sc.y_element[0], rev[3])) && 
	  (close_float_p(sc.y_element[1], rev[4])) && 
	  (close_float_p(sc.y_element[2], rev[5]))) {
	 if ((close_float_p(sc.z_element[0], rev[6])) && 
	     (close_float_p(sc.z_element[1], rev[7])) && 
	     (close_float_p(sc.z_element[2], rev[8]))) {
	    if ((close_float_p(sc.trans_frac(0), tev[0])) && 
		(close_float_p(sc.trans_frac(1), tev[1])) && 
		(close_float_p(sc.trans_frac(2), tev[2]))) {
	       r = 1;
	    }
	 }
      }
   }

   return r;
}

int test_coot_atom_tree() {

   std::cout << "Atom tree test" << std::endl;
   int r = 0;
   coot::dictionary_residue_restraints_t rest("XXX",0);
   mmdb::Residue *res = 0;

   // test that the exception is thrown by setting b = 1 there, don't continue if not
   bool b = 0;
   try { 
      coot::atom_tree_t tree(rest, res, "");
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
      b = 1;
   } 
   if (b == 0) { 
      std::cout << "No throw on null res" << std::endl;
      return 0;
   } 

   // Now test that the exception is thrown when there is no atom tree
   // (this test might not be necessary later).
   //
   b = 0;
   res = new mmdb::Residue;
   try { 
      coot::atom_tree_t tree(rest, res, "");
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
      b = 1;
   } 
   if (b == 0) { 
      std::cout << "No throw on no tree" << std::endl;
      return 0;
   } 

   delete res; // clear up first
   
   // OK give it something correct then

   // setup the restraints (p.second). Reading from a file in this
   // directory.  That should be fixed at some stage.
   
   std::string cif_file_name = "libcheck_ASP.cif";
   coot::protein_geometry geom;
   coot::read_refmac_mon_lib_info_t rmit = geom.init_refmac_mon_lib(greg_test(cif_file_name), 0);
   std::pair<short int, coot::dictionary_residue_restraints_t> p = 
      geom.get_monomer_restraints("ASP", 0);

   if (!p.first) {
      std::cout << "No ASP in dictionary" << std::endl;
      return 0;
   } 

   // Now get a residue
   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);
   bool ifound = 0;

   int imod = 1;
   res = test_get_residue(atom_sel.mol, "B", 1);

   if (0) {
      try {
	 r = test_tree_rotation(p.second, res, " CB ", " CG ", 0);
	 if (r)
	    r = test_tree_rotation(p.second, res, " CB ", " CG ", 1);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << rte.what() << std::endl;
      }
   }

   if (1) {
      try {
	 filename = "monomer-3GP.pdb";
	 atom_selection_container_t atom_sel = get_atom_selection(greg_test(filename), true, true);
	 if (!atom_sel.read_success) {
	    std::cout << "monomer-3GP.pdb not read successfully." << std::endl;
	 } else { 
	    mmdb::Residue *res = test_get_residue(atom_sel.mol, "A", 1);
	    if (res) {
	       coot::read_refmac_mon_lib_info_t rmit =
		  geom.init_refmac_mon_lib(greg_test("libcheck_3GP.cif"), 0);
	       std::pair<short int, coot::dictionary_residue_restraints_t> p = 
		  geom.get_monomer_restraints("3GP", 0);
	       if (p.first) { 
		  bool test1 = test_tree_rotation(p.second, res, " N9 ", " C1*", 0);
		  atom_sel.mol->WritePDBASCII("3GP-test1.pdb");
		  bool test2 = test_tree_rotation(p.second, res, " N9 ", " C1*", 1);
		  atom_sel.mol->WritePDBASCII("3GP-test2.pdb");
		  if (test1 && test2)
		     r = 1;
	       } else { 
		  std::cout << "Getting restraints for 3GP failed" << std::endl;
	       }
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << rte.what() << std::endl;
	 r = 0;
      }
   }

   return r;
}



bool
test_atom_tree_t::test_atom_vec(const std::vector<std::vector<int> > &contact_indices) const {

   bool r = 0;

   for (unsigned int iav=0; iav<atom_vertex_vec.size(); iav++) {
      std::cout << " vertex number: " << iav << " back [";
      for (unsigned int ib=0; ib<atom_vertex_vec[iav].backward.size(); ib++) {
	 std::cout << atom_vertex_vec[iav].backward[ib] << ","; 
      }
      std::cout << "] "; // end of backatom list

      std::cout << "forward ["; 
      for (unsigned int ifo=0; ifo<atom_vertex_vec[iav].forward.size(); ifo++) {
	 std::cout << atom_vertex_vec[iav].forward[ifo] << ","; 
      }
      std::cout << "] " << std::endl; // end of forward atom list
   }
   return r;
}

int
test_coot_atom_tree_2() {

   //                    o 4 
   //                   /
   //     0          1 /
   //     o --------- o                 0: NE1 (at origin)
   //     |           |                 1: CD
   //     |           |                 2: CE2
   //     |           |                 3: CZ
   //     o --------- o                 4: CG
   //     3           2
   // 
   int r = 0;

   std::vector<std::pair<std::string, clipper::Coord_orth> > names;
   names.push_back(std::pair<std::string, clipper::Coord_orth> (" NE1", clipper::Coord_orth(0,0,0)));
   names.push_back(std::pair<std::string, clipper::Coord_orth> (" CD ", clipper::Coord_orth(1,0,0)));
   names.push_back(std::pair<std::string, clipper::Coord_orth> (" CE2", clipper::Coord_orth(1,-1,0)));
   names.push_back(std::pair<std::string, clipper::Coord_orth> (" CZ ", clipper::Coord_orth(0,-1,0)));
   names.push_back(std::pair<std::string, clipper::Coord_orth> (" CG ", clipper::Coord_orth(0.5,1.5,0)));
		   
   mmdb::Residue *residue_p = new mmdb::Residue;
   for (int i=0; i<5; i++) {
      mmdb::Atom *at = new mmdb::Atom();
      at->SetAtomName(names[i].first.c_str());
      at->SetCoordinates(names[i].second.x(), names[i].second.y(), names[i].second.z(), 1.0, 20.0);
      residue_p->AddAtom(at);
   } 
   
   std::vector<std::vector<int> > contact_indices(5);
   contact_indices[0].push_back(1);
   contact_indices[0].push_back(3);
   contact_indices[1].push_back(0);
   contact_indices[1].push_back(2);
   contact_indices[1].push_back(4);
   contact_indices[2].push_back(1);
   contact_indices[2].push_back(3);
   contact_indices[3].push_back(0);
   contact_indices[3].push_back(2);
   contact_indices[4].push_back(1);

   test_atom_tree_t tat(contact_indices, 2, residue_p, "");
   bool tat_result = tat.test_atom_vec(contact_indices);

   double test_angle = 30.0; // degress
   double tar = (M_PI/180.0)* test_angle;
   bool reverse_flag = 0;
   tat.rotate_about(" CZ ", " CE2", tar, reverse_flag);
   reverse_flag = 1;
   tat.rotate_about(" CZ ", " CE2", tar, reverse_flag);
   
   delete residue_p;
   return r;
}

int
test_coot_atom_tree_proline() {

   int r = 0; 
   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename, true, true);
   mmdb::Residue *res_pro = test_get_residue(atom_sel.mol, "A", 12);
   if (res_pro) {
      coot::protein_geometry geom;
      geom.init_standard();
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      res_pro->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<clipper::Coord_orth> before_pos(n_residue_atoms);
      std::vector<clipper::Coord_orth>  after_pos(n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++)
	 before_pos[iat]=clipper::Coord_orth(residue_atoms[iat]->x,
					     residue_atoms[iat]->y,
					     residue_atoms[iat]->z);
      for (int iat=0; iat<n_residue_atoms; iat++) 
	 std::cout << "Atom Table " << iat << " " << residue_atoms[iat]->name << std::endl;
      std::vector<std::vector<int> > contact_indices =
	 coot::util::get_contact_indices_for_PRO_residue(residue_atoms,
							 n_residue_atoms,
							 &geom);

      int base_atom_index = 0;
      test_atom_tree_t tat(contact_indices, base_atom_index, res_pro, "");
      bool reverse_flag = 0;
      double test_angle = 30.0; // degress
      double tar = (M_PI/180.0)* test_angle;
      //      tat.rotate_about(" N  ", " CA ", tar, reverse_flag);
      tat.rotate_about(" CA ", " CB ", tar, reverse_flag);
      //      tat.rotate_about(" CB ", " CG ", tar, reverse_flag);

      for (int iat=0; iat<n_residue_atoms; iat++)
	 after_pos[iat]=clipper::Coord_orth(residue_atoms[iat]->x,
					    residue_atoms[iat]->y,
					    residue_atoms[iat]->z);
      for (int i=0; i<n_residue_atoms; i++) { 
	 double d = clipper::Coord_orth::length(before_pos[i],after_pos[i]);
	 if (d > 0.0001)
	    std::cout << "test: atom " << i << " " << residue_atoms[i]->name << " moved" << std::endl;
	 else 
	    std::cout << "test: atom " << i << " " << residue_atoms[i]->name << " static" << std::endl;
      }
   }
   return r;
}


// can return null;
mmdb::Residue *test_get_residue(mmdb::Manager *mol, const std::string &chain_id_ref, int resno_ref) {

   mmdb::Residue *residue_p = 0;
   
   int imod = 1;
   mmdb::Residue *res = 0;

   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      std::string chain_id = chain_p->GetChainID();
      if (chain_id == chain_id_ref) {
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) { 
	    res = chain_p->GetResidue(ires);
	    int resno = res->GetSeqNum();
	    if (resno == resno_ref) {
	       residue_p = res;
	       break;
	    }
	 }
      }
      if (residue_p)
	 break;
   }
   return residue_p;
} 


bool test_tree_rotation(const coot::dictionary_residue_restraints_t &rest,
			mmdb::Residue *res,
			const std::string &rotate_atom_1,
			const std::string &rotate_atom_2,
			bool reverse_flag) {


   
   bool r = 0;
   coot::atom_tree_t tree(rest, res, "");
   mmdb::PPAtom residue_atoms;
   int n_residue_atoms;
   res->GetAtomTable(residue_atoms, n_residue_atoms);
   std::vector<clipper::Coord_orth> before_pos(n_residue_atoms);
   std::vector<clipper::Coord_orth>  after_pos(n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++)
      before_pos[iat]=clipper::Coord_orth(residue_atoms[iat]->x,
					  residue_atoms[iat]->y,
					  residue_atoms[iat]->z);
   if (0) 
      for (int i=0; i<n_residue_atoms; i++)
	 std::cout << "   Test Before atom " << residue_atoms[i] << std::endl;

   double test_angle = 3.0; // degress
   double tar = (M_PI/180.0)* test_angle;
   
   tree.rotate_about(rotate_atom_1, rotate_atom_2, tar, reverse_flag);
   std::cout << std::endl; // separator
   for (int iat=0; iat<n_residue_atoms; iat++)
      after_pos[iat]=clipper::Coord_orth(residue_atoms[iat]->x,
					 residue_atoms[iat]->y,
					 residue_atoms[iat]->z);
   if (0) 
      for (int i=0; i<n_residue_atoms; i++)
	 std::cout << "    Test After atom " << residue_atoms[i] << std::endl;
      
   double dist = 0.0;
   for (int i=0; i<n_residue_atoms; i++) { 
      double d = clipper::Coord_orth::length(before_pos[i],after_pos[i]);
      if (d > 0.0001)
	 std::cout << "test: atom " << i << " " << residue_atoms[i]->name << " moved" << std::endl;
      else 
	 std::cout << "test: atom " << i << " " << residue_atoms[i]->name << " static" << std::endl;
      dist += d;
   }

   // Test that the moving atoms rotated the correct amount.
   //
   // First we need to find the rotate positions for the give atoms names.
   //

   clipper::Coord_orth r_pt_1;
   clipper::Coord_orth r_pt_2;
   for (int iat=0; iat<n_residue_atoms; iat++) {
      std::string at_name(residue_atoms[iat]->name);
      if (at_name == rotate_atom_1)
	 r_pt_1 = clipper::Coord_orth(residue_atoms[iat]->x,
				      residue_atoms[iat]->y,
				      residue_atoms[iat]->z);
      if (at_name == rotate_atom_2)
	 r_pt_2 = clipper::Coord_orth(residue_atoms[iat]->x,
				      residue_atoms[iat]->y,
				      residue_atoms[iat]->z);
   }

   r = 1; // intially success
   for (int i=0; i<n_residue_atoms; i++) { 
      double d = clipper::Coord_orth::length(before_pos[i], after_pos[i]);
      if (d > 0.0001) {
	 std::string at_name(residue_atoms[i]->name);
	 bool v = test_rotate_atom_angle(at_name,
					 r_pt_1, r_pt_2, before_pos[i], after_pos[i], test_angle);
	 if (v == 0) {
	    std::cout << " fail in test_rotate_atom_angle " << i << " "
		      << residue_atoms[i]->name << std::endl;
	    r = 0;
	    break;
	 }
      }
   }

   return r;

}


int 
test_rotate_around_vector() {

   int r = 0;
   
   std::string filename = "monomer-3GP.pdb";
   atom_selection_container_t atom_sel = get_atom_selection(greg_test(filename), true, true);

   std::string rotate_atom_1 = " N9 ";
   std::string rotate_atom_2 = " C1*";

   mmdb::Residue *residue_p = test_get_residue(atom_sel.mol, "A", 1);

   if (! residue_p) {
      std::cout << "residue not found for test_rotate_round_vector()" << std::endl;
   } else {
      int found_n_rotate_pts = 0;
      std::vector<int> exclude_atoms;
      clipper::Coord_orth rotate_pt_1;
      clipper::Coord_orth rotate_pt_2;
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 std::string at_name(residue_atoms[iat]->name);
	 if (at_name == rotate_atom_1) {
	    found_n_rotate_pts++;
	    rotate_pt_1 = clipper::Coord_orth(residue_atoms[iat]->x,
					      residue_atoms[iat]->y,
					      residue_atoms[iat]->z);
	    exclude_atoms.push_back(iat);
	 } 
	 if (at_name == rotate_atom_2) {
	    found_n_rotate_pts++;
	    rotate_pt_2 = clipper::Coord_orth(residue_atoms[iat]->x,
					      residue_atoms[iat]->y,
					      residue_atoms[iat]->z);
	    exclude_atoms.push_back(iat);
	 } 
      }

      if (found_n_rotate_pts != 2) {
	 std::cout << "rotate atoms not found for test_rotate_round_vector()" << std::endl;
      } else {
	 bool correct_rotation = 1; // success initially.
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    bool exclude = 0;
	    for (unsigned int iex=0; iex<exclude_atoms.size(); iex++) {
	       if (iat == exclude_atoms[iex]) {
		  exclude = 1;
		  break;
	       }
	    }

	    // We want point D.  Then we can measure CDC' for all
	    // rotated atoms. (C' is where C moves to after rotation).
	    //
	    //  \A
	    //   \\                                  ~
	    //    \  \                               ~
	    //     \   \                             ~
	    //      \    \                           ~
	    //       \     \                         ~
	    //        \      \                       ~
	    //         \ b     \                     ~
	    //        B \-------- C
	    //           .      /
	    //            .   /
	    //              /
	    //            D .         Angle BDC is 90 degrees.
	    //               .        BC is the hypotenuse of the BCD triangle.
	    //                .

	    if (0) {
	       rotate_pt_1 = clipper::Coord_orth(-1.0, 0.0, 0.0);
	       rotate_pt_2 = clipper::Coord_orth( 0.0, 0.0, 0.0);
	    }
	    
	    if (! exclude) {
	       clipper::Coord_orth C_pt(residue_atoms[iat]->x,
					residue_atoms[iat]->y,
					residue_atoms[iat]->z);

	       clipper::Coord_orth ab(rotate_pt_2-rotate_pt_1);
	       clipper::Coord_orth ba(rotate_pt_1-rotate_pt_2);
	       clipper::Coord_orth bc(C_pt-rotate_pt_2);
	       clipper::Coord_orth ab_unit(ab.unit());
	       double bclen = clipper::Coord_orth::length(rotate_pt_2, C_pt);
	       double ablen = clipper::Coord_orth::length(rotate_pt_1, rotate_pt_2);
	       double cos_b = clipper::Coord_orth::dot(ba, bc)/(ablen*bclen);
	       double bdlen = bclen * cos(M_PI-acos(cos_b));
	       clipper::Coord_orth D_pt = rotate_pt_2 + bdlen * ab_unit;

	       // Print angle at BDC:
	       clipper::Coord_orth cd(D_pt - C_pt);
	       clipper::Coord_orth bd(D_pt - rotate_pt_2);
	       double bdlen2 = clipper::Coord_orth::length(D_pt, rotate_pt_2);
	       double cdlen = clipper::Coord_orth::length(D_pt, C_pt);

	       if (0) { 
		  std::cout << "   D: " << D_pt.format() << " BD "<< bd.format() << " bdlen: "
			    << bdlen << " " << bdlen2 << " " << bc.format() << " " << cdlen
			    << " " << clipper::Coord_orth::dot(cd, bd) << " -> "
			    << (180.0/M_PI)*acos(clipper::Coord_orth::dot(cd, bd)/(bdlen*cdlen))
			    << " degrees " << std::endl;
		  std::cout << "   ab " << ab.format() << std::endl;
		  std::cout << "   ab_unit " << ab_unit.format() << std::endl;
		  std::cout << "   bclen " << bclen << std::endl;
		  std::cout << "   ablen " << ablen << std::endl;
		  std::cout << "   dot prod " << clipper::Coord_orth::dot(ab, bc) << std::endl;
		  std::cout << "   cos_b " << cos_b << std::endl;
		  std::cout << "   b     " << acos(cos_b) << std::endl;
		  std::cout << "   bdlen " << bdlen << std::endl;
		  std::cout << "   D_pt " << D_pt.format() << std::endl;
		  // Make sure that D_pt does not move when rotated:
		  for (double a=0; a<7.0; a+=1.0) {
		     clipper::Coord_orth D_pt_r =  coot::util::rotate_around_vector(ab, D_pt, rotate_pt_2, a);
		     std::cout << "   " << a << " " << D_pt_r.format() << std::endl;
		  }
	       }

	       double test_angle = 20.0; // degrees

	       clipper::Coord_orth C_prime_pt =
		  coot::util::rotate_around_vector(ab, C_pt, rotate_pt_2, (M_PI/180.0)*test_angle);

	       clipper::Coord_orth dc(C_pt-D_pt);
	       clipper::Coord_orth dc_prime(C_prime_pt-D_pt);
	       double dc_len = clipper::Coord_orth::length(C_pt,D_pt);
	       double dc_prime_len = clipper::Coord_orth::length(C_prime_pt,D_pt);

	       if (0) { 
		  std::cout << "   dc       " << dc.format() << std::endl;
		  std::cout << "   dc_prime " << dc_prime.format() << std::endl;
		  std::cout << "   dc_len   " << dc_len << std::endl;
		  std::cout << "   dc_prime_len " << dc_prime_len << std::endl;
		  std::cout << "   dot_prod " << clipper::Coord_orth::dot(dc, dc_prime) << std::endl;
	       }
	       

	       double cos_theta = clipper::Coord_orth::dot(dc, dc_prime)/(dc_len * dc_prime_len);

	       if (0) { 

		  std::cout << "   AB.DC " << clipper::Coord_orth::dot(ab,dc);
		  std::cout << "   AB.DC'" << clipper::Coord_orth::dot(ab,dc_prime) << std::endl;
		  
		  std::cout << "   AB.DC  angle " << (180.0/M_PI)*clipper::Coord_orth::dot(ab,dc)/(ablen * dc_len);
		  std::cout << "   AB.DC' angle " << (180.0/M_PI)*clipper::Coord_orth::dot(ab,dc_prime)/(ablen * dc_prime_len);
		  std::cout << std::endl;
	       }
	       
	       std::cout << "   " << iat << " " << residue_atoms[iat]->name << " "
			 << cos_theta << " -> " << acos(cos_theta)*180.0/M_PI << " degrees" << std::endl;

	       double real_rotation_angle = acos(cos_theta)*180.0/M_PI;
	       if (! close_float_p(test_angle, real_rotation_angle))
		  correct_rotation = 0;
	       
	       residue_atoms[iat]->x = C_prime_pt.x();
	       residue_atoms[iat]->y = C_prime_pt.y();
	       residue_atoms[iat]->z = C_prime_pt.z();
	    }
	 }
	 atom_sel.mol->WritePDBASCII("3gp-rotated.pdb");
	 r = correct_rotation;
      }
   }

   return r;
}

// Is the angle that after_pos is rotated round the r_pt_1 and r_pt_2
// vector the same as test_angle?
// 
bool
test_rotate_atom_angle(const std::string &atom_name,
		       const clipper::Coord_orth &rotate_pt_1,
		       const clipper::Coord_orth &rotate_pt_2,
		       const clipper::Coord_orth &C_pt,
		       const clipper::Coord_orth &C_prime_pt,
		       double test_angle) {
   bool r = 0;

   clipper::Coord_orth ab(rotate_pt_2-rotate_pt_1);
   clipper::Coord_orth bc(C_pt-rotate_pt_2);
   clipper::Coord_orth ab_unit(ab.unit());
   double bclen = clipper::Coord_orth::length(rotate_pt_2, C_pt);
   double ablen = clipper::Coord_orth::length(rotate_pt_1, rotate_pt_2);
   double cos_b = clipper::Coord_orth::dot(ab, bc)/(ablen*bclen);
   double bdlen = bclen * cos(M_PI-acos(cos_b));
   clipper::Coord_orth D_pt = rotate_pt_2 - bdlen * ab_unit;

   // Print angle at BDC:
   clipper::Coord_orth cd(D_pt - C_pt);
   clipper::Coord_orth bd(D_pt - rotate_pt_2);
   double bdlen2 = clipper::Coord_orth::length(D_pt, rotate_pt_2);
   double cdlen = clipper::Coord_orth::length(D_pt, C_pt);
   clipper::Coord_orth dc(C_pt-D_pt);
   clipper::Coord_orth dc_prime(C_prime_pt-D_pt);
   double dc_len = clipper::Coord_orth::length(C_pt,D_pt);
   double dc_prime_len = clipper::Coord_orth::length(C_prime_pt,D_pt);
   double cos_theta = clipper::Coord_orth::dot(dc, dc_prime)/(dc_len * dc_prime_len);
   double real_angle = acos(cos_theta)*180.0/M_PI;
   std::cout << "  " << atom_name << " " << cos_theta << " -> "
	     << real_angle << " degrees" << std::endl;
   if (close_float_p(real_angle, test_angle)) { 
      r = 1;
   } else {
      std::cout << "   Ooops " << real_angle << " not close to " << test_angle << std::endl;
   } 
   return r;
} 

int test_ssm_sequence_formatting() {

#ifdef HAVE_SSMLIB   
   graphics_info_t g;
   std::pair<std::string, std::string> aligned_sequences;

   std::string s;
   std::string t;
   aligned_sequences.first = s;
   aligned_sequences.second = t;
   g.print_horizontal_ssm_sequence_alignment(aligned_sequences);
   std::cout << "--" << std::endl;

   s = "DVSGTVCLSALPPEATDTLNLIASDGPFPYSQDGVVFQNR--ESVLPTQSYG";
   s += "YYHEYTVITP--GARTRGTRRI.ICGEATQEDY..YTGDHYATFSLIDQTC";
   t = "---SGTVCLSALPPEATDTLNLIASDGPFPYSQDG";
   aligned_sequences.first = s;
   aligned_sequences.second = t;
   g.print_horizontal_ssm_sequence_alignment(aligned_sequences);
   std::cout << "--" << std::endl;

   s = "D";
   t = "--SGTVCLSALPPEATDTLNLIASDGPFPYSQDG";
   aligned_sequences.first = s;
   aligned_sequences.second = t;
   g.print_horizontal_ssm_sequence_alignment(aligned_sequences);
   std::cout << "--" << std::endl;
   
   s = "DVSGTVCLSALPPEATDTLNIASDGPFPYSQDGVVFQNR--ESVLPQSYG";
   t = "--SGTVCLSALPPEATDTLNIASDGPFPYSQDXXxxxxxxxxxxxxxxxG";
   aligned_sequences.first = s;
   aligned_sequences.second = t;
   g.print_horizontal_ssm_sequence_alignment(aligned_sequences);
   std::cout << "--" << std::endl;
#endif // HAVE_SSMLIB
   
   return 1;
} 

int test_OXT_in_restraints() {

   int r = 0; // initially fail.
   coot::protein_geometry geom;
   geom.init_standard();
   std::string cif_file_name = greg_test("libcheck_BCS.cif");
   coot::read_refmac_mon_lib_info_t rmit = geom.init_refmac_mon_lib(cif_file_name, 0);
   if (! rmit.success) {
      std::cout << "Fail to get good status from reading " << cif_file_name << std::endl;
   } else { 
      bool v1 = geom.OXT_in_residue_restraints_p("TRP");
      bool v2 = geom.OXT_in_residue_restraints_p("BCS");

// This doesn't work for new restraints (which *do* have OXTs)
//       if (v1 == 0)
// 	 if (v2 == 1)
// 	    r = 1;
// 	 else
// 	    std::cout << "fail to find OXT in BSC" << std::endl;
//       else
// 	 std::cout << "Fail to not find OXT in TRP" << std::endl;

      if (v2 == 1)
	 r = 1;
      else
	 std::cout << "fail to find OXT in BSC" << std::endl;

   }
   return r; 
} 

int test_relativise_file_name () {

   std::string f1 = "/a/b";
   std::string f2 = "/c/a/b";
   std::string f3 = "/c/b";
   std::string f4 = "/a";
   std::string  c = "/a";

   std::string r1 = coot::util::relativise_file_name(f1, c);
   if (r1 != "b") {
      std::cout << "FAIL: relativise_file_name(" << f1 << ", " << c << ") gives " << r1 << "\n";
      return 0;
   }
   std::string r2 = coot::util::relativise_file_name(f2, c);
   if (r2 != f2) {
      std::cout << "FAIL: relativise_file_name(" << f2 << ", " << c << ") gives " << r2 << "\n";
      return 0;
   }
   std::string r3 = coot::util::relativise_file_name(f3, c);
   if (r3 != f3) {
      std::cout << "FAIL: relativise_file_name(" << f3 << ", " << c << ") gives " << r3 << "\n";
      return 0;
   }
   std::string r4 = coot::util::relativise_file_name(f4, c);
   if (r4 != f4) {
      std::cout << "FAIL: relativise_file_name(" << f4 << ", " << c << ") gives " << r4 << "\n";
      return 0;
   }

   return 1;
}

int test_previous_water() {

   coot::protein_geometry geom;
   geom.init_standard();
   int status = 0;
   molecule_class_info_t mci;
   mci.handle_read_draw_molecule(1,
				 greg_test("pathological-water-test.pdb"),
				 coot::util::current_working_dir(),
				 &geom,
				 0, 0, true, true, 1, coot::NORMAL_BONDS, false);
   mci.delete_atom("D", 162, "", " O  ", "");
   coot::Cartesian rc(0,0,0); // hack?
   int iprev = mci.intelligent_previous_atom("D", 162, " O  ", "", rc);
   mmdb::Atom *at = mci.atom_sel.atom_selection[iprev];
   std::cout << "previous atom: " << at << std::endl;
   if (std::string(at->GetChainID()) == "D")
      if (at->GetSeqNum() == 161)
	 status = 1;
   
   std::cout << "returning " << status << std::endl;
   return status;

}

int test_ccp4srs() {

#ifdef HAVE_CCP4SRS   
   int status = 0;
   testing_data t;
   const char *sbase_dir = getenv("COOT_SBASE_DIR");
   if (sbase_dir) {
      std::string sbase_monomer_dir = sbase_dir;
      int init_status = t.geom.init_ccp4srs(sbase_monomer_dir);
      std::cout << "   SBase init status: " << init_status << std::endl;
      if (init_status != ccp4srs::CCP4SRS_Ok) {
	 std::cout << "   WARNING:: Trouble initialising SBase" << std::endl;
      } else { 
	 // t.geom.read_sbase_residues();
	 status = 1;

	 std::string test_name= "AMP";
	 std::vector<std::pair<std::string,std::string> > v = t.geom.matching_ccp4srs_residues_names(test_name);
	 std::cout << "INFO:: " << v.size() << " matching residue names" << std::endl;
	 for (unsigned int i=0; i<v.size(); i++) { 
	    std::cout << "    " << i << " of " << v.size() << " "
		      << v[i].first << std::endl;
	 }
      }
   }
   return status;
#else
   return 1;
#endif // HAVE_CCP4SRS   
}


int test_coordinated_waters() {

   int status = 0;
   double water_limit = 2.9; // contacts need to be less than this to be written out.

   int data_type = COOT_COORDS_FILE_SELECTION;
   std::vector<std::string> file_names = filtered_by_glob("coot-download", data_type);
   for (unsigned int i=0; i<file_names.size(); i++) {
      atom_selection_container_t atom_sel = get_atom_selection(file_names[i], true, true);
      if (atom_sel.mol) { 
	 coot::util::water_coordination_t wc(atom_sel.mol, 3.3);
	 std::vector<coot::util::contact_atoms_info_t> water_contacts = 
	    wc.get_highly_coordinated_waters(5, water_limit);
	 if (water_contacts.size() > 0) {
	    std::cout << "    " << water_contacts.size() << std::endl;
	    for (unsigned int j=0; j<water_contacts.size(); j++) {
	       std::cout << "       ";
	       std::cout << water_contacts[j].central_atom();
	       std::cout << "\n";
	       for (unsigned int k=0; k<water_contacts[j].size(); k++) {
		  coot::util::contact_atoms_info_t::contact_atom_t at = water_contacts[j][k];
		  if (at.dist < water_limit) 
		     std::cout << "              " << at.dist << "  " << at.at << std::endl;
	       }
	    }
	 }
      }
   }

   return status;
} 

int test_geometry_distortion_info_type() {

   int status = 0;
   coot::simple_restraint rest;
   coot::residue_spec_t rs("A", 1);
   
   coot::geometry_distortion_info_t gdi_1(6, rest, rs);
   coot::geometry_distortion_info_t gdi_2(7, rest, rs);
   coot::geometry_distortion_info_t gdi_3; // undefined, no distortion
   

   if (gdi_1 < gdi_2) {
      if (gdi_2 > gdi_1) {
	 bool cont = 1;
	 try {
	    if (gdi_3 < gdi_2) // just use the test - does nothing useful
	       int x = 2; 
	    cont = 0;
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "    Good gdi < exception thrown" << std::endl;
	 }
	 if (cont) { 
	    try {
	       if (gdi_3 < gdi_2)
		  int x = 2;
	       cont = 0;
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "    Good gdi > exception thrown" << std::endl;
	       status = 1;
	    }
	 }
      } else {
	 std::cout << "test geometry_distortion_info_t > fails" << std::endl;
      }
   } else { 
      std::cout << "test geometry_distortion_info_t < fails" << std::endl;
   }

   return status;
}


int test_translate_close_to_origin() {

   int status = 0;

   clipper::Coord_orth origin(0,0,0);
   std::vector<clipper::Coord_orth> pts;
   pts.push_back(clipper::Coord_orth(99.9, 100.1, 100.0));
   mmdb::Manager *mol = coot::util::create_mmdbmanager_from_points(pts);
   clipper::Cell_descr cell_descr(100,100,100,
				  clipper::Util::d2rad(90.0), 
				  clipper::Util::d2rad(90.0), 
				  clipper::Util::d2rad(90.0));
   clipper::Cell cell(cell_descr);
   bool cell_status = coot::util::set_mol_cell(mol, cell);
   if (! cell_status) {
      std::cout << "failure to set cell" << std::endl;
   } else { 
      mol->SetSpaceGroup("P 21 21 21");
      coot::util::translate_close_to_origin(mol);
      std::pair<bool, clipper::Coord_orth> c = coot::centre_of_molecule(mol);
      if (c.first) { 
	 double len = clipper::Coord_orth::length(c.second, origin);
	 std::cout << "    Got length " << len << std::endl;
	 if (len < 0.2)
	    status = 1;
      }
   }
   return status;
}

int test_flev_aromatics() {

   int status = 0;
   // std::string filename = "test-with-5GP.pdb";
   std::string filename = "test-with-5GP-with-ideal-A37-PHE.pdb";
   // std::string filename = "coot-download/1x8b.pdb";
   atom_selection_container_t atom_sel = get_atom_selection(greg_test(filename), true, true);
   mmdb::Residue *res_ref = coot::util::get_residue("C", 1, "", atom_sel.mol);
   // mmdb::Residue *res_ref = coot::util::get_residue("A", 901, "", atom_sel.mol);
   if (! res_ref) {
      std::cout << "failed to get reference residue in test_flev_aromatics()" << std::endl;
      return 0;
   } else {

      testing_data t;
      int dynamic_add_status = t.geom.try_dynamic_add("5GP", 1);
      // read_cif_dictionary("coot-ccp4/.coot-to-lbg-lib");
      std::cout << "DEBUG:: dynamic_add_status " << dynamic_add_status << std::endl;
      float residues_near_radius = 4.0;
      std::vector<mmdb::Residue *> residues =
	 coot::residues_near_residue(res_ref, atom_sel.mol, residues_near_radius);
      std::pair<bool, coot::dictionary_residue_restraints_t> p = 
	 t.geom.get_monomer_restraints("5GP", 0);
      coot::pi_stacking_container_t pi_stack_info(p.second, residues, res_ref);

      if (pi_stack_info.stackings.size() > 0)
	 status = 1;

   }
   return status;
}

int test_map_segmentation() {

   std::string filename = "emd_1661.map";
   clipper::CCP4MAPfile file;
   try { 
      file.open_read(filename);
      clipper::Xmap<float> xmap;
      file.import_xmap(xmap);
      float low_level = 0.0524; // 0.075; // 0.02; // 0.005;
      coot::util::segment_map s;
      std::pair<int, clipper::Xmap<int> > segmented_map = s.segment(xmap, low_level);

      clipper::CCP4MAPfile mapout;
      mapout.open_write(std::string("segmented.map"));
      mapout.export_xmap(segmented_map.second);
      mapout.close_write();
   }
   catch (const clipper::Message_base &exc) {
      std::cout <<  "WARNING:: failed to open " << filename << std::endl;
   }

   return 1;

}

int test_lsq_plane() {

   int r = 0;
   
   std::vector<clipper::Coord_orth> v;
   clipper::Coord_orth pt(0.5, 0.5, 0.1); // in the plane, hopefully.

   // The eigenvalues should not be equal, else nans.
   // 
   v.push_back(clipper::Coord_orth(0.0, 0.0, 0.0));
   v.push_back(clipper::Coord_orth(1.0, 0.0, 0.2));
   v.push_back(clipper::Coord_orth(1.0, 1.1, 0.2));
   v.push_back(clipper::Coord_orth(0.0, 1.0, 0.0));
   
   std::pair<double, double> d = coot::lsq_plane_deviation(v, pt);
   std::cout << "LSQ deviations: " << d.first << " " << d.second << std::endl;
   if (close_float_p(d.first, 0.0)) {
      r = 1;
   } 
   return r;
}


// This test fails.  I can't get PutPDBString() to work.
// BL says:: PutPDB works now (used the wrong CRYST format)
// However test still fails as Cryst.Vol is not calculated I think
// 
int test_copy_cell_symm_orig_scale_headers() {

   int r = 0;

   mmdb::Manager *m1 = new mmdb::Manager;
   mmdb::Manager *m2 = new mmdb::Manager;

   int set1 = m1->PutPDBString("CRYST1   69.782   69.782  157.017  90.00  90.00  90.00 P 41 21 2     8");
   int set2 = m1->PutPDBString("ORIGX1      1.000000  0.000000  0.000000        0.00000");                         
   int set3 = m1->PutPDBString("ORIGX2      0.000000  1.000000  0.000000        0.00000");                         
   int set4 = m1->PutPDBString("ORIGX3      0.000000  0.000000  1.000000        0.00000");                         
   int set5 = m1->PutPDBString("SCALE1      0.014330  0.000000  0.000000        0.00000");                         
   int set6 = m1->PutPDBString("SCALE2      0.000000  0.014330  0.000000        0.00000");                         
   int set7 = m1->PutPDBString("SCALE3      0.000000  0.000000  0.006369        0.00000");
   m1->PutPDBString("ATOM      1  N   MET A 291     -11.787  76.079  32.455  1.00 46.95           N  ");
   m1->PutPDBString("ATOM      2  CA  MET A 291     -10.759  74.985  32.450  1.00 46.65           C  ");
   m1->PutPDBString("ATOM      3  C   MET A 291      -9.337  75.415  32.821  1.00 45.29           C  ");
   m1->PutPDBString("ATOM      4  O   MET A 291      -9.056  75.720  33.979  1.00 45.72           O  ");


   std::cout << "sets: "
	     << set1 << " " << set2 << " " << set3 << " "
	     << set4 << " " << set5 << " " << set6 << " "
	     << set7 << " " << std::endl;

   // m1->ReadPDBASCII("test.pdb");

   const char *spc_1 = m1->GetSpaceGroup();
   if (spc_1 == 0)
      throw std::runtime_error("fail to set spacegroup with PutPDBString");

   std::cout << "m1 spacegroup " << spc_1 << std::endl;

   mmdb::realtype a[6];
   mmdb::realtype vol;
   int orthcode;
   m1->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
   std::cout << "PutPDBString: cell "
	     << a[0] << " " << a[1] << " " << a[2] << " "
	     << a[3] << " " << a[4] << " " << a[5] << " " << vol  << " " << orthcode
	     << std::endl;

   
   bool r1 = coot::util::copy_cell_and_symm_headers(m1, m2);

   const char *spc_2 = m2->GetSpaceGroup();
   if (!spc_2)
      throw std::runtime_error("fail to convert spacegroup (NULL)");

   std::cout << "debug spacegroup " << spc_2 << std::endl;
   std::string sp = spc_2;
   if (sp != "P 41 21 2")
      throw std::runtime_error("failed to set correct space group");

   m2->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
   std::cout << "Copied cell "
	     << a[0] << " " << a[1] << " " << a[2] << " "
	     << a[3] << " " << a[4] << " " << a[5] << " " << vol  << " " << orthcode
	     << std::endl;

   if (vol < 70*70*150)
      throw std::runtime_error("failed to set correct cell");

   
   delete m1;
   delete m2;

   if (! r1) {
      std::cout << "coot::util::copy_cell_and_symm_headers() fails" << std::endl;
      return 0;
   } 
   return 1;
}

int test_residue_atom_renaming() {

   // std::string f1 = "coot-ccp4/store-1/prodrg-DRG.pdb";
   // std::string f2 = "coot-ccp4/store-2/prodrg-DRG.pdb";

   std::string f1 = "coot-ccp4/store-1/prodrg-DRG-with-H.pdb";
   std::string f2 = "coot-ccp4/store-2/prodrg-DRG-with-H.pdb";

   atom_selection_container_t atom_sel_ref = get_atom_selection(f1, true, true);
   atom_selection_container_t atom_sel_mov = get_atom_selection(f2, true, true);

   std::vector<std::string> orig_atom_names;
   std::vector<std::string> curr_atom_names;

   mmdb::Residue *res_ref = atom_sel_ref.mol->GetResidue(1, "", 1, "");
   mmdb::Residue *res_mov = atom_sel_mov.mol->GetResidue(1, "", 1, "");

   if (! res_ref || ! res_mov) {
      std::cout << "test_residue_atom_renaming(): failed to get residues "
		<< res_ref << " " << res_mov << std::endl;
   } else {
      // normal path
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      res_mov->GetAtomTable(residue_atoms, n_residue_atoms);
      orig_atom_names.resize(n_residue_atoms);
      curr_atom_names.resize(n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++)
	 orig_atom_names[i] = residue_atoms[i]->name;

      bool match_hydrogens_also = 1;
      coot::graph_match_info_t gm = coot::graph_match(res_mov, res_ref, 0, match_hydrogens_also);
      gm.match_names(res_mov);

      for (int i=0; i<n_residue_atoms; i++)
	 curr_atom_names[i] = residue_atoms[i]->name;

      //  CAL ->  CAD
      //  HAE ->  HAM
      //  CAJ ->  CAE
      // CLAN -> CLAH
      //  CAI ->  CAA
      //  HAA ->  HAL
      
      for (int i=0; i<n_residue_atoms; i++)
	 std::cout << " test ---  was :" << orig_atom_names[i] << ":  now :"
		   << curr_atom_names[i] << ":" << std::endl;

      // Test that there are no atoms in the residue that occur more
      // than once.
      std::map<std::string, int> atom_name_map;
      for (int i=0; i<n_residue_atoms; i++) {
	 atom_name_map[residue_atoms[i]->name]++;
	 if (atom_name_map[residue_atoms[i]->name] != 1) {
	    std::string m = "   fail - more than 1 of :";
	    m += residue_atoms[i]->name;
	    m+= ":";
	    throw (std::runtime_error(m));
	 }
      }
      
      
   }
   return 0;
}

int test_mcd_and_thornton_h_bonds() {

   int r = 0;

   testing_data t;
   t.geom.init_refmac_mon_lib(greg_test("SGP-modified.cif"), 0);
   atom_selection_container_t asc = get_atom_selection(greg_test("test-hydrogenated-region.pdb"), true, false);
   if (asc.read_success) {

      int SelHnd_all = asc.mol->NewSelection();
      int SelHnd_lig = asc.mol->NewSelection();
      asc.mol->SelectAtoms(SelHnd_all, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      asc.mol->SelectAtoms(SelHnd_lig, 0, "A", 97, "", 97, "", "*", "*", "*", "*");

      coot::h_bonds hb;
      std::vector<coot::h_bond> hbonds =
	 hb.get_mcdonald_and_thornton(SelHnd_lig, SelHnd_all, asc.mol, t.geom);

      std::cout << "Returned H-bonds:" << std::endl;
      for (unsigned int i=0; i<hbonds.size(); i++) { 
	 std::cout << "   " << i << "  " << hbonds[i] << std::endl;
      }
   } 
   return r;
}

int test_COO_mod() {

   graphics_info_t g;
   testing_data t;
   int status = 0;
   std::string file_name = greg_test("hideous-OXT.pdb");
   atom_selection_container_t asc = get_atom_selection(file_name, true, true);

   
   if (!  asc.read_success) {
      std::cout << "failed to correctly read hideous-OXT.pdb " << std::endl;
   } else { 
      std::cout << "read " << asc.n_selected_atoms << " atom " << std::endl;

      // refine ...
      
      mmdb::PResidue *SelResidues = new mmdb::PResidue[1];
      SelResidues[0] = asc.atom_selection[0]->residue;

      residue_selection_t result =
	 testing_func_probabilities_refine_fragment(asc, SelResidues,
						    1, "A", 93, t.geom, 0, 0, 0, 0);
      delete [] SelResidues;

      // now test...
   
      std::vector<int> atom_index(3);
      atom_index[0] = 1; // CA
      atom_index[1] = 6; // C
      atom_index[2] = 7; // O
      for (unsigned int i=0; i<3; i++) {
	 clipper::Coord_orth pos(result.SelResidues[0]->GetAtom(atom_index[i])->x,
				 result.SelResidues[0]->GetAtom(atom_index[i])->y,
				 result.SelResidues[0]->GetAtom(atom_index[i])->z);

	 g.lsq_plane_atom_positions->push_back(pos);
      }
      clipper::Coord_orth oxt_pos(result.SelResidues[0]->GetAtom(8)->x,
				  result.SelResidues[0]->GetAtom(8)->y,
				  result.SelResidues[0]->GetAtom(8)->z);
      clipper::Coord_orth o_pos(result.SelResidues[0]->GetAtom(7)->x,
				result.SelResidues[0]->GetAtom(7)->y,
				result.SelResidues[0]->GetAtom(7)->z);
      result.clear_up();
   
      std::pair<float,float> d_pair =
	 coot::lsq_plane_deviation(*g.lsq_plane_atom_positions, oxt_pos);
      float d = fabs(d_pair.first);
      std::cout << "OXT out of plane distance: " << d << std::endl;
      double oxt_dist = clipper::Coord_orth::length(o_pos, oxt_pos);
      std::cout << "OXT->O distance: " << oxt_dist << std::endl;
   
      if (d < 0.02)
	 if (oxt_dist > 2.0)
	    status = 1;

   }
   return status;
}

int test_new_comp_id() {

   int status = 1;

   std::vector<std::pair<std::string, std::string> > comp_ids;
   comp_ids.push_back(std::pair<std::string, std::string> ("L19", "L20"));
   comp_ids.push_back(std::pair<std::string, std::string> ("LIG", "LI2"));
   comp_ids.push_back(std::pair<std::string, std::string> ("L01", "L02"));
   comp_ids.push_back(std::pair<std::string, std::string> ("119", "120"));
   comp_ids.push_back(std::pair<std::string, std::string> ("120", "121"));
   comp_ids.push_back(std::pair<std::string, std::string> ("D99", "")); // failure case

   for (unsigned int i=0; i<comp_ids.size(); i++) {
      std::string n = coot::suggest_new_comp_id(comp_ids[i].first);
      if (n != comp_ids[i].second) {
	 std::cout << "New comp_id fail on " << comp_ids[i].first << " wanted " << comp_ids[i].second
		   << " but got \"" << n << "\"" << std::endl;
	 status = 0; // fail
	 break;
      }
   }
   return status;
}

int test_trailing_slash() {

   int status = 1; // OK
   std::string s = "x/";
   if (coot::util::remove_trailing_slash(s) != "x") {
      status = 0;
   }
   s = "/";
   if (coot::util::remove_trailing_slash(s) != "") {
      status = 0;
   }
   s = "ss";
   if (coot::util::remove_trailing_slash(s) != "ss") {
      status = 0;
   }
   s = "\\"; // single
   if (coot::util::remove_trailing_slash(s) != "") {
      status = 0;
   }
   s = "";
   if (coot::util::remove_trailing_slash(s) != "") {
      status = 0;
   }
   return status;
}

 

int test_remove_whitespace() {

   std::string s = "";
   if (coot::util::remove_trailing_whitespace(s) != "") {
      std::cout << "fail on 1" << std::endl;
      return 0;
   } 
   s = "zz";
   if (coot::util::remove_trailing_whitespace(s) != "zz") {
      std::cout << "fail on 2" << std::endl;
      return 0;
   } 
   s = "  zz";
   if (coot::util::remove_trailing_whitespace(s) != "  zz") { 
      std::cout << "fail on 3" << std::endl;
      return 0;
   } 
   s = "  zz x";
   if (coot::util::remove_trailing_whitespace(s) != "  zz x") { 
      std::cout << "fail on 4" << std::endl;
      return 0;
   } 
   s = "  zz  xx   ";
   if (coot::util::remove_trailing_whitespace(s) != "  zz  xx") { 
      std::cout << "fail on 5" << std::endl;
      return 0;
   }
   return 1;
}

int test_position_residue_by_internal_coords() {

#ifdef HAVE_CCP4SRS   

   // This test doesn't do proper test (as yet) and the class it tests
   // is not complete either.

   int status = 0;

   testing_data t;
   int idealised_flag = 1;

   const char *sbase_dir = getenv("COOT_CCP4SRS_DIR");
   if (sbase_dir) {
      std::string sbase_monomer_dir = sbase_dir;
      int init_status = t.geom.init_ccp4srs(sbase_monomer_dir);
      if (init_status != ccp4srs::CCP4SRS_Ok) {
	 std::cout << "   WARNING:: Trouble initialising SBase" << std::endl;
      } else { 

	 mmdb::Manager *r_mol = new mmdb::Manager;
	 r_mol->ReadPDBASCII("coot-ccp4/monomer-ASN.pdb");
   
	 mmdb::Residue *r = coot::util::get_first_residue(r_mol);
	 mmdb::Residue *m = t.geom.get_ccp4srs_residue("NAG");

	 if (! r_mol) {
	    std::cout << "Failed to get ASN molecule " << std::endl;
	 } else {
	    if (!r) {
	       std::cout << "Failed to get ASN residue " << std::endl;
	    } else { 
	       if (!m) {
		  std::cout << "Failed to get NAG residue " << std::endl;
	       } else {

		  coot::atom_name_quad quad(" CG ", " ND2", " C1 ", " O5 ");
		  quad.set_atom_residue_index(2,2);
		  quad.set_atom_residue_index(3,2);
		  coot::position_residue_by_internal_coordinates p(r, m, quad, 1.45, 108.9, -100.4);
		  status = p.move_moving_residue();
		  std::cout << "DEBUG:: move_moving_residue status: " << status << std::endl;
	       } 
	    }
	 }
      }
   }
   return status;
#else
   return 1;
#endif // HAVE_CCP4SRS   
} 


int test_beam_in_residue() {

   testing_data t;
   int status = 0;
   mmdb::Manager *r_mol = new mmdb::Manager;
   r_mol->ReadPDBASCII("coot-ccp4/monomer-ASN.pdb");
   mmdb::Residue *r = coot::util::get_first_residue(r_mol);
   if (r) {
      coot::beam_in_linked_residue lr(r, "NAG-ASN", "NAG", &t.geom);
      mmdb::Residue *result = lr.get_residue();
      if (result) {
	 status = 1;
      } 
   } 

   return status;
}

int test_multi_residue_torsion() {

   int status = 0;
   testing_data t;
   t.geom.try_dynamic_add("NAG", 1);
   mmdb::Manager *mol = new mmdb::Manager;
   mol->ReadPDBASCII("ASN-NAG-pair.pdb");
   mmdb::Residue *res_1 = coot::util::get_first_residue(mol);
   if (res_1) {
      coot::residue_spec_t specs[2];
      int selhnd = mol->NewSelection();
      specs[0] = coot::residue_spec_t("A", 131, "");
      specs[1] = coot::residue_spec_t("A", 361, "");
      for (unsigned int i=0; i<2; i++) {
	 specs[i].select_atoms(mol, selhnd, mmdb::SKEY_OR); 
      }

      mmdb::PPAtom atom_selection;
      int n_selected_atoms;
      mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);

      // now we need to add the link, we need a mmdb::Residue * for the
      // second residue to do that.
      //
      int selhnd_res2 = mol->NewSelection();
      mol->Select(selhnd_res2, mmdb::STYPE_RESIDUE, 1, "A",
		  361, "", 361, "", "*", "*", "*", "*", mmdb::SKEY_NEW);
      int nSelResidues;
      mmdb::PPResidue SelResidues;
      mol->GetSelIndex(selhnd_res2, SelResidues, nSelResidues);
      
      if (nSelResidues != 1) {
	 std::cout << "problem in test_multi_residue_torsion" << std::endl;
      } else {
	 mmdb::Residue *res_2 = SelResidues[0];
	 coot::bonded_pair_t bp(res_1, res_2, 0, 0, "NAG-ASN");
	 coot::bonded_pair_container_t bpc;
	 bpc.try_add(bp);
	 atom_selection_container_t asc;
	 asc.mol = mol;
	 asc.atom_selection = atom_selection;
	 asc.n_selected_atoms = n_selected_atoms;
	 asc.SelectionHandle = selhnd;

	 for (int i=0; i<n_selected_atoms; i++) { 
	    std::cout << "selected atom: " << atom_selection[i] << std::endl;
	 }
	 int imol = 0; // dummy
	 coot::contact_info contacts(asc, imol, &t.geom, bpc);
 	 std::vector<std::vector<int> > contact_indices =
 	    contacts.get_contact_indices_with_reverse_contacts();
 	 coot::atom_tree_t tree(contact_indices, 0, mol, selhnd);
 	 tree.rotate_about(1, 4, 3, 0); // CA-CB
	 mol->WritePDBASCII("rotated.pdb");
      } 

      // clean up
      // 
      mol->DeleteSelection(selhnd_res2);
      mol->DeleteSelection(selhnd);
      
   }
   return status;
}

int
test_torsions_from_residue_selection() {

   int status = 1;
   graphics_info_t g;
   testing_data t;
   t.geom.try_dynamic_add("NAG", 1);
   mmdb::Manager *mol = new mmdb::Manager;
   mol->ReadPDBASCII("frank.pdb"); // frankenstein, standard PDB file
				   // with added refmac LINKRs (refmac
				   // deletes HETATMs so we can't use
				   // the unmodified refmac results.
   mmdb::Residue *res_1 = coot::util::get_first_residue(mol);
   if (! res_1) {
      std::cout << "no res_1" << std::endl;
   } else { 
      coot::residue_spec_t specs[2];
      int selhnd = mol->NewSelection();
      specs[0] = coot::residue_spec_t("A", 121, "");
      specs[1] = coot::residue_spec_t("A", 200, "");
      for (unsigned int i=0; i<2; i++) {
	 specs[i].select_atoms(mol, selhnd, mmdb::SKEY_OR); 
      }

      mmdb::PPAtom atom_selection;
      int n_selected_atoms;
      mol->GetSelIndex(selhnd, atom_selection, n_selected_atoms);
      int imol = 0; // dummy
      std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v = 
	 coot::torsionable_bonds(imol, mol, atom_selection, n_selected_atoms, &t.geom);

      // tidy up
      mol->DeleteSelection(selhnd);
   }
   delete mol;

   return status;
}



int
test_read_prosmart_distance_restraints() {

   int status = 1;
   std::string file_name("ProSMART_Output/tutorial-modern.txt");

   int imol = read_pdb("test.pdb");
   add_refmac_extra_restraints(imol, file_name.c_str());

   return status;
}


#endif // BUILT_IN_TESTING

