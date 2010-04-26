/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
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

#include <string.h> // for strcmp


// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL

#include <algorithm> // for sort
#include <stdexcept>

#include "simple-restraint.hh"

//
#include "coot-coord-extras.hh"  // is_nucleotide_by_dict

// #include "mmdb.h" // for printing of CAtom pointers as info not raw
                     // pointers.  Removed. Too much (linking issues in)
                     // Makefile pain.

#include "coot-sysdep.h"

// #include "mmdb.h"

// iend_res is inclusive, so that 17,17 selects just residue 17.
//   have_disulfide_residues: other residues are included in the
//				residues_mol for disphide restraints.
// 
coot::restraints_container_t::restraints_container_t(int istart_res_in, int iend_res_in,
						     short int have_flanking_residue_at_start,
						     short int have_flanking_residue_at_end,
						     short int have_disulfide_residues,
						     const std::string &altloc,
						     const char *chain_id,
						     CMMDBManager *mol_in, 
						     const std::vector<coot::atom_spec_t> &fixed_atom_specs) {

   from_residue_vector = 0;
   lograma.init(LogRamachandran::All, 2.0, true);
   include_map_terms_flag = 0;
   init_from_mol(istart_res_in, iend_res_in, 
		 have_flanking_residue_at_start, 
		 have_flanking_residue_at_end,
		 have_disulfide_residues,
		 altloc,
		 chain_id, mol_in, fixed_atom_specs);
}

// Used in omega distortion graph
// 
coot::restraints_container_t::restraints_container_t(atom_selection_container_t asc_in,
						     const std::string &chain_id) { 
   from_residue_vector = 0;
   include_map_terms_flag = 0;
   lograma.init(LogRamachandran::All, 2.0, true);
   verbose_geometry_reporting = 0;
   mol = asc_in.mol;
   have_oxt_flag = 0;
   do_numerical_gradients_flag = 0;

   istart_res = 999999;
   iend_res = -9999999;

   PCResidue *SelResidues = NULL;
   int nSelResidues;

   // -------- Find the max and min res no -----------------------------
   int selHnd = mol->NewSelection();
   mol->Select(selHnd, STYPE_RESIDUE, 1,
	       (char *) chain_id.c_str(),
	       ANY_RES, "*",
	       ANY_RES, "*",
	       "*",  // residue name
	       "*",  // Residue must contain this atom name?
	       "*",  // Residue must contain this Element?
	       "*",  // altLocs
	       SKEY_NEW // selection key
	       );
   mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

   int resno;
   for (int ires=0; ires<nSelResidues; ires++) {
      resno = SelResidues[ires]->GetSeqNum();
      if (resno < istart_res)
	 istart_res = resno;
      if (resno > iend_res)
	 iend_res = resno;
   }
   mol->DeleteSelection(selHnd);
   // 
   // -------- Found the max and min res no -----------------------------

   // -------------------------------------------------------------------
   // Set class variables atom to the selection that includes the
   // chain (which is not the same as the input atom selection)
   //
   int SelHnd = mol->NewSelection();
   atom = NULL;
   mol->SelectAtoms(SelHnd, 0,
		    (char *) chain_id.c_str(),
		    ANY_RES, // starting resno, an int
		    "*", // any insertion code
		    ANY_RES, // ending resno
		    "*", // ending insertion code
		    "*", // any residue name
		    "*", // atom name
		    "*", // elements
		    "*"  // alt loc.
		    );
   mol->GetSelIndex(SelHnd, atom, n_atoms);

   // -------------------------------------------------------------------

   initial_position_params_vec.resize(3*n_atoms);
   for (int i=0; i<n_atoms; i++) {
      initial_position_params_vec[3*i  ] = atom[i]->x; 
      initial_position_params_vec[3*i+1] = atom[i]->y; 
      initial_position_params_vec[3*i+2] = atom[i]->z;
      // std::cout << "    " << i << "  " << coot::atom_spec_t(atom[i]) << "\n";
   }
}

coot::restraints_container_t::restraints_container_t(PCResidue *SelResidues, int nSelResidues,
						     const std::string &chain_id,
						     CMMDBManager *mol_in) { 
   
   include_map_terms_flag = 0;
   from_residue_vector = 0;
   lograma.init(LogRamachandran::All, 2.0, true);
   std::vector<coot::atom_spec_t> fixed_atoms_dummy;
   int istart_res = 999999;
   int iend_res = -9999999;
   int resno;
   
   for (int i=0; i<nSelResidues; i++) { 
      resno = SelResidues[i]->seqNum;
      if (resno < istart_res)
	 istart_res = resno;
      if (resno > iend_res)
	 iend_res = resno;
   }
   
   short int have_flanking_residue_at_start = 0;
   short int have_flanking_residue_at_end = 0;
   short int have_disulfide_residues = 0;
   char *chn = (char *) chain_id.c_str();

   std::cout << "DEBUG:  ==== istart_res iend_res " << istart_res << " " << iend_res << std::endl; 

   init_from_mol(istart_res, iend_res, 
		 have_flanking_residue_at_start,
		 have_flanking_residue_at_end,
		 have_disulfide_residues, 
		 std::string(""), chn, mol_in, fixed_atoms_dummy);

}

coot::restraints_container_t::restraints_container_t(int istart_res_in, int iend_res_in,
						     short int have_flanking_residue_at_start,
						     short int have_flanking_residue_at_end,
						     short int have_disulfide_residues,
						     const std::string &altloc,
						     const char *chain_id,
						     CMMDBManager *mol, // const in an ideal world
						     const std::vector<coot::atom_spec_t> &fixed_atom_specs,
						     const clipper::Xmap<float> &map_in,
						     float map_weight_in) {

   from_residue_vector = 0;
   lograma.init(LogRamachandran::All, 2.0, true);
   init_from_mol(istart_res_in, iend_res_in, 		 
		 have_flanking_residue_at_start, 
		 have_flanking_residue_at_end,
		 have_disulfide_residues,
		 altloc,
		 chain_id, mol, fixed_atom_specs);
   map = map_in;
   map_weight = map_weight_in;
   include_map_terms_flag = 1;

}

// 20081106 construct from a vector of residues, each of which
// has a flag attached that denotes whether or not it is a fixed
// residue (it would be set, for example in the case of flanking
// residues).
coot::restraints_container_t::restraints_container_t(const std::vector<std::pair<bool,CResidue *> > &residues,
						     const coot::protein_geometry &geom,
						     CMMDBManager *mol,
						     const std::vector<atom_spec_t> &fixed_atom_specs) {
   from_residue_vector = 1;
   lograma.init(LogRamachandran::All, 2.0, true);
   init_from_residue_vec(residues, geom, mol, fixed_atom_specs);
   include_map_terms_flag = 0;
}



// What are the rules for dealing with alt conf in flanking residues?
// 
// First we want to try to select only atom that have the same alt
// conf, if that fails to select atoms in flanking residues (and there
// are atoms in flanking residues (which we know from
// have_flanking_residue_at_end/start), then we should try an mmdb construction [""|"A"]
// (i.e. blank or "A").  If that fails, give up - it's a badly formed pdb file (it 
// seems to me).
// 
void
coot::restraints_container_t::init_from_mol(int istart_res_in, int iend_res_in,
					    short int have_flanking_residue_at_start,
					    short int have_flanking_residue_at_end,
					    short int have_disulfide_residues,
					    const std::string &altloc,
					    const char *chain_id,
					    CMMDBManager *mol_in, 
					    const std::vector<coot::atom_spec_t> &fixed_atom_specs) {

   init_shared_pre(mol_in);

   istart_res = istart_res_in;
   iend_res   = iend_res_in;
   chain_id_save = chain_id;

   // internal flags that mirror passed have_flanking_residue_at_* variables
   istart_minus_flag = have_flanking_residue_at_start;
   iend_plus_flag    = have_flanking_residue_at_end;

   int iselection_start_res = istart_res;
   int iselection_end_res   = iend_res;
   // std::cout << "start res range: " << istart_res << " " << iend_res << " " << chain_id << "\n";

   // Are the flanking atoms available in mol_in?  (mol_in was
   // constrcted outside so the mol_in constructing routine know if
   // they were there or not.
   // 
   if (have_flanking_residue_at_start) iselection_start_res--;
   if (have_flanking_residue_at_end)   iselection_end_res++;

   SelHnd_atom = mol->NewSelection();
   mol->SelectAtoms(SelHnd_atom,
		    0,
		    chain_id,
		    iselection_start_res, "*",
		    iselection_end_res, "*",
		    "*", // rnames
		    "*", // anames
		    "*", // elements
		    "*"  // altLocs 
		    );

   // set the PPCAtom atom (class variable) and n_atoms:
   // 
   mol->GetSelIndex(SelHnd_atom, atom, n_atoms);

   bool debug = 0;
   if (debug) { 
      std::cout << "DEBUG:: Selecting residues in chain " << chain_id << " gives " << n_atoms
		<< " atoms " << std::endl;
      for (int iat=0; iat<n_atoms; iat++) {
	 std::cout << "   " << iat << " " << atom[iat]->name << " "  << atom[iat]->GetSeqNum()
		   << " " << atom[iat]->GetChainID() << std::endl;
      }
   }

   // debugging, you need to uncomment mmdb.h at the top and add coords to the linking
   // atom_selection_container_t tmp_res_asc = make_asc(mol_in);
   //    std::cout << "There are " << tmp_res_asc.n_selected_atoms
   //              << " atoms in tmp_res_asc\n";
//    for (int kk=0; kk<tmp_res_asc.n_selected_atoms; kk++) { 
//       std::cout << "In simple rest " << kk << " "
//                 << tmp_res_asc.atom_selection[kk] << "\n";
//    } 

//    std::cout << "INFO::" << n_atoms << " atoms selected from molecule for refinement" ;
//    std::cout << " (this includes fixed and flanking atoms)." << std::endl;

   if (n_atoms == 0) { 
      std::cout << "ERROR:: atom selection disaster:" << std::endl;
      std::cout << "this should not happen" << std::endl;
      std::cout << "res range: " << iselection_start_res << " " 
		<< iselection_end_res << " " << chain_id << " " 
		<< "flanking flags: " << have_flanking_residue_at_start 
		<< " " << have_flanking_residue_at_end << std::endl;
   }

   init_shared_post(fixed_atom_specs);
}

void
coot::restraints_container_t::init_shared_pre(CMMDBManager *mol_in) {

   do_numerical_gradients_flag = 0;
   verbose_geometry_reporting = 0;
   have_oxt_flag = 0; // set in mark_OXT()
   mol = mol_in;
} 

void
coot::restraints_container_t::init_shared_post(const std::vector<atom_spec_t> &fixed_atom_specs) {


   bonded_atom_indices.resize(n_atoms);

   initial_position_params_vec.resize(3*n_atoms); 
   for (int i=0; i<n_atoms; i++) {
      initial_position_params_vec[3*i  ] = atom[i]->x; 
      initial_position_params_vec[3*i+1] = atom[i]->y; 
      initial_position_params_vec[3*i+2] = atom[i]->z; 
   }

   // Set the UDD have_bond_or_angle to initally all "not".  They get
   // set to "have" (1) in make_restraints (and functions thereof).
   // 
   // udd_handle becomes member data so that it can be used in
   // make_restraints() without passing it back (this function is part
   // of a constructor, don't forget).
   // 
   udd_bond_angle = mol->RegisterUDInteger (UDR_ATOM, "bond or angle");
   if (udd_bond_angle < 0) { 
     std::cout << "ERROR:: can't make udd_handle in init_from_mol\n";
   } else { 
     for (int i=0; i<n_atoms; i++) {
       atom[i]->PutUDData(udd_bond_angle,0);
     }
   }

   // Set the UDD of the indices in the atom array (i.e. the thing
   // that get_asc_index returns)
   udd_atom_index_handle = mol->RegisterUDInteger ( UDR_ATOM, "atom_array_index");
   if (udd_atom_index_handle < 0) { 
     std::cout << "ERROR:: can't make udd_handle in init_from_mol\n";
   } else { 
     for (int i=0; i<n_atoms; i++) {
       atom[i]->PutUDData(udd_atom_index_handle,i);
     }
   }

   use_map_gradient_for_atom.resize(n_atoms,0);
   if (! from_residue_vector) {
      // convential way
      for (int i=0; i<n_atoms; i++) {
	 if (atom[i]->residue->seqNum >= istart_res &&
	     atom[i]->residue->seqNum <= iend_res) {
	    use_map_gradient_for_atom[i] = 1;
	 } else {
	    use_map_gradient_for_atom[i] = 0;
	 }
      }
   } else {
      // blank out the non moving atoms (i.e. flanking residues)
      for (int i=0; i<n_atoms; i++) {
	 CResidue *res_p = atom[i]->residue;
	 if (is_a_moving_residue_p(res_p)) {
	    use_map_gradient_for_atom[i] = 1;
	 } else { 
	    use_map_gradient_for_atom[i] = 0;
	 }
      }
   }

   // z weights
   atom_z_weight.resize(n_atoms);
   std::vector<std::pair<std::string, int> > atom_list = coot::util::atomic_number_atom_list();
   for (int i=0; i<n_atoms; i++) {
      double z = coot::util::atomic_number(atom[i]->element, atom_list);
      if (z < 0.0) {
	 std::cout << "Unknown element :" << atom[i]->element << ": " << std::endl;
	 z = 6.0; // as for carbon
      } 
      atom_z_weight[i] = z;
   }
   
   // similarly for the fixed atoms:   
   // 
   assign_fixed_atom_indices(fixed_atom_specs); // convert from std::vector<CAtom *>
   				                // to std::vector<int> fixed_atom_indices;

   // blank out those atoms from seeing electron density map gradients
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      use_map_gradient_for_atom[fixed_atom_indices[ifixed]] = 0;
   } 
   
   if (verbose_geometry_reporting)
      for (int i=0; i<n_atoms; i++)
	 std::cout << atom[i]->name << " " << atom[i]->residue->seqNum << " "
		   << use_map_gradient_for_atom[i] << std::endl;

} 

void
coot::restraints_container_t::init_from_residue_vec(const std::vector<std::pair<bool,CResidue *> > &residues,
						    const coot::protein_geometry &geom,
						    CMMDBManager *mol,
						    const std::vector<atom_spec_t> &fixed_atom_specs) {

   init_shared_pre(mol);
   residues_vec = residues;

   // Need to set class members PPCAtom atom and int n_atoms.
   // ...
   // 20090620: or do we?

   // debug:
   if (0) { 
      for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
	 std::cout << "debug:: =============== in init_from_residue_vec() residue "
		   << ir << " of " << residues_vec.size() << " "
		   << residues_vec[ir].second << std::endl;
	 PCAtom *res_atom_selection = NULL;
	 int n_res_atoms;
	 residues_vec[ir].second->GetAtomTable(res_atom_selection, n_res_atoms);
	 std::cout << "debug:: =============== in init_from_residue_vec() residue "
		   << ir << " of " << residues_vec.size() << " has : "
		   << n_res_atoms << " atom " << std::endl;
	 std::cout << "debug:: =============== in init_from_residue_vec() residue "
		   << ir << " of " << residues_vec.size() << " seqnum: "
		   << residues_vec[ir].second->GetSeqNum() << " chainid: "
		   << residues_vec[ir].second->GetChainID() << std::endl;
	 for (int iat=0; iat<n_res_atoms; iat++) {
	    CAtom *at =  res_atom_selection[iat];
	    std::cout << "DEBUG:: in init_from_residue_vec: atom "
		      << iat << " of " << n_res_atoms << " \"" 
		      << at->name << "\" \"" << at->altLoc << "\" " 
		      << at->GetSeqNum() << " \"" << at->GetInsCode() << "\" \"" 
		      << at->GetChainID() << "\"" << std::endl;
	 }
      }
   }
   

   // what about adding the flanking residues?  How does the atom
   // indexing of that work when (say) adding a bond?
   coot::bonded_pair_container_t bpc = bonded_flanking_residues_by_residue_vector(geom);

   //    std::cout << "   DEBUG:: made " << bpc.size() << " bonded flanking pairs " << std::endl;

   // passed and flanking
   std::vector<CResidue *> all_residues;
   for (unsigned int i=0; i<residues.size(); i++) {
      all_residues.push_back(residues[i].second);
   }

   // Include only the fixed residues, because they are the flankers,
   // the other residues are the ones in the passed residues vector.
   // We don't have members of bonded_pair_container_t that are both
   // fixed.
   int n_bonded_flankers_in_total = 0; // debug/info counter
   for (unsigned int i=0; i<bpc.size(); i++) {
      if (bpc[i].is_fixed_first) { 
	 all_residues.push_back(bpc[i].res_1);
	 n_bonded_flankers_in_total++;
      } 
      if (bpc[i].is_fixed_second) { 
	 all_residues.push_back(bpc[i].res_2);
	 n_bonded_flankers_in_total++;
      } 
   }

//    std::cout << "   DEBUG:: There are " << residues.size() << " passed residues and "
// 	     << all_residues.size() << " residues total (including flankers)"
// 	     << std::endl;

   n_atoms = 0;
   for (unsigned int i=0; i<all_residues.size(); i++) {
      n_atoms += all_residues[i]->GetNumberOfAtoms();
   }
   std::cout << "   DEUBG:: There are " << n_atoms << " atoms total (including flankers)"
	     << std::endl;

   atom = new PCAtom[n_atoms];
   int atom_index = 0;
   for (unsigned int i=0; i<all_residues.size(); i++) {
      PPCAtom residue_atoms = 0;
      int n_res_atoms;
      all_residues[i]->GetAtomTable(residue_atoms, n_res_atoms);
      for (int iat=0; iat<n_res_atoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 atom[atom_index] = at;
	 atom_index++;
      }
   }

   init_shared_post(fixed_atom_specs); // use n_atoms
} 



void
coot::restraints_container_t::assign_fixed_atom_indices(const std::vector<coot::atom_spec_t> &fixed_atom_specs) {

   fixed_atom_indices.clear();
//    std::cout << "Finding atom indices for " << fixed_atom_specs.size()
// 	     << " fixed atoms " << std::endl;
   for (unsigned int i=0; i<fixed_atom_specs.size(); i++) {
      for (int iat=0; iat<n_atoms; iat++) {
	 if (fixed_atom_specs[i].matches_spec(atom[iat])) {
	    fixed_atom_indices.push_back(iat);
	 }
      }
   }
   //    std::cout << "Found indices for " << fixed_atom_indices.size()
   // << " fixed atoms" << std::endl;
}

// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
// 
coot::refinement_results_t
coot::restraints_container_t::minimize(restraint_usage_Flags usage_flags) {

   short int print_chi_sq_flag = 1;
   return minimize(usage_flags, 1000, print_chi_sq_flag);

}

 
// return success: GSL_ENOPROG, GSL_CONTINUE, GSL_ENOPROG (no progress)
// 
coot::refinement_results_t
coot::restraints_container_t::minimize(restraint_usage_Flags usage_flags, 
				       int nsteps_max,
				       short int print_initial_chi_sq_flag) {

   restraints_usage_flag = usage_flags;
   
   const gsl_multimin_fdfminimizer_type *T;
   gsl_multimin_fdfminimizer *s;


   // check that we have restraints before we start to minimize:
   if (restraints_vec.size() == 0) {
      if (restraints_usage_flag != NO_GEOMETRY_RESTRAINTS) {
      std::cout << "SPECIFICATION ERROR:  There are no restraints. ";
      std::cout << "No minimization will happen" << std::endl;
      return coot::refinement_results_t(0, 0, "No Restraints!");
      }
   } 
   
   setup_gsl_vector_variables();  //initial positions 


   setup_multimin_func(); // provide functions for f, df, fdf

   // T = gsl_multimin_fdfminimizer_conjugate_fr; // not as good as pr
   // T = gsl_multimin_fdfminimizer_steepest_descent; // pathetic
   // T = gsl_multimin_fminimizer_nmsimplex; // you can't just drop this in, 
                                             // because simplex is a 
                                             // non-gradient method. 
   // T = gsl_multimin_fdfminimizer_vector_bfgs;
   T = gsl_multimin_fdfminimizer_conjugate_pr;

   //   This is the Polak-Ribiere conjugate gradient algorithm. It is
   //   similar to the Fletcher-Reeves method, differing only in the
   //   choice of the coefficient \beta. Both methods work well when
   //   the evaluation point is close enough to the minimum of the
   //   objective function that it is well approximated by a quadratic
   //   hypersurface.

   s = gsl_multimin_fdfminimizer_alloc (T, n_variables());

   // x and multimin_func are class variables
   // 
   // if (print_initial_chi_sq_flag) 
      // std::cout << "sizes: n_variables() " << n_variables() 
		// << " s: " << s->x->size 
		// << " x: " << x->size
		// << " and fdf->n: " << multimin_func.n << std::endl;
   
   // info(); 
   
   // restraints_usage_flag = BONDS; 
   // restraints_usage_flag = BONDS_AND_ANGLES; 
   // restraints_usage_flag = BONDS_ANGLES_AND_TORSIONS;
   // restraints_usage_flag = BONDS_ANGLES_TORSIONS_AND_PLANES;
   // restraints_usage_flag = BONDS_ANGLES_AND_PLANES;

   // restraints_usage_flag = BONDS_MASK; 
   // restraints_usage_flag = ANGLES_MASK; 
   // restraints_usage_flag = TORSIONS_MASK;
   // restraints_usage_flag = PLANES_MASK;
   
   
   // std::cout << "DEBUG:: in minimize, restraints_usage_flag is " 
   // << restraints_usage_flag << std::endl;

   // We get ~1ms/residue with bond and angle terms and no density terms.
   // 
   // for (int i=0; i<100; i++) { // time testing

   gsl_multimin_fdfminimizer_set (s, &multimin_func, x, 0.02, 0.03);
   // stepsize 0.01, tolerance 1e-4


   if (print_initial_chi_sq_flag) { 
      double d = coot::distortion_score(x,(double *)this); 
      std::cout << "initial distortion_score: " << d << std::endl; 
      chi_squareds("Initial Chi Squareds", s->x);
   }

//      std::cout << "pre minimization atom positions\n";
//      for (int i=0; i<n_atoms*3; i+=3)
//         std::cout << i/3 << " ("
// 		  << gsl_vector_get(x, i) << ", "
// 		  << gsl_vector_get(x, i+1) << ", "
// 		  << gsl_vector_get(x, i+2) << ")"
// 		  << std::endl;

   // fix_chiral_atoms_maybe(s->x);
   // fix_chiral_atoms_maybe(s->x);

//     std::cout << "Post chiral fixing\n";
//     for (int i=0; i<n_atoms*3; i++)
//        std::cout << gsl_vector_get(x, i) << std::endl;

   size_t iter = 0; 
   int status;
   std::vector<coot::refinement_lights_info_t> lights_vec;
   do
      {
	 iter++;
	 // std::cout << "debug:: iteration number " << iter << std::endl;
	 status = gsl_multimin_fdfminimizer_iterate (s);

	 if (status) { 
	    cout << "unexpected error from gsl_multimin_fdfminimizer_iterate"
		 << endl;
	    if (status == GSL_ENOPROG) {
	       cout << "Error in gsl_multimin_fdfminimizer_iterate was GSL_ENOPROG"
		    << endl; 
	       chi_squareds("Final Estimated Average Z scores", s->x);
	    }
	    break;
	 } 

	 // back of envelope calculation suggests g_crit = 0.1 for
	 // coordinate shift of 0.001:  So let's choose 0.05
	 status = gsl_multimin_test_gradient (s->gradient, 0.2);

	 if (status == GSL_SUCCESS) { 
	    std::cout << "Minimum found (iteration number " << iter << ") at ";
	    std::cout << s->f << "\n";
	    std::vector <coot::refinement_lights_info_t> results = 
	       chi_squareds("Final Estimated Average Z Scores:", s->x);
	    lights_vec = results;
	 }

	 if (verbose_geometry_reporting)
	    cout << "iteration number " << iter << " " << s->f << endl; 

      }
   while ((status == GSL_CONTINUE) && (int(iter) < nsteps_max));

   // std::cout << "Post refinement fixing\n";
   // fix_chiral_atoms_maybe(s->x);

   // } time testing
   update_atoms(s->x); // do OXT here

   gsl_multimin_fdfminimizer_free (s);
   gsl_vector_free (x);
   // (we don't get here unless restraints were found)
   coot::refinement_results_t rr(1, status, lights_vec);
//    std::cout << "DEBUG:: returning from minimize() :" << results_string
// 	     << ":" << std::endl;
   
   return rr;
}


coot::geometry_distortion_info_container_t
coot::restraints_container_t::geometric_distortions(coot::restraint_usage_Flags flags) {

   restraints_usage_flag = flags;
   setup_gsl_vector_variables();  //initial positions in x array
   coot::geometry_distortion_info_container_t dv = distortion_vector(x);
   
   return dv;
}

std::ostream&
coot::operator<<(std::ostream &s, geometry_distortion_info_container_t gdic) {

   s << "[ chain :" << gdic.chain_id << ": residues " << gdic.min_resno << " to " << gdic.max_resno
     << " residues: \n" ;
   for (unsigned int ires=0; ires<gdic.geometry_distortion.size(); ires++)
      s << "      " << gdic.geometry_distortion[ires] << "\n";
   s << "]\n";
   return s;
} 

std::ostream&
coot::operator<<(std::ostream &s, geometry_distortion_info_t gdi) {

   if (gdi.set) {
      s << gdi.restraint << " " << gdi.residue_spec << " distortion: " << gdi.distortion_score;
   } else {
      s << "{geometry_distortion_info-unset}";
   } 
   return s;
}


std::ostream &
coot::operator<<(std::ostream &s, simple_restraint r) {

   s << "{restraint: ";
   if (r.restraint_type == coot::BOND_RESTRAINT)
      s << "Bond   ";
   if (r.restraint_type == coot::ANGLE_RESTRAINT)
      s << "Angle  ";
   if (r.restraint_type == coot::TORSION_RESTRAINT)
      s << "Torsion";
   if (r.restraint_type == coot::PLANE_RESTRAINT)
      s << "Plane  ";
   if (r.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT)
      s << "NBC    ";
   if (r.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT)
      s << "Chiral ";
   if (r.restraint_type == coot::RAMACHANDRAN_RESTRAINT)
      s << "Rama   ";
   s << "}";
   return s;
}
coot::omega_distortion_info_container_t
coot::restraints_container_t::omega_trans_distortions(int mark_cis_peptides_as_bad_flag) {

   // restraints_usage_flag = flags;  // not used?
   setup_gsl_vector_variables();  //initial positions in x array
   std::string chain_id("");

   if (n_atoms > 0)
      chain_id = atom[0]->GetChainID();
   else
      chain_id = "blank"; // shouldn't happen.


   CChain *chain_p = atom[0]->GetChain();
   // I think there will be need for some sort of scaling thing here.
   std::pair<short int, int> minr = coot::util::min_resno_in_chain(chain_p);
   std::pair<short int, int> maxr = coot::util::max_resno_in_chain(chain_p);
   int min_resno = 0;
   int max_resno = 0;
   if (minr.first && maxr.first) {
      min_resno = minr.second;
      max_resno = maxr.second;
   }
   double scale = 1.2; // how much to scale the omega difference by so that
		       // it looks good on the graph.

   if (!mark_cis_peptides_as_bad_flag)
      scale = 2.5;

   coot::omega_distortion_info_container_t dc(chain_id, min_resno, max_resno);

   // we need to pick out the CA and C of this and N and Ca of next.

   // now add data to dc.omega_distortions.
   PCResidue *first = NULL;
   PCResidue *second = NULL;
   int nfirst, nnext;
//    int ifirst;
//    int inext;

   for (int i=istart_res; i<iend_res; i++) {
      int selHnd1 = mol->NewSelection();
      mol->Select(selHnd1, STYPE_RESIDUE, 1,
		  (char *) chain_id.c_str(),
		  i, "*",
		  i, "*",
		  "*",  // residue name
		  "*",  // Residue must contain this atom name?
		  "*",  // Residue must contain this Element?
		  "*",  // altLocs
		  SKEY_NEW // selection key
		  );
      mol->GetSelIndex(selHnd1, first, nfirst);

      int selHnd2 = mol->NewSelection();
      mol->Select(selHnd2, STYPE_RESIDUE, 1,
		  (char *) chain_id.c_str(),
		  i+1, "*",
		  i+1, "*",
		  "*",  // residue name
		  "*",  // Residue must contain this atom name?
		  "*",  // Residue must contain this Element?
		  "*",  // altLocs
		  SKEY_NEW // selection key
		  );
      mol->GetSelIndex(selHnd2, second, nnext);


      if ((nfirst > 0) && (nnext > 0)) {

	 if (! first[0]->chain->isSolventChain()) { 

	    // So we have atoms selected in both residues, lets look for those atoms:

	    CAtom *at;
	    clipper::Coord_orth ca_first, c_first, n_next, ca_next;
	    short int got_ca_first = 0, got_c_first = 0, got_n_next = 0, got_ca_next = 0;
	    PPCAtom res_selection = NULL;
	    int i_no_res_atoms;

	    first[0]->GetAtomTable(res_selection, i_no_res_atoms);
	    if (i_no_res_atoms > 0) {
	       for (int iresatom=0; iresatom<i_no_res_atoms; iresatom++) {
		  at = res_selection[iresatom];
		  std::string atom_name(at->name);
		  if (atom_name == " CA ") {
		     ca_first = clipper::Coord_orth(at->x, at->y, at->z);
		     got_ca_first = 1;
		  }
		  if (atom_name == " C  ") {
		     c_first = clipper::Coord_orth(at->x, at->y, at->z);
		     got_c_first = 1;
		  }
	       }
	    }
	    second[0]->GetAtomTable(res_selection, i_no_res_atoms);
	    if (i_no_res_atoms > 0) {
	       for (int iresatom=0; iresatom<i_no_res_atoms; iresatom++) {
		  at = res_selection[iresatom];
		  std::string atom_name(at->name);
		  if (atom_name == " CA ") {
		     ca_next = clipper::Coord_orth(at->x, at->y, at->z);
		     got_ca_next = 1;
		  }
		  if (atom_name == " N  ") {
		     n_next = clipper::Coord_orth(at->x, at->y, at->z);
		     got_n_next = 1;
		  }
	       }
	    }

	    if (got_ca_first && got_c_first && got_n_next && got_ca_next) {
	       double tors = clipper::Coord_orth::torsion(ca_first, c_first, n_next, ca_next);
	       double torsion = clipper::Util::rad2d(tors);
	       torsion = (torsion > 0.0) ? torsion : 360.0 + torsion;
	       std::string info = coot::util::int_to_string(i);
	       info += chain_id;
	       info += " ";
	       info += first[0]->name;
	       info += " Omega: ";
	       info += coot::util::float_to_string(torsion);
	       double distortion = fabs(180.0 - torsion);
	       // distortion = (distortion < 60.0) ? distortion : 60.0;

	       if (!mark_cis_peptides_as_bad_flag)
		  // consider the cases: torsion=1 and torsion=359
		  if (distortion > 90.0) {
		     distortion = torsion;
		     if (distortion > 180.0) {
			distortion -= 360;
			distortion = fabs(distortion);
		     }
		  }
	       distortion *= scale;
	       dc.omega_distortions.push_back(coot::omega_distortion_info_t(i, distortion, info));
	    } else {
	       std::cout << "INFO:: failed to get all atoms for omega torsion "
			 << "chain " << first[0]->GetChainID() << " residues "
			 << first[0]->GetSeqNum() << " to " << second[0]->GetSeqNum()
			 << std::endl;
	    }
	 }
      }
      mol->DeleteSelection(selHnd1);
      mol->DeleteSelection(selHnd2);
   }
   return dc;
} 

void
coot::restraints_container_t::adjust_variables(const atom_selection_container_t &asc) { 



}

// 
double
starting_structure_diff_score(const gsl_vector *v, void *params) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   double d;
   double dist = 0; 

   // if (v->size != restraints->init_positions_size() ) {

   for (int i=0; i<restraints->init_positions_size(); i++) { 
      d = restraints->initial_position(i) - gsl_vector_get(v, i);
      dist += 0.01*d*d;
   }
   cout << "starting_structure_diff_score: " << dist << endl; 
   return dist; 
}

// Return the distortion score.
// 
double
coot::distortion_score(const gsl_vector *v, void *params) { 
   
   // so we are comparing the geometry of the value in the gsl_vector
   // v and the ideal values.
   // 

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   double distortion = 0; 


   // distortion += starting_structure_diff_score(v, params); 

   // tmp debugging stuff
   double nbc_diff = 0.0;
   double d;

   for (int i=0; i<restraints->size(); i++) {
      if (restraints->restraints_usage_flag & coot::BONDS_MASK) { // 1: bonds
	 if ( (*restraints)[i].restraint_type == coot::BOND_RESTRAINT) {
  	    // cout << "adding an bond restraint " << i << endl;
	    distortion += coot::distortion_score_bond((*restraints)[i], v);
	 }
      }

      if (restraints->restraints_usage_flag & coot::ANGLES_MASK) { // 2: angles
	 if ( (*restraints)[i].restraint_type == coot::ANGLE_RESTRAINT) {
  	    // cout << "adding an angle restraint " << i << endl;
	    distortion += coot::distortion_score_angle((*restraints)[i], v);
	 }
      }

      if (restraints->restraints_usage_flag & coot::TORSIONS_MASK) { // 4: torsions
	 if ( (*restraints)[i].restraint_type == coot::TORSION_RESTRAINT) {
	    // cout << "adding an torsion restraint: number " << i << endl;  
	    // std::cout << "distortion sum pre-adding a torsion: " << distortion << std::endl;
	    distortion += coot::distortion_score_torsion((*restraints)[i], v); 
	    // std::cout << "distortion sum post-adding a torsion: " << distortion << std::endl;
	 }
      }

      if (restraints->restraints_usage_flag & coot::PLANES_MASK) { // 8: planes
	 if ( (*restraints)[i].restraint_type == coot::PLANE_RESTRAINT) {
	    // 	    std::cout << "adding an plane restraint " << i << std::endl;  
	    // std::cout << "distortion sum pre-adding a torsion: " << distortion << std::endl;
	    distortion += coot::distortion_score_plane((*restraints)[i], v); 
	    // std::cout << "distortion sum post-adding a torsion: " << distortion << std::endl;
	 }
      }

      if (restraints->restraints_usage_flag & coot::NON_BONDED_MASK) { // 16: 
	 if ( (*restraints)[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) { 
	    // std::cout << "adding a non_bonded restraint " << i << std::endl;
	    // std::cout << "distortion sum pre-adding  " << distortion << std::endl;
	    d = coot::distortion_score_non_bonded_contact( (*restraints)[i], v);
	    distortion += d;
	    // 	    std::cout << "distortion sum post-adding " << distortion << std::endl;
	    nbc_diff += d;
	 }
      }

      if (restraints->restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) { 
	 // if (0) { 
	 
   	 if ( (*restraints)[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) { 
   	    d = coot::distortion_score_chiral_volume( (*restraints)[i], v);
	    // std::cout << "DEBUG:: distortion for chiral: " << d << std::endl;
   	    distortion += d;
   	 }
      }

      if (restraints->restraints_usage_flag & coot::RAMA_PLOT_MASK) {
   	 if ( (*restraints)[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) { 
   	    d = coot::distortion_score_rama( (*restraints)[i], v, restraints->LogRama());
	    // std::cout << "DEBUG:: distortion for rama: " << d << std::endl;
   	    distortion += d; // positive is bad...  negative is good.
   	 }
      }
   }

//     std::cout << "nbc_diff   distortion: " << nbc_diff << std::endl;
//     std::cout << "post-terms distortion: " << distortion << std::endl;

   if ( restraints->include_map_terms() )
      distortion += coot::electron_density_score(v, params); // good map fit: low score
   
   // cout << "distortion (in distortion_score): " << distortion << endl; 
   return distortion; 
}

coot::geometry_distortion_info_container_t
coot::restraints_container_t::distortion_vector(const gsl_vector *v) const {

   std::string chainid;
   if (n_atoms > 0)
      chainid = atom[0]->GetChainID();
   else
      chainid = "blank";

   coot::geometry_distortion_info_container_t distortion_vec_container(atom, n_atoms, chainid);
   double distortion = 0.0;
   int atom_index = -1; // initial unset, the atom index used to spec the residue.

   for (unsigned int i=0; i<restraints_vec.size(); i++) {
      if (restraints_usage_flag & coot::BONDS_MASK) 
	 if (restraints_vec[i].restraint_type == coot::BOND_RESTRAINT) { 
	    distortion = coot::distortion_score_bond(restraints_vec[i], v);
	    atom_index = restraints_vec[i].atom_index_1;
	 } 

      if (restraints_usage_flag & coot::ANGLES_MASK) 
	 if (restraints_vec[i].restraint_type == coot::ANGLE_RESTRAINT) { 
	    distortion = coot::distortion_score_angle(restraints_vec[i], v);
	    atom_index = restraints_vec[i].atom_index_1;
	 } 

      if (restraints_usage_flag & coot::TORSIONS_MASK)
	 if (restraints_vec[i].restraint_type == coot::TORSION_RESTRAINT) { 
	    distortion = coot::distortion_score_torsion(restraints_vec[i], v);
	    atom_index = restraints_vec[i].atom_index_1;
	 } 

      if (restraints_usage_flag & coot::PLANES_MASK) 
	 if (restraints_vec[i].restraint_type == coot::PLANE_RESTRAINT) { 
	    distortion = coot::distortion_score_plane(restraints_vec[i], v);
	    atom_index = restraints_vec[i].atom_index[0];
	 } 

      if (restraints_usage_flag & coot::NON_BONDED_MASK)  
	 if (restraints_vec[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) { 
	    distortion = coot::distortion_score_non_bonded_contact(restraints_vec[i], v);
	    atom_index = restraints_vec[i].atom_index_1;
	 }

      if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK)
   	 if (restraints_vec[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
   	    distortion = coot::distortion_score_chiral_volume(restraints_vec[i], v);
	    atom_index = restraints_vec[i].atom_index_centre;
	 } 

      if (restraints_usage_flag & coot::RAMA_PLOT_MASK) 
    	 if (restraints_vec[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) { 
	    distortion = coot::distortion_score_rama(restraints_vec[i], v, lograma);
	    atom_index = restraints_vec[i].atom_index_1;
	 } 

      if (atom_index != -1) {
	 coot::residue_spec_t rs(atom[atom_index]->GetResidue());
	 coot::geometry_distortion_info_t gdi(distortion, restraints_vec[i], rs);
	 distortion_vec_container.geometry_distortion.push_back(gdi);
      }
   }


   // Find max_resno and min_resno.
   // 
   // Notice we only use BOND_RESTRAINT to do this (atom_index_1 is
   // not defined in plane restraint).
   int max_resno = -9999999;
   int min_resno = 999999999;
   int idx1, idx2;
   int this_resno1, this_resno2;
   for (unsigned int i=0; i<distortion_vec_container.geometry_distortion.size(); i++) {
      if (restraints_usage_flag & coot::BONDS_MASK) 
	 if (restraints_vec[i].restraint_type == coot::BOND_RESTRAINT) {
	    idx1 = distortion_vec_container.geometry_distortion[i].restraint.atom_index_1;
	    idx2 = distortion_vec_container.geometry_distortion[i].restraint.atom_index_2;
	    
	    this_resno1 = distortion_vec_container.atom[idx1]->GetSeqNum();
	    this_resno2 = distortion_vec_container.atom[idx2]->GetSeqNum();
	    if (this_resno1 < min_resno)
	       min_resno = this_resno1;
	    if (this_resno2 < min_resno)
	       min_resno = this_resno2;
	    if (this_resno1 > max_resno) 
	       max_resno = this_resno1;
	    if (this_resno2 > max_resno) 
	       max_resno = this_resno2;
	 }
   }
   distortion_vec_container.set_min_max(min_resno, max_resno);

   return distortion_vec_container;
}

void
coot::restraints_container_t::fix_chiral_atoms_maybe(gsl_vector *s) {

   if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) {
      for(unsigned int i=0; i<restraints_vec.size(); i++) {
	 if ( restraints_vec[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	    fix_chiral_atom_maybe(restraints_vec[i], s);
	 }
      }
   }
}

double
coot::distortion_score_bond(const coot::simple_restraint &bond_restraint,
			    const gsl_vector *v) {

   // cout << "adding a bond restraint" << endl;
   int idx = 3*(bond_restraint.atom_index_1 - 0); 
   clipper::Coord_orth a1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(bond_restraint.atom_index_2 - 0); 
   clipper::Coord_orth a2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   
   double weight = 1/(bond_restraint.sigma * bond_restraint.sigma);
   // std::cout << "   weight is " << weight << std::endl; 
   double bit = 
      (clipper::Coord_orth::length(a1,a2) - bond_restraint.target_value);
   
   return weight * bit *bit;
}

double
coot::distortion_score_angle(const coot::simple_restraint &angle_restraint,
			     const gsl_vector *v) {
   

   int idx = 3*(angle_restraint.atom_index_1); 
   clipper::Coord_orth a1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(angle_restraint.atom_index_2); 
   clipper::Coord_orth a2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(angle_restraint.atom_index_3); 
   clipper::Coord_orth a3(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   clipper::Coord_orth d1 = a1 - a2; 
   clipper::Coord_orth d2 = a3 - a2; 
   double len1 = clipper::Coord_orth::length(a1,a2);
   double len2 = clipper::Coord_orth::length(a3,a2);
   len1 = len1 > 0.01 ? len1 : 0.01; 
   len2 = len2 > 0.01 ? len2 : 0.01; 
   double theta = acos(clipper::Coord_orth::dot(d1,d2)/(len1*len2));
   double bit = 
      clipper::Util::rad2d(theta) - angle_restraint.target_value;
   double weight = 1/(angle_restraint.sigma * angle_restraint.sigma);
   // 	    std::cout << "actual: " << clipper::Util::rad2d(theta) << " target: "
   // 		      << (*restraints)[i].target_value
   // 		      << " adding distortion score: "
   // 		      << weight * pow(bit, 2.0) << std::endl;
   return weight * bit * bit;
}


//
// Return the distortion score from a single torsion restraint.
double
coot::distortion_score_torsion(const coot::simple_restraint &torsion_restraint,
			       const gsl_vector *v) {

   // First calculate the torsion:
   // theta = arctan(E/G); 
   // where E = a.(bxc) and G = -a.c + (a.b)(b.c)


   int idx; 

   idx = 3*(torsion_restraint.atom_index_1);
   clipper::Coord_orth P1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(torsion_restraint.atom_index_2); 
   clipper::Coord_orth P2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(torsion_restraint.atom_index_3); 
   clipper::Coord_orth P3(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(torsion_restraint.atom_index_4); 
   clipper::Coord_orth P4(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));

//    P1 = clipper::Coord_orth(1.0, 0.0, 1.0); 
//    P2 = clipper::Coord_orth(0.0, -1.0, 1.0); 
//    P3 = clipper::Coord_orth(0.0, 0.0, 0.0); 
//    P4 = clipper::Coord_orth(-1.0, -1.0, 1.0); 
//    P4 = clipper::Coord_orth(1.0, 1.0, 1.0); 

   clipper::Coord_orth a = P2 - P1; 
   clipper::Coord_orth b = P3 - P2; 
   clipper::Coord_orth c = P4 - P3;

   // b*b * [ a.(bxc)/b ]
   double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
      sqrt( b.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
      + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

   double theta = clipper::Util::rad2d(atan2(E,G));
//    double clipper_theta = 
//       clipper::Util::rad2d(clipper::Coord_orth::torsion(P1, P2, P3, P4));

   if ( clipper::Util::isnan(theta) ) {
      std::cout << "WARNING: observed torsion theta is a NAN!" << std::endl;
   } 

   if (theta < 0.0) theta += 360.0; 

   //if (torsion_restraint.periodicity == 1) {
      // }
   // use period here
   // 
//    double diff = theta - torsion_restraint.target_value;   

   double diff = 99999.9; 
   double tdiff; 
   double trial_target; 
   int per = torsion_restraint.periodicity;
   for(int i=0; i<per; i++) { 
      // trial_target = torsion_restraint.target_value + double(i)*360.0/double(per);  ??
      trial_target = torsion_restraint.target_value + double(i)*360.0/double(per); 
      if (trial_target >= 360.0) trial_target -= 360.0; 
      tdiff = theta - trial_target; 
      if (abs(tdiff) < abs(diff)) { 
	 diff = tdiff;
      }
   }
//    std::cout << "DEBUG:: atom index: " << torsion_restraint.atom_index_1
// 	     << " target " << torsion_restraint.target_value
// 	     << ", theta: " << theta << " \tdiff: " << diff << std::endl;

   if (diff >= 99999.0) { 
      std::cout << "Error in periodicity (" << per << ") check" << std::endl;
      std::cout << "target_value: " << torsion_restraint.target_value
		<< ", theta: " << theta << std::endl;
   }
   if (diff < -180.0) { 
      diff += 360.; 
   } else { 
      if (diff > 180.0) { 
	 diff -= 360.0; 
      }
   }

//    std::cout << "distortion score torsion "
// 	     << diff*diff/(torsion_restraint.sigma * torsion_restraint.sigma) << " ";

//      cout << "distortion_torsion theta (calc): " << theta 
//   	<< " periodicity " << torsion_restraint.periodicity
//   	<< " target "      << torsion_restraint.target_value
//   	<< " diff: " << diff << endl ;

//     std::cout << "in distortion_torsion: sigma = " << torsion_restraint.sigma
//  	     << ", weight=" << pow(torsion_restraint.sigma,-2.0)
//  	     << " and diff is " << diff << std::endl;
   
      
   return diff*diff/(torsion_restraint.sigma * torsion_restraint.sigma);

}

double
coot::distortion_score_plane(const coot::simple_restraint &plane_restraint,
			     const gsl_vector *v) {

   coot::plane_distortion_info_t info =
      distortion_score_plane_internal(plane_restraint, v);

   return info.distortion_score;

}

double 
coot::distortion_score_chiral_volume(const coot::simple_restraint &chiral_restraint,
				     const gsl_vector *v) { 

   int idx = 3*(chiral_restraint.atom_index_centre);
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   idx = 3*(chiral_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));

   clipper::Coord_orth a = a1 - centre;
   clipper::Coord_orth b = a2 - centre;
   clipper::Coord_orth c = a3 - centre;

   double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
   // double volume_sign = chiral_restraint.chiral_volume_sign;
   double distortion = cv - chiral_restraint.target_chiral_volume;

   if (0) {
      std::cout << "atom indices: "
		<< chiral_restraint.atom_index_centre << " "
		<< chiral_restraint.atom_index_1 << " " 
		<< chiral_restraint.atom_index_2 << " " 
		<< chiral_restraint.atom_index_3 << " ";
      std::cout << "DEBUG:: (distortion) chiral volume target "
		<< chiral_restraint.target_chiral_volume
		<< " chiral actual " << cv << " diff: " << distortion;
   }
 	     

   distortion *= distortion;
   distortion /= chiral_restraint.sigma * chiral_restraint.sigma;

   if (0)
      std::cout << " score: " << distortion << "\n";

   return distortion;
}

double 
coot::distortion_score_rama(const coot::simple_restraint &rama_restraint,
			    const gsl_vector *v,
			    const LogRamachandran &lograma) {

   double distortion = 0;
   // First calculate the torsions:
   // theta = arctan(E/G); 
   // where E = a.(bxc) and G = -a.c + (a.b)(b.c)

   int idx; 

   idx = 3*(rama_restraint.atom_index_1);
   clipper::Coord_orth P1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_2); 
   clipper::Coord_orth P2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_3); 
   clipper::Coord_orth P3(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_4); 
   clipper::Coord_orth P4(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_5); 
   clipper::Coord_orth P5(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));

//    P1 = clipper::Coord_orth(1.0, 0.0, 1.0); 
//    P2 = clipper::Coord_orth(0.0, -1.0, 1.0); 
//    P3 = clipper::Coord_orth(0.0, 0.0, 0.0); 
//    P4 = clipper::Coord_orth(-1.0, -1.0, 1.0); 
//    P4 = clipper::Coord_orth(1.0, 1.0, 1.0); 

   clipper::Coord_orth a = P2 - P1; 
   clipper::Coord_orth b = P3 - P2; 
   clipper::Coord_orth c = P4 - P3;
   clipper::Coord_orth d = P5 - P4;

   // Old (6 atom) wrong:
   // TRANS    psi      1 N (P1)    1 CA (P2)     1 C  (P3)    2 N (P4)
   // TRANS    phi      1 C (P3)    2  N (P4)     2 CA (P5)    2 C (P6)
   //
   // New assignements:
   // TRANS    phi    (1st C) (2nd N ) (2nd CA) (2nd C) 
   // TRANS    psi    (2nd N) (2nd CA) (2nd C ) (3nd N)
   // 
   // So Rama_atoms in this order:
   //   0       1        2      3         4
   //  P1      P2       P3     P4        P5
   // (1st C) (2nd N) (2nd CA) (2nd C) (3rd N)

   // ---------- phi ------------------
   // b*b * [ a.(bxc)/b ]
   double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
      sqrt( b.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
      + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

   double phi = clipper::Util::rad2d(atan2(E,G));
   if (phi < 180.0)
      phi += 360.0;
   if (phi > 180.0)
      phi -= 360.0;

   // ---------- psi ------------------
   // b*b * [ a.(bxc)/b ]
   double H = clipper::Coord_orth::dot(b, clipper::Coord_orth::cross(c,d)) *
      sqrt( c.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double I = - clipper::Coord_orth::dot(b,d)*c.lengthsq()
      + clipper::Coord_orth::dot(b,c)*clipper::Coord_orth::dot(c,d);

   double psi = clipper::Util::rad2d(atan2(H,I));
   if (psi < 180.0)
      psi += 360.0;
   if (psi > 180.0)
      psi -= 360.0;

   double lr = lograma.interp(clipper::Util::d2rad(phi), clipper::Util::d2rad(psi));
   double R = 10.0 * lr;
   // std::cout << "rama distortion for " << phi << " " << psi << " is " << R << std::endl;

   if ( clipper::Util::isnan(phi) ) {
      std::cout << "WARNING: observed torsion phi is a NAN!" << std::endl;
   } 
   if ( clipper::Util::isnan(psi) ) {
      std::cout << "WARNING: observed torsion phi is a NAN!" << std::endl;
   } 

   return R;
}


// Before we ask this question, we need to match residues which chiral
// restraints.
//
// i.e we need a function that takes a CMMDBManager * and returns
// pairs of simple_restraint, CResidue, that is done in
//  molecule_class_info-other.cc
// 
std::vector<std::pair<short int, coot::atom_spec_t> >
coot::is_bad_chiral_atom_p(const coot::dict_chiral_restraint_t &chiral_restraint,
			   CResidue *res) {

   short int ibad=0;
   coot::atom_spec_t chiral_atom;
   std::vector<std::pair<short int, coot::atom_spec_t> > v;

   int i_no_res_atoms = res->GetNumberOfAtoms();
   
   for (int iat1=0; iat1<i_no_res_atoms; iat1++) {
      std::string pdb_atom_name1(res->GetAtom(iat1)->name);
      if (pdb_atom_name1 == chiral_restraint.atom_id_1_4c()) {
	 
	 for (int iat2=0; iat2<i_no_res_atoms; iat2++) {
	    std::string pdb_atom_name2(res->GetAtom(iat2)->name);
	    if (pdb_atom_name2 == chiral_restraint.atom_id_2_4c()) {
	       
	       for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
		  std::string pdb_atom_name3(res->GetAtom(iat3)->name);
		  if (pdb_atom_name3 == chiral_restraint.atom_id_3_4c()) {
		     
		     for (int iatc=0; iatc<i_no_res_atoms; iatc++) {
			std::string pdb_atom_namec(res->GetAtom(iatc)->name);
			if (pdb_atom_namec == chiral_restraint.atom_id_c_4c()) {

			   // Now, do they have corresponding
			   // altconfs?  The bonded atoms must have
			   // the same altconf as each other (or blank).
			   
			   // 
			   // e.g. for Chiral atom CB in THR: CA   CB,A OG1,A CG2,A (1)
			   //                                 CA   CB   OG1,A CG2,A (2)
			   //                 but not:        CA   CB   OG1,A CG2,B (3)

			   std::string chiral_alt_conf(res->GetAtom(iatc)->altLoc);
			   std::string altLoc1(res->GetAtom(iat1)->altLoc);
			   std::string altLoc2(res->GetAtom(iat2)->altLoc);
			   std::string altLoc3(res->GetAtom(iat3)->altLoc);
			   short int matching_altlocs = 0;

			   // These should catch most cases, 
			   // but there may be some that they don't.

			   if ( chiral_alt_conf != "" &&  
			       (altLoc1 == "" || altLoc1 == chiral_alt_conf) &&
			       (altLoc2 == "" || altLoc2 == chiral_alt_conf) &&
			       (altLoc3 == "" || altLoc3 == chiral_alt_conf) ) { 
			      matching_altlocs = 1;
			   } 
			   if ( chiral_alt_conf == "" &&  
			       (altLoc1 == "" ) &&
			       (altLoc2 == "" ) &&
			       (altLoc3 == "" ) ) { 
			      matching_altlocs = 1;
			   }

			   // We need to check like this to get rid of case (3):
			   if ( chiral_alt_conf == "" &&  
			       (altLoc1 == altLoc2 || altLoc1 == "" || altLoc2 == "") &&
			       (altLoc1 == altLoc3 || altLoc1 == "" || altLoc3 == "") &&
			       (altLoc2 == altLoc3 || altLoc2 == "" || altLoc3 == "") ) { 
			      matching_altlocs = 1;
			   }

			   if (matching_altlocs) { 

			      clipper::Coord_orth centre(res->GetAtom(iatc)->x,
							 res->GetAtom(iatc)->y,
							 res->GetAtom(iatc)->z);
			      clipper::Coord_orth a1(res->GetAtom(iat1)->x,
						     res->GetAtom(iat1)->y,
						     res->GetAtom(iat1)->z);
			      clipper::Coord_orth a2(res->GetAtom(iat2)->x,
						     res->GetAtom(iat2)->y,
						     res->GetAtom(iat2)->z);
			      clipper::Coord_orth a3(res->GetAtom(iat3)->x,
						     res->GetAtom(iat3)->y,
						     res->GetAtom(iat3)->z);
			   
			      clipper::Coord_orth a = a1 - centre;
			      clipper::Coord_orth b = a2 - centre;
			      clipper::Coord_orth c = a3 - centre;
			      double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));

			      chiral_atom = coot::atom_spec_t(res->GetChainID(),
							      res->GetSeqNum(),
							      res->GetInsCode(),
							      res->GetAtom(iatc)->name,
							      res->GetAtom(iatc)->altLoc);

			      if (cv*chiral_restraint.volume_sign < 0) {
// 				 std::cout << "DEBUG:: " << res->name << " "
// 					   << res->GetSeqNum() << " " << pdb_atom_namec
// 					   << " has bad chiral volume\n";
				 ibad = 1;
			      } else { 
// 				 std::cout << "DEBUG:: " << res->name << " "
// 					   << res->GetSeqNum() << " " << pdb_atom_namec
// 					   << " has good chiral volume\n";
				 ibad = 0;
			      }
			      v.push_back(std::pair<short int, coot::atom_spec_t> (ibad, chiral_atom));
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   
   // return std::pair<short int, coot::atom_spec_t> (ibad, chiral_atom);
   return v;
}


// Return a list of bad chiral volumes for this molecule:
// 
// Return also a flag for the status of this test, were there any
// residues for which we didn't find restraints?  The flag is the
// number of residue names in first part of the returned pair.
// 0 is good.
// 
std::pair<std::vector<std::string> , std::vector <coot::atom_spec_t> >
coot::bad_chiral_volumes(CMMDBManager *mol, protein_geometry *geom_p,
			 int cif_dictionary_read_number) {

   std::vector <coot::atom_spec_t> v;
   int restraints_status = 1;
   std::vector<std::string> unknown_types_vec; 

#ifdef HAVE_GSL
   
   int n_models = mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) { 
      
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      if (nchains <= 0) { 
	 std::cout << "bad nchains in molecule " << nchains
		   << std::endl;
      } else { 
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    if (chain_p == NULL) {  
	       // This should not be necessary. It seem to be a
	       // result of mmdb corruption elsewhere - possibly
	       // DeleteChain in update_molecule_to().
	       std::cout << "NULL chain in ... " << std::endl;
	    } else { 
	       int nres = chain_p->GetNumberOfResidues();
	       PCResidue residue_p;
	       // CAtom *at;
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  if (n_atoms > 3) {
		     std::string residue_type(residue_p->name);
		     if (residue_type == "UNK")
			residue_type = "ALA";
		     if (! geom_p->have_dictionary_for_residue_type(residue_type,
								  cif_dictionary_read_number)) {
			std::cout << "WARNING::! Failed to find restraint for residue type "
				  << residue_type << std::endl;
			restraints_status = 0;
			// only add the residue_type if it is not already in the vector
			short int irfound = 0;
			for (unsigned int ir=0; ir<unknown_types_vec.size(); ir++)
			   if (unknown_types_vec[ir] == residue_type) {
			      irfound = 1;
			      break;
			   }
			if (irfound == 0) 
			   unknown_types_vec.push_back(residue_type);
		     } else { 
			std::vector<coot::dict_chiral_restraint_t> chiral_restraints = 
			   geom_p->get_monomer_chiral_volumes(std::string(residue_p->name));
			coot::dict_chiral_restraint_t chiral_restraint;
			for (unsigned int irestr=0; irestr<chiral_restraints.size(); irestr++) { 
			   chiral_restraint = chiral_restraints[irestr];
			   if (! chiral_restraint.is_a_both_restraint()) {
			      std::vector<std::pair<short int, coot::atom_spec_t> > c = 
				 coot::is_bad_chiral_atom_p(chiral_restraint, residue_p);
			      for (unsigned int ibad=0; ibad<c.size(); ibad++) { 
				 if (c[ibad].first) {
				    std::cout << "INFO:: found bad chiral atom: " 
					      << chain_p->GetChainID() << " " 
					      << residue_p->GetSeqNum() << " "
					      << residue_p->GetInsCode() << " "
					      << chiral_restraint.atom_id_c_4c() << " "
					      << c[ibad].second.alt_conf << std::endl;
				    
				    v.push_back(coot::atom_spec_t(chain_p->GetChainID(),
								  residue_p->GetSeqNum(),
								  residue_p->GetInsCode(),
								  chiral_restraint.atom_id_c_4c(),
								  c[ibad].second.alt_conf));
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
   }
#endif //  HAVE_GSL
   std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t> > pair(unknown_types_vec, v);
   return pair;
}


			     
void
coot::fix_chiral_atom_maybe (const simple_restraint &chiral_restraint,
			     gsl_vector *v) {
   int idx = 3*(chiral_restraint.atom_index_centre);
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   idx = 3*(chiral_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   clipper::Coord_orth a = a1 - centre;
   clipper::Coord_orth b = a2 - centre;
   clipper::Coord_orth c = a3 - centre;

   double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
   double volume_sign = chiral_restraint.chiral_volume_sign;

   std::cout << "DEBUG:::::::: Fix chiral maybe :::::: "
	     << volume_sign*cv << std::endl;

   if (volume_sign*cv < 0.0) {
      // we have a mismatch between restraint and actual, needs fixing
      std::cout << "Atom index " << chiral_restraint.atom_index_centre
		<< " is undergoing chiral centre inversion\n";
      fix_chiral_atom_internal(chiral_restraint, v);
   }
}

void
coot::fix_chiral_atom_internal(const simple_restraint &chiral_restraint,
			       gsl_vector *v) {

   // we need to find the chiral plane (i.e. the plane of the
   // non-centre atoms of the chiral restraint).
   //
   // We need to find the distance of the chiral atom from the chiral
   // plane (the vector is normal to the plane - and the plane
   // equation defines the normal).
   //
   // We will then move the atom to the other side of the plane.
   //

   //Recall d = Ax+By+Cz-D, where d is the distance from the
   //plane of the chiral atom.
   //

   // So we have 3 points, a b and c.  The cross-product b-a x c-a
   // gives us the normal the the plane.

   int idx = 3*(chiral_restraint.atom_index_centre);
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   idx = 3*(chiral_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));

   clipper::Coord_orth a1_2 = a1 - a2;
   clipper::Coord_orth a2_3 = a2 - a3;
   clipper::Coord_orth a3_1 = a3 - a1;

//     clipper::Coord_orth normal =
//        clipper::Coord_orth(b_a.y()*c_a.z() - b_a.z()*c_a.y(),
//  			  b_a.z()*c_a.x() - b_a.x()*c_a.z(),
//  			  b_a.x()*c_a.y() - b_a.y()*c_a.x());

   clipper::Coord_orth scaled_normal =
      clipper::Coord_orth(a1.y()*a2_3.z() + a2.y()*a3_1.z() + a3.y()*a1_2.z(),
			  a1.z()*a2_3.x() + a2.z()*a3_1.x() + a3.z()*a1_2.x(),
			  a1.x()*a2_3.y() + a2.x()*a3_1.y() + a3.x()*a1_2.y());
   double scaled_D =   (
			a1.x()*(a2.y()*a3.z()-a3.y()*a2.z()) + 
		        a2.x()*(a3.y()*a1.z()-a1.y()*a3.z()) +
			a3.x()*(a1.y()*a2.z()-a2.y()*a1.z()));

   // normalize:
   double len_s_n = sqrt(scaled_normal.lengthsq());
   double one_over_len_sn = 1.0/len_s_n;

   clipper::Coord_orth normal = one_over_len_sn * scaled_normal;
   double D = one_over_len_sn * scaled_D;

   std::cout << "normal now: " << normal.format() << "D: " << D << "\n";
   
   //d = Ax+By+Cz-D
   double d = normal.x()*centre.x() + normal.y()*centre.y() + normal.z()*centre.z() - D;

   std::cout << "d is " << d << " for atom index "
	     << chiral_restraint.atom_index_centre << "\n";

   // Cool, now lets move the centre atom to the other side of the plane
   // by 0.5A

   double factor;
   if (d<0)
      factor = d - 0.5;
   else
      factor = d + 0.5;
   clipper::Coord_orth shift = - factor*normal;

   std::cout  << "DEBUG::  CHIRAL: shifting atom index "
	      << chiral_restraint.atom_index_centre << " by "
	      << shift.format() << "\n";

   idx = 3*(chiral_restraint.atom_index_centre);
   gsl_vector_set(v, idx,   gsl_vector_get(v, idx  ) + shift.x());
   gsl_vector_set(v, idx+1, gsl_vector_get(v, idx+1) + shift.y());
   gsl_vector_set(v, idx+2, gsl_vector_get(v, idx+2) + shift.z());
}


double
coot::distortion_score_non_bonded_contact(const coot::simple_restraint &bond_restraint,
					  const gsl_vector *v) {

   int idx = 3*(bond_restraint.atom_index_1 - 0); 
   clipper::Coord_orth a1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(bond_restraint.atom_index_2 - 0); 
   clipper::Coord_orth a2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));

   double weight = 1/(bond_restraint.sigma * bond_restraint.sigma);

   double dist = clipper::Coord_orth::length(a1,a2);
   double bit; 

   double r = 0.0;
   if (dist < bond_restraint.target_value) {
      bit = dist - bond_restraint.target_value;
      r = weight * bit * bit;
   }

   return r;
}

coot::plane_distortion_info_t
coot::distortion_score_plane_internal(const coot::simple_restraint &plane_restraint,
				      const gsl_vector *v) {

   coot::plane_distortion_info_t info;
   double sum_devi = 0; // return this (after weighting)

   // Recall the equation of a plane: ax + by + cz + d = 0:
   // 
   // First find the centres of the sets of atoms: x_cen, y_cen,
   // z_cen.  We move the plane down so that it crosses the origin and
   // there for d=0 (we'll add it back later).  This makes the problem
   // into 3 equations, 3 unknowns, an eigenvalue problem, with the
   // smallest eigenvalue corresponding to the best fit plane.
   //
   // Google for: least squares plane:
   // http://www.infogoaround.org/JBook/LSQ_Plane.html
   //
   //
   //
   // So, first get the centres:
   double sum_x = 0, sum_y = 0, sum_z = 0;
   int idx;
   int n_atoms = plane_restraint.atom_index.size();
   double dn_atoms = double(n_atoms); 

   if (n_atoms > 0) { 
      for (int i=0; i<n_atoms; i++) {
	 idx = 3*(plane_restraint.atom_index[i]);
// 	 std::cout << "atom_index vs n_atom " << plane_restraint.atom_index[i]
// 		   << " " << n_atoms << std::endl;
	 if (plane_restraint.atom_index[i] < 0) {
	    std::cout << "trapped bad plane restraint! " << plane_restraint.atom_index[i]
		      << std::endl;
	    // return info;
	 } else {
	    sum_x += gsl_vector_get(v,idx);
	    sum_y += gsl_vector_get(v,idx+1);
	    sum_z += gsl_vector_get(v,idx+2);
	 }
      }
      double x_cen = sum_x/dn_atoms;
      double y_cen = sum_y/dn_atoms;
      double z_cen = sum_z/dn_atoms;

      clipper::Matrix<double> mat(3,3);
      
      for (int i=0; i<n_atoms; i++) {
	 idx = 3*(plane_restraint.atom_index[i]);
	 if (plane_restraint.atom_index[i] < 0) {
	 } else { 
	    mat(0,0) += (gsl_vector_get(v,idx  ) - x_cen) * (gsl_vector_get(v,idx  ) - x_cen);
	    mat(1,1) += (gsl_vector_get(v,idx+1) - y_cen) * (gsl_vector_get(v,idx+1) - y_cen);
	    mat(2,2) += (gsl_vector_get(v,idx+2) - z_cen) * (gsl_vector_get(v,idx+2) - z_cen);
	    mat(0,1) += (gsl_vector_get(v,idx  ) - x_cen) * (gsl_vector_get(v,idx+1) - y_cen);
	    mat(0,2) += (gsl_vector_get(v,idx  ) - x_cen) * (gsl_vector_get(v,idx+2) - z_cen);
	    mat(1,2) += (gsl_vector_get(v,idx+1) - y_cen) * (gsl_vector_get(v,idx+2) - z_cen);
	 }
      }
      mat(1,0) = mat(0,1);
      mat(2,0) = mat(0,2);
      mat(2,1) = mat(1,2);
      std::vector<double> eigens = mat.eigen(true);

//       std::cout << "we get eigen values: "
// 		<< eigens[0] << "  "
// 		<< eigens[1] << "  "
// 		<< eigens[2] << std::endl;

      // Let's now extract the values of a,b,c normalize them
      std::vector<double> abcd;
      abcd.push_back(mat(0,0));
      abcd.push_back(mat(1,0));
      abcd.push_back(mat(2,0));
      double sqsum = 1e-20;

      for (int i=0; i<3; i++)
	 sqsum += abcd[i] * abcd[i];
      for (int i=0; i<3; i++)
	 abcd[i] /= sqsum;


      //set d, recall di = Axi+Byi+Czi-D, so xi = x_cen, yi = y_cen, zi = z_cen:
      abcd.push_back( abcd[0]*x_cen + abcd[1]*y_cen + abcd[2]*z_cen );
      info.abcd = abcd;

      double val;
      for (int i=0; i<n_atoms; i++) {
	 idx = 3*(plane_restraint.atom_index[i]);
	 if (idx < 0) {
	 } else { 
	    val = 
	       abcd[0]*gsl_vector_get(v,idx  ) +
	       abcd[1]*gsl_vector_get(v,idx+1) +
	       abcd[2]*gsl_vector_get(v,idx+2) -
	       abcd[3];
	    sum_devi += val * val;
	 }
      }
      // std::cout << "sum_devi is : " << sum_devi << std::endl;

   }
//    std::cout << "distortion_score_plane returning "
// 	     <<  plane_restraint.sigma << "^-2 * "
// 	     << sum_devi << " = "
//  	     << pow(plane_restraint.sigma,-2.0) * sum_devi << std::endl;
   info.distortion_score = sum_devi / (plane_restraint.sigma * plane_restraint.sigma);
   
   return info;
}


std::vector<coot::refinement_lights_info_t>
coot::restraints_container_t::chi_squareds(std::string title, const gsl_vector *v) const {

   std::vector<coot::refinement_lights_info_t> lights_vec;
   int n_bond_restraints = 0; 
   int n_angle_restraints = 0; 
   int n_torsion_restraints = 0; 
   int n_plane_restraints = 0; 
   int n_non_bonded_restraints = 0;
   int n_chiral_volumes = 0;
   int n_rama_restraints = 0;

   double bond_distortion = 0; 
   double angle_distortion = 0; 
   double torsion_distortion = 0; 
   double plane_distortion = 0; 
   double non_bonded_distortion = 0;
   double chiral_vol_distortion = 0;
   double rama_distortion = 0;

   for (int i=0; i<size(); i++) {
      if (restraints_usage_flag & coot::BONDS_MASK) { // 1: bonds
	 if ( (*this)[i].restraint_type == coot::BOND_RESTRAINT) {
	    n_bond_restraints++;
	    bond_distortion += coot::distortion_score_bond((*this)[i], v);
	 }
      }

      if (restraints_usage_flag & coot::ANGLES_MASK) { // 2: angles
	 if ( (*this)[i].restraint_type == coot::ANGLE_RESTRAINT) {
	    n_angle_restraints++;
	    angle_distortion += coot::distortion_score_angle((*this)[i], v);
	 }
      }

//       std::cout << "DEBUG:: ---- restraints_usage_flag & coot::TORSIONS_MASK: "
// 		<< (restraints_usage_flag & coot::TORSIONS_MASK) << std::endl;
      if (restraints_usage_flag & coot::TORSIONS_MASK) { // 4: torsions
	 if ( (*this)[i].restraint_type == coot::TORSION_RESTRAINT) {
	    n_torsion_restraints++;
	    torsion_distortion += coot::distortion_score_torsion((*this)[i], v); 
	 }
      }

      if (restraints_usage_flag & coot::PLANES_MASK) { // 8: planes
	 if ( (*this)[i].restraint_type == coot::PLANE_RESTRAINT) {
	    n_plane_restraints++;
	    plane_distortion += coot::distortion_score_plane((*this)[i], v); 
	 }
      }

      if (restraints_usage_flag & coot::NON_BONDED_MASK) { 
	 if ( (*this)[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) { 
	    n_non_bonded_restraints++;
	    non_bonded_distortion += coot::distortion_score_non_bonded_contact((*this)[i], v);
	 }
      }

      if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) { 
  	 if ( (*this)[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) { 
  	    n_chiral_volumes++;
  	    chiral_vol_distortion += coot::distortion_score_chiral_volume((*this)[i], v);
  	 }
      }

      if (restraints_usage_flag & coot::RAMA_PLOT_MASK) {
  	 if ( (*this)[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) { 
  	    n_rama_restraints++;
  	    rama_distortion += coot::distortion_score_rama((*this)[i], v, lograma);
  	 }
      }
   }

   std::string r = "";

   r += title;
   r += "\n";
   std::cout << "    " << title << std::endl;
   if (n_bond_restraints == 0) {
      std::cout << "bonds:      N/A " << std::endl;
   } else {
      std::cout << "bonds:      " << bond_distortion/double(n_bond_restraints)
		<< std::endl;
      double bd = bond_distortion/double(n_bond_restraints);
      r += "   bonds:  ";
      r += coot::util::float_to_string_using_dec_pl(bd, 3);
      r += "\n";
      std::string s = "Bonds:  ";
      s += coot::util::float_to_string_using_dec_pl(bd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Bonds", s, bd));
   } 
   if (n_angle_restraints == 0) {
      std::cout << "angles:     N/A " << std::endl;
   } else {
      double ad = angle_distortion/double(n_angle_restraints);
      std::cout << "angles:     " << angle_distortion/double(n_angle_restraints)
		<< std::endl;
      r += "   angles: ";
      r += coot::util::float_to_string_using_dec_pl(angle_distortion/double(n_angle_restraints), 3);
      r += "\n";
      std::string s = "Angles: ";
      s += coot::util::float_to_string_using_dec_pl(ad, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Angles", s, ad));

   } 
   if (n_torsion_restraints == 0) {
      std::cout << "torsions:   N/A " << std::endl;
   } else {
      double td = torsion_distortion/double(n_torsion_restraints);
      std::cout << "torsions:   " << td << std::endl;
      r += "   torsions: ";
      r += coot::util::float_to_string_using_dec_pl(td, 3);
      r += "\n";
      std::string s = "Torsions: ";
      s += coot::util::float_to_string_using_dec_pl(td, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Torsions", s, td));
   } 
   if (n_plane_restraints == 0) {
      std::cout << "planes:     N/A " << std::endl;
   } else {
      double pd = plane_distortion/double(n_plane_restraints);
      std::cout << "planes:     " << pd << std::endl;
      r += "   planes: ";
      r += coot::util::float_to_string_using_dec_pl(pd, 3);
      r += "\n";
      std::string s = "Planes: ";
      s += coot::util::float_to_string_using_dec_pl(pd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Planes", s, pd));
   }
   if (n_non_bonded_restraints == 0) {
      std::cout << "non-bonded: N/A " << std::endl;
   } else {
      double nbd = non_bonded_distortion/double(n_non_bonded_restraints);
      std::cout << "non-bonded: " << nbd
		<< std::endl;
      r += "   non-bonded: ";
      r += coot::util::float_to_string_using_dec_pl(nbd, 3);
      r += "\n";
      std::string s = "Non-bonded: ";
      s += coot::util::float_to_string_using_dec_pl(nbd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Non-bonded", s, nbd));
   }
   if (n_chiral_volumes == 0) { 
      std::cout << "chiral vol: N/A " << std::endl;
   } else {
      double cd = chiral_vol_distortion/double(n_chiral_volumes);
      std::cout << "chiral vol: " << cd << std::endl;
      r += "   chirals: ";
      r += coot::util::float_to_string_using_dec_pl(cd, 3);
      std::string s = "Chirals: ";
      s += coot::util::float_to_string_using_dec_pl(cd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Chirals", s, cd));
   }
   if (n_rama_restraints == 0) { 
      std::cout << "rama plot:  N/A " << std::endl;
   } else {
      double rd = rama_distortion/double(n_rama_restraints);
         std::cout << "rama plot:  " << rd << std::endl;
      r += "   rama plot: ";
      r += coot::util::float_to_string_using_dec_pl(rd, 3);
      std::string s = "Rama Plot: ";
      s += coot::util::float_to_string_using_dec_pl(rd, 3);
      lights_vec.push_back(coot::refinement_lights_info_t("Rama", s, rd)); // correct name?
   }
   return lights_vec;
} 

void 
coot::numerical_gradients(gsl_vector *v, 
			  void *params, 
			  gsl_vector *df) {

//    clipper::Coord_orth a(0,0,2); 
//    cout << "length test: " << a.lengthsq() << endl; 

   // double S = coot::distortion_score(v, params); 
   double tmp; 
   double new_S_plus; 
   double new_S_minu; 
   double val;
   double micro_step = 0.0001;  // the difference between the gradients
			        // seems not to depend on the
			        // micro_step size (0.0001 vs 0.001)

   std::cout << "analytical_gradients" << std::endl; 
   for (unsigned int i=0; i<df->size; i++) {
      tmp = gsl_vector_get(df, i); 
      cout << tmp << "  " << i << endl; 
   }
   
   std::cout << "numerical_gradients" << std::endl; 
   for (unsigned int i=0; i<v->size; i++) { 
      
      tmp = gsl_vector_get(v, i); 
      gsl_vector_set(v, i, tmp+micro_step); 
      new_S_plus = coot::distortion_score(v, params); 
      gsl_vector_set(v, i, tmp-micro_step); 
      new_S_minu = coot::distortion_score(v, params);
      // new_S_minu = 2*tmp - new_S_plus; 

      // now put v[i] back to tmp
      gsl_vector_set(v, i, tmp);
      
      val = (new_S_plus - new_S_minu)/(2*micro_step); 
      cout << val << "  " << i << endl; 
      //

      // overwrite the analytical gradients with numerical ones:
      // gsl_vector_set(df, i, val);
   } 
} // Note to self: try 5 atoms and doctor the .rst file if necessary.
  // Comment out the bond gradients.
  // Try getting a perfect structure (from refmac idealise) and
  // distorting an atom by adding 0.5 to one of its coordinates.  We
  // should be able to trace the gradients that way.


void coot::my_df(const gsl_vector *v, 
		 void *params, 
		 gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   int n_var = restraints->n_variables();

   // first initialize the derivative vector:
   for (int i=0; i<n_var; i++) {
      gsl_vector_set(df,i,0);
   }

   coot::my_df_bonds     (v, params, df); 
   coot::my_df_angles    (v, params, df);
   coot::my_df_torsions  (v, params, df);
   coot::my_df_rama      (v, params, df);
   coot::my_df_planes    (v, params, df);
   coot::my_df_non_bonded(v, params, df);
   coot::my_df_chiral_vol(v, params, df);
   
   // my_df_start_pos(v,params,df); // If you add these back, don't
                                    // forget to create the corresponding 
                                    // distortion_score code.
   
   if (restraints->include_map_terms()) {
      // std::cout << "Using map terms " << std::endl;
      coot::my_df_electron_density((gsl_vector *)v,params,df);
   } 

   if (restraints->do_numerical_gradients_status())
      coot::numerical_gradients((gsl_vector *)v,params,df); 

}
   
/* The gradients of f, df = (df/dx(k), df/dy(k) .. df/dx(l) .. ). */
void coot::my_df_bonds (const gsl_vector *v, 
			void *params,
			gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 
   
   // the length of gsl_vector should be equal to n_var: 
   // 
   // int n_var = restraints->n_variables();
   //    float derivative_value; 
   int idx; 
   int n_bond_restr = 0; // debugging counter
   
   //     for (int i=0; i<n_var; i++) { 
   //       gsl_vector_set(df, i, derivative_value); 
   //     } 
   
    // Now run over the bonds
    // and add the contribution from this bond/restraint to 
    // dS/dx_k dS/dy_k dS/dz_k dS/dx_l dS/dy_l dS/dz_l for each bond
    // 

   if (restraints->restraints_usage_flag & coot::BONDS_MASK) {

      double target_val;
      double b_i;
      double weight;
      double x_k_contrib;
      double y_k_contrib;
      double z_k_contrib;
      
      double x_l_contrib;
      double y_l_contrib;
      double z_l_contrib;

      for (int i=0; i<restraints->size(); i++) {
       
	 if ( (*restraints)[i].restraint_type == coot::BOND_RESTRAINT) { 

// 	    std::cout << "DEBUG bond restraint fixed flags: "
// 		      << (*restraints)[i].fixed_atom_flags[0] << " "
// 		      << (*restraints)[i].fixed_atom_flags[1] << " "
// 		      << restraints->get_atom((*restraints)[i].atom_index_1)->GetSeqNum() << " "
// 		      << restraints->get_atom((*restraints)[i].atom_index_1)->name
// 		      << " to " 
// 		      << restraints->get_atom((*restraints)[i].atom_index_2)->GetSeqNum() << " "
// 		      << restraints->get_atom((*restraints)[i].atom_index_2)->name
// 		      << std::endl;
	    
// 	    int n_fixed=0;
// 	    if  ((*restraints)[i].fixed_atom_flags[0])
// 	       n_fixed++;
// 	    if  ((*restraints)[i].fixed_atom_flags[1])
// 	       n_fixed++;
	    
	    n_bond_restr++; 

	    target_val = (*restraints)[i].target_value;

	    // what is the index of x_k?
	    idx = 3*((*restraints)[i].atom_index_1); 
	    clipper::Coord_orth a1(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_2); 
	    clipper::Coord_orth a2(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));

	    // what is b_i?
	    b_i = clipper::Coord_orth::length(a1,a2); 
	    b_i = b_i > 0.1 ? b_i : 0.1;  // Garib's stabilization

	    weight = 1/((*restraints)[i].sigma * (*restraints)[i].sigma);
	    
	    // weight = 1.0;
	    // weight = pow(0.021, -2.0);
	    // std::cout << "df weight is " << weight << std::endl;

	    double constant_part = 2.0*weight*(b_i - target_val)/b_i;
	    
	    x_k_contrib = constant_part*(a1.x()-a2.x());
	    y_k_contrib = constant_part*(a1.y()-a2.y());
	    z_k_contrib = constant_part*(a1.z()-a2.z());
	    
	    x_l_contrib = constant_part*(a2.x()-a1.x());
	    y_l_contrib = constant_part*(a2.y()-a1.y());
	    z_l_contrib = constant_part*(a2.z()-a1.z());

	    if (!(*restraints)[i].fixed_atom_flags[0]) { 
	       idx = 3*((*restraints)[i].atom_index_1 - 0);  
	       // cout << "first  idx is " << idx << endl; 
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_k_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_k_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_k_contrib); 
	    } else {
	       idx = 3*((*restraints)[i].atom_index_1 - 0);  
// 	       std::cout << "Fixed atom[0] "
// 			 << restraints->get_atom((*restraints)[i].atom_index_1)->GetSeqNum() << " " 
// 			 << restraints->get_atom((*restraints)[i].atom_index_1)->name << " " 
// 			 << ", Not adding " << x_k_contrib << " "
// 			 << y_k_contrib << " "
// 			 << z_k_contrib << " to " << gsl_vector_get(df, idx) << " "
// 			 << gsl_vector_get(df, idx+1) << " "
// 			 << gsl_vector_get(df, idx+2) << std::endl;
	    } 

	    if (!(*restraints)[i].fixed_atom_flags[1]) { 
	       idx = 3*((*restraints)[i].atom_index_2 - 0); 
	       // cout << "second idx is " << idx << endl; 
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_l_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_l_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_l_contrib); 
	    } else {
	       idx = 3*((*restraints)[i].atom_index_2 - 0);  
// 	       std::cout << "Fixed atom[1] "
// 			 << restraints->get_atom((*restraints)[i].atom_index_2)->GetSeqNum() << " " 
// 			 << restraints->get_atom((*restraints)[i].atom_index_2)->name << " " 
// 			 << ", Not adding " << x_k_contrib << " "
// 			 << y_k_contrib << " "
// 			 << z_k_contrib << " to "
// 			 << gsl_vector_get(df, idx) << " "
// 			 << gsl_vector_get(df, idx+1) << " "
// 			 << gsl_vector_get(df, idx+2) << std::endl;
	    } 
	 }
      }
   }
}

void
coot::my_df_non_bonded(const  gsl_vector *v, 
			void *params, 
			gsl_vector *df) {
   
   
   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 
   
   // the length of gsl_vector should be equal to n_var: 
   // 
   // int n_var = restraints->n_variables();
   // float derivative_value; 
   int idx; 
   int n_non_bonded_restr = 0; // debugging counter
   

   if (restraints->restraints_usage_flag & coot::NON_BONDED_MASK) { 

      double target_val;
      double b_i;
      double weight;
      double x_k_contrib;
      double y_k_contrib;
      double z_k_contrib;
      
      double x_l_contrib;
      double y_l_contrib;
      double z_l_contrib;

      for (int i=0; i<restraints->size(); i++) { 
      
	 if ( (*restraints)[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) { 

	    n_non_bonded_restr++; 
	    
	    target_val = (*restraints)[i].target_value;
	    
	    // what is the index of x_k?
	    idx = 3*( (*restraints)[i].atom_index_1 );

	    clipper::Coord_orth a1(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_2); 
	    clipper::Coord_orth a2(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));

	    // what is b_i?
	    b_i = clipper::Coord_orth::length(a1,a2);

	    // Just a bit of debugging	    
// 	    if (b_i > 0.1 ) 
// 	       // nothing
// 	       float jj=0;
// 	    else
// 	       std::cout << "Garib stabilization in play!" << std::endl;
	    
	    b_i = b_i > 0.1 ? b_i : 0.1;  // Garib's stabilization

	    weight = 1/( (*restraints)[i].sigma * (*restraints)[i].sigma );

	    if (b_i < (*restraints)[i].target_value) {

	       double constant_part = 2.0*weight*(b_i - target_val)/b_i;
	       x_k_contrib = constant_part*(a1.x()-a2.x());
	       y_k_contrib = constant_part*(a1.y()-a2.y());
	       z_k_contrib = constant_part*(a1.z()-a2.z());

	       x_l_contrib = constant_part*(a2.x()-a1.x());
	       y_l_contrib = constant_part*(a2.y()-a1.y());
	       z_l_contrib = constant_part*(a2.z()-a1.z());

	       if (!(*restraints)[i].fixed_atom_flags[0]) { 
		  idx = 3*((*restraints)[i].atom_index_1 - 0); 
		  // cout << "first  idx is " << idx << endl; 
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_k_contrib); 
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_k_contrib); 
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_k_contrib); 
	       }

	       if (!(*restraints)[i].fixed_atom_flags[1]) { 
		  idx = 3*((*restraints)[i].atom_index_2 - 0); 
		  // cout << "second idx is " << idx << endl; 
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_l_contrib); 
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_l_contrib); 
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_l_contrib); 
	       }
	    }
	 }
      }
   }
}

// Add in the angle gradients
//
void coot::my_df_angles(const gsl_vector *v, 
			void *params, 
			gsl_vector *df) {

   int n_angle_restr = 0; 
   int idx; 

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   if (restraints->restraints_usage_flag & coot::ANGLES_MASK) { // 2: angles

      double a;
      double b;
      double l_over_a_sqd;
      double l_over_b_sqd;
      double l_ab;

      double a_dot_b;
      double cos_theta;
      double theta;
      double prem;

      double x_k_contrib;
      double y_k_contrib;
      double z_k_contrib;

      double x_m_contrib;
      double y_m_contrib;
      double z_m_contrib;

      double term1x;
      double term1y;
      double term1z;

      double term2x;
      double term2y;
      double term2z;

      double x_l_mid_contrib;
      double y_l_mid_contrib;
      double z_l_mid_contrib;

      double weight;
      double ds_dth;
      double w_ds_dth;

      for (int i=0; i<restraints->size(); i++) {
      
	 if ( (*restraints)[i].restraint_type == coot::ANGLE_RESTRAINT) {

	    n_angle_restr++;

	    double target_value = (*restraints)[i].target_value*DEGTORAD;

	    idx = 3*((*restraints)[i].atom_index_1); 
	    clipper::Coord_orth k(gsl_vector_get(v,idx), 
			gsl_vector_get(v,idx+1), 
			gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_2); 
	    clipper::Coord_orth l(gsl_vector_get(v,idx), 
			gsl_vector_get(v,idx+1), 
			gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_3); 
	    clipper::Coord_orth m(gsl_vector_get(v,idx), 
			gsl_vector_get(v,idx+1), 
			gsl_vector_get(v,idx+2));

	    clipper::Coord_orth a_vec = (k - l); 
	    clipper::Coord_orth b_vec = (m - l);  
	    a = sqrt(a_vec.lengthsq());
	    b = sqrt(b_vec.lengthsq()); 

	    a = a > 0.1 ? a : 0.1;  // Garib's stabilization
	    b = b > 0.1 ? b : 0.1;  // 
	    
	    l_over_a_sqd = 1/(a*a);
	    l_over_b_sqd = 1/(b*b);
	    l_ab = 1/(a*b); 

	    // for the end atoms: 
	    // \frac{\partial \theta}{\partial x_k} =
	    //    -\frac{1}{sin\theta} [(x_l-x_k)cos\theta + \frac{x_m-x_l}{ab}]
	 
	    // we need to stabilize $\theta$ too.
	    a_dot_b = clipper::Coord_orth::dot(a_vec,b_vec); 
	    cos_theta = a_dot_b/(a*b); 
	    theta = acos(cos_theta); 
	    // theta = theta > 0.001 ? theta : 0.001;
	    if (theta < 0.001) theta = 0.001; // it was never -ve.

	    prem = -1/sin(theta); 
	 
	    // The end atoms:
	    x_k_contrib = prem*(cos_theta*(l.x()-k.x())*l_over_a_sqd + l_ab*(m.x()-l.x()));
	    y_k_contrib = prem*(cos_theta*(l.y()-k.y())*l_over_a_sqd + l_ab*(m.y()-l.y()));
	    z_k_contrib = prem*(cos_theta*(l.z()-k.z())*l_over_a_sqd + l_ab*(m.z()-l.z()));
	    
	    x_m_contrib = prem*(cos_theta*(l.x()-m.x())*l_over_b_sqd + l_ab*(k.x()-l.x()));
	    y_m_contrib = prem*(cos_theta*(l.y()-m.y())*l_over_b_sqd + l_ab*(k.y()-l.y()));
	    z_m_contrib = prem*(cos_theta*(l.z()-m.z())*l_over_b_sqd + l_ab*(k.z()-l.z()));

	    // For the middle atom, we have more cross terms in 
	    // the derivatives of ab and a_dot_b.
	    // 
	    // I will split it up so that it is easier to read: 
	    // 
	    term1x = (-cos_theta*(l.x()-k.x())*l_over_a_sqd) -cos_theta*(l.x()-m.x())*l_over_b_sqd;
	    term1y = (-cos_theta*(l.y()-k.y())*l_over_a_sqd) -cos_theta*(l.y()-m.y())*l_over_b_sqd;
	    term1z = (-cos_theta*(l.z()-k.z())*l_over_a_sqd) -cos_theta*(l.z()-m.z())*l_over_b_sqd;

	    term2x = (-(k.x()-l.x())-(m.x()-l.x()))*l_ab;
	    term2y = (-(k.y()-l.y())-(m.y()-l.y()))*l_ab; 
	    term2z = (-(k.z()-l.z())-(m.z()-l.z()))*l_ab; 

	    x_l_mid_contrib = prem*(term1x + term2x); 
	    y_l_mid_contrib = prem*(term1y + term2y); 
	    z_l_mid_contrib = prem*(term1z + term2z);

	    // and finally the term that is common to all, $\frac{\partial S}{\partial \theta}
	    // dS/d(th).
	    //
	    weight = 1/((*restraints)[i].sigma * (*restraints)[i].sigma);
	    ds_dth = 2*(theta - target_value)*RADTODEG*RADTODEG;
	    w_ds_dth = weight * ds_dth; 

	    if (!(*restraints)[i].fixed_atom_flags[0]) { 
	       idx = 3*((*restraints)[i].atom_index_1);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_k_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_k_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_k_contrib*w_ds_dth); 
	    }
	    if (!(*restraints)[i].fixed_atom_flags[2]) { 
	       idx = 3*((*restraints)[i].atom_index_3);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_m_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_m_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_m_contrib*w_ds_dth); 
	    }

	    // and mid atom
	    if (!(*restraints)[i].fixed_atom_flags[1]) { 
	       idx = 3*((*restraints)[i].atom_index_2);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_l_mid_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_l_mid_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_l_mid_contrib*w_ds_dth); 
	    }
	 }
      }
   }
   // cout << "added " << n_angle_restr << " angle restraint gradients" << endl; 
} 

// Add in the torsion gradients
//
void coot::my_df_torsions(const gsl_vector *v, 
			  void *params, 
			  gsl_vector *df) {

   my_df_torsions_internal(v, params, df, 0);
}


coot::distortion_torsion_gradients_t
coot::fill_distortion_torsion_gradients(const clipper::Coord_orth &P1,
					const clipper::Coord_orth &P2,
					const clipper::Coord_orth &P3,
					const clipper::Coord_orth &P4) {

   coot::distortion_torsion_gradients_t dtg; 
   clipper::Coord_orth a = P2 - P1;
   clipper::Coord_orth b = P3 - P2; 
   clipper::Coord_orth c = P4 - P3;

   double b_lengthsq = b.lengthsq();
   double b_length = sqrt(b_lengthsq); 
   if (b_length < 0.01) { 
      b_length = 0.01; // Garib's stabilization
      b_lengthsq = b_length * b_length; 
   }

   double H = -clipper::Coord_orth::dot(a,c);
   double J =  clipper::Coord_orth::dot(a,b); 
   double K =  clipper::Coord_orth::dot(b,c); 
   double L = 1/b_lengthsq;
   double one_over_b = 1/b_length;

   // a.(bxc)/b
   double E = one_over_b*clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)); 
   // -a.c+(a.b)(b.c)/(b*b)
   double G = H+J*K*L;
   double F = 1/G;
   if (G == 0.0) F = 999999999.9;

   dtg.theta = clipper::Util::rad2d(atan2(E,G));
   // 	    double clipper_theta = 
   // 	       clipper::Util::rad2d(clipper::Coord_orth::torsion(P1, P2, P3, P4));

   // x
   double dH_dxP1 =  c.x(); 
   double dH_dxP2 = -c.x(); 
   double dH_dxP3 =  a.x(); 
   double dH_dxP4 = -a.x(); 

   double dK_dxP1 = 0; 
   double dK_dxP2 = -c.x(); 
   double dK_dxP3 =  c.x() - b.x(); 
   double dK_dxP4 =  b.x(); 

   double dJ_dxP1 = -b.x();
   double dJ_dxP2 =  b.x() - a.x();
   double dJ_dxP3 =  a.x();
   double dJ_dxP4 =  0;

   double dL_dxP1 = 0;
   double dL_dxP2 =  2*(P3.x()-P2.x())*L*L; // check sign
   double dL_dxP3 = -2*(P3.x()-P2.x())*L*L;
   double dL_dxP4 = 0;

   // y 
   double dH_dyP1 =  c.y(); 
   double dH_dyP2 = -c.y(); 
   double dH_dyP3 =  a.y(); 
   double dH_dyP4 = -a.y(); 

   double dK_dyP1 = 0; 
   double dK_dyP2 = -c.y(); 
   double dK_dyP3 =  c.y() - b.y(); 
   double dK_dyP4 =  b.y(); 

   double dJ_dyP1 = -b.y();
   double dJ_dyP2 =  b.y() - a.y();
   double dJ_dyP3 =  a.y();
   double dJ_dyP4 =  0;

   double dL_dyP1 = 0;
   double dL_dyP2 =  2*(P3.y()-P2.y())*L*L; // check sign
   double dL_dyP3 = -2*(P3.y()-P2.y())*L*L;
   double dL_dyP4 = 0;

   // z 
   double dH_dzP1 =  c.z(); 
   double dH_dzP2 = -c.z(); 
   double dH_dzP3 =  a.z(); 
   double dH_dzP4 = -a.z(); 

   double dK_dzP1 = 0; 
   double dK_dzP2 = -c.z(); 
   double dK_dzP3 =  c.z() - b.z(); 
   double dK_dzP4 =  b.z(); 

   double dJ_dzP1 = -b.z();
   double dJ_dzP2 =  b.z() - a.z();
   double dJ_dzP3 =  a.z();
   double dJ_dzP4 =  0;

   double dL_dzP1 = 0;
   double dL_dzP2 =  2*(P3.z()-P2.z())*L*L;
   double dL_dzP3 = -2*(P3.z()-P2.z())*L*L;
   double dL_dzP4 = 0;

   // M
   double dM_dxP1 = -(b.y()*c.z() - b.z()*c.y());
   double dM_dxP2 =  (b.y()*c.z() - b.z()*c.y()) + (a.y()*c.z() - a.z()*c.y());
   double dM_dxP3 =  (b.y()*a.z() - b.z()*a.y()) - (a.y()*c.z() - a.z()*c.y());
   double dM_dxP4 = -(b.y()*a.z() - b.z()*a.y());

   double dM_dyP1 = -(b.z()*c.x() - b.x()*c.z());
   double dM_dyP2 =  (b.z()*c.x() - b.x()*c.z()) + (a.z()*c.x() - a.x()*c.z());
   double dM_dyP3 = -(a.z()*c.x() - a.x()*c.z()) + (b.z()*a.x() - b.x()*a.z());
   double dM_dyP4 = -(b.z()*a.x() - b.x()*a.z());

   double dM_dzP1 = -(b.x()*c.y() - b.y()*c.x());
   double dM_dzP2 =  (b.x()*c.y() - b.y()*c.x()) + (a.x()*c.y() - a.y()*c.x());
   double dM_dzP3 = -(a.x()*c.y() - a.y()*c.x()) + (a.y()*b.x() - a.x()*b.y());
   double dM_dzP4 = -(a.y()*b.x() - a.x()*b.y());

   // E
   double dE_dxP1 = dM_dxP1*one_over_b;
   double dE_dyP1 = dM_dyP1*one_over_b;
   double dE_dzP1 = dM_dzP1*one_over_b; 

   // M = Eb
   double dE_dxP2 = dM_dxP2*one_over_b + E*(P3.x() - P2.x())*L;
   double dE_dyP2 = dM_dyP2*one_over_b + E*(P3.y() - P2.y())*L;
   double dE_dzP2 = dM_dzP2*one_over_b + E*(P3.z() - P2.z())*L;
	    
   double dE_dxP3 = dM_dxP3*one_over_b - E*(P3.x() - P2.x())*L;
   double dE_dyP3 = dM_dyP3*one_over_b - E*(P3.y() - P2.y())*L;
   double dE_dzP3 = dM_dzP3*one_over_b - E*(P3.z() - P2.z())*L;
	    
   double dE_dxP4 = dM_dxP4*one_over_b;
   double dE_dyP4 = dM_dyP4*one_over_b;
   double dE_dzP4 = dM_dzP4*one_over_b;

   double EFF = E*F*F;
   double JL = J*L;
   double KL = K*L;
   double JK = J*K;

   // x
   dtg.dD_dxP1 = F*dE_dxP1 - EFF*(dH_dxP1 + JL*dK_dxP1 + KL*dJ_dxP1 + JK*dL_dxP1);
   dtg.dD_dxP2 = F*dE_dxP2 - EFF*(dH_dxP2 + JL*dK_dxP2 + KL*dJ_dxP2 + JK*dL_dxP2);
   dtg.dD_dxP3 = F*dE_dxP3 - EFF*(dH_dxP3 + JL*dK_dxP3 + KL*dJ_dxP3 + JK*dL_dxP3);
   dtg.dD_dxP4 = F*dE_dxP4 - EFF*(dH_dxP4 + JL*dK_dxP4 + KL*dJ_dxP4 + JK*dL_dxP4);

   // y
   dtg.dD_dyP1 = F*dE_dyP1 - EFF*(dH_dyP1 + JL*dK_dyP1 + KL*dJ_dyP1 + JK*dL_dyP1);
   dtg.dD_dyP2 = F*dE_dyP2 - EFF*(dH_dyP2 + JL*dK_dyP2 + KL*dJ_dyP2 + JK*dL_dyP2);
   dtg.dD_dyP3 = F*dE_dyP3 - EFF*(dH_dyP3 + JL*dK_dyP3 + KL*dJ_dyP3 + JK*dL_dyP3);
   dtg.dD_dyP4 = F*dE_dyP4 - EFF*(dH_dyP4 + JL*dK_dyP4 + KL*dJ_dyP4 + JK*dL_dyP4);

   // z
   dtg.dD_dzP1 = F*dE_dzP1 - EFF*(dH_dzP1 + JL*dK_dzP1 + KL*dJ_dzP1 + JK*dL_dzP1);
   dtg.dD_dzP2 = F*dE_dzP2 - EFF*(dH_dzP2 + JL*dK_dzP2 + KL*dJ_dzP2 + JK*dL_dzP2);
   dtg.dD_dzP3 = F*dE_dzP3 - EFF*(dH_dzP3 + JL*dK_dzP3 + KL*dJ_dzP3 + JK*dL_dzP3);
   dtg.dD_dzP4 = F*dE_dzP4 - EFF*(dH_dzP4 + JL*dK_dzP4 + KL*dJ_dzP4 + JK*dL_dzP4);
  
   return dtg;
} 

// Add in the torsion gradients
//
void coot::my_df_torsions_internal(const gsl_vector *v, 
				   void *params, 
				   gsl_vector *df,
				   bool do_rama_torsions) {
   
   int n_torsion_restr = 0; 
   int idx; 

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   if (restraints->restraints_usage_flag & coot::TORSIONS_MASK) { 
     
      for (int i=0; i<restraints->size(); i++) {
      
	 if ( (*restraints)[i].restraint_type == coot::TORSION_RESTRAINT) {

	    n_torsion_restr++;

	    // double target_value = (*restraints)[i].target_value*DEGTORAD;

	    idx = 3*((*restraints)[i].atom_index_1); 
	    clipper::Coord_orth P1(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_2); 
	    clipper::Coord_orth P2(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_3); 
	    clipper::Coord_orth P3(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_4); 
	    clipper::Coord_orth P4(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));

	    coot::distortion_torsion_gradients_t dtg =
	       fill_distortion_torsion_gradients(P1, P2, P3, P4);

	    if (! do_rama_torsions) { 
	       // 
	       // use period 
	       // 
	       double diff = 99999.9; 
	       double tdiff; 
	       double trial_target; 
	       int per = (*restraints)[i].periodicity;

	       if (dtg.theta < 0.0) dtg.theta += 360.0; 

	       for(int iper=0; iper<per; iper++) { 
		  trial_target = (*restraints)[i].target_value + double(iper)*360.0/double(per); 
		  if (trial_target >= 360.0) trial_target -= 360.0; 
		  tdiff = dtg.theta - trial_target; 
		  if (abs(tdiff) < abs(diff)) { 
		     diff = tdiff;
		  }
	       }
	       if (diff < -180.0) { 
		  diff += 360.; 
	       } else { 
		  if (diff > 180.0) { 
		     diff -= 360.0; 
		  }
	       }
	       //  	    cout << "in df_torsion: theta is " << theta 
	       //  		 <<  " and target is " << (*restraints)[i].target_value 
	       //  		 << " and diff is " << diff 
	       //  		 << " and periodicity: " << (*restraints)[i].periodicity <<  endl;

	       double torsion_scale = (1.0/(1+pow(tan(clipper::Util::d2rad(dtg.theta)),2))) *
		  clipper::Util::rad2d(1.0);

	       double weight = 1/((*restraints)[i].sigma * (*restraints)[i].sigma);

	       // 	    std::cout << "torsion weight: " << weight << std::endl;
	       // 	    std::cout << "torsion_scale : " << torsion_scale << std::endl; 
	       // 	    std::cout << "diff          : " << torsion_scale << std::endl; 

	       // 	    std::cout << "E: " << E << ", F: " << F << ", G: " << G
	       // 		      << ", H: " << H << ", J: "
	       // 		      << J << ", K: " << K << ", L: " << L << " " << std::endl;

	       // 	    std::cout << "gradients: " << std::endl
	       // 		      <<  dD_dxP1 << " " << dD_dyP1 << " " << dD_dzP1 << std::endl
	       // 		      <<  dD_dxP2 << " " << dD_dyP2 << " " << dD_dzP2 << std::endl
	       // 		      <<  dD_dxP3 << " " << dD_dyP3 << " " << dD_dzP3 << std::endl
	       // 		      <<  dD_dxP4 << " " << dD_dyP4 << " " << dD_dzP4 << std::endl
	       // 		      << std::endl;

	       double xP1_contrib = 2.0*diff*dtg.dD_dxP1*torsion_scale * weight;
	       double xP2_contrib = 2.0*diff*dtg.dD_dxP2*torsion_scale * weight;
	       double xP3_contrib = 2.0*diff*dtg.dD_dxP3*torsion_scale * weight;
	       double xP4_contrib = 2.0*diff*dtg.dD_dxP4*torsion_scale * weight;

	       double yP1_contrib = 2.0*diff*dtg.dD_dyP1*torsion_scale * weight;
	       double yP2_contrib = 2.0*diff*dtg.dD_dyP2*torsion_scale * weight;
	       double yP3_contrib = 2.0*diff*dtg.dD_dyP3*torsion_scale * weight;
	       double yP4_contrib = 2.0*diff*dtg.dD_dyP4*torsion_scale * weight;

	       double zP1_contrib = 2.0*diff*dtg.dD_dzP1*torsion_scale * weight;
	       double zP2_contrib = 2.0*diff*dtg.dD_dzP2*torsion_scale * weight;
	       double zP3_contrib = 2.0*diff*dtg.dD_dzP3*torsion_scale * weight;
	       double zP4_contrib = 2.0*diff*dtg.dD_dzP4*torsion_scale * weight;
	    
	       if (! (*restraints)[i].fixed_atom_flags[0]) { 
		  idx = 3*((*restraints)[i].atom_index_1);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP1_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP1_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP1_contrib);
	       }

	       if (! (*restraints)[i].fixed_atom_flags[1]) { 
		  idx = 3*((*restraints)[i].atom_index_2);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP2_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP2_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP2_contrib);
	       }

	       if (! (*restraints)[i].fixed_atom_flags[2]) { 
		  idx = 3*((*restraints)[i].atom_index_3);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP3_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP3_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP3_contrib);
	       }

	       if (! (*restraints)[i].fixed_atom_flags[3]) { 
		  idx = 3*((*restraints)[i].atom_index_4);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP4_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP4_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP4_contrib);
	       }
	    }
	 } 
      }
   }
}

// Add in the torsion gradients
//
void coot::my_df_rama(const gsl_vector *v, 
		      void *params, 
		      gsl_vector *df) {

   // First calculate the torsions:
   // theta = arctan(E/G); 
   // where E = a.(bxc) and G = -a.c + (a.b)(b.c)

   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   if (restraints->restraints_usage_flag & coot::RAMA_PLOT_MASK) { 
     
      for (int i=0; i<restraints->size(); i++) {
      
	 if ( (*restraints)[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) {

	    int idx;
	    coot::simple_restraint rama_restraint = (*restraints)[i];

	    idx = 3*(rama_restraint.atom_index_1);
	    clipper::Coord_orth P1(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*(rama_restraint.atom_index_2); 
	    clipper::Coord_orth P2(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*(rama_restraint.atom_index_3); 
	    clipper::Coord_orth P3(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*(rama_restraint.atom_index_4); 
	    clipper::Coord_orth P4(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*(rama_restraint.atom_index_5); 
	    clipper::Coord_orth P5(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));

	    clipper::Coord_orth a = P2 - P1; 
	    clipper::Coord_orth b = P3 - P2; 
	    clipper::Coord_orth c = P4 - P3;
	    clipper::Coord_orth d = P5 - P4;

	    // New assignements:
	    // TRANS    psi    (2nd N) (2nd CA) (2nd C ) (3nd N)
	    // TRANS    phi    (1st C) (2nd N ) (2nd CA) (2nd C) 
	    // 
	    // So Rama_atoms in this order:
	    //   0       1        2      3         4
	    //  P1      P2       P3     P4        P5
	    // (1st C) (2nd N) (2nd CA) (2nd C) (3rd N)

	    // ---------- phi ------------------
	    // b*b * [ a.(bxc)/b ]
	    double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
	       sqrt( b.lengthsq() );

	    // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
	    double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
	       + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

	    double phi = clipper::Util::rad2d(atan2(E,G));
	    if (phi < 180.0)
	       phi += 360.0;
	    if (phi > 180.0)
	       phi -= 360.0;

	    // ---------- psi ------------------
	    // b*b * [ a.(bxc)/b ]
	    double H = clipper::Coord_orth::dot(b, clipper::Coord_orth::cross(c,d)) *
	       sqrt( c.lengthsq() );

	    // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
	    double I = - clipper::Coord_orth::dot(b,d)*c.lengthsq()
	       + clipper::Coord_orth::dot(b,c)*clipper::Coord_orth::dot(c,d);

	    double psi = clipper::Util::rad2d(atan2(H,I));
	    if (psi < 180.0)
	       psi += 360.0;
	    if (psi > 180.0)
	       psi -= 360.0;


	    if ( clipper::Util::isnan(phi) ) {
	       std::cout << "WARNING: observed torsion phi is a NAN!" << std::endl;
	       // throw an exception
	    } 
	    if ( clipper::Util::isnan(psi) ) {
	       std::cout << "WARNING: observed torsion phi is a NAN!" << std::endl;
	       // throw an exception
	    }

	    double phir = clipper::Util::d2rad(phi);
	    double psir = clipper::Util::d2rad(psi);
	    double R = restraints->rama_prob(phir, psir);
	    
	    // std::cout << "df rama distortion for " << phi << " " << psi << " is "
	    // << R << std::endl;

	    coot::distortion_torsion_gradients_t dtg_phi =
	       fill_distortion_torsion_gradients(P1, P2, P3, P4);

	    coot::distortion_torsion_gradients_t dtg_psi =
	       fill_distortion_torsion_gradients(P2, P3, P4, P5);

	    // Faster to use these, not calculate them above?
	    // 
	    // phir = clipper::Util::d2rad(dtg_phi.theta);
	    // psir = clipper::Util::d2rad(dtg_psi.theta);
	    LogRamachandran::Lgrad lgrd = restraints->rama_grad(phir, psir);

	    double tan_phir = tan(phir);
	    double tan_psir = tan(psir);

	    double multiplier_phi = 10.0/(1.0 + tan_phir*tan_phir) * lgrd.DlogpDphi;
	    double multiplier_psi = 10.0/(1.0 + tan_psir*tan_psir) * lgrd.DlogpDpsi;
	    
	    double xP1_contrib = multiplier_phi*dtg_phi.dD_dxP1;
	    double yP1_contrib = multiplier_phi*dtg_phi.dD_dyP1;
	    double zP1_contrib = multiplier_phi*dtg_phi.dD_dzP1;

	    double xP2_contrib = multiplier_phi*dtg_phi.dD_dxP2;
	    double yP2_contrib = multiplier_phi*dtg_phi.dD_dyP2;
	    double zP2_contrib = multiplier_phi*dtg_phi.dD_dzP2;

	    double xP3_contrib = multiplier_phi*dtg_phi.dD_dxP3;
	    double yP3_contrib = multiplier_phi*dtg_phi.dD_dyP3;
	    double zP3_contrib = multiplier_phi*dtg_phi.dD_dzP3;

	    double xP4_contrib = multiplier_phi*dtg_phi.dD_dxP4;
	    double yP4_contrib = multiplier_phi*dtg_phi.dD_dyP4;
	    double zP4_contrib = multiplier_phi*dtg_phi.dD_dzP4;

	    // The class variable gives a misleading name here for the
	    // follwing blocks. P2 is in postion 1 for dtg_phi, P3 is
	    // in position 2, P4 is called in the 3rd position (and
	    // P5 in 4th).

	    xP2_contrib += multiplier_psi * dtg_psi.dD_dxP1;
	    yP2_contrib += multiplier_psi * dtg_psi.dD_dyP1;
	    zP2_contrib += multiplier_psi * dtg_psi.dD_dzP1;
	    
	    xP3_contrib += multiplier_psi * dtg_psi.dD_dxP2;
	    yP3_contrib += multiplier_psi * dtg_psi.dD_dyP2;
	    zP3_contrib += multiplier_psi * dtg_psi.dD_dzP2;
	    
	    xP4_contrib += multiplier_psi * dtg_psi.dD_dxP3;
	    yP4_contrib += multiplier_psi * dtg_psi.dD_dyP3;
	    zP4_contrib += multiplier_psi * dtg_psi.dD_dzP3;

	    if (0) { 
	       xP2_contrib = 0.0;
	       yP2_contrib = 0.0;
	       zP2_contrib = 0.0;
	       
	       xP3_contrib = 0.0;
	       yP3_contrib = 0.0;
	       zP3_contrib = 0.0;
	       
	       xP4_contrib = 0.0;
	       yP4_contrib = 0.0;
	       zP4_contrib = 0.0;
	    }
	    
	    double xP5_contrib = multiplier_psi*dtg_psi.dD_dxP4;
	    double yP5_contrib = multiplier_psi*dtg_psi.dD_dyP4;
	    double zP5_contrib = multiplier_psi*dtg_psi.dD_dzP4;


	    if (! (*restraints)[i].fixed_atom_flags[0]) { 
	       idx = 3*((*restraints)[i].atom_index_1);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP1_contrib);
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP1_contrib);
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP1_contrib);
	    }

	    if (! (*restraints)[i].fixed_atom_flags[1]) { 
	       idx = 3*((*restraints)[i].atom_index_2);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP2_contrib);
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP2_contrib);
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP2_contrib);
	    }

	    if (! (*restraints)[i].fixed_atom_flags[2]) { 
	       idx = 3*((*restraints)[i].atom_index_3);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP3_contrib);
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP3_contrib);
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP3_contrib);
	    }

	    if (! (*restraints)[i].fixed_atom_flags[3]) { 
	       idx = 3*((*restraints)[i].atom_index_4);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP4_contrib);
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP4_contrib);
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP4_contrib);
	    }

	    if (! (*restraints)[i].fixed_atom_flags[4]) { 
	       idx = 3*((*restraints)[i].atom_index_5);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP5_contrib);
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP5_contrib);
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP5_contrib);
	    }
	 }
      }
   }
}

//  the chiral volumes
void 
coot::my_df_chiral_vol(const gsl_vector *v, void *params, gsl_vector *df) { 

   int n_chiral_vol_restr = 0;
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   int idx;
   double cv;
   double distortion;
   
   if (restraints->restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) {
      // if (0) {
      
      for (int i=0; i<restraints->size(); i++) {
	 
	 if ( (*restraints)[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {

	    n_chiral_vol_restr++;

	    idx = 3*(*restraints)[i].atom_index_centre;
	    clipper::Coord_orth centre(gsl_vector_get(v, idx),
				       gsl_vector_get(v, idx+1),
				       gsl_vector_get(v, idx+2));

	    idx = 3*( (*restraints)[i].atom_index_1);
	    clipper::Coord_orth a1(gsl_vector_get(v, idx),
				   gsl_vector_get(v, idx+1),
				   gsl_vector_get(v, idx+2));
	    idx = 3*( (*restraints)[i].atom_index_2);
	    clipper::Coord_orth a2(gsl_vector_get(v, idx),
				   gsl_vector_get(v, idx+1),
				   gsl_vector_get(v, idx+2));
	    idx = 3*( (*restraints)[i].atom_index_3);
	    clipper::Coord_orth a3(gsl_vector_get(v, idx),
				   gsl_vector_get(v, idx+1),
				   gsl_vector_get(v, idx+2));

	    clipper::Coord_orth a = a1 - centre;
	    clipper::Coord_orth b = a2 - centre;
	    clipper::Coord_orth c = a3 - centre;

	    cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));

	    distortion = cv - (*restraints)[i].target_chiral_volume;
	    
// 	    std::cout << "---- xxx ---- DEBUG:: chiral volume deriv: " 
// 		      << cv << " chiral distortion " 
// 		      << distortion << "\n";
	    // distortion /= ((*restraints)[i].sigma * (*restraints)[i].sigma);

	    double P0_x_contrib =
	       - (b.y()*c.z() - b.z()*c.y())
	       - (a.z()*c.y() - a.y()*c.z())
	       - (a.y()*b.z() - a.z()*b.y());
		  
	    double P0_y_contrib = 
	       - (b.z()*c.x() - b.x()*c.z())
	       - (a.x()*c.z() - a.z()*c.x())
	       - (a.z()*b.x() - a.x()*b.z());

	    double P0_z_contrib = 
	       - (b.x()*c.y() - b.y()*c.x())
	       - (a.y()*c.x() - a.x()*c.y())
	       - (a.x()*b.y() - a.y()*b.x());

	    double P1_x_contrib = b.y()*c.z() - b.z()*c.y();
	    double P1_y_contrib = b.z()*c.x() - b.x()*c.z();
	    double P1_z_contrib = b.x()*c.y() - b.y()*c.x();

	    double P2_x_contrib = a.z()*c.y() - a.y()*c.z();
	    double P2_y_contrib = a.x()*c.z() - a.z()*c.x();
	    double P2_z_contrib = a.y()*c.x() - a.x()*c.y();

	    double P3_x_contrib = a.y()*b.z() - a.z()*b.y();
	    double P3_y_contrib = a.z()*b.x() - a.x()*b.z();
	    double P3_z_contrib = a.x()*b.y() - a.y()*b.x();

	    double s = 2*distortion/((*restraints)[i].sigma * (*restraints)[i].sigma);

	    if (!(*restraints)[i].fixed_atom_flags[0]) { 
	       idx = 3*( (*restraints)[i].atom_index_centre);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P0_x_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P0_y_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P0_z_contrib);
	    }
	       
	    if (!(*restraints)[i].fixed_atom_flags[1]) { 
	       idx = 3*( (*restraints)[i].atom_index_1);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P1_x_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P1_y_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P1_z_contrib);
	    }

	    if (!(*restraints)[i].fixed_atom_flags[2]) { 
	       idx = 3*( (*restraints)[i].atom_index_2);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P2_x_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P2_y_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P2_z_contrib);
	    }

	    if (!(*restraints)[i].fixed_atom_flags[3]) { 
	       idx = 3*( (*restraints)[i].atom_index_3);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P3_x_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P3_y_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P3_z_contrib);
	    }
	 }
      }
   }
} 

// manipulate the gsl_vector_get *df
// 
void
coot::my_df_planes(const gsl_vector *v, 
		   void *params, 
		   gsl_vector *df) {

   
   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   


   int idx; 
   int n_plane_restr = 0; // debugging counter

   coot::plane_distortion_info_t plane_info;

   
   if (restraints->restraints_usage_flag & coot::PLANES_MASK) {
      
      int n_atoms;
      double devi_len;
      double weight;

      for (int i=0; i<restraints->size(); i++) {
       
	 if ( (*restraints)[i].restraint_type == coot::PLANE_RESTRAINT) {

	    simple_restraint plane_restraint = (*restraints)[i];
	    plane_info = distortion_score_plane_internal(plane_restraint, v);
	    n_atoms = plane_restraint.atom_index.size();
	    for (int j=0; j<n_atoms; j++) {
	       if (! (*restraints)[i].fixed_atom_flags[j] ) { 
		  n_plane_restr++;
		  idx = 3*(plane_restraint.atom_index[j]);
		  devi_len =
		     plane_info.abcd[0]*gsl_vector_get(v,idx  ) + 
		     plane_info.abcd[1]*gsl_vector_get(v,idx+1) +
		     plane_info.abcd[2]*gsl_vector_get(v,idx+2) -
		     plane_info.abcd[3];
		  weight = 1/((*restraints)[i].sigma * (*restraints)[i].sigma);

		  clipper::Grad_orth<double> d(2.0 * weight * devi_len * plane_info.abcd[0],
					       2.0 * weight * devi_len * plane_info.abcd[1],
					       2.0 * weight * devi_len * plane_info.abcd[2]);

		  if (!(*restraints)[i].fixed_atom_flags[j]) { 
		     gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + d.dx());
		     gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + d.dy());
		     gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + d.dz());
// 		     if (n_plane_restr == 1) {
// 			clipper::Coord_orth normal(plane_info.abcd[0],
// 						   plane_info.abcd[1],
// 						   plane_info.abcd[2]);
// 		     }
		  }
	       }
	    }
	 }
      }
   }
}

double
coot::restraints_container_t::electron_density_score_at_point(const clipper::Coord_orth &ao) const {
      double dv; 
      
      clipper::Coord_frac af = ao.coord_frac(map.cell()); 
      clipper::Coord_map  am = af.coord_map(map.grid_sampling()); 
      // clipper::Interp_linear::interp(map, am, dv); 
      clipper::Interp_cubic::interp(map, am, dv); 
      
      return dv;  
}

clipper::Grad_orth<double>
coot::restraints_container_t::electron_density_gradient_at_point(const clipper::Coord_orth &ao) const {
   
   clipper::Grad_map<double> grad;
   double dv;
   
   clipper::Coord_frac af = ao.coord_frac(map.cell()); 
   clipper::Coord_map  am = af.coord_map(map.grid_sampling()); 
   clipper::Interp_cubic::interp_grad(map, am, dv, grad);
   clipper::Grad_frac<double> grad_frac = grad.grad_frac(map.grid_sampling());
   return grad_frac.grad_orth(map.cell());
} 


// Ah, but (c.f. distortion) we want to return a low value for a good
// fit and a high one for a bad.
double
coot::electron_density_score(const gsl_vector *v, void *params) { 

   // We sum to the score and negate.  That will do?
   // 
   double score = 0; 
   // double e = 2.718281828; 
   
   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 

   if (restraints->include_map_terms() == 1) { 
      
      // convert from variables to coord_orths of where the atoms are
      
      for (unsigned int i=0; i< v->size; i += 3) { 
	 int iat = i/3;
	 if (restraints->use_map_gradient_for_atom[iat] == 1) {
	    bool use_it = 1;
// 	    for (unsigned int ifixed=0; ifixed<restraints->fixed_atom_indices.size(); ifixed++) {
// 	       if (restraints->fixed_atom_indices[ifixed] == iat) { 
// 		  std::cout << "ignoring density term for atom " << iat << std::endl;
// 		  use_it = 0;
// 		  break;
// 	       }
// 	    }
	    if (use_it) { 
	       clipper::Coord_orth ao(gsl_vector_get(v,i), 
				      gsl_vector_get(v,i+1), 
				      gsl_vector_get(v,i+2));
	       
	       score += restraints->Map_weight() *
		  restraints->atom_z_weight[iat] *
		  restraints->electron_density_score_at_point(ao);
	    }
	 }
      }
   }
   
   // return pow(e,-score*0.01);
   return -score;

}

// Note that the gradient for the electron density is opposite to that
// of the gradient for the geometry (consider a short bond on the edge
// of a peak - in that case the geometry gradient will be negative as
// the bond is lengthened and the electron density gradient will be
// positive).
//
// So we want to change that positive gradient for a low score when
// the atoms cooinside with the density - hence the contributions that
// we add are negated.
// 
void coot::my_df_electron_density (gsl_vector *v, 
				   void *params, 
				   gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 

   if (restraints->include_map_terms() == 1) { 


      clipper::Grad_orth<double> grad_orth;
      float scale = restraints->Map_weight();
      float zs;
      
      for (unsigned int i=0; i<v->size; i+=3) {

	 int iat = i/3;
// 	 std::cout << "restraints->use_map_gradient_for_atom[" << iat << "] == "
// 		   << restraints->use_map_gradient_for_atom[iat] << "\n";
	 
	 if (restraints->use_map_gradient_for_atom[iat] == 1) {

	    clipper::Coord_orth ao(gsl_vector_get(v,i), 
				   gsl_vector_get(v,i+1), 
				   gsl_vector_get(v,i+2));
	    
	    // std::cout << "gradients: " << grad_orth.format() << std::endl;
	    // std::cout << "adding to " << gsl_vector_get(df, i  ) << std::endl;
	    // add this density term to the gradient
	    //
	    // 
	    grad_orth = restraints->electron_density_gradient_at_point(ao);
	    zs = scale * restraints->atom_z_weight[iat];

	    if (0) { 
	       std::cout << "electron density df: adding "
			 <<  - zs * grad_orth.dx() << " "
			 <<  - zs * grad_orth.dy() << " "
			 <<  - zs * grad_orth.dz() << " to "
			 <<  gsl_vector_get(df, i  ) << " "
			 <<  gsl_vector_get(df, i+1) << " "
			 <<  gsl_vector_get(df, i+2) << "\n";
	    }
	    
	    gsl_vector_set(df, i,   gsl_vector_get(df, i  ) - zs * grad_orth.dx());
	    gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) - zs * grad_orth.dy());
	    gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) - zs * grad_orth.dz());
	 } else {
	    // atom is private	    
// 	    std::cout << "  Not adding elecron density for atom "
// 		      << restraints->atom[iat]->GetChainID() << " "
// 		      << restraints->atom[iat]->GetSeqNum() << " "
// 		      << restraints->atom[iat]->GetAtom() << std::endl;
	 } 
      }
   }
}

void coot::my_df_electron_density_old (gsl_vector *v, 
				   void *params, 
				   gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 

   if (restraints->include_map_terms() == 1) { 

      double new_S_minu, new_S_plus, tmp, val; 

      cout << "density_gradients" << endl; 
      for (unsigned int i=0; i<v->size; i++) { 
      
	 tmp = gsl_vector_get(v, i); 
	 gsl_vector_set(v, i, tmp+0.01); 
	 new_S_plus = coot::electron_density_score(v, params); 
	 gsl_vector_set(v, i, tmp-0.01); 
	 new_S_minu = coot::electron_density_score(v, params); 
	 // new_S_minu = 2*tmp - new_S_plus; 

	 // restore the initial value: 
	 gsl_vector_set(v, i, tmp);

	 val = (new_S_plus - new_S_minu)/(2*0.01); 
	 cout << "density gradient: " << i << " " << val << endl;

	 // add this density term to the gradient
	 gsl_vector_set(df, i, gsl_vector_get(df, i) + val);
      } 
   }
}

// Include the gradients of the terms that describe the deviation from
// the starting position.
// 
void my_df_start_pos (const gsl_vector *v, 
		      void *params, 
		      gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   
   if (int(v->size) != int(restraints->init_positions_size()) ) {
      cout << "very worry. A bug. " << v->size << " "
	   << restraints->init_positions_size() << endl;
   } else {
      // These are similar to the bond terms, however, the d, distance
      // between the atoms is tiny when we are at the beginning of the
      // refinement.
      //
      // We do x, y, z for each round
      //
      double val; 
      for (unsigned int i=0; i<v->size; i += 3) {
	 //
	 clipper::Coord_orth current_pos(gsl_vector_get(v,i),
					 gsl_vector_get(v,i+1),
					 gsl_vector_get(v,i+2));
	 clipper::Coord_orth start_pos(restraints->initial_position(i),
				       restraints->initial_position(i+1),
				       restraints->initial_position(i+2));

	 double b = clipper::Coord_orth::length(start_pos,current_pos); 

	 b = b < 0.01 ? 0.01 : b; // Stabilization

	 val = 0.01*(gsl_vector_get(v, i) - restraints->initial_position(i));
	 gsl_vector_set(df, i, gsl_vector_get(df, i) + val);
	 // cout << "added in my_df_start_pos: " << val << endl; 

	 val = 0.1*(gsl_vector_get(v, i+1) - restraints->initial_position(i+1));
	 gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) + val);

	 val = 0.1*(gsl_vector_get(v, i+2) - restraints->initial_position(i+2));
	 gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) + val);
	 
      }
   }
}

// Compute both f and df together.
void coot::my_fdf(const gsl_vector *x, void *params, 
		  double *f, gsl_vector *df) { 

   *f = coot::distortion_score(x, params); 
    coot::my_df(x, params, df); 
}
  

// We need to fill restraints_vec (which is a vector of
// simple_restraint) using the coordinates () and the dictionary of
// restraints, protein_geometry geom.
//
// The plan is to get a list of residues, and for each of those
// residues, look in geom for a residue of that type and if found,
// fill restraints_vec appropriately.
int
coot::restraints_container_t::make_restraints(const coot::protein_geometry &geom,
					      coot::restraint_usage_Flags flags_in, 
					      short int do_residue_internal_torsions,
					      float rama_plot_target_weight,
					      bool do_rama_plot_restraints, 
					      coot::pseudo_restraint_bond_type sec_struct_pseudo_bonds) {

//   std::cout << "----- make restraints called with link_torsions_restraints_type: "
// 	    << link_torsions_restraints_type << std::endl;

   int iret = 0;
   restraints_usage_flag = flags_in; // also set in minimize() and geometric_distortions()
   mark_OXT(geom);
   iret += make_monomer_restraints(geom, do_residue_internal_torsions);

   iret += make_link_restraints(geom, do_rama_plot_restraints);
   //  don't do torsions, ramas maybe.   
   iret += make_flanking_atoms_restraints(geom, do_rama_plot_restraints);
   int iret_prev = restraints_vec.size();

   if (sec_struct_pseudo_bonds == coot::HELIX_PSEUDO_BONDS) {
      make_helix_pseudo_bond_restraints();
   } 
   if (sec_struct_pseudo_bonds == coot::STRAND_PSEUDO_BONDS) {
      make_strand_pseudo_bond_restraints();
   }

   if (restraints_usage_flag & coot::NON_BONDED_MASK) {
      if (iret_prev > 0) {
	 int n_nbcr = make_non_bonded_contact_restraints();
	 std::cout << "INFO:: made " << n_nbcr << " non-bonded restraints\n";
      }
   }

   return restraints_vec.size();
}

// This only marks the first OXT atom we find that has all its
// reference atoms (which is reasonable (I hope)).
// 
void 
coot::restraints_container_t::mark_OXT(const coot::protein_geometry &geom) { 

   std::string oxt(" OXT");
   for (int i=0; i<n_atoms; i++) { 
      if (std::string(atom[i]->name) == oxt) { 

	 // OK, we need to find the reference positions for this OXT,
	 // the N, CA, C and O in that order.
	 
	 CResidue *residue = atom[i]->residue;
	 CAtom *res_atom = NULL;
	 if (residue) {

	    std::string residue_type(residue->GetResName());
	    if (! geom.OXT_in_residue_restraints_p(residue_type)) { 

	    
	       res_atom = residue->GetAtom(" N  ");
	       if (res_atom) 
		  oxt_reference_atom_pos.push_back(clipper::Coord_orth(res_atom->x, res_atom->y, res_atom->z));
	       res_atom = residue->GetAtom(" CA ");
	       if (res_atom) 
		  oxt_reference_atom_pos.push_back(clipper::Coord_orth(res_atom->x, res_atom->y, res_atom->z));
	       res_atom = residue->GetAtom(" C  ");
	       if (res_atom) 
		  oxt_reference_atom_pos.push_back(clipper::Coord_orth(res_atom->x, res_atom->y, res_atom->z));
	       res_atom = residue->GetAtom(" O  ");
	       if (res_atom) 
		  oxt_reference_atom_pos.push_back(clipper::Coord_orth(res_atom->x, res_atom->y, res_atom->z));

	       have_oxt_flag = 1;
	       oxt_index = i;
	       break;
	    }
	 }
      }
   }
}


std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2) const {

   std::vector<bool> r(2,0);
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index1 == fixed_atom_indices[ifixed])
	 r[0] = 1;
      if (index2 == fixed_atom_indices[ifixed])
	 r[1] = 1;
   }
   return r;
} 

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2, int index3) const {

   std::vector<bool> r(3,0);
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index1 == fixed_atom_indices[ifixed])
	 r[0] = 1;
      if (index2 == fixed_atom_indices[ifixed])
	 r[1] = 1;
      if (index3 == fixed_atom_indices[ifixed])
	 r[2] = 1;
   }
   return r;
} 

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(int index1, int index2, int index3, int index4) const {

   std::vector<bool> r(4,0);
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      if (index1 == fixed_atom_indices[ifixed])
	 r[0] = 1;
      if (index2 == fixed_atom_indices[ifixed])
	 r[1] = 1;
      if (index3 == fixed_atom_indices[ifixed])
	 r[2] = 1;
      if (index4 == fixed_atom_indices[ifixed])
	 r[3] = 1;
   }
   return r;
}

std::vector<bool>
coot::restraints_container_t::make_fixed_flags(const std::vector<int> &indices) const {

   std::vector<bool> r(indices.size(), 0);
   for (unsigned int ifixed=0; ifixed<fixed_atom_indices.size(); ifixed++) {
      for (unsigned int i_index=0; i_index<indices.size(); i_index++) {
	 if (indices[i_index] == fixed_atom_indices[ifixed])
	    r[i_index] = 1;
      }
   }
   return r;
} 


void
coot::restraints_container_t::make_helix_pseudo_bond_restraints() {

   // This method of making pseudo bonds relies on the residue range
   // being continuous in sequence number (seqNum) and no insertion
   // codes messing up the number scheme.  If these are not the case
   // then we will make bonds between the wrongs atoms (of the wrong
   // residues)... a more sophisticated algorithm would be needed.
   // 
   // The method ignores residues with alt confs.
   //
   float pseudo_bond_esd = 0.04; // just a guess
   int selHnd = mol->NewSelection();
   int nSelResidues;
   PPCResidue SelResidue;
   PPCAtom res_1_atoms = NULL;
   PPCAtom res_2_atoms = NULL;
   int n_res_1_atoms;
   int n_res_2_atoms;
   int index1 = -1; 
   int index2 = -1; 
   mol->Select (selHnd, STYPE_RESIDUE, 1, // .. TYPE, iModel
		(char *) chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"",  // altLocs
		SKEY_NEW // selection key
		);
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 0) {
      for (int i=4; i<nSelResidues; i++) {
	 // nN -> (n-4)O 2.91
         // nN -> (n-3)O 3.18
         // nO -> (n+3)N 3.18  // symmetric.  No need to specify both forwards
         // nO -> (n+4)N 2.91  // and backwards directions.
	 SelResidue[i]->GetAtomTable(res_1_atoms, n_res_1_atoms);
	 for (int iat1=0; iat1<n_res_1_atoms; iat1++) {
	    std::string at_1_name(res_1_atoms[iat1]->name);
	    
	    if (at_1_name == " N  ") {
	       CResidue *contact_res = SelResidue[i-4];
	       if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 4)) {
		  contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		  for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
		     std::string at_2_name(res_2_atoms[iat2]->name);
		     if (at_2_name == " O  ") {
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			add(BOND_RESTRAINT, index1, index2, fixed_flags,
			    2.91, pseudo_bond_esd, 1.2);
			std::cout << "Helix Bond restraint (" << res_1_atoms[iat1]->name << " "
				  << res_1_atoms[iat1]->GetSeqNum() << ") to ("
				  << res_2_atoms[iat2]->name << " "
				  << res_2_atoms[iat2]->GetSeqNum() << ") 2.91" << std::endl;
		     }
		  }
	       }

	       contact_res = SelResidue[i-3];
	       if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 3)) {
		  contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		  for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
		     std::string at_2_name(res_2_atoms[iat2]->name);
		     if (at_2_name == " O  ") {
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			add(BOND_RESTRAINT, index1, index2, fixed_flags,
			    3.18, pseudo_bond_esd, 1.2);
			std::cout << "Helix Bond restraint (" << res_1_atoms[iat1]->name << " "
				  << res_1_atoms[iat1]->GetSeqNum() << ") to ("
				  << res_2_atoms[iat2]->name << " "
				  << res_2_atoms[iat2]->GetSeqNum() << ") 3.18" << std::endl;
		     }
		  }
	       }
	       
	    } 
	 } 
      }
   }
   mol->DeleteSelection(selHnd);
}


void
coot::restraints_container_t::make_strand_pseudo_bond_restraints() { 

   // This method of making pseudo bonds relies on the residue range
   // being continuous in sequence number (seqNum) and no insertion
   // codes messing up the number scheme.  If these are not the case
   // then we will make bonds between the wrongs atoms (of the wrong
   // residues)... a more sophisticated algorithm would be needed.
   // 
   // The method ignores residues with alt confs.
   //
   float pseudo_bond_esd = 0.08; // just a guess
   int selHnd = mol->NewSelection();
   int nSelResidues;
   PPCResidue SelResidue;
   PPCAtom res_1_atoms = NULL;
   PPCAtom res_2_atoms = NULL;
   PPCAtom res_3_atoms = NULL;
   int n_res_1_atoms;
   int n_res_2_atoms;
   int n_res_3_atoms;
   int index1 = -1; 
   int index2 = -1; 
   int index3 = -1; 
   mol->Select (selHnd, STYPE_RESIDUE, 1, // .. TYPE, iModel
		(char *) chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"",  // altLocs
		SKEY_NEW // selection key
		);
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 0) {
      for (int i=1; i<nSelResidues; i++) {
         // nO -> (n-1)O 4.64  // symmetric.  No need to specify both forwards
	 // Angle nO-(n+1)O-(n+2)O: 
	 SelResidue[i]->GetAtomTable(res_1_atoms, n_res_1_atoms);
	 for (int iat1=0; iat1<n_res_1_atoms; iat1++) {
	    std::string at_1_name(res_1_atoms[iat1]->name);
	    // O Pseudo bonds and angles
	    if (at_1_name == " O  ") {
	       CResidue *contact_res = SelResidue[i-1];
	       if (SelResidue[i]->GetSeqNum() == (contact_res->GetSeqNum() + 1)) {
		  contact_res->GetAtomTable(res_2_atoms, n_res_2_atoms);
		  for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
		     std::string at_2_name(res_2_atoms[iat2]->name);
		     if (at_2_name == " O  ") {
			std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
			res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
			res_2_atoms[iat2]->GetUDData(udd_atom_index_handle, index2);
			add(BOND_RESTRAINT, index1, index2, fixed_flags,
			    4.64, pseudo_bond_esd, 1.2);
			std::cout << "Strand Bond restraint ("
				  << res_1_atoms[iat1]->name << " "
				  << res_1_atoms[iat1]->GetSeqNum() << ") to ("
				  << res_2_atoms[iat2]->name << " "
				  << res_2_atoms[iat2]->GetSeqNum() << ") 4.64" << std::endl;

			// now the pseudo angle
			if (i<(nSelResidues-1)) { 
			   CResidue *contact_res_2 = SelResidue[i+1];
			   if (SelResidue[i]->GetSeqNum() == (contact_res_2->GetSeqNum() - 1)) {
			      contact_res_2->GetAtomTable(res_3_atoms, n_res_3_atoms);
			      for (int iat3=0; iat3<n_res_3_atoms; iat3++) {
				 std::string at_3_name(res_3_atoms[iat3]->name);
				 if (at_3_name == " O  ") {
				    std::vector<bool> fixed_flag =
				       make_fixed_flags(index2, index1, index3);
				    res_3_atoms[iat3]->GetUDData(udd_atom_index_handle, index3);
				    // 98.0 degrees
				    add(ANGLE_RESTRAINT, index2, index1, index3,
					fixed_flag, 98.0, 0.5, 1.2);
				    std::cout << "Strand Angle restraint ("
					      << res_1_atoms[iat1]->name << " "
					      << res_1_atoms[iat1]->GetSeqNum() << ") to ("
					      << res_2_atoms[iat2]->name << " "
					      << res_2_atoms[iat2]->GetSeqNum()
					      << ") to ("
					      << res_3_atoms[iat3]->name << " "
					      << res_3_atoms[iat3]->GetSeqNum()
					      << ") 98.0 " << std::endl;
				    break;
				 }
			      }
			   }
			}
			break;
		     }
		  }
	       }
	    } // end of O

	    // Now make a CA-CA-CA pseudo angle of 120 degrees
	    if (at_1_name == " CA ") {
	       CResidue *contact_res_2 = SelResidue[i-1];
	       if (SelResidue[i]->GetSeqNum() == (contact_res_2->GetSeqNum() + 1)) {
		  contact_res_2->GetAtomTable(res_2_atoms, n_res_2_atoms);
		  for (int iat2=0; iat2<n_res_2_atoms; iat2++) {
		     std::string at_2_name(res_2_atoms[iat2]->name);
		     if (at_2_name == " CA ") {
			if (i<(nSelResidues-1)) {
			   CResidue *contact_res_3 = SelResidue[i+1];
			   if (SelResidue[i]->GetSeqNum() == (contact_res_3->GetSeqNum() - 1)) {
			      contact_res_3->GetAtomTable(res_3_atoms, n_res_3_atoms);
			      for (int iat3=0; iat3<n_res_3_atoms; iat3++) {
				 std::string at_3_name(res_3_atoms[iat3]->name);
				 if (at_3_name == " CA ") {
				    std::vector<bool> fixed_flag =
				       make_fixed_flags(index1, index2, index3);
				    res_1_atoms[iat1]->GetUDData(udd_atom_index_handle, index1);
				    res_2_atoms[iat3]->GetUDData(udd_atom_index_handle, index2);
				    res_3_atoms[iat3]->GetUDData(udd_atom_index_handle, index3);
				    add(ANGLE_RESTRAINT, index2, index1, index3,
					fixed_flag, 120.0, 0.5, 1.2);
				    std::cout << "Strand Angle restraint ("
					      << res_1_atoms[iat1]->name << " "
					      << res_1_atoms[iat1]->GetSeqNum() << ") to ("
					      << res_2_atoms[iat2]->name << " "
					      << res_2_atoms[iat2]->GetSeqNum()
					      << ") to ("
					      << res_3_atoms[iat3]->name << " "
					      << res_3_atoms[iat3]->GetSeqNum()
					      << ") 120.0 " << std::endl;
				    break;
				 }
			      }
			   }
			}
			break;
		     }
		  }
	       }
	    }
	 } 
      }
   }
   mol->DeleteSelection(selHnd);
}


int
coot::restraints_container_t::make_monomer_restraints(const coot::protein_geometry &geom,
						      short int do_residue_internal_torsions) {

   if (from_residue_vector)
      return make_monomer_restraints_from_res_vec(geom, do_residue_internal_torsions);
   else
      return make_monomer_restraints_by_linear(geom, do_residue_internal_torsions);

}

int
coot::restraints_container_t::make_monomer_restraints_by_linear(const coot::protein_geometry &geom,
								short int do_residue_internal_torsions) {
   
   int iret = 0;
   
   int selHnd = mol->NewSelection();
   int nSelResidues;
   coot::restraints_container_t::restraint_counts_t sum;

   mol->Select (selHnd, STYPE_RESIDUE, 1, // .. TYPE, iModel
		(char *) chain_id_save.c_str(), // Chain(s)
		istart_res, "*", // starting res
		iend_res,   "*", // ending   res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",  // altLocs
		SKEY_NEW // selection key
		);
   SelResidue_active = NULL;
   mol->GetSelIndex (selHnd, SelResidue_active, nSelResidues);
//    std::cout << "INFO:: GetSelIndex returned " << nSelResidues
// 	     << " residues (monomer restraints) " << std::endl;
   // save the (new (7Nov2003)) class variables (used in non_bonded
   // stuff) that keep the "active" (as opposed to "flanking") residues:
   nSelResidues_active = nSelResidues; 
   if (nSelResidues > 0) { 
      for (int i=0; i<nSelResidues; i++) {
	 if (SelResidue_active[i]) {
	    coot::restraints_container_t::restraint_counts_t local = 
	       make_monomer_restraints_by_residue(SelResidue_active[i], geom,
						  do_residue_internal_torsions);
	    sum += local;
	 }
      }
   } else {
      std::cout << "get_monomer_restraints: There were no residues selected!? "
		<< std::endl;
   }
   // mol->DeleteSelection(selHnd); // -> makes crash! 6-feb-2004.
   //
   // This is because SelResidue_active is used elsewhere.
   // 
   // This should go into the destructor, I guess.

   sum.report(do_residue_internal_torsions);
   std::cout << "created " << restraints_vec.size() << " restraints" << std::endl;
   std::cout << std::endl; 
   return iret; // return 1 on success.  Hmm... how is this set? (and subsequently used?)
}

int
coot::restraints_container_t::make_monomer_restraints_from_res_vec(const coot::protein_geometry &geom,
								   short int do_residue_internal_torsions) {
   
   int iret = 0;

   coot::restraints_container_t::restraint_counts_t sum;

   for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
      coot::restraints_container_t::restraint_counts_t local = 
	 make_monomer_restraints_by_residue(residues_vec[ir].second, geom,
					    do_residue_internal_torsions);
      sum += local;
   } 
   sum.report(do_residue_internal_torsions);
   std::cout << "created " << restraints_vec.size() << " restraints" << std::endl;
   std::cout << std::endl;
   return iret;
} 



coot::restraints_container_t::restraint_counts_t
coot::restraints_container_t::make_monomer_restraints_by_residue(CResidue *residue_p,
								 const protein_geometry &geom,
								 short int do_residue_internal_torsions) {

   coot::restraints_container_t::restraint_counts_t local;
   int i_no_res_atoms;
   PPCAtom res_selection = NULL;
   
   // idr: index dictionary residue 
   for (int idr=0; idr<geom.size(); idr++) {
      std::string pdb_resname(residue_p->name);
      if (pdb_resname == "UNK") pdb_resname = "ALA";
      // if (geom[idr].comp_id == pdb_resname) {
      // old style comp_id usage
      // if (dictionary_name_matches_coords_resname(geom[idr].comp_id,pdb_resname)) {

      // OK, we need the 3 letter code for carbohydrates, the
      // comp_id for nucleotides:
      //
      // comp_id 3-letter-code name group
      // Ar         A        'Adenosine                    ' RNA                33  22 .
      // GAL-b-D    GAL      'beta_D_galactose             ' D-pyranose         24  12 .

      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr),
						 pdb_resname) ||
	  dictionary_name_matches_coords_resname(geom[idr].comp_id, pdb_resname)) {

	 
// 	 std::cout << "DEBUG:: ------------- dict/pdb name matches " << pdb_resname
// 		   << " --------------- " << std::endl; 

	 // now get a list of atoms in that residue
	 // (SelResidue[i]) and compare them to the atoms in
	 // geom[idr].bond_restraint[ib].

	 residue_p->GetAtomTable(res_selection, i_no_res_atoms);
		  
	 if (i_no_res_atoms > 0) {

	    // 		     std::cout << "   bonds... " << std::endl;
	    local.n_bond_restraints += add_bonds(idr, res_selection, i_no_res_atoms,
						 residue_p, geom);
	    // 		     std::cout << "   angles... " << std::endl;
	    local.n_angle_restraints += add_angles(idr, res_selection, i_no_res_atoms,
						   residue_p, geom);
	    if (do_residue_internal_torsions) { 
	       // 	 	        std::cout << "   torsions... " << std::endl;
	       std::string residue_type = residue_p->GetResName();
	       if (residue_type != "PRO") 
		  local.n_torsion_restr += add_torsions(idr, res_selection, i_no_res_atoms,
							residue_p, geom);
	    }
	    // 		     std::cout << "   planes... " << std::endl;
	    local.n_plane_restraints += add_planes(idr, res_selection, i_no_res_atoms,
						   residue_p, geom);

	    local.n_chiral_restr += add_chirals(idr, res_selection, i_no_res_atoms, 
						residue_p, geom);

	 }
      }
   }
   return local;
} 

int
coot::restraints_container_t::make_link_restraints(const coot::protein_geometry &geom,
						   bool do_rama_plot_restraints) {

   if (from_residue_vector)
      return make_link_restraints_from_res_vec(geom, do_rama_plot_restraints);
   else 
      return make_link_restraints_by_linear(geom, do_rama_plot_restraints); // conventional

   int ir = 0;
   return ir; 
}

int
coot::restraints_container_t::make_link_restraints_by_linear(const coot::protein_geometry &geom,
							     bool do_rama_plot_restraints) {


   // Last time (for monomer geometry), we got a residue and added
   // bonds, angles and torsions by checking the atoms of that single
   // residue.
   //
   // This time, we need 2 consecutive residues, checking the atom
   // types of each - note that the atoms have to correspond to the
   // correct comp_id for that atom.
   //
   int selHnd = mol->NewSelection();
   PPCResidue     SelResidue;
   int nSelResidues;


//    timeval start_time;
//    timeval current_time;
//    double td;


   mol->Select ( selHnd,STYPE_RESIDUE, 1, // .. TYPE, iModel
		 chain_id_save.c_str(), // Chain(s)
		 istart_res, "*",  // starting res
		 iend_res,   "*",  // ending res
		 "*",  // residue name
		 "*",  // Residue must contain this atom name?
		 "*",  // Residue must contain this Element?
		 "*",  // altLocs
		 SKEY_NEW // selection key
		 );
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   std::cout << "INFO:: GetSelIndex (make_link_restraints) returned " << nSelResidues
	     << " residues (for link restraints) between (and including) residues "
	     << istart_res << " and " << iend_res << " of chain " << chain_id_save
	     << std::endl;
   
   coot::bonded_pair_container_t bonded_residue_pairs =
      bonded_residues_conventional(selHnd, geom);

   int iv = make_link_restraints_by_pairs(geom, bonded_residue_pairs, "Link");

   if (do_rama_plot_restraints)
      add_rama_links(selHnd, geom);

   mol->DeleteSelection(selHnd);
   return iv;
}


int
coot::restraints_container_t::make_link_restraints_from_res_vec(const coot::protein_geometry &geom,
								bool do_rama_plot_restraints) {

   // this determines the link type
   coot::bonded_pair_container_t bonded_residue_pairs = bonded_residues_from_res_vec(geom);
//     std::cout << "   DEBUG:: in make_link_restraints_from_res_vec() found "
//  	     << bonded_residue_pairs.size() << " bonded residues " << std::endl;
   int iv = make_link_restraints_by_pairs(geom, bonded_residue_pairs, "Link");

   if (do_rama_plot_restraints)
      add_rama_links_from_res_vec(geom);

   return iv;
}

int
coot::restraints_container_t::make_link_restraints_by_pairs(const coot::protein_geometry &geom,
							    const coot::bonded_pair_container_t &bonded_residue_pairs,
							    std::string link_flank_link_string) {

   int iret = 0;
   int n_link_bond_restr = 0;
   int n_link_angle_restr = 0;
   int n_link_torsion_restr = 0;
   int n_link_plane_restr = 0;

   for (unsigned int ibonded_residue=0;
	ibonded_residue<bonded_residue_pairs.size();
	ibonded_residue++) {
      
      // now select pairs of residues: but not for the last one, which
      // doesn't have a peptide link on its C atom.
      //
      // It has a modification there, I suppose, but I am ignoring
      // that for now.
      //
      //
      std::string link_type = bonded_residue_pairs[ibonded_residue].link_type;

      CResidue *sel_res_1 = bonded_residue_pairs[ibonded_residue].res_1;
      CResidue *sel_res_2 = bonded_residue_pairs[ibonded_residue].res_2;

      if (0) { 
	 std::cout << " ------- looking for link :" << link_type
		   << ": restraints etc. between residues " 
		   << sel_res_1->GetChainID() << " " << sel_res_1->seqNum << " - " 
		   << sel_res_2->GetChainID() << " " << sel_res_2->seqNum << std::endl;
      }
	 
// 	    gettimeofday(&start_time, NULL);

      // Are these residues neighbours?  We can add some sort of
      // ins code test here in the future.
      
      //if ( abs(sel_res_1->GetSeqNum() - sel_res_2->GetSeqNum()) <= 1
      // || abs(sel_res_1->index - sel_res_2->index) <= 1) {

      if (1) { 

	 bool is_fixed_first_residue = bonded_residue_pairs[ibonded_residue].is_fixed_first;
	 bool is_fixed_second_residue = bonded_residue_pairs[ibonded_residue].is_fixed_second;

	 // link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);

	 n_link_bond_restr += add_link_bond(link_type,
					    sel_res_1, sel_res_2,
					    is_fixed_first_residue,
					    is_fixed_second_residue,
					    geom);
	 // 	    gettimeofday(&current_time, NULL);
	 // td = time_diff(current_time, start_time);
	 // t0 = td;

	 n_link_angle_restr += add_link_angle(link_type,
					      sel_res_1, sel_res_2,
					      is_fixed_first_residue,
					      is_fixed_second_residue,
					      geom);
	 // 	    gettimeofday(&current_time, NULL);
	 // td = time_diff(current_time, start_time);
	 // t1 = td;

	 n_link_plane_restr += add_link_plane(link_type,
					      sel_res_1, sel_res_2,
					      is_fixed_first_residue,
					      is_fixed_second_residue,
					      geom);
	 // 	    gettimeofday(&current_time, NULL);
	 // td = time_diff(current_time, start_time);
	 // t2 = td;

	 // link_torsions_type_flag is the type of phi/psi restraints
	 // std::cout << "---------------link_torsions_type_flag: "
	 // << link_torsions_type_flag
	 // << std::endl;

	 // gettimeofday(&current_time, NULL);
	 // td = time_diff(current_time, start_time);
	 // t3 = td;

	 // 	    std::cout << "after bonds    " << t0 << std::endl;
	 // 	    std::cout << "after angles   " << t1 << std::endl;
	 // 	    std::cout << "after planes   " << t2 << std::endl;
	 // 	    std::cout << "after torsions " << t3 << std::endl;
      }
   }

   std::cout << link_flank_link_string << " restraints: " << std::endl;
   std::cout << "   " << n_link_bond_restr    << " bond    links" << std::endl;
   std::cout << "   " << n_link_angle_restr   << " angle   links" << std::endl;
   std::cout << "   " << n_link_plane_restr   << " plane   links" << std::endl;
   return iret; 
}

// return in millisecs
// double 
// coot::restraints_container_t::time_diff(const timeval &current, const timeval &start) const {

//    double d = current.tv_sec - start.tv_sec;
//    d *= 1000.0;
//    d += double(current.tv_usec - start.tv_usec)/1000.0;
//    return d;


// Need to add test that residues are linked with trans
//
std::vector<coot::rama_triple_t>
coot::restraints_container_t::make_rama_triples(int SelResHnd,
						const coot::protein_geometry &geom) const {
   std::vector<coot::rama_triple_t> v;
   PPCResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(SelResHnd, SelResidue, nSelResidues);

   for (int i=0; i<(nSelResidues-2); i++) {
      if (SelResidue[i] && SelResidue[i+1] && SelResidue[i+2]) {
	 coot::rama_triple_t t(SelResidue[i], SelResidue[i+1], SelResidue[i+2]);
	 v.push_back(t);
      }
   }
   return v;
}

void
coot::restraints_container_t::add_rama_links(int selHnd, const coot::protein_geometry &geom) {
   // 
   // Note do_rama_plot_restraints is true before we come here.
   
   int n_link_torsion_restr = 0;
   std::vector<coot::rama_triple_t> v = make_rama_triples(selHnd, geom);
   for (unsigned int ir=0; ir<v.size(); ir++) {
      std::string link_type = "TRANS";
      add_rama(link_type,
	       v[ir].r_1, v[ir].r_2, v[ir].r_3,
		  0,0,0, geom); // no residues are fixed.
      n_link_torsion_restr++;
   }
   std::cout << "   " << n_link_torsion_restr << " torsion/rama links" << std::endl;
}

void
coot::restraints_container_t::add_rama_links_from_res_vec(const coot::protein_geometry &geom) {

   std::cout << "::::::::: FIX ME :::::: residue vec rama triples " << std::endl;
   // need to get rama triples
   //
   // make_rama_triples(, geom);
}

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_conventional(int selHnd,
					      const coot::protein_geometry &geom) const {

   float dist_crit = 3.0; // if atoms in different residues are closer
			  // than this, then they are considered
			  // bonded (potentially).
   
   // First add "linear" links
   // 
   coot::bonded_pair_container_t c = bonded_residues_by_linear(selHnd, geom);

   // Now, are there any other links?
   //
   PPCResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   if (nSelResidues > 1) {
      for (int ii=0; ii<nSelResidues; ii++) {
	 for (int jj=0; jj<nSelResidues; jj++) {
	    if (jj>ii) {
	       if (! c.linked_already_p(SelResidue[ii], SelResidue[jj])) {
		  std::pair<bool, float> d = closest_approach(SelResidue[ii], SelResidue[jj]);
		  if (d.first) {
		     if (d.second < dist_crit) {
			std::pair<std::string, bool> l =
			   find_link_type_rigourous(SelResidue[ii], SelResidue[jj], geom);
			if (l.first != "") {
			} 
		     } 
		  }
	       }
	    }
	 }
      }
   }
   return c;
}

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_by_linear(int SelResHnd,
							const coot::protein_geometry &geom) const {
   
   coot::bonded_pair_container_t c;
   PPCResidue SelResidue;
   int nSelResidues;
   mol->GetSelIndex(SelResHnd, SelResidue, nSelResidues);
   if (nSelResidues > 1) {
   
      std::string link_type("TRANS");
      if (coot::util::is_nucleotide_by_dict(SelResidue[0], geom))
	 link_type = "p"; // phosphodiester linkage

      for (int i=0; i<(nSelResidues-1); i++) {
	 if (SelResidue[i] && SelResidue[i+1]) {


	    // There are a couple of ways we allow links.  First, that
	    // the residue numbers are consecutive.
	    //
	    // If there is a gap in residue numbering, the C and N
	    // have to be within 3A.
	    
	    // Are these residues neighbours?  We can add some sort of
	    // ins code test here in the future.
	    if ( abs(SelResidue[i]->GetSeqNum() - SelResidue[i+1]->GetSeqNum()) <= 1) { 
	       link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);
	       if (link_type != "") {
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  coot::bonded_pair_t p(SelResidue[i], SelResidue[i+1],
					whole_first_residue_is_fixed,
					whole_second_residue_is_fixed, link_type);
		  c.try_add(p);
	       }
	    }

	    // distance check this one.... it could be the opposite of
	    // an insertion code, or simply a gap - and we don't want
	    // to make a bond for a gap.
	    // 
	    if (abs(SelResidue[i]->index - SelResidue[i+1]->index) <= 1) {
	       // link_type = find_link_type(SelResidue[i], SelResidue[i+1], geom);
	       std::pair<std::string, bool> link_info =
		  find_link_type_rigourous(SelResidue[i], SelResidue[i+1], geom);
	       if (link_info.first != "") {
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  if (link_info.second == 0) { 
		     coot::bonded_pair_t p(SelResidue[i], SelResidue[i+1],
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_type);
		     c.try_add(p);
		  } else {
		     // order switch
		     coot::bonded_pair_t p(SelResidue[i+1], SelResidue[i],
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_type);
		     c.try_add(p);
		  }
	       } 
	    }
	 }
      }
   }
   return c;
} 

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_residues_from_res_vec(const coot::protein_geometry &geom) const {

   bool debug = 1;
   coot::bonded_pair_container_t bpc;
   float dist_crit = 3.0;
   // std::cout << "  debug:: residues_vec.size() " << residues_vec.size() << std::endl;

   int nres = residues_vec.size();
   for (int ii=0; ii<residues_vec.size(); ii++) {
      CResidue *res_f = residues_vec[ii].second;
      for (int jj=ii+1; jj<residues_vec.size(); jj++) {
	 CResidue *res_s = residues_vec[jj].second;
	 std::pair<bool, float> d = closest_approach(res_f, res_s);
// 	 std::cout << " closest approach given " << coot::residue_spec_t(res_f)
// 		   << " and " << coot::residue_spec_t(res_s) << std::endl;
// 	 std::cout << " closest approach d " << d.first << " " << d.second << std::endl;
	 if (d.first) {
	    if (d.second < dist_crit) {
	       std::pair<std::string, bool> l = find_link_type_rigourous(res_f, res_s, geom);
	       std::string link_type = l.first;
	       if (link_type != "") {

		  // too verbose?
		  if (1) 
		     std::cout << "   INFO:: "
			       << coot::residue_spec_t(res_f) << " " << coot::residue_spec_t(res_s)
			       << " link_type :" << link_type << ":" << std::endl;
		  bool whole_first_residue_is_fixed = 0;
		  bool whole_second_residue_is_fixed = 0;
		  bool order_switch_flag = l.second;
		  if (!order_switch_flag) { 
		     coot::bonded_pair_t p(res_f, res_s,
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_type);
		     bpc.try_add(p);
		  } else {
		     coot::bonded_pair_t p(res_s, res_f,
					   whole_first_residue_is_fixed,
					   whole_second_residue_is_fixed, link_type);
		     bpc.try_add(p);
		  }
	       } 
	    }
	 }
      }
   }
   return bpc;
}

bool
coot::bonded_pair_container_t::linked_already_p(CResidue *r1, CResidue *r2) const {

   bool r = 0;
   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      if (((bonded_residues[i].res_1 == r1) &&
	   (bonded_residues[i].res_2 == r2)) ||
	  ((bonded_residues[i].res_1 == r2) &&
	   (bonded_residues[i].res_2 == r1))) {
	 r = 1;
	 break;
      }
   }
   return r;
}

bool
coot::bonded_pair_container_t::try_add(const coot::bonded_pair_t &bp) {

   bool found = 0;
   for (unsigned int i=0; i<bonded_residues.size(); i++) {
      if ( (bonded_residues[i].res_1 == bp.res_1 &&
	    bonded_residues[i].res_2 == bp.res_2) ||
	   (bonded_residues[i].res_1 == bp.res_2 &&
	    bonded_residues[i].res_2 == bp.res_1) ) {
	 found = 1;
	 break;
      }
   }
   
   if (! found) {
      bonded_residues.push_back(bp);
   }
   return found;
}

std::ostream&
coot::operator<<(std::ostream &s, coot::bonded_pair_container_t bpc) {

   s << "Bonded Pair Container contains " << bpc.bonded_residues.size() << " bonded residues"
     << "\n";

   for (unsigned int i=0; i<bpc.bonded_residues.size(); i++)
      s << "   " << i << "  [:"
	<< bpc[i].link_type << ": "
	<< bpc[i].res_1->GetChainID() << " " << bpc[i].res_1->GetSeqNum() << " "
	<< bpc[i].res_1->GetInsCode() << " to " << bpc[i].res_2->GetChainID() << " "
	<< bpc[i].res_2->GetSeqNum() << " " << bpc[i].res_2->GetInsCode() << "]"
	<< "\n";

   return s; 
}


// Simple, residue-name and restraints group type link finder
// 
std::string
coot::restraints_container_t::find_link_type(CResidue *first, CResidue *second,
					     const coot::protein_geometry &geom) const {

   // Should return TRANS, PTRANS (PRO-TRANS), CIS, PCIS, p, BETA1-2, BETA1-4 etc.
   std::string link_type(""); // unset

   std::string residue_type_1 = first->name;
   std::string residue_type_2 = second->name;
   if (residue_type_1 == "UNK") residue_type_1 = "ALA"; // hack for KDC.
   if (residue_type_2 == "UNK") residue_type_2 = "ALA";

   std::string t1="";
   std::string t2="";

   for (int idr=0; idr<geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_1)) {
	 t1 = geom[idr].residue_info.group;
	 break;
      }
   }
   for (int idr=0; idr<geom.size(); idr++) {
      if (dictionary_name_matches_coords_resname(geom.three_letter_code(idr), residue_type_2)) {
	 t2 = geom[idr].residue_info.group;
	 break;
      }
   }

   if (t1 == "L-peptide" || t1 == "D-peptide")
      if (t2 == "L-peptide" || t2 == "D-peptide")
	 link_type = "TRANS";
   
   if (coot::util::is_nucleotide_by_dict(first, geom))
      link_type = "p"; // phosphodiester linkage

   if (t1 == "D-pyranose" || t1 == "D-furanose" || t1 == "L-pyranose" || t1 == "L-furanose") { 
      if (t2 == "D-pyranose" || t2 == "D-furanose" || t2 == "L-pyranose" || t2 == "L-furanose") {
	 link_type = find_glycosidic_linkage_type(first, second, t1, t2, geom);
	 std::cout << "INFO:: glycosidic_linkage type :" << link_type << ":\n";
      } 
   }

   return link_type; 
}

// Return "" on failure to find link.  
// 
// If this returns "", consider calling this function again, with
// reversed arguments.
// 
std::string
coot::restraints_container_t::find_glycosidic_linkage_type(CResidue *first, CResidue *second,
							   const std::string &group1,
							   const std::string &group2,
							   const protein_geometry &geom) const {

   double critical_dist = 3.0; // A, less than that and Coot should
			       // try to make the bond.
   PPCAtom res_selection_1 = NULL;
   PPCAtom res_selection_2 = NULL;
   int i_no_res_atoms_1;
   int i_no_res_atoms_2;
   double d;
   std::vector<coot::glycosidic_distance> close;
   
   first->GetAtomTable( res_selection_1, i_no_res_atoms_1);
   second->GetAtomTable(res_selection_2, i_no_res_atoms_2);

   for (int i1=0; i1<i_no_res_atoms_1; i1++) { 
      clipper::Coord_orth a1(res_selection_1[i1]->x,
			     res_selection_1[i1]->y,
			     res_selection_1[i1]->z);
      for (int i2=0; i2<i_no_res_atoms_2; i2++) {
	 clipper::Coord_orth a2(res_selection_2[i2]->x,
				res_selection_2[i2]->y,
				res_selection_2[i2]->z);
	 d = (a1-a2).lengthsq();
	 if (d < critical_dist*critical_dist) {
	    close.push_back(coot::glycosidic_distance(res_selection_1[i1],
						      res_selection_2[i2],
						      sqrt(d)));
	 }
      }
   }

   std::sort(close.begin(), close.end());

   // if you consider to uncomment this to debug a repulsion instead
   // the forming of a glycosidic bond, consider the residue numbering
   // of the residues involved: the "residue 1" should have the O4 and
   // the "residue 2" (+1 residue number) should have the C1.
   // 
   if (0) { 
      std::cout << "DEBUG:: number of sorted distances in glycosidic_linkage: "
		<< close.size() << std::endl;
      for (unsigned int i=0; i<close.size(); i++) {
	 std::cout << "#### glyco close: " << close[i].distance << "  "
		   << close[i].at1->GetChainID() << " " 
		   << close[i].at1->GetSeqNum() << " " 
		   << close[i].at1->GetAtomName() << " " 
		   << " to "
		   << close[i].at2->GetChainID() << " " 
		   << close[i].at2->GetSeqNum() << " " 
		   << close[i].at2->GetAtomName() << " " 
		   << std::endl;
      }
   }

   std::string link_type("");
   
   // short int found_link = 0;
   float smallest_link_dist = 99999.9;
   for (unsigned int i=0; i<close.size(); i++) {
      std::string name_1(close[i].at1->name);
      std::string name_2(close[i].at2->name);
      if (name_1 == " O4 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "BETA1-4";
	    }
      
      if (name_1 == " O2 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "BETA1-2";
	    }
      
      if (name_1 == " O3 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "BETA1-3";
	    }
      
      if (name_1 == " O3 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "BETA1-3";
	    }
	       
      if (name_1 == " C2 " )
	 if (name_2 == " O3 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "BETA2-3";
	    }
	       
      if (name_1 == " O6 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "BETA1-6";
	    }
	       
      if (name_1 == " O2 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "ALPHA1-2";
	    }
	       
      if (name_1 == " O3 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "ALPHA1-3";
	    }
      
      if (name_1 == " C2 " )
	 if (name_2 == " O3 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "ALPHA2-3";
	    }
      
      if (name_1 == " O4 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "ALPHA1-4";
	    }
      
      if (name_1 == " O6 " )
	 if (name_2 == " C1 ")
	    if (close[i].distance < smallest_link_dist) {
	       smallest_link_dist = close[i].distance;
	       link_type = "ALPHA1-6";
	    }
   }
   return link_type;
}

bool
coot::restraints_container_t::link_infos_are_glycosidic_p(const std::vector<std::pair<coot::chem_link, bool> > &link_infos) const {

   bool is_sweet = 0;
   for (unsigned int i=0; i<link_infos.size(); i++) {
      std::string id = link_infos[i].first.Id();
      if (id.length() > 4) {
	 if ((id.substr(0,5) == "ALPHA") ||
	     (id.substr(0,4) == "BETA")) {
	    is_sweet = 1;
	    break;
	 } 
      } 
   } 
   return is_sweet;
} 

// Return the link type and a residue order switch flag.
// Return link_type as "" if not found.
// 
std::pair<std::string, bool>
coot::restraints_container_t::find_link_type_rigourous(CResidue *first, CResidue *second,
						       const coot::protein_geometry &geom) const {

   bool debug = 0;
   std::string link_type = "";
   bool order_switch_flag = 0;
   std::string comp_id_1 = first->GetResName();
   std::string comp_id_2 = second->GetResName();

   try {
      std::string group_1 = geom.get_group(first);
      std::string group_2 = geom.get_group(second);
      if (debug)
	 std::cout << "====== debug:: find_link_type_rigourous() from "
		   << first->GetChainID() << " " << first->GetSeqNum() << " " << first->GetResName()
		   << " <--> " 
		   << second->GetChainID() << " " << second->GetSeqNum() << " " << second->GetResName()
		   << " " << std::endl;
      try {
	 std::vector<std::pair<coot::chem_link, bool> > link_infos =
	    geom.matching_chem_link(comp_id_1, group_1, comp_id_2, group_2);

	 if (debug) { 
	    std::cout << "     DEBUG:: find_link_type_rigourous: "
		      << link_infos.size() 
		      << " possible links: (link_infos):\n";
	    for (unsigned int il=0; il<link_infos.size(); il++)
	       std::cout << "            possible links: (link_infos): "
			 << il << " " << link_infos[il].first << " "
			 << link_infos[il].second << std::endl;
	 }

	 // Now, if link is a TRANS (default-peptide-link), then
	 // make sure that the C and N (or N and C) atoms of the
	 // first and second residue are within dist_crit (2.0A) of
	 // each other.  If not, then we should fail to find a link
	 // between these 2 residues.
	 // 
	 // unsigned int ilink = 0; // the first link

	 for (unsigned int ilink=0; ilink<link_infos.size(); ilink++) { 

	    if (link_infos[ilink].first.is_peptide_link_p()) {
	       std::pair<bool, bool> close_info = peptide_C_and_N_are_close_p(first, second);
	       if (close_info.first) {
		  order_switch_flag = close_info.second;
		  // link_type = link_infos[ilink].first.Id();
		  link_type = "TRANS"; // 200100415 for now, we force
				       // all peptide links to be
				       // TRANS.  (We don't yet (as of
				       // today) know if this link was
				       // CIS or TRANS). TRANS has
				       // 5-atom (plane3) plane
				       // restraints, CIS does not.
		  if (debug) 
		     std::cout << "   ==== TRANS or CIS pass NvsC dist test with order switch "
			       << order_switch_flag << std::endl;
	       } else {
		  // std::cout << "   TRANS or CIS FAIL NvsC dist test " << std::endl;
		  std::vector<std::pair<coot::chem_link, bool> > link_infos_non_peptide =
		     geom.matching_chem_link_non_peptide(comp_id_1, group_1, comp_id_2, group_2);

		  // debug::
		  if (debug)
		     for (unsigned int il=0; il<link_infos_non_peptide.size(); il++)
			std::cout << "   debug:: non-peptide link: " << link_infos_non_peptide[il].first.Id()
				  << std::endl;
	       
		  // 20100330 eh?  is something missing here?  What shall
		  // we do with link_infos_non_peptide?  did I mean that,
		  // now the peptide test has failed, perhaps we should
		  // distance check potential other (non-peptide) links?

		  // returns "" if no link found
		  // second is the order switch flag;
		  // 
		  std::pair<std::string, bool> non_peptide_close_link_info = 
		     general_link_find_close_link(link_infos_non_peptide, first, second, geom);
		  if (non_peptide_close_link_info.first != "") { 
		     link_type = non_peptide_close_link_info.first;
		     order_switch_flag = non_peptide_close_link_info.second;
		  }

	       } 
	    } else {

	       if (debug)
		  std::cout << "   find_link_type_rigourous() not a peptide link... " << std::endl;
	       //
	       order_switch_flag = link_infos[ilink].second;
	    
	       if (link_infos_are_glycosidic_p(link_infos)) {
		  link_type = find_glycosidic_linkage_type(first, second, group_1, group_2, geom);
		  if (link_type == "") {
		     link_type = find_glycosidic_linkage_type(second, first, group_2, group_1, geom);
		     if (link_type != "")
			order_switch_flag = 1;
		  }
	       } else {

		  if (debug)
		     std::cout << "   find_link_type_rigourous() not a glycosidic_linkage... " << std::endl;

		  std::vector<std::pair<coot::chem_link, bool> > link_infos_non_peptide =
		     geom.matching_chem_link_non_peptide(comp_id_1, group_1, comp_id_2, group_2);

		  // debug::
		  if (debug)
		     for (unsigned int il=0; il<link_infos_non_peptide.size(); il++)
			std::cout << "   debug:: " << il << " of " << link_infos_non_peptide.size()
				  << " non-peptide link: " << link_infos_non_peptide[il].first.Id()
				  << std::endl;
	       
		  // returns "" if no link found
		  // second is the order switch flag;
		  // 
		  std::pair<std::string, bool> non_peptide_close_link_info = 
		     general_link_find_close_link(link_infos_non_peptide, first, second, geom);
		  if (non_peptide_close_link_info.first != "") { 
		     link_type = non_peptide_close_link_info.first;
		     order_switch_flag = non_peptide_close_link_info.second;
		  }
	       }
	    }
	 } 
      }
      catch (std::runtime_error mess_in) {
	 if (debug) { 
	    // didn't find a chem_link for this pair, that's OK sometimes.
	    std::cout << "CAUGHT exception: " << mess_in.what() << std::endl;
	    // geom.print_chem_links();
	 }
      } 
   }
   catch (std::runtime_error mess) {
      // Failed to get group.  We don't want to hear about not getting
      // the group of thousands of HOHs.
      // std::cout << mess.what() << std::endl;
   }

   if (debug)
      std::cout << "   DEBUG:: given " << coot::residue_spec_t(first) << " and "
		<< coot::residue_spec_t(second)
		<< " find_link_type_rigourous returns link type :"
		<< link_type << ": and order-switch flag " << order_switch_flag << std::endl;
   return std::pair<std::string, bool> (link_type, order_switch_flag);
}

// a pair, first is if C and N are close and second if and order
// switch is needed to make it so.
std::pair<bool, bool>
coot::restraints_container_t::peptide_C_and_N_are_close_p(CResidue *r1, CResidue *r2) const {

   float dist_crit = 2.0; // 2.0 A for a peptide link - so that we
			  // don't find unintentional peptides - which
			  // would be a disaster.

   CAtom *at_c_1 = NULL;
   CAtom *at_n_1 = NULL;
   CAtom *at_c_2 = NULL;
   CAtom *at_n_2 = NULL;

   PPCAtom residue_atoms_1 = NULL;
   PPCAtom residue_atoms_2 = NULL;
   int n_residue_atoms_1;
   int n_residue_atoms_2;
   r1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   r2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

   for (int iat=0; iat<n_residue_atoms_1; iat++) {
      std::string atom_name(residue_atoms_1[iat]->name);
      if (atom_name == " C  ") {
	 at_c_1 = residue_atoms_1[iat];
      } 
      if (atom_name == " N  ") {
	 at_n_1 = residue_atoms_1[iat];
      } 
   }

   for (int iat=0; iat<n_residue_atoms_2; iat++) {
      std::string atom_name(residue_atoms_2[iat]->name);
      if (atom_name == " C  ") {
	 at_c_2 = residue_atoms_2[iat];
      } 
      if (atom_name == " N  ") {
	 at_n_2 = residue_atoms_2[iat];
      } 
   }

   if (at_c_1 && at_n_2) {
      clipper::Coord_orth c1(at_c_1->x, at_c_1->y, at_c_1->z);
      clipper::Coord_orth n2(at_n_2->x, at_n_2->y, at_n_2->z);
      float d = clipper::Coord_orth::length(c1, n2);
      // std::cout << "   C1->N2 dist " << d << std::endl;
      if (d < dist_crit)
	 return std::pair<bool, bool> (1,0);
   } 

   if (at_n_1 && at_c_2) {
      clipper::Coord_orth n1(at_n_1->x, at_n_1->y, at_n_1->z);
      clipper::Coord_orth c2(at_c_2->x, at_c_2->y, at_c_2->z);
      float d = clipper::Coord_orth::length(n1, c2);
      // std::cout << "   N1->C2 dist " << d << std::endl;
      if (d < dist_crit)
	 return std::pair<bool, bool> (1,1);
   } 

   return std::pair<bool, bool> (0, 0);

}

// a pair, first is if C and N are close and second is if an order
// switch is needed to make it so.
//
// return "" as first if no close link found.
//
// alt confs are ignored when finding the atoms for the potential
// link.  Perhaps they should not be (but the calling function does
// not have them).
// 
std::pair<std::string, bool>
coot::restraints_container_t::general_link_find_close_link(std::vector<std::pair<coot::chem_link, bool> > &li,
							   CResidue *r1, CResidue *r2,
							   const coot::protein_geometry &geom) const {

   
   std::pair<std::string, bool> r("", 0);
   std::string rs = general_link_find_close_link_inner(li, r1, r2, geom);
   if (rs != "") {
      r.first = rs;
   } else { 
      rs = general_link_find_close_link_inner(li, r2, r1, geom);
      if (rs  != "") {
	 r.first = rs;
	 r.second = 1;
      }
   }

   return r;

}

std::string
coot::restraints_container_t::general_link_find_close_link_inner(std::vector<std::pair<coot::chem_link, bool> > &li,
								 CResidue *r1, CResidue *r2,
								 const coot::protein_geometry &geom) const {

   float dist_crit = 3.0; // Angstroms.
   bool debug = 0;
   
   std::string r("");
   std::pair<bool,float> close = closest_approach(r1, r2);
   if (close.first) {
      if (close.second < dist_crit) {
	 for (unsigned int ilink=0; ilink<li.size(); ilink++) {
	    coot::chem_link link = li[ilink].first;
	    // now find link with id in  dict_link_res_restraints;
	    dictionary_residue_link_restraints_t lr = geom.link(link.Id());
	    if (lr.link_id != "") { // found link with lind_id link.Id() { 
	       for (unsigned int ib=0; ib<lr.link_bond_restraint.size(); ib++) {
		  std::string atom_id_1 = lr.link_bond_restraint[ib].atom_id_1_4c();
		  std::string atom_id_2 = lr.link_bond_restraint[ib].atom_id_2_4c();
		  CAtom *at_1 = r1->GetAtom(atom_id_1.c_str());
		  CAtom *at_2 = r2->GetAtom(atom_id_2.c_str());
		  if (at_1 && at_2) {
		     clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		     clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		     double d = clipper::Coord_orth::length(p1,p2);
		     if (debug) 
			std::cout << "  dist check " << " link-bond-number: "
				  << ib << " " << coot::atom_spec_t(at_1) << " -to- "
				  <<  coot::atom_spec_t(at_2) << "   " << d << " for link type "
				  <<  lr.link_id << std::endl;
		     if (d < dist_crit) {
			r = link.Id();
			break;
		     } 
		  } else {
		     std::cout << "WARNING:: dictionary vs model problem? "
			       << "We didn't find the bonded linking atoms :"
			       << lr.link_bond_restraint[ib].atom_id_1_4c() << ": :"
			       << lr.link_bond_restraint[ib].atom_id_2_4c() << ":" << std::endl;
		  } 
	       }
	    }
	    if (r != "")
	       break;
	 }
      } else {
	 std::cout << "FAIL close enough with closer than dist_crit " << close.second << std::endl;
      }
   }
   return r;
}



std::pair<bool,float>
coot::restraints_container_t::closest_approach(CResidue *r1, CResidue *r2) const {

   return coot::closest_approach(mol, r1, r2);
} 


int
coot::restraints_container_t::make_flanking_atoms_restraints(const coot::protein_geometry &geom,
							     bool do_rama_plot_restraints) {

   coot::bonded_pair_container_t bonded_residue_pairs = bonded_flanking_residues(geom);
   int iv = make_link_restraints_by_pairs(geom, bonded_residue_pairs, "Flanking residue");

   int n_rama_restraints = -1; // unset, don't output an info line if
			       // do_rama_plot_restraints is not set.
   if (do_rama_plot_restraints) {
      // e.g 1 free 2 free 3 flanking (fixed).
      n_rama_restraints = make_flanking_atoms_rama_restraints(geom);  // returns 0 or something.
   }
   return iv;
}


int coot::restraints_container_t::make_flanking_atoms_rama_restraints(const protein_geometry &geom) {
   int n_rama_restraints = 0;

   if (istart_minus_flag && iend_plus_flag) {  // have flanking residues

      std::vector<coot::ramachandran_restraint_flanking_residues_helper_t> vrrfr;
      coot::ramachandran_restraint_flanking_residues_helper_t rrfr_1;
      rrfr_1.resno_first = istart_res-1;
      rrfr_1.resno_third = istart_res+1;
      rrfr_1.is_fixed[0] = 1;
      if (istart_res == iend_res) // i.e. just one moving residue
	 rrfr_1.is_fixed[2] = 1;
      vrrfr.push_back(rrfr_1);

      // we don't want to add 2 sets of flanking ramas for when
      // refining just one residue (with 2 flanking residues)
      if (istart_res != iend_res) { 
	 coot::ramachandran_restraint_flanking_residues_helper_t rrfr_2;
	 rrfr_2.resno_first = iend_res-1;
	 rrfr_2.resno_third = iend_res+1;
	 rrfr_2.is_fixed[2] = 1;
	 vrrfr.push_back(rrfr_2);
      }

      for (unsigned int iround=0; iround<vrrfr.size(); iround++) { 
      
	 int selHnd = mol->NewSelection();
	 PPCResidue SelResidue = NULL;
	 int nSelResidues;
	 mol->Select (selHnd,STYPE_RESIDUE, 1, // .. TYPE, iModel
		      (char *) chain_id_save.c_str(), // Chain(s)
		      vrrfr[iround].resno_first,   "*",  // starting res
		      vrrfr[iround].resno_third,   "*",  // ending res
		      "*",  // residue name
		      "*",  // Residue must contain this atom name?
		      "*",  // Residue must contain this Element?
		      "*",  // altLocs
		      SKEY_NEW); // selection key 
	 mol->GetSelIndex ( selHnd, SelResidue,nSelResidues );
	 std::cout << "DEBUG:: GetSelIndex (make_flanking_atoms_rama_restraints) returned " 
		   << nSelResidues << " residues (for flanking rama restraints)" << std::endl;
      
	 if (nSelResidues == 3) {
	    // super careful would mean that we check the link type of
	    // both pairs before calling this function:

	    if (0) { // debugging fixed atoms
	       for (int i=0; i<3; i++)
		  std::cout << "   calling add_rama with resno "
			    << SelResidue[i]->GetSeqNum() << " Fixed: "
			    << vrrfr[iround].is_fixed[i] << std::endl;
	    }
	       
	    add_rama("TRANS",
		     SelResidue[0], SelResidue[1], SelResidue[2],
		     vrrfr[iround].is_fixed[0],
		     vrrfr[iround].is_fixed[1],
		     vrrfr[iround].is_fixed[2], geom);
	 }
      
	 mol->DeleteSelection(selHnd);
      }
   }
      
   return n_rama_restraints;
}

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_flanking_residues(const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t bpc;

   // residue n is at the end of the active selection.  What is
   // residue n+1 from the mol? We will make a bonded pair of residues
   // n and n+1, and make the is_fixed_residue true for the n+1
   // (flanking) residue (and False for residue n, obviously).
   //
   // We need to ignore (don't add) a [n,n+1] pair if residue n+1 is
   // in the vector of residues that we pass.
   // 
   // We need to do this for [n,n+1] pairs and [n,n-1] pairs.
   //
   // First, what residue n?  I suppose that that would be a member of
   // the residue vector passed to the constructor.  But what if we
   // are not using the residues vector constructor?
   // 

   if (from_residue_vector)
      bpc = bonded_flanking_residues_by_residue_vector(geom);
   else
      bpc = bonded_flanking_residues_by_linear(geom);
   
   return bpc;
}



coot::bonded_pair_container_t
coot::restraints_container_t::bonded_flanking_residues_by_linear(const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t bpc;
   std::string link_type = "TRANS";
   PPCResidue SelResidue = NULL;
   int nSelResidues;
   int selHnd = mol->NewSelection();
   mol->Select (selHnd,STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		istart_res-1,   "*",  // starting res
		istart_res,     "*",  // ending res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",  // altLocs
		SKEY_NEW); // selection key 
   mol->GetSelIndex (selHnd, SelResidue, nSelResidues);
   std::cout << "INFO:: GetSelIndex (make_flanking_atoms_restraints) returned " 
	     << nSelResidues << " residues (flanking restraints)" << std::endl;
   if (nSelResidues > 1) {
      link_type = find_link_type(SelResidue[0], SelResidue[1], geom);
      if (coot::util::is_nucleotide_by_dict(SelResidue[0], geom))
	 link_type = "p"; // phosphodiester linkage
   
      coot::bonded_pair_t bp(SelResidue[0], SelResidue[1], 1, 0, link_type);
      bpc.try_add(bp);
   }
   mol->DeleteSelection(selHnd);
   
   // And now again for the C-terminal flanking residue:
   // 
   selHnd = mol->NewSelection();
   mol->Select (selHnd,STYPE_RESIDUE, 1, // .. TYPE, iModel
		chain_id_save.c_str(), // Chain(s)
		iend_res,   "*",  // starting res
		iend_res+1, "*",  // ending res
		"*",  // residue name
		"*",  // Residue must contain this atom name?
		"*",  // Residue must contain this Element?
		"*",  // altLocs
		SKEY_NEW); // selection key 
   mol->GetSelIndex (selHnd, SelResidue, nSelResidues);
   std::cout << "INFO:: GetSelIndex (make_flanking_atoms_restraints) returned " 
	     << nSelResidues << " residues (flanking restraints)" << std::endl;
   if (nSelResidues > 1) {
      link_type = find_link_type(SelResidue[0], SelResidue[1], geom);
      if (coot::util::is_nucleotide_by_dict(SelResidue[0], geom))
	 link_type = "p"; // phosphodiester linkage
      coot::bonded_pair_t bp(SelResidue[0], SelResidue[1], 0, 1, link_type);
      bpc.try_add(bp);
   }
   mol->DeleteSelection(selHnd);

   // std::cout << "DEBUG:: bonded_flanking_residues_by_linear() reutrns " << bpc;
   return bpc;
}

coot::bonded_pair_container_t
coot::restraints_container_t::bonded_flanking_residues_by_residue_vector(const coot::protein_geometry &geom) const {

   coot::bonded_pair_container_t bpc;
   float dist_crit = 3.0;

   // Don't forget to consider the case were we refine one residue,
   // and that is the ASN for a glycosylation.  So the neighbouring
   // residues and the GLC/NAG (say) will be flanking residues.
   //
   // But that is a really hard thing - isn't it? To find a GLC that
   // is connected to this residue?

   
   for (unsigned int ir=0; ir<residues_vec.size(); ir++) {

      std::vector<CResidue *> neighbours = coot::residues_near_residue(residues_vec[ir].second,
								       mol, dist_crit);
//       std::cout << " DEBUG:: bonded_flanking_residues_by_residue_vector using reference res "
//  		<< ir << " of " << residues_vec.size() << " " 
//  		<< residues_vec[ir].second << " with " << neighbours.size() << " neighbours"
// 		<< std::endl;
      
      for (unsigned int ineighb=0; ineighb<neighbours.size(); ineighb++) {

	 bool found = 0;
	 for (unsigned int ires=0; ires<residues_vec.size(); ires++) {
	    // pointer comparison
	    if (neighbours[ineighb] == residues_vec[ires].second) {
	       found = 1;
	       break;
	    } 
	 }

	 if (! found) {
	    // OK, so this neighbour was not in the passed set of
	    // moving residues, it can be a flanking residue then...
	    std::pair<bool, float> d = closest_approach(neighbours[ineighb], residues_vec[ir].second);
// 	    std::cout << " not in residue vec... " << d.first << " " << d.second
// 		      << " for " << coot::residue_spec_t(neighbours[ineighb]) << " to "
// 		      << coot::residue_spec_t(residues_vec[ir].second) << std::endl;
	    if (d.first) {
	       if (d.second < dist_crit) {
		  std::pair<std::string, bool> l =
		     find_link_type_rigourous(neighbours[ineighb], residues_vec[ir].second, geom);
		  std::string link_type = l.first;
		  if (link_type != "") {
		     bool order_switch_flag = l.second;
		     if (! order_switch_flag) {
			coot::bonded_pair_t bp(neighbours[ineighb], residues_vec[ir].second, 1, 0, link_type);
			bpc.try_add(bp);
		     } else {
			coot::bonded_pair_t bp(residues_vec[ir].second, neighbours[ineighb], 0, 1, link_type);
			bpc.try_add(bp);
		     }
		  }
	       }
	    }
	 }
      } 
   } 
   return bpc;
}


// Atoms that are not involved in bonds or angles, but are in the
// residue selection should be at least 2.7A away from each other.
// 
// Here are my anti-bumping notes:
//
//
//     Anti-bumping restraints in regularization:
//
//     Considering totally screwed-up geometry: We should add a strong
//     repulsion for atoms that are not bonded so that they go away from
//     each other.  
//
//     Something like a triangle function between 0->2A and 0 beyond that.
//
//     Each atom has to check the distance to each other atom in the
//     selection: if they are not bonded, get a repulsion score for that
//     distance.
//
//     Derivative of that should be not too tricky, similar to bonds, but
//     not the same.
//
//     Instead of 500, use 10*matrix?  Doesn't really matter, I think.
//
//     Instead of using 2.0 as the critical distance, let's instead use
//     d_crit:
//
//     Infact,
//
//     f = 1000-d*(1000/d_crit)   for d<d_crit
//     df/dd = -1000/d_crit
//     df/dx = df/dd dd/dx
//           = -1000/d_crit
//
//     It's like bonds:
//     d = sqrt[ (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 ]
//     => dd/dx = (xi-xj)/d
//
//     So df/dx = -1000/d_crit * (xi-xj)/d
//
//     Need to keep a list of repulsing atom pairs so that we don't have
//     to calculate them once each for distortion_score and derivates..?
//


int 
coot::restraints_container_t::make_non_bonded_contact_restraints() { 

   construct_non_bonded_contact_list();
   
   // so now filtered_non_bonded_atom_indices is filled.
   //    std::vector<bool> fixed_atom_flags(2);  // 2 atoms in this restraint.


//   std::cout << "non-bonded list:" << std::endl;
//   std::cout << "--------------------------------------------------\n";
//   for (int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
//     std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
//     for (int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) { 
//       std::cout << filtered_non_bonded_atom_indices[i][j] << " ";
//     } 
//     std::cout << std::endl;
//   } 
//  std::cout << "--------------------------------------------------\n";

   int n_nbc_r = 0;
   for (unsigned int i=0; i<filtered_non_bonded_atom_indices.size(); i++) { 
      for (unsigned int j=0; j<filtered_non_bonded_atom_indices[i].size(); j++) {
// 	 fixed_atom_flags[0] = is_fixed_first;
// 	 fixed_atom_flags[1] = is_fixed_second;

	 if (0) { 
	    std::cout << "adding non-bonded contact restraint " 
		      << i << ":  " << filtered_non_bonded_atom_indices[i][j] << " ["
		      << atom[i]->GetSeqNum() << " " << atom[i]->name << "] to [" 
		      << atom[filtered_non_bonded_atom_indices[i][j]]->GetSeqNum() << " " 
		      << atom[filtered_non_bonded_atom_indices[i][j]]->name << "]" << std::endl;
	 }
	 add_non_bonded(i, filtered_non_bonded_atom_indices[i][j]);
	 n_nbc_r++;
      }
   }
   return n_nbc_r;
}

// fill the member data filtered_non_bonded_atom_indices
// 
void
coot::restraints_container_t::construct_non_bonded_contact_list() {

   if (from_residue_vector)
      construct_non_bonded_contact_list_by_res_vec();
   else 
      construct_non_bonded_contact_list_conventional();

}

   
void
coot::restraints_container_t::construct_non_bonded_contact_list_conventional() {

   // So first, I need a method/list to determine what is not-bonded
   // to what.
   // 

//    std::cout << "bonded list:" << std::endl;
//    std::cout << "--------------------------------------------------\n";
//    for (int i=0; i<bonded_atom_indices.size(); i++) { 
//       std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
//       for (int j=0; j<bonded_atom_indices[i].size(); j++) { 
// 	 std::cout << bonded_atom_indices[i][j] << " ";
//       } 
//       std::cout << std::endl;
//    } 
//    std::cout << "--------------------------------------------------\n";

  // Now we need to know which indices into the PPCAtom atoms are in
  // the moving set (rather than the flanking atoms).
  // 
  std::vector<std::vector<int> > non_bonded_atom_indices;
  //
  short int was_bonded_flag;
  // Note: bonded_atom_indices is sized to n_atoms in init_shared_post().
  non_bonded_atom_indices.resize(bonded_atom_indices.size());

  // Set up some PPCAtom things needed in the loop:
  PPCAtom res_selection_local;
  int n_res_atoms;
  PPCAtom res_selection_local_inner;
  int n_res_atoms_inner;
  int atom_index, atom_index_inner;
  int ierr;

  // iar = i active residue, nSelResidues_active is class variable

  // std::cout << "INFO:: There are " << nSelResidues_active << " active residues\n";
  for (int iar=0; iar<nSelResidues_active; iar++) { 

     SelResidue_active[iar]->GetAtomTable(res_selection_local, n_res_atoms);
     // std::cout << "There are " << n_res_atoms << " active atoms in this active residue\n";

     for (int iat=0; iat<n_res_atoms; iat++) { 

	// set atom_index
	ierr = res_selection_local[iat]->GetUDData(udd_atom_index_handle, atom_index);
	if (ierr != UDDATA_Ok) { 
	   std::cout << "ERROR:: in getting UDDATA res_selection_local, ierr=" 
		     << ierr << " "
		     << res_selection_local[iat]->GetSeqNum() << " " 
		     << res_selection_local[iat]->GetAtomName() << " \n";
	}
	
	short int matched_oxt = 0;
	if (have_oxt_flag) { 
	   if (std::string(res_selection_local[iat]->name) == " OXT") { 
	      matched_oxt = 1;
	   } else { 
	      matched_oxt = 0;
	   }
	}

	if (! matched_oxt) { 

	   // For each of the bonds of atom with atom_index index
	   // we need to check if bonded_atom_indices[atom_index][j]
	   // matches any atom index of the active atoms:

	   for (int jar=0; jar<nSelResidues_active; jar++) { 
	   
	      SelResidue_active[jar]->GetAtomTable(res_selection_local_inner, 
						   n_res_atoms_inner);
	   
	      for (int jat=0; jat<n_res_atoms_inner; jat++) { 
	      
		 // set atom_index_inner
		 ierr =  res_selection_local_inner[jat]->GetUDData(udd_atom_index_handle, 
								   atom_index_inner);

		 if (atom_index == atom_index_inner) { 
		    // std::cout << "skipping same index " << std::endl;
		 } else {
		    
// 		    std::cout << "DEBUG:: checking bond pair " << atom_index << " " 
// 			      << atom_index_inner << " " 
// 			      << atom[atom_index]->name << " " << atom[atom_index]->GetSeqNum() << "    " 
// 			      << atom[atom_index_inner]->name << " " << atom[atom_index_inner]->GetSeqNum()
//          		      << std::endl;
	      
		    was_bonded_flag = 0;

		    if (have_oxt_flag) 
		       if (! strcmp(res_selection_local_inner[jat]->name, " OXT")) // matched
			  matched_oxt = 1;

		    if (! matched_oxt) { 

		       for (unsigned int j=0; j<bonded_atom_indices[atom_index].size(); j++) { 
		 
			  if (bonded_atom_indices[atom_index][j] == atom_index_inner) { 
			     was_bonded_flag = 1;
			     break;
			  } 
		       }
		 
		       if (was_bonded_flag == 0) { 
			  non_bonded_atom_indices[atom_index].push_back(atom_index_inner);
		       }
		    }
		 }
	      }
	   }
	}
     }
  }

//   std::cout << "non-bonded list (unfiltered):" << std::endl;
//   std::cout << "--------------------------------------------------\n";
//   for (int i=0; i<non_bonded_atom_indices.size(); i++) { 
//      std::cout << i << "  " << atom[i]->GetSeqNum() << " " << atom[i]->name << " : "; 
//      for (int j=0; j<non_bonded_atom_indices[i].size(); j++) { 
// 	std::cout << non_bonded_atom_indices[i][j] << " ";
//      } 
//      std::cout << std::endl;
//   } 
//   std::cout << "--------------------------------------------------\n";

  filter_non_bonded_by_distance(non_bonded_atom_indices, 8.0);
}

void
coot::restraints_container_t::construct_non_bonded_contact_list_by_res_vec() {

   double dist_crit = 8.0;
   filtered_non_bonded_atom_indices.resize(bonded_atom_indices.size());

   for (int i=0; i<bonded_atom_indices.size(); i++) {
      for (int j=0; j<bonded_atom_indices.size(); j++) {
	 if (i != j) {
	    if (!coot::is_member_p(bonded_atom_indices[i], j)) {

	       // We don't want NCBs for fixed residues (i.e. the
	       // flanking residues), so if both atoms are in residues
	       // that are not in residue_vec, then we don't add a NCB
	       // for that atom pair.

	       // atom j is not bonded to atom i, is it close? (i.e. within dist_crit?)
	       clipper::Coord_orth pt1(atom[i]->x, atom[i]->y, atom[i]->z);
	       clipper::Coord_orth pt2(atom[j]->x, atom[j]->y, atom[j]->z);
	       double d = clipper::Coord_orth::length(pt1, pt2);
	       if (d < dist_crit) {
		  CResidue *r1 = atom[i]->residue;
		  CResidue *r2 = atom[j]->residue;
		  
		  if (is_a_moving_residue_p(r1) && is_a_moving_residue_p(r2)) { 
		     filtered_non_bonded_atom_indices[i].push_back(j);
		  }
	       } 
	    }
	 } 
      }
   } 
} 



void 
coot::restraints_container_t::filter_non_bonded_by_distance(const std::vector<std::vector<int> > &non_bonded_atom_indices, double dist) { 

   filtered_non_bonded_atom_indices.resize(non_bonded_atom_indices.size());

   CAtom *atom_1;
   CAtom *atom_2;
   double dist2;
   double dist_lim2 = dist*dist;
   int i_at_ind;
   
   for (unsigned int i=0; i<non_bonded_atom_indices.size(); i++) { 
      for (unsigned int j=0; j<non_bonded_atom_indices[i].size(); j++) {
	 
	 atom_1 = atom[i];
	 atom_2 = atom[non_bonded_atom_indices[i][j]];

// 	 dist2 = clipper::Coord_orth::lengthsq(clipper::Coord_orth(atom_1->x, atom_1->y, atom_1->z), 
// 					       clipper::Coord_orth(atom_2->x, atom_2->y, atom_2->z));

	 dist2 = (clipper::Coord_orth(atom_1->x, atom_1->y, atom_1->z) -
		  clipper::Coord_orth(atom_2->x, atom_2->y, atom_2->z)).lengthsq();

	 if (dist2 < dist_lim2) { 
// 	    std::cout << "accepting non-bonded contact between " << atom_1->GetSeqNum() 
// 		      << " " << atom_1->name << " and " << atom_2->GetSeqNum() 
// 		      << " " << atom_2->name  << "\n";
	    atom_2->GetUDData(udd_atom_index_handle, i_at_ind); // sets i_at_ind.
	    filtered_non_bonded_atom_indices[i].push_back(i_at_ind);
	 } else { 
// 	    std::cout << "          non-bonded contact between " << atom_1->GetSeqNum() 
// 		      << " " << atom_1->name  << " and " << atom_2->GetSeqNum() 
// 		      << " " << atom_2->name << " rejected by distance\n";
	 } 
      }
   }
}

bool
coot::restraints_container_t::is_a_moving_residue_p(CResidue *r) const {

   bool ret = 0;
   for (unsigned int i=0; i<residues_vec.size(); i++) {
      if (residues_vec[i].second == r) {
	 ret = 1;
	 break;
      } 
   }
   return ret;
}


int
coot::restraints_container_t::add_bonds(int idr, PPCAtom res_selection,
					int i_no_res_atoms,
					PCResidue SelRes,
					const coot::protein_geometry &geom) {

   int n_bond_restr = 0;
   int index1, index2;
   bool debug = 0;

   for (unsigned int ib=0; ib<geom[idr].bond_restraint.size(); ib++) {
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name1(res_selection[iat]->name);

	 if (debug)
	    std::cout << "comparing first (pdb) :" << pdb_atom_name1
		      << ": with (dict) :"
		      << geom[idr].bond_restraint[ib].atom_id_1_4c()
		      << ":" << std::endl; 

	 if (pdb_atom_name1 == geom[idr].bond_restraint[ib].atom_id_1_4c()) {
	    for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

	       std::string pdb_atom_name2(res_selection[iat2]->name);

	       if (debug)
		  std::cout << "comparing second (pdb) :" << pdb_atom_name2
			    << ": with (dict) :"
			    << geom[idr].bond_restraint[ib].atom_id_2_4c()
			    << ":" << std::endl;
	       
	       if (pdb_atom_name2 == geom[idr].bond_restraint[ib].atom_id_2_4c()) {

		  // check that the alt confs aren't different
		  std::string alt_1(res_selection[iat ]->altLoc);
		  std::string alt_2(res_selection[iat2]->altLoc);
		  if (alt_1 == "" || alt_2 == "" || alt_1 == alt_2) { 

		     if (debug) { 
			std::cout << "atom match 1 " << pdb_atom_name1;
			std::cout << " atom match 2 " << pdb_atom_name2
				  << std::endl;
		     }

		     // now we need the indices of
		     // pdb_atom_name1 and
		     // pdb_atom_name2 in asc.atom_selection:

		     //  		  int index1_old = get_asc_index(pdb_atom_name1,
		     //  					     SelRes->seqNum,
		     //  					     SelRes->GetChainID());
		     //  		  int index2_old = get_asc_index(pdb_atom_name2,
		     //  					     SelRes->seqNum,
		     //  					     SelRes->GetChainID());

		     res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
		     res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);

		     // set the UDD flag for this residue being bonded/angle with 
		     // the other
		  
		     bonded_atom_indices[index1].push_back(index2);
		     bonded_atom_indices[index2].push_back(index1);

		     // 		  std::cout << "add_bond: " << index1_old << " " << index1 << std::endl;
		     // 		  std::cout << "add_bond: " << index2_old << " " << index2 << std::endl;

		     // this needs to be fixed for fixed atom (rather
		     // than just knowing that these are not flanking
		     // atoms).
		     // 
		     std::vector<bool> fixed_flags = make_fixed_flags(index1, index2);
// 		     std::cout << "creating (monomer) bond restraint with fixed flags "
// 			       << fixed_flags[0] << " " << fixed_flags[1] << " "
// 			       << atom[index1]->GetSeqNum() << " "
// 			       << atom[index1]->name << " to "
// 			       << atom[index2]->GetSeqNum() << " "
// 			       << atom[index2]->name
// 			       << " restraint index " << n_bond_restr << "\n";

		     add(BOND_RESTRAINT, index1, index2,
			 fixed_flags,
			 geom[idr].bond_restraint[ib].dist(),
			 geom[idr].bond_restraint[ib].esd(),
			 1.2);  // junk value
		     n_bond_restr++;
		  }
	       }
	    }
	 }
      }
   }

   return n_bond_restr;

}

int
coot::restraints_container_t::add_angles(int idr, PPCAtom res_selection,
					 int i_no_res_atoms,
					 PCResidue SelRes,
					 const coot::protein_geometry &geom) {

   int n_angle_restr = 0;
   int index1, index2, index3;

//    std::cout << "There are " << geom[idr].angle_restraint.size()
// 	     << " angle restraints for this residue type" << std::endl; 

   for (unsigned int ib=0; ib<geom[idr].angle_restraint.size(); ib++) {
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name1(res_selection[iat]->name);

//  	 std::cout << "angle:  comparing :" << pdb_atom_name1 << ": with :"
//  		   << geom[idr].angle_restraint[ib].atom_id_1_4c()
//  		   << ":" << std::endl;
	 
	 if (pdb_atom_name1 == geom[idr].angle_restraint[ib].atom_id_1_4c()) {
	    for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

	       std::string pdb_atom_name2(res_selection[iat2]->name);
	       if (pdb_atom_name2 == geom[idr].angle_restraint[ib].atom_id_2_4c()) {
				    
// 		  std::cout << "angle: atom match 1 " << pdb_atom_name1;
// 		  std::cout << " atom match 2 " << pdb_atom_name2
// 			    << std::endl;

		  for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
		     
		     std::string pdb_atom_name3(res_selection[iat3]->name);
		     if (pdb_atom_name3 == geom[idr].angle_restraint[ib].atom_id_3_4c()) {
		  

			// now we need the indices of
			// pdb_atom_name1 and
			// pdb_atom_name2 in asc.atom_selection:

//  			int index1_old = get_asc_index(pdb_atom_name1,
//  						   SelRes->seqNum,
// 						   SelRes->GetChainID());
//  			int index2_old = get_asc_index(pdb_atom_name2,
//  						   SelRes->seqNum,
//  						   SelRes->GetChainID());
//  			int index3_old = get_asc_index(pdb_atom_name3,
//  						   SelRes->seqNum,
//  						   SelRes->GetChainID());

			std::string alt_1(res_selection[iat ]->altLoc);
			std::string alt_2(res_selection[iat2]->altLoc);
			std::string alt_3(res_selection[iat3]->altLoc);

			if (((alt_1 == alt_2) && (alt_1 == alt_3)) ||
			    ((alt_1 == ""   ) && (alt_2 == alt_3)) ||
			    ((alt_2 == ""   ) && (alt_1 == alt_3)) ||
			    ((alt_3 == ""   ) && (alt_1 == alt_2)))
			   {
			
			   res_selection[iat ]->GetUDData(udd_atom_index_handle, index1);
			   res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
			   res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);

			   // std::cout << "add_angles: " << index1_old << " " << index1 << std::endl;
			   // std::cout << "add_angles: " << index2_old << " " << index2 << std::endl;
			   // std::cout << "add_angles: " << index3_old << " " << index3 << std::endl;
			
			   // set the UDD flag for this residue being bonded/angle with 
			   // the other
			
			   bonded_atom_indices[index1].push_back(index3);
			   bonded_atom_indices[index3].push_back(index1);
		  
			   // this needs to be fixed for fixed atom (rather
			   // than just knowing that these are not flanking
			   // atoms).
			   // 
			   std::vector<bool> fixed_flag = make_fixed_flags(index1, index2, index3);

			   add(ANGLE_RESTRAINT, index1, index2, index3,
			       fixed_flag,
			       geom[idr].angle_restraint[ib].angle(),
			       geom[idr].angle_restraint[ib].esd(),
			       1.2);  // junk value
			   n_angle_restr++;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return n_angle_restr;
}

int
coot::restraints_container_t::add_torsions(int idr, PPCAtom res_selection,
					   int i_no_res_atoms,
					   PCResidue SelRes,
					   const coot::protein_geometry &geom) {

   int n_torsion_restr = 0; 

   for (unsigned int ib=0; ib<geom[idr].torsion_restraint.size(); ib++) {

      // Joel Bard fix: Don't add torsion restraints for torsion that
      // have either s.d. or period 0

      if (geom[idr].torsion_restraint[ib].periodicity() > 0) { // we had this test most inner
	 if (geom[idr].torsion_restraint[ib].esd() > 0.000001) { // new test
	 
	    // now find the atoms
	    for (int iat=0; iat<i_no_res_atoms; iat++) {
	       std::string pdb_atom_name1(res_selection[iat]->name);

	       if (pdb_atom_name1 == geom[idr].torsion_restraint[ib].atom_id_1_4c()) {
		  for (int iat2=0; iat2<i_no_res_atoms; iat2++) {

		     std::string pdb_atom_name2(res_selection[iat2]->name);
		     if (pdb_atom_name2 == geom[idr].torsion_restraint[ib].atom_id_2_4c()) {
				    
			// 		  std::cout << "atom match 1 " << pdb_atom_name1;
			// 		  std::cout << " atom match 2 " << pdb_atom_name2
			// 			    << std::endl;

			for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
		     
			   std::string pdb_atom_name3(res_selection[iat3]->name);
			   if (pdb_atom_name3 == geom[idr].torsion_restraint[ib].atom_id_3_4c()) {
		  
			      for (int iat4=0; iat4<i_no_res_atoms; iat4++) {
		     
				 std::string pdb_atom_name4(res_selection[iat4]->name);
				 if (pdb_atom_name4 == geom[idr].torsion_restraint[ib].atom_id_4_4c()) {
		  

				    // now we need the indices of
				    // pdb_atom_name1 and
				    // pdb_atom_name2 in asc.atom_selection:

				    int index1 = get_asc_index(res_selection[iat]->name,
							       res_selection[iat]->altLoc,
							       SelRes->seqNum,
							       SelRes->GetInsCode(),
							       SelRes->GetChainID());
				    int index2 = get_asc_index(res_selection[iat2]->name,
							       res_selection[iat2]->altLoc,
							       SelRes->seqNum,
							       SelRes->GetInsCode(),
							       SelRes->GetChainID());
				    int index3 = get_asc_index(res_selection[iat3]->name,
							       res_selection[iat3]->altLoc,
							       SelRes->seqNum,
							       SelRes->GetInsCode(),
							       SelRes->GetChainID());
				    int index4 = get_asc_index(res_selection[iat4]->name,
							       res_selection[iat4]->altLoc,
							       SelRes->seqNum,
							       SelRes->GetInsCode(),
							       SelRes->GetChainID());

				    std::vector<bool> fixed_flags =
				       make_fixed_flags(index1, index2, index3, index4);
				    add(TORSION_RESTRAINT, index1, index2, index3, index4,
					fixed_flags,
					geom[idr].torsion_restraint[ib].angle(),
					geom[idr].torsion_restraint[ib].esd(),
					1.2,  // junk value
					geom[idr].torsion_restraint[ib].periodicity());
				    n_torsion_restr++;
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
   }

   return n_torsion_restr;
}


int
coot::restraints_container_t::add_chirals(int idr, PPCAtom res_selection,
					  int i_no_res_atoms,
					  PCResidue SelRes,
					  const coot::protein_geometry &geom) { 

   int n_chiral_restr = 0;
   int index1, index2, index3, indexc;
   
   //   std::cout << "DEBUG:: trying to add chirals for this residue..." << std::endl;
   
   for (unsigned int ic=0; ic<geom[idr].chiral_restraint.size(); ic++) {
      // for now, let's just reject restraints that are a "both",
      // better would be to check the geometry and refine to the one
      // that is closest.
      if (!geom[idr].chiral_restraint[ic].is_a_both_restraint()) { 
	 for (int iat1=0; iat1<i_no_res_atoms; iat1++) {
	    std::string pdb_atom_name1(res_selection[iat1]->name);
	    if (pdb_atom_name1 == geom[idr].chiral_restraint[ic].atom_id_1_4c()) {
	       
	       for (int iat2=0; iat2<i_no_res_atoms; iat2++) {
		  std::string pdb_atom_name2(res_selection[iat2]->name);
		  if (pdb_atom_name2 == geom[idr].chiral_restraint[ic].atom_id_2_4c()) {
		     
		     for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
			std::string pdb_atom_name3(res_selection[iat3]->name);
			if (pdb_atom_name3 == geom[idr].chiral_restraint[ic].atom_id_3_4c()) {
			   
			   for (int iatc=0; iatc<i_no_res_atoms; iatc++) {
			      std::string pdb_atom_namec(res_selection[iatc]->name);
			      if (pdb_atom_namec == geom[idr].chiral_restraint[ic].atom_id_c_4c()) {
				 
//   			      std::cout << "DEBUG:: adding chiral number " << ic << " for " 
//   					<< res_selection[iatc]->GetSeqNum() << " "
//   					<< res_selection[0]->GetResName()
//   					<< pdb_atom_namec << " bonds to "
//   					<< pdb_atom_name1 << " "
//   					<< pdb_atom_name2 << " "
//   					<< pdb_atom_name3 << " "
//   					<< std::endl;

				 std::string alt_conf_c = res_selection[iatc]->altLoc;
				 std::string alt_conf_1 = res_selection[iat1]->altLoc;
				 std::string alt_conf_2 = res_selection[iat2]->altLoc;
				 std::string alt_conf_3 = res_selection[iat3]->altLoc;

				 if (((alt_conf_1 == alt_conf_c) || (alt_conf_1 == "")) &&
				     ((alt_conf_2 == alt_conf_c) || (alt_conf_2 == "")) && 
				     ((alt_conf_3 == alt_conf_c) || (alt_conf_3 == ""))) { 

				    res_selection[iat1]->GetUDData(udd_atom_index_handle, index1);
				    res_selection[iat2]->GetUDData(udd_atom_index_handle, index2);
				    res_selection[iat3]->GetUDData(udd_atom_index_handle, index3);
				    res_selection[iatc]->GetUDData(udd_atom_index_handle, indexc);

				    if (0) 
				       std::cout << "   Adding chiral restraint for " << res_selection[iatc]->name
						 << " " << res_selection[iatc]->GetSeqNum() <<  " "
						 << res_selection[iatc]->GetChainID() 
						 << " with target volume " << geom[idr].chiral_restraint[ic].target_volume()
						 << " with volume sigma " << geom[idr].chiral_restraint[ic].volume_sigma()
						 << " with volume sign " << geom[idr].chiral_restraint[ic].volume_sign
						 << " idr index: " << idr << " ic index: " << ic << std::endl;

				    std::vector<bool> fixed_flags =
				       make_fixed_flags(indexc, index1, index2, index3);
				    restraints_vec.push_back(simple_restraint(CHIRAL_VOLUME_RESTRAINT, indexc,
									      index1, index2, index3,
									      geom[idr].chiral_restraint[ic].volume_sign,
									      geom[idr].chiral_restraint[ic].target_volume(),
									      geom[idr].chiral_restraint[ic].volume_sigma(),
									      fixed_flags));
				    n_chiral_restr++;
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
   }
   return n_chiral_restr;
}


// Creates any number of simple_restraints for this monomer and adds
// them to restraints_vec.
// 
// idr provides the index of the comp_id (e.g. "ALA") match in geom.
// 
int
coot::restraints_container_t::add_planes(int idr, PPCAtom res_selection,
					 int i_no_res_atoms,
					 PCResidue SelRes,
					 const coot::protein_geometry &geom) {

//    std::cout << "There are " << geom[idr].plane_restraint.size()
// 	     << " plane restraints for " << SelRes->seqNum << " "
// 	     << geom[idr].comp_id << std::endl; 
   int n_plane_restr = 0;
   std::vector<std::string> altconfs;
   for (unsigned int ip=0; ip<geom[idr].plane_restraint.size(); ip++) {
      std::vector <int> pos; 
      for (int iat=0; iat<i_no_res_atoms; iat++) {
	 std::string pdb_atom_name(res_selection[iat]->name);
	 for (int irest_at=0; irest_at<geom[idr].plane_restraint[ip].n_atoms(); irest_at++) {
	    if (pdb_atom_name == geom[idr].plane_restraint[ip].atom_id(irest_at)) {
	       int idx = get_asc_index(res_selection[iat]->name,
				       res_selection[iat]->altLoc,
				       SelRes->seqNum,
				       SelRes->GetInsCode(),
				       SelRes->GetChainID());
	       if (idx >= 0) { 
		  pos.push_back(idx);
		  altconfs.push_back(res_selection[iat]->altLoc);
	       }
	    }
	 }
      }
      if (pos.size() > 3 ) {
	 if (check_altconfs_for_plane_restraint(altconfs)) { 
	 // Hoorah, sufficient number of plane restraint atoms found a match
	    //  	 std::cout << "adding a plane for " << SelRes->seqNum << " "
	    //  		   << geom[idr].comp_id << " esd: "
	    // 		   << geom[idr].plane_restraint[ip].dist_esd() << std::endl; 
	    // add_plane is in the .h file
	    std::vector<bool> fixed_flags = make_fixed_flags(pos);
// 	    std::cout << "DEBUG:: add_monomer plane restraint for plane-id: comp-id: "
// 		      << geom[idr].comp_id << " "
// 		      << geom[idr].plane_restraint[ip].plane_id << " "
// 		      << SelRes->GetChainID()
// 		      << " " << SelRes->GetSeqNum()
// 		      << " :" << SelRes->GetInsCode() << ":"
// 		      << std::endl;   
// 	    std::cout << "DEBUG:: adding monomer plane with pos indexes ";
// 	    for (int ipos=0; ipos<pos.size(); ipos++)
// 	       std::cout << " " << pos[ipos];
// 	    std::cout << "\n";
	    add_plane(pos, fixed_flags, geom[idr].plane_restraint[ip].dist_esd());
	    n_plane_restr++;
	 }
      } else {
	 if (verbose_geometry_reporting == 1) {  
	    std::cout << "in add_planes: sadly not every restraint atom had a match"
		      << std::endl;
	    std::cout << "   needed " << geom[idr].plane_restraint[ip].n_atoms()
		      << " got " << pos.size()<< std::endl;
	    std::cout << "   residue type: " << geom[idr].comp_id << " "
		      << "   plane id: " << ip << std::endl;
	 }
      } 
   }
   return n_plane_restr; 
}

// Add a trans bond linkage
//
// Residue 1 of the link is the first atom of the link
int
coot::restraints_container_t::add_link_bond(std::string link_type,
					    PCResidue first, PCResidue second,
					    short int is_fixed_first,
					    short int is_fixed_second,
					    const coot::protein_geometry &geom) {

   PPCAtom first_sel;
   PPCAtom second_sel;
   int n_first_res_atoms, n_second_res_atoms;
   std::vector<bool> fixed_atom_flags(2);  // 2 atoms in this restraint.

   first->GetAtomTable(first_sel,   n_first_res_atoms); 
   second->GetAtomTable(second_sel, n_second_res_atoms);
   short int found_link_type = 0;
   int index1, index2;

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

//    std::cout << "INFO:: geom.link_size() is " << geom.link_size() << std::endl;
//    std::cout << "first residue:\n";
//    for (int i=0; i<n_first_res_atoms; i++)
//       std::cout << "    " << first_sel[i]->name  << " " << first_sel[i]->GetSeqNum() << "\n";
//    std::cout << "second residue:\n";
//    for (int i=0; i<n_second_res_atoms; i++)
//       std::cout << "    " << second_sel[i]->name  << " " << second_sel[i]->GetSeqNum() << "\n";

   int nbond = 0;
   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) { // typically "TRANS"
	 found_link_type = 1;
	 for (unsigned int j=0; j<geom.link(i).link_bond_restraint.size(); j++) {
	    if (geom.link(i).link_bond_restraint[j].atom_1_comp_id == 1 && 
		geom.link(i).link_bond_restraint[j].atom_2_comp_id == 2) {
	       // as expected.
	    } else {
	       std::cout << "PROGRAMMER ERROR (shortsighted fool)" << std::endl;
	       std::cout << "bad things will now happen..." << std::endl; 
	    }
	    for (int ifat=0; ifat<n_first_res_atoms; ifat++) { 
	       std::string pdb_atom_name_1(first_sel[ifat]->name);

	       if (pdb_atom_name_1 == geom.link(i).link_bond_restraint[j].atom_id_1_4c()) {
		  for (int isat=0; isat<n_second_res_atoms; isat++) { 
		     std::string pdb_atom_name_2(second_sel[isat]->name);

		     if (pdb_atom_name_2 == geom.link(i).link_bond_restraint[j].atom_id_2_4c()) {

//   			std::cout << "DEBUG::  adding " << link_type << " bond for "
// 				  << first->seqNum
//    				  << " -> " << second->seqNum << " atoms " 
// 				  << first_sel [ifat]->GetAtomName() << " to "
// 				  << second_sel[isat]->GetAtomName() 
// 				  << std::endl;

			// Now, do the alt confs match?
			//
			std::string alt_conf_1 =  first_sel[ifat]->altLoc;
			std::string alt_conf_2 = second_sel[isat]->altLoc;
			if ((alt_conf_1 == alt_conf_2) || (alt_conf_1 == "") || (alt_conf_2 == "")) {

			   first_sel [ifat]->GetUDData(udd_atom_index_handle, index1);
			   second_sel[isat]->GetUDData(udd_atom_index_handle, index2);

			   // set the UDD flag for this residue being bonded/angle with 
			   // the other
		  
			   bonded_atom_indices[index1].push_back(index2);
			   bonded_atom_indices[index2].push_back(index1);

			   fixed_atom_flags[0] = is_fixed_first;
			   fixed_atom_flags[1] = is_fixed_second;

			   std::vector<bool> other_fixed_flags = make_fixed_flags(index1, index2);
			   for (int ii=0; ii<2; ii++)
			      if (other_fixed_flags[ii])
				 fixed_atom_flags[ii] = 1;

			   add(BOND_RESTRAINT, index1, index2,
			       fixed_atom_flags,
			       geom.link(i).link_bond_restraint[j].dist(),
			       geom.link(i).link_bond_restraint[j].esd(),
			       1.2); // junk value
			   nbond++;
			}
		     }
		  }
	       }
	    }
	 }
      } 
   } 

   if (found_link_type == 0)
      std::cout << "link type: " << link_type << " not found in dictionary!!\n";

   // std::cout << "add link bond returns " <<  nbond << std::endl;
   return nbond;
}


// Note that we have to convert between fixed residues and fixed
// atoms.  It was easier with bonds of course where there was a 1:1
// relationship.  Residue 1 of the link was the first atom of the link
int
coot::restraints_container_t::add_link_angle(std::string link_type,
					     PCResidue first, PCResidue second,
					     short int is_fixed_first,  // residue
					     short int is_fixed_second, // residue
					     const coot::protein_geometry &geom) {

   int nangle = 0;
   
   PPCAtom first_sel;
   PPCAtom second_sel;
   int n_first_res_atoms, n_second_res_atoms;
   int n_atom_1, n_atom_2, n_atom_3;
   int index1, index2, index3;

   first->GetAtomTable(first_sel,   n_first_res_atoms); 
   second->GetAtomTable(second_sel, n_second_res_atoms);

   PPCAtom atom_1_sel, atom_2_sel, atom_3_sel; // assigned to either
					       // first_sel or
					       // second_sel when
					       // atom_1_comp_id
					       // (etc.) are known.

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   std::vector<bool> fixed_flag(3);
   fixed_flag[0] = 0;  // not fixed
   fixed_flag[1] = 0;
   fixed_flag[2] = 0;

   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) { // typically TRANS
	 for (unsigned int j=0; j<geom.link(i).link_angle_restraint.size(); j++) {

	    if (geom.link(i).link_angle_restraint[j].atom_1_comp_id == 1) {
	       atom_1_sel = first_sel;
	       n_atom_1   = n_first_res_atoms;
	       fixed_flag[0] = is_fixed_first;
	    } else {
	       atom_1_sel = second_sel;
	       n_atom_1   = n_second_res_atoms;
	       fixed_flag[0] = is_fixed_second;
	    }
	    if (geom.link(i).link_angle_restraint[j].atom_2_comp_id == 1) {
	       atom_2_sel = first_sel;
	       n_atom_2   = n_first_res_atoms;
	       fixed_flag[1] = is_fixed_first;
	    } else {
	       atom_2_sel = second_sel;
	       n_atom_2   = n_second_res_atoms;
	       fixed_flag[1] = is_fixed_second;
	    }
	    if (geom.link(i).link_angle_restraint[j].atom_3_comp_id == 1) {
	       atom_3_sel = first_sel;
	       n_atom_3   = n_first_res_atoms;
	       fixed_flag[2] = is_fixed_first;
	    } else {
	       atom_3_sel = second_sel;
	       n_atom_3   = n_second_res_atoms;
	       fixed_flag[2] = is_fixed_second;
	    }
	    
	    for (int ifat=0; ifat<n_atom_1; ifat++) { 
	       std::string pdb_atom_name_1(atom_1_sel[ifat]->name);

	       if (pdb_atom_name_1 == geom.link(i).link_angle_restraint[j].atom_id_1_4c()) {
		  for (int isat=0; isat<n_atom_2; isat++) { 
		     std::string pdb_atom_name_2(atom_2_sel[isat]->name);

		     if (pdb_atom_name_2 == geom.link(i).link_angle_restraint[j].atom_id_2_4c()) {
			for (int itat=0; itat<n_atom_3; itat++) { 
			   std::string pdb_atom_name_3(atom_3_sel[itat]->name);
			   
			   if (pdb_atom_name_3 == geom.link(i).link_angle_restraint[j].atom_id_3_4c()) {

// 			      std::cout << "INFO: res "
// 					<< atom_1_sel[ifat]->residue->seqNum << " "
// 					<< atom_1_sel[ifat]->name << " to "
// 					<< atom_2_sel[isat]->residue->seqNum << " "
// 					<< atom_2_sel[isat]->name << " to "
// 					<< atom_3_sel[itat]->residue->seqNum << " "
// 					<< atom_3_sel[itat]->name << std::endl;
				 
//  			      int index1_old = get_asc_index(pdb_atom_name_1,
//  							 atom_1_sel[ifat]->residue->seqNum,
//  							 atom_1_sel[ifat]->residue->GetChainID());
			
//  			      int index2_old = get_asc_index(pdb_atom_name_2,
//  							 atom_2_sel[isat]->residue->seqNum,
//  							 atom_2_sel[isat]->residue->GetChainID());

//  			      int index3_old = get_asc_index(pdb_atom_name_3,
//  							 atom_3_sel[itat]->residue->seqNum,
//  							 atom_3_sel[itat]->residue->GetChainID());

			      std::string alt_conf_1 = atom_1_sel[ifat]->altLoc;
			      std::string alt_conf_2 = atom_2_sel[isat]->altLoc;
			      std::string alt_conf_3 = atom_3_sel[itat]->altLoc;
			      
			      if (((alt_conf_1 == alt_conf_2) && (alt_conf_1 == alt_conf_3)) ||
				  ((alt_conf_1 == alt_conf_2) && (alt_conf_3 == "")) ||
				  ((alt_conf_2 == alt_conf_3) && (alt_conf_1 == ""))) { 
			      
				 atom_1_sel[ifat]->GetUDData(udd_atom_index_handle, index1);
				 atom_2_sel[isat]->GetUDData(udd_atom_index_handle, index2);
				 atom_3_sel[itat]->GetUDData(udd_atom_index_handle, index3);
			      
				 if (0) { 
				    std::cout << "bonded_atom_indices.size(): "
					      <<  bonded_atom_indices.size() << std::endl;
				    std::cout << "   add_link_angle: "  << " " << index1 << std::endl;
				    std::cout << "   add_link_angle: "  << " " << index2 << std::endl;
				    std::cout << "   add_link_angle: "  << " " << index3 << std::endl;
				 }

				 // set the UDD flag for this residue being bonded/angle with 
				 // the other
			     
				 bonded_atom_indices[index1].push_back(index3);
				 bonded_atom_indices[index3].push_back(index1);

				 std::vector<bool> other_fixed_flags = make_fixed_flags(index1,
											index2,
											index3);
				 for (int ii=0; ii<other_fixed_flags.size(); ii++)
				    if (other_fixed_flags[ii])
				       fixed_flag[ii] = 1;

				 add(ANGLE_RESTRAINT, index1, index2, index3,
				     fixed_flag,
				     geom.link(i).link_angle_restraint[j].angle(),
				     geom.link(i).link_angle_restraint[j].angle_esd(),
				     1.2); // junk value
				 nangle++;
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
   return nangle;
}


int
coot::restraints_container_t::add_link_torsion(std::string link_type,
					       int phi_psi_restraints_type,
					       PCResidue first, PCResidue second,
					       short int is_fixed_first, short int is_fixed_second,
					       const coot::protein_geometry &geom) {

   // link_type is "p", "TRANS" etc.

  std::cout << "--------- :: Adding link torsion phi_psi_restraints_type: " 
	    << phi_psi_restraints_type << std::endl;
   
   int n_torsion = 0;
      
   PPCAtom first_sel;
   PPCAtom second_sel;
   int n_first_res_atoms, n_second_res_atoms;
   int n_atom_1, n_atom_2, n_atom_3, n_atom_4;

   
   first->GetAtomTable(first_sel,   n_first_res_atoms); 
   second->GetAtomTable(second_sel, n_second_res_atoms);

   // assigned to either first_sel or second_sel when atom_1_comp_id
   // (etc.) are known.
   PPCAtom atom_1_sel, atom_2_sel, atom_3_sel, atom_4_sel;

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   std::vector<bool> fixed_flag(4);
   fixed_flag[0] = 0;  // not fixed
   fixed_flag[1] = 0;
   fixed_flag[2] = 0;
   fixed_flag[3] = 0;

   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) {
	 for (unsigned int j=0; j<geom.link(i).link_torsion_restraint.size(); j++) {

	    // This could have been more compact if we were using a
	    // vector of atom ids... (or lisp)... heyho.
	    // 
	    if (geom.link(i).link_torsion_restraint[j].atom_1_comp_id == 1) {
	       atom_1_sel = first_sel;
	       n_atom_1 = n_first_res_atoms; 
	       fixed_flag[0] = is_fixed_first;
	    } else {
	       atom_1_sel = second_sel;
	       n_atom_1 = n_second_res_atoms; 
	       fixed_flag[0] = is_fixed_second;
	    }
	    if (geom.link(i).link_torsion_restraint[j].atom_2_comp_id == 1) {
	       atom_2_sel = first_sel;
	       n_atom_2 = n_first_res_atoms; 
	       fixed_flag[1] = is_fixed_first;
	    } else {
	       atom_2_sel = second_sel;
	       n_atom_2 = n_second_res_atoms; 
	       fixed_flag[1] = is_fixed_second;
	    }
	    if (geom.link(i).link_torsion_restraint[j].atom_3_comp_id == 1) {
	       atom_3_sel = first_sel;
	       n_atom_3 = n_first_res_atoms; 
	       fixed_flag[2] = is_fixed_first;
	    } else {
	       atom_3_sel = second_sel;
	       n_atom_3 = n_second_res_atoms; 
	       fixed_flag[2] = is_fixed_second;
	    }
	    if (geom.link(i).link_torsion_restraint[j].atom_4_comp_id == 1) {
	       n_atom_4 = n_first_res_atoms; 
	       atom_4_sel = first_sel;
	       fixed_flag[3] = is_fixed_first;
	    } else {
	       atom_4_sel = second_sel;
	       n_atom_4 = n_second_res_atoms; 
	       fixed_flag[3] = is_fixed_second;
	    }
	    for (int ifat=0; ifat<n_atom_1; ifat++) { 
	       std::string pdb_atom_name_1(atom_1_sel[ifat]->name);
	       
	       if (pdb_atom_name_1 == geom.link(i).link_torsion_restraint[j].atom_id_1_4c()) {
		  for (int isat=0; isat<n_atom_2; isat++) { 
		     std::string pdb_atom_name_2(atom_2_sel[isat]->name);
		     
		     if (pdb_atom_name_2 == geom.link(i).link_torsion_restraint[j].atom_id_2_4c()) {
			for (int itat=0; itat<n_atom_3; itat++) { 
			   std::string pdb_atom_name_3(atom_3_sel[itat]->name);
			   
			   if (pdb_atom_name_3 == geom.link(i).link_torsion_restraint[j].atom_id_3_4c()) {
			      for (int iffat=0; iffat<n_atom_4; iffat++) {
				 std::string pdb_atom_name_4(atom_4_sel[iffat]->name);
				 			   
				 if (pdb_atom_name_4 == geom.link(i).link_torsion_restraint[j].atom_id_4_4c()) {
				    
				    int index1 = get_asc_index(atom_1_sel[ifat]->name,
							       atom_1_sel[ifat]->altLoc,
							       atom_1_sel[ifat]->residue->seqNum,
							       atom_1_sel[ifat]->GetInsCode(),
							       atom_1_sel[ifat]->GetChainID());
			
				    int index2 = get_asc_index(atom_2_sel[isat]->name,
							       atom_2_sel[isat]->altLoc,
							       atom_2_sel[isat]->residue->seqNum,
							       atom_2_sel[ifat]->GetInsCode(),
							       atom_2_sel[isat]->GetChainID());
				    
				    int index3 = get_asc_index(atom_3_sel[itat]->name,
							       atom_3_sel[itat]->altLoc,
							       atom_3_sel[itat]->residue->seqNum,
							       atom_3_sel[itat]->GetInsCode(),
							       atom_3_sel[itat]->GetChainID());

				    int index4 = get_asc_index(atom_4_sel[iffat]->name,
							       atom_4_sel[iffat]->altLoc,
							       atom_4_sel[iffat]->residue->seqNum,
							       atom_4_sel[iffat]->GetInsCode(),
							       atom_4_sel[iffat]->residue->GetChainID());

// 				    std::cout << "making torsion " << geom.link(i).link_torsion_restraint[j].id()
// 					      << " from atoms \n    "
// 					      << atom_1_sel[ifat]->name << " " 
// 					      << atom_1_sel[ifat]->GetSeqNum() << "\n    " 
// 					      << atom_2_sel[isat]->name << " " 
// 					      << atom_2_sel[isat]->GetSeqNum() << "\n    " 
// 					      << atom_3_sel[itat]->name << " " 
// 					      << atom_3_sel[itat]->GetSeqNum() << "\n    " 
// 					      << atom_4_sel[iffat]->name << " " 
// 					      << atom_4_sel[iffat]->GetSeqNum() << "\n";
				       

				    double target_phi = -57.82 + 360.0;
				    double target_psi = -47.0  + 360.0;
				    double esd = 5.01;
				    if (phi_psi_restraints_type == coot::restraints_container_t::LINK_TORSION_ALPHA_HELIX) {
				       // theortical helix CA-N bond: -57.82 phi
				       //                  CA-C bond: -47    psi
				       target_phi = -57.82 + 360.0;
				       target_psi = -47.00 + 360.0;
				    }
				       
				    if (phi_psi_restraints_type == coot::restraints_container_t::LINK_TORSION_BETA_STRAND) {
				       // beta strand
				       target_phi = -110.0 + 360.0;
				       target_psi =  120.0; // approx values
				       esd = 15.0; // guess
				    }
				    
				    if (geom.link(i).link_torsion_restraint[j].id() == "phi") { 
				       add(TORSION_RESTRAINT, index1, index2, index3, index4,
					   fixed_flag,
					   target_phi,
					   esd,
					   1.2, // junk value
					   1);
				       n_torsion++;
				    }
				    if (geom.link(i).link_torsion_restraint[j].id() == "psi") { 
				       add(TORSION_RESTRAINT, index1, index2, index3, index4,
					   fixed_flag,
					   target_psi,
					   esd, 
					   1.2, // junk value (obs)
					   1);
				       n_torsion++;
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
      }
   }
   return n_torsion; 
}

// make RAMACHANDRAN_RESTRAINTs, not TORSION_RESTRAINTs these days.
int
coot::restraints_container_t::add_rama(std::string link_type,
				       PCResidue prev_res,
				       PCResidue this_res,
				       PCResidue post_res,
				       bool is_fixed_first,
				       bool is_fixed_second,
				       bool is_fixed_third,
				       const coot::protein_geometry &geom) {

   // Old notes:
   // TRANS    psi      1 N      1 CA     1 C      2 N   
   // TRANS    phi      1 C      2 N      2 CA     2 C   
   // TRANS    omega    1 CA     1 C      2 N      2 CA
   //
   // New assignements:
   // TRANS    psi    (2nd N) (2nd CA) (2nd C ) (3nd N)
   // TRANS    phi    (1st C) (2nd N ) (2nd CA) (2nd C) 
   // 
   // So Rama_atoms in this order:
   //   0       1        2      3         4
   // (1st C) (2nd N) (2nd CA) (2nd C) (3rd N)

   
   // std::cout << "DEBUG:: --------- :: Adding RAMA phi_psi_restraints_type" << std::endl;
   
   int n_rama = 0;
      
   PPCAtom prev_sel;
   PPCAtom this_sel;
   PPCAtom post_sel;
   int n_first_res_atoms, n_second_res_atoms, n_third_res_atoms;

   prev_res->GetAtomTable(prev_sel,  n_first_res_atoms); 
   this_res->GetAtomTable(this_sel, n_second_res_atoms);
   post_res->GetAtomTable(post_sel,  n_third_res_atoms);

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
      // throw 
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
      // throw
   }
   if (n_third_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
      // throw
   }

   std::vector<bool> fixed_flag(5, 0);
   if (is_fixed_first) {
      fixed_flag[0] = 1;
   }
   if (is_fixed_second) {
      fixed_flag[1] = 1;
      fixed_flag[2] = 1;
      fixed_flag[3] = 1;
   }
   if (is_fixed_third) {
      fixed_flag[4] = 1;
   }
   std::vector<CAtom *> rama_atoms(5);
   for (int ir=0; ir<5; ir++)
      rama_atoms[ir] = 0;
   
   for (int i=0; i<n_first_res_atoms; i++) {
      std::string atom_name(prev_sel[i]->name);
      if (atom_name == " C  ")
	 rama_atoms[0] = prev_sel[i];
   }
   for (int i=0; i<n_second_res_atoms; i++) {
      std::string atom_name(this_sel[i]->name);
      if (atom_name == " N  ")
	 rama_atoms[1] = this_sel[i];
      if (atom_name == " CA ")
	 rama_atoms[2] = this_sel[i];
      if (atom_name == " C  ")
	 rama_atoms[3] = this_sel[i];
   }
   for (int i=0; i<n_third_res_atoms; i++) {
      std::string atom_name(post_sel[i]->name);
      if (atom_name == " N  ")
	 rama_atoms[4] = post_sel[i];
   }

   if (rama_atoms[0] && rama_atoms[1] && rama_atoms[2] && 
       rama_atoms[3] && rama_atoms[4]) {

      std::vector<int> atom_indices(5, -1);
      for (int i=0; i<5; i++) {
	 atom_indices[i] = get_asc_index(rama_atoms[i]->name,
					 rama_atoms[i]->altLoc,
					 rama_atoms[i]->residue->seqNum,
					 rama_atoms[i]->GetInsCode(),
					 rama_atoms[i]->GetChainID());
      }

      if ( (atom_indices[0] != -1) && (atom_indices[1] != -1) && (atom_indices[2] != -1) && 
	   (atom_indices[3] != -1) && (atom_indices[4] != -1)) { 

	 // std::cout << "Adding RAMACHANDRAN_RESTRAINT " << std::endl;
	 
	 add(RAMACHANDRAN_RESTRAINT,
	     atom_indices[0], atom_indices[1], atom_indices[2],
	     atom_indices[3], atom_indices[4], fixed_flag);
	 n_rama++;
      }
   }
   // std::cout << "returning..." << n_rama << std::endl;
   return n_rama; 
}

int coot::restraints_container_t::add_link_plane(std::string link_type,
						 PCResidue first, PCResidue second,
						 short int is_fixed_first_res,
						 short int is_fixed_second_res,
						 const coot::protein_geometry &geom) {

//    std::cout << "DEBUG:: add_link_plane for " << first->GetChainID() << " " << first->GetSeqNum()
// 	     << " :" << first->GetInsCode() << ":"
// 	     << " -> " << second->GetChainID() << " " << second->GetSeqNum()
// 	     << " :" << second->GetInsCode() << ":" << std::endl;
   
   int n_plane = 0;

   PPCAtom first_sel;
   PPCAtom second_sel;
   PPCAtom atom_sel;  // gets assigned to either first or second
   int n_first_res_atoms, n_second_res_atoms;
   int link_res_n_atoms; // gets assigned to either n_first_res_atoms, n_second_res_atoms
   PCResidue res; 

   first->GetAtomTable(first_sel,   n_first_res_atoms); 
   second->GetAtomTable(second_sel, n_second_res_atoms);

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   // 20100424 there can be any number of plane restraints made here (not just one for each
   // link_plane_restraint) because there can be any number of alt confs passed in the
   // residues.  So we need to find all the alt conf in the 2 residues passed and make a
   // mapping between the alt conf and the plane restraint's vector of atom indices.
   //
   std::map<std::string, std::vector<int> > atom_indices_map;
   for (unsigned int iat=0; iat<n_first_res_atoms; iat++) {
      std::string alt_loc(first_sel[iat]->altLoc);
      std::map<std::string, std::vector<int> >::const_iterator it = atom_indices_map.find(alt_loc);
      if (it == atom_indices_map.end()) {
	 std::vector<int> v;
	 atom_indices_map[alt_loc] = v;
      }
   }
   for (unsigned int iat=0; iat<n_second_res_atoms; iat++) {
      std::string alt_loc(second_sel[iat]->altLoc);
      std::map<std::string, std::vector<int> >::const_iterator it = atom_indices_map.find(alt_loc);
      if (it == atom_indices_map.end()) {
	 std::vector<int> v;
	 atom_indices_map[alt_loc] = v;
      }
   }
   
   

   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) { // typically "TRANS"
	 for (unsigned int ip=0; ip<geom.link(i).link_plane_restraint.size(); ip++) {
	    std::vector<bool> fixed_flag(geom.link(i).link_plane_restraint[ip].n_atoms(), 0);

	    std::map<std::string, std::vector<int> >::iterator it;
	    for (it = atom_indices_map.begin(); it != atom_indices_map.end(); it++)
	       it->second.clear();

	    for (int irest_at=0; irest_at<geom.link(i).link_plane_restraint[ip].n_atoms(); irest_at++) {
	       if (geom.link(i).link_plane_restraint[ip].atom_comp_ids[irest_at] == 1) {
		  // std::cout << "comp_id: 1 " << std::endl;
		  link_res_n_atoms = n_first_res_atoms;
		  atom_sel = first_sel;
		  res = first;
		  fixed_flag[irest_at] = is_fixed_first_res;
	       } else {
		  // std::cout << "comp_id: 2 " << std::endl;
		  link_res_n_atoms = n_second_res_atoms;
		  atom_sel = second_sel;
		  res = second; 
		  fixed_flag[irest_at] = is_fixed_second_res;
	       }
	       for (int iat=0; iat<link_res_n_atoms; iat++) { 
		  std::string pdb_atom_name(atom_sel[iat]->name);
		  if (geom.link(i).link_plane_restraint[ip].atom_id(irest_at) == pdb_atom_name) {
		     if (0)
			std::cout << "     pushing back to :" << atom_sel[iat]->altLoc << ": vector "
				  << res->GetChainID() << " "
				  << res->seqNum << " :"
				  << atom_sel[iat]->name << ": :"
				  << atom_sel[iat]->altLoc << ":" << std::endl;

		     atom_indices_map[atom_sel[iat]->altLoc].push_back(get_asc_index(atom_sel[iat]->name,
										     atom_sel[iat]->altLoc,
										     res->seqNum,
										     res->GetInsCode(),
										     res->GetChainID()));
		  }
	       }
	    }


	    for (it = atom_indices_map.begin(); it != atom_indices_map.end(); it++) {
	       
	       if (it->second.size() > 3) {
		  
		  if (0) {  // debugging.
		     std::cout << "  atom indices: for plane restraint :" << it->first
			       << ":" << std::endl;
		     for (int indx=0; indx<it->second.size(); indx++) {
			std::cout << "           " << it->second[indx] << std::endl;
		     }
		     for (int ind=0; ind<it->second.size(); ind++) {
			std::cout << ind << " " << atom[it->second[ind]]->GetChainID() << " "
				  << atom[it->second[ind]]->GetSeqNum() << " :"
				  << atom[it->second[ind]]->name << ": :"
				  << atom[it->second[ind]]->altLoc << ":\n";
		     }
		     std::cout << "DEBUG:: adding link plane with pos indexes ";
		     for (int ipos=0; ipos<it->second.size(); ipos++)
			std::cout << " " << it->second[ipos];
		     std::cout << "\n";
		  }
		  
		  std::vector<bool> other_fixed_flags = make_fixed_flags(it->second);
		  for (int ii=0; ii<2; ii++)
		     if (other_fixed_flags[ii])
			fixed_flag[ii] = 1;
		  
		  add_plane(it->second, fixed_flag, geom.link(i).link_plane_restraint[ip].dist_esd());
		  n_plane++;
	       }
	    }
	 }
      }
   }
   return n_plane;
}

int coot::restraints_container_t::add_link_plane_tmp(std::string link_type,
						 PCResidue first, PCResidue second,
						 short int is_fixed_first_res,
						 short int is_fixed_second_res,
						 const coot::protein_geometry &geom) {

   return 0;

}



int
coot::restraints_container_t::get_asc_index(const char *at_name,
					    const char *alt_loc,
					    int resno,
					    const char *ins_code,
					    const char *chain_id) const {

   return get_asc_index_new(at_name, alt_loc, resno, ins_code, chain_id);
   
}

int
coot::restraints_container_t::get_asc_index_new(const char *at_name,
						const char *alt_loc,
						int resno,
						const char *ins_code,
						const char *chain_id) const {

   int index = -1;
   int SelHnd = mol->NewSelection();
   
   mol->SelectAtoms(SelHnd,
		    0,
		    chain_id,
		    resno, ins_code,
		    resno, ins_code,
		    "*",      // resnames
		    at_name,  // anames
		    "*",      // elements
		    alt_loc  // altLocs 
		    );

   int nSelAtoms;
   PPCAtom SelAtom = NULL;
   mol->GetSelIndex(SelHnd, SelAtom, nSelAtoms);

   if (nSelAtoms > 0) {
     if (udd_atom_index_handle >= 0) { 
       SelAtom[0]->GetUDData(udd_atom_index_handle, index); // sets index
     } else { 
       index = get_asc_index_old(at_name, resno, chain_id);
     } 
   }
   mol->DeleteSelection(SelHnd);

   return index;
}

int
coot::restraints_container_t::get_asc_index_old(const std::string &at_name,
						int resno,
						const char *chain_id) const {

   int index = -1;
   int SelHnd = mol->NewSelection();
   
   mol->SelectAtoms(SelHnd,
			0,
			(char *)chain_id,
			resno, "*",
			resno, "*",
			"*", // rnames
			(char *) at_name.c_str(), // anames
			"*", // elements
			"*" // altLocs 
			);

   int nSelAtoms;
   PPCAtom SelAtom;
   mol->GetSelIndex(SelHnd, SelAtom, nSelAtoms);

   if (nSelAtoms > 0) {
      // now check indices.
      // 
      // Sigh. Is mmdb really this shit or am I missing something?
      //
      // It's not shit, you are not using it as it is supposed to be
      // used, I think. Instead of passing around the index to an
      // atom selection, you should simply be passing a pointer to an
      // atom.
      for (int i=0; i<n_atoms; i++) {
	 if (atom[i] == SelAtom[0]) {
	    index = i;
	    break;
	 }
      }
   }
   mol->DeleteSelection(SelHnd);

   if (index == -1 ) { 
      std::cout << "ERROR:: failed to find atom index for "
		<< at_name << " " << resno << " " << chain_id
		<< std::endl;
   } 
   return index; 
}
					    

// now setup the gsl_vector with initial values
//
// We presume that the atoms in PPCAtom are exactly the same 
// order as they are in the pdb file that refmac/libcheck uses
// to generate the restraints. 
//  
void 
coot::restraints_container_t::setup_gsl_vector_variables() {

   int idx; 

   // recall that x is a class variable, 
   // (so are n_atoms and atom, which were set in the constructor)
   //  

   x = gsl_vector_alloc(3*n_atoms);

   // If atom is going out of date, check how atom is handled when the
   // restraints (this object) goes out of scope.  the destructor.
   
//    std::cout << "DEBUG:: using atom array pointer " << atom << std::endl;
//    std::cout << "DEBUG:: Top few atoms" << std::endl;
//    for (int i=0; i<10; i++) {
//       std::cout << "Top atom " << i << " "
// 		<< atom[i] << "   "
// 		<< atom[i]->x << " "
// 		<< atom[i]->y << " "
// 		<< atom[i]->z << " "
// 		<< std::endl;
//    } 

   for (int i=0; i<n_atoms; i++) {
      idx = 3*i; 
      gsl_vector_set(x, idx,   atom[i]->x);
      gsl_vector_set(x, idx+1, atom[i]->y);
      gsl_vector_set(x, idx+2, atom[i]->z);
   }

}


void 
coot::restraints_container_t::update_atoms(gsl_vector *s) { 

   int idx; 

   for (int i=0; i<n_atoms; i++) { 

      idx = 3*i; 
//       atom[i]->SetCoordinates(gsl_vector_get(s,idx), 
// 			      gsl_vector_get(s,idx+1),
// 			      gsl_vector_get(s,idx+2),
// 			      atom[i]->occupancy,
// 			      atom[i]->tempFactor+1); 
      atom[i]->x = gsl_vector_get(s,idx);
      atom[i]->y = gsl_vector_get(s,idx+1);
      atom[i]->z = gsl_vector_get(s,idx+2);
   }

   if (have_oxt_flag) { 
      position_OXT();
   } 
} 


void
coot::restraints_container_t::position_OXT() { 

   if (oxt_reference_atom_pos.size()== 4) { 
      // std::cout << "DEBUG:: Positioning OXT by dictionary" << std::endl;
      double tors_o = 
	 clipper::Coord_orth::torsion(oxt_reference_atom_pos[0], 
				      oxt_reference_atom_pos[1], 
				      oxt_reference_atom_pos[2], 
				      oxt_reference_atom_pos[3]);
      double angl_o = clipper::Util::d2rad(120.8);
      clipper::Coord_orth oxt_pos(oxt_reference_atom_pos[0], 
				  oxt_reference_atom_pos[1], 
				  oxt_reference_atom_pos[2], 
				  1.231, angl_o, tors_o + M_PI);
      atom[oxt_index]->x = oxt_pos.x();
      atom[oxt_index]->y = oxt_pos.y();
      atom[oxt_index]->z = oxt_pos.z();
   }
} 

int
coot::restraints_container_t::write_new_atoms(std::string pdb_file_name) { 

   //
   int status = -1;
   if (mol != NULL) {
      // return 0 on success, non-zero on failure.
      status = mol->WritePDBASCII((char *)pdb_file_name.c_str());
      if (status == 0)
	 std::cout << "INFO:: output file: " << pdb_file_name
		   << " written." << std::endl;
      else
	 std::cout << "WARNING:: output file: " << pdb_file_name
		   << " not written." << std::endl;
   } else { 
      cout << "not constructed from asc, not writing coords" << endl; 
   }
   return status;
}

void
coot::restraints_container_t::info() const {

   std::cout << "There are " << restraints_vec.size() << " restraints" << std::endl;

   for (unsigned int i=0; i< restraints_vec.size(); i++) {
      if (restraints_vec[i].restraint_type == coot::TORSION_RESTRAINT) {
	 std::cout << "restraint " << i << " is of type "
		   << restraints_vec[i].restraint_type << std::endl;

	 std::cout << restraints_vec[i].atom_index_1 << " "
		   << restraints_vec[i].atom_index_2 << " "
		   << restraints_vec[i].atom_index_3 << " "
		   << restraints_vec[i].atom_index_4 << " "
		   << restraints_vec[i].target_value << " "
		   << restraints_vec[i].sigma << " " << std::endl
		   << " with "
  		   << restraints_vec[i].atom_index.size() << " vector atoms " << std::endl
		   << " with periodicity "
  		   << restraints_vec[i].periodicity << std::endl;
      }

      std::cout << "restraint number " << i << " is restraint_type " <<
	 restraints_vec[i].restraint_type << std::endl;
   }
} 




#endif // HAVE_GSL
