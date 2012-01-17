/* src/graphics-info-superpose.cc
 * 
 * Copyright 2004, 2005 by The University of York
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
 * Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#if defined _MSC_VER
#include <windows.h>
#endif

#include <iostream>
#include "graphics-info.h"

// ----------------------------------------------------------------------------
//                              superposing
// ----------------------------------------------------------------------------
// 
// This is moved to graphics_info_t because the
// superpose_optionmenu_activate_mol1 were here (I don't recall why now...)
// 

// superposing:
// static 
void 
graphics_info_t::superpose_optionmenu_activate_mol1(GtkWidget *item, GtkPositionType pos) { 
   //    std::cout << "DEBUG:: superpose imol1 now " <<  pos << std::endl;
   graphics_info_t::superpose_imol1 = pos;
   GtkWidget *chain_optionmenu = lookup_widget(GTK_WIDGET(item), 
					       "superpose_reference_chain_optionmenu");
   GtkWidget *checkbutton = lookup_widget(GTK_WIDGET(item),
					  "superpose_reference_chain_checkbutton");
   if (GTK_TOGGLE_BUTTON(checkbutton)->active)
      graphics_info_t::fill_superpose_option_menu_with_chain_options(chain_optionmenu, 1);
} 

// static 
void
graphics_info_t::superpose_optionmenu_activate_mol2(GtkWidget *item, GtkPositionType pos) { 
   //    std::cout << "DEBUG:: superpose imol2 now " <<  pos << std::endl;
   graphics_info_t::superpose_imol2 = pos;
   GtkWidget *chain_optionmenu = lookup_widget(GTK_WIDGET(item), 
					       "superpose_moving_chain_optionmenu");
   GtkWidget *checkbutton = lookup_widget(GTK_WIDGET(item),
					  "superpose_moving_chain_checkbutton");
   if (GTK_TOGGLE_BUTTON(checkbutton)->active)
      graphics_info_t::fill_superpose_option_menu_with_chain_options(chain_optionmenu, 0);
} 
 

// static
void graphics_info_t::fill_superpose_option_menu_with_chain_options(GtkWidget *chain_optionmenu, 
								    int is_reference_structure_flag) {
   GtkWidget *mol_optionmenu = NULL;

   if (is_reference_structure_flag)
      mol_optionmenu = lookup_widget(chain_optionmenu,
				     "superpose_dialog_reference_mol_optionmenu");
   else 
      mol_optionmenu = lookup_widget(chain_optionmenu,
				     "superpose_dialog_moving_mol_optionmenu");

   GtkSignalFunc callback_func;
   int imol;
   if (is_reference_structure_flag) { 
      imol = graphics_info_t::superpose_imol1;
      callback_func =
	 GTK_SIGNAL_FUNC(graphics_info_t::superpose_reference_chain_option_menu_item_activate);
   } else {
      imol = graphics_info_t::superpose_imol2;
       callback_func =
	 GTK_SIGNAL_FUNC(graphics_info_t::superpose_moving_chain_option_menu_item_activate);
   }
   
   if (imol >=0 && imol < n_molecules()) { 
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_optionmenu,
								      imol, callback_func);
      if (is_reference_structure_flag) {
	 graphics_info_t::superpose_imol1_chain = set_chain;
      } else {
	 graphics_info_t::superpose_imol2_chain = set_chain;
      }
      
   } else {
      std::cout << "ERROR in imol in fill_superpose_option_menu_with_chain_options "
		<< std::endl;
   }
}

// static
void
graphics_info_t::superpose_reference_chain_option_menu_item_activate (GtkWidget *item,
								      GtkPositionType pos) {
   std::string s = graphics_info_t::menu_item_label(item);
   graphics_info_t::superpose_imol1_chain = s;
}



// This does a redraw
// 
// static
void
graphics_info_t::superpose_moving_chain_option_menu_item_activate (GtkWidget *item,
								   GtkPositionType pos) { 

   graphics_info_t::superpose_imol2_chain = menu_item_label(item);
}


// Note that the 1 and 2 meanings get swapped.  Here they mean: 1 for
// moving and 2 for reference (except imol2 (which is the number of the
// moving molecule)).
// 
//  For example, let's look at 1sar:A (reference) onto 1py3:A (moving)

// MSDfold gives:

// >q|PDB:1sar:A|PDB:1py3:A RIBONUCLEASE SA (E.C.3.1.4.8) 1SAR 3
// dvSGTVCLSALPPEATDTLNLIASDGPFPYSQDGVVFQNRESVLPTQSYGYYHEYTVITPGARTRGTRRI
// ICGEATqEDYYTGDHYATFSLIDQTC

// >t|PDB:1py3:A|PDB:1sar:A CRYSTAL STRUCTURE OF RIBONUCLEASE SA2
// -aLADVCRTKLPSQAQDTLALIAKNGPYPYNRDGVVFENRESRLPKKGNGYYHEFTVVTPGSNDRGTRRV
// VTGGYG-EQYWSPDHYATFQEIDPRC

// I do the same using SSMlib and write out the 
// SSMAlign->Ca1[i1] 0<=i1<nsel1
// and 
// SSMAlign->Ca2[i2] 0<=i2<nsel2

// alignment of 94 residues in moving molecule:
// -1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95

// alignment of 96 residues in reference molecule:
// -1 -1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 -1 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93

// How do these indexes convert to MSDfold alignment?

// OK, the numbers correspond to the residue index (selindex) of the *other* molecule
// 
int
graphics_info_t::superpose_with_atom_selection(atom_selection_container_t asc_ref,
					       atom_selection_container_t asc_mov,
					       int imol2,
					       std::string moving_mol_name,
					       std::string reference_mol_name,
					       short int move_copy_of_imol2_flag) {

   int imodel_return = -1;

#ifdef HAVE_SSMLIB

   int precision = SSMP_Normal;
   int connectivity = CSSC_Flexible;

   // probably not necessary, not sure:
   SetSSConnectivityCheck ( CSSC_Flexible );
   SetSSMatchPrecision    ( SSMP_Normal );

   if (asc_ref.n_selected_atoms > 0) { 
      if (asc_mov.n_selected_atoms > 0) {
	 

	 CMMDBManager *mol1 = asc_ref.mol;
	 CMMDBManager *mol2 = asc_mov.mol;

	 // debug atom selections:
	 PCAtom *debug_atom_selection1 = NULL;
	 PCAtom *debug_atom_selection2 = NULL;
	 int n_selected_atoms_1_debug, n_selected_atoms_2_debug;
	 asc_mov.mol->GetSelIndex(asc_mov.SelectionHandle,
				  debug_atom_selection1,
				  n_selected_atoms_1_debug);
	 asc_ref.mol->GetSelIndex(asc_ref.SelectionHandle,
				  debug_atom_selection2,
				  n_selected_atoms_2_debug);

	 std::string name1 = reference_mol_name;
	 std::string name2 = moving_mol_name;
	 std::cout << "superposing..." << std::endl;

	 // Remove the pointer one day.
	 // 
	 CSSMAlign *SSMAlign = new CSSMAlign();
	 int rc = SSMAlign->Align(mol2, mol1, precision, connectivity,
				  asc_mov.SelectionHandle,
				  asc_ref.SelectionHandle);

	 if (rc)  {
	    std::string ws;
	    switch (rc)  {
	    case SSM_noHits :
	       std::cout << " *** secondary structure does not match.\n";
	       ws = "secondary structure does not match";
	       break;
	    case SSM_noSPSN :
	       std::cout << " *** structures are too remote.\n";
	       ws = "structures are too remote";
	       break;
	    case SSM_noGraph :
	       std::cout << " *** can't make graph for " << name1 << "\n";
	       ws = "can't make graph for " + name1;
	       ws += " structure";
	       break;
	    case SSM_noVertices :
	       std::cout << " *** empty graph for " << name1 << "\n";
	       ws = "empty graph for " + name1;
	       break;
	    case SSM_noGraph2 :
	       std::cout << " *** can't make graph for " << name2 << "\n";
	       ws = "can't make graph for " + name2;
	       break;
	    case SSM_noVertices2 :
	       std::cout << " *** empty graph for " << name2 << "\n";
	       ws = "empty graph for " + name2;
	       break;
	    default :
	       std::cout << " *** undocumented return code: " << rc << "\n";
	    }
	    GtkWidget *w = wrapped_nothing_bad_dialog(ws);
	    gtk_widget_show(w);
	 } else  {


	    if (move_copy_of_imol2_flag) { 
	       mol2 = new CMMDBManager;
	       mol2->Copy(asc_mov.mol, MMDBFCM_All);
	       std::string name = "Copy_of_";
	       name += moving_mol_name;
	       int imol2_new = graphics_info_t::create_molecule();
	       graphics_info_t::molecules[imol2_new].install_model(imol2_new, make_asc(mol2), name, 1);
	       imol2 = imol2_new;
	    }
	    
	    imodel_return = imol2;
	    
	    // OK, let's get a consistent naming system:  1 is moving: 2 is reference

	    PCAtom *atom_selection1 = NULL;
	    PCAtom *atom_selection2 = NULL;
	    int n_selected_atoms_1, n_selected_atoms_2;
	    asc_mov.mol->GetSelIndex(SSMAlign->selHndCa1,
				     atom_selection1,
				     n_selected_atoms_1);
	    asc_ref.mol->GetSelIndex(SSMAlign->selHndCa2,
				     atom_selection2,
				     n_selected_atoms_2);

	    std::cout << "number of Ca atoms in selections: " << n_selected_atoms_1
		      << " (moving) and " << n_selected_atoms_2 << " (reference)"
		      << std::endl;
	    
	    // realtype *rv = SSMAlign->dist1;
	    std::string res_type1;
	    std::string res_type2;
	    std::string res_no_str_1;
	    std::string res_no_str_2;
	    std::string reference_seq, moving_seq;

	    // another go at the sequence alignment used in the superposition.
	    // 
	    // correct: number of residues in reference structure: " << SSMAlign->nres2

	    print_ssm_sequence_alignment(SSMAlign,
					 asc_ref, asc_mov,
					 atom_selection1, atom_selection2,
					 n_selected_atoms_1, n_selected_atoms_2,
					 move_copy_of_imol2_flag);

	    make_and_print_horizontal_ssm_sequence_alignment(SSMAlign, asc_ref, asc_mov,
							     atom_selection1, atom_selection2,
							     n_selected_atoms_1, n_selected_atoms_2);

	    graphics_info_t::molecules[imol2].transform_by(SSMAlign->TMatrix);
	    graphics_info_t::molecules[imol2].set_show_symmetry(0); 
	    // 
	    // ivector Ca1, Ca2;
	    // rvector dist1;
	    // int nCa1, nCa2;
	    // mat44 tmat;
	    // realtype rmsdAchieved, seqIdentity, nCombs;
	    // int nAligned, nGaps, nMisD;
	    // CSuperpose Superpose;  type in class CSSMAlign
	    // 
	    // Superpose->GetSuperposition(Ca1, dist1, nCa1, Ca2,
	    // nCa2, tmat, rmsdAchieved, nAligned, nGaps, seqIdentity,
	    // nMisD, nCombs);
	    //
	    std::cout << "INFO: core rmsd achieved: " << SSMAlign->rmsd << " Angstroems\n"
		      << "      number of residues in reference structure: " << SSMAlign->nres2 << "\n"
		      << "      number of residues in moving structure:    " << SSMAlign->nres1 << "\n"
		      << "      number of residues in aligned sections (reference):  " << SSMAlign->nsel2 << "\n"
		      << "      number of residues in aligned sections (moving):     " << SSMAlign->nsel1 << "\n"
		      << "      number of aligned residues:  " << SSMAlign->nalgn << "\n"
		      << "      number of gaps:              " << SSMAlign->ngaps << "\n"
		      << "      number of misdirections:     " << SSMAlign->nmd << "\n"
		      << "      number of SSE combinations:  " << SSMAlign->ncombs << "\n"
		      << "      sequence identity:           " << SSMAlign->seqIdentity*100.0 << "%\n";

	    graphics_draw();
	 }
	 delete SSMAlign;
      } else {
	 std::cout << "WARNING:: Molecule moving has no model atoms\n";
      }
   } else {
      std::cout << "WARNING:: Molecule reference has no model atoms\n";
   }
#endif // HAVE_SSMLIB
   return imodel_return;
}


#ifdef HAVE_SSMLIB
void
graphics_info_t::make_and_print_horizontal_ssm_sequence_alignment(CSSMAlign *SSMAlign,
							 atom_selection_container_t asc_ref,
							 atom_selection_container_t asc_mov,
							 PCAtom *atom_selection1, PCAtom *atom_selection2,
							 int n_selected_atoms_1, int n_selected_atoms_2) const {

   std::pair<std::string, std::string> aligned_sequences =
      get_horizontal_ssm_sequence_alignment(SSMAlign, asc_ref, asc_mov,
					    atom_selection1, atom_selection2,
					    n_selected_atoms_1, n_selected_atoms_2);

   print_horizontal_ssm_sequence_alignment(aligned_sequences);
}
#endif // HAVE_SSMLIB

#ifdef HAVE_SSMLIB
// 
// To make a GUI dialog with the alignment, we need this function to
// generate a string, rather than print to the screen.  Easily converted.
//
void
graphics_info_t::print_horizontal_ssm_sequence_alignment(std::pair<std::string, std::string> aligned_sequences) const {

   bool debug = 0;
   
   if (debug) { 
      std::cout << "DEBUG:: moving :" << aligned_sequences.first  << ":" << std::endl;
      std::cout << "DEBUG:: target :" << aligned_sequences.second << ":" << std::endl;
   }

   int chars_per_line = 70;
   int lf = aligned_sequences.first.length();
   int ls = aligned_sequences.second.length();
   int l = lf;
   if (ls > l)
      l = ls;
   int n_lines = 1 + l/chars_per_line;
   // std::cout << "DEUBG:: n_lines: " << n_lines << " " << lf << " " << ls << std::endl;

   for (unsigned int i=0; i<n_lines; i++) {
      int f_start = i*chars_per_line;
      int f_end = chars_per_line;
      if (f_end > lf)
	 f_end = lf - f_start;
      // std::cout << "DEBUG:: comparing  first " << f_start << " with " << lf << std::endl;
      if (f_start < lf)
	 std::cout << " Moving: " << aligned_sequences.first.substr(f_start, f_end) << std::endl;
      
      int s_start = i*chars_per_line;
      int s_end = chars_per_line;
      if (s_end > ls)
	 s_end = ls - s_start;
      // std::cout << "DEBUG:: comparing second " << s_start << " with " << ls << std::endl;
      if (s_start < ls)
	 std::cout << " Target: " << aligned_sequences.second.substr(s_start, s_end) << std::endl;
      std::cout << std::endl; // for neatness
   } 
}
#endif // HAVE_SSMLIB


#ifdef HAVE_SSMLIB
std::pair<std::string, std::string>
graphics_info_t::get_horizontal_ssm_sequence_alignment(CSSMAlign *SSMAlign,
						       atom_selection_container_t asc_ref,
						       atom_selection_container_t asc_mov,
						       PCAtom *atom_selection1, PCAtom *atom_selection2,
						       int n_selected_atoms_1, int n_selected_atoms_2) const {

   std::string s;
   std::string t;
   int previous_t_index = -1;
   int previous_s_index = -1;
   bool debug = 0;

   // Talk to Eugene about this: turn on debugging and run Alice
   // Dawson test.
   //
   if (debug) {
      std::cout << "DEBUG:: t_indexes: ";
      for (unsigned int i1=0; i1<SSMAlign->nsel1; i1++) {
	 std::cout << "[" << i1 << " " << SSMAlign->Ca1[i1] << "] " ;
      }
      std::cout << std::endl;
      std::cout << "DEBUG:: s_indexes: ";
      for (unsigned int i2=0; i2<SSMAlign->nsel2; i2++) {
	 std::cout << "[" << i2 << " " << SSMAlign->Ca2[i2] << "] " ;
      }
      std::cout << std::endl;
   }

   for (unsigned int i1=0; i1<SSMAlign->nsel1; i1++) {
      int t_index = SSMAlign->Ca1[i1];
      if (0)
	 std::cout << "debug i1: " << i1 << " t_index is " << t_index << " and range is "
		   << SSMAlign->nsel2 << std::endl;
      if (t_index == -1) {
	 s += coot::util::three_letter_to_one_letter(atom_selection1[i1]->GetResName());
	 t += "-";
      } else {
	 // was t_index sensible?
	 if (t_index < SSMAlign->nsel2) {
	    int s_index = SSMAlign->Ca2[t_index];
	    if (t_index != (previous_t_index + 1)) {
	       // this path rarely happens?
	       int atom_sel2_index = previous_t_index + 1;
	       if (atom_sel2_index < SSMAlign->nsel2) { 
		  CAtom *at = atom_selection2[t_index];
		  CAtom *ti_at = atom_selection2[atom_sel2_index];
		  if (debug) { 
		     std::cout << "t_index: " << t_index << std::endl;
		     std::cout << "atom_sel2_index: " << atom_sel2_index << std::endl;
		     std::cout << "at:    " << at << std::endl;
		     std::cout << "ti_at: " << ti_at << std::endl;
		  }
		  t += coot::util::three_letter_to_one_letter(ti_at->GetResName());
		  s += "."; // there was (at least) one extra residue in t sequence
	       } else {
		  // 20091221 confusion over what to do here
		  // atom_sel2_index is out of bounds (adding the test
		  // for it fixes a crash). Let's add nothing.
		  t += "";
		  s += ""; // there was (at least) one extra residue in t sequence
	       } 
	    }
	    if (s_index == i1) { 
	       s += coot::util::three_letter_to_one_letter(atom_selection1[i1]->GetResName());
	       t += coot::util::three_letter_to_one_letter(atom_selection2[t_index]->GetResName());
	    }
	 } else {
	    t += "^";
	    s += "^";
	 }

	 // for next round
	 previous_t_index = t_index;
      }
   }
   std::cout << std::endl;
   return std::pair<std::string, std::string> (s,t);
}
#endif // HAVE_SSMLIB

#ifdef HAVE_SSMLIB
void
graphics_info_t::print_ssm_sequence_alignment(CSSMAlign *SSMAlign,
					      atom_selection_container_t asc_ref,
					      atom_selection_container_t asc_mov,
					      PCAtom *atom_selection1, PCAtom *atom_selection2,
					      int n_selected_atoms_1, int n_selected_atoms_2,
					      short int move_copy_of_imol2_flag) {

   std::cout << "Another Go...\n\n";

   CChain *moving_chain_p = 0;
   CChain *reference_chain_p = 0;
   std::string mov_chain_id = std::string(atom_selection1[0]->GetChainID());
   std::string ref_chain_id = std::string(atom_selection2[0]->GetChainID());
   std::string slc_1, slc_2; // single letter code.

   int nchains_ref = asc_ref.mol->GetNumberOfChains(1);
   for (int ich=0; ich<nchains_ref; ich++) {
      CChain *chain_p = asc_ref.mol->GetChain(1, ich);
      std::string mol_chain_id(chain_p->GetChainID());
      if (mol_chain_id == std::string(ref_chain_id)) {
	 reference_chain_p = chain_p;
	 break;
      }
   }
   int nchains_mov = asc_mov.mol->GetNumberOfChains(1);
   for (int ich=0; ich<nchains_mov; ich++) {
      CChain *chain_p = asc_mov.mol->GetChain(1, ich);
      std::string mol_chain_id(chain_p->GetChainID());
      if (mol_chain_id == std::string(mov_chain_id)) {
	 moving_chain_p = chain_p;
	 break;
      }
   }

   if (moving_chain_p && reference_chain_p) {

      // print_alignment_table (not sequence)
      // 
      if (n_selected_atoms_1 > 0) {
 	 clipper::RTop_orth ssm_matrix = coot::util::matrix_convert(SSMAlign->TMatrix);
//       this is output later, no need to reproduce it here.
// 	 std::cout << "ssm_matrix: \n" << ssm_matrix.format() << std::endl;
	 std::cout << "      Moving  Reference   Distance" << std::endl;
	 for (int ires=0; ires<n_selected_atoms_1; ires++) {
	    if (ires < SSMAlign->nalgn) { 
	       CAtom *mov_at = atom_selection1[ires];
	    
	       int mov_index = SSMAlign->Ca1[ires];
	       std::cout << "      " << mov_at->GetChainID() << " " << mov_at->GetSeqNum();
	       if ((mov_index > -1) && (mov_index < n_selected_atoms_1)) { 
		  CAtom *ref_at = atom_selection2[mov_index];
		  if (ref_at) {
		     clipper::Coord_orth pos1(mov_at->x, mov_at->y, mov_at->z);
		     clipper::Coord_orth pos2(ref_at->x, ref_at->y, ref_at->z);
		     clipper::Coord_orth pos3 = pos1.transform(ssm_matrix);
		     double d = clipper::Coord_orth::length(pos3, pos2);
		     std::cout << " <---> " << ref_at->GetChainID() << " "
			       << ref_at->GetSeqNum() << "  : " << d << "  "
			       << " A\n";
		  }
	       } else {
		  std::cout << "\n";
	       }
	    }
	 }
      }
   } else {
      std::cout << "ERROR:: Failed to get moving or reference_chain pointer\n";
   }
}
#endif
