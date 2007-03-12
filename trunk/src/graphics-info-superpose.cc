/* src/graphics-info-superpose.cc
 * 
 * Copyright 2004, 2005 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#if defined _MSC_VER
#include <windows.h>
#endif

#include <iostream>
#include "graphics-info.h"
#ifdef HAVE_SSMLIB
#include "ssm_align.h"
#endif

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
   
   if (imol >=0 && imol < n_molecules) { 
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

   char *data = NULL;
   data = (char *)pos;
   // this can fail when more than one sequence mutate is used at the same time:
   if (data) 
      graphics_info_t::superpose_imol1_chain = data;
}

// This does a redraw
// 
// static
void
graphics_info_t::superpose_moving_chain_option_menu_item_activate (GtkWidget *item,
								   GtkPositionType pos) { 

   char *data = NULL;
   data = (char *)pos;
//    std::cout << "INFO:: superpose_moving_chain_option_menu_item_activate "
// 	     << " got data: " << data << std::endl;
   // this can fail when more than one sequence mutate is used at the same time:
   if (data) 
      graphics_info_t::superpose_imol2_chain = data;
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
void
graphics_info_t::superpose_with_atom_selection(atom_selection_container_t asc_ref,
					       atom_selection_container_t asc_mov,
					       int imol2,
					       std::string moving_mol_name,
					       std::string reference_mol_name,
					       short int move_copy_of_imol2_flag) {

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
	       std::string name = "Copy of ";
	       name += moving_mol_name;
	       graphics_info_t::molecules[graphics_info_t::n_molecules].install_model(make_asc(mol2), name, 1);
	       imol2 = graphics_info_t::n_molecules;
	       graphics_info_t::n_molecules++;
	    }
	    
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
	    std::string slc_1, slc_2; // single letter code.
	    int rn1, rn2;

	    // another go at the sequence alignment used in the superposition.
	    // 
	    // correct: number of residues in reference structure: " << SSMAlign->nres2

	    std::cout << "Another Go...\n\n";

	    CChain *moving_chain_p = 0;
	    CChain *reference_chain_p = 0;
	    std::string mov_chain_id = std::string(atom_selection1[0]->GetChainID());
	    std::string ref_chain_id = std::string(atom_selection2[0]->GetChainID());

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

// 	       std::cout << "Reference Chain: "<< reference_chain_p->GetChainID() << " has "
// 			 << reference_chain_p->GetNumberOfResidues() << " residues\n";
// 	       std::cout << "Moving    Chain: "<< moving_chain_p->GetChainID() << " has "
// 			 << moving_chain_p->GetNumberOfResidues() << " residues\n";

	       std::vector<std::string> ref_seq(n_selected_atoms_2);
	       std::vector<std::string> mov_seq(n_selected_atoms_1);

	       for (int i=0; i<n_selected_atoms_2; i++)
		  ref_seq[i] = coot::util::downcase(coot::util::three_letter_to_one_letter(atom_selection2[i]->GetResName()));
	       for (int i=0; i<n_selected_atoms_1; i++)
		  mov_seq[i] = coot::util::downcase(coot::util::three_letter_to_one_letter(atom_selection1[i]->GetResName()));

// 	       std::cout << "Ca1 indexes: \n";
// 	       for (int i=0; i<n_selected_atoms_1; i++)
// 		  std::cout << SSMAlign->Ca1[i] << " ";
// 	       std::cout << "\n";
// 	       std::cout << "Ca2 indexes: \n";
// 	       for (int i=0; i<n_selected_atoms_2; i++)
// 		  std::cout << SSMAlign->Ca2[i] << " ";
// 	       std::cout << "\n";
		     
// 	       std::cout << "Simple mov sequence:\n";
// 	       for (int i=0; i<n_selected_atoms_1; i++)
// 		  std::cout << mov_seq[i];
// 	       std::cout << "\n";
// 	       std::cout << "Simple ref sequence:\n";
// 	       for (int i=0; i<n_selected_atoms_2; i++)
// 		  std::cout << ref_seq[i];
// 	       std::cout << "\n";

	       std::string ral = "";
	       std::string mal = "";

	       std::string push_ral = "";
	       std::string push_mal = "";

	       std::string pushed_ral_char = " ";
	       std::string pushed_mal_char = " ";

	       std::string sum_push_ral = "";
	       std::string sum_push_mal = "";

	       // reference 2:
	       // moving    1:

	       int n_max_residues = SSMAlign->nres1;
	       if (SSMAlign->nres2 > n_max_residues)
		  n_max_residues = SSMAlign->nres2;

// 	       slc_1 = "-";
// 	       slc_2 = "-";

	       int prev_diff_1 = 0;
	       int prev_diff_2 = 0;
	       int diff_1 = 1; 
	       int diff_2 = 1; 
	       for (int ires=0; ires<n_max_residues; ires++) {

		  slc_1 = " ";
		  slc_2 = " ";
		  if (ires<SSMAlign->nres1) { 
		     int mov_index = SSMAlign->Ca1[ires];
		     if ((SSMAlign->Ca1[ires]>=0) && (SSMAlign->Ca1[ires-1]>=0)) { 
			diff_1 =  SSMAlign->Ca1[ires] - SSMAlign->Ca1[ires-1];
			if (diff_1 > prev_diff_1) { 
			   std::cout << "dashing mal for " << ires << " " << diff_1 << " "
				     << prev_diff_1 << " " << std::endl;
			   for (int igap=0; igap<(diff_1-prev_diff_1); igap++)
			      mal += "-";
			}
		     }
		     if (mov_index > -1) { 
			slc_1 = coot::util::three_letter_to_one_letter(atom_selection1[ires]->GetResName());
			pushed_ral_char = "";
		     } else { 
			slc_1 = coot::util::downcase(coot::util::three_letter_to_one_letter(atom_selection1[ires]->GetResName()));
		     }

		     // push along the reference chain (maybe)
		     if (mov_index == -1)
			push_ral = pushed_ral_char;
		     else
			push_ral = "";
			
		  }
		  if (ires<SSMAlign->nres2) { 
		     int ref_index = SSMAlign->Ca2[ires];
		     if ((SSMAlign->Ca2[ires]>=0) && (SSMAlign->Ca2[ires-1]>=0)) { 
			diff_2 =  SSMAlign->Ca2[ires] - SSMAlign->Ca2[ires-1];
			if (diff_2 > prev_diff_2) { 
			   std::cout << "dashing ral for " << ires << " " << diff_1 << " "
				     << prev_diff_1 << " " << std::endl;
			   for (int igap=0; igap<(diff_1-prev_diff_1); igap++)
			      ral += "-";
			}
		     }
		     if (ref_index > -1) {
			// for neatness a function to do this.
			slc_2 = coot::util::three_letter_to_one_letter(atom_selection2[ires]->GetResName());
			pushed_mal_char = "";
		     } else {
			slc_2 = coot::util::downcase(coot::util::three_letter_to_one_letter(atom_selection2[ires]->GetResName()));
		     }

		     // push along the moving chain (maybe)
		     if (ref_index == -1) { 
			push_mal = pushed_mal_char;
		     } else {
			push_mal = "";
		     }
		  }

// 		  std::cout << "DEbug for ires " << ires << " adding :" << slc_1 << ": to mal"
// 			    << " and " << slc_2 << " to ral" << std::endl;

// 		  std::cout << "DEbug for ires " << ires << " adding :" << push_mal << ": to sum_push_mal"
// 			    << " and " << push_ral << " to sum_push_ral" << std::endl;

		  if (ral.length() == 0) 
		     sum_push_mal += push_mal;
		  sum_push_ral += push_ral;
		  mal += slc_1;
		  ral += slc_2;
		  prev_diff_1 = diff_1;
		  prev_diff_2 = diff_2;
	       }
	       std::cout << "INFO:: alignment (leading spacing might be wrong): " << std::endl;
	       std::cout << "mov: " << sum_push_mal << mal << std::endl;
	       std::cout << "ref: " << sum_push_ral << ral << std::endl;

	       // print_alignment_table : just too painful to put an CSSMAlign into graphics-info.h
	       // 
	       if (SSMAlign->nres1 > 0) {
		  std::cout << "      Moving  Reference   Distance" << std::endl;
		  for (int ires=0; ires<SSMAlign->nres1; ires++) {
		     CAtom *mov_at = atom_selection1[ires];
		     std::cout << "      " << mov_at->GetChainID() << " " << mov_at->GetSeqNum();
		     
		     int mov_index = SSMAlign->Ca1[ires];
		     if (mov_index > -1) { 
			CAtom *ref_at = atom_selection2[mov_index];
			clipper::Coord_orth pos1(mov_at->x, mov_at->y, mov_at->z);
			clipper::Coord_orth pos2(ref_at->x, ref_at->y, ref_at->z);
			double d = clipper::Coord_orth::length(pos1, pos2);
			if (move_copy_of_imol2_flag) { 
			   clipper::Coord_orth pos3 =
			      pos1.transform(coot::util::matrix_convert(SSMAlign->TMatrix));
			   d = clipper::Coord_orth::length(pos3, pos2);
			}
			std::cout << " <---> " << ref_at->GetChainID() << " "
				  << ref_at->GetSeqNum() << "  : " << d << "  "
// 				  << pos1.format() << "   "
// 				  << pos2.format() << "   "
// 				  << pos3.format() << "   "
				  << " A\n";
		     } else {
			std::cout << "\n";
		     } 
		  }
	       }

	    } else {
	       std::cout << "ERROR:: Failed to get moving or reference_chain pointer\n";
	    }
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
}

