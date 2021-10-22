/* src/positioned-widgets.h
 * 
 * Copyright 2005 by The University of York
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, 
 * MA 02110-1301, USA.
 */



// don't overlap with those defined in c-inteface.h
// These have statics int describing the position:
#define COOT_ROTAMER_SELECTION_DIALOG   1001
#define COOT_EDIT_CHI_DIALOG            1002

#define COOT_UNDEFINED_WINDOW               0
#define COOT_DELETE_WINDOW                  1
#define COOT_MUTATE_RESIDUE_RANGE_WINDOW    2
#define COOT_ACCEPT_REJECT_WINDOW           3
#define COOT_DISTANCES_ANGLES_WINDOW        4
#define COOT_ROTATE_TRANSLATE_DIALOG        5
#define COOT_DISPLAY_CONTROL_WINDOW         6
#define COOT_DISPLAY_CONTROL_MAPS_VBOX      7 
#define COOT_DISPLAY_CONTROL_MOLECULES_VBOX 8
#define COOT_DISPLAY_CONTROL_PANE           9
#define COOT_ACCESSION_CODE_WINDOW_OCA     10
#define COOT_ACCESSION_CODE_WINDOW_OCA_WITH_SF 11
#define COOT_ACCESSION_CODE_WINDOW_EDS     12
#define COOT_ACCESSION_CODE_WINDOW_PDB_REDO 112
#define COOT_MODEL_REFINE_DIALOG           13
#define COOT_GO_TO_ATOM_WINDOW             14
#define COOT_RAMACHANDRAN_PLOT_WINDOW      15
#define COOT_FILESELECTION_DIALOG          16
#define COOT_REFINEMENT_CONTROL_DIALOG     17

#define COOT_CHECKED_WATERS_BADDIES_DIALOG 18
#define COOT_CHECK_WATERS_DIALOG           29
#define COOT_MTZ_COLUMN_SELECTOR_DIALOG    19
#define COOT_DIFF_MAPS_PEAK_DIALOG         30
#define COOT_EDIT_RESTRAINTS_DIALOG        31

#define COOT_SCREENDUMP_SIMPLE   20
#define COOT_SCREENDUMP_RASTER3D 21
#define COOT_SCREENDUMP_POVRAY   22

#define COOT_FIXED_ATOM_DIALOG   23
