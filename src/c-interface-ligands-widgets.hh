/*
 * src/c-interface-ligands-widgets.hh
 *
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */


// widgets for ligand-fitting, specifically.

#ifndef C_INTERFACE_LIGANDS_WIDGETS_HH
#define C_INTERFACE_LIGANDS_WIDGETS_HH

class ligand_wiggly_ligand_data_t {
   void init() {
      finish = false;
      immediate_execute_ligand_search = true; // unless we have wiggly ligands
      progress_bar = 0;
      progress_bar_label = 0;
      progress_bar_window = 0;
      wlig = 0;
   }
public:
   int imol_ligand;
   coot::wligand *wlig;
   GtkWidget *progress_bar;
   GtkWidget *progress_bar_window;
   GtkWidget *progress_bar_label;
   bool finish;
   bool immediate_execute_ligand_search;
   ligand_wiggly_ligand_data_t() {
      init();
   }
   ligand_wiggly_ligand_data_t(coot::wligand *wlig_in) {
      init();
      wlig = wlig_in;
   }
};

ligand_wiggly_ligand_data_t setup_ligands_progress_bar(); // return the progress bar
void setup_ligands_progress_bar_idle(coot::wligand *wlig,
				     int imol_ligand,
				     ligand_wiggly_ligand_data_t ld);

gboolean install_simple_wiggly_ligand_idle_fn(gpointer data);

#endif // C_INTERFACE_LIGANDS_WIDGETS_HH
