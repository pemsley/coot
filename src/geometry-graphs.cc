/* src/graphics-info.cc
 * 
 * Copyright 2004, 2005 by The University of York
 * Author Paul Emsley, Bernhard Lohkamp
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

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "coot-utils/coot-coord-utils.hh" // compilation error on MacOSX if this
                               // doesn't come before the next 3 lines
                               // (for gtk/gtk-canvas) (I suspect
                               // canvas is the problem).

#include "utils/coot-utils.hh"

#include <gtk/gtk.h>
// #include <gdk_imlib.h>
#include <goocanvas.h>

#include "geometry-graphs.hh"
#include "interface.h"
#include "support.h"  // for lookup

#include "c-interface.h" // for the coot graphics callbacks in button_press()
#include "cc-interface.hh"
#include "graphics-info.h"


gint
coot::on_geometry_graph_block_clicked(GooCanvasItem *item, GooCanvasItem *target, GdkEvent *event, gpointer user_data) {

   coot::geometry_graph_block_info *block_info_p = static_cast<coot::geometry_graph_block_info*>(user_data);
   if (block_info_p) {
      const coot::geometry_graph_block_info &block_info(*block_info_p);

      if (false)
         std::cout << "INFO:: trying to centre on: " << block_info.atom_spec.chain_id << " "
                   << block_info.atom_spec.res_no << " "
                   << block_info.atom_spec.atom_name << std::endl;

      set_go_to_atom_molecule(block_info.imol);
      set_go_to_atom_chain_residue_atom_name(block_info.atom_spec.chain_id.c_str(),
                                             block_info.atom_spec.res_no,
                                             block_info.atom_spec.atom_name.c_str());
   }
   return TRUE; // it's been handled.
}

// this could be a static?
//
void
coot::geometry_graphs::button_press(GooCanvasItem *item, GdkEvent *event,
                                    const coot::geometry_graph_block_info &binfo) {

   // "button_press" event was for the old gnomecanvas. We probably don't need this now.

   if (false)
      std::cout << "INFO:: trying to centre on: " << binfo.atom_spec.chain_id << " "
                << binfo.atom_spec.res_no << " "
                << binfo.atom_spec.atom_name << std::endl;

   set_go_to_atom_molecule(binfo.imol);
   set_go_to_atom_chain_residue_atom_name(binfo.atom_spec.chain_id.c_str(),
                                          binfo.atom_spec.res_no,
                                          binfo.atom_spec.atom_name.c_str());

}

// not with a tooltip, it can't.
//
void
coot::geometry_graphs::mouse_over(GooCanvasItem *item, GdkEvent *event,
                                  const coot::geometry_graph_block_info &block_info) {

   // do a tooltip
   //std::cout << block_info.resno << std::endl;

   // tooltip_like_box(block_info, event);

   std::cout << "Don't do tooltips like this" << std::endl;


}

// "Molecule 2: " (for example) gets prepended to the graph_label_in
// before becoming at GtkLabel.
//
coot::geometry_graphs::geometry_graphs(coot::geometry_graph_type graph_type_in,
                                       int imol_in,
                                       std::string graph_label_in,
                                       int nchains_in, int max_chain_length_in) {

   graph_type = graph_type_in;
   imol = imol_in;
   n_chains = nchains_in;
   max_chain_length = max_chain_length_in;
   setup_internal();
   setup_canvas(nchains_in, max_chain_length_in);
   // std::cout << "DEBUG:: in geometry_graphs constuctor resizing geometry_graphs block to "
   // << nchains_in << std::endl;
   blocks.resize(nchains_in);
   offsets.resize(nchains_in);
   chain_index.resize(n_chains);
   graph_label = "Molecule ";
   graph_label += coot::util::int_to_string(imol);
   graph_label += ": ";
   graph_label += graph_label_in;
   if (graph_type == GEOMETRY_GRAPH_DENSITY_FIT) {
      graph_label += " fit vs. Map ";
      graph_label += coot::util::int_to_string(imol_refinement_map());
   }

   // set the window title:
   std::string title("Graph");
   switch(graph_type) {
        case coot::GEOMETRY_GRAPH_GEOMETRY:
                title = _("Geometry Graphs");
                break;
////B
        case coot::GEOMETRY_GRAPH_CALC_B_FACTOR:
                title = _(" < B > Factor Graphs ");
                break;
////E
        case coot::GEOMETRY_GRAPH_B_FACTOR:
                title = _(" B Factor Variance Graphs ");
                break;
        case coot::GEOMETRY_GRAPH_DENSITY_FIT:
                title = _("Density Fit Graphs");
                break;
        case coot::GEOMETRY_GRAPH_OMEGA_DISTORTION:
                title = _("Peptide Omega Distortion Graphs");
                break;
        case coot::GEOMETRY_GRAPH_ROTAMER:
                title = _("Unusual Rotamer Graphs");
                break;
        case coot::GEOMETRY_GRAPH_NCS_DIFFS:
                title = _("NCS Differences");
                break;
        default:
                break;
   }

   // adjust distortion_max depending on graph type
   if (graph_type == coot::GEOMETRY_GRAPH_OMEGA_DISTORTION)
      distortion_max = 90.0; // degrees

   gtk_window_set_title(GTK_WINDOW(dialog()), title.c_str());

   // And now the graph label for EJD:
   GtkWidget *label = lookup_widget(dialog(), "geometry_graphs_label");
   gtk_label_set_text(GTK_LABEL(label), graph_label.c_str());
}

void
coot::geometry_graphs::render_to_canvas(const coot::geometry_distortion_info_container_t &dc,
                                        int chain_number) {

   // First we get min and max from the distortion info container.
   // This tells us how long to draw the x axis line.  We want to do
   // this the first time only.
   //
   // In the case of an update, the x axis line is already drawn.  In
   // that case we just need to draw the block in the right place.
   // The right place is determined by the residue number and the
   // minimum residue number, which is given(determined) in the render
   // to canvas call (and we need to look up the chain id of the
   // updated residues in the chain-id->chain number conversion
   // vector).
   //
   // Let's call the place for the residue after account has been
   // taken of the minimum residue number in the chain: let's call it
   // offset_residue_place.

//    std::cout << "INFO:: there are " << dc.geometry_distortion.size()
//              << " distortions in container " << chain_number << std::endl;

   if (chain_number < static_cast<int>(chain_index.size()))
      chain_index[chain_number] = dc.chain_id;

   // dc can contain "unset" values for max_resno and min_resno, so check them here
   //
   int max_resno = dc.max_resno;
   int min_resno = dc.min_resno;
   if (max_resno >= 0) {
      int nres = max_resno - min_resno + 1;
      offsets[chain_number] = min_resno - 1;
      if (true) // debug
         std::cout << "::::::::::: in render_to_canvas() offsets[" << chain_number << "] is set to "
                   << offsets[chain_number] << std::endl;

      draw_chain_axis(nres, chain_number);
      draw_chain_axis_tick_and_tick_labels(min_resno, max_resno, chain_number);

      if (true) {
         std::cout << "DEBUG:: blocks.size(): " << blocks.size() << std::endl;
         std::cout << "DEBUG:: resizing blocks[" << chain_number << "] to " << nres << std::endl;
      }

      int n_blocks = blocks.size();
      if (chain_number < n_blocks) {
         if (nres > 0) {
            blocks[chain_number].resize(nres + 1); // needs to index max_resno
            render_geometry_distortion_blocks_internal_linear(dc, min_resno, max_resno);
            label_chain(dc.chain_id, chain_number); // labels last so they are on top.
         }
      }
   }


}

void
coot::geometry_graphs::render_geometry_distortion_blocks_internal(const coot::geometry_distortion_info_container_t &dc) {

   std::map<coot::residue_spec_t, std::pair<double, std::string> > residue_distortions; // for block height
   std::map<coot::residue_spec_t, std::pair<coot::geometry_distortion_info_t, std::string> > worst_distortions;
   int chain_number = chain_id_to_chain_index(dc.chain_id);

   int idx_1, idx_2, idx_3, idx_4;
   mmdb::realtype occ_1, occ_2, occ_3;
   for (unsigned int i=0; i<dc.geometry_distortion.size(); i++) {

      if (false)
         std::cout << "now examining restraint number " << i << " type "
                   << dc.geometry_distortion[i].restraint.restraint_type << std::endl;

      if (dc.geometry_distortion[i].restraint.restraint_type == coot::BOND_RESTRAINT) {
         idx_1 = dc.geometry_distortion[i].restraint.atom_index_1;
         idx_2 = dc.geometry_distortion[i].restraint.atom_index_2;
         coot::residue_spec_t rs_1(dc.atom[idx_1]->GetResidue());
         coot::residue_spec_t rs_2(dc.atom[idx_2]->GetResidue());
         occ_1 = sane_occupancy(dc.atom[idx_1]->occupancy);
         occ_2 = sane_occupancy(dc.atom[idx_2]->occupancy);
         coot::geometry_distortion_info_t extra_distortion = dc.geometry_distortion[i];
         // we try here because if worst_distortions[rs_] is undefined
         // (as it will be the first time we use it then we are doing
         // a comparison against an undefined distorition
         //
         try {
            if (extra_distortion < worst_distortions[rs_1].first) {
            } else {
               worst_distortions[rs_1].first = extra_distortion;
               worst_distortions[rs_1].second = make_distortion_string(dc.geometry_distortion[i], dc);
            }
         }
         catch (const std::runtime_error &rte) {
            if (extra_distortion.initialised_p()) {
               worst_distortions[rs_1].first = extra_distortion;
               worst_distortions[rs_1].second = make_distortion_string(dc.geometry_distortion[i], dc);
            }
         }
         residue_distortions[rs_1].first += 0.5 * extra_distortion.distortion_score * occ_1;
         residue_distortions[rs_2].first += 0.5 * extra_distortion.distortion_score * occ_2;
      }

      if (dc.geometry_distortion[i].restraint.restraint_type == coot::ANGLE_RESTRAINT) {
         idx_1 = dc.geometry_distortion[i].restraint.atom_index_1;
         idx_2 = dc.geometry_distortion[i].restraint.atom_index_2;
         idx_3 = dc.geometry_distortion[i].restraint.atom_index_2;
         coot::residue_spec_t rs_1(dc.atom[idx_1]->GetResidue());
         coot::residue_spec_t rs_2(dc.atom[idx_2]->GetResidue());
         coot::residue_spec_t rs_3(dc.atom[idx_3]->GetResidue());
         occ_1 = sane_occupancy(dc.atom[idx_1]->occupancy);
         occ_2 = sane_occupancy(dc.atom[idx_2]->occupancy);
         occ_3 = sane_occupancy(dc.atom[idx_3]->occupancy);
         coot::geometry_distortion_info_t extra_distortion = dc.geometry_distortion[i];

         // we try here because if worst_distortions[rs_] is undefined
         // (as it will be the first time we use it then we are doing
         // a comparison against an undefined distorition
         //
         try {
            if (extra_distortion < worst_distortions[rs_1].first) {
            } else {
               worst_distortions[rs_1].first = extra_distortion;
               worst_distortions[rs_1].second = make_distortion_string(dc.geometry_distortion[i], dc);
            }
         }
         catch (const std::runtime_error &rte) {
            if (extra_distortion.initialised_p()) {
               worst_distortions[rs_1].first = extra_distortion;
               worst_distortions[rs_1].second = make_distortion_string(dc.geometry_distortion[i], dc);
            }
         }
         residue_distortions[rs_1].first += 0.333 * extra_distortion.distortion_score * occ_1;
         residue_distortions[rs_2].first += 0.333 * extra_distortion.distortion_score * occ_2;
         residue_distortions[rs_3].first += 0.333 * extra_distortion.distortion_score * occ_2;
      }

      if (dc.geometry_distortion[i].restraint.restraint_type == coot::PLANE_RESTRAINT) {
         int n_atoms = dc.geometry_distortion[i].restraint.plane_atom_index.size();
         if (n_atoms > 0) { 
            double factor = 1/double(n_atoms);
            for (unsigned int iat=0; iat<dc.geometry_distortion[i].restraint.plane_atom_index.size(); iat++) {
               idx_1 = dc.geometry_distortion[i].restraint.plane_atom_index[iat].first;
               occ_1 = sane_occupancy(dc.atom[idx_1]->occupancy);
               coot::geometry_distortion_info_t extra_distortion = dc.geometry_distortion[i];
               coot::residue_spec_t rs_1(dc.atom[idx_1]->GetResidue());
               residue_distortions[rs_1].first += factor * extra_distortion.distortion_score * occ_1;

               // we try here because if worst_distortions[rs_] is undefined
               // (as it will be the first time we use it then we are doing
               // a comparison against an undefined distorition
               // 
               try { 
                  if (extra_distortion < worst_distortions[rs_1].first) {
                  } else {
                     worst_distortions[rs_1].first = extra_distortion;
                     worst_distortions[rs_1].second = make_distortion_string(dc.geometry_distortion[i], dc);
                  }
               }
               catch (const std::runtime_error &rte) {
                  if (extra_distortion.initialised_p()) { 
                     worst_distortions[rs_1].first = extra_distortion;
                     worst_distortions[rs_1].second = make_distortion_string(dc.geometry_distortion[i], dc);
                  }
               }

            }
         }
      }
   }
   std::map<coot::residue_spec_t, std::pair<double, std::string> >::iterator it;
   for (it=residue_distortions.begin(); it!=residue_distortions.end(); it++) {
      it->second.second = worst_distortions[it->first].second;
   }
   plot_blocks(residue_distortions, chain_number);
}

std::string
coot::geometry_graphs::make_distortion_string(const coot::geometry_distortion_info_t &geometry_distortion,
                                              const coot::geometry_distortion_info_container_t &dc) const {

   std::string s;
   if (geometry_distortion.restraint.restraint_type == coot::BOND_RESTRAINT) {
      int idx_1 = geometry_distortion.restraint.atom_index_1;
      int idx_2 = geometry_distortion.restraint.atom_index_2;
      coot::residue_spec_t rs_1(dc.atom[idx_1]->GetResidue());
      s += dc.atom[idx_1]->GetChainID();
      s += " ";
      s += coot::util::int_to_string(dc.atom[idx_1]->GetSeqNum());
      s += " ";
      s += dc.atom[idx_1]->GetResName();
      s += " ";
      s  = "Bond: ";
      s += dc.atom[idx_1]->name;
      s += "-";
      s += dc.atom[idx_2]->name;
      s += " z score: ";
      s += coot::util::float_to_string(sqrt(geometry_distortion.distortion_score));
   }
   
   if (geometry_distortion.restraint.restraint_type == coot::ANGLE_RESTRAINT) {
      int idx_1 = geometry_distortion.restraint.atom_index_1;
      int idx_2 = geometry_distortion.restraint.atom_index_2;
      int idx_3 = geometry_distortion.restraint.atom_index_3;
      coot::residue_spec_t rs_1(dc.atom[idx_1]->GetResidue());
      s += dc.atom[idx_1]->GetChainID();
      s += " ";
      s += coot::util::int_to_string(dc.atom[idx_1]->GetSeqNum());
      s += " ";
      s += dc.atom[idx_1]->GetResName();
      s += " ";
      s  = "Angle: ";
      s += dc.atom[idx_1]->name;
      s += "-";
      s += dc.atom[idx_2]->name;
      s += "-";
      s += dc.atom[idx_3]->name;
      s += " z score: ";
      s += coot::util::float_to_string(sqrt(geometry_distortion.distortion_score));
   }
   if (geometry_distortion.restraint.restraint_type == coot::PLANE_RESTRAINT) {
      if (geometry_distortion.restraint.plane_atom_index.size() > 0) {
         int idx = geometry_distortion.restraint.plane_atom_index[0].first;
         s += dc.atom[idx]->GetChainID();
         s += " ";
         s += coot::util::int_to_string(dc.atom[idx]->GetSeqNum());
         s += " ";
         s += dc.atom[idx]->GetResName();
         s += " ";
      }
      s = "Plane: ";
      for (unsigned int iat=0; iat<geometry_distortion.restraint.plane_atom_index.size(); iat++) {
         int idx = geometry_distortion.restraint.plane_atom_index[iat].first;
         s += dc.atom[idx]->name;
         s += " ";
      }
      s += " z score: ";
      s += coot::util::float_to_string(sqrt(geometry_distortion.distortion_score));
   }
   // std::cout << "make_distortion_string() " << s << std::endl;
   return s;
}


void
coot::geometry_graphs::plot_blocks(const std::map<coot::residue_spec_t, std::pair<double, std::string> > &residue_distortions,
                                   int chain_number) {

   float magic_scale = 0.6;
   std::map<coot::residue_spec_t, std::pair<double, std::string> >::const_iterator it;
   std::string distortion_string;
   
   for (it=residue_distortions.begin(); it != residue_distortions.end(); it++) {
      // std::cout << "plot residue " << it->first << " with height " << it->second << std::endl;

      // make an atom spec (graphics_info is not allowed here (well, not included anyway)).
      // 
      coot::atom_spec_t try_atom_spec(it->first.chain_id, it->first.res_no, it->first.ins_code,
                                      " CA ", "");

      geometry_graph_block_info block_info(imol,
                                           it->first.res_no,
                                           try_atom_spec,
                                           magic_scale * it->second.first,
                                           it->second.second, canvas);
      plot_block(block_info, offsets[chain_number], chain_number);
   }

}



void
coot::geometry_graphs::render_geometry_distortion_blocks_internal_linear(const coot::geometry_distortion_info_container_t &dc,
                                                                  int min_resno, int max_resno) {

   int this_resno1;
   int this_resno2;
   int this_resno3;
   int idx1, idx2, idx3;
   double occ1, occ2, occ3;
   int NO_DISTORTION_IN_THIS_RESIDUE = -9999;

   std::vector<double> distortion_sum  (max_resno - min_resno + 1, 0);
   std::vector<double> distortion_worst(max_resno - min_resno + 1, 0);
   std::vector<std::string> atom_of_distortion(max_resno - min_resno + 1, " CA ");

   // resi_of_distortion gets set by bond restraints (only).
   std::vector<int>         resi_of_distortion(max_resno - min_resno + 1,
                                               NO_DISTORTION_IN_THIS_RESIDUE); // unassigned.
   std::vector<std::string> distortion_string(distortion_sum.size());

   // std::cout << "DEBUG:: dc has size " << dc.size() << " in render_blocks_internal\n";
   
   for (unsigned int i=0; i<dc.geometry_distortion.size(); i++) {
      std::string info_stub;
      std::string info;
      if (dc.geometry_distortion[i].restraint.restraint_type == coot::BOND_RESTRAINT) {
         idx1 = dc.geometry_distortion[i].restraint.atom_index_1;
         idx2 = dc.geometry_distortion[i].restraint.atom_index_2;
         this_resno1 = dc.atom[idx1]->GetSeqNum();
         this_resno2 = dc.atom[idx2]->GetSeqNum();
         occ1 = sane_occupancy(dc.atom[idx1]->occupancy);
         occ2 = sane_occupancy(dc.atom[idx2]->occupancy);
         info_stub  = coot::util::int_to_string(this_resno1);
         info_stub += dc.atom[idx1]->GetChainID();
         info_stub += " ";
         info_stub += dc.atom[idx1]->GetResName();
//          std::cout << "restraint " <<  i  << " is a bond restraint between " 
//                    << dc.atom[idx1]->GetSeqNum()
//                    << " to residue " 
//                    << dc.atom[idx2]->GetSeqNum()
//                    << " " 
//                    << dc.geometry_distortion[i].distortion_score << "\n";
         
         double extra = dc.geometry_distortion[i].distortion_score;
         // std::cout << "Bond restraint extra " << extra << std::endl;
         if (extra > distortion_worst[this_resno1 - min_resno]) {
            info  = info_stub;
            info += " Bond: ";
            info += dc.atom[idx1]->name;
            info += " ";
            info += dc.atom[idx2]->name;
            info += " z score: ";
            info += coot::util::float_to_string(sqrt(extra));
            //             std::cout << "new worst " << info << std::endl;
            atom_of_distortion[this_resno1 - min_resno] = dc.atom[idx1]->name;
//             std::cout << "DEBUG::  BOND_RESTRAINT setting resi_of_distortion["
//                       << this_resno1 - min_resno << "] to " << dc.atom[idx1]->GetSeqNum()
//                       << std::endl;
            resi_of_distortion[this_resno1 - min_resno] = dc.atom[idx1]->GetSeqNum();
            distortion_string[this_resno1 - min_resno] = info;
            distortion_worst[this_resno1 - min_resno] = extra;
         }

         distortion_sum[this_resno1 - min_resno] += 0.5 * extra * occ1;
         distortion_sum[this_resno2 - min_resno] += 0.5 * extra * occ2;
      }
      if (dc.geometry_distortion[i].restraint.restraint_type == coot::ANGLE_RESTRAINT) {
         idx1 = dc.geometry_distortion[i].restraint.atom_index_1;
         idx2 = dc.geometry_distortion[i].restraint.atom_index_2;
         idx3 = dc.geometry_distortion[i].restraint.atom_index_3;
         this_resno1 = dc.atom[idx1]->GetSeqNum();
         this_resno2 = dc.atom[idx2]->GetSeqNum();
         this_resno3 = dc.atom[idx3]->GetSeqNum();
         occ1 = sane_occupancy(dc.atom[idx1]->occupancy);
         occ2 = sane_occupancy(dc.atom[idx2]->occupancy);
         occ3 = sane_occupancy(dc.atom[idx3]->occupancy);
         info_stub  = coot::util::int_to_string(this_resno1);
         info_stub += dc.atom[idx1]->GetChainID();
         info_stub += " ";
         info_stub += dc.atom[idx1]->GetResName();
         double extra = dc.geometry_distortion[i].distortion_score;
         // std::cout << "Angle restraint extra " << extra << std::endl;
         if (extra > distortion_worst[this_resno1 - min_resno]) {
            info  = info_stub;
            info += " Angle: ";
            info += coot::util::remove_whitespace(dc.atom[idx1]->name);
            info += "-";
            info += coot::util::remove_whitespace(dc.atom[idx2]->name);
            info += "-";
            info += coot::util::remove_whitespace(dc.atom[idx3]->name);
            info += " z score: ";
            info += coot::util::float_to_string(sqrt(extra));
            // std::cout << "new worst " << info << std::endl;
            distortion_string[this_resno1 - min_resno] = info;
            atom_of_distortion[this_resno1 - min_resno] = dc.atom[idx2]->name;
            distortion_worst[this_resno1 - min_resno] = extra;
         }
         distortion_sum[this_resno1 - min_resno] += 0.333 * dc.geometry_distortion[i].distortion_score * occ1;
         distortion_sum[this_resno2 - min_resno] += 0.333 * dc.geometry_distortion[i].distortion_score * occ2;
         distortion_sum[this_resno3 - min_resno] += 0.333 * dc.geometry_distortion[i].distortion_score * occ3;
      }

      if (dc.geometry_distortion[i].restraint.restraint_type == coot::PLANE_RESTRAINT) {
         double factor = 1/double(dc.geometry_distortion[i].restraint.plane_atom_index.size());
         for (unsigned int iat=0; iat<dc.geometry_distortion[i].restraint.plane_atom_index.size(); iat++) {
            idx1 = dc.geometry_distortion[i].restraint.plane_atom_index[iat].first;
            if (idx1 >= 0) { 
               this_resno1 = dc.atom[idx1]->GetSeqNum();
               occ1 = sane_occupancy(dc.atom[idx1]->occupancy);
               info_stub  = coot::util::int_to_string(this_resno1);
               info_stub += dc.atom[idx1]->GetChainID();
               info_stub += " ";
               info_stub += dc.atom[idx1]->GetResName();
               double extra = dc.geometry_distortion[i].distortion_score;
               // std::cout << "Plane restraint extra " << extra << std::endl;
               if (extra > distortion_worst[this_resno1 - min_resno]) {
                  info = info_stub;
                  info += " Plane distortion at: ";
                  info += dc.atom[idx1]->GetChainID();
                  info += coot::util::int_to_string(dc.atom[idx1]->GetSeqNum());
                  info += " z score: ";
                  info += coot::util::float_to_string(sqrt(extra));
                  // std::cout << "new worst " << info << std::endl;
                  distortion_string[this_resno1 - min_resno] = info;
                  atom_of_distortion[this_resno1 - min_resno] = dc.atom[idx1]->name;
                  distortion_worst[this_resno1 - min_resno] = extra;
               }
               distortion_sum[this_resno1 - min_resno] +=
                  factor * dc.geometry_distortion[i].distortion_score * occ1;
            }
         }
      }
   }

   // Make the blocks and plot them.
   //
   //
   int chain_number = chain_id_to_chain_index(dc.chain_id);
   if (chain_number >= 0) {

      for (unsigned int ires=0; ires<distortion_sum.size(); ires++) {
         this_resno1 = resi_of_distortion[ires];
         if (this_resno1 != NO_DISTORTION_IN_THIS_RESIDUE) {
            std::string inscode("");
            std::string at_name = atom_of_distortion[ires];
            std::string altconf("");
            coot::atom_spec_t atom_spec(dc.chain_id, this_resno1, inscode, at_name, altconf);
//             std::cout << "DEBUG:: making block_info with this_resno1 " << this_resno1
//                       << " from ires " << ires << std::endl;
            geometry_graph_block_info block_info(imol, this_resno1, atom_spec,
                                                 0.6 * distortion_sum[ires],
                                                 distortion_string[ires], canvas);


            if (false)
               std::cout << "DEBUG:: plot_block (geom dist) " << block_info.resno
                         << " " << min_resno << " " << offsets[chain_number] << " "
                         << chain_number << std::endl;
            plot_block(block_info, offsets[chain_number], chain_number);

            //          } else {
            //             std::cout << "DEBUG:: Found NO_DISTORTION_IN_THIS_RESIDUE for "
            //                       << ires << std::endl;
         }
      }
   } else {
      // Failing here?  Then check that the chain of these residues
      // was set by a call to render_to_canvas()
      std::cout << "ERROR: failed to get chain in render_blocks_internal" << std::endl;
   }
}

// Note: offset is 0 for a chain starting at 1.
void
coot::geometry_graphs::render_b_factor_blocks(int imol,
                                              int chain_number,
                                              const std::string  &chain_id,
                                              int offset, 
                                              const std::vector<coot::b_factor_block_info_t> &biv) {

   int this_resno;

   if (chain_number < int(chain_index.size()))
      chain_index[chain_number] = chain_id;


   int min_resno = offset + 1;
   int max_resno = offset + biv.size();
   int nres = biv.size();
   if (biv.size() > 0) {
      for (unsigned int i=0; i<biv.size(); i++)
         if (biv[i].resno > max_resno)
            max_resno = biv[i].resno;

      nres = max_resno - min_resno + 1;
      offsets[chain_number] = offset;
      draw_chain_axis(nres, chain_number);
      draw_chain_axis_tick_and_tick_labels(min_resno, max_resno, chain_number);
      blocks[chain_number].resize(nres+1); // needs to index max_resno

      std::string inscode("");
      std::string at_name(" CA ");
      std::string altconf("");
      for (unsigned int i=0; i<biv.size(); i++) {
         this_resno = i + offset;
         coot::atom_spec_t atom_spec(chain_id, biv[i].resno, inscode, biv[i].atom_name, altconf);
         geometry_graph_block_info block_info(imol, biv[i].resno, atom_spec, biv[i].b_factor_var,
                                              biv[i].info_string, canvas);
//          std::cout << "DEBUG:: plot_block (b-factor) " << block_info.resno
//                    << " " << min_resno <<  " " << offset << " " << chain_number << std::endl;
         plot_block(block_info, offset, chain_number);
      }
      label_chain(chain_id, chain_number); // come back to this...
   }
}

// return -1 on "no such chain"
int
coot::geometry_graphs::chain_id_to_chain_index(const std::string &chain_id) const {

   int chain_number = -1;

   // get the chain_number
   for (unsigned int ich=0; ich<chain_index.size(); ich++) {
      if (chain_id == chain_index[ich]) {
         chain_number = ich;
         break;
      }
   }
   return chain_number;
}


void
coot::geometry_graphs::update_residue_blocks(const coot::geometry_distortion_info_container_t &dc) {

   // This can be called with any arrangement of residues, ie. with gaps.
   //

   // OK, so dc has an entry for each restraint, each residue is
   // represented multiple times, we don't want to keep deleting the
   // same residue, so let's make a vector of deleted residues.
   std::vector<coot::residue_spec_t> deleted_block_res_specs;

   int chain_number = chain_id_to_chain_index(dc.chain_id);
   for (unsigned int iblock=0; iblock<dc.geometry_distortion.size(); iblock++) {
      coot::residue_spec_t rs(dc.geometry_distortion[iblock].residue_spec);
      bool ifound = false;
      for (unsigned int i=0; i<deleted_block_res_specs.size(); i++) {
         if (deleted_block_res_specs[i] == rs) {
            ifound = true;
            break;
         }
      }
      if (! ifound) {
         delete_block(chain_number, rs);
         deleted_block_res_specs.push_back(rs);
      }
   }
   render_geometry_distortion_blocks_internal(dc);

}

// This is only to be used on a contiguous range of residues in the
// same chain.
//
void
coot::geometry_graphs::update_residue_blocks_linear(const coot::geometry_distortion_info_container_t &dc) {

   // First, we need to delete the old blocks.
   //
   // We do that by getting the chain id of the residue range then
   // from that find the chain id, and the chain index is part of the
   // indexing of the residue, the other index being the "offsetted"
   // residue number.

   int chain_number = chain_id_to_chain_index(dc.chain_id);

   if (chain_number >= 0) {
      for (int ires=dc.min_resno; ires<=dc.max_resno; ires++) {
         std::cout << "DEBUG:: DELETING block " << ires << std::endl;
         delete_block(chain_number, ires); // delete_block does the offsetting.
      }
   }

   // now render to canvas the new blocks:
   // std::cout << "render_blocks_internal " << dc.min_resno << " " << dc.max_resno << std::endl;
   render_geometry_distortion_blocks_internal_linear(dc, dc.min_resno, dc.max_resno);
}


void
coot::geometry_graphs::update_residue_blocks(const std::vector<coot::geometry_graph_block_info_generic> &dv) {

   // std::cout << "DEBUG:: updating " << dv.size() << " blocks to canvas: " << canvas << std::endl;
   // for (int ires=0; ires<dv.size(); ires++) {
   // std::cout << "DEBUG:: dv[" << ires << "] has resno " << dv[ires].resno << std::endl;
   // }
   for (unsigned int ires=0; ires<dv.size(); ires++) {
      int chain_number = chain_id_to_chain_index(dv[ires].atom_spec.chain_id);
      int resno = dv[ires].resno;
      // std::cout << "DEBUG:: deleting block: " << dv[ires].atom_spec.resno << " "
      // << dv[ires].atom_spec.chain << std::endl;
      delete_block(chain_number, resno); // delete_block does the offsetting.
   }

   // render generic blocks
   for(unsigned int i=0; i<dv.size(); i++) {
      int chain_number = chain_id_to_chain_index(dv[i].atom_spec.chain_id);
      geometry_graph_block_info block_info(imol, dv[i].resno,
                                           dv[i].atom_spec,
                                           dv[i].distortion,
                                           dv[i].distortion_info_string, canvas);
      if (chain_number > -1)
         plot_block(block_info, offsets[chain_number], chain_number);
   }
}

void
coot::geometry_graphs::update_omega_blocks(const coot::omega_distortion_info_container_t &om_dist,
                                           int chain_number,
                                           const std::string &chain_id) {

   int chain_index = chain_id_to_chain_index(chain_id);
   if (chain_index != -1) {

      // delete the omega blocks in om_dist
      for (unsigned int iblock=0; iblock<om_dist.omega_distortions.size(); iblock++) {
         int raw_resno = om_dist.omega_distortions[iblock].resno;
         int offsetted_residue_number = raw_resno - offsets[chain_index];
         if (offsetted_residue_number < int(blocks[chain_index].size())) {

            // add a test here that offsetted_residue_number is
            // sensible for blocks[chain_index].

            if (offsetted_residue_number >=0 &&
                offsetted_residue_number < int(blocks[chain_index].size())) {
               if (blocks[chain_index][offsetted_residue_number]) {
                  if (true) // debug
                     std::cout << ":::: in update_omega_blocks() blocks[" << chain_index
                               << "] is of size(): " << blocks[chain_index].size()
                               << " and deleting block " << offsetted_residue_number
                               << " which is (" << raw_resno << " - " << "offsets["
                               << chain_index << "]=" << offsets[chain_index] << ")"
                               << std::endl;

                  // FIX_THIS
                  // gtk_object_destroy(G_OBJECT(blocks[chain_index][offsetted_residue_number]));

                  blocks[chain_index][offsetted_residue_number] = NULL;
               }
            }
         }
      }

      // a new one
      render_omega_blocks(om_dist, chain_index, chain_id, offsets[chain_index]);
   }
}


void
coot::geometry_graphs::render_omega_blocks(const coot::omega_distortion_info_container_t &om_dist,
                                           int chain_number,
                                           const std::string &chain_id,
                                           int offset_in) {

   int this_resno;

   if (chain_number < int(chain_index.size()))
      chain_index[chain_number] = chain_id;

   int min_resno = offset_in + 1;
   int max_resno = -9999;
   if (om_dist.omega_distortions.size() > 0) {
      for (unsigned int i=0; i<om_dist.omega_distortions.size(); i++)
         if (om_dist.omega_distortions[i].resno > max_resno)
            max_resno = om_dist.omega_distortions[i].resno;
      int nres = max_resno - min_resno + 1;

      offsets[chain_number] = offset_in;
      draw_chain_axis(nres, chain_number);
      draw_chain_axis_tick_and_tick_labels(min_resno, max_resno, chain_number);
      blocks[chain_number].resize(nres+1);

      std::string inscode("");
      std::string at_name(" CA ");
      std::string altconf("");
      for (unsigned int i=0; i<om_dist.omega_distortions.size(); i++) {
         this_resno = om_dist.omega_distortions[i].resno;
         coot::atom_spec_t atom_spec(chain_id, this_resno, inscode, at_name, altconf);
     if (om_dist.omega_distortions[i].info_string.find("Cis") != std::string::npos) {
        if (residue_name(imol, chain_id, this_resno+1, inscode) == "PRO")
           atom_spec.string_user_data = "Cis Peptide Pro";
        else
           atom_spec.string_user_data = "Cis Peptide";
     }
         geometry_graph_block_info block_info(imol,
                                              om_dist.omega_distortions[i].resno,
                                              atom_spec,
                                              om_dist.omega_distortions[i].distortion,
                                              om_dist.omega_distortions[i].info_string,
                                              canvas);
//           std::cout << "DEBUG:: plot_block (omega) resno " << block_info.resno
//                     << " min_resno: " << min_resno <<  " offset " << offset_in
//                    << " chain_number: " << chain_number << std::endl;
         plot_block(block_info, offset_in, chain_number);
      }
      label_chain(chain_id, chain_number);
   }
}

void
coot::geometry_graphs::delete_block(int chain_number, const coot::residue_spec_t &rs) {

   delete_block(chain_number, rs.res_no);

}


void
coot::geometry_graphs::delete_block(int chain_number, int raw_resno) {

   bool debug = 0;
   if (chain_number > -1) {
      if (chain_number < int(blocks.size())) {
         int offsetted_residue_number = raw_resno - offsets[chain_number];
         if (offsetted_residue_number < int(blocks[chain_number].size())) {
            if (debug)
               std::cout << "DEBUB:: destroying block chain_number = " << chain_number << " "
                         << "raw_resno = " << raw_resno << " offsets[chain_number] = "
                         << offsets[chain_number] << " offsetted_residue_number = "
                         << offsetted_residue_number << std::endl;
            if (offsetted_residue_number >= 1) {
               if (offsetted_residue_number < int(blocks[chain_number].size())) {
                  if (blocks[chain_number][offsetted_residue_number]) {
                     // FIX_THIS - Done - keeping for reference.
                     // gtk_object_destroy(GTK_OBJECT(blocks[chain_number][offsetted_residue_number]));
                     GooCanvasItem *item = blocks[chain_number][offsetted_residue_number];
                     GooCanvasItem *parent = goo_canvas_item_get_parent(item);
                     int n_children = goo_canvas_item_get_n_children(parent);
                     for (int i=0; i<n_children; i++) {
                        GooCanvasItem *child = goo_canvas_item_get_child(parent, i);
                        if (item == child) {
                           goo_canvas_item_remove_child(parent, i);
                           break;
                        }
                     }
                  }
               }
            }
         } else {
            std::cout << "ERROR:: Attempt to delete non-existant residue block raw: "
                      << raw_resno <<  " after-offset: " << offsetted_residue_number << " chain-num: "
                      << chain_number << std::endl;
         }
      } else {
         std::cout << "ERROR:: Attempt to delete non-existant residue block in "
                   << "non-existant chain " << chain_number << " " << blocks.size()
                   << std::endl;
      }
   }
}

void
coot::geometry_graphs::delete_block(const std::string &chain_id, int resno) {

   int chain_number = chain_id_to_chain_index(chain_id);
   if (chain_number != -1) {
      delete_block(chain_number, resno);
   }
}

// make this a static of coot::geometry_graphs
//
// static
void
coot::geometry_graphs::density_fit_rescale_button_callback(GtkButton *button, gpointer user_data) {

   coot::geometry_graphs *graph = static_cast<coot::geometry_graphs *>(user_data);
   GtkWidget *entry = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "rescale_entry"));
   const char *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      std::string t(txt);
      float scale = coot::util::string_to_float(t);
      graphics_info_t::residue_density_fit_scale_factor = scale;
      graphics_info_t g;
      std::vector<coot::geometry_graph_block_info_generic> block_set =
         g.density_fit_from_mol(g.molecules[graph->get_imol()].atom_sel,
                                graph->get_imol(),
                                g.Imol_Refinement_Map());
      graph->update_residue_blocks(block_set);
   }
};


void
coot::geometry_graphs::setup_canvas(int n_chains, int max_chain_length) {

   // Fixes: could not find argument "points" in the `GnomeCanvasLine' class ancestry
   // goo_canvas_init();

   GtkWidget *dialog = create_geometry_graphs_dialog();
   canvas = GOO_CANVAS(goo_canvas_new());
   g_object_set(G_OBJECT(canvas), "has-tooltip", TRUE, NULL);

#ifdef WINDOWS_MINGW
   int canvas_usize_x = max_chain_length*10 + 325; // add a bit more to get (may be
                                                   // general for GNOME_CANVAS)
#else
   int canvas_usize_x = max_chain_length*10 + 200; // add a bit to get
                                                  // the label at the
                                                  // right hand side
#endif //MINGW

   // If the canvas is 3500+ or so, then it can't be resized width-wise (Doug Kuntz)
   if (canvas_usize_x > 32100) {
      std::cout << "WARNING:: truncating canvas width! " << std::endl;
      canvas_usize_x = 32100;
   }
   int canvas_usize_y = 80 * n_chains + 30; 
   double scroll_width  = canvas_usize_x + 20.0;
   double scroll_height = canvas_usize_y + 20.0;

   gtk_widget_set_size_request(GTK_WIDGET(canvas), canvas_usize_x, canvas_usize_y);
   gtk_widget_set_size_request(dialog, 600, 400);
   GtkWidget *scrolled_window = lookup_widget(dialog, "geometry_graphs_scrolledwindow");

//    std::cout << "INFO:: canvas size based on " << n_chains << " chains with max"
//              << " length " << max_chain_length << std::endl;
//    std::cout << "INFO:: canvas size: " << canvas_usize_x << " " << canvas_usize_y
//              << std::endl;

   gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
                                         GTK_WIDGET(canvas));

   // this doesn't work: perhaps because of this? gtk_scrolled_window_add: assertion 'child_widget == NULL' failed
   // gtk_container_add(GTK_CONTAINER(scrolled_window), GTK_WIDGET(canvas));

   
   double left_limit    = 0;
   double upper_limit   = 0;

   std::cout << "set scroll region here " << std::endl;
   // goo_canvas_set_scroll_region(canvas, left_limit, upper_limit, scroll_width, scroll_height);

   // gtk_canvas_set_pixels_per_unit(GTK_CANVAS(canvas),zoom);


   gtk_widget_show(GTK_WIDGET(canvas));
   gtk_widget_show(dialog);

   // gtk_widget_ref(GTK_WIDGET(canvas));
   g_object_set_data(G_OBJECT(dialog), "geometry_graph_canvas", canvas);
   g_object_set_data(G_OBJECT(canvas), "geometry-graph", this); // looked up by geometry_graph_dialog_to_object()

   if (graph_type == GEOMETRY_GRAPH_DENSITY_FIT) {
      GtkWidget *dialog_vbox = lookup_widget(dialog, "geometry_graphs_dialog_vbox");
      if (dialog_vbox) {

         GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
         GtkWidget *button = gtk_button_new_with_label("Rescale");
         GtkWidget *entry = gtk_entry_new();
         gtk_entry_set_text(GTK_ENTRY(entry), "1.0");
         gtk_box_pack_start (GTK_BOX (vbox), entry,  FALSE, FALSE, 3);
         gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 3);
         gtk_widget_show(entry);
         gtk_widget_show(button);
         gtk_widget_show(vbox);
         gtk_widget_set_size_request(entry, 90, -1);
         gtk_box_pack_start(GTK_BOX(dialog_vbox), vbox, FALSE, FALSE, 3);
         g_object_set_data(G_OBJECT(button), "rescale_entry", entry);
         g_signal_connect(G_OBJECT(button), "clicked",
                          G_CALLBACK(density_fit_rescale_button_callback),
                          this);
      }
   }

}

//  0,0 -------> +x
//   |
//   |
//   |
//   |
//   |
//   \/
//   +y
//
// y2 needs to be greater than y1 (and same for the xs of course)
//
void
coot::geometry_graphs::plot_block(const coot::geometry_graph_block_info &block_info,
                                  int offset_in, int chain_number) {

   // when updating, it would be better to set g_object_set_data to resize it
   // rather than delete an create a new one.

   std::string outline_color = "grey10";
   int resno = block_info.resno;
   GooCanvasItem *item;
   double scale = 10.0;
   double chain_scale = 80.0;
   double distortion = block_info.distortion;
   int offset = offset_in;

   // x2 and y2 are height and width
   double x1 = (resno-offset)*scale;
   double x2 = 9.0;
   double y1 = 50 + chain_number*chain_scale;
   double y2 = -distortion * 0.5;

   std::string colour = distortion_to_colour(distortion);

   // Lets make the blocks with missing atoms a different colour.
   //

   if (block_info.atom_spec.string_user_data == "Missing Atoms")
      colour = "slate blue";

   if (block_info.atom_spec.string_user_data == "Cis Peptide Pro"
       && graphics_info_t::mark_cis_peptides_as_bad_flag)
      colour = "slate blue";


   if (false) {
      std::cout << " distortion: " << distortion << std::endl;
      std::cout << " block: " << x1 << " " << x2 << " " << y1 << " " << y2 << " fill " << colour << std::endl;
   }

   std::string tooltip = block_info.distortion_info_string;
   item = goo_canvas_rect_new(goo_canvas_get_root_item(canvas),
                              x1, y1, x2, y2,
                              "fill_color", colour.c_str(),
                              "stroke-color", outline_color.c_str(),
                              "tooltip", tooltip.c_str(),
                              NULL);

//    std::cout << "plotting block[" << chain_number << "][" << resno-offset<< "]" << std::endl;
//    std::cout << "plotting " << resno << " " << offset << " " << chain_index.size() << std::endl;
//    std::cout << "plotting " << blocks[chain_number].size() << std::endl;

   if (resno-offset < int(blocks[chain_number].size())) {
      blocks[chain_number][resno-offset] = item;
   } else {
      std::cout << "ERROR:: storing plot block[" << chain_number << "][" << resno-offset<< "]" << std::endl;
      std::cout << "ERROR:: storing plot resno: " << resno << " offset: " << offset << " chain_index.size(): "
                << chain_index.size() << std::endl;
      std::cout << "ERROR:: storing plot: blocks[chain_number].size(): " << blocks[chain_number].size()
                << std::endl;
   }

   coot::geometry_graph_block_info *local_block_info_p = new coot::geometry_graph_block_info(block_info);

   g_object_set_data(G_OBJECT(item), "geometry-graph", this);

   g_signal_connect(G_OBJECT(item), "button_press_event",
                    G_CALLBACK(coot::on_geometry_graph_block_clicked),
                    gpointer(local_block_info_p));

}

void
coot::geometry_graphs::draw_chain_axis(int nres, int ichain) const {

   double chain_scale = 80.0;
   double scale = 10.0;

   //   x2,y2
   //  |
   //  |x0,y0
   //  ------------------------------- x1,y1
   double x0 = 7;
   double y0 = 55 + ichain*chain_scale;

   double x1 = x0 + scale * nres + 10;
   double y1 = y0;

   double x2 = x0;
   double y2 = y0 - 50;

//    std::cout << "chain axis: " << ichain << " "
//              << "(" << x0 << "," << y0 << ")"
//              << "(" << x1 << "," << y1 << ")"
//              << "(" << x2 << "," << y2 << ")"
//              << std::endl;

   GooCanvasItem *item;
   GooCanvasPoints *points = goo_canvas_points_new(3);

   points->coords[0] = x1;
   points->coords[1] = y1;

   points->coords[2] = x0;
   points->coords[3] = y0;

   points->coords[4] = x2;
   points->coords[5] = y2;

   item = goo_canvas_polyline_new(goo_canvas_get_root_item(canvas),
                                  TRUE, 0,
                                  "width_pixels", 3.0,
                                  "points", points,
                                  "fill_color", "black",
                                  NULL);
   goo_canvas_points_unref(points);

}


void
coot::geometry_graphs::draw_chain_axis_tick_and_tick_labels(int min_resno,
                                                            int max_resno,
                                                            int chain_number) const {

   GooCanvasItem *item;
   GooCanvasPoints *points;

   double res_scale = 10.0;
   double chain_scale = 80.0;

   std::string tick_colour = "grey80";

   int tick_start_res_no = min_resno + 9; // make this more clever.

   for (int i=tick_start_res_no; i<max_resno; i+= 10) {

      points = goo_canvas_points_new(2);
      double x0 = (i - min_resno + 1) * res_scale + 5.0;
      double y0 = 55.0 + chain_number * chain_scale;
      double x1 = x0;
      double y1 = y0 + 5.0;

      points->coords[0] = x0;
      points->coords[1] = y0;

      points->coords[2] = x1;
      points->coords[3] = y1;

      item = goo_canvas_polyline_new(goo_canvas_get_root_item(canvas),
                                     TRUE, 0,
                                     "width_pixels", 2.0,
                                     "points", points,
                                     "fill_color", tick_colour.c_str(),
                                     NULL);
      goo_canvas_points_unref(points);

      item = goo_canvas_text_new(goo_canvas_get_root_item(canvas),
                                 std::to_string(i).c_str(),
                                 x0 - 6.0,
                                 y1 + 4.0,
                                 -1,
                                 GOO_CANVAS_ANCHOR_WEST,
                                 "font", fixed_font_str.c_str(),
                                 "fill_color", tick_colour.c_str(),
                                 NULL);
   }

}



void
coot::geometry_graphs::label_chain(const std::string &label, int ichain) const {

   double chain_scale = 80.0;

   double x = 20;
   double y = 10 + ichain*chain_scale;
   std::string text = "Chain ";
   text += label;
   GooCanvasItem *item;

   item = goo_canvas_text_new(goo_canvas_get_root_item(canvas),
                              text.c_str(),
                               x, y,
                              -1,
                              GOO_CANVAS_ANCHOR_WEST,
                              "font", fixed_font_str.c_str(),
                              "fill_color", "grey90",
                              NULL);
}


void
coot::geometry_graphs::setup_internal() {

   distortion_max = 100.0; // guess

   colour_list.push_back("green1");
   colour_list.push_back("chartreuse1");
   colour_list.push_back("OliveDrab1");
   colour_list.push_back("DarkOliveGreen1");
   colour_list.push_back("yellow2");
   colour_list.push_back("gold1");
   colour_list.push_back("goldenrod1");
   colour_list.push_back("tan1");
   colour_list.push_back("orange1");
   colour_list.push_back("coral1");
   colour_list.push_back("tomato1");
   colour_list.push_back("OrangeRed1");
   colour_list.push_back("firebrick1");
   colour_list.push_back("red");

   fixed_font_str = "fixed";
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   fixed_font_str = "monospace";
#endif

   tooltip_item = NULL;
   tooltip_item_text = NULL;


}


std::string
coot::geometry_graphs::distortion_to_colour(const double &distortion) const {

   if (distortion >= distortion_max)
      return colour_list.back();

   if (distortion < 0) // shouldn't be any more - now that I've remove
                       // negative occs (which multiplied the distortions).
      return colour_list[0];

   int ind = int(double(colour_list.size()) * distortion/distortion_max);
   if (ind == int(colour_list.size()))
      ind = colour_list.size() - 1;

   return colour_list[ind];

}


GtkWidget *
coot::geometry_graphs::dialog() const {
   return lookup_widget(GTK_WIDGET(canvas), "geometry_graphs_dialog");
}



void
coot::geometry_graphs::tooltip_like_box(const geometry_graph_block_info &bi,
                                        const GdkEvent *event) {

   clear_tooltip_box();

   std::string label = bi.distortion_info_string;
// BL says:: think the box should be larger
// think that is for Gnome_canvas build
#ifdef HAVE_GNOME_CANVAS
   double tw = label.size() * 8.0 + 10.0;
#else
   double tw = label.size() * 6.0 + 10.0;
#endif // GNOME_CANVAS

   double x1 = 0.0, y1 = 0.0;
   // int x_as_int, y_as_int;
   // GdkModifierType state;

   // gdk_window_get_pointer(GDK_WINDOW(dialog()), &x_as_int, &y_as_int, &state);

   if (event->type == GDK_MOTION_NOTIFY) {
      GdkEventMotion mymotion = event->motion;
//       std::cout << "my motion: " <<  event->type << " "
//     << mymotion.x << "  " << mymotion.y << std::endl;
      x1 = 10.0;
      y1 = 10.0;

      // 0.21 is good
      // 0.41 is good
      // 0.61 is good
      // 0.66 is bad
      // 0.71 is bad
      // 0.81 is bad
      //
      // Argh, I give up.  Time to Ask Kevin.
      //
      double xtmp = double(mymotion.x);
      double ytmp = double(mymotion.y);

      // magic voodoo: who knows why this works... some sort of
      // integer thing?
      int ix = int(xtmp/10.0);
      int iy = int(ytmp/10.0);
      x1 += 10.0 * ix;
      y1 += 10.0 * iy;

      // std::cout << "at " << x1 << " " << y1 << "   " << ix << " " << iy << std::endl;

      double xt = x1 + 5.0;
      double yt = y1 + 7.0;
      double x2 = x1 + tw;
      double y2 = y1 + 16.0;

      tooltip_item = goo_canvas_rect_new(goo_canvas_get_root_item(canvas),
                                         x1,
                                         y1,
                                         x2,
                                         y2,
                                         "fill_color", "PaleGreen",
                                         "outline_color", "black",
                                         NULL); 

      tooltip_item_text = goo_canvas_text_new(goo_canvas_get_root_item(canvas),
                                              label.c_str(),
                                              xt, yt,
                                              -1,
                                              GOO_CANVAS_ANCHOR_WEST,
                                              "font", fixed_font_str.c_str(),
                                              "fill_color", "black",
                                              NULL);
   }
}

void
coot::geometry_graphs::clear_tooltip_box() {

   
   // if (tooltip_item)
   //    g_object_destroy(G_OBJECT(tooltip_item));
   // if (tooltip_item_text)
   //    g_object_destroy(G_OBJECT(tooltip_item_text));
   // tooltip_item = NULL;
   // tooltip_item_text = NULL;

   std::cout << "This is not the way to do tooltips " << std::endl;
}



void
coot::geometry_graphs::render_to_canvas(const std::vector<coot::geometry_graph_block_info_generic> &gbi,
                                        int chain_number,
                                        const std::string &chain_id,
                                        int max_resno,
                                        int min_resno,
                                        int offset) {

//    std::cout << "render_to_canvas with offset "  << offset
//              << " max_resno " << max_resno << " min_resno " << min_resno << std::endl;

   // int this_resno;

   if (chain_number < int(chain_index.size()))
      chain_index[chain_number] = chain_id;

   int nres = max_resno - min_resno + 1;
   offsets[chain_number] = offset;
   draw_chain_axis(nres, chain_number);
//    std::cout << "DEBUG:: min_resno: " << min_resno << " max_resno: " << max_resno
//              << " chain_number " << chain_number << std::endl;
   draw_chain_axis_tick_and_tick_labels(min_resno, max_resno, chain_number);
   blocks[chain_number].resize(nres + 1); // needs to index max_resno

   for(unsigned int i=0; i<gbi.size(); i++) {
      geometry_graph_block_info block_info(imol, gbi[i].resno, gbi[i].atom_spec, gbi[i].distortion,
                                           gbi[i].distortion_info_string, canvas);
//        std::cout << "DEBUG:: plot_block (generic) " << block_info.resno
//                   << " " << min_resno <<  " " << offset << " " << chain_number << std::endl;
      plot_block(block_info, offset, chain_number);
   }
   label_chain(chain_id, chain_number);
}


void
coot::geometry_graphs::close_yourself() {

   // free up blocks (and ticks?)

   for (unsigned int ichain=0; ichain<blocks.size(); ichain++) {
      for (unsigned int ires=0; ires<blocks[ichain].size(); ires++) {
         if (blocks[ichain][ires]) {
            // gtk_object_destroy(GTK_OBJECT(blocks[ichain][ires]));
            // g_object_destroy(G_OBJECT(blocks[ichain][ires]));
            std::cout << "destroy block here A " << blocks[ichain][ires] << std::endl;
            blocks[ichain][ires] = NULL;
         }
      }
   }
   // The actual destroy is happening in the OK button callback.  This
   // is part of the destroy callback of the window which is called by
   // that callback.
}
