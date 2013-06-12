/* src/test-goograph.cc
 * 
 * Copyright 2011, 2012 by The University of Oxford
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


#include "goograph.hh"

int main (int argc, char **argv) {

   gtk_init (&argc, &argv);
   std::vector<std::pair<double, double> > data;

   data.push_back(std::pair<double, double> ( 104.5, 4));
   data.push_back(std::pair<double, double> ( 104.75, 1));
   data.push_back(std::pair<double, double> ( 105, 2));
   data.push_back(std::pair<double, double> ( 105.25, 1));
   data.push_back(std::pair<double, double> ( 105.5, 0));
   data.push_back(std::pair<double, double> ( 105.75, 0));
   data.push_back(std::pair<double, double> ( 106, 0));
   data.push_back(std::pair<double, double> ( 106.25, 0));
   data.push_back(std::pair<double, double> ( 106.5, 0));
   data.push_back(std::pair<double, double> ( 106.75, 0));
   data.push_back(std::pair<double, double> ( 107, 0));
   data.push_back(std::pair<double, double> ( 107.25, 0));
   data.push_back(std::pair<double, double> ( 107.5, 0));
   data.push_back(std::pair<double, double> ( 107.75, 0));
   data.push_back(std::pair<double, double> ( 108, 0));
   data.push_back(std::pair<double, double> ( 108.25, 0));
   data.push_back(std::pair<double, double> ( 108.5, 0));
   data.push_back(std::pair<double, double> ( 108.75, 0));
   data.push_back(std::pair<double, double> ( 109, 0));
   data.push_back(std::pair<double, double> ( 109.25, 0));
   data.push_back(std::pair<double, double> ( 109.5, 1));
   data.push_back(std::pair<double, double> ( 109.75, 1));
   data.push_back(std::pair<double, double> ( 110, 0));
   data.push_back(std::pair<double, double> ( 110.25, 1));
   data.push_back(std::pair<double, double> ( 110.5, 3));
   data.push_back(std::pair<double, double> ( 110.75, 2));


//     for (unsigned int i=0; i<data.size(); i++) { 
//        data[i].first *= 87;
//        data[i].second *= 83;
//        data[i].second += 0.6;
//     }
   
   coot::goograph g;
   int trace = g.trace_new();
   g.set_plot_title("Angle distribution");
   g.set_data(trace, data);
   
   g.set_axis_label(coot::goograph::X_AXIS, "Angle");
   g.set_axis_label(coot::goograph::Y_AXIS, "Counts");
   g.set_trace_type(trace, coot::graph_trace_info_t::PLOT_TYPE_BAR);

   double f = 1;
   double fy = 1;

//    for (unsigned int i=0; i<data.size(); i++) { 
//       data[i].first += 80;
//       data[i].second *= 0.6;
//       data[i].second += 0.5;
//    }
   
   bool dashed = true;
   trace = g.trace_new();
   g.set_data(trace, data);
   std::string colour = "blue";
   g.set_trace_type(trace, coot::graph_trace_info_t::PLOT_TYPE_SMOOTHED_LINE, dashed);
   g.set_trace_colour(trace, colour);

   bool do_annotations = false;
   do_annotations = true; 
   if (do_annotations) {
      // red lines
      lig_build::pos_t p1(150*f, 0*fy);
      lig_build::pos_t p2(150*f, 40*fy);
      lig_build::pos_t p3(195*f, 35*fy);
      lig_build::pos_t p4(152*f, 35*fy);
      // red text
      lig_build::pos_t p5(230*f, 35*fy);

      // black lines
      lig_build::pos_t p6(211*f, 30*fy);
      lig_build::pos_t p7(211*f, 0*fy);
      lig_build::pos_t p8(213*f, 28*fy);
      lig_build::pos_t p9(250*f, 28*fy);
      // black text
      lig_build::pos_t p10(280*f, 28*fy);
      
      g.add_annotation_line(p1, p2, "#aa0000", 3, dashed, false, false);
      g.add_annotation_line(p3, p4, "#aa0000", 2, dashed, false, true);
      g.add_annotation_text("Median 4.52", p5, "#aa0000", "");
      g.add_annotation_line(p6, p7, "black", 2, dashed, false, false);
      g.add_annotation_line(p8, p9, "black", 2, dashed, true, false);
      g.add_annotation_text("Model 5.3", p10, "black", "");
   }
   g.show_dialog();

   gtk_main();
   return 0;
}
