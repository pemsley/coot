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
   data.push_back(std::pair<double, double> (1.0, 0.1));
   data.push_back(std::pair<double, double> (1.1, 0.2));
   data.push_back(std::pair<double, double> (1.2, 0.25));
   data.push_back(std::pair<double, double> (1.3, 0.32));
   data.push_back(std::pair<double, double> (1.4, 0.37));
   data.push_back(std::pair<double, double> (1.5, 0.44));
   data.push_back(std::pair<double, double> (1.6, 0.42));
   data.push_back(std::pair<double, double> (0.5, 0.01));
   data.push_back(std::pair<double, double> (0.434, 0.01));
   data.push_back(std::pair<double, double> (1.7, 0.32));
   data.push_back(std::pair<double, double> (1.8, 0.3));
   data.push_back(std::pair<double, double> (1.9, 0.2));
   data.push_back(std::pair<double, double> (2.0, 0.12));
   data.push_back(std::pair<double, double> (2.1, 0.1));
   data.push_back(std::pair<double, double> (2.2, 0.08));
   data.push_back(std::pair<double, double> (2.3, 0.05));
   data.push_back(std::pair<double, double> (2.4, 0.01));
   data.push_back(std::pair<double, double> (2.5, 0.0 ));
   data.push_back(std::pair<double, double> (2.9, 0.0 ));
   data.push_back(std::pair<double, double> (3.0, 0.01));
   data.push_back(std::pair<double, double> (3.1, 0.03));
   data.push_back(std::pair<double, double> (3.2, 0.06));
   data.push_back(std::pair<double, double> (3.3, 0.05));
   data.push_back(std::pair<double, double> (3.4, 0.02));
   data.push_back(std::pair<double, double> (3.5, 0.01));

    for (unsigned int i=0; i<data.size(); i++) { 
       data[i].first *= 87;
       data[i].second *= 83;
       data[i].second += 0.6;
    }
   
   coot::goograph g;
   int trace = g.trace_new();
   g.set_plot_title("Bond length distribution vs. database");
   g.set_data(trace, data);
   
   g.set_axis_label(coot::goograph::X_AXIS, "Bond length");
   g.set_axis_label(coot::goograph::Y_AXIS, "Counts");
   g.plot(trace, coot::goograph::PLOT_TYPE_BAR, "");

   double f = 1;
   double fy = 1;

   for (unsigned int i=0; i<data.size(); i++) { 
      data[i].first += 80;
      data[i].second *= 0.6;
      data[i].second += 0.5;
   }
   bool dashed = true;
   trace = g.trace_new();
   g.set_data(trace, data);
   std::string colour = "blue";
   g.plot(trace, coot::goograph::PLOT_TYPE_SMOOTHED_LINE, colour, dashed);

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
