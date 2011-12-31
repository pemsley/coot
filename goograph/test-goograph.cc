

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
       data[i].first *= 60;
       data[i].second *= 4;
       data[i].second += 0.3;
    }
   
   coot::goograph g;
   int trace = g.trace_new();
   g.set_plot_title("Test graph");
   g.set_data(trace, data);
   // g.set_ticks(  coot::goograph::X_AXIS, 50,  10);
   // g.set_ticks(  coot::goograph::Y_AXIS, 10,   2);
   // g.set_extents(coot::goograph::X_AXIS, 0.0,   4.0);
   // g.set_extents(coot::goograph::Y_AXIS, 0.0, 0.6);

   
   g.set_axis_label(coot::goograph::X_AXIS, "Bond length");
   g.set_axis_label(coot::goograph::Y_AXIS, "Counts");
   g.plot(trace, coot::goograph::PLOT_TYPE_BAR);
   // g.plot(trace, coot::goograph::PLOT_TYPE_LINE);
   // g.plot(trace, coot::goograph::PLOT_TYPE_SMOOTHED_LINE);

   double f = 1;
   double fy = 0.04;

   for (unsigned int i=0; i<data.size(); i++) { 
      data[i].first += 50;
      data[i].second *= 0.6;
      data[i].second += 0.5;
   }
    trace = g.trace_new();
    g.set_data(trace, data);
    g.plot(trace, coot::goograph::PLOT_TYPE_SMOOTHED_LINE);

   bool do_annotations = false;
   // do_annotations = true; 
   if (do_annotations) {
      // red lines
      lig_build::pos_t p1(150*f, 12*fy);
      lig_build::pos_t p2(150*f, 36*fy);
      lig_build::pos_t p3(195*f, 35*fy);
      lig_build::pos_t p4(152*f, 35*fy);
      // red text
      lig_build::pos_t p5(230*f, 35*fy);

      // black lines
      lig_build::pos_t p6(200*f, 33*fy);
      lig_build::pos_t p7(200*f, 12*fy);
      lig_build::pos_t p8(202*f, 28*fy);
      lig_build::pos_t p9(250*f, 28*fy);
      // black text
      lig_build::pos_t p10(280*f, 28*fy);
      
      bool dashed = true;
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
