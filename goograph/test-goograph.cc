

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
   data.push_back(std::pair<double, double> (0.5, 0.01));

    for (unsigned int i=0; i<data.size(); i++) { 
       data[i].first *= 100;
       data[i].second *= 10;
       data[i].second += 5.5;
    }
   
   coot::goograph g;
   int trace = g.trace_new();
   g.set_plot_title("Test graph");
   g.set_data(trace, data);
   g.set_ticks(coot::goograph::X_AXIS, 0.5, 0.1);
   g.set_ticks(coot::goograph::Y_AXIS, 0.1, 0.02);
   // g.set_extents(coot::goograph::X_AXIS, 0.0,   4.0);
   // g.set_extents(coot::goograph::Y_AXIS, 0.0, 0.6);

   if (1) {
      g.set_ticks(  coot::goograph::X_AXIS, 50,  10);
      g.set_ticks(  coot::goograph::Y_AXIS, 1,   0.2);
      // g.set_extents(coot::goograph::X_AXIS, 40.0, 400.0);
      // g.set_extents(coot::goograph::Y_AXIS, 0.0, 60);
   } 
   
   g.set_axis_label(coot::goograph::X_AXIS, "Bond length");
   g.set_axis_label(coot::goograph::Y_AXIS, "Counts");
   g.plot(trace, coot::goograph::PLOT_TYPE_BAR);
   // g.plot(trace, coot::goograph::PLOT_TYPE_LINE);
   // g.plot(trace, coot::goograph::PLOT_TYPE_SMOOTHED_LINE);

   double f = 0.01;

   bool do_annotations = false; 
   if (do_annotations) { 
      lig_build::pos_t p1(150*f,  0*f);
      lig_build::pos_t p2(150*f, 50*f);
      lig_build::pos_t p3(195*f, 47*f);
      lig_build::pos_t p4(152*f, 47*f);
      lig_build::pos_t p5(230*f, 47*f);
      lig_build::pos_t p6(200*f, 35*f);
      lig_build::pos_t p7(200*f,  0*f);
      lig_build::pos_t p8(202*f, 30*f);
      lig_build::pos_t p9(240*f, 30*f);
      lig_build::pos_t p10(270*f, 30*f);
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
