

#include "goograph/goograph.hh"

int main (int argc, char **argv) {

   gtk_init (&argc, &argv);
   std::vector<std::pair<double, double> > data;
   data.push_back(std::pair<double, double> (1.0, 0.1));
   data.push_back(std::pair<double, double> (1.1, 0.1));
   data.push_back(std::pair<double, double> (1.2, 0.25));
   data.push_back(std::pair<double, double> (1.3, 0.32));
   data.push_back(std::pair<double, double> (1.4, 0.37));
   data.push_back(std::pair<double, double> (1.5, 0.47));
   data.push_back(std::pair<double, double> (1.6, 0.42));
   data.push_back(std::pair<double, double> (1.7, 0.32));
   data.push_back(std::pair<double, double> (1.8, 0.3));
   data.push_back(std::pair<double, double> (1.9, 0.2));
   data.push_back(std::pair<double, double> (2.0, 0.12));
   data.push_back(std::pair<double, double> (2.1, 0.1));
   data.push_back(std::pair<double, double> (2.2, 0.08));
   data.push_back(std::pair<double, double> (2.3, 0.05));
   data.push_back(std::pair<double, double> (2.4, 0.01));
   data.push_back(std::pair<double, double> (3.0, 0.01));
   data.push_back(std::pair<double, double> (3.1, 0.03));
   data.push_back(std::pair<double, double> (3.2, 0.06));
   data.push_back(std::pair<double, double> (3.3, 0.05));
   data.push_back(std::pair<double, double> (3.4, 0.02));
   data.push_back(std::pair<double, double> (3.5, 0.01));

   for (unsigned int i=0; i<data.size(); i++) { 
      data[i].first *= 100;
      data[i].second *= 100;
   }
   
   coot::goograph g;
   int trace = g.trace_new();
   g.set_plot_title("Test graph");
   g.set_data(trace, data);
   g.set_ticks(coot::goograph::X_AXIS,  0.5, 0.1);
   g.set_ticks(coot::goograph::Y_AXIS, 0.1, 0.02);
   // g.set_extents(coot::goograph::X_AXIS, 0.0,   4.0);
   // g.set_extents(coot::goograph::Y_AXIS, 0.0, 0.6);

   if (1) {
      g.set_ticks(  coot::goograph::X_AXIS, 50,  10);
      g.set_ticks(  coot::goograph::Y_AXIS, 10,   2);
      g.set_extents(coot::goograph::X_AXIS, 40.0, 400.0);
      g.set_extents(coot::goograph::Y_AXIS, 0.0, 60);
   } 
   
   g.set_axis_label(coot::goograph::X_AXIS, "Bond length");
   g.set_axis_label(coot::goograph::Y_AXIS, "Counts");
   // g.plot(trace, coot::goograph::PLOT_TYPE_BAR);
   g.plot(trace, coot::goograph::PLOT_TYPE_LINE);
   g.show_dialog();

   gtk_main();
   return 0;
}
