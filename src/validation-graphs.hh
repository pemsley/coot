
#ifndef VALIDATION_GRAPHS_HH
#define VALIDATION_GRAPHS_HH


#include <gtk/gtk.h>

namespace coot {

   class validation_graphs_t {
   public:
      GtkWidget *dynarama_is_displayed;
      GtkWidget *sequence_view_is_displayed;
      GtkWidget *geometry_graph;
      GtkWidget *b_factor_variance_graph;
////B B GRAPH
      GtkWidget *b_factor_graph;
////E B GRAPH
      GtkWidget *residue_density_fit_graph;
      GtkWidget *omega_distortion_graph;
      GtkWidget *rotamer_graph; 
      GtkWidget *ncs_diffs_graph;
      validation_graphs_t() {
	 init();
      }
      
      void init() {
	 dynarama_is_displayed      = NULL;
	 sequence_view_is_displayed = NULL;
	 geometry_graph             = NULL;
	 b_factor_variance_graph    = NULL;
////B B GRAPH
	 b_factor_graph  	    = NULL;
////E B GRAPH
	 residue_density_fit_graph  = NULL;
	 omega_distortion_graph     = NULL;
	 rotamer_graph              = NULL;
	 ncs_diffs_graph            = NULL;
      }
   };
}

#endif // VALIDATION_GRAPHS_HH

