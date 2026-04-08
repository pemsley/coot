#ifdef DO_GEOMETRY_GRAPHS
   coot::geometry_graphs *geometry_graph_dialog_to_object(GtkWidget *w) const {
      coot::geometry_graphs *gr = NULL;
      if (!w) {
         std::cout << "geometry_graph_dialog_to_object case A" << std::endl;
	 std::cout << "ERROR:: null w in geometry_graph_dialog_to_object" << std::endl;
      } else {
         std::cout << "geometry_graph_dialog_to_object case B" << std::endl;
	 GtkWidget *local_graph_canvas = lookup_widget(w, "geometry_graph_canvas");
	 if (local_graph_canvas) {
	    // gr = (coot::geometry_graphs *) (gtk_object_get_user_data(GTK_OBJECT(local_graph_dialog)));
            // this need a corresponding change to set the g_object data - whereever that is.
            gr = static_cast<coot::geometry_graphs *> (g_object_get_data(G_OBJECT(local_graph_canvas), "geometry-graph"));
            if (! gr)
               std::cout << "geometry_graph_dialog_to_object case C - bad news" << std::endl;
         } else {
            std::cout << "geometry_graph_dialog_to_object case D - bad news" << std::endl;
         }
      }
      return gr;
   }
#endif
