

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
      wlig = wlig_in;
   }
};

ligand_wiggly_ligand_data_t setup_ligands_progress_bar(); // return the progress bar
void setup_ligands_progress_bar_idle(coot::wligand *wlig,
				     int imol_ligand,
				     ligand_wiggly_ligand_data_t ld);

gboolean install_simple_wiggly_ligand_idle_fn(gpointer data);

#endif // C_INTERFACE_LIGANDS_WIDGETS_HH
