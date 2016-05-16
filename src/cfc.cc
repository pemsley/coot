
// header-here

#ifdef USE_PYTHON


// needed to parse cc-interface.hh
#ifdef USE_GUILE
#include <libguile.h> 
#endif // USE_GUILE


#include "cfc.hh"
#include "cfc-widgets.hh"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/chemical-feature-clusters.hh"
#include "cc-interface.hh"
#include "graphics-info.h"
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#include "interface.h" // because we use create_cfc_dialog();
#include "c-interface-python.hh" // because we use display_python().

// return a Python object that contains (with indexing)
//   water positions around ligand
//   chemical features of ligands
//
// environment_residues_py is a list of residue specs
// solvated_ligand_info_py is a list of
//   list mol_no ligand_spec
// which, with radius will be used to find the waters
// 
PyObject *chemical_feature_clusters_py(PyObject *environment_residues_py,
				       PyObject *solvated_ligand_info_py,
				       float radius) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   PyObject *r = Py_False;

   if (PyList_Check(environment_residues_py)) {
      if (PyList_Check(solvated_ligand_info_py)) {

	 std::vector<coot::residue_spec_t> neighbs =
	    py_to_residue_specs(environment_residues_py);

	 std::vector<coot::chem_feat_solvated_ligand_spec> ligands; // fill this from
         	                                           // solvated_ligand_info_py

	 int n = PyObject_Length(solvated_ligand_info_py);
	 for(int i=0; i<n; i++) {

	    PyObject *o = PyList_GetItem(solvated_ligand_info_py, i);

	    // o should be a list of molecule-idx, residue-spec
	    if (PyList_Check(o)) {
	       int no = PyObject_Length(o);
	       if (no == 2) {
		  PyObject *mol_idx_py  = PyList_GetItem(o, 0);
		  PyObject *lig_spec_py = PyList_GetItem(o, 1);

		  if (PyInt_Check(mol_idx_py)) {
		     int imol = PyInt_AsLong(mol_idx_py);
		     coot::residue_spec_t ligand_spec = residue_spec_from_py(lig_spec_py);

		     if (graphics_info_t::is_valid_model_molecule(imol)) {

			mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;

			std::vector<coot::residue_spec_t> neighbs_waters; // fill this
			std::vector<coot::residue_spec_t> neighbs_raw =
			   coot::residues_near_residue(ligand_spec, mol, radius);
			for (unsigned int i=0; i<neighbs_raw.size(); i++) { 
			   mmdb::Residue *res = coot::util::get_residue(neighbs_raw[i], mol);
			   if (res) {
			      std::string res_name = res->GetResName();
			      if (res_name == "HOH")
				 neighbs_waters.push_back(neighbs_raw[i]);
			   }
			}

			if (false) { // debug
			   std::cout << "molecule imol: " << imol << " has "
				     << neighbs_waters.size() << " waters" << std::endl;
			   for (unsigned int iwat=0; iwat<neighbs_waters.size(); iwat++)
			      std::cout << "   " << iwat << " " << neighbs_waters[iwat]
					<< std::endl;
			}

			coot::chem_feat_solvated_ligand_spec lig(ligand_spec, neighbs_waters, mol);
			ligands.push_back(lig);
			
		     } else {
			std::cout << "ERROR:: invalid model molecule " << imol << std::endl;
		     }
		     
		  } else {
		     std::cout << "ERROR:: mol_idx is not a int " << std::endl;
		     break;
		  }
	       }
	    }
	 }

	 coot::chem_feat_clust cl(neighbs, ligands, graphics_info_t::Geom_p());

	 std::vector<coot::chem_feat_clust::water_attribs> water_positions =
	    cl.get_water_positions();

	 std::cout << "INFO:: Got " << water_positions.size() << " waters" << std::endl;

	 PyObject *water_attribs_py = PyList_New(water_positions.size());
	 for (unsigned int iw=0; iw<water_positions.size(); iw++) {
	    PyObject *o = PyList_New(3);
	    PyObject *pos_py = PyList_New(3);
	    PyList_SetItem(pos_py, 0, PyFloat_FromDouble(water_positions[iw].pos.x()));
	    PyList_SetItem(pos_py, 1, PyFloat_FromDouble(water_positions[iw].pos.y()));
	    PyList_SetItem(pos_py, 2, PyFloat_FromDouble(water_positions[iw].pos.z()));
	    PyList_SetItem(o, 0, PyInt_FromLong(water_positions[iw].ligand_idx));
	    PyList_SetItem(o, 1, py_residue(water_positions[iw].residue_spec()));
	    PyList_SetItem(o, 2, pos_py);
	    PyList_SetItem(water_attribs_py, iw, o);
	 }

	 r = water_attribs_py;

      } else {
	 std::cout << "ERROR:: chemical_feature_clusters_py() arg 2 is not a list"  << std::endl;
      }
   } else {
      std::cout << "ERROR:: chemical_feature_clusters_py() arg 1 is not a list"  << std::endl;
   } 

   
#endif // MAKE_ENHANCED_LIGAND_TOOLS   

   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;
}

void chemical_feature_clusters_accept_info_py(PyObject *env_residue,
					      PyObject *mol_ligand_specs,
					      PyObject *cluster_info) {

   if (graphics_info_t::use_graphics_interface_flag) {


      // I want a vector of things that contain
      // the molecule number
      // the water residue spec
      // the cluster number
      //
      // and a vector of cluster specs
      //
      cfc::extracted_cluster_info_from_python eci(cluster_info);

      GtkWidget *w = wrapped_create_cfc_dialog(eci);

      if (w) {
	 gtk_widget_show(w);
      } else {
	 std::cout << "Null w in chemical_feature_clusters_accept_info_py()" << std::endl;
      } 

   }
}


// should this be in a cfc widgets-specific file?
// 
GtkWidget *cfc::wrapped_create_cfc_dialog(const cfc::extracted_cluster_info_from_python &extracted_cluster_info) {

   GtkWidget *cfc_dialog = create_cfc_dialog();

   // we want a vector of (water) clusters.  For each cluster we need
   // to know what structures have waters in that cluster
   //
   std::map<int, std::vector<int> > cluster_map;
   std::map<int, std::vector<int> >::const_iterator it;
   for (unsigned int i=0; i<extracted_cluster_info.cw.size(); i++) {
      const clustered_water_info_from_python &cw = extracted_cluster_info.cw[i];

      // add imol if it is not already a member of this cluster
      if (std::find(cluster_map[cw.cluster_number].begin(),
		    cluster_map[cw.cluster_number].end(),
		    cw.imol) == cluster_map[cw.cluster_number].end())
	 cluster_map[cw.cluster_number].push_back(cw.imol);
   }
   // how many structures are there?
   unsigned int n_structures = extracted_cluster_info.n_structures();

   unsigned int cluster_idx = 0;
   double inv_n = 1.0/double(n_structures);

   // we need to convert the map to a vector, because we can sort a
   // vector using the number of structures, so that when we make the
   // buttons, the cluster with the most residues will appear at the
   // top
   // 
   std::vector<std::pair<std::vector<int>, water_cluster_info_from_python> > cluster_vec(extracted_cluster_info.cluster_idx_max()+1);
   for (it=cluster_map.begin(); it!=cluster_map.end(); it++) {

      unsigned int idx = it->first;
      if (idx < extracted_cluster_info.wc.size()) {
	 std::pair<std::vector<int>, water_cluster_info_from_python> p(it->second, extracted_cluster_info.wc[idx]);
	 cluster_vec[it->first] = p;
      } else {
	 std::cout << "ERROR::::::::::::::: indexing fail " << idx << " " << extracted_cluster_info.wc.size() << std::endl;
      }
   }
   std::sort(cluster_vec.begin(), cluster_vec.end(), extracted_cluster_info_from_python::cluster_vector_sorter);

   GtkWidget *waters_table = lookup_widget(cfc_dialog, "cfc_waters_table");

   if (! waters_table) {
      std::cout << "no waters_table" << std::endl;
   } else { 

      gtk_table_resize(GTK_TABLE(waters_table), cluster_vec.size(), 2);
   
      for (unsigned int i=0; i<cluster_vec.size(); i++) {
	 unsigned int n_this = cluster_vec[i].first.size();
	 double f = inv_n * n_this;
	 if (false)
	    std::cout << "   Water cluster " << i << " is present in " << n_this << " = "
		      << f*100 << " % of structures" << std::endl;

	 // button time!
	 //
	 // the left-hand water cluster button
	 // on clicked: centre on this cluster
	 //
	 std::string lb_label = "Water ";
	 lb_label += coot::util::int_to_string(i);
	 lb_label += ": ";
	 lb_label += coot::util::float_to_string_using_dec_pl(f*100, 1);
	 lb_label += " % conserved";
	 GtkWidget *left_button = gtk_button_new_with_label(lb_label.c_str());
	 // so that we can look up the widget names
	 unsigned int site_idx=0;
	 std::string button_name = "cfc_site_";
	 button_name += coot::util::int_to_string(site_idx);
	 button_name += "_water_cluster_";
	 button_name += coot::util::int_to_string(i);
	 button_name += "_button";
	 gtk_object_set_data_full (GTK_OBJECT (cfc_dialog), button_name.c_str(), 
				   left_button,
				   (GtkDestroyNotify) gtk_widget_unref);

	 water_cluster_info_from_python *wat_clust_p =
	    new water_cluster_info_from_python(cluster_vec[i].second);

	 std::cout << "creating a button with a index " << i << " and cluster centre " << cluster_vec[i].second.pos.format() << std::endl;

	 gtk_signal_connect(GTK_OBJECT(left_button), "clicked",
			    GTK_SIGNAL_FUNC(on_cfc_water_cluster_button_clicked),
			    (gpointer) wat_clust_p);
	 
	 gtk_table_attach(GTK_TABLE(waters_table), left_button,
			  0, 1, i, i+1,
			  (GtkAttachOptions) (GTK_FILL),
			  (GtkAttachOptions) (0), 0, 0);
	 gtk_widget_show(left_button);

	 // right hand buttons - one for every structure
	 // name: cfc_site_xx_water_cluster_yy_hbox
	 //
	 // we need an hbox to put them in
	 //
	 GtkWidget *hbox = gtk_hbox_new(FALSE, 0);

	 for (unsigned int j=0; j<n_structures; j++) {
	    std::string struct_button_name = "cfc_site_";
	    struct_button_name += coot::util::int_to_string(site_idx);
	    struct_button_name += "_water_cluster_";
	    struct_button_name += coot::util::int_to_string(i);
	    struct_button_name += "_structure_";
	    struct_button_name += coot::util::int_to_string(j);
	    struct_button_name += "_button";
	    std::string struct_button_label = " ";
	    GtkWidget *button = gtk_button_new_with_label(struct_button_label.c_str());
	    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	    gtk_widget_show(button);


	    // change the color to green when the structure (j) is a
	    // member of the vector cluster_vec[i];

	    if (std::find(cluster_vec[i].first.begin(), cluster_vec[i].first.end(), j) == cluster_vec[i].first.end()) {
	       // wasn't there (that's OK)
	    } else { 
	       // happy path
	       GdkColor color;
	       gdk_color_parse ("lightgreen", &color);
	       gtk_widget_modify_bg(GTK_WIDGET(button), GTK_STATE_NORMAL, &color);
	    }
	    
	 }

	 gtk_table_attach(GTK_TABLE(waters_table), hbox,
			  1,2,i,i+1,
			  (GtkAttachOptions) (GTK_FILL),
			  (GtkAttachOptions) (0), 0, 0);
	 gtk_widget_show(hbox);
	 
      }
   }
   
   return cfc_dialog;
}

void
cfc::on_cfc_water_cluster_button_clicked(GtkButton *button,
					 gpointer user_data) {

   water_cluster_info_from_python *water_clust_p = static_cast<water_cluster_info_from_python *> (user_data);

   std::cout << "go here " << water_clust_p->pos.format() << std::endl;
   coot::Cartesian c(water_clust_p->pos[0], water_clust_p->pos[1],    water_clust_p->pos[2]);
   graphics_info_t g;
   g.setRotationCentre(c);
   g.graphics_draw();
}



cfc::extracted_cluster_info_from_python::extracted_cluster_info_from_python(PyObject *cluster_info_py) {

   std::vector<water_cluster_info_from_python> v;
   std::vector<clustered_water_info_from_python> v_cw;

   if (! PyList_Check(cluster_info_py)) {
      std::cout << "ERROR:: not a list in cfc_extract_cluster_info()" << std::endl;
   } else {
      int list_size = PyObject_Length(cluster_info_py);
      if (list_size > 0) {
	 PyObject *water_cluster_info_py = PyList_GetItem(cluster_info_py, 0);

	 if (! PyList_Check(water_cluster_info_py)) {
	    std::cout << "ERROR:: not a list for water_cluster_info_py in cfc_extract_cluster_info()"
		      << std::endl;
	 } else {

	    int n_clusters = PyObject_Length(water_cluster_info_py);

	    std::cout << "found a water cluster list of length " << n_clusters << std::endl;

	    for (int iclust=0; iclust<n_clusters; iclust++) {
	       PyObject *cluster_py = PyList_GetItem(water_cluster_info_py, iclust);

	       // should be a list of length 3: with position, weight and cluster-sphere size/length/radius

	       if (! PyTuple_Check(cluster_py)) {

		  PyObject *dp = display_python(cluster_py);
		  if (dp == NULL) {
		     std::cout << "ERROR:: not a list for water_cluster item in "
			       << "cfc_extract_cluster_info() (null dp)" << std::endl;
		  } else { 
		     std::cout << "ERROR:: not a list for water_cluster item in cfc_extract_cluster_info()"
			       << PyString_AsString(dp) << std::endl;
		  }

	       } else {
		  
		  int n_items = PyObject_Length(cluster_py);

		  if (n_items != 3) {
		     std::cout << "strange cluster info " << n_items << std::endl;
		  } else {
		     
		     PyObject *pos_py    = PyTuple_GetItem(cluster_py, 0);
		     PyObject *weight_py = PyTuple_GetItem(cluster_py, 1);
		     PyObject *radius_py = PyTuple_GetItem(cluster_py, 2);

		     if (! PyList_Check(pos_py)) {
			std::cout << "ERROR:: position is not a list " << std::endl;
		     } else {
			int n_xyz = PyObject_Length(pos_py);
			if (n_xyz != 3) {
			   std::cout << "strange position list " << n_xyz << std::endl;
			} else {
			   PyObject *x_py = PyList_GetItem(pos_py, 0);
			   PyObject *y_py = PyList_GetItem(pos_py, 1);
			   PyObject *z_py = PyList_GetItem(pos_py, 2);

			   if (PyFloat_Check(x_py)) { // fatigue...
			      if (PyFloat_Check(y_py)) {
				 if (PyFloat_Check(z_py)) {

				    double x = PyFloat_AsDouble(x_py);
				    double y = PyFloat_AsDouble(y_py);
				    double z = PyFloat_AsDouble(z_py);

				    if (PyFloat_Check(weight_py)) {
				       if (PyFloat_Check(radius_py)) {

					  double weight = PyFloat_AsDouble(weight_py);
					  double radius = PyFloat_AsDouble(radius_py);

					  clipper::Coord_orth pos(x,y,z);
					  v.push_back(cfc::water_cluster_info_from_python(pos, weight, radius));
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }

      if (list_size > 1) { // it should be
	 PyObject *cluster_assignments_py = PyList_GetItem(cluster_info_py, 1);

	 // contains a list of [imol, water-residue-spec, cluster-number]

	 if (! PyList_Check(cluster_assignments_py)) {
	    std::cout << "ERROR:: cluster_assignments_py is not a list " << std::endl;
	 } else {
	    int n = PyObject_Length(cluster_assignments_py);
	    std::cout << "found " << n << " cluster assignments" << std::endl;
	    for (int i=0; i<n; i++) {
	       PyObject *item_py = PyList_GetItem(cluster_assignments_py, i);
	       if (! PyList_Check(item_py)) {
		  std::cout << "ERROR:: item in cluster_assignments_py is not a list " << std::endl;
	       } else {
		  int n_in_item = PyObject_Length(item_py);
		  if (n_in_item != 3) {
		     std::cout << "ERROR:: item in cluster_assignments_py is not of length 3 " << std::endl;
		  } else {
		     PyObject *imol_py   = PyList_GetItem(item_py, 0);
		     PyObject *spec_py   = PyList_GetItem(item_py, 1);
		     PyObject *iclust_py = PyList_GetItem(item_py, 2);

		     coot::residue_spec_t water_spec = residue_spec_from_py(spec_py);
		     int imol_water = PyInt_AsLong(imol_py);
		     int iclust     = PyInt_AsLong(iclust_py);

		     clustered_water_info_from_python cw(imol_water, water_spec, iclust);
		     v_cw.push_back(cw);
		  }
	       }
	    }
	 }
      }
   }

   std::cout << "------------ v.size() " << v.size() << std::endl;

   std::cout << "-------------v_cw.size() " << v_cw.size() << std::endl;

   wc = v;
   cw = v_cw;
}

unsigned int
cfc::extracted_cluster_info_from_python::n_structures() const {

   std::map<int, int> imol_map;

   for (unsigned int i=0; i<cw.size(); i++) {
      // is zero the default value?
      imol_map[cw[i].imol]++;
   }
   return imol_map.size();
}


unsigned int
cfc::extracted_cluster_info_from_python::cluster_idx_max() const {

   int idx_max = 0;
   for (unsigned int i=0; i<cw.size(); i++) {
      if (cw[i].cluster_number > idx_max)
	 idx_max = cw[i].cluster_number;
   }
   return idx_max;
}


#endif // USE_PYTHON
