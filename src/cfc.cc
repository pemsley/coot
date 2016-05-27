
// header-here

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#ifdef USE_PYTHON

// needed to parse cc-interface.hh
#ifdef USE_GUILE
#include <libguile.h> 
#endif // USE_GUILE

// For drawing balls
#if __APPLE__
#   include <OpenGL/gl.h>
#   include <OpenGL/glu.h>
#else
#   include <GL/gl.h>
#   include <GL/glu.h>
#endif

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
				       double radius_1, double radius_2) {

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
			   coot::residues_near_residue(ligand_spec, mol, radius_1);
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

			coot::chem_feat_solvated_ligand_spec lig(ligand_spec, neighbs_waters, mol, imol);
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

	 coot::chem_feat_clust cl(neighbs, ligands, radius_2, graphics_info_t::Geom_p());

	 std::vector<coot::chem_feat_clust::water_attribs> water_positions =
	    cl.get_water_positions();

	 std::vector<coot::simple_chemical_feature_attributes> chemical_features =
	    cl.get_chemical_features();

	 std::cout << "INFO:: Got " << water_positions.size() << " waters" << std::endl;
	 std::cout << "INFO:: Got " << chemical_features.size() << " chemical features"
		   << std::endl;

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

	 PyObject *chemical_feature_attribs_py = PyList_New(chemical_features.size());
	 for (unsigned int i=0; i<chemical_features.size(); i++) {
	    // need: type, position, imol, res_spec
	    PyObject *o = PyList_New(4);
	    PyObject *pos_py = PyList_New(3);
	    PyList_SetItem(pos_py, 0, PyFloat_FromDouble(chemical_features[i].pos.x()));
	    PyList_SetItem(pos_py, 1, PyFloat_FromDouble(chemical_features[i].pos.y()));
	    PyList_SetItem(pos_py, 2, PyFloat_FromDouble(chemical_features[i].pos.z()));

	    PyList_SetItem(o, 0, PyString_FromString(chemical_features[i].type.c_str()));
	    PyList_SetItem(o, 1, pos_py);
	    PyList_SetItem(o, 2, PyInt_FromLong(chemical_features[i].imol));
	    PyList_SetItem(o, 3, py_residue(chemical_features[i].residue_spec));
	    PyList_SetItem(chemical_feature_attribs_py, i, o);
	 }

	 // ------------------------------------------------------------------------
	 //            side chains 
	 // ------------------------------------------------------------------------

	 PyObject *residue_sidechain_attribs_py = PyList_New(1);
	 PyObject *x_py = PyInt_FromLong(12);

	 // PyList_SetItem(residue_sidechain_attribs_py, 0, x_py);

	 // ------------------------------------------------------------------------
	 //            combine
	 // ------------------------------------------------------------------------

	 r = PyList_New(2); // 3 with residues
	 
	 PyList_SetItem(r, 0, water_attribs_py);
	 PyList_SetItem(r, 1, chemical_feature_attribs_py);

      } else {
	 std::cout << "ERROR:: chemical_feature_clusters_py() arg 2 is not a list"
		   << std::endl;
      }
   } else {
      std::cout << "ERROR:: chemical_feature_clusters_py() arg 1 is not a list"
		<< std::endl;
   }

   
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

      int idx_gdo = eci.show_water_balls();

      if (w) {
	 gtk_widget_show(w);
      } else {
	 std::cout << "Null w in chemical_feature_clusters_accept_info_py()" << std::endl;
      } 

   }
}


// should this be in a cfc widgets-specific file?
// 
GtkWidget *cfc::wrapped_create_cfc_dialog(cfc::extracted_cluster_info_from_python &extracted_cluster_info) {

   GtkWidget *cfc_dialog = create_cfc_dialog();

   wrapped_create_cfc_dialog_add_waters(extracted_cluster_info, cfc_dialog);
   wrapped_create_cfc_dialog_add_pharmacophores(extracted_cluster_info, cfc_dialog);

   return cfc_dialog;
}

void cfc::wrapped_create_cfc_dialog_add_waters(cfc::extracted_cluster_info_from_python &extracted_cluster_info, GtkWidget *cfc_dialog) {

   // we want a vector of (water) clusters.  For each cluster we need
   // to know what structures have waters in that cluster
   //
   std::map<int, std::vector<int> > cluster_map;
   std::map<int, std::vector<int> >::const_iterator it;
   for (unsigned int i=0; i<extracted_cluster_info.cw.size(); i++) {
      const clustered_feature_info_from_python &cw = extracted_cluster_info.cw[i];

      // add imol if it is not already a member of this cluster
      if (std::find(cluster_map[cw.cluster_number].begin(),
		    cluster_map[cw.cluster_number].end(),
		    cw.imol) == cluster_map[cw.cluster_number].end())
	 cluster_map[cw.cluster_number].push_back(cw.imol);
   }
   // how many structures are there?
   unsigned int n_structures = extracted_cluster_info.n_water_structures();
   // if we iterate through n_structures, j=0; j<n_structures, to get
   // to the molecule number (imol) of the jth structure, index into
   // the structures_vec: imol = structures_vec[j];
   //
   std::vector<int> structures_vec = extracted_cluster_info.water_structures_vec();

   unsigned int cluster_idx = 0;
   double inv_n = 1.0/double(n_structures);

   // we need to convert the map to a vector, because we can sort a
   // vector using the number of structures, so that when we make the
   // buttons, the cluster with the most residues will appear at the
   // top
   // 
   std::vector<std::pair<std::vector<int>, water_cluster_info_from_python> > cluster_vec(extracted_cluster_info.water_cluster_idx_max()+1);
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

	 // we should do a better job at clearing up the memory when
	 // these buttons are destroyed
	 // 
	 gtk_object_set_data_full (GTK_OBJECT (cfc_dialog), button_name.c_str(), 
				   left_button,
				   (GtkDestroyNotify) gtk_widget_unref);

	 water_cluster_info_from_python *wat_clust_p =
	    new water_cluster_info_from_python(cluster_vec[i].second);

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

	    if (j < structures_vec.size()) {

	       int imol_local = structures_vec[j]; // should be some function of j...
	       // the question is: what is
	       // the imol (molecule number) of the jth
	       // structure?
	    
	       std::string struct_button_name = "cfc_site_";
	       struct_button_name += coot::util::int_to_string(site_idx);
	       struct_button_name += "_water_cluster_";
	       struct_button_name += coot::util::int_to_string(i);
	       struct_button_name += "_structure_";
	       struct_button_name += coot::util::int_to_string(imol_local);
	       struct_button_name += "_button";
	       std::string struct_button_label = " ";
	       GtkWidget *button = gtk_button_new_with_label(struct_button_label.c_str());
	       gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	       gtk_signal_connect(GTK_OBJECT(button), "clicked",
				  GTK_SIGNAL_FUNC(on_cfc_water_cluster_structure_button_clicked),
				  GINT_TO_POINTER(imol_local));
	    
	       // change the color to green when the structure (j) is a
	       // member of the vector cluster_vec[i];

	       if (std::find(cluster_vec[i].first.begin(), cluster_vec[i].first.end(), imol_local) == cluster_vec[i].first.end()) {
		  // wasn't there (that's OK)
	       } else {
		  // happy path
		  GdkColor color;
		  gdk_color_parse ("#AACCAA", &color);
		  gtk_widget_modify_bg(GTK_WIDGET(button), GTK_STATE_NORMAL, &color);
	       }
	       gtk_widget_show(button);

	    } else {
	       std::cout << "ERRROR:: out of index range jth mol " << j
			 << " " << structures_vec.size() << std::endl;
	    }
	    
	 }

	 gtk_table_attach(GTK_TABLE(waters_table), hbox,
			  1,2,i,i+1,
			  (GtkAttachOptions) (GTK_FILL),
			  (GtkAttachOptions) (0), 0, 0);
	 gtk_widget_show(hbox);
	 
      }
   }
}


void
cfc::wrapped_create_cfc_dialog_add_pharmacophores(cfc::extracted_cluster_info_from_python &extracted_cluster_info, GtkWidget *cfc_dialog) {

   GtkWidget *ligands_table = lookup_widget(cfc_dialog, "cfc_ligands_table");

   // for each cluster, by which structures are they contributed? And what is the
   // position of the cluster (needed for callback info construction).
   //
   std::map<std::string, std::vector<std::pair<std::vector<int>, clipper::Coord_orth> > > cluster_structure_vector;
   std::map<std::string, std::vector<clipper::Coord_orth> >::const_iterator it;

   for (it  = extracted_cluster_info.pharmacophore_model_cluster_means.begin();
	it != extracted_cluster_info.pharmacophore_model_cluster_means.end();
	it++) {

      const std::string &type = it->first; // e.g. Donor
      const std::vector<clipper::Coord_orth> &means = it->second;

      cluster_structure_vector[type].resize(means.size());

      unsigned int n_means = means.size();
      for (unsigned int i=0; i<n_means; i++) { // types match cluster_number
	 for (unsigned int j=0; j<extracted_cluster_info.pharmacophore[type].size(); j++) {

	    if (extracted_cluster_info.pharmacophore[type][j].cluster_number == i) {
	       if (std::find(cluster_structure_vector[type][i].first.begin(),
			     cluster_structure_vector[type][i].first.end(),
			     extracted_cluster_info.pharmacophore[type][j].imol) ==
		   cluster_structure_vector[type][i].first.end()) {
		  cluster_structure_vector[type][i].first.push_back(extracted_cluster_info.pharmacophore[type][j].imol);
		  cluster_structure_vector[type][i].second = means[i];
	       }
	    }
	 }
      }
   }

   // what did we make?
   //
   std::map<std::string, std::vector<std::pair<std::vector<int>, clipper::Coord_orth> > >::const_iterator it_2;
   unsigned int n_pharacophores = 0;

   for (it_2  = cluster_structure_vector.begin();
	it_2 != cluster_structure_vector.end();
	it_2++) {
      for (unsigned int i=0; i<it_2->second.size(); i++)
	 n_pharacophores++;
   }


   if (false) { 
      for (it_2  = cluster_structure_vector.begin();
	   it_2 != cluster_structure_vector.end();
	   it_2++) {
	 for (unsigned int i=0; i<it_2->second.size(); i++) {
	    if (false) { // debug
	       std::cout << "DEBUG:: cluster_structure_vector[" << it_2->first
			 << "][" << i << "] : ";
	       for (unsigned int j=0; j<it_2->second[i].first.size(); j++)
		  std::cout << " " << it_2->second[i].first[j];
	       std::cout << std::endl;
	    }
	 }
      }
      std::cout << "DEUBG:: found " << n_pharacophores << " pharmacophores " << std::endl;
   }


   // this table has a header, so indexing is +1 vertical
   gtk_table_resize(GTK_TABLE(ligands_table), n_pharacophores+1, 2);

   // OK, now we can make some buttons
   //
   // do we want to "flatten out" cluster_structure_vector so that the
   // pharmacophores appear in order of conservedness?  No need - we
   // can make a simple container class and sort at the end.
   //
   unsigned int n_structures = extracted_cluster_info.n_pharmacophore_structures();
   double inv_n = 1.0/double(n_structures);

   // we want to sort the pharmacophore buttons.
   // store them here and sort them when filled.
   //
   std::vector<pharm_button_set> pharm_button_set_collection; // fill then then sort it.
   
   for (it_2  = cluster_structure_vector.begin();
	it_2 != cluster_structure_vector.end();
	it_2++) {

      const std::string &type = it_2->first;

      for (unsigned int i=0; i<it_2->second.size(); i++) {

	 // imol_residue_spec_vec should contain the imol,res_spec of ligands
	 // that contribute to this pharmacophore, they are not in structure order and
	 // there can be any number of them.
	 //
	 std::vector<std::pair<int, coot::residue_spec_t> > imol_residue_spec_vec =
	    extracted_cluster_info.pharmacophore_cluster_imol_residue_spec_vec(type, i);

	 std::string lhb_label = it_2->first;
	 unsigned int n_this = it_2->second[i].first.size();
	 double f = inv_n * n_this;
	 lhb_label += " ";
	 lhb_label += coot::util::int_to_string(i);
	 lhb_label += ": ";
	 lhb_label += coot::util::float_to_string_using_dec_pl(f*100, 1);
	 lhb_label += " % conserved";
	 GtkWidget *left_button = gtk_button_new_with_label(lhb_label.c_str());

	 pharm_button_set pbs(it_2->first, left_button, f);
	 
	 // on clicked, go a cluster mean, show all structures that
	 // contribute to this cluster, and the chemical features of
	 // those ligands that match the type of the cluster.
	 // 
	 // We need cluster mean, and the contributors to this cluster.

	 GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
	 pbs.add_structure_buttons_hbox(hbox);
	 std::vector<std::pair<int, coot::residue_spec_t> > contributing_specs;

	 for (unsigned int j=0; j<n_structures; j++) {

	    if (true) {

	       // we should store this vector, not keep evaluating it
	       //
	       // pharmacophore_structures_vec()[j] yeilds an imol

	       int imol_this_structure = extracted_cluster_info.pharmacophore_structures_vec()[j];
	       std::pair<int, coot::residue_spec_t> imol_res_spec_this_structure =
		  extracted_cluster_info.pharmacophore_structures_and_specs_vec()[j];

	       // for contributing structures, used in the structure button callback.
	       //
	       // std::pair<int, coot::residue_spec_t> imol_res_spec_local =
	       // imol_residue_spec_vec[j];

	       // what is imol for the jth structure?
	       // Not imol_res_spec_local.first

	       // I don't think that it does - delete
	       // std::cout << "structure " << j << " has imol " << imol_res_spec_local.first
	       // << std::endl;
	       
	       std::string struct_button_label = " ";
	       GtkWidget *button = gtk_button_new_with_label(struct_button_label.c_str());

	       // on clicked, go to this pharmacophore position,
	       // show only this structure and this ligand and the pharmacophores
	       // for this ligand. We need imol, res_spec, and pharmacophore pos.
	       //
	       const clipper::Coord_orth &pos = it_2->second[i].second;
	       coot::residue_spec_t res_spec; // what is this now? Do we have access to the
	                                      // ligand spec of imol_this_structure?
	       on_pharmacophore_structure_click_info_t ci(pos, imol_this_structure, res_spec);

	       on_pharmacophore_structure_click_info_t *ci_p = new on_pharmacophore_structure_click_info_t(ci);
	       
	       gtk_signal_connect(GTK_OBJECT(button), "clicked",
				  GTK_SIGNAL_FUNC(on_cfc_pharmacophore_cluster_structure_button_clicked),
				  (gpointer)ci_p);

	       if (std::find(it_2->second[i].first.begin(),
			     it_2->second[i].first.end(),
			     imol_this_structure) == it_2->second[i].first.end()) {
		  
		  // wasn't there (that's OK)
		  
	       } else {
		  // happy path
		  GdkColor color;
		  gdk_color_parse("#AACCAA", &color);
		  gtk_widget_modify_bg(GTK_WIDGET(button), GTK_STATE_NORMAL, &color);

		  if (std::find(contributing_specs.begin(),
				contributing_specs.end(),
				imol_res_spec_this_structure) == contributing_specs.end()) {
		     contributing_specs.push_back(imol_res_spec_this_structure);
		  }
	       }
	       gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	       gtk_widget_show(button);
	    }
	 }
	 pharm_button_set_collection.push_back(pbs);
	 on_pharmacophore_click_info_t ci(it_2->second[i].second, contributing_specs);
	 // make a pointer to a copy of that that we can use in the callback user-data.
	 //
	 on_pharmacophore_click_info_t *ci_p = new on_pharmacophore_click_info_t(ci);
	 gtk_signal_connect(GTK_OBJECT(left_button), "clicked",
			    GTK_SIGNAL_FUNC(on_cfc_pharmacophore_cluster_button_clicked),
			    (gpointer) ci_p);

      }
   }

   std::sort(pharm_button_set_collection.begin(),
	     pharm_button_set_collection.end(),
	     pharm_button_set::sorter);

   int i_row = 1; // start after the header
   
   for (unsigned int i=0; i<pharm_button_set_collection.size(); i++) {
      const pharm_button_set &p = pharm_button_set_collection[i];
      gtk_table_attach(GTK_TABLE(ligands_table), p.left_hand_button,
		       0,1, i_row+i,i_row+i+1,
		       (GtkAttachOptions) (GTK_FILL),
		       (GtkAttachOptions) (0), 0, 0);
      gtk_table_attach(GTK_TABLE(ligands_table), p.structure_buttons_hbox,
		       1,2, i_row+i,i_row+i+1,
		       (GtkAttachOptions) (GTK_FILL),
		       (GtkAttachOptions) (0), 0, 0);
      gtk_widget_show(p.left_hand_button);
      gtk_widget_show(p.structure_buttons_hbox);
   }
}

void
cfc::on_cfc_water_cluster_button_clicked(GtkButton *button,
					 gpointer user_data) {

   water_cluster_info_from_python *water_clust_p =
      static_cast<water_cluster_info_from_python *> (user_data);

   coot::Cartesian c(water_clust_p->pos[0],
		     water_clust_p->pos[1],
		     water_clust_p->pos[2]);
   graphics_info_t g;
   g.setRotationCentre(c);
   g.display_all_model_molecules(); // no redraw
   g.graphics_draw();

}

void
cfc::on_cfc_water_cluster_structure_button_clicked(GtkButton *button,
						   gpointer user_data) {

   int imol = GPOINTER_TO_INT(user_data);

   // std::cout << "undisplay all except " << imol << std::endl;

   graphics_info_t g;
   g.undisplay_all_model_molecules_except(imol);  // no redraw
   g.graphics_draw();
   
}

// on clicked, go a cluster mean, show all structures that
// contribute to this cluster, and the chemical features of
// those ligands that match the type of the cluster.
// 
// We need cluster mean, and the contributors to this cluster.
//
void
cfc::on_cfc_pharmacophore_cluster_button_clicked(GtkButton *button,
						 gpointer user_data) {
   on_pharmacophore_click_info_t *ci_p =
      static_cast<on_pharmacophore_click_info_t *> (user_data);

   graphics_info_t g;
   coot::Cartesian c(ci_p->pos[0],
		     ci_p->pos[1],
		     ci_p->pos[2]);
   std::cout << "set rotation centre to " << ci_p->pos.format() << std::endl;
   g.setRotationCentre(c);
   
   std::vector<int> keep_these; // imols
   for (unsigned int i=0; i<ci_p->imol_residue_specs.size(); i++) {

      // imol if it is not already in the keep_these vector
      if (std::find(keep_these.begin(),
		    keep_these.end(), ci_p->imol_residue_specs[i].first) ==
	  keep_these.end())
	 keep_these.push_back(ci_p->imol_residue_specs[i].first);
   }

   if (false) { 
      std::cout << "undisplay all except ";
      for (unsigned int i=0; i<keep_these.size(); i++) { 
	 std::cout << " " << keep_these[i];
      }
      std::cout << std::endl;
   }
   
   g.undisplay_all_model_molecules_except(keep_these);  // no redraw
   g.graphics_draw();

}

void
cfc::on_cfc_pharmacophore_cluster_structure_button_clicked(GtkButton *button,
							   gpointer user_data) {

   on_pharmacophore_structure_click_info_t *ci_p =
      static_cast<on_pharmacophore_structure_click_info_t *> (user_data);

   graphics_info_t g;
   coot::Cartesian c(ci_p->pos[0],
		     ci_p->pos[1],
		     ci_p->pos[2]);
   g.setRotationCentre(c);
   // std::cout << "undisplay all except " << ci_p->imol << std::endl;

   g.undisplay_all_model_molecules_except(ci_p->imol);  // no redraw
   g.graphics_draw();

}


void
cfc::extracted_cluster_info_from_python::extract_water_info(PyObject *cluster_info_py) {

   std::vector<water_cluster_info_from_python> v;
   std::vector<clustered_feature_info_from_python> v_cw;

   int list_size = PyObject_Length(cluster_info_py);
   if (list_size > 0) {
      PyObject *water_cluster_info_py = PyList_GetItem(cluster_info_py, 0);

      if (! PyList_Check(water_cluster_info_py)) {
	 std::cout << "ERROR:: not a list for water_cluster_info_py in cfc_extract_cluster_info()"
		   << std::endl;
      } else {

	 int n_clusters = PyObject_Length(water_cluster_info_py);

	 for (int iclust=0; iclust<n_clusters; iclust++) {
	    PyObject *cluster_py = PyList_GetItem(water_cluster_info_py, iclust);

	    // should be a list of length 3: with position, weight and
	    // cluster-sphere size/length/radius

	    if (! PyTuple_Check(cluster_py)) {

	       PyObject *dp = display_python(cluster_py);
	       if (dp == NULL) {
		  std::cout << "ERROR:: not a list for water_cluster item in "
			    << "cfc_extract_cluster_info() (null dp)" << std::endl;
	       } else { 
		  std::cout << "ERROR:: not a list for water_cluster item in "
			    << "cfc_extract_cluster_info()" << PyString_AsString(dp)
			    << std::endl;
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

		  clustered_feature_info_from_python cw(imol_water, water_spec, iclust);
		  v_cw.push_back(cw);
	       }
	    }
	 }
      }
   }
   wc = v;
   cw = v_cw;

}

std::vector<clipper::Coord_orth>
cfc::extracted_cluster_info_from_python::extract_cluster_means(PyObject *means_py) {

   std::vector<clipper::Coord_orth> v;
   
   if (! PyList_Check(means_py)) {
      std::cout << "ERROR:: means_py not a list in extract_cluster_means()" << std::endl;
   } else {

      int n_means = PyObject_Length(means_py);
      for (int i=0; i<n_means; i++) { 
	 PyObject *mean_pos_py = PyList_GetItem(means_py, i);

	 if (! PyList_Check(mean_pos_py)) {
	    std::cout << "ERROR:: mean_pos_py not a list in extract_cluster_means()"
		      << std::endl;
	 } else {
	    int n_pos = PyObject_Length(mean_pos_py);
	    if (n_pos == 3) {
	       
	       PyObject *x_py = PyList_GetItem(mean_pos_py, 0);
	       PyObject *y_py = PyList_GetItem(mean_pos_py, 1);
	       PyObject *z_py = PyList_GetItem(mean_pos_py, 2);

	       if (PyFloat_Check(x_py)) {
		  if (PyFloat_Check(y_py)) {
		     if (PyFloat_Check(z_py)) {

			double x = PyFloat_AsDouble(x_py);
			double y = PyFloat_AsDouble(y_py);
			double z = PyFloat_AsDouble(z_py);
			clipper::Coord_orth pos(x,y,z);
			v.push_back(pos);
		     }
		  }
	       }
	    }
	 }
      }
   }

   return v;
}


void
cfc::extracted_cluster_info_from_python::extract_chemical_feature_info(PyObject *cf_py) {

   // cf_py is [type, features-annotated-by-cluster-number, cluster_means]

   // key is type
   // std::map<std::string, std::vector<clustered_feature_info_from_python> > pharmacophore;
   // std::map<std::string, std::vector<clipper::Coord_orth> > pharmacophore_model_cluster_means;

   if (! PyList_Check(cf_py)) {
      std::cout << "ERROR:: not a list 0 in extract_chemical_feature_info()"
		<< std::endl;
   } else {

      int n_cf = PyObject_Length(cf_py);
      if (n_cf == 3) {

	 PyObject *type_py  = PyList_GetItem(cf_py, 0);
	 PyObject *facn_py  = PyList_GetItem(cf_py, 1);
	 PyObject *means_py = PyList_GetItem(cf_py, 2);

	 int n = PyObject_Length(facn_py);

	 std::string type;
	 if (PyString_Check(type_py))
	    type = PyString_AsString(type_py);

	 pharmacophore_model_cluster_means[type] = extract_cluster_means(means_py);

	 for (int i=0; i<n; i++) {
	    PyObject *pharm_py = PyList_GetItem(facn_py, i);

	    // pharm_py [[[-21.515, 3.507, -8.572], 0, [True, 'A', 1501, '']], 4]
	    
	    if (! PyList_Check(pharm_py)) {
	       std::cout << "ERROR:: pharm_py - Not a list " << std::endl;
	    } else {

	       int n_py = PyObject_Length(pharm_py);
	       // std::cout << "pharm_py is a list of length " << n_py << std::endl;

	       if (n_py == 2) {
		  PyObject *pos_imol_spec_list_py = PyList_GetItem(pharm_py, 0);
		  PyObject *cluster_number_py     = PyList_GetItem(pharm_py, 1);

		  if (PyList_Check(pos_imol_spec_list_py)) { // fatigue
		     if (PyInt_Check(cluster_number_py)) {

			int cluster_number = PyInt_AsLong(cluster_number_py);

			int n_pos_imol_spec_list =
			   PyObject_Length(pos_imol_spec_list_py);

			if (n_pos_imol_spec_list == 3) {

			   PyObject *pos_py  = PyList_GetItem(pos_imol_spec_list_py, 0);
			   PyObject *imol_py = PyList_GetItem(pos_imol_spec_list_py, 1);
			   PyObject *spec_py = PyList_GetItem(pos_imol_spec_list_py, 2);

			   if (PyList_Check(pos_py)) {
			      if (PyInt_Check(imol_py)) {
				 if (PyList_Check(spec_py)) {

				    PyObject *x_py = PyList_GetItem(pos_py, 0);
				    PyObject *y_py = PyList_GetItem(pos_py, 1);
				    PyObject *z_py = PyList_GetItem(pos_py, 2);

				    int imol = PyInt_AsLong(imol_py);

				    coot::residue_spec_t res_spec =
				       residue_spec_from_py(spec_py);

				    // it's not a water of course
				    clustered_feature_info_from_python cwi(imol,
									   res_spec,
									   cluster_number);
				    // store cwi somewhere

				    pharmacophore[type].push_back(cwi);

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

cfc::extracted_cluster_info_from_python::extracted_cluster_info_from_python(PyObject *cluster_info_py) {

   if (! PyList_Check(cluster_info_py)) {
      std::cout << "ERROR:: not a list in cfc_extract_cluster_info()" << std::endl;
   } else {

      int n = PyObject_Length(cluster_info_py);

      if (n == 2) { // water then chemical_features

	 PyObject *o0 = PyList_GetItem(cluster_info_py, 0);
	 PyObject *o1 = PyList_GetItem(cluster_info_py, 1);
	 extract_water_info(o0);

	 if (PyList_Check(o1)) {
	    int n_fc = PyObject_Length(o1);
	    for (int i=0; i<n_fc; i++) {
	       PyObject *o_list_item_py = PyList_GetItem(o1, i);
	       extract_chemical_feature_info(o_list_item_py);
	    }
	 }
      }
   }
}

unsigned int
cfc::extracted_cluster_info_from_python::n_water_structures() const {

   std::map<int, int> imol_map;

   for (unsigned int i=0; i<cw.size(); i++) {
      // is zero the default value?
      imol_map[cw[i].imol]++;
   }
   return imol_map.size();

}

std::vector<int>
cfc::extracted_cluster_info_from_python::water_structures_vec() const {

   // I don't like using a list, I'd rather do it the usual way, with
   // a vector ans std::find().
   //
   std::list<int> imol_list;
   std::vector<int> imol_vec;

   for (unsigned int i=0; i<cw.size(); i++) {
      imol_list.push_back(cw[i].imol);
   }
   imol_list.sort();
   imol_list.unique();
   std::list<int>::const_iterator it;
   for (it=imol_list.begin(); it!=imol_list.end(); it++) {
      imol_vec.push_back(*it);
   }
   return imol_vec;
}


// how many structures have pharmacophores?
unsigned int
cfc::extracted_cluster_info_from_python::n_pharmacophore_structures() const {

   return pharmacophore_structures_vec().size();
}

// running through the indexing of pharmachore_structures_vec by
// i=0->n_pharacophores_structures() give us a set of imols that
// have pharmocophores.  Used for pharmacophe structure buttons.
// 
std::vector<int> 
cfc::extracted_cluster_info_from_python::pharmacophore_structures_vec() const {

   std::vector<int> imols_collection;

   std::map<std::string, std::vector<clustered_feature_info_from_python> >::const_iterator it;   
   for (it=pharmacophore.begin(); it!=pharmacophore.end(); it++) {
      for (unsigned int i=0; i<it->second.size(); i++) {

	 if (std::find(imols_collection.begin(),
		       imols_collection.end(),
		       it->second[i].imol) == imols_collection.end())
	    imols_collection.push_back(it->second[i].imol);
      }
   }

   std::sort(imols_collection.begin(), imols_collection.end());
   
   return imols_collection;
   
}

// as above, but we return the spec too
//
std::vector<std::pair<int, coot::residue_spec_t> >
cfc::extracted_cluster_info_from_python::pharmacophore_structures_and_specs_vec() const {

   std::vector<std::pair<int, coot::residue_spec_t> > imols_collection;

   std::map<std::string, std::vector<clustered_feature_info_from_python> >::const_iterator it;   

   for (it=pharmacophore.begin(); it!=pharmacophore.end(); it++) {
      for (unsigned int i=0; i<it->second.size(); i++) {

	 std::pair<int, coot::residue_spec_t> p(it->second[i].imol, it->second[i].residue_spec);
	 if (std::find(imols_collection.begin(),
		       imols_collection.end(),
		       p) == imols_collection.end()) {
	    imols_collection.push_back(p);
	 }
      }
   }

   std::sort(imols_collection.begin(), imols_collection.end());
   
   return imols_collection;
}


std::vector<std::pair<int, coot::residue_spec_t> >
cfc::extracted_cluster_info_from_python::water_cluster_imol_residue_spec_vec() const {

   // this needs checking (doing a == test on a pair in a vector)

   if (false) {
      std::cout << "in imol_residue_spec_vec() checking cw:" << std::endl;
      for (unsigned int i=0; i<cw.size(); i++) {
	 std::cout << "cw[" << i << "] " << cw[i].imol << " " << cw[i].residue_spec
		   << std::endl;
      }
   }
   
   std::vector<std::pair<int, coot::residue_spec_t> > v;
   for (unsigned int i=0; i<cw.size(); i++) {
      std::pair<int, coot::residue_spec_t> p(cw[i].imol, cw[i].residue_spec);
      if (std::find(v.begin(), v.end(), p) == v.end())
	 v.push_back(p);
   }
   return v;
}


// which imol,residue_specs are contributing to the pharmacophore[type][idx] ? 
// 
std::vector<std::pair<int, coot::residue_spec_t> >
cfc::extracted_cluster_info_from_python::pharmacophore_cluster_imol_residue_spec_vec(const std::string &type, unsigned int cluster_idx) {

   std::vector<std::pair<int, coot::residue_spec_t> > v;

   // what about if two ligand pharmacophores contribute to the same cluster (model)?

   for (unsigned int i=0; i<pharmacophore_model_cluster_means[type].size(); i++) {
      for (unsigned int j=0; j<pharmacophore[type].size(); j++) { 
	 if (pharmacophore[type][j].cluster_number == i) {
	    std::pair<int, coot::residue_spec_t> p(pharmacophore[type][j].imol,
						   pharmacophore[type][j].residue_spec);
	    v.push_back(p);
	 }
      }
   }
   return v;
}




unsigned int
cfc::extracted_cluster_info_from_python::water_cluster_idx_max() const {

   unsigned int idx_max = 0;
   for (unsigned int i=0; i<cw.size(); i++) {
      if (cw[i].cluster_number > idx_max)
	 idx_max = cw[i].cluster_number;
   }
   return idx_max;
}

// return the generic display object index
int
cfc::extracted_cluster_info_from_python::show_water_balls() const {

   coot::generic_display_object_t obj("CFC conserved water balls");
   
   if (graphics_info_t::use_graphics_interface_flag) {
      int n = n_water_structures();
      unsigned int n_wc = wc.size();
   
      for (unsigned int i=0; i<n_wc; i++) {
	 unsigned int n_this = 0;
	 // count how many have this cluster
	 for (unsigned int j=0; j<cw.size(); j++) { 
	    if (cw[j].cluster_number == i) {
	       n_this++;
	    }
	 }

	 double f = double(n_this)/double(n);

	 if (f > 0.01) {

	    double radius = f * 1.1;
	    coot::generic_display_object_t::sphere_t sphere(wc[i].pos, radius);
	    sphere.col = coot::colour_t(0.9, 0.2, 0.2);
	    obj.spheres.push_back(sphere);
	 }
      }
   }

   obj.is_displayed_flag = true;
   graphics_info_t::generic_objects_p->push_back(obj);
   graphics_info_t::graphics_draw();

   return (graphics_info_t::generic_objects_p->size() -1);
   
}


#endif // MAKE_ENHANCED_LIGAND_TOOLS


#endif // USE_PYTHON
