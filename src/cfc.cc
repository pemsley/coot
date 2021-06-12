/* src/cfc.cc
 * 
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
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


#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#ifdef USE_PYTHON

#include "Python.h" // _XOPEN_SOURCE definition - python one before /usr/include/features.h

#include <cstddef>

#include <epoxy/gl.h>

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
				       PyObject *solvated_ligand_info_py, // [imol, ligand-specs]s
				       double radius_1, double radius_2) {

   std::cout << "debug:: ################ chemical_feature_clusters_py() start" << std::endl;

   PyObject *r = Py_False;

   if (PyList_Check(environment_residues_py)) {
      if (PyList_Check(solvated_ligand_info_py)) {

	 std::vector<coot::residue_spec_t> neighbs =
	    py_to_residue_specs(environment_residues_py);

	 std::vector<coot::chem_feat_solvated_ligand_spec> ligands; // fill this from
         	                                           // solvated_ligand_info_py

	 int n = PyObject_Length(solvated_ligand_info_py);
         std::cout << "debug:: ################ chemical_feature_clusters_py() here with size of imol_ligand_specs "
                   << n << std::endl;
	 for(int i=0; i<n; i++) {

            // std::cout << "debug:: ################ chemical_feature_clusters_py() here with i " << i << std::endl;
	    PyObject *o = PyList_GetItem(solvated_ligand_info_py, i);

	    // o should be a list of molecule-idx, residue-spec
	    if (PyList_Check(o)) {
	       int no = PyObject_Length(o);
               // std::cout << "debug:: ################ chemical_feature_clusters_py() here with no " << no << std::endl;
	       if (no == 2) {
		  PyObject *mol_idx_py  = PyList_GetItem(o, 0);
		  PyObject *lig_spec_py = PyList_GetItem(o, 1);

		  if (PyLong_Check(mol_idx_py)) {
		     int imol = PyLong_AsLong(mol_idx_py);
		     coot::residue_spec_t ligand_spec = residue_spec_from_py(lig_spec_py);

		     if (graphics_info_t::is_valid_model_molecule(imol)) {

			mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;

			std::vector<coot::residue_spec_t> neighbs_waters; // fill this
			std::vector<coot::residue_spec_t> neighbs_raw =
			   coot::residues_near_residue(ligand_spec, mol, radius_1);
                        if (false)
                           std::cout << "debug:: ##### in chemical_feature_clusters_py() neighbs_raw size is "
                                     << neighbs_raw.size() << std::endl;

			for (unsigned int i_neighb=0; i_neighb<neighbs_raw.size(); i_neighb++) { 
			   mmdb::Residue *res = coot::util::get_residue(neighbs_raw[i_neighb], mol);
			   if (res) {
			      std::string res_name = res->GetResName();
			      if (res_name == "HOH")
				 neighbs_waters.push_back(neighbs_raw[i_neighb]);
			   } else {
                              std::cout << "debug:: ##### in chemical_feature_clusters_py() got a NULL residue"
                                        << std::endl;
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
	    PyList_SetItem(o, 0, PyLong_FromLong(water_positions[iw].ligand_idx));
	    PyList_SetItem(o, 1, residue_spec_to_py(water_positions[iw].residue_spec()));
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

	    PyList_SetItem(o, 0, PyUnicode_FromString(chemical_features[i].type.c_str()));
	    PyList_SetItem(o, 1, pos_py);
	    PyList_SetItem(o, 2, PyLong_FromLong(chemical_features[i].imol));
	    PyList_SetItem(o, 3, residue_spec_to_py(chemical_features[i].residue_spec));
	    PyList_SetItem(chemical_feature_attribs_py, i, o);
	 }

	 // ------------------------------------------------------------------------
	 //            side chains 
	 // ------------------------------------------------------------------------

	 PyObject *residue_sidechain_attribs_py = PyList_New(1);
	 PyObject *x_py = PyLong_FromLong(12);

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

void chemical_feature_clusters_accept_info_py(unsigned int site_number,
					      PyObject *env_residue,
					      PyObject *mol_ligand_specs,
					      PyObject *cluster_info) {

   std::cout << "debug:: ################################## chemical_feature_clusters_accept_info_py()"
             << std::endl;

   if (graphics_info_t::use_graphics_interface_flag) {
      cfc::extracted_cluster_info_from_python eci(cluster_info);
      std::cout << "::::::::::::::::::::::::::: in chemical_feature_clusters_accept_info_py() "
		<< site_number << " graphics_info_t::cfc_dialog test " << std::endl;
      if (graphics_info_t::cfc_dialog)
	 chemical_feature_clusters_add_site_info(site_number, eci, graphics_info_t::cfc_dialog);
      else
	 std::cout << "::::::::::::::::::::::::::: in chemical_feature_clusters_accept_info_py() "
		   << site_number << " graphics_info_t::cfc_dialog is null " << std::endl;
      int idx_gdo = eci.show_water_balls(site_number);
   }
}

void chemical_feature_clusters_setup_dialog() {

   GtkWidget *w = create_cfc_dialog();
   if (w) {
      graphics_info_t::cfc_dialog = w;
   } else {
      std::cout << "Null w in chemical_feature_clusters_accept_info_py()" << std::endl;
   }
}


// where are the ligand binding sites?
//
// site_info_py looks like this:
// 
// [(0, [0, [True, 'A', 1501, '']]),
//  (1, [0, [True, 'A', 1501, '']]),
//  (0, [0, [True, 'A', 1501, '']]),
//  (0, [0, [True, 'A', 1501, '']]),
//  (0, [0, [True, 'A', 1501, '']]),
//  (0, [0, [True, 'A', 1501, '']]),
//  (0, [0, [True, 'A', 1501, '']]),
//  (0, [0, [True, 'A', 1501, '']]),  .. ]
// for evey atom in every ligand:
// [(site_idx, [imol, residue_spec])]
// 
PyObject *chemical_feature_clusters_accept_site_clusters_info_py(PyObject *site_info_py) {

   // we need to make a list of ligand sites and the
   // (imol, residue_spec) pairs the contribute to each site.
   // 
   // We are passed a list of
   // [cluster_idx, [imol, residue_spec]] for every atom of every ligand, so
   // when constructing the vector of pairs, we will check that that pair
   // is not already there.
   //
   // the key is the cluster index
   //
   std::map<int, std::vector<std::pair<int, coot::residue_spec_t> > > ligand_sites;

   if (! PyList_Check(site_info_py)) {
      std::cout << "chemical_feature_clusters_accept_site_clusters_info_py "
		<< "site_info_py is not a list" << std::endl;
   } else {
      int l = PyObject_Length(site_info_py);
      std::cout << "chemical_feature_clusters_accept_site_clusters_info_py"
		<< " site_info_py length " << l << std::endl;
      for (int i=0; i<l; i++) {
	 PyObject *tup = PyList_GetItem(site_info_py, i);

	 if (PyTuple_Check(tup)) {
	    PyObject *site_idx_py          = PyTuple_GetItem(tup, 0);
	    PyObject *imol_residue_spec_py = PyTuple_GetItem(tup, 1);
	    
	    if (PyLong_Check(site_idx_py) || PyLong_Check(site_idx_py)) { // Argh! PyLong_Check() passes on Mac, fails on PC.
	                                                                 // (it was actually a numpy int64, but keep
	                                                                 //  this extra test now that we have it)
	       if (PyList_Check(imol_residue_spec_py)) {
		  int ll = PyObject_Length(imol_residue_spec_py);
		  if (ll == 2) {
		     int site_idx = -1;
		     if (PyLong_Check(site_idx_py))
			 site_idx = PyLong_AsLong(site_idx_py);
		     if (PyLong_Check(site_idx_py))
			 site_idx = PyLong_AsLong(site_idx_py);
		     PyObject *imol_py     = PyList_GetItem(imol_residue_spec_py, 0);
		     PyObject *res_spec_py = PyList_GetItem(imol_residue_spec_py, 1);
		     int imol = PyLong_AsLong(imol_py);
		     std::pair<bool, coot::residue_spec_t> spec =
			make_residue_spec_py(res_spec_py);

		     if (spec.first) { 
			if (false)
			   std::cout << "site: " << site_idx << " " << imol
				     << " " << spec.second << std::endl;
			
			std::pair<int, coot::residue_spec_t> p(imol, spec.second);
			if (std::find(ligand_sites[site_idx].begin(),
				      ligand_sites[site_idx].end(),
				      p) == ligand_sites[site_idx].end())
			   ligand_sites[site_idx].push_back(p);
		     }
		  }
	       }
	    } else {
	       std::cout << "site_idx_py was not a PyInt or a PyLong" << std::endl;
	       PyObject *o = PyObject_Type(site_idx_py);
	       PyObject *dp2 = display_python(o);
	       if (dp2 == NULL) {
		  std::cout << "ERROR:: chemical_feature_clusters_accept_site_clusters_info_py (null dp)" << std::endl;
	       } else { 
		  std::cout << "ERROR:: chemical_feature_clusters_accept_site_clusters_info_py() site_idx_py type: "
			    << PyUnicode_AsUTF8String(dp2) << std::endl;
	    }
	    } 
	 }
      }
   }

   // return a list of pairs of imol,residue_specs
   //
   std::cout << "---------------------------- debug creating ligand_sites_py with size " << ligand_sites.size() << std::endl;
   PyObject *ligand_sites_py = PyList_New(ligand_sites.size());
   
   std::map<int, std::vector<std::pair<int, coot::residue_spec_t> > >::const_iterator it;
   unsigned int list_idx = 0; // because we are using iterator
   for (it = ligand_sites.begin(); it != ligand_sites.end(); ++it) {
      PyObject *li = PyList_New(it->second.size());
      for (unsigned int i=0; i<it->second.size(); i++) { 
	 PyObject *l = PyList_New(2);
	 PyList_SetItem(l, 0, PyLong_FromLong(it->second[i].first));
	 PyList_SetItem(l, 1, residue_spec_to_py(it->second[i].second));
	 PyList_SetItem(li, i, l);
      }
      PyList_SetItem(ligand_sites_py, list_idx, li);
      list_idx++; // for next round
   }

   return ligand_sites_py;
}


void
cfc::cfc_dialog_add_waters(unsigned int site_number,
			   cfc::extracted_cluster_info_from_python extracted_cluster_info,
			   GtkWidget *cfc_dialog) {

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

   double inv_n = 1.0/double(n_structures);

   // we need to convert the map to a vector, because we can sort a
   // vector using the number of structures, so that when we make the
   // buttons, the cluster with the most residues will appear at the
   // top
   //

   std::cout << "--------- cfc_dialog_add_waters() for site " << site_number
	     << " with water_cluster_idx_max() "
	     << extracted_cluster_info.water_cluster_idx_max() << std::endl;
   
   std::vector<std::pair<std::vector<int>, water_cluster_info_from_python> > cluster_vec(extracted_cluster_info.water_cluster_idx_max()+1);
   for (it=cluster_map.begin(); it!=cluster_map.end(); ++it) {

      unsigned int idx = it->first;
      if (idx < extracted_cluster_info.wc.size()) {
	 std::pair<std::vector<int>, water_cluster_info_from_python> p(it->second, extracted_cluster_info.wc[idx]);
	 cluster_vec[it->first] = p;
      } else {
	 std::cout << "ERROR::::::::::::::: indexing fail " << idx << " " << extracted_cluster_info.wc.size()
		   << std::endl;
      }
   }
   std::sort(cluster_vec.begin(), cluster_vec.end(), extracted_cluster_info_from_python::cluster_vector_sorter);

   // Now waters_table is not made by glade, we make a new one of them for each site
   // and add the table to the vbox
   
   // GtkWidget *waters_table = lookup_widget(cfc_dialog, "cfc_waters_table");

   GtkWidget *waters_vbox = lookup_widget(cfc_dialog, "cfc_waters_vbox");
   GtkWidget *waters_table = gtk_table_new(cluster_vec.size(), 2, FALSE);

   std::string waters_table_name = "cfc_waters_table_site_";
   waters_table_name += coot::util::int_to_string(site_number);

   // we want to be able to look up the widget and undisplay it.
   // 
   // old
   // gtk_object_set_data_full (GTK_OBJECT (cfc_dialog), waters_table_name.c_str(),
   //                           waters_table,
   // 			        (GtkDestroyNotify) gtk_widget_unref);
   g_object_set_data(G_OBJECT(cfc_dialog), waters_table_name.c_str(), waters_table);

   // perhaps we want to do this actually?
   std::string *wtn_p = new std::string(waters_table_name);
   gpointer gp = (gpointer) wtn_p;
   g_object_set_data(G_OBJECT(waters_table), "widget_name", gp);
   
   gtk_box_pack_start(GTK_BOX(waters_vbox), waters_table, FALSE, FALSE, 0);
   if (site_number == 0)
      gtk_widget_show(waters_table);


   std::cout << ":::::::::::::::::: cfc_waters_table: " << waters_table << std::endl;

   if (! waters_table) {
      std::cout << "no waters_table" << std::endl;
   } else {

      // gtk_table_resize(GTK_TABLE(waters_table), cluster_vec.size(), 2);

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
	 std::string button_name = "cfc_site_";
	 button_name += coot::util::int_to_string(site_number);
	 button_name += "_water_cluster_";
	 button_name += coot::util::int_to_string(i);
	 button_name += "_button";

	 // we should do a better job at clearing up the memory when
	 // these buttons are destroyed
	 //
	 // std::cout << "debug:: gtk_object_set_data_full() on button with name " << button_name
	 // << " and label: " << lb_label << std::endl;
	 // gtk_object_set_data_full (GTK_OBJECT (cfc_dialog), button_name.c_str(),
         // left_button, (GtkDestroyNotify) gtk_widget_unref);
         //
         g_object_set_data(G_OBJECT(cfc_dialog), button_name.c_str(), left_button);

	 water_cluster_info_from_python *wat_clust_p =
	    new water_cluster_info_from_python(cluster_vec[i].second);

	 g_signal_connect(G_OBJECT(left_button), "clicked",
                          G_CALLBACK(on_cfc_water_cluster_button_clicked),
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
	       struct_button_name += coot::util::int_to_string(site_number);
	       struct_button_name += "_water_cluster_";
	       struct_button_name += coot::util::int_to_string(i);
	       struct_button_name += "_structure_";
	       struct_button_name += coot::util::int_to_string(imol_local);
	       struct_button_name += "_button";
	       std::string struct_button_label = " ";

	       GtkWidget *button = gtk_button_new_with_label(struct_button_label.c_str());
	       gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
	       g_signal_connect(G_OBJECT(button), "clicked",
                                G_CALLBACK(on_cfc_water_cluster_structure_button_clicked),
                                GINT_TO_POINTER(imol_local));
	    
	       // change the color to green when the structure (j) is a
	       // member of the vector cluster_vec[i];

	       if (std::find(cluster_vec[i].first.begin(), cluster_vec[i].first.end(), imol_local) == cluster_vec[i].first.end()) {
		  // wasn't there (that's OK)
	       } else {
		  // happy path
		  GdkColor color;
		  gdk_color_parse("#AACCAA", &color);
		  gtk_widget_modify_bg(GTK_WIDGET(button), GTK_STATE_NORMAL, &color);
                  gtk_widget_set_name(GTK_WIDGET(button), "cfc-green-button");
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
cfc::cfc_dialog_add_pharmacophores(unsigned int site_number,
				   cfc::extracted_cluster_info_from_python extracted_cluster_info,
				   GtkWidget *cfc_dialog) {

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
	++it_2) {
      for (unsigned int i=0; i<it_2->second.size(); i++)
	 n_pharacophores++;
   }


   if (true) { 
      for (it_2  = cluster_structure_vector.begin();
	   it_2 != cluster_structure_vector.end();
	   ++it_2) {
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
      std::cout << "DEBUG:: found " << n_pharacophores << " pharmacophores " << std::endl;
   }

   GtkWidget *cfc_ligands_vbox = lookup_widget(cfc_dialog, "cfc_ligands_vbox");
   GtkWidget *ligands_table = gtk_table_new(n_pharacophores+1, 2, FALSE);

   std::string ligands_table_name = "cfc_ligands_table_site_";
   ligands_table_name += coot::util::int_to_string(site_number);
   
   // we want to be able to look up the widget and undisplay it.
   // 
   // gtk_object_set_data_full (GTK_OBJECT (cfc_dialog), ligands_table_name.c_str(),
   //   ligands_table, (GtkDestroyNotify) gtk_widget_unref);
   g_object_set_data(G_OBJECT(cfc_dialog), ligands_table_name.c_str(), ligands_table);

   // perhaps we want to do this actually?
   std::string *ltn_p = new std::string(ligands_table_name);
   gpointer gp = (gpointer) ltn_p;
   g_object_set_data(G_OBJECT(ligands_table), "widget_name", gp);

   gtk_box_pack_start(GTK_BOX(cfc_ligands_vbox), ligands_table, FALSE, FALSE, 0);
   if (site_number == 0)
      gtk_widget_show(ligands_table);

   // this table has a header, so indexing is +1 vertical
   // gtk_table_resize(GTK_TABLE(ligands_table), n_pharacophores+1, 2);

   // OK, now we can make some buttons
   //
   // do we want to "flatten out" cluster_structure_vector so that the
   // pharmacophores appear in order of conservedness?  No need - we
   // can make a simple container class and sort at the end.
   //
   unsigned int n_structures = extracted_cluster_info.n_pharmacophore_structures();

   std::cout << "debug:: in cfc_dialog_add_pharmacophores() n_structures: " << n_structures << std::endl;

   if (n_structures == 0) {

      std::cout << "WARNING:: Ooops! in cfc_dialog_add_pharmacophores() n_structures is 0" << std::endl;

   } else {

      // Happy Path
      
      double inv_n = 1.0/double(n_structures);

      // we want to sort the pharmacophore buttons.
      // store them here and sort them when filled.
      //
      std::vector<pharm_button_set> pharm_button_set_collection; // fill then then sort it.
   
      for (it_2  = cluster_structure_vector.begin();
	   it_2 != cluster_structure_vector.end();
	   ++it_2) {

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
	       
		  g_signal_connect(G_OBJECT(button), "clicked",
                                   G_CALLBACK(on_cfc_pharmacophore_cluster_structure_button_clicked),
                                   (gpointer)ci_p);

		  if (std::find(it_2->second[i].first.begin(),
				it_2->second[i].first.end(),
				imol_this_structure) == it_2->second[i].first.end()) {
		  
		     // wasn't there (that's OK)
                     gtk_button_set_label(GTK_BUTTON(button), "_"); // space doesn't match the size of an I
		  
		  } else {
		     // happy path, let's change the button colour and add it to contributing specs
		     GdkColor color;
		     gdk_color_parse("#AACCAA", &color);

                     if (false) // we know now that this is working right - it's the window manager that
                                // obscures the colour of the background
                        std::cout << "debug:: cfc_dialog_add_pharmacophores() setting the button bg to #AACCAA " << std::endl;
		     gtk_widget_modify_bg(GTK_WIDGET(button), GTK_STATE_NORMAL, &color);
                     gtk_widget_set_name(GTK_WIDGET(button), "cfc-green-button");
                     gtk_button_set_label(GTK_BUTTON(button), "I");

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
	    g_signal_connect(G_OBJECT(left_button), "clicked",
                             G_CALLBACK(on_cfc_pharmacophore_cluster_button_clicked),
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
}

void
cfc::on_cfc_site_button_clicked(GtkButton *button,
				gpointer user_data) {

   // we also need to look up the vboxes of ligands, residues and
   // waters and display/undisplay the tables

   if (user_data) {

      std::pair<int, clipper::Coord_orth> *pos_p =
	 static_cast<std::pair<int, clipper::Coord_orth> *> (user_data);

      int site_number = pos_p->first;
      // undisplay other tables
      GtkWidget *ligands_vbox  = lookup_widget(GTK_WIDGET(button), "cfc_ligands_vbox");
      GtkWidget *waters_vbox   = lookup_widget(GTK_WIDGET(button), "cfc_waters_vbox");
      GtkWidget *residues_vbox = lookup_widget(GTK_WIDGET(button), "cfc_residues_vbox");

      if (ligands_vbox) {
	 std::string show_this_one_name = "cfc_ligands_table_site_";
	 show_this_one_name += coot::util::int_to_string(site_number);
	 cfc_table_show_hide(show_this_one_name, ligands_vbox);
	 
      }
      if (waters_vbox) {
	 std::string show_this_one_name = "cfc_waters_table_site_";
	 show_this_one_name += coot::util::int_to_string(site_number);
	 cfc_table_show_hide(show_this_one_name, waters_vbox);
      }
      if (residues_vbox) {
	 std::string show_this_one_name = "cfc_residues_table_site_";
	 show_this_one_name += coot::util::int_to_string(site_number);
	 cfc_table_show_hide(show_this_one_name, residues_vbox);
      }
      
      coot::Cartesian c(pos_p->second.x(), pos_p->second.y(), pos_p->second.z());
      graphics_info_t g;
      g.setRotationCentre(c);
      g.graphics_draw();
   }
}


void
cfc::cfc_table_show_hide(std::string show_this_one_name, GtkWidget *vbox) {
   
   GList *dlist = gtk_container_get_children(GTK_CONTAINER(vbox));
   GList *free_list = dlist;
   
   while (dlist) {
      GtkWidget *list_item = (GtkWidget *) (dlist->data);
      // get and test the name
      gpointer gp = g_object_get_data(G_OBJECT(list_item), "widget_name");
      if (gp) {
	 std::string *name = static_cast<std::string *> (gp);
	 if (*name == show_this_one_name) {
	    gtk_widget_show(GTK_WIDGET(list_item));
	 } else {
	    gtk_widget_hide(GTK_WIDGET(list_item));
	 }
      }
      dlist = dlist->next;
   }
   g_list_free(free_list);
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

std::pair<bool, clipper::Coord_orth>
cfc::extracted_cluster_info_from_python::pharmacophores_centre() const {

   std::pair<bool, clipper::Coord_orth> pos_pair(false, clipper::Coord_orth(0,0,0));
   clipper::Coord_orth sum(0,0,0);
   unsigned int n = 0;

   std::map<std::string, std::vector<clipper::Coord_orth> >::const_iterator it;

   for (it =  pharmacophore_model_cluster_means.begin();
	it != pharmacophore_model_cluster_means.end();
	it++) {
      for (unsigned int i=0; i<it->second.size(); i++) { 
	 sum += it->second[i];
	 n++;
      }
   }

   if (n > 0) {
      pos_pair.first = true;
      double f(n);
      pos_pair.second = clipper::Coord_orth(sum.x()/f, sum.y()/f, sum.z()/f);
   }

   return pos_pair;
}



void
cfc::extracted_cluster_info_from_python::extract_water_info(PyObject *cluster_info_py) {

   std::vector<water_cluster_info_from_python> v;
   std::vector<clustered_feature_info_from_python> v_cw;

   int list_size = PyObject_Length(cluster_info_py);

   std::cout << "extracted_cluster_info_from_python::extract_water_info() list_size is " << list_size
             << std::endl;

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
			    << "cfc_extract_cluster_info()" << PyUnicode_AsUTF8String(dp)
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
		  int imol_water = PyLong_AsLong(imol_py);
		  int iclust     = PyLong_AsLong(iclust_py);

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
	 if (PyUnicode_Check(type_py))
	    type = PyBytes_AS_STRING(PyUnicode_AsUTF8String(type_py));

	 pharmacophore_model_cluster_means[type] = extract_cluster_means(means_py);

	 for (int i=0; i<n; i++) {
	    PyObject *pharm_py = PyList_GetItem(facn_py, i);

	    // pharm_py [[[-21.515, 3.507, -8.572], 0, [True, 'A', 1501, '']], 4]
	    
	    if (! PyList_Check(pharm_py)) {
	       std::cout << "ERROR:: pharm_py - Not a list " << std::endl;
	    } else {

	       int n_py = PyObject_Length(pharm_py);

	       if (n_py == 2) {
		  PyObject *pos_imol_spec_list_py = PyList_GetItem(pharm_py, 0);
		  PyObject *cluster_number_py     = PyList_GetItem(pharm_py, 1);

		  if (PyList_Check(pos_imol_spec_list_py)) {
		     if (PyLong_Check(cluster_number_py)) {

			int cluster_number = PyLong_AsLong(cluster_number_py);

			int n_pos_imol_spec_list =
			   PyObject_Length(pos_imol_spec_list_py);

			if (n_pos_imol_spec_list == 3) {

			   PyObject *pos_py  = PyList_GetItem(pos_imol_spec_list_py, 0);
			   PyObject *imol_py = PyList_GetItem(pos_imol_spec_list_py, 1);
			   PyObject *spec_py = PyList_GetItem(pos_imol_spec_list_py, 2);

			   if (PyList_Check(pos_py)) {
			      if (PyLong_Check(imol_py)) {
				 if (PyList_Check(spec_py)) {

				    PyObject *x_py = PyList_GetItem(pos_py, 0);
				    PyObject *y_py = PyList_GetItem(pos_py, 1);
				    PyObject *z_py = PyList_GetItem(pos_py, 2);

				    int imol = PyLong_AsLong(imol_py);

				    coot::residue_spec_t res_spec =
				       residue_spec_from_py(spec_py);

				    clustered_feature_info_from_python cwi(imol, res_spec, cluster_number);
				    // store cwi 
				    pharmacophore[type].push_back(cwi);

				 } else {
				    std::cout << "ERROR:: spec_py is not a list" << std::endl;
				 }
			      } else {
				 std::cout << "ERROR:: imol_py is not an int" << std::endl;
			      }
			   } else {
			      std::cout << "ERROR:: pos_py is not a list" << std::endl;
			   } 
			} else {
			   std::cout << "ERROR:: n_pos_imol_spec_list is not 3" << std::endl;
			}
		     } else {
			std::cout << "ERROR:: cluster_number_py is not an int" << std::endl;
		     } 
		  } else {
		     std::cout << "ERROR:: pos_imol_spec_list_py is not a list " << std::endl;
		  } 
	       } else {
		  std::cout << "ERROR:: pharm_py is a list of length " << n_py << std::endl;
	       } 
	    }
	 }
      }
   }
}

// constructor
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
   // a vector and std::find().
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

#include "c-interface-generic-objects.h"

// return the generic display object index
int
cfc::extracted_cluster_info_from_python::show_water_balls(unsigned int site_number) const {

   graphics_info_t g;
   std::string s = "CFC Site " + coot::util::int_to_string(site_number) + " conserved waters";

   int water_balls_object = new_generic_object_number(s);
   meshed_generic_display_object &obj = g.generic_display_objects[water_balls_object];

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
	    meshed_generic_display_object::sphere_t sphere(wc[i].pos, radius);
	    sphere.col = glm::vec4(0.9, 0.2, 0.2, 1.0);
	    obj.add(sphere);
	 }
      }
   }

   // obj.is_displayed_flag = true;
   // graphics_info_t::generic_objects_p->push_back(obj);

   Material material;
   obj.mesh.setup(&g.shader_for_moleculestotriangles, material); // fast return if already donen
   set_display_generic_object(water_balls_object, 1);

   return water_balls_object;
   
}

void
cfc::cfc_dialog_add_site_info(unsigned int site_number,
			      const extracted_cluster_info_from_python &eci) {

   GtkWidget *sites_table = lookup_widget(graphics_info_t::cfc_dialog,
					  "cfc_sites_table");
   if (sites_table) {
      gtk_table_resize(GTK_TABLE(sites_table), site_number+1, 4);
      unsigned int n_structures = eci.n_pharmacophore_structures();
      std::string structures_word = " structures";
      if (n_structures == 1)
	 structures_word = " structure";
      std::string button_label = "Site ";
      button_label += coot::util::int_to_string(site_number+1);
      // info for callback
      std::pair<int, clipper::Coord_orth> *site_no_pos_p = NULL;
      std::pair<bool, clipper::Coord_orth> pos = eci.pharmacophores_centre();
      if (pos.first) {
	 site_no_pos_p = new std::pair<int, clipper::Coord_orth> (site_number, pos.second);
      }
      GtkWidget *site_button = gtk_button_new_with_label(button_label.c_str());
      GtkWidget *label_1 = gtk_label_new(" contributed to by ");
      GtkWidget *label_2 = gtk_label_new(coot::util::int_to_string(n_structures).c_str());
      GtkWidget *label_3 = gtk_label_new(structures_word.c_str());
      g_signal_connect(G_OBJECT(site_button), "clicked",
                       G_CALLBACK(on_cfc_site_button_clicked),
                       (gpointer) site_no_pos_p);
      gtk_table_attach(GTK_TABLE(sites_table), site_button,
		       0, 1, site_number, site_number+1,
		       (GtkAttachOptions) (GTK_FILL),
		       (GtkAttachOptions) (0), 0, 0);
      gtk_table_attach(GTK_TABLE(sites_table), label_1,
		       1, 2, site_number, site_number+1,
		       (GtkAttachOptions) (GTK_FILL),
		       (GtkAttachOptions) (0), 0, 0);
      gtk_table_attach(GTK_TABLE(sites_table), label_2,
		       2, 3, site_number, site_number+1,
		       (GtkAttachOptions) (GTK_FILL),
		       (GtkAttachOptions) (0), 0, 0);
      gtk_table_attach(GTK_TABLE(sites_table), label_3,
		       3, 4, site_number, site_number+1,
		       (GtkAttachOptions) (GTK_FILL),
		       (GtkAttachOptions) (0), 0, 0);
      gtk_widget_show(site_button);
      gtk_widget_show(label_1);
      gtk_widget_show(label_2);
      gtk_widget_show(label_3);
   }
}


void
cfc::chemical_feature_clusters_add_site_info(unsigned int site_number,
					     const extracted_cluster_info_from_python &eci,
					     GtkWidget *cfc_dialog) {

   cfc_dialog = graphics_info_t::cfc_dialog;
   cfc_dialog_add_waters(site_number, eci, cfc_dialog);
   cfc_dialog_add_pharmacophores(site_number, eci, cfc_dialog);
   cfc_dialog_add_site_info(site_number, eci);

   gtk_window_set_default_size(GTK_WINDOW(cfc_dialog), 600, 400);
   gtk_widget_show(cfc_dialog);
}


#endif // MAKE_ENHANCED_LIGAND_TOOLS


#endif // USE_PYTHON
