/* src/graphics-info.cc
 * 
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Copyright 2008  by The University of Oxford
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
#include <GL/glut.h>  // for some reason...  // Eh?

#include <iostream>

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coords/mmdb-crystal.h"
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "clipper/core/map_utils.h" // Map_stats
#include "skeleton/graphical_skel.h"


#include "interface.h"

#include "molecule-class-info.h"
// #include "rama_plot.hh"
#include "skeleton/BuildCas.h"


#include "gl-matrix.h" // for baton rotation
#include "trackball.h" // for baton rotation

// #include "rottrans-buttons.hh"  old and deletable.

#include "analysis/bfkurt.hh"

#include "globjects.h"
#ifdef USE_DUNBRACK_ROTAMERS
#include "ligand/dunbrack.hh"
#else 
#include "ligand/richardson-rotamer.hh"
#endif 
#include "ligand/ligand.hh"
#include "graphics-info.h"

#include "coot-utils/coot-map-utils.hh"
#include "geometry-graphs.hh"

// Validation stuff	    //



#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

void
coot::set_validation_graph(int imol, coot::geometry_graph_type type, GtkWidget *dialog) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {

      bool found = 0; 
      if (type == coot::GEOMETRY_GRAPH_GEOMETRY) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.geometry_graph = dialog;
      }
      if (type == coot::GEOMETRY_GRAPH_B_FACTOR) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.b_factor_variance_graph = dialog;
      }
////B
      if (type == coot::GEOMETRY_GRAPH_CALC_B_FACTOR) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.b_factor_graph = dialog;
      }
////E
      if (type == coot::GEOMETRY_GRAPH_DENSITY_FIT) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.residue_density_fit_graph = dialog;
      }
      if (type == coot::GEOMETRY_GRAPH_OMEGA_DISTORTION) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.omega_distortion_graph = dialog;
      }
      if (type == coot::GEOMETRY_GRAPH_ROTAMER) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.rotamer_graph = dialog;
      }
      if (type == coot::GEOMETRY_GRAPH_NCS_DIFFS) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.ncs_diffs_graph = dialog;
      }
      if (type == coot::SEQUENCE_VIEW) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.sequence_view_is_displayed = dialog;
      }
      if (type == coot::RAMACHANDRAN_PLOT) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.dynarama_is_displayed = dialog;
      }


      if (!found) {
	 std::cout << "ERROR:: graph type " << type << " not found " << std::endl;
      }
      
   } else {
      std::cout << "WARNING:: set_validation_graph no valid molecule for imol = "
		<< imol << std::endl;
   } 
} 
#endif // defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)


#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

GtkWidget *
coot::get_validation_graph(int imol, coot::geometry_graph_type type) {

   GtkWidget *w = 0;
   if (graphics_info_t::is_valid_model_molecule(imol)) {
      bool found = 1; 
	switch(type){
	case coot::GEOMETRY_GRAPH_GEOMETRY:
	   w = graphics_info_t::molecules[imol].validation_graphs.geometry_graph;
	   break;
	case coot::GEOMETRY_GRAPH_B_FACTOR:
	   w = graphics_info_t::molecules[imol].validation_graphs.b_factor_variance_graph;
	   break;
	case coot::GEOMETRY_GRAPH_CALC_B_FACTOR:
	   w = graphics_info_t::molecules[imol].validation_graphs.b_factor_graph;
	   break;
	case coot::GEOMETRY_GRAPH_DENSITY_FIT:
	   w = graphics_info_t::molecules[imol].validation_graphs.residue_density_fit_graph;
	   break;
	case coot::GEOMETRY_GRAPH_OMEGA_DISTORTION:
	   w = graphics_info_t::molecules[imol].validation_graphs.omega_distortion_graph;
	   break;
	case coot::GEOMETRY_GRAPH_ROTAMER:
	   w = graphics_info_t::molecules[imol].validation_graphs.rotamer_graph;
	   break;
	case coot::GEOMETRY_GRAPH_NCS_DIFFS:
	   w = graphics_info_t::molecules[imol].validation_graphs.ncs_diffs_graph;
	   break;
	case coot::SEQUENCE_VIEW:
	   w = graphics_info_t::molecules[imol].validation_graphs.sequence_view_is_displayed;
	   break;
	case coot::RAMACHANDRAN_PLOT:
	   w = graphics_info_t::molecules[imol].validation_graphs.dynarama_is_displayed;
	   break;
	default:
	   found=0;
	   break;
	}

   }
   return w;
} 
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)


// convenience function
void
graphics_info_t::update_geometry_graphs(int imol) {
   update_geometry_graphs(molecules[imol].atom_sel, imol);
}


// imol map is passed in case that the density fit graph was displayed.
// If there is no imol_map then the geometry graph could not have been displayed.  You can
// pass -1 for the map in that case.
void
graphics_info_t::update_geometry_graphs(mmdb::PResidue *SelResidues, int nSelResidues, int imol, int imol_map) { // searching for update_validation_graphs? Check the next function also

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   
   GtkWidget *graph = coot::get_validation_graph(imol, coot::GEOMETRY_GRAPH_ROTAMER);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get rotamer_graph from dialog\n";
      } else {
	 std::vector<coot::geometry_graph_block_info_generic> dv =
	    rotamers_from_residue_selection(SelResidues, nSelResidues, imol);
	 gr->update_residue_blocks(dv);
      }
   }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL   
}

#include "nsv.hh"

// The molecule-based version of the above.
void
graphics_info_t::update_geometry_graphs(const atom_selection_container_t &moving_atoms_asc_local,  // searching for update_validation_graphs?
					int imol_moving_atoms) {

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

   GtkWidget *graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_GEOMETRY);
   if (graph) {
      // get deviations and replace those positions in the graph:
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get geometry_graph from dialog\n";
      } else {
	 bool with_nbcs = false;
	 std::vector<coot::geometry_distortion_info_container_t> dv =
	    geometric_distortions_from_mol(imol_moving_atoms, moving_atoms_asc_local, with_nbcs);
	 for(unsigned int ich=0; ich<dv.size(); ich++)
// 	    std::cout << "       ich " << ich << " residue blocks for updating:\n"
// 		      << dv[ich] << std::endl;
	 for(unsigned int ich=0; ich<dv.size(); ich++) 
	    gr->update_residue_blocks(dv[ich]);
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_DENSITY_FIT);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get residue_density_fit_graph from dialog\n";
      } else {
	 // Imol_Refinement_Map should be set by the time we get
	 // here.  There may be a pathological case where it has
	 // been closed when we get here.  density_fit_from_mol checks that.
	 std::vector<coot::geometry_graph_block_info_generic> dv =
	    density_fit_from_mol(moving_atoms_asc_local,
				 imol_moving_atoms,
				 Imol_Refinement_Map());
	 gr->update_residue_blocks(dv);
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_ROTAMER);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get rotamer_graph from dialog\n";
      } else {
	 std::vector<coot::geometry_graph_block_info_generic> dv =
	    rotamers_from_mol(moving_atoms_asc_local, imol_moving_atoms);

	 gr->update_residue_blocks(dv);
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_NCS_DIFFS);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get rotamer_graph from dialog\n";
      } else {
	 std::vector<coot::geometry_graph_block_info_generic> dv =
	    ncs_diffs_from_mol(imol_moving_atoms); // update everything
	 gr->update_residue_blocks(dv);
      }
   }
   
   graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_OMEGA_DISTORTION);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get omega_graph from dialog\n";
      } else {

	 if (! moving_atoms_asc_local.empty()) {

	    // We do this long handedly (c.f. above) because here we use
	    // render_omega_blocks() which needs the offset (which is a
	    // per-chain variable:
	    //
	    int n_models = moving_atoms_asc_local.mol->GetNumberOfModels();
	    for (int imodel = 1; imodel <= n_models; imodel++) { 
	       mmdb::Model *model_p = moving_atoms_asc_local.mol->GetModel(imodel);
	       mmdb::Chain *chain_p;
	       const char *chain_id;
	       int n_chains = model_p->GetNumberOfChains();

	       for (int ich=0; ich<n_chains; ich++) {
		  chain_p = model_p->GetChain(ich);
		  chain_id = chain_p->GetChainID();
		  std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
		  if (m.first) {
		     // not used:
		     // int offset = m.second - 1; // min resno = 1 -> offset = 0

		     coot::omega_distortion_info_container_t om_dist = 
			omega_distortions_from_mol(moving_atoms_asc_local, chain_id);	

		     if (0)
			std::cout << "DEBUG:: update omega dist graph chain "
				  << om_dist.chain_id << " " << om_dist.omega_distortions.size()
				  << " blocks" << std::endl;

		     gr->update_omega_blocks(om_dist, ich, std::string(chain_id));
		  }
	       }
	    }
	 }
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::SEQUENCE_VIEW);
   if (graph) {

      exptl::nsv *sequence_view = static_cast<exptl::nsv *>(g_object_get_data(G_OBJECT(graph), "nsv"));

      if (sequence_view) {
	 mmdb::Manager *mol = molecules[imol_moving_atoms].atom_sel.mol;
	 sequence_view->regenerate(mol);
      }
   }

   // and now ramachandran also


#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL
}

void
graphics_info_t::update_validation_graphs(int imol) {

#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   GtkWidget *w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
   if (w) {
      coot::rama_plot *plot = reinterpret_cast<coot::rama_plot *>(gtk_object_get_user_data(GTK_OBJECT(w)));
      std::cout << "doing handle_rama_plot_update() " << std::endl;
      handle_rama_plot_update(plot);
   }
   // now update the geometry graphs, so get the asc
   atom_selection_container_t u_asc = molecules[imol].atom_sel;
   update_geometry_graphs(u_asc, imol);
#endif // HAVE_GTK_CANVAS
}



void
graphics_info_t::delete_residue_from_geometry_graphs(int imol,
						     coot::residue_spec_t res_spec) {
#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

   std::vector<coot::geometry_graph_type> graph_types;
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_GEOMETRY);
   graph_types.push_back(coot::GEOMETRY_GRAPH_B_FACTOR);
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_OMEGA_DISTORTION);
   graph_types.push_back(coot::GEOMETRY_GRAPH_ROTAMER);
   graph_types.push_back(coot::GEOMETRY_GRAPH_NCS_DIFFS);

   for (unsigned int igt=0; igt<graph_types.size(); igt++) { 
      GtkWidget *graph =
	 coot::get_validation_graph(imol_moving_atoms, graph_types[igt]);
      if (graph) {
	 coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
	 if (gr) {
	    gr->delete_block(res_spec.chain_id, res_spec.res_no);
	 }
      }
   }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL
}

void
graphics_info_t::delete_residues_from_geometry_graphs(int imol,
						      const std::vector<coot::residue_spec_t> &res_specs) {

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

   std::vector<coot::geometry_graph_type> graph_types;
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_GEOMETRY);
   graph_types.push_back(coot::GEOMETRY_GRAPH_B_FACTOR);
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_OMEGA_DISTORTION);
   graph_types.push_back(coot::GEOMETRY_GRAPH_ROTAMER);
   graph_types.push_back(coot::GEOMETRY_GRAPH_NCS_DIFFS);

   for (unsigned int igt=0; igt<graph_types.size(); igt++) {
      GtkWidget *graph =
	 coot::get_validation_graph(imol_moving_atoms, graph_types[igt]);
      if (graph) {
	 coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
	 if (gr) {

	    for (std::size_t ires=0; ires<res_specs.size(); ires++) {
	       const coot::residue_spec_t &res_spec = res_specs[ires];
	       gr->delete_block(res_spec.chain_id, res_spec.res_no);
	    }
	 }
      }
   }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL
}

void
graphics_info_t::delete_chain_from_geometry_graphs(int imol, const std::string &chain_id) {

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

   std::vector<coot::geometry_graph_type> graph_types;
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_GEOMETRY);
   graph_types.push_back(coot::GEOMETRY_GRAPH_B_FACTOR);
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_OMEGA_DISTORTION);
   graph_types.push_back(coot::GEOMETRY_GRAPH_ROTAMER);
   graph_types.push_back(coot::GEOMETRY_GRAPH_NCS_DIFFS);

   for (unsigned int igt=0; igt<graph_types.size(); igt++) { 
      GtkWidget *graph =
	 coot::get_validation_graph(imol_moving_atoms, graph_types[igt]);
      if (graph) {
	 coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
	 if (gr) {
	    // for res_spec in residues of chain that has bee deleted...
	    // gr->delete_block(res_spec.chain_id, res_spec.res_no);
	 }
      }
   }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL
}


void
graphics_info_t::geometric_distortion(int imol) {

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

   // we need to assign these
//    int resno_1;
//    int resno_2;
   std::string chain_id_1;
   
   //    short int irest = 0;  // make 1 if restraints were found

   // make the selection and build a new molecule inside restraints.

   // short int have_flanking_residue_at_start = 0; 
   // short int have_flanking_residue_at_end = 0;
   // short int have_disulfide_residues = 0;  // other residues are included in the
   // residues_mol for disphide
   // restraints.
   
   // 9 Sept 2003: The atom selection goes mad if residue with seqnum
   // iend_res+1 does not exist, but is not at the end of the chain.

   // Therefore we will set 2 flags, which tell us if istart_res-1 and
   // iend_res+1 exist.  And we do that by trying to select atoms from
   // them - if they exist, the number of selected atoms will be more
   // than 0.

//    istart_minus_flag = 0;  // from simple restraint code
//    iend_plus_flag    = 0;

   mmdb::Manager *mol = molecules[imol].atom_sel.mol; // short-hand usage

   // This dcv used to be inside the mol test, but that tickled what I
   // believe a compiler bug (on deconstructing this vector).  So now
   // it is here and coot doesn't crash when doing geometry analysis
   // of NMR model(s).  (In minimal testing).
   
   std::vector<coot::geometry_distortion_info_container_t> dcv;

   if (mol) {

      // big copy:
      bool with_nbcs = false;
      dcv = geometric_distortions_from_mol(imol, molecules[imol].atom_sel, with_nbcs);

      int max_chain_length = coot::util::max_min_max_residue_range(mol);
      if (max_chain_length <= 0) {
	 std::cout << "WARNING:: Funny coords - no graphs for this molecule" << std::endl;
      } else { 
	 int nchains = coot::util::number_of_chains(mol);

	 std::string name = graphics_info_t::molecules[imol].name_for_display_manager();
	 coot::geometry_graphs *graphs = new coot::geometry_graphs(coot::GEOMETRY_GRAPH_GEOMETRY,
								   imol, name, nchains, max_chain_length);
	 coot::set_validation_graph(imol, coot::GEOMETRY_GRAPH_GEOMETRY, graphs->dialog()); // store for potential updates

	 // debugging.  Problem was negative occupancies.
// 	 for(unsigned int i=0; i<dcv.size(); i++) {
// 	    std::cout << i << " chain: " << dcv[i].chain_id << std::endl;
// 	    for(unsigned int j=0; j<dcv[i].geometry_distortion.size(); j++) {
// 	       std::cout << j << " " << dcv[i].geometry_distortion[j].distortion_score
// 			 << std::endl;
// 	    }
// 	 }

	 for(unsigned int i=0; i<dcv.size(); i++) { 
	    graphs->render_to_canvas(dcv[i], i);
	 }
      }
   }
   
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif
}

#ifdef HAVE_GSL
coot::geometry_distortion_info_container_t
graphics_info_t::geometric_distortions(int imol, mmdb::Residue *residue_p, bool with_nbcs) {

   coot::geometry_distortion_info_container_t gdc(NULL, 0, "");
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

   if (residue_p) {
      mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
      if (mol) {
	 atom_selection_container_t asc = make_asc(mol);
	 std::vector<coot::geometry_distortion_info_container_t> v =
	    geometric_distortions_from_mol(imol, asc, with_nbcs);
	 if (v.size() == 1) {
	    if (v[0].geometry_distortion.size() > 1) {
	       gdc = v[0];
	    }
	 } 
	 asc.clear_up();
      }
   }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   return gdc;
}
#endif // HAVE_GSL


#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
std::vector<coot::geometry_distortion_info_container_t>
graphics_info_t::geometric_distortions_from_mol(int imol, const atom_selection_container_t &asc,
						bool with_nbcs) {

   std::vector<coot::geometry_distortion_info_container_t> dcv;
   std::string altconf("");  // use this (e.g. "A") or "".

   if (! asc.mol)
      return dcv;
   
   int n_models = asc.mol->GetNumberOfModels();
   
   if (n_models > 0) { 

      // 20100629, crash!  Ooops, we can't run over many models
      // because geometry_graphs (and its data member `blocks' are
      // sized to the number of chains).  If we run over all models,
      // then there are too many chains for the indexing of `blocks`
      // -> crash.  So we just use the first model.
      
      // for (int imod=1; imod<=n_models; imod++) {
      int imod=1;
      {
	 
	 mmdb::Model *model_p = asc.mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 const char *chain_id;
	 
	 
	 for (int ichain=0; ichain<nchains; ichain++) {
	    
	    chain_p = model_p->GetChain(ichain);

	    if (! chain_p->isSolventChain()) { 
	       chain_id = chain_p->GetChainID();
	    
	       // First make an atom selection of the residues selected to regularize.
	       // 
	       int selHnd = asc.mol->NewSelection(); // yes, it's deleted.
	       int nSelResidues;
	       mmdb::PResidue *SelResidues = NULL;
	    
	       // Consider as the altconf the altconf of one of the residues (we
	       // must test that the altlocs of the selected atoms to be either
	       // the same as each other (A = A) or one of them is "".  We need to
	       // know the mmdb syntax for "either".  Well, now I know that's ",A"
	       // (for either blank or "A").
	       // 
	       // 
	       // 
	       asc.mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
			       chain_id,
			       mmdb::ANY_RES, "*",
			       mmdb::ANY_RES, "*",
			       "*",  // residue name
			       "*",  // Residue must contain this atom name?
			       "*",  // Residue must contain this Element?
			       "*",  // altLocs
			       mmdb::SKEY_NEW // selection key
			       );
	       asc.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
	       std::pair<int, std::vector<std::string> > icheck = 
		  check_dictionary_for_residue_restraints(imol, SelResidues, nSelResidues);
	    
	       if (icheck.first == 0) { 
		  for (unsigned int icheck_res=0; icheck_res<icheck.second.size(); icheck_res++) { 
		     std::cout << "WARNING:: Failed to find restraints for " 
			       << icheck.second[icheck_res] << std::endl;
		  }
	       }
	    
	       std::cout << "INFO:: " << nSelResidues 
			 << " residues selected for geometry checking object" << std::endl;
	    
	       if (nSelResidues <= 0) {

		  std::cout << "ERROR:: No Residues!!   This should never happen:" << std::endl;
		  std::cout << "  in create_regularized_graphical_object" << std::endl;

	       } else { // normal
	       
		  std::vector<mmdb::Atom *> fixed_atoms;
		  std::vector<coot::atom_spec_t> fixed_atom_specs;
	       
		  // Notice that we have to make 2 atom selections, one, which includes
		  // flanking (and disulphide eventually) residues that is used for the
		  // restraints (restraints_container_t constructor) and one that is the
		  // moving atoms (which does not have flanking atoms).
		  // 
		  // The restraints_container_t moves the atom of the mol that is passes to
		  // it.  This must be the same mol as the moving atoms mol so that the
		  // changed atom positions can be seen.  However (as I said) the moving
		  // atom mol should not have flanking residues shown.  So we make an asc
		  // that has the same mol as that passed to the restraints but a different
		  // atom selection (it is the atom selection that is used in the bond
		  // generation).
		  // 

//	          20100210 try vector
// 		  coot::restraints_container_t restraints(SelResidues, nSelResidues,
// 							  std::string(chain_id),
// 							  asc.mol);
		  std::vector<std::pair<bool,mmdb::Residue *> > residue_vec;
		  for (int ires=0; ires<nSelResidues; ires++)
		     residue_vec.push_back(std::pair<bool, mmdb::Residue *> (0, SelResidues[ires]));

		  std::vector<mmdb::Link> links;
		  clipper::Xmap<float> dummy_xmap;

		  coot::restraints_container_t restraints(residue_vec,
							  links,
							  *Geom_p(),
							  asc.mol,
							  fixed_atom_specs,
							  &dummy_xmap);
	       
		  // coot::restraint_usage_Flags flags = coot::BONDS;
		  // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES;
		  // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
		  // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES; 
		  coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
		  flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
		  flags = coot::BONDS_ANGLES_AND_PLANES;
		  flags = coot::BONDS_ANGLES_PLANES_AND_CHIRALS;

		  if (with_nbcs)
		     flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
		  
		  unsigned int n_threads = coot::get_max_number_of_threads();
		  if (n_threads > 0)
		     restraints.thread_pool(&static_thread_pool, n_threads);
		  short int do_residue_internal_torsions = 0;
	       
		  // 	       if (do_torsion_restraints) { 
		  // 		  do_residue_internal_torsions = 1;
		  // 		  flags = coot::BONDS_ANGLES_TORSIONS_PLANES_AND_NON_BONDED;
		  // 	       } 
	       
		  // 	       if (do_peptide_torsion_restraints)
		  // 		  do_link_torsions = 1;

		  coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
		  bool do_trans_peptide_restraints = false;
		  int nrestraints = 
		     restraints.make_restraints(imol, *geom_p,
						flags,
						do_residue_internal_torsions,
						do_trans_peptide_restraints,
						0.0, 0, false, false, false,
						pseudos);

		  if (nrestraints > 0) {

// 		     std::cout << "DEBUG:: model " << imod << " pushing back " << nrestraints
// 			       << " restraints" << std::endl;

		     dcv.push_back(restraints.geometric_distortions());
		  
		  } else {

		     // don't give this annoying dialog if restraints
		     // have been read for the residues in this chain.
		     // e.g. a single CLs residues in a chain.
		     std::vector<std::string> res_types = coot::util::residue_types_in_chain(chain_p);
		     bool hd = geom_p->have_dictionary_for_residue_types(res_types, imol,
									 cif_dictionary_read_number);
		     cif_dictionary_read_number += res_types.size();
		     if (! hd) {
			if (use_graphics_interface_flag) { 
			   GtkWidget *widget = create_no_restraints_info_dialog();
			   gtk_widget_show(widget);
			} else {
			   std::cout << "WARNING:: No dictionary for some residue types " << std::endl;
			} 
		     }
		  }
	       }
	       asc.mol->DeleteSelection(selHnd);
	    }
	 }
      }
   }
   // print_geometry_distortion(dcv);
   return dcv;
}
#endif // HAVE_GSL
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

#ifdef HAVE_GSL   
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
void
graphics_info_t::print_geometry_distortion(const std::vector<coot::geometry_distortion_info_container_t> &v) const {
   for (unsigned int i=0; i<v.size(); i++) { 
      std::cout << v[i] << "\n";
   }
}
#endif // HAVE_GSL
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)


////B
void
graphics_info_t::calc_b_factor_graphs(int imol) {

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   
   if (imol<n_molecules())
      if (imol >= 0)
	 if (molecules[imol].has_model()) {
	    mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	    bool is_shelx_mol = molecules[imol].is_from_shelx_ins();

	    coot_extras::b_factor_analysis bfa(mol, is_shelx_mol);
	    std::vector<coot_extras::my_chain_of_stats_t> bfa_chain_info =
	       bfa.chain_details();

	    int n_models = mol->GetNumberOfModels();
	    int max_chain_length = coot::util::max_min_max_residue_range(mol);
	    if (max_chain_length <= 0) {
	       std::cout << "WARNING:: Funny coords - no graphs for this molecule"
			 << std::endl;
	    } else { 
	       for (int imodel = 1; imodel <= n_models; imodel++) { 
		  mmdb::Model *model_p = mol->GetModel(imodel);
		  mmdb::Chain *chain_p;
		  unsigned int n_chains = model_p->GetNumberOfChains();
		  coot::geometry_graphs *graphs =
		     new coot::geometry_graphs(coot::GEOMETRY_GRAPH_CALC_B_FACTOR,
					       imol,
					       graphics_info_t::molecules[imol].name_for_display_manager(), 
					       n_chains, max_chain_length);
		  // b_factor_variance_graph[imol] = graphs->dialog();
		  set_validation_graph(imol, coot::GEOMETRY_GRAPH_CALC_B_FACTOR, graphs->dialog());

		  coot::b_factor_block_info_t bfi[3];
		  float std_dev;
		  int offset;

		  for (unsigned int ich=0; ich<n_chains; ich++) {

		     if (ich < bfa_chain_info.size()) { 
			chain_p = model_p->GetChain(ich);
			std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
			
			if (m.first) { 
			   std::vector<coot::b_factor_block_info_t> bfiv;
			   offset = m.second - 1; 
			  for (unsigned int ires=0; ires<bfa_chain_info[ich].residue_properties.size(); ires++) {

				double variance[2]={0.0,0.0};
////	B FACTOR CALC FOR MAIN AND SIDECHAIN
				mmdb::PPAtom residue_atoms;
				mmdb::PResidue residue_p = chain_p->GetResidue(ires);
   				int nResidueAtoms=0;
				double running_sum[4]	= {0.0,0.0,0.0,0.0}; 
 				double mean[3]		= {0.0,0.0,0.0};
				double std_dev[2]	= {0.0,0.0};
				double bf  	= 0.0;
				double occ 	= 0.0;
				double bfo 	= 0.0;
 				double div[2] 	= {0.0,0.0};
   				int    bMC 	= 0;

				residue_p->GetAtomTable(residue_atoms, nResidueAtoms); 
				if (nResidueAtoms > 0) { 
					for (int i=0; i<nResidueAtoms; i++) {
						bMC = ( coot::is_main_chain_p(residue_atoms[i]) )?1:0;
						std::string ele = residue_atoms[i]->element;
						if ((ele != " H") && (ele != " D")) {
							bf  = residue_atoms[i]->tempFactor;
							occ = residue_atoms[i]->occupancy;
							if ( ((bf > 0.0) && (occ >= 0.0) && (occ <= 1.0)) ||
							(is_shelx_mol && (occ < 11.001) && (occ > 10.999))) {
								if (is_shelx_mol)
									occ = 1.0; 
								div[bMC] 		+= occ;
								bfo			 = bf*occ;
								running_sum[0+2*bMC] 	+= bfo;
								running_sum[1+2*bMC] 	+= bfo*bfo;
							}
						}
					}
					if ( div[0] > 0 || div[1] > 0 ) { 
						mean[0]     = (div[0]>0)?(running_sum[0]/div[0]):(0.0);		// notMC
						mean[1]     = (running_sum[0]+running_sum[2])/(div[0]+div[1]);
					}
				}else{ //SO THIS SHOULD NEVER HAPPEN WILL JUST YIELD ZERO PLOTS
					std::cout << "ERROR::  IN B FACTOR CALCULATION, EMPTY RESIDUE" << std::endl;
				}
////	END AVERAGE B FACTOR
////	SIDECHAIN <B>_sc
			      bfi[0].resno = bfa_chain_info[ich].residue_properties[ires].resno;
			      bfi[0].b_factor_var = mean[0]+mean[1]; // IN ORDER TO PLOT SC ON THE TOP
			      bfi[0].info_string  = int_to_string(bfi[0].resno);
			      bfi[0].info_string += chain_p->GetChainID();
			      bfi[0].info_string += " ";
			      bfi[0].info_string += bfa_chain_info[ich].residue_properties[ires].resname;
			      bfi[0].info_string += ":SC: ";
			      bfi[0].info_string += float_to_string(bfi[0].b_factor_var - mean[1]); // USER SEES AVERAGE SC B FACTOR
			      bfi[0].atom_name = bfa_chain_info[ich].residue_properties[ires].atom_name;
////	TOTAL <B>_tot
		 	      bfi[1].resno = bfa_chain_info[ich].residue_properties[ires].resno;
			      bfi[1].b_factor_var = mean[1];
			      bfi[1].info_string  = int_to_string(bfi[1].resno);
			      bfi[1].info_string += chain_p->GetChainID();
			      bfi[1].info_string += " ";
			      bfi[1].info_string += bfa_chain_info[ich].residue_properties[ires].resname;
			      bfi[1].info_string += ":TOT: ";
			      bfi[1].info_string += float_to_string(bfi[1].b_factor_var); // USER SEES THE ACTUAL VALUE SINCE IT IS IN THE FRONT
			      bfi[1].atom_name = bfa_chain_info[ich].residue_properties[ires].atom_name;

			      bfiv.push_back(bfi[0]); //SC
			      bfiv.push_back(bfi[1]); //TOT
			 }
			 graphs->render_b_factor_blocks(imol, ich, bfa_chain_info[ich].chain_id,
							  offset, bfiv);
			}
		     }
		  }
	       }
	    }
	 }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL  
}
////E

void
graphics_info_t::b_factor_graphs(int imol) {

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   
   if (imol<n_molecules())
      if (imol >= 0)
	 if (molecules[imol].has_model()) {
	    mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	    bool is_shelx_mol = molecules[imol].is_from_shelx_ins();

	    coot_extras::b_factor_analysis bfa(mol, is_shelx_mol);
	    std::vector<coot_extras::my_chain_of_stats_t> bfa_chain_info =
	       bfa.chain_details();

	    int n_models = mol->GetNumberOfModels();
	    int max_chain_length = coot::util::max_min_max_residue_range(mol);
	    if (max_chain_length <= 0) {
	       std::cout << "WARNING:: Funny coords - no graphs for this molecule"
			 << std::endl;
	    } else { 
	       for (int imodel = 1; imodel <= n_models; imodel++) { 
		  mmdb::Model *model_p = mol->GetModel(imodel);
		  mmdb::Chain *chain_p;
		  unsigned int n_chains = model_p->GetNumberOfChains();
		  coot::geometry_graphs *graphs =
		     new coot::geometry_graphs(coot::GEOMETRY_GRAPH_B_FACTOR,
					       imol,
					       graphics_info_t::molecules[imol].name_for_display_manager(), 
					       n_chains, max_chain_length);
		  // b_factor_variance_graph[imol] = graphs->dialog();
		  set_validation_graph(imol, coot::GEOMETRY_GRAPH_B_FACTOR, graphs->dialog());

		  coot::b_factor_block_info_t bfi;
		  float std_dev;
		  int offset;

		  for (unsigned int ich=0; ich<n_chains; ich++) {

		     if (ich < bfa_chain_info.size()) { 
			chain_p = model_p->GetChain(ich);
			std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
			
			if (m.first) { 
			   std::vector<coot::b_factor_block_info_t> bfiv;
			   offset = m.second - 1;
			   
			   for (unsigned int ires=0; ires<bfa_chain_info[ich].residue_properties.size(); ires++) {
			      bfi.resno = bfa_chain_info[ich].residue_properties[ires].resno;
			      std_dev = bfa_chain_info[ich].residue_properties[ires].std_dev;
			      bfi.b_factor_var = std_dev * std_dev;
			      bfi.info_string  = int_to_string(bfi.resno);
			      bfi.info_string += chain_p->GetChainID();
			      bfi.info_string += " ";
			      bfi.info_string += bfa_chain_info[ich].residue_properties[ires].resname;
			      bfi.info_string += ": ";
			      bfi.info_string += float_to_string(bfi.b_factor_var);
			      bfi.atom_name = bfa_chain_info[ich].residue_properties[ires].atom_name;
			      bfiv.push_back(bfi);
			   }
			   graphs->render_b_factor_blocks(imol, ich, bfa_chain_info[ich].chain_id,
							  offset, bfiv);
			}
		     }
		  }
	       }
	    }
	 }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL   
}

void
graphics_info_t::omega_graphs(int imol) {

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   
   if (imol >= 0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	    int n_models = mol->GetNumberOfModels();
	    int max_chain_length = coot::util::max_min_max_residue_range(mol);
	    if (max_chain_length <= 0) {
	       std::cout << "WARNING:: Funny coords - no graphs for this molecule" << std::endl;
	    } else { 
	       for (int imodel = 1; imodel <= n_models; imodel++) { 
		  mmdb::Model *model_p = mol->GetModel(imodel);
		  mmdb::Chain *chain_p;
		  const char *chain_id;
		  int n_chains = model_p->GetNumberOfChains();
		  coot::geometry_graphs *graphs =
		     new coot::geometry_graphs(coot::GEOMETRY_GRAPH_OMEGA_DISTORTION,
					       imol,
					       graphics_info_t::molecules[imol].name_for_display_manager(), 
					       n_chains, max_chain_length);
		  // omega_distortion_graph[imol] = graphs->dialog();
		  coot::set_validation_graph(imol, coot::GEOMETRY_GRAPH_OMEGA_DISTORTION, graphs->dialog());

		  for (int ich=0; ich<n_chains; ich++) {
		     chain_p = model_p->GetChain(ich);
		     chain_id = chain_p->GetChainID();
		     std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
		     if (m.first) {
			int offset = m.second - 1; // min resno = 1 -> offset = 0
			int selHnd = mol->NewSelection();
			mmdb::PResidue *SelResidues = NULL;
			int nSelResidues;

			mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
				    chain_id,
				    mmdb::ANY_RES, "*",
				    mmdb::ANY_RES, "*",
				    "*",  // residue name
				    "*",  // Residue must contain this atom name?
				    "*",  // Residue must contain this Element?
				    "*",  // altLocs
				    mmdb::SKEY_NEW // selection key
				    );
			mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
			clipper::Xmap<float> dummy_xmap;

			if (nSelResidues > 0) { 
			   coot::restraints_container_t restraints(molecules[imol].atom_sel,
								   std::string(chain_id),
								   &dummy_xmap);

			   coot::omega_distortion_info_container_t om_dist = 
			      restraints.omega_trans_distortions(*geom_p,
								 mark_cis_peptides_as_bad_flag);
			   // std::cout << "DEBUG: got om_dist." << std::endl;

			   graphs->render_omega_blocks(om_dist, ich, std::string(chain_id),
						       offset);
			}
			mol->DeleteSelection(selHnd);
		     }
		  }
	       }
	    }
	 }
      }
   }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL   
}

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
coot::omega_distortion_info_container_t 
graphics_info_t::omega_distortions_from_mol(const atom_selection_container_t &asc,
					    const std::string &chain_id) {

   clipper::Xmap<float> dummy_xmap;
   coot::restraints_container_t restraints(asc, chain_id, &dummy_xmap);
   coot::omega_distortion_info_container_t om_dist =
      restraints.omega_trans_distortions(*geom_p, mark_cis_peptides_as_bad_flag);
   return om_dist;
}
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL   

coot::rotamer_graphs_info_t 
graphics_info_t::rotamer_graphs(int imol) {

   coot::rotamer_graphs_info_t info;
      
#ifdef HAVE_GSL    
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   if (imol >= 0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	    int n_models = mol->GetNumberOfModels();
	    int max_chain_length = coot::util::max_min_max_residue_range(mol);
	    if (max_chain_length <= 0) {
	       std::cout << "WARNING:: Funny coords - no graphs for this molecule" << std::endl;
	    } else {
	       for (int imodel = 1; imodel <= n_models; imodel++) { 
		  mmdb::Model *model_p = mol->GetModel(imodel);
		  mmdb::Chain *chain_p;
		  const char *chain_id;
		  int n_chains = model_p->GetNumberOfChains();
		  coot::geometry_graphs *graphs = 0;
		  if (use_graphics_interface_flag) {
		     std::string mol_name = graphics_info_t::molecules[imol].name_for_display_manager();
		     graphs = new coot::geometry_graphs(coot::GEOMETRY_GRAPH_ROTAMER,
							imol, mol_name,
							n_chains, max_chain_length);
		     
		     // rotamer_graph[imol] = graphs->dialog();
		     coot::set_validation_graph(imol, coot::GEOMETRY_GRAPH_ROTAMER, graphs->dialog());
		  }
		  for (int ich=0; ich<n_chains; ich++) {
		     chain_p = model_p->GetChain(ich);
		     if (! chain_p->isSolventChain()) { 
			chain_id = chain_p->GetChainID();
			std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
			if (m.first) {
			   int offset = m.second - 1;
			   int selHnd = mol->NewSelection();
			   mmdb::PResidue *SelResidues = NULL;
			   int nSelResidues;

			   mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
				       chain_id,
				       mmdb::ANY_RES, "*",
				       mmdb::ANY_RES, "*",
				       "*",  // residue name
				       "*",  // Residue must contain this atom name?
				       "*",  // Residue must contain this Element?
				       "*",  // altLocs
				       mmdb::SKEY_NEW // selection key
				       );
			   mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

			   if (nSelResidues > 0) {
			      int max_resno = -9999;

			      std::string altconf = ""; // fixme, I guess
			      std::vector<coot::geometry_graph_block_info_generic> v;
			      for (int ir=0; ir<nSelResidues; ir++) {
				 int this_resno = SelResidues[ir]->GetSeqNum();
				 std::string res_name = SelResidues[ir]->GetResName();
				 if (coot::util::is_standard_amino_acid_name(res_name)) { 
				 // if (res_name != "HOH") { 
				    std::string this_inscode = SelResidues[ir]->GetInsCode();
				    if (this_resno > max_resno)
				       max_resno = this_resno;
				    coot::rotamer_probability_info_t d_score = 
				       get_rotamer_probability(SelResidues[ir], altconf, mol,
							       rotamer_lowest_probability, 1);
				 
				    double distortion = 0.0;
				    std::string str = int_to_string(this_resno);
				    str += chain_id;
				    str += " ";
				    str += SelResidues[ir]->name;
				    str += " ";
				    coot::atom_spec_t atom_spec(chain_id, this_resno,
								"", " CA ", "");
				    coot::graph_rotamer_info_t ri("", -9999, "", -99.9, "unset");
				    switch (d_score.state) {

				    case 1:
				       // On reflection we don't want found
				       // rotamers with low probabilities to be
				       // marked (nearly as or more) bad than not
				       // found rotamers:
				       distortion = 100.0 * rotamer_distortion_scale/d_score.probability;
				       distortion = distortion > 100.0 ? 100.0 : distortion;
				       str += " ";
				       str += d_score.rotamer_name;
				       str += " ";
				       str += "Probability: ";
				       str += float_to_string(d_score.probability);
				       str += "%";
				       ri = coot::graph_rotamer_info_t (chain_id, this_resno,
									this_inscode,
									d_score.probability,
									d_score.rotamer_name);
				       info.info.push_back(ri);
				       break;
				 
				    case 0:
				       distortion = 100.0;
				       str += "Missing Atoms";
				       atom_spec.string_user_data = "Missing Atoms";
				       ri = coot::graph_rotamer_info_t (chain_id, this_resno,
									this_inscode,
									0.0, "Missing Atoms");
				       info.info.push_back(ri);
				       break;
				 
				    case -1:
				       distortion = 100.0;
				       ri = coot::graph_rotamer_info_t (chain_id, this_resno,
									this_inscode,
									0.0, "Rotamer not recognised");
				       info.info.push_back(ri);
				       str += "Rotamer not recognised";
				       break;

				    case -2:   // don't plot a block for this one.
				       distortion = 0.0;
				       break;
				    }
				    if (d_score.state != -2) {
				       v.push_back(coot::geometry_graph_block_info_generic(imol, this_resno, atom_spec, distortion, str));
				    }
				 }
			      } // end residue for loop
			      // done residue loop:
// 			      std::cout << "render_to_canvas: (rotamer) chain: " << ich << " min_resno: "
// 					<< m.second << " max_resno: " << max_resno
// 					<< " offset: " << offset << std::endl;
			      if (use_graphics_interface_flag && max_resno > 0)
				 graphs->render_to_canvas(v, ich, std::string(chain_id),
							  max_resno, m.second, offset);

			   }

			   mol->DeleteSelection(selHnd);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL
   return info;
}

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::rotamers_from_mol(const atom_selection_container_t &asc,
				  int imol_moving_atoms) {

   // this does not use the provided atom_selection_container_t asc
   
   std::vector<coot::geometry_graph_block_info_generic> dv;

   mmdb::Manager *mol = molecules[imol_moving_atoms].atom_sel.mol;
   // int n_models = mol->GetNumberOfModels();
   // for (int imodel = 1; imodel <= n_models; imodel++) { 
   int imodel = 1;
   mmdb::Model *model_p = mol->GetModel(imodel);
   mmdb::Chain *chain_p;
   const char *chain_id;
   int n_chains = model_p->GetNumberOfChains();

   for (int ich=0; ich<n_chains; ich++) {
      chain_p = model_p->GetChain(ich);
      if (! chain_p->isSolventChain()) { 
	 chain_id = chain_p->GetChainID();
	 std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
	 if (m.first) {
	    // int offset = m.second - 1;
	    int selHnd = mol->NewSelection();
	    mmdb::PResidue *SelResidues = NULL;
	    int nSelResidues;

	    mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
			chain_id,
			mmdb::ANY_RES, "*",
			mmdb::ANY_RES, "*",
			"*",  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			mmdb::SKEY_NEW // selection key 
			);
	    mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

	    if (nSelResidues > 0) {
	       int max_resno = -9999;

	       for (int ir=0; ir<nSelResidues; ir++) {
		  int this_resno = SelResidues[ir]->GetSeqNum();
		  std::string res_name = SelResidues[ir]->GetResName();
		  if (res_name != "HOH") { 
		     if (this_resno > max_resno)
			max_resno = this_resno;
		     std::string altconf = "";
		     coot::rotamer_probability_info_t d_score = 
			get_rotamer_probability(SelResidues[ir], altconf, mol,
						rotamer_lowest_probability, 1);
		     double distortion = 0.0;
		     std::string str = int_to_string(this_resno);
		     str += chain_id;
		     str += " ";
		     str += SelResidues[ir]->name;
		     str += " ";
		     coot::atom_spec_t atom_spec(chain_id, this_resno,
						 "", " CA ", "");
		     switch (d_score.state) {

		     case 1:
			// On reflection we don't want found
			// rotamers with low probabilities to be
			// marked (nearly as or more) bad than not
			// found rotamers:
			distortion = 100.0 * rotamer_distortion_scale/d_score.probability;
			distortion = distortion > 100.0 ? 100.0 : distortion;
			str += " ";
			str += d_score.rotamer_name;
			str += " ";
			str += "Probability: ";
			str += float_to_string(d_score.probability);
			str += "%";
			break;
				 
		     case 0:
			distortion = 100.0;
			str += "Missing Atoms";
			atom_spec.string_user_data = "Missing Atoms";
			break;
				 
		     case -1:
			distortion = 100.0;
			str += "Rotamer not recognised";
			break;

		     case -2:   // don't plot a block for this one.
			distortion = 0.0;
			break;
		     }
		     if (d_score.state != -2) {
			dv.push_back(coot::geometry_graph_block_info_generic(imol_moving_atoms, this_resno, atom_spec, distortion, str));
		     }
		  }
	       }
	    }
	    mol->DeleteSelection(selHnd);
	 }
      }
   }
   return dv;
}
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL


#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::rotamers_from_residue_selection(mmdb::PResidue *SelResidues,
						 int nSelResidues, int imol) {
   
   std::vector<coot::geometry_graph_block_info_generic> v;
   for (int ires=0; ires<nSelResidues; ires++) {
      std::string res_name = SelResidues[ires]->GetResName();
      if (res_name != "HOH") {
	 int this_resno = SelResidues[ires]->GetSeqNum();
	 std::string chain_id = SelResidues[ires]->GetChainID();
	 std::string alt_conf = ""; // fixme?
	 coot::rotamer_probability_info_t d_score =
	    get_rotamer_probability(SelResidues[ires], alt_conf, 0, // 0? woo..
				    rotamer_lowest_probability, 1);
	 double distortion = 0.0;
	 std::string str = int_to_string(this_resno);
	 str += chain_id;
	 str += " ";
	 str += SelResidues[ires]->name;
	 str += " ";
	 switch (d_score.state) {

	 case 1:
	    // On reflection we don't want found
	    // rotamers with low probabilities to be
	    // marked (nearly as or more) bad than not
	    // found rotamers:
	    distortion = 100.0 * rotamer_distortion_scale/d_score.probability;
	    distortion = distortion > 100.0 ? 100.0 : distortion;
	    str += " ";
	    str += d_score.rotamer_name;
	    str += " ";
	    str += "Probability: ";
	    str += float_to_string(d_score.probability);
	    str += "%";
	    break;
				 
	 case 0:
	    distortion = 100.0;
	    str += "Missing Atoms";
	    break;
				 
	 case -1:
	    distortion = 100.0;
	    str += "Rotamer not recognised";
	    break;

	 case -2:   // don't plot a block for this one.
	    distortion = 0.0;
	    break;
	 }
	 if (d_score.state != -2) {
	    coot::atom_spec_t atom_spec(chain_id, this_resno,
					"", " CA ", "");
	    v.push_back(coot::geometry_graph_block_info_generic(imol_moving_atoms,
								this_resno,
								atom_spec,
								distortion,
								str));
	 }
      }
   }
   return v;
}
#endif //  defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GS


void
graphics_info_t::density_fit_graphs(int imol) {

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
   if (imol >= 0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    
	    int imol_for_map = Imol_Refinement_Map();
	    if (imol_for_map == -1)
	       show_select_map_dialog();
            // maybe we have it now?!
            imol_for_map = Imol_Refinement_Map();
	    if (imol_for_map > -1) {
	       // std::cout << "DEBUG:: starting" << std::endl;
	       mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	       // double max_grid_factor = coot::util::max_gridding(molecules[imol_for_map].xmap_list[0]);
	       // std::cout << "DEBUG:: max_grid_factor: " << max_grid_factor << std::endl;
	       // double sq_max_grid_fac = max_grid_factor * max_grid_factor;
	       // 1.8A data -> max_grid_factor = 0.6
	       // 3.0A data -> max_grid_factor = 1.0
	       // squared max_grid_factor: 0.36: 4.0 is good
	       // squared max_grid_factor: 1.0 : 4.0 is 4 times too much
	       // so we want 1/squared(max_grid_factor) not 4.0
	       int n_models = mol->GetNumberOfModels();
	       int max_chain_length = coot::util::max_min_max_residue_range(mol);
	       if (max_chain_length <= 0) {
		  std::cout << "WARNING:: Funny coords - no graphs\n";
	       } else {
		  for (int imodel = 1; imodel <= n_models; imodel++) { 
		     mmdb::Model *model_p = mol->GetModel(imodel);
                     if (! model_p) continue;
		     int n_chains = model_p->GetNumberOfChains();
		     coot::geometry_graphs *graphs =
			new coot::geometry_graphs(coot::GEOMETRY_GRAPH_DENSITY_FIT,
						  imol,
						  graphics_info_t::molecules[imol].name_for_display_manager(), 
						  n_chains, max_chain_length);
		     
		     // residue_density_fit_graph[imol] = graphs->dialog();
		     set_validation_graph(imol, coot::GEOMETRY_GRAPH_DENSITY_FIT, graphs->dialog());
		     
		     for (int ich=0; ich<n_chains; ich++) {
                        mmdb::Chain *chain_p = model_p->GetChain(ich);
			const char *chain_id = chain_p->GetChainID();
			std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
			if (m.first) {
			   int offset = m.second - 1;
			   int selHnd = mol->NewSelection(); // d
			   mmdb::PResidue *SelResidues = NULL;
			   int nSelResidues;
			   mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
				       chain_id,
				       mmdb::ANY_RES, "*",
				       mmdb::ANY_RES, "*",
				       "*",  // residue name
				       "*",  // Residue must contain this atom name?
				       "*",  // Residue must contain this Element?
				       "*",  // altLocs
				       mmdb::SKEY_NEW // selection key
				       );
			   mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

			   std::vector<coot::geometry_graph_block_info_generic> v = 
			      graphics_info_t::density_fit_from_residues(SelResidues, nSelResidues,
									 imol, 
									 imol_for_map);
			   
			   if (nSelResidues > 0) {
			      int max_resno = -9999;
			      for (int ires=0; ires<nSelResidues; ires++)
				 if (SelResidues[ires]->GetSeqNum() > max_resno)
				    max_resno = SelResidues[ires]->GetSeqNum();

			      graphs->render_to_canvas(v, ich, std::string(chain_id),
						       max_resno, m.second, offset);
			   }
			   mol->DeleteSelection(selHnd);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL   
}

// Use this to update the molecule only.  (Because we call a graph
// function that doesn't setup the canvas and draw the ticks (it only
// updates the blocks that it has been given).
//
// We pass imol_moving_atoms because we will be updating a graph and
// we want to know which graph to update.
// 
#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::density_fit_from_mol(const atom_selection_container_t &asc,
				      int imol_moving_atoms,
				      int imol_map) {


   std::vector<coot::geometry_graph_block_info_generic> drv;
   std::string altconf("");  // use this (e.g. "A") or "".

   if (! asc.mol)
      return drv;
   
   if (imol_map < n_molecules() && graphics_info_t::molecules[imol_map].has_xmap()) { 
      int n_models = asc.mol->GetNumberOfModels();
   
      if (n_models > 0) { 
      
	 for (int imod=1; imod<=n_models; imod++) { 
	 
	    mmdb::Model *model_p = asc.mol->GetModel(imod);
	    mmdb::Chain *chain_p;
	    int nchains = model_p->GetNumberOfChains();
	 
	    for (int ichain=0; ichain<nchains; ichain++) {
	    
	       chain_p = model_p->GetChain(ichain);
               const char *chain_id = chain_p->GetChainID();

	       // Maybe we could do a chain->GetResidueTable() here
	       // instead of a selection.
	    
	       int selHnd = asc.mol->NewSelection();
	       int nSelResidues;
	       mmdb::PResidue *SelResidues = NULL;
	    
	       asc.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
			       chain_id,
			       mmdb::ANY_RES, "*",
			       mmdb::ANY_RES, "*",
			       "*",  // residue name
			       "*",  // Residue must contain this atom name?
			       "*",  // Residue must contain this Element?
			       "*",  // altLocs
			       mmdb::SKEY_NEW // selection key
			       );
	       asc.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

	       std::vector<coot::geometry_graph_block_info_generic> v = 
		  density_fit_from_residues(SelResidues, nSelResidues, imol_moving_atoms, imol_map);

	       for (unsigned int i=0; i<v.size(); i++)
		  drv.push_back(v[i]);

	       // the graph update is done in the function that calls this one.

	       asc.mol->DeleteSelection(selHnd);
	    
	    }
	 }
      }
   }
   return drv;
}
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL   



// To be called for each chain in the molecule (or atom selection).
// 
#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::density_fit_from_residues(mmdb::PResidue *SelResidues, int nSelResidues,
					   int imol,
					   int imol_for_map) const {

   std::vector<coot::geometry_graph_block_info_generic> v;
   if (nSelResidues > 0) {
      int max_resno = -9999;
      double max_grid_factor = coot::util::max_gridding(molecules[imol_for_map].xmap);

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
	 int this_resno = residue_p->GetSeqNum();
	 if (this_resno > max_resno)
	    max_resno = this_resno;

	 mmdb::PAtom *residue_atoms=0;
	 int n_residue_atoms;

	 SelResidues[ir]->GetAtomTable(residue_atoms, n_residue_atoms);

	 double residue_density_score =
	    coot::util::map_score(residue_atoms,
				  n_residue_atoms,
				  molecules[imol_for_map].xmap, 1);
	 double occ_sum = coot::util::occupancy_sum(residue_atoms, n_residue_atoms);
	 if (occ_sum > 0) {
            float distortion_max_abs = 132.0;
            float distortion_max = distortion_max_abs;
	    residue_density_score /= occ_sum;
	    std::string str = int_to_string(this_resno);
	    str += residue_p->GetChainID();
	    str += " ";
	    str += residue_p->name;
	    str += " ";
	    str += float_to_string(residue_density_score);

	    if (residue_density_score < 0.0001)
	       residue_density_score = 0.0001;

	    // std::cout << "DEBUG::          max_grid_factor " << max_grid_factor
	    // << " score " << residue_density_score << std::endl;
	    double sf = residue_density_fit_scale_factor * 1.25;
	    // high resolution maps have high grid factors (say 0.5) and high
	    // residue_density_ scores (say 2.0)
	    // double distortion = sf/(pow(max_grid_factor,3) * residue_density_score);
	    double distortion = sf/(pow(max_grid_factor,4) * residue_density_score);  // Hmm.

            std::cout << " residue " << coot::residue_spec_t(residue_p) << " distortion: " << distortion
                      << " vs distortion_max " << distortion_max
                      << " and distortion_max_abs " << distortion_max_abs
                      << std::endl;

	    // distortion *= distortion; // non-linear, provides distinction.

	    if (distortion > distortion_max)
	       distortion = distortion_max;
	    // use intelligent atom name here
	    std::string chain_id = residue_p->GetChainID();
            std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
            coot::atom_spec_t atom_spec(chain_id, this_resno, "", atom_name, "");
            // std::cout << "creating block with distortion " << distortion << std::endl;
	    v.push_back(coot::geometry_graph_block_info_generic(imol, this_resno, atom_spec, distortion, str));
	 }
      }
      //      graphs->render_to_canvas(v, ich, std::string(chain_id),
      //                               max_resno, m.second, offset);
   }
   return v;
}
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::ncs_diffs(int imol, const coot::ncs_chain_difference_t &d) {
   std::vector<coot::geometry_graph_block_info_generic> v;

   std::cout << "peer chain id in ncs_diffs: " << d.peer_chain_id << std::endl;
   for (unsigned int ires=0; ires<d.residue_info.size(); ires++) {
      if (d.residue_info[ires].filled) { 
	 float distance = d.residue_info[ires].mean_diff;
	 mmdb::Atom *at = molecules[imol].atom_intelligent(d.peer_chain_id,
						      d.residue_info[ires].resno,
						      d.residue_info[ires].inscode);
	 std::string atom_name = " CA ";
	 std::string altconf = "";
	 if (at) {
	    atom_name = at->name;
	    altconf = at->altLoc;
	 }
	 
	 coot::atom_spec_t as(d.peer_chain_id,
			      d.residue_info[ires].resno,
			      d.residue_info[ires].inscode,
			      atom_name, altconf);
	 std::string str = coot::util::int_to_string(d.residue_info[ires].resno);
	 str += d.peer_chain_id;
	 str += " mean d = ";
	 str += coot::util::float_to_string(distance);
	 str += "A";
	 coot::geometry_graph_block_info_generic block(imol, d.residue_info[ires].resno, as,
						       30.0*distance, str);
	 v.push_back(block);
      }
   }
   return v;
}
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL

#ifdef HAVE_GSL
#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::ncs_diffs_from_mol(int imol) {

   int imodel = 1;
   std::vector<coot::geometry_graph_block_info_generic> drv;
   std::string altconf("");  // use this (e.g. "A") or "".

   if (is_valid_model_molecule(imol)) { 
      std::pair<bool, std::string> master_info = graphics_info_t::molecules[imol].first_ncs_master_chain_id();
      if (master_info.first) {
	 std::string master = master_info.second; // a target_chain_id
	 float w = 1.0; // main chain weight
	 coot::ncs_differences_t diff = graphics_info_t::molecules[imol].ncs_chain_differences(master, w);
   
	 mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	 mmdb::Model *model_p = mol->GetModel(imodel);
	 int n_chains = diff.diffs.size();
	 int max_chain_length = coot::util::max_min_max_residue_range(mol);
	 coot::geometry_graphs *graphs =
	    new coot::geometry_graphs(coot::GEOMETRY_GRAPH_NCS_DIFFS, imol,
				      graphics_info_t::molecules[imol].name_for_display_manager(), 
				      n_chains, max_chain_length);

	 // ncs_diffs_graph[imol] = graphs->dialog();
	 set_validation_graph(imol, coot::GEOMETRY_GRAPH_NCS_DIFFS, graphs->dialog());
	 for (unsigned int incs_set=0; incs_set<diff.diffs.size(); incs_set++) {

	    // do this for each chain
	    int min_resno =  99999; 
	    int max_resno = -99999; 
	    int offset = 0;

	    // diffs.diffs is a vector of chain differences (vector of ncs_chain_differences_t)
	    // A ncs_chain_differences_t contains a vector of residue difference infos.

	    for (unsigned int ires=0; ires<diff.diffs[incs_set].residue_info.size(); ires++) {
	       if (0) 
		  std::cout << "DEBUG:: resno for diffs: "
			    << diff.diffs[incs_set].residue_info[ires].resno
			    << std::endl;
	    
	       if (diff.diffs[incs_set].residue_info[ires].resno < min_resno) { 
		  min_resno = diff.diffs[incs_set].residue_info[ires].resno;
	       }
	       if (diff.diffs[incs_set].residue_info[ires].resno > max_resno) { 
		  max_resno = diff.diffs[incs_set].residue_info[ires].resno;
	       }
	    }
	    offset = min_resno - 1; 
	    //       std::cout << "max_resno, min_resno " << max_resno << " "
	    // 		<< min_resno << std::endl;
      
	    std::vector<coot::geometry_graph_block_info_generic> v = 
	       graphics_info_t::ncs_diffs(imol, diff.diffs[incs_set]);
	    graphs->render_to_canvas(v, incs_set, diff.diffs[incs_set].peer_chain_id,
				     max_resno, min_resno, offset);
	 }
      }
   } 
   return drv;
}
#endif // defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
#endif // HAVE_GSL
