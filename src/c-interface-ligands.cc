/* src/c-interface-ligands.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2008, 2009, 2010, 2011, 2012 by The University of Oxford
 * Copyright 2012, 2013, 2014, 2015 by Medical Research Council
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,  02110-1301, USA
 */

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"


#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#if defined _MSC_VER
#include <windows.h>
#endif

#include "globjects.h" //includes gtk/gtk.h

#include <vector>
#include <string>

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"

#include "coot-utils/read-sm-cif.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/peak-search.hh"
#include "coot-utils/coot-h-bonds.hh"

#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-generic-objects.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"

// #include "coot-compare-residues.hh"

#include "ligand/wligand.hh"

#include "c-interface-ligands-swig.hh"

// #include "guile-fixups.h"

#ifdef USE_GUILE
#include <cstdio> /* for std::FILE in gmp.h for libguile.h */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#include <libguile.h>
#pragma GCC diagnostic pop
#endif

#include "utils/logging.hh"
extern logging logger;

#include "get-monomer.hh"


#ifdef HAVE_GOOCANVAS
#include <goocanvas.h>
#include "lbg/wmolecule.hh"
#endif

#include "c-interface-bonds.hh"

#include "widget-headers.hh"
#include "widget-from-builder.hh"

PyObject *go_to_ligand_py() {

   clipper::Coord_orth new_pos = go_to_ligand_inner();
   PyObject *r = PyList_New(3);
   for (std::size_t i=0; i<3; i++)
      PyList_SetItem(r, i, PyFloat_FromDouble(new_pos[i]));
   return r;
}

clipper::Coord_orth
go_to_ligand_inner() {

   // note to self: this function *starts* the view interpolation and finishes quickly

   clipper::Coord_orth new_rotation_centre;

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      if (is_valid_model_molecule(imol)) {

	 graphics_info_t g;
	 clipper::Coord_orth rc(graphics_info_t::RotationCentre_x(),
				graphics_info_t::RotationCentre_y(),
				graphics_info_t::RotationCentre_z());
	 coot::new_centre_info_t new_centre =
	    graphics_info_t::molecules[pp.second.first].new_ligand_centre(rc, graphics_info_t::go_to_ligand_n_atoms_limit);
         if (false)
            std::cout << "debug:: go_to_ligand_inner(): new_centre: "
                      << new_centre.position.format() << " info \""
                      << new_centre.info_string << "\" "
                      << new_centre.residue_spec << " type "
                      << new_centre.type << " "
                      << std::endl;
	 new_rotation_centre = new_centre.position;
	 if (new_centre.type == coot::NORMAL_CASE) {
	    // g.setRotationCentre(new_centre.position); // the target position is (should be) reached by the tick animation
                                                         // in coot::view_info_t::interpolate()
	    g.perpendicular_ligand_view(imol, new_centre.residue_spec);
            std::cout << "::::::::::::::::::::::::::: go_to_ligand_inner() C " << std::endl;
            std::string name = g.molecules[imol].get_residue_name(new_centre.residue_spec);

            // 20220303-PE  No, we are still at the starting point when we get here.
            // view_info_t::interpolate() sets up it's own animation function (glarea_tick_function() is not used)
	    // g.update_things_on_move_and_redraw(); // now not done in perpendicular_ligand_view()

	    std::string s = "INFO:: Centred on residue ";
	    s += new_centre.residue_spec.chain_id;
	    s+= " ";
	    s += coot::util::int_to_string(new_centre.residue_spec.res_no);
	    s += new_centre.residue_spec.ins_code;
	    s+= " ";
	    s+= name;
	    s+= " ";
	    s += " in molecule #";
	    s += coot::util::int_to_string(pp.second.first);
	    s += ".";
	    add_status_bar_text(s.c_str());
            std::cout << "INFO:: status bar text: " << s << std::endl;
	 } else {
	    if (new_centre.type == coot::NO_LIGANDS) {
	       std::string s = "No ligand (hetgroup) found in this molecule (#";
	       s += coot::util::int_to_string(pp.second.first);
	       s += ").";
	       add_status_bar_text(s.c_str());
	    }
	    if (new_centre.type == coot::SINGLE_LIGAND_NO_MOVEMENT) {
	       g.perpendicular_ligand_view(pp.second.first, new_centre.residue_spec);
	       std::string s = "This ligand (";
	       s += coot::util::int_to_string(new_centre.residue_spec.res_no);
	       s += new_centre.residue_spec.ins_code;
	       s += new_centre.residue_spec.chain_id;
	       s += ") ";
	       s += " is the only ligand in this molecule (#";
	       s += coot::util::int_to_string(pp.second.first);
	       s += ").";
	       add_status_bar_text(s.c_str());
	    }
	 }
      }
   }
   return new_rotation_centre;
}


/*! \brief centre on the ligand of the "active molecule", if we are
  already there, centre on the next hetgroup (etc) */
void go_to_ligand() {

   go_to_ligand_inner();
}

void set_go_to_ligand_n_atoms_limit(int n_atoms_min) {
   graphics_info_t::go_to_ligand_n_atoms_limit = n_atoms_min;
}




/*  ----------------------------------------------------------------------- */
/*                  ligand overlay                                          */
/*  ----------------------------------------------------------------------- */
/*! \brief "Template"-based matching.  Overlap the first residue in
  imol_ligand onto the residue specified by the reference parameters.
  Use graph matching, not atom names.  */

#ifdef USE_GUILE
SCM
overlap_ligands(int imol_ligand, int imol_ref, const char *chain_id_ref,
		int resno_ref) {

   SCM scm_status = SCM_BOOL_F;
   coot::graph_match_info_t rtop_info =
      overlap_ligands_internal(imol_ligand, imol_ref, chain_id_ref, resno_ref, 1);

   if (rtop_info.success) {
      SCM match_info = scm_cons(scm_from_int(rtop_info.n_match), SCM_EOL);
      match_info = scm_cons(scm_from_double(rtop_info.dist_score), match_info);
      SCM s = scm_cons(match_info, SCM_EOL);
      scm_status = scm_cons(rtop_to_scm(rtop_info.rtop), s);
   }
   return scm_status;
}
#endif

// should be in utils?
void
match_ligand_torsions(int imol_ligand, int imol_ref, const char *chain_id_ref, int resno_ref) {

   // 20110608 I don't think that this function should make
   // restraints.  Let's get restraints outside and test for their
   // existance before we run this function.

   if (!is_valid_model_molecule(imol_ligand)) {
      std::cout << "WARNING molecule number " << imol_ligand << " is not a valid model molecule"
		<< std::endl;
   } else {
      if (! is_valid_model_molecule(imol_ref)) {
	 std::cout << "WARNING molecule number " << imol_ref << " is not a valid model molecule"
		   << std::endl;
      } else {

	 if (is_valid_model_molecule(imol_ref)) {
	    graphics_info_t g;
	    mmdb::Manager *mol_ref = g.molecules[imol_ref].atom_sel.mol;
	    mmdb::Residue *res_ref = coot::util::get_residue(chain_id_ref,
							resno_ref, "",
							mol_ref);
	    if (res_ref) {

	       std::string res_name_ref_res(res_ref->GetResName());
	       std::pair<bool, coot::dictionary_residue_restraints_t> restraints_info =
		  g.Geom_p()->get_monomer_restraints(res_name_ref_res, imol_ref);
	       if (restraints_info.first) {
		  std::vector <coot::dict_torsion_restraint_t> tr =
		     g.Geom_p()->get_monomer_torsions_from_geometry(res_name_ref_res, 0);

		  if (tr.size() == 0) {
		     std::cout << "WARNING:: No torsion restraints from PRODRG"
			       << std::endl;
		  } else {
		     // normal case
		     int n_rotated = g.molecules[imol_ligand].match_torsions(res_ref, tr, *g.Geom_p());
		     std::cout << "INFO:: rotated " << n_rotated << " torsions in matching torsions"
			       << std::endl;
		  }
	       }
	       graphics_draw();
	    }
	 }
      }
   }
}


#ifdef USE_PYTHON
PyObject *overlap_ligands_py(int imol_ligand, int imol_ref, const char *chain_id_ref,
		int resno_ref) {

   PyObject *python_status = Py_False;
   coot::graph_match_info_t rtop_info =
      overlap_ligands_internal(imol_ligand, imol_ref, chain_id_ref, resno_ref, 1);

   if (rtop_info.success) {
      PyObject *match_info = PyList_New(2);
      PyList_SetItem(match_info, 0, PyFloat_FromDouble(rtop_info.dist_score));
      PyList_SetItem(match_info, 1, PyLong_FromLong(rtop_info.n_match));
      python_status = PyList_New(2);
      PyList_SetItem(python_status, 0, rtop_to_python(rtop_info.rtop));
      PyList_SetItem(python_status, 1, match_info);
   }
   if (PyBool_Check(python_status)) {
     Py_INCREF(python_status);
   }
   return python_status;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
SCM
analyse_ligand_differences(int imol_ligand, int imol_ref, const char *chain_id_ref,
			   int resno_ref) {

   SCM scm_status = SCM_BOOL_F;
   coot::graph_match_info_t rtop_info =
      overlap_ligands_internal(imol_ligand, imol_ref, chain_id_ref, resno_ref, 0);

   std::cout << "analyse_ligand_differences: success       " << rtop_info.success << std::endl;
   std::cout << "analyse_ligand_differences: n_match       " << rtop_info.n_match << std::endl;
   std::cout << "analyse_ligand_differences: dist_score    " << rtop_info.dist_score << std::endl;
   std::cout << "analyse_ligand_differences: atoms matched " << rtop_info.matching_atom_names.size() << std::endl;
   std::cout << "analyse_ligand_differences: rtop: \n" << rtop_info.rtop.format() << std::endl;

   if (rtop_info.success) {
      SCM match_info = scm_cons(scm_from_int(rtop_info.n_match), SCM_EOL);
      match_info = scm_cons(scm_from_double(rtop_info.dist_score), match_info);
      SCM s = scm_cons(match_info, SCM_EOL);
      scm_status = scm_cons(rtop_to_scm(rtop_info.rtop), s);
   }
   return scm_status;
}
#endif

#ifdef USE_PYTHON
PyObject *analyse_ligand_differences_py(int imol_ligand, int imol_ref, const char *chain_id_ref,
			   int resno_ref) {

   PyObject *python_status;
   python_status = Py_False;
   coot::graph_match_info_t rtop_info =
      overlap_ligands_internal(imol_ligand, imol_ref, chain_id_ref, resno_ref, 0);

   std::cout << "analyse_ligand_differences: success       " << rtop_info.success << std::endl;
   std::cout << "analyse_ligand_differences: n_match       " << rtop_info.n_match << std::endl;
   std::cout << "analyse_ligand_differences: dist_score    " << rtop_info.dist_score << std::endl;
   std::cout << "analyse_ligand_differences: atoms matched " << rtop_info.matching_atom_names.size() << std::endl;
   std::cout << "analyse_ligand_differences: rtop: \n" << rtop_info.rtop.format() << std::endl;

   if (rtop_info.success) {
      PyObject *match_info;
      match_info = PyList_New(2);
      PyList_SetItem(match_info, 0, PyFloat_FromDouble(rtop_info.dist_score));
      PyList_SetItem(match_info, 1, PyLong_FromLong(rtop_info.n_match));
      python_status = PyList_New(2);
      PyList_SetItem(python_status, 0, rtop_to_python(rtop_info.rtop));
      PyList_SetItem(python_status, 1, match_info);
   }
   if (PyBool_Check(python_status)) {
     Py_INCREF(python_status);
   }
   return python_status;
}
#endif //USE_PYTHON



#ifdef USE_GUILE
// For testing that pyrogen generates consistent atom types.  Not
// useful for general public.
//
SCM
compare_ligand_atom_types_scm(int imol_ligand, int imol_ref, const char *chain_id_ref,
			      int resno_ref) {

   SCM scm_status = SCM_BOOL_F;
   bool match_hydrogens_also = false;
   bool apply_rtop_flag = true;

   if (! is_valid_model_molecule(imol_ligand)) {
      std::cout << "WARNING:: not a valid model molecule (ligand) " << imol_ligand << std::endl;
   } else {
      if (! is_valid_model_molecule(imol_ref)) {
	 std::cout << "WARNING:: not a valid model molecule (ref) " << imol_ligand << std::endl;
      } else {

	 graphics_info_t g;
	 mmdb::Residue *res_ref = g.molecules[imol_ref].get_residue(chain_id_ref, resno_ref, "");
	 mmdb::Residue *res_mov = g.molecules[imol_ligand].get_first_residue();

	 if (! res_ref) {
	    std::cout << "WARNING failed to find reference residue" << std::endl;
	 } else {
	    if (! res_mov) {
	       std::cout << "WARNING failed to find moving residue" << std::endl;
	    } else {

	       coot::graph_match_info_t rtop_info =
		  coot::graph_match(res_mov, res_ref, apply_rtop_flag, match_hydrogens_also);

	       // work atom-names then ref atom-names
	       std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > matches = rtop_info.matching_atom_names;

	       std::cout << "found " << matches.size() << " graph matched atoms" << std::endl;
	       std::string rn1 = res_mov->GetResName();
	       std::string rn2 = res_ref->GetResName();

	       if (matches.size()) {
		  int n_fail = 0;

		  for (unsigned int i=0; i<matches.size(); i++) {
		     std::string energy_type_1 = g.Geom_p()->get_type_energy(matches[i].first.first,  rn1, imol_ligand);
		     std::string energy_type_2 = g.Geom_p()->get_type_energy(matches[i].second.first, rn2, imol_ref);

		     std::string extra_space = "";
		     if (i<10) extra_space = " ";
		     std::cout << "   " << extra_space << i << " names: \""
			       << matches[i].second.first << "\" \""
			       << matches[i].first.first << "\" ->  "
			       << "\"" << energy_type_2 << "\"  and  \"" << energy_type_1 << "\"";
		     if (energy_type_1 != energy_type_2) {
			std::cout << "   #### fail " << std::endl;
			n_fail++;
		     } else {
			std::cout << std::endl;
		     }
		  }
		  scm_status = scm_from_int(n_fail);
	       }
	    }
	 }
      }
   }
   return scm_status;

}
#endif

#ifdef USE_PYTHON
// For testing that pyrogen generates consistent atom types.  Not
// useful for general public.
//
PyObject *
compare_ligand_atom_types_py(int imol_ligand, int imol_ref, const char *chain_id_ref,
                             int resno_ref) {

   PyObject *py_status = Py_False;
   bool match_hydrogens_also = false;
   bool apply_rtop_flag = true;

   if (! is_valid_model_molecule(imol_ligand)) {
      std::cout << "WARNING:: not a valid model molecule (ligand) " << imol_ligand << std::endl;
   } else {
      if (! is_valid_model_molecule(imol_ref)) {
	 std::cout << "WARNING:: not a valid model molecule (ref) " << imol_ligand << std::endl;
      } else {

	 graphics_info_t g;
	 mmdb::Residue *res_ref = g.molecules[imol_ref].get_residue(chain_id_ref, resno_ref, "");
	 mmdb::Residue *res_mov = g.molecules[imol_ligand].get_first_residue();

     if (! res_ref) {
        std::cout << "WARNING failed to find reference residue" << std::endl;
     } else {
        if (! res_mov) {
           std::cout << "WARNING failed to find moving residue" << std::endl;
        } else {

           coot::graph_match_info_t rtop_info =
              coot::graph_match(res_mov, res_ref, apply_rtop_flag, match_hydrogens_also);

	       // work atom-names then ref atom-names
	       std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > matches = rtop_info.matching_atom_names;

	       std::cout << "found " << matches.size() << " graph matched atoms" << std::endl;
	       std::string rn1 = res_mov->GetResName();
	       std::string rn2 = res_ref->GetResName();

	       if (matches.size()) {
              int n_fail = 0;

              for (unsigned int i=0; i<matches.size(); i++) {
                 std::string energy_type_1 = g.Geom_p()->get_type_energy(matches[i].first.first,  rn1, imol_ligand);
                 std::string energy_type_2 = g.Geom_p()->get_type_energy(matches[i].second.first, rn2, imol_ref);

                 std::string extra_space = "";
                 if (i<10) extra_space = " ";
                 std::cout << "   " << extra_space << i << " names: \""
                           << matches[i].second.first << "\" \""
                           << matches[i].first.first << "\" ->  "
                           << "\"" << energy_type_2 << "\"  and  \"" << energy_type_1 << "\"";
                 if (energy_type_1 != energy_type_2) {
                    std::cout << "   #### fail " << std::endl;
                    n_fail++;
                 } else {
                    std::cout << std::endl;
                 }
              }
              py_status = PyLong_FromLong(n_fail);
	       }
	    }
	 }
      }
   }
   if (PyBool_Check(py_status)) {
      Py_XINCREF(py_status);
   }
   return py_status;

}
#endif // USE_PYTHON

coot::graph_match_info_t
overlap_ligands_internal(int imol_ligand, int imol_ref, const char *chain_id_ref,
			 int resno_ref, bool apply_rtop_flag) {

   coot::graph_match_info_t graph_info;

   mmdb::Residue *residue_moving = 0;
   mmdb::Residue *residue_reference = 0;

   // running best ligands:
   mmdb::Residue *best_residue_moving = NULL;
   double best_score = 99999999.9; // low score good.
   clipper::RTop_orth best_rtop;


   if (! is_valid_model_molecule(imol_ligand))
      return graph_info;

   if (! is_valid_model_molecule(imol_ref))
      return graph_info;

   mmdb::Manager *mol_moving = graphics_info_t::molecules[imol_ligand].atom_sel.mol;
   mmdb::Manager *mol_ref    = graphics_info_t::molecules[imol_ref].atom_sel.mol;

   for (int imod=1; imod<=mol_moving->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol_moving->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::PResidue residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p) {
	       int n_atoms = residue_p->GetNumberOfAtoms();
	       if (n_atoms > 0) {
		  residue_moving = residue_p;
		  break;
	       }
	    }
	 }
	 if (residue_moving)
	    break;
      }

      if (! residue_moving) {
	 std::cout << "Oops.  Failed to find moving residue" << std::endl;
      } else {
	 int imodel_ref = 1;
	 mmdb::Model *model_ref_p = mol_ref->GetModel(imodel_ref);
	 mmdb::Chain *chain_p;
	 int nchains = model_ref_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_ref_p->GetChain(ichain);
	    if (std::string(chain_p->GetChainID()) == std::string(chain_id_ref)) {
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::PResidue residue_p;
	       for (int ires=0; ires<nres; ires++) {
		  residue_p = chain_p->GetResidue(ires);
		  if (residue_p) {
		     int seqnum = residue_p->GetSeqNum();
		     if (seqnum == resno_ref) {
			residue_reference = residue_p;
			break;
		     }
		  }
	       }
	       if (residue_reference)
		  break;
	    }
	 }

	 if (!residue_reference) {
	    std::cout << "Oops.  Failed to find reference residue" << std::endl;
	 } else {
	    bool match_hydrogens_also = 0;
	    coot::graph_match_info_t rtop_info =
	       coot::graph_match(residue_moving, residue_reference, apply_rtop_flag, match_hydrogens_also);
	    if (rtop_info.success) {
	       if (rtop_info.dist_score < best_score) { // low score good.
		  best_score = rtop_info.dist_score;
		  best_residue_moving = residue_moving;
		  best_rtop = rtop_info.rtop;
		  graph_info = rtop_info;
	       }
	    } else {
	       // std::cout << "Oops.  Match failed somehow" << std::endl;
	    }
	 }
      }
   }
   if (apply_rtop_flag) {
      if (best_residue_moving) {
	 // move just the best ligand:
	 graphics_info_t::molecules[imol_ligand].transform_by(best_rtop, best_residue_moving);
	 // delete everything except the best ligand:
	 graphics_info_t::molecules[imol_ligand].delete_all_except_res(best_residue_moving);
	 graphics_draw();
      }
   }
   return graph_info;
}

// Read cif_dict_in, match the atoms there-in to those of the dictionary of reference_comp_id.
// Write a new dictionary to cif_dict_out.
//
// Return 1 if was successful in doing the atom matching and writing the cif file.
int
match_residue_and_dictionary(int imol, std::string chain_id, int res_no, std::string ins_code,
			     std::string cif_dict_in,
			     std::string cif_dict_out,
			     std::string cif_dict_comp_id,
			     std::string reference_comp_id,
			     std::string output_comp_id,
			     std::string output_compound_name) {

   int result = 0;

   if (coot::file_exists(cif_dict_in)) {
      coot::protein_geometry geom_local;
      geom_local.try_dynamic_add(reference_comp_id, 0);
      std::pair<short int, coot::dictionary_residue_restraints_t> rp_1 =
	 geom_local.get_monomer_restraints(reference_comp_id, imol);
      if (rp_1.first) {
	 coot::protein_geometry geom_matcher;
	 coot::read_refmac_mon_lib_info_t rmit = geom_matcher.init_refmac_mon_lib(cif_dict_in, 0);
	 if (rmit.n_bonds == 0) {
	    std::cout << "No bonds from " << cif_dict_in << std::endl;
	 } else {
	    std::pair<short int, coot::dictionary_residue_restraints_t> rp_2 =
	       geom_matcher.get_monomer_restraints(cif_dict_comp_id, imol);
	    if (rp_2.first) {
	       mmdb::Residue *residue_p = NULL; // for the moment

	       std::cout << "DEBUG:: ------ about to match "
			 << rp_2.second.residue_info.comp_id << " to "
			 << rp_1.second.residue_info.comp_id << " names" << std::endl;
	       coot::dictionary_match_info_t dmi =
		  rp_2.second.match_to_reference(rp_1.second, residue_p, output_comp_id,
						 output_compound_name);
	       if (dmi.n_matches > 0) {
		  dmi.dict.residue_info.comp_id = output_comp_id;
		  dmi.dict.residue_info.name = output_compound_name;
		  dmi.dict.write_cif(cif_dict_out);
	       } else {
		  std::cout << "WARNING:: not similar enough, n_matches = "
			    << dmi.n_matches << std::endl;
	       }
	    } else {
	       std::cout << "WARNING:: no monomer restraints for " << cif_dict_comp_id
			 << " in " << cif_dict_in << std::endl;
	    }
	 }
      } else {
	 std::cout << "WARNING:: No restraints for " << reference_comp_id << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << cif_dict_in << " file not found" << std::endl;
   }
   return result;
}


// This is the GUI interface: User sits over on top of the residue
// they want to change and provides an output dictionary cif file
// name, a reference-comp-id and a compd-id for the new residue (type).
int
match_this_residue_and_dictionary(int imol, std::string chain_id, int res_no, std::string ins_code,
				  std::string cif_dict_out,
				  std::string reference_comp_id,
				  std::string output_comp_id) {

   int result = 0;

   if (false)
      std::cout << "here in match_this_residue_and_dictionary() "
		<< " cif_dict_out " << cif_dict_out
		<< " reference_comp_id " << reference_comp_id
		<< " output_comp_id " << output_comp_id
		<< std::endl;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Residue *this_residue = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (this_residue) {
	 std::string this_residue_type = this_residue->GetResName();
	 std::pair<short int, coot::dictionary_residue_restraints_t> dict_1 =
	    g.Geom_p()->get_monomer_restraints(this_residue_type, imol);
	 if (dict_1.first) {

	    g.Geom_p()->try_dynamic_add(reference_comp_id, g.cif_dictionary_read_number++);
	    int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
	    std::pair<short int, coot::dictionary_residue_restraints_t> dict_2 =
	       g.Geom_p()->get_monomer_restraints(reference_comp_id, imol_enc);

	    if (dict_2.first) {

	       coot::dictionary_match_info_t dmi =
		  dict_1.second.match_to_reference(dict_2.second, this_residue,
						   output_comp_id, output_comp_id); // placeholder for name

	       if (dmi.n_matches > 0) {
		  dmi.dict.residue_info.comp_id = output_comp_id;
		  dmi.dict.residue_info.name =  ".";
		  dmi.dict.write_cif(cif_dict_out);
	       } else {
		  std::cout << "No matches " << std::endl;
	       }

	    } else {
	       std::cout << "WARNING:: match_this_residue_and_dictionary, no dictionary "
			 << " for reference type " << output_comp_id << std::endl;
	    }
	 } else {
	    std::cout << "WARNING:: match_this_residue_and_dictionary, no dictionary for type "
	       << this_residue_type << std::endl;
	 }
      } else {
	 std::cout << "WARNING:: match_this_residue_and_dictionary, no such residue " <<
	    coot::residue_spec_t(chain_id, res_no, ins_code) << std::endl;
      }
   } else {
      std::cout << "WARNING:: match_this_residue_and_dictionary, no such molecule " << imol
		<< std::endl;
   }
   return result;
}








void ligand_expert() {  /* sets the flag to have expert option ligand entries in the GUI */

   graphics_info_t::ligand_expert_flag = 1;
}

/*! \brief set the protein molecule for ligand searching */
void set_ligand_search_protein_molecule(int imol) {
   graphics_info_t::set_ligand_protein_mol(imol);

}

/*! \brief set the map molecule for ligand searching */
void set_ligand_search_map_molecule(int imol_map) {
   graphics_info_t::set_ligand_map_mol(imol_map);
}

/*! \brief add a ligand molecule to the list of ligands to search for
  in ligand searching */
void add_ligand_search_wiggly_ligand_molecule(int imol_ligand) {
   graphics_info_t::find_ligand_add_flexible_ligand(imol_ligand);
}

/*! \brief add a ligand molecule to the list of ligands to search for
  in ligand searching */
void add_ligand_search_ligand_molecule(int imol_ligand) {
   if (is_valid_model_molecule(imol_ligand))
      graphics_info_t::find_ligand_add_rigid_ligand(imol_ligand);
}

void add_ligand_clear_ligands() {

   graphics_info_t g;
   g.find_ligand_clear_ligand_mols();
}

#include "c-interface-ligand-search.hh"

// Called from callbacks.c
void execute_ligand_search() {

   // 20161031: If you come to look at this in the future, scrap all ligand search interface
   // and re-write from scratch.
   //
   // It should have 2 paths: one with a info/progress-bar dialog and the (traditional)
   // silent one (say, that is used in no-graphics mode).

   if (graphics_info_t::use_graphics_interface_flag) {
      ligand_wiggly_ligand_data_t lwld = ligand_search_install_wiggly_ligands();
      // if we have wiggly/conformer ligands, execute_ligand_search_internal(wlig_p)
      // is run in install_simple_wiggly_ligand_idle_fn()

      if (lwld.immediate_execute_ligand_search)
	      execute_ligand_search_internal(lwld.wlig);

   } else {
      ligand_wiggly_ligand_data_t lwld = ligand_search_install_wiggly_ligands();
      execute_ligand_search_internal(lwld.wlig);
   }
}


// execute_find_ligands_real, you might say
//
#ifdef USE_GUILE
SCM execute_ligand_search_scm() {

   ligand_wiggly_ligand_data_t lwld = ligand_search_install_wiggly_ligands();
   std::vector<int> solutions = execute_ligand_search_internal(lwld.wlig);
   return generic_int_vector_to_list_internal(solutions);
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *execute_ligand_search_py() {

   std::vector<int> solutions;
   PyObject *r = generic_int_vector_to_list_internal_py(solutions);
   ligand_wiggly_ligand_data_t lwld = ligand_search_install_wiggly_ligands();
   if (lwld.immediate_execute_ligand_search) {
      std::vector<int> solutions = execute_ligand_search_internal(lwld.wlig);
      // now sort solutions (because they (probably) have been real-space refined now)
      r = generic_int_vector_to_list_internal_py(solutions);
   }
   return r;
}
#endif // USE_PYTHON

#include "utils/ctpl.h"

/*! \brief  Allow the user a scripting means to find ligand at the rotation centre */
void set_find_ligand_here_cluster(int state) {
   graphics_info_t g;
   g.find_ligand_here_cluster_flag = state;
}

#include "c-interface-ligands-widgets.hh"

ligand_wiggly_ligand_data_t
ligand_search_install_wiggly_ligands() {

   graphics_info_t g;
   // the bool is for flag is-flexible?
   std::vector<std::pair<int, bool> > ligands = g.find_ligand_ligand_mols();
   coot::wligand *wlig_p = new coot::wligand;
   ligand_wiggly_ligand_data_t lfwd(wlig_p);

   for(unsigned int i=0; i<ligands.size(); i++) {

      std::cout << "INFO:: ligand number " << i << " is molecule number "
		          << g.find_ligand_ligand_mols()[i].first << "  "
		          << " with wiggly flag: "
		          << g.find_ligand_ligand_mols()[i].second << std::endl;

      if (ligands[i].second) {
	      // argh (i).
	      coot::minimol::molecule mmol(g.molecules[ligands[i].first].atom_sel.mol);

//	 for(unsigned int ifrag=0; ifrag<mmol.fragments.size(); ifrag++) {
// 	    for (int ires=mmol[ifrag].min_res_no(); ires<=mmol[ifrag].max_residue_number();
//		 ires++) {
// 	       if (mmol[ifrag][ires].n_atoms() > 0) {
// 		  std::cout << "DEBUG:: in execute_ligand_search:  mmol["
// 			    << ifrag << "][" << ires << "].name :"
// 			    <<  mmol[ifrag][ires].name << ":" << std::endl;
// 	       }
//	    }
//	 }

	      // std::pair<short int, std::string> istat_pair =
	      try {
	         bool optim_geom = true;
	         bool fill_vec = false;

	         if (graphics_info_t::use_graphics_interface_flag) {
	            ligand_wiggly_ligand_data_t lfwd_local = setup_ligands_progress_bar();
               lfwd.imol_ligand = ligands[i].first;
               lfwd.progress_bar        = lfwd_local.progress_bar;
               lfwd.progress_bar_window = lfwd_local.progress_bar_window;
               lfwd.progress_bar_label  = lfwd_local.progress_bar_label;
               lfwd.immediate_execute_ligand_search = false;

               setup_ligands_progress_bar_idle(wlig_p, ligands[i].first, lfwd);

               // this GtkFunction returns a gint and takes a widget
               // gint idle = gtk_idle_add((GtkFunction) install_simple_wiggly_ligand_idle_fn,
               // progress_bar);

               // wlig.install_simple_wiggly_ligands(g.Geom_p(), mmol, ligands[i].first,
               // g.ligand_wiggly_ligand_n_samples,
               // optim_geom, fill_vec);

	         } else {
               unsigned int n_threads = coot::get_max_number_of_threads();
               ctpl::thread_pool thread_pool(n_threads);
	            wlig_p->install_simple_wiggly_ligands(g.Geom_p(), mmol, ligands[i].first,
						     g.ligand_wiggly_ligand_n_samples,
						     optim_geom, fill_vec, &thread_pool, n_threads);
	         }
	      }
	      catch (const std::runtime_error &mess) {
	         std::cout << "ERROR:: failure in flexible ligand definition.\n";
	         std::cout << mess.what() << std::endl;
	         if (graphics_info_t::use_graphics_interface_flag) {
	            GtkWidget *w = wrapped_nothing_bad_dialog(mess.what());
	            gtk_widget_set_visible(w, TRUE);
	         }
	         // return solutions;
	      }
      } else {
	      // argh (ii).
	      wlig_p->install_ligand(g.molecules[ligands[i].first].atom_sel.mol);
      }
   }
   return lfwd;
}

#include "c-interface-ligand-search.hh"

std::vector<int>
execute_ligand_search_internal(coot::wligand *wlig_p) {

   std::vector<int> solutions;

   std::cout << "in execute_ligand_search_internal() find_ligand_here_cluster_flag "
             << graphics_info_t::find_ligand_here_cluster_flag << std::endl;

   graphics_info_t g;

   if (! is_valid_model_molecule(g.find_ligand_protein_mol())) {
      std::cout << "Protein molecule for ligand search not set" << std::endl;
      std::cout << "Aborting ligand search" << std::endl;
      return solutions;
   }
   if (! is_valid_map_molecule(g.find_ligand_map_mol())) {
      std::cout << "Map molecule for ligand search not set" << std::endl;
      std::cout << "Aborting ligand search" << std::endl;
      return solutions;
   }
   if (g.find_ligand_ligand_mols().size() == 0) {
      std::cout << "No defined ligand molecules" << std::endl;
      std::cout << "Aborting ligand search" << std::endl;
      return solutions;
   }

   std::cout << "::::::: Debug in execute_ligand_search" << g.find_ligand_ligand_mols().size() << std::endl;
   for (size_t i=0; i < g.find_ligand_ligand_mols().size(); i++) {
      std::cout << i << " " << g.find_ligand_ligand_mols()[i].first << " " << g.find_ligand_ligand_mols()[i].second << std::endl;
   }


   mmdb::Manager *protein_mol =
      g.molecules[g.find_ligand_protein_mol()].atom_sel.mol;

   if (g.ligand_verbose_reporting_flag) {
      wlig_p->set_verbose_reporting();
      // debugging, output the post-conformer generation ligands wligand-*.pdb
      // (but pre-idealized).
      wlig_p->set_debug_wiggly_ligands();
   }

   std::cout << ":::::::::::: in execute_ligand_search_internal() import maps from mol "
	          << g.find_ligand_map_mol() << std::endl;
   wlig_p->import_map_from(g.molecules[g.find_ligand_map_mol()].xmap);
   std::vector<std::pair<int, bool> > ligands = g.find_ligand_ligand_mols();

   short int mask_waters_flag; // treat waters like other atoms?
   mask_waters_flag = g.find_ligand_mask_waters_flag;
   if (! g.find_ligand_here_cluster_flag) {
      int imol = graphics_info_t::create_molecule();
      if (graphics_info_t::map_mask_atom_radius > 0) {
	      // only do this if it was set by the user.
	      wlig_p->set_map_atom_mask_radius(graphics_info_t::map_mask_atom_radius);
      } else {
	      wlig_p->set_map_atom_mask_radius(2.0);  // Angstroms
      }

      std::string name("ligand masked map");
      // std::cout << "DEBUG:: calling mask_map\n";
      wlig_p->mask_map(protein_mol, mask_waters_flag); // mask by protein
      // std::cout << "DEBUG:: done mask_map\n";

      bool is_em_flag = graphics_info_t::molecules[g.find_ligand_map_mol()].is_EM_map();
      g.molecules[imol].install_new_map(wlig_p->masked_map(), wlig_p->masked_map_name(), is_em_flag);

      // This should not be be scroll map?
      if (0) {
         g.scroll_wheel_map = imol;  // change the current scrollable map to
                                     // the masked map.
      }
      wlig_p->find_clusters(g.ligand_cluster_sigma_level);  // trashes the xmap
      wlig_p->set_acceptable_fit_fraction(g.ligand_acceptable_fit_fraction);
      wlig_p->fit_ligands_to_clusters(g.find_ligand_n_top_ligands); // 10 clusters

   } else {

      // don't search the map, just use the peak/cluser at the screen
      // centre.

      wlig_p->mask_map(protein_mol, mask_waters_flag);
      // wlig.set_acceptable_fit_fraction(g.ligand_acceptable_fit_fraction);

      clipper::Coord_orth pt(g.X(), g.Y(), g.Z()); // close to 3GP peak (not in it).
      float n_sigma = g.ligand_cluster_sigma_level; // cluster points must be more than this.
      wlig_p->cluster_from_point(pt, n_sigma);
      wlig_p->fit_ligands_to_clusters(1); // just this cluster.

   }
   wlig_p->make_pseudo_atoms(); // put anisotropic atoms at the ligand sites

   // now add in the solution ligands:
   int n_clusters = wlig_p->n_clusters_final();

   int n_new_ligand = 0;
   coot::minimol::molecule m;
   for (int iclust=0; iclust<n_clusters; iclust++) {

      // frac_lim is the fraction of the score of the best solutions
      // that we should consider as solutions. 0.5 is generous, I
      // think.
      float frac_lim = g.find_ligand_score_by_correl_frac_limit; // 0.7;
      float correl_frac_lim = 0.9;
      // nino-mode
      unsigned int nlc = wlig_p->n_ligands_for_cluster(iclust, frac_lim);
      wlig_p->score_and_resort_using_correlation(iclust, nlc);

      // false is the default case
      if (g.find_ligand_multiple_solutions_per_cluster_flag == false) {
	      nlc = 1;
	      correl_frac_lim = 0.975;
      }

      if (nlc > 12) nlc = 12; // arbitrary limit of max 12 solutions per cluster
      float tolerance = 20.0;
      // limit_solutions should be run only after a post-correlation sort.
      //
      wlig_p->limit_solutions(iclust, correl_frac_lim, nlc, tolerance, true);

      for (unsigned int isol=0; isol<nlc; isol++) {

	 m = wlig_p->get_solution(isol, iclust);
	 if (! m.is_empty()) {
	    float bf = graphics_info_t::default_new_atoms_b_factor;
	    mmdb::Manager *ligand_mol = m.pcmmdbmanager();
	    coot::hetify_residues_as_needed(ligand_mol);
	    atom_selection_container_t asc = make_asc(ligand_mol);
	    int g_mol = graphics_info_t::create_molecule();
	    std::string label = "Fitted ligand #";
	    label += g.int_to_string(iclust);
	    label += "-";
	    label += g.int_to_string(isol);
	    g.molecules[g_mol].install_model(g_mol, asc, g.Geom_p(), label, 1);
	    g.molecules[g_mol].assign_hetatms();
	    if (g.find_ligand_do_real_space_refine_state()) {
	       set_imol_refinement_map(g.find_ligand_map_mol());
 	       int previous_state = refinement_immediate_replacement_state();
 	       g.refinement_immediate_replacement_flag = 1;
 	       g.refine_residue_range(g_mol, "A", "A", 1, "", 1, "", "", 0);
 	       accept_regularizement();
 	       g.refinement_immediate_replacement_flag = previous_state;
	    }
	    solutions.push_back(g_mol);
	    n_new_ligand++;
	    if (g.go_to_atom_window){
	       g.update_go_to_atom_window_on_new_mol();
	       g.update_go_to_atom_window_on_changed_mol(g_mol);
	    }
	 }
      }
   }

   if (graphics_info_t::use_graphics_interface_flag) {
      if (n_new_ligand) {

   // We need some python code here to match post-ligand-fit-gui
#if defined USE_GUILE && !defined WINDOWS_MINGW
	 safe_scheme_command("(post-ligand-fit-gui)");
#else

#ifdef USE_PYGTK
	 safe_python_command("post_ligand_fit_gui()");
#endif // USE_PYGTK
#endif // USE_GUILE

	 // GtkWidget *w = create_new_ligands_info_dialog();
	 GtkWidget *w = widget_from_builder("new_ligands_info_dialog");
	 GtkWidget *label = widget_from_builder("new_ligands_info_dialog_label");
	 std::string label_str("  Found ");
	 label_str += graphics_info_t::int_to_string(n_new_ligand);
	 if (n_new_ligand == 1)
	    label_str += " acceptable ligand  ";
	 else
	    label_str += " acceptable ligands  ";
	 gtk_label_set_text(GTK_LABEL(label), label_str.c_str());
	 gtk_widget_set_visible(w, TRUE);
      } else {
	 // GtkWidget *w = create_no_new_ligands_info_dialog();
	 GtkWidget *w = widget_from_builder("no_new_ligands_info_dialog");
	 gtk_widget_set_visible(w, TRUE);
      }
   }

   graphics_draw();
   return solutions;
}


std::vector<int> ligand_search_make_conformers_internal() {

   std::vector<int> mol_list;
   graphics_info_t g;
   coot::wligand wlig;
   if (g.ligand_verbose_reporting_flag)
      wlig.set_verbose_reporting();
   std::vector<std::pair<int, bool> > ligands = g.find_ligand_ligand_mols();
   for(unsigned int i=0; i<ligands.size(); i++) {
      if (ligands[i].second) {
	 try {
	    bool optim_geom = true;
	    bool fill_return_vec = false; // would give input molecules (not conformers)
	    coot::minimol::molecule mmol(g.molecules[ligands[i].first].atom_sel.mol);
            unsigned int n_threads = coot::get_max_number_of_threads();
            ctpl::thread_pool thread_pool(n_threads);
	    wlig.install_simple_wiggly_ligands(g.Geom_p(), mmol, ligands[i].first,
					       g.ligand_wiggly_ligand_n_samples,
					       optim_geom, fill_return_vec, &thread_pool, n_threads);
	 }
	 catch (const std::runtime_error &mess) {
	    std::cout << "Error in flexible ligand definition.\n";
	    std::cout << mess.what() << std::endl;
	    if (graphics_info_t::use_graphics_interface_flag) {
	       GtkWidget *w = wrapped_nothing_bad_dialog(mess.what());
	       gtk_widget_set_visible(w, TRUE);
	    }
	 }
      }
   }

   std::vector<coot::minimol::molecule> conformers = wlig.get_conformers();
   for (unsigned int iconf=0; iconf<conformers.size(); iconf++) {
      atom_selection_container_t asc = make_asc(conformers[iconf].pcmmdbmanager());
      int g_mol = graphics_info_t::create_molecule();
      std::string label = "conformer-";
      label+= coot::util::int_to_string(iconf);
      graphics_info_t::molecules[g_mol].install_model(g_mol, asc, g.Geom_p(), label, 1);
      mol_list.push_back(g_mol);
   }

   return mol_list;
}

/*! \brief make conformers of the ligand search molecules, each in its
  own molecule.

  Don't search the density. */
#ifdef USE_GUILE
SCM
ligand_search_make_conformers_scm() {

   std::vector<int> v = ligand_search_make_conformers_internal();
   return generic_int_vector_to_list_internal(v);
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *ligand_search_make_conformers_py() {

   std::vector<int> v = ligand_search_make_conformers_internal();
   return generic_int_vector_to_list_internal_py(v);
}
#endif // USE_PYTHON

void set_ligand_cluster_sigma_level(float f) { /* default 2.2 */
   graphics_info_t::ligand_cluster_sigma_level = f;
}

void set_find_ligand_mask_waters(int i) {

   graphics_info_t::find_ligand_mask_waters_flag = i;

}

/*! \brief set the atom radius for map masking */
void set_map_mask_atom_radius(float rad) {
   graphics_info_t::map_mask_atom_radius = rad;
}

/*! \brief get the atom radius for map masking, -99 is initial "unset" */
float map_mask_atom_radius() {
   return graphics_info_t::map_mask_atom_radius;
}


void set_ligand_flexible_ligand_n_samples(int i) {
   graphics_info_t::ligand_wiggly_ligand_n_samples = i;
}


void set_ligand_acceptable_fit_fraction(float f) {
   if (f >= 0.0 && f<= 1.0) {
      graphics_info_t::ligand_acceptable_fit_fraction = f;
      // std::cout << "ligand acceptable fit fraction set to " << f << "\n";
   } else {
      // std::cout << "ligand acceptable fit fraction " << f << " ignored\n";
   }
}

void set_find_ligand_n_top_ligands(int n) { /* fit the top n ligands,
					      not all of them, default
					      10. */
   graphics_info_t g;
   g.find_ligand_n_top_ligands = n;

}

void set_find_ligand_multi_solutions_per_cluster(float lim_1, float lim_2) {
   graphics_info_t g;
   g.find_ligand_multiple_solutions_per_cluster_flag = true;
   g.find_ligand_score_by_correl_frac_limit = lim_1;
   g.find_ligand_score_correl_frac_interesting_limit = lim_2;
}


/* flip the ligand (usually active residue) around its eigen vectoros
   to the next flip number.  Immediate replacement (like flip
   peptide). */
void flip_ligand(int imol, const char *chain_id, int resno) {

   // note that if the ligand's current flip_number is not 0, then we
   // need to undo a the current flip number before applying the next
   // one.  If the new flip_number is 0, we just undo of course.
   if (is_valid_model_molecule(imol)) {
      coot::minimol::molecule m =
	 graphics_info_t::molecules[imol].eigen_flip_residue(chain_id, resno);

      if (0) {
	 atom_selection_container_t asc = make_asc(m.pcmmdbmanager());
	 int g_mol = graphics_info_t::create_molecule();
	 std::string label = "flipping residue";
	 graphics_info_t g;
	 g.molecules[g_mol].install_model(g_mol, asc, g.Geom_p(), label, 1);
      }
   }
   graphics_draw();
}


//! \brief JED-flip intermediate atoms
//
// @ return a success status - 0 for "intermediate atoms were not not shown"
// 1 for intermediate atoms were shown (but not necessarily flipped
//
int jed_flip_intermediate_atoms() {
   return graphics_info_t::jed_flip_intermediate_atoms(false);
}

int reverse_jed_flip_intermediate_atoms() {
   return graphics_info_t::jed_flip_intermediate_atoms(true);
}

void jed_flip(int imol, const char *chain_id, int res_no, const char *ins_code,
	      const char *atom_name, const char *alt_conf, short int invert_selection) {

   if (is_valid_model_molecule(imol)) {
      bool invert_selection_flag = invert_selection;
      graphics_info_t g;
      std::string alt_conf_str(alt_conf);
      std::string atom_name_str(atom_name);
      coot::residue_spec_t spec(chain_id, res_no, ins_code);
      std::cout << "jed-flipping static atoms " << atom_name_str << std::endl;
      std::string problem_string = g.molecules[imol].jed_flip(spec, atom_name_str, alt_conf_str,
							      invert_selection_flag, g.Geom_p());
      if (! problem_string.empty()) {
	 add_status_bar_text(problem_string.c_str());
      }
   }
   graphics_draw();

}


//! \brief side-chain 180 flip on the active atom
int side_chain_flip_180_intermediate_atoms() {
   graphics_info_t g;
   return g.side_chain_flip_180_intermediate_atoms();
}


/*  ----------------------------------------------------------------------- */
/*                  Mask                                                    */
/*  ----------------------------------------------------------------------- */
/* The idea is to generate a new map that has been masked by some
   coordinates. */
/*!     (mask-map-by-molecule map-no mol-no invert?)  creates and
        displays a masked map, cuts down density where the coordinates
        are not.  If invert? is #t, cut the density down where the
        atoms are.  */
int mask_map_by_molecule(int map_mol_no, int coord_mol_no, short int invert_flag) {

   int imol_new_map = -1;

   // create a new map
   //
   // c.f. a new map that is created "Masked by Protein"
   // we use molecule_class_info_t's new_map()
   //

   // Where should the bulk of this function be?  Let's put it here for now.
   coot::ligand lig;
   graphics_info_t g;
   if (map_mol_no >= g.n_molecules()) {
      std::cout << "No such molecule (no map) at molecule number " << map_mol_no << std::endl;
   } else {

      if (coord_mol_no >= g.n_molecules()) {
	 std::cout << "No such molecule (no coords) at molecule number " << map_mol_no << std::endl;
      } else {

	 if (!g.molecules[map_mol_no].has_xmap()) {
	    std::cout << "No xmap in molecule number " << map_mol_no << std::endl;
	 } else {

	    if (! g.molecules[coord_mol_no].has_model()) {
	       std::cout << "No model in molecule number " << map_mol_no << std::endl;
	    } else {
	       short int mask_waters_flag; // treat the waters like protein atoms?
	       mask_waters_flag = g.find_ligand_mask_waters_flag;
	       lig.import_map_from(g.molecules[map_mol_no].xmap);
	       int selectionhandle = g.molecules[coord_mol_no].atom_sel.mol->NewSelection();

	       if (g.map_mask_atom_radius > 0) {
		  lig.set_map_atom_mask_radius(g.map_mask_atom_radius);
	       }

	       // make a selection:
	       std::string rnames = "*";
	       if (!mask_waters_flag)
		  rnames = "!HOH,WAT"; // treat waters differently to regular atoms.
	       g.molecules[coord_mol_no].atom_sel.mol->SelectAtoms(selectionhandle, 0, "*",
								   mmdb::ANY_RES, "*",
								   mmdb::ANY_RES, "*",
								   (char *) rnames.c_str(),
								   "*", "*", "*");

	       lig.mask_map(g.molecules[coord_mol_no].atom_sel.mol, selectionhandle, invert_flag);
	       g.molecules[coord_mol_no].atom_sel.mol->DeleteSelection(selectionhandle);
	       imol_new_map = g.create_molecule();
	       std::cout << "INFO:: Creating masked  map in molecule number " << imol_new_map << std::endl;
	       bool is_em_map_flag = g.molecules[map_mol_no].is_EM_map();
	       std::string old_name = g.molecules[map_mol_no].get_name();
	       std::string new_name = "Masked Map from " + old_name;
	       g.molecules[imol_new_map].install_new_map(lig.masked_map(), new_name, is_em_map_flag);
	       graphics_draw();
	    }
	 }
      }
   }
   return imol_new_map;
}

int
mask_map_by_atom_selection(int map_mol_no, int coords_mol_no, const char *mmdb_atom_selection, short int invert_flag) {

   int imol_new_map = -1;
   graphics_info_t g;
   if (is_valid_map_molecule(map_mol_no)) {
      if (is_valid_model_molecule(coords_mol_no)) {
	 // std::cout << "making lig..." << std::endl;
	 coot::ligand lig;
	 lig.import_map_from(g.molecules[map_mol_no].xmap);

	 if (g.map_mask_atom_radius > 0) {
	    lig.set_map_atom_mask_radius(g.map_mask_atom_radius);
	 }
	 int selectionhandle = g.molecules[coords_mol_no].atom_sel.mol->NewSelection();
	 g.molecules[coords_mol_no].atom_sel.mol->Select(selectionhandle, mmdb::STYPE_ATOM,
							 mmdb_atom_selection, mmdb::SKEY_NEW);
	 lig.mask_map(g.molecules[coords_mol_no].atom_sel.mol, selectionhandle, invert_flag);
	 imol_new_map = g.create_molecule();
	 std::string name = g.molecules[map_mol_no].name_;
	 std::string new_name = name + " Masked Map";
	 bool is_em_map_flag = g.molecules[map_mol_no].is_EM_map();
         float contour_level = 0.99f * g.molecules[map_mol_no].get_contour_level();
         g.molecules[imol_new_map].install_new_map_with_contour_level(lig.masked_map(), new_name, contour_level, is_em_map_flag);
	 graphics_draw();
      } else {
	 std::cout << "No model molecule in " << coords_mol_no << std::endl;
      }
   } else {
      std::cout << "No map molecule in " << map_mol_no << std::endl;
   }
   return imol_new_map;
}


#include "c-interface-gui.hh"

void do_find_ligands_dialog() {

   GtkWidget *dialog = widget_from_builder("find_ligand_dialog");
   int istate = fill_ligands_dialog(dialog); /* return OK, we have map(s), ligand(s), masking(s) */

   if (istate == 0) {
      gtk_widget_set_visible(dialog, FALSE);
      std::string s("Problem finding maps, coords or ligands!");
      graphics_info_t g;
      g.add_status_bar_text(s);
      std::cout << s << std::endl;
   } else {

      // std::cout << "do_find_ligands_dialog()  showing dialog " << dialog << std::endl;
      set_transient_for_main_window(dialog);

      gtk_widget_set_visible(dialog, TRUE);
   }

}

/*  ----------------------------------------------------------------------- */
/*                  monomer lib                                             */
/*  ----------------------------------------------------------------------- */
std::vector<std::pair<std::string, std::string> > monomer_lib_3_letter_codes_matching(const std::string &search_string, short int allow_minimal_descriptions_flag) {

   graphics_info_t g;
   std::vector<std::pair<std::string, std::string > > v = g.Geom_p()->matching_names(search_string, allow_minimal_descriptions_flag);

   // search the list of monomers for text search_text:

   return v;
}

#include "c-interface-image-widget.hh"

// we allocate new memory here, without ever giving it back.  The
// memory should be freed when the dialog is destroyed.
//
int
handle_make_monomer_search(const char *text, GtkWidget *viewport) {

   unsigned int image_size = 128;
#if __APPLE__
   // image_size = 450;
#endif

   int stat = 0;
   std::string t(text);
   GtkWidget *vbox_current = widget_from_builder("monomer_search_results_vbox");
   bool allow_minimal_descriptions_flag = false;
   GtkWidget *dialog = widget_from_builder("monomer_search_dialog");
   gtk_widget_set_size_request(dialog, 500, 300);

   graphics_info_t g;
   std::vector<std::pair<std::string, std::string> > v =
      monomer_lib_3_letter_codes_matching(t, allow_minimal_descriptions_flag);

   if (false)
      std::cout << "DEBUG::  " << " found " << v.size() << " matching molecules "
		<< " using string :" << t << ":"
		<< std::endl;

   GtkWidget *vbox = vbox_current;
   g.clear_out_container(vbox);

   GtkWidget *grid = gtk_grid_new();
   gtk_grid_set_row_spacing(GTK_GRID(grid), 2);
   gtk_grid_set_column_spacing(GTK_GRID(grid), 2);
   gtk_box_append(GTK_BOX(vbox), grid);
   for (unsigned int i=0; i<v.size(); i++) {
      const std::string &tlc           = v[i].first;
      const std::string &molecule_name = v[i].second;
      std::string *s = new std::string(tlc); // memory leak
      GtkWidget *wp = get_image_widget_for_comp_id(tlc, coot::protein_geometry::IMOL_ENC_ANY, image_size);
      if (wp) {
         GtkWidget *button_for_image = gtk_button_new();
         gtk_button_set_child(GTK_BUTTON(button_for_image), wp);
         // gtk_widget_set_size_request(wp, image_size, image_size);
         // gtk_widget_set_hexpand(wp, TRUE);
         // gtk_widget_set_vexpand(wp, TRUE);
         gtk_image_set_pixel_size(GTK_IMAGE(wp), image_size); // Example: request 64 logical pixels
         gtk_grid_attach(GTK_GRID(grid), button_for_image, 0, i, 1, 1);
         g_signal_connect(G_OBJECT(button_for_image), "clicked", G_CALLBACK(on_monomer_lib_search_results_button_press), s);
      } else {
         // 20251211-PE too noisy. Many tlc don't have dictionaries
         // logger.log(log_t::WARNING, "failed to get image widget for", tlc);
      }
      std::string button_label = tlc + std::string(": ") + molecule_name;
      GtkWidget *button_for_label = gtk_button_new();
      GtkWidget *label = gtk_label_new(button_label.c_str());
      gtk_button_set_child(GTK_BUTTON(button_for_label), label);
      gtk_widget_set_halign(label, GTK_ALIGN_START);
      gtk_grid_attach(GTK_GRID(grid), button_for_label, 1, i, 1, 1);
      g_signal_connect(G_OBJECT(button_for_label), "clicked", G_CALLBACK(on_monomer_lib_search_results_button_press), s);
   }

   return stat;
}


void
on_monomer_lib_sbase_molecule_button_press (GtkButton *button,
					    gpointer user_data) {
   std::string *s = (std::string *) user_data;
   get_sbase_monomer(s->c_str());

}

void
on_monomer_lib_search_results_button_press (GtkButton *button,
					    gpointer user_data) {

   std::string *s = (std::string *) user_data;
   get_monomer(s->c_str());

}

// Sigh.  This is wrong. We don't want to replace the coordinates - we
// do that in molecule_class_info_t::replace_coords()
//
// But this may be useful in future
// (not not yet compile tested)
//
// atom_selection_container_t rigid_body_asc;

//       std::vector<mmdb::Atom *> atoms;
//       mmdb::Atom *at;

//       for (int n=0; n<g.molecules[imol].atom_sel.n_selected_atoms; n++) {
// 	 at = g.molecules[imol].atom_sel.atom_selection[n];
// 	 std::string mmdb_chain(at->residue-GetChainID());
// 	 for(int ifrag=0; ifrag<mol.fragments.size(); ifrag++) {
// 	    for (int ires=1; ires<=mol[ifrag].n_residues(); ires++) {
// 	       for (int iatom=0; iatom<mol[ifrag][ires].atoms.size(); iatom++) {
// 		  if (mmdb_chain == mol[ifrag].fragment_id) {
// 		     if (ires == at->residue->seqNum) {
// 			std::string mmdb_name(at->name);
// 			if (mmdb_name == mol[ifrag][ires][iatom].name) {
// 			   atoms.push_back(at);
// 			}
// 		     }
// 		  }
// 	       }
// 	    }
// 	 }
//       }

//       mmdb::PPAtom atom_selection = new mmdb::Atom *[atoms.size()];
//       for (int i=0; i<atoms.size(); i++)
// 	 atom_selection[i] = atoms[i];

//       rigid_body_asc.atom_selection = atom_selection;
//       rigid_body_asc.n_selected_atoms = atoms.size();
//       std::cout << "INFO:: make_atom_selection found " << atoms.size()
// 		<< " atoms" << std::endl;

//       rigid_body_asc.mol = g.molecules[g.imol_rigid_body_refine].atom_sel.mol;


// c.f. molecule_class_info_t::atom_index
//
// This is not a member class of molecule_class_info_t because we dont
// want minimol dependency in molecule_class_info_t and because we
// want the number of atoms *and* the pointer to be returned... messy
// messy - lets do that here.
//
// We only use the n_selected_atoms and atom_selection fields of this
// atom_selection_container_t.
//
// atom_selection_container_t
// make_atom_selection(int imol, const coot::minimol::molecule &mol) {

// }

//       std::vector<coot::minimol::atom *> atoms = range_mol.select_atoms_serial();
//       for (int i=0; i<atoms.size(); i++) {
// 	 std::cout << "range mol atom: " << atoms[i]->pos.format() << std::endl;
//       }


#ifdef USE_GUILE
SCM add_dipole_scm(int imol, const char* chain_id, int res_no, const char *ins_code) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> res_specs;
      coot::residue_spec_t rs(chain_id, res_no, ins_code);
      res_specs.push_back(rs);
      graphics_info_t g;
      std::pair<coot::dipole, int> dp =
	 g.molecules[imol].add_dipole(res_specs, *g.Geom_p());
      r = dipole_to_scm(dp);
   }
   graphics_draw();
   return r;
}
#endif

#ifdef USE_GUILE
SCM dipole_to_scm(std::pair<coot::dipole, int> dp) {

   SCM r = SCM_EOL;

   clipper::Coord_orth co = dp.first.get_dipole();
   SCM co_scm = SCM_EOL;
   co_scm = scm_cons(scm_from_double(co.z()), co_scm);
   co_scm = scm_cons(scm_from_double(co.y()), co_scm);
   co_scm = scm_cons(scm_from_double(co.x()), co_scm);
   r = scm_cons(co_scm, r);
   r = scm_cons(scm_from_int(dp.second), r);
   return r;
}
#endif


#ifdef USE_GUILE
/*! \brief generate a dipole from all atoms in the given residues. */
SCM add_dipole_for_residues_scm(int imol, SCM residue_specs) {

   // int r = -1;
   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> res_specs;
      SCM residue_spec_length = scm_length(residue_specs);
      int l = scm_to_int(residue_spec_length);
      for (int ispec=0; ispec<l; ispec++) {
	 SCM residue_spec_scm = scm_list_ref(residue_specs, scm_from_int(ispec));
	 coot::residue_spec_t rs = residue_spec_from_scm(residue_spec_scm);
	 res_specs.push_back(rs);
      }
      graphics_info_t g;
      std::pair<coot::dipole, int> d =
	 g.molecules[imol].add_dipole(res_specs, *g.Geom_p());
      r = dipole_to_scm(d);
   }
   graphics_draw();
   return r;
}
#endif /* USE_GUILE */

#ifdef USE_PYTHON
PyObject *add_dipole_py(int imol, const char* chain_id, int res_no, const char *ins_code) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> res_specs;
      coot::residue_spec_t rs(chain_id, res_no, ins_code);
      res_specs.push_back(rs);
      graphics_info_t g;
      std::pair<coot::dipole, int> dp =
	 g.molecules[imol].add_dipole(res_specs, *g.Geom_p());
      // set r
      r = dipole_to_py(dp);
   }
   graphics_draw();
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *dipole_to_py(std::pair<coot::dipole, int> dp) {

  PyObject *r = PyList_New(2);

  clipper::Coord_orth co = dp.first.get_dipole();
  PyObject *co_py = PyList_New(3);
  PyList_SetItem(co_py, 0, PyFloat_FromDouble(co.x()));
  PyList_SetItem(co_py, 1, PyFloat_FromDouble(co.y()));
  PyList_SetItem(co_py, 2, PyFloat_FromDouble(co.z()));
  PyList_SetItem(r, 0, PyLong_FromLong(dp.second));
  PyList_SetItem(r, 1, co_py);
  return r;
}
#endif


#ifdef USE_PYTHON
/*! \brief generate a dipole from all atoms in the given residues. */
PyObject *add_dipole_for_residues_py(int imol, PyObject *residue_specs) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> res_specs;
      unsigned int l = PyObject_Length(residue_specs);
      for (unsigned int ispec=0; ispec<l; ispec++) {
	 PyObject *residue_spec_py = PyList_GetItem(residue_specs, ispec);
	 coot::residue_spec_t rs = residue_spec_from_py(residue_spec_py);
	 res_specs.push_back(rs);
      }
      graphics_info_t g;
      std::pair<coot::dipole, int> dp =
	 g.molecules[imol].add_dipole(res_specs, *g.Geom_p());
      // set r as a python object here.
      r = dipole_to_py(dp);
   }
   graphics_draw();
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif /* USE_PYTHON */


void delete_dipole(int imol, int dipole_number) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].delete_dipole(dipole_number);
   }
   graphics_draw();
}

#ifdef USE_GUILE
SCM non_standard_residue_names_scm(int imol) {

   SCM r = SCM_EOL;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      std::vector<std::string> resnames =
	 coot::util::non_standard_residue_types_in_molecule(mol);

      // remove water if it is there
      std::vector<std::string>::iterator it =
	 std::find(resnames.begin(), resnames.end(), "HOH");
      if (it != resnames.end())
	 resnames.erase(it);

      r = generic_string_vector_to_list_internal(resnames);
   }

   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *non_standard_residue_names_py(int imol) {

  PyObject *r = PyList_New(0);
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      std::vector<std::string> resnames =  coot::util::non_standard_residue_types_in_molecule(mol);

      // remove water if it is there
      std::vector<std::string>::iterator it =
	 std::find(resnames.begin(), resnames.end(), "HOH");
      if (it != resnames.end())
	 resnames.erase(it);

      r = generic_string_vector_to_list_internal_py(resnames);
   }

   return r;
}
#endif



#ifdef USE_GUILE
/*! \brief return a list of compoundIDs of in SBase of which the
  given string is a substring of the compound name */
SCM matching_compound_names_from_sbase_scm(const char *compound_name_fragment) {

   graphics_info_t g;
   std::vector<std::pair<std::string, std::string> > matching_comp_ids =
      g.Geom_p()->matching_ccp4srs_residues_names(compound_name_fragment);

   std::vector<std::string> rv;
   for (unsigned int i=0; i<matching_comp_ids.size(); i++)
      rv.push_back(matching_comp_ids[i].first);

   SCM r = generic_string_vector_to_list_internal(rv);
   return r;
}
#endif

#ifdef USE_PYTHON
/*! \brief return a list of compoundIDs of in SBase of which the
  given string is a substring of the compound name */
PyObject *matching_compound_names_from_sbase_py(const char *compound_name_fragment) {

   graphics_info_t g;
   std::vector<std::pair<std::string, std::string> > matching_comp_ids =
      g.Geom_p()->matching_ccp4srs_residues_names(compound_name_fragment);

   std::vector<std::string> rv;
   for (unsigned int i=0; i<matching_comp_ids.size(); i++)
      rv.push_back(matching_comp_ids[i].first);

   PyObject *r = generic_string_vector_to_list_internal_py(rv);
   return r;
}
#endif /* USE_PYTHON */

#ifdef USE_GUILE
/*! \brief return a list of compoundIDs in the dictionary of which the
  given string is a substring of the compound name */
SCM matching_compound_names_from_dictionary_scm(const char *compound_name_fragment,
						short int allow_minimal_descriptions_flag) {

   graphics_info_t g;
   std::vector<std::pair<std::string, std::string> > matching_comp_ids =
      g.Geom_p()->matching_names(compound_name_fragment, allow_minimal_descriptions_flag);

   if (false) { // this was needed to see the bug in mmdb parser on long molecule names.
      for (unsigned int i=0; i<matching_comp_ids.size(); i++) {
         std::cout << " " << i << " " << matching_comp_ids[i].first << " " << matching_comp_ids[i].second
                   << std::endl;
      }
   }
   std::vector<std::string> rv;
   for (unsigned int i=0; i<matching_comp_ids.size(); i++)
      rv.push_back(matching_comp_ids[i].first);

   SCM r = generic_string_vector_to_list_internal(rv);
   return r;

}

#endif /* USE_GUILE */


/*! \brief add residue name to the list of residue names that don't
  get auto-loaded from the Refmac dictionary. */
void add_non_auto_load_residue_name(const char *s) {
   graphics_info_t g;
   g.Geom_p()->add_non_auto_load_residue_name(s);
}

/*! \brief remove residue name from the list of residue names that don't
  get auto-loaded from the Refmac dictionary. */
void remove_non_auto_load_residue_name(const char *s) {
   graphics_info_t g;
   g.Geom_p()->remove_non_auto_load_residue_name(s);
}

/*! \brief try to auto-load the dictionary for comp_id from the refmac monomer library.

   return 0 on failure.
   return 1 on successful auto-load.
   return 2 on already-read.
   */
int auto_load_dictionary(const char *comp_id) {

   graphics_info_t g;
   int r = 0;
   if (comp_id) {
      std::string s = comp_id;
      int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
      if (g.Geom_p()->have_dictionary_for_residue_type_no_dynamic_add(s, imol_enc)) {
	 r = 2;
      } else {
	 int status = g.Geom_p()->try_dynamic_add(s, g.cif_dictionary_read_number);
	 if (status)
	    r = 1;
      }
   }
   return r;
}

/*! \brief as above, but dictionary is cleared and re-read if it already exists */
int reload_dictionary(const char *comp_id) {

   graphics_info_t g;
   int r = 0;
   if (comp_id) {
      std::string s = comp_id;
      int status = g.Geom_p()->try_dynamic_add(s, g.cif_dictionary_read_number);
      if (status)
	 r = 1;
   }
   return r;

}




#ifdef USE_PYTHON
/*! \brief return a list of compoundIDs in the dictionary which the
  given string is a substring of the compound name */
PyObject *matching_compound_names_from_dictionary_py(const char *compound_name_fragment,
						     short int allow_minimal_descriptions_flag) {

   graphics_info_t g;
   std::vector<std::pair<std::string, std::string> > matching_comp_ids =
      g.Geom_p()->matching_names(compound_name_fragment,
				 allow_minimal_descriptions_flag);

   std::vector<std::string> rv;
   for (unsigned int i=0; i<matching_comp_ids.size(); i++)
      rv.push_back(matching_comp_ids[i].first);

   PyObject *r = generic_string_vector_to_list_internal_py(rv);
   return r;
}
#endif /* USE_PYTHON */



int get_sbase_monomer(const char *comp_id) {

   return get_ccp4srs_monomer_and_dictionary(comp_id);

}

int get_ccp4srs_monomer_and_dictionary(const char *comp_id) {

   int imol = -1;

   graphics_info_t g;
   mmdb::Residue *residue_p = g.Geom_p()->get_ccp4srs_residue(comp_id);
   if (residue_p) {
      mmdb::Manager *mol = new mmdb::Manager;
      mmdb::Model *model_p = new mmdb::Model;
      mmdb::Chain *chain_p = new mmdb::Chain;
      residue_p->SetResID(comp_id, 1, "");
      chain_p->AddResidue(residue_p);
      chain_p->SetChainID("A");
      model_p->AddChain(chain_p);
      mol->AddModel(model_p);
      imol = g.create_molecule();
      std::string name = "Monomer ";
      name += coot::util::upcase(comp_id);
      graphics_info_t::molecules[imol].install_model(imol, make_asc(mol), g.Geom_p(), name, 1, 0);
      move_molecule_to_screen_centre_internal(imol);
      g.Geom_p()->fill_using_ccp4srs(comp_id); // if you can (should be able to)
      graphics_draw();
   }

   return imol;
}


#ifdef USE_GUILE
/*! \brief return the monomer name
  return scheme false if not found */
SCM comp_id_to_name_scm(const char *comp_id) {
   SCM r = SCM_BOOL_F;
   graphics_info_t g;
   int imol = coot::protein_geometry::IMOL_ENC_ANY;
   std::pair<bool, std::string> name = g.Geom_p()->get_monomer_name(comp_id, imol);

   if (! name.first) {
      g.Geom_p()->try_dynamic_add(comp_id, g.cif_dictionary_read_number);
      g.cif_dictionary_read_number++;
      name = g.Geom_p()->get_monomer_name(comp_id, imol);
   }
   if (name.first) {
      r = scm_from_locale_string(name.second.c_str());
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
/*! \brief return the monomer name

  return python false if not found */
PyObject *comp_id_to_name_py(const char *comp_id) {

   PyObject *r = Py_False;

   graphics_info_t g;
   int imol = coot::protein_geometry::IMOL_ENC_ANY;
   std::pair<bool, std::string> name = g.Geom_p()->get_monomer_name(comp_id, imol);
   if (name.first) {
      r = myPyString_FromString(name.second.c_str());
   }
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;
}

#endif // USE_PYTHON



/*  ----------------------------------------------------------------------- */
//                  animated ligand interactions
/*  ----------------------------------------------------------------------- */
void add_animated_ligand_interaction(int imol, const pli::fle_ligand_bond_t &lb) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].add_animated_ligand_interaction(lb);
   }
}




/* Just pondering, for use with exporting ligands from the 2D sketcher
   to the main coot window. It should find the residue that this
   residue is sitting on top of that is in a molecule that has lot of
   atoms (i.e. is a protein) and create a new molecule that is a copy
   of the molecule without the residue/ligand that this (given)
   ligand overlays - and a copy of this given ligand.  */
int exchange_ligand(int imol_lig, const char *chain_id_lig, int resno_lig, const char *ins_code_lig) {

   int imol = -1;

   return imol;

}

/*! \brief Match ligand atom names

  By using graph matching, make the names of the atoms of the
  given ligand/residue match those of the reference residue/ligand as
  closely as possible - where there would be a atom name clash, invent
  a new atom name.
 */
void match_ligand_atom_names(int imol_ligand, const char *chain_id_ligand, int resno_ligand, const char *ins_code_ligand,
			     int imol_reference, const char *chain_id_reference, int resno_reference, const char *ins_code_reference) {

   if (! is_valid_model_molecule(imol_ligand)) {
      std::cout << "Not a valid model number " << imol_ligand << std::endl;
   } else {
      if (! is_valid_model_molecule(imol_reference)) {
	      std::cout << "Not a valid model number " << imol_reference << std::endl;
      } else {
	 graphics_info_t g;
	 mmdb::Residue *res_ref =
	    g.molecules[imol_reference].get_residue(chain_id_reference, resno_reference, ins_code_reference);
	 if (! res_ref) {
	    std::cout << "No reference residue " << chain_id_reference << " " << resno_reference
		      << " "  << ins_code_reference << std::endl;
	 } else {
	    // now lock res_ref - when multi-threaded
	    g.molecules[imol_ligand].match_ligand_atom_names(chain_id_ligand, resno_ligand, ins_code_ligand, res_ref);
	    graphics_draw();
	 }
      }
   }
}

/*! \brief Match ligand atom names to a reference ligand type (comp_id)

  By using graph matching, make the names of the atoms of the
  given ligand/residue match those of the reference ligand from the
  geometry store as closely as possible. Where there would be an
  atom name clash, invent a new atom name.
 */
void match_ligand_atom_names_to_comp_id(int imol_ligand, const char *chain_id_ligand, int resno_ligand, const char *ins_code_ligand, const char *comp_id_ref) {

   if (! is_valid_model_molecule(imol_ligand)) {
      std::cout << "Not a valid model number " << imol_ligand << std::endl;
   } else {
      graphics_info_t g;
      mmdb::Residue *res_ref = NULL;

      g.Geom_p()->try_dynamic_add(comp_id_ref, g.cif_dictionary_read_number++);
      std::pair<bool, coot::dictionary_residue_restraints_t> p =
	 g.Geom_p()->get_monomer_restraints(comp_id_ref, imol_ligand);

      if (p.first)
	 res_ref = p.second.GetResidue(true, g.default_new_atoms_b_factor);

      if (! res_ref) {
	 std::cout << "No reference residue for comp_id " << comp_id_ref << std::endl;
      } else {

	 // now lock res_ref - when multi-threaded
	 g.molecules[imol_ligand].match_ligand_atom_names(chain_id_ligand, resno_ligand, ins_code_ligand, res_ref);
	 graphics_draw();
      }
   }
}



/*  return a list of chiral centre ids as determined from topological
    equivalence analysis based on the bond info (and element names).

    A better OO design would put this inside protein-geometry.

*/
std::vector<std::string>
topological_equivalence_chiral_centres(const std::string &residue_type) {

   std::vector<std::string> centres;
#ifdef HAVE_GOOCANVAS

   graphics_info_t g;

   int imol = 0; // dummy
   std::pair<bool, coot::dictionary_residue_restraints_t> p =
      g.Geom_p()->get_monomer_restraints(residue_type, imol);

   if (p.first) {

      const coot::dictionary_residue_restraints_t &restraints = p.second;
      lig_build::molfile_molecule_t mm(restraints); // makes dummy atoms
      widgeted_molecule_t wm(mm, NULL); // we have the atom names already.
      topological_equivalence_t top_eq(wm.atoms, wm.bonds);
      centres = top_eq.chiral_centres();

      std::cout << "-------- chiral centres by topology analysis -----------"
		<< std::endl;
      for (unsigned int ic=0; ic<centres.size(); ic++) {
	 std::cout << "     " << ic << "  " << centres[ic] << std::endl;
      }
      std::cout << "-------------------" << std::endl;
   }

#endif // HAVE_GOOCANVAS
   return centres;
}


// This doesn't reset the view, perhaps it should (normally, we do not
// reset the view when using install_model).
//
int read_small_molecule_cif(const char *file_name) {

   int imol = -1;
   coot::smcif smcif;
   mmdb::Manager *mol = smcif.read_sm_cif(file_name);

   if (mol) {
      graphics_info_t g;
      imol = g.create_molecule();
      // do the "Show Symmetry" and Bond (ball) smoothness for the user:
      if (g.bond_smoothness_factor == 1)
         g.bond_smoothness_factor = 2;
      short int display_in_diplay_control = 1;
      bool is_from_shelx_ins_flag = false;
      g.show_symmetry = 1;
      g.molecules[imol].install_model(imol, mol, g.Geom_p(), file_name, display_in_diplay_control, is_from_shelx_ins_flag);
      // g.molecules[imol].update_symmetry();
      update_go_to_atom_window_on_new_mol();
      graphics_draw();
   }

   return imol;
}



// This doesn't reset the view, perhaps it should (normally, we do not
// reset the view when using install_model).
//
int read_small_molecule_data_cif(const char *file_name) {

   int imol = -1;
   coot::smcif smcif;
   bool state = smcif.read_data_sm_cif(file_name);
   if (state) {
      graphics_info_t g;
      imol = g.create_molecule();

      bool is_em_map_flag = false;
      std::pair<clipper::Xmap<float>, clipper::Xmap<float> > maps = smcif.sigmaa_maps();
      if (not (maps.first.is_null())) {
	 std::string map_name = file_name;
	 map_name += " SigmaA";
	 g.molecules[imol].install_new_map(maps.first, map_name, is_em_map_flag);
	 g.scroll_wheel_map = imol;
	 int imol_diff = g.create_molecule();
	 map_name = file_name;
	 map_name += " Diff-SigmaA";
	 g.molecules[imol_diff].install_new_map(maps.second, map_name, is_em_map_flag);
	 g.molecules[imol_diff].set_map_is_difference_map(true);
      }
      graphics_draw();
   }
   return imol;
}

int read_small_molecule_data_cif_and_make_map_using_coords(const char *file_name,
							   int imol_coords) {

   graphics_info_t g;
   int imol_map = -1;
   if (is_valid_model_molecule(imol_coords)) {

      // we make a new mol because somewhere in read_data_sm_cif(), the atom
      // selections get trashed.
      //
      mmdb::Manager *new_mol = new mmdb::Manager;
      new_mol->Copy(g.molecules[imol_coords].atom_sel.mol, mmdb::MMDBFCM_All);
      atom_selection_container_t asc = make_asc(new_mol);
      mmdb::Atom** atom_selection = asc.atom_selection;
      int n_selected_atoms = asc.n_selected_atoms;
      coot::smcif smcif;
      bool state = smcif.read_data_sm_cif(file_name);

      std::pair<clipper::Xmap<float>, clipper::Xmap<float> > maps =
	 smcif.sigmaa_maps_by_calc_sfs(atom_selection, n_selected_atoms);
      if (not (maps.first.is_null())) {
	 imol_map = g.create_molecule();
	 std::string map_name = file_name;
	 map_name += " SigmaA";
	 g.molecules[imol_map].install_new_map(maps.first, map_name, false);
	 g.scroll_wheel_map = imol_map;
	 int imol_diff = g.create_molecule();
	 map_name = file_name;
	 map_name += " Diff-SigmaA";
	 g.molecules[imol_diff].install_new_map(maps.second, map_name, false);
	 g.molecules[imol_diff].set_map_is_difference_map(true);
      }
   }
   return imol_map;
}




#ifdef USE_GUILE
void
multi_residue_torsion_scm(int imol, SCM residues_specs_scm) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::vector<coot::residue_spec_t> residue_specs =
	 scm_to_residue_specs(residues_specs_scm);
      g.multi_torsion_residues(imol, residue_specs);

      graphics_draw();
   }

}
#endif

// /*! \brief set the torsion of the given residues. The
// torsion_residues_specs_scm is a list of: (list  torsion-angle atom-spec-1 atom-spec-2
// atom-spec-2 atom-spec-3).
// */
// #ifdef USE_GUILE
// SCM set_multi_torsion_scm(int imol, SCM torsion_residues_specs_scm) {

//    SCM r = SCM_BOOL_F;
//    if (is_valid_model_molecule(imol)) {
//       graphics_info_t g;
//       if (scm_is_true(scm_list_p(torsion_residues_specs_scm))) {
// 	 std::vector<coot::torsion_atom_specs_t> tas;
// 	 SCM specs_length_scm = scm_length(torsion_residues_specs_scm);
// 	 unsigned int specs_length = scm_to_int(specs_length_scm);
// 	 for (unsigned int i=0; i<specs_length; i++) {
// 	    SCM torsion_spec_scm = scm_list_ref(torsion_residues_specs_scm, scm_from_int(i));
// 	    if (! scm_is_true(scm_list_p(torsion_spec_scm))) {
// 	       break;
// 	    } else {
// 	       // happy path
// 	       SCM torsion_spec_length_scm = scm_length(torsion_spec_scm);
// 	       unsigned int torsion_spec_length = scm_to_int(torsion_spec_length_scm);
// 	       if (torsion_spec_length == 5) {
// 		  SCM torsion_angle_scm = scm_car(torsion_spec_scm);
// 		  if (scm_is_true(scm_number_p(torsion_angle_scm))) {
// 		     double ta = scm_to_double(torsion_angle_scm);
// 		     SCM atom_spec_1_scm = scm_list_ref(torsion_spec_scm, scm_from_int(1));
// 		     SCM atom_spec_2_scm = scm_list_ref(torsion_spec_scm, scm_from_int(2));
// 		     SCM atom_spec_3_scm = scm_list_ref(torsion_spec_scm, scm_from_int(3));
// 		     SCM atom_spec_4_scm = scm_list_ref(torsion_spec_scm, scm_from_int(4));
// 		     coot::atom_spec_t atom_spec_1 = atom_spec_from_scm_expression(atom_spec_1_scm);
// 		     coot::atom_spec_t atom_spec_2 = atom_spec_from_scm_expression(atom_spec_2_scm);
// 		     coot::atom_spec_t atom_spec_3 = atom_spec_from_scm_expression(atom_spec_3_scm);
// 		     coot::atom_spec_t atom_spec_4 = atom_spec_from_scm_expression(atom_spec_4_scm);
// 		     coot::torsion_atom_specs_t t(atom_spec_1, atom_spec_2, atom_spec_3, atom_spec_4, ta);
// 		     tas.push_back(t);
// 		  }
// 	       }
// 	    }
// 	 }
// 	 if (specs_length == tas.size()) {
// 	    g.molecules[imol].set_multi_torsion(tas);
// 	 }
//       }
//    }
//    return r;
// }
// #endif // USE_GUILE


#ifdef USE_PYTHON
void
multi_residue_torsion_py(int imol, PyObject *residues_specs_py) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::vector<coot::residue_spec_t> residue_specs = py_to_residue_specs(residues_specs_py);
      g.multi_torsion_residues(imol, residue_specs);
      graphics_draw();
   }
}

#endif // USE_PYTHON


void
multi_residue_torsion_fit(int imol, const std::vector<coot::residue_spec_t> &residue_specs, int n_trials) {

   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map())) {
	 graphics_info_t g;
	 const clipper::Xmap<float> &xmap =
	    g.molecules[imol_refinement_map()].xmap;
	 g.molecules[imol].multi_residue_torsion_fit(residue_specs, xmap, n_trials, g.Geom_p());
	 graphics_draw();
      }
   }
}


/*! \brief fit residues

(note: fit to the current-refinement map)
*/
#ifdef USE_GUILE
SCM multi_residue_torsion_fit_scm(int imol, SCM residues_specs_scm, int n_trials) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map())) {
	 graphics_info_t g;
	 std::vector<coot::residue_spec_t> residue_specs =
	    scm_to_residue_specs(residues_specs_scm);
	 const clipper::Xmap<float> &xmap =
	    g.molecules[imol_refinement_map()].xmap;
	 g.molecules[imol].multi_residue_torsion_fit(residue_specs, xmap, n_trials, g.Geom_p());
	 graphics_draw();
	 r = SCM_BOOL_T; // we did something at least.
      }
   }
   return r;
}
#endif


/*! \brief fit residues

(note: fit to the current-refinement map)
*/
#ifdef USE_PYTHON
PyObject *multi_residue_torsion_fit_py(int imol, PyObject *residues_specs_py, int n_trials) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_refinement_map())) {
	 graphics_info_t g;
	 std::vector<coot::residue_spec_t> residue_specs = py_to_residue_specs(residues_specs_py);
	 const clipper::Xmap<float> &xmap = g.molecules[imol_refinement_map()].xmap;
	 g.molecules[imol].multi_residue_torsion_fit(residue_specs, xmap, n_trials, g.Geom_p());
	 graphics_draw();
	 r = Py_True; // we did something at least.
      }
   }
   Py_INCREF(r);
   return r;
}
#endif // PYTHON


void
setup_multi_residue_torsion() {

   graphics_info_t g;
   g.in_multi_residue_torsion_define = true; // we are about to pick atoms
   g.multi_residue_torsion_picked_residue_specs.clear(); // clear previous picked residues
   pick_cursor_maybe(); // depends on ctrl key for rotate

   g.multi_residue_torsion_reverse_fragment_mode = false;
   // GtkWidget *w = create_multi_residue_torsion_pick_dialog();
   GtkWidget *w = widget_from_builder("multi_residue_torsion_pick_dialog");
   gtk_widget_set_visible(w, TRUE);
}


/* show the rotatable bonds dialog */
void
show_multi_residue_torsion_dialog() {

   graphics_info_t g;
   if (g.multi_residue_torsion_picked_residue_specs.size()) {
      g.multi_torsion_residues(g.multi_residue_torsion_picked_residues_imol,
			       g.multi_residue_torsion_picked_residue_specs);
      g.in_multi_residue_torsion_mode = true;
   }
   graphics_draw();
}

void clear_multi_residue_torsion_mode() {

   graphics_info_t g;
   g.in_multi_residue_torsion_mode = false;
}

void set_multi_residue_torsion_reverse_mode(short int mode) {

   graphics_info_t g;
   g.multi_residue_torsion_reverse_fragment_mode = mode;
}


/* ------------------------------------------------------------------------- */
/*                      prodrg import function                               */
/* ------------------------------------------------------------------------- */
// the function passed to lbg, which is called when a new
// prodrg-in.mdl file has been made.  We no longer have a timeout
// function waiting for prodrg-in.mdl to be updated/written.
//
void prodrg_import_function(std::string file_name, std::string comp_id) {

   std::string func_name = "import-from-3d-generator-from-mdl";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(file_name));
   args.push_back(single_quote(comp_id));
   coot::scripting_function(func_name, args);
}


/* ------------------------------------------------------------------------- */
/*                       SBase import function                               */
/* ------------------------------------------------------------------------- */
// the function passed to lbg, so that it calls it when a new
// SBase comp_id is required.  We no longer have a timeout
// function waiting for prodrg-in.mdl to be updated/written.
//
void sbase_import_function(std::string comp_id) {

   bool done = false;
#ifdef USE_PYTHON
   if (graphics_info_t::prefer_python) {
      std::string s = "get_sbase_monomer_and_overlay(";
      s += single_quote(comp_id);
      s += ")";
      safe_python_command(s);
      done = true;
   }
#endif

#ifdef USE_GUILE
   if (! done) {
      std::string s = "(get-ccp4srs-monomer-and-overlay ";
      s += single_quote(comp_id);
      s += ")";
      safe_scheme_command(s);
   }
#endif

}


// return a spec for the first residue with the given type.
// test the returned spec for unset_p().
//
coot::residue_spec_t
get_residue_by_type(int imol, const std::string &residue_type) {

   coot::residue_spec_t spec;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      spec = g.molecules[imol].get_residue_by_type(residue_type);
   }
   return spec;
}


std::vector<coot::residue_spec_t> get_residue_specs_in_mol(int imol, const std::string &residue_type) {

   std::vector<coot::residue_spec_t> v;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      v = g.molecules[imol].get_residues_by_type(residue_type);
   }
   return v;
}

#ifdef USE_PYTHON
PyObject *get_residue_specs_in_mol_py(int imol,
				      const std::string &residue_type) {

   std::vector<coot::residue_spec_t> v =
      get_residue_specs_in_mol(imol, residue_type);

   PyObject *r = PyList_New(v.size());
   for (unsigned int i=0; i<v.size(); i++)
      PyList_SetItem(r, i, residue_spec_to_py(v[i]));

   return r;
}
#endif



#ifdef USE_GUILE
SCM get_residue_by_type_scm(int imol, const std::string &residue_type) {

   SCM r = SCM_BOOL_F;
   coot::residue_spec_t spec = get_residue_by_type(imol, residue_type);
   if (! spec.unset_p())
      r = residue_spec_to_scm(spec);
   return r;
}
#endif


#ifdef USE_PYTHON
PyObject *get_residue_by_type_py(int imol, const std::string &residue_type) {

   PyObject *r = Py_False;
   coot::residue_spec_t spec = get_residue_by_type(imol, residue_type);
   if (! spec.unset_p())
      r = residue_spec_to_py(spec);
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }

   return r;
}
#endif


#ifdef USE_GUILE
SCM het_group_residues_scm(int imol) {

   SCM r = SCM_EOL;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> specs = graphics_info_t::molecules[imol].het_groups();
      for (unsigned int i=0; i<specs.size(); i++) {
	 SCM s = residue_spec_to_scm(specs[i]);
	 r = scm_cons(s, r);
      }
      r = scm_reverse(r);
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *het_group_residues_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> specs = graphics_info_t::molecules[imol].het_groups();
      r = PyList_New(specs.size());
      for (unsigned int i=0; i<specs.size(); i++) {
         PyList_SetItem(r, i, residue_spec_to_py(specs[i]));
      }
   }
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }

   return r;
}
#endif

/*! \brief return the number of non-hydrogen atoms in the given
  het-group (comp-id).

Return -1 on comp-id not found in dictionary.  */
int het_group_n_atoms(const char *comp_id) {

   graphics_info_t g;
   int r = g.Geom_p()->n_non_hydrogen_atoms(comp_id);
   return r;
}

#ifdef USE_GUILE
SCM new_molecule_sans_biggest_ligand_scm(int imol) {

   SCM r = SCM_BOOL_F;
   std::pair<mmdb::Residue *, int> res = new_molecule_sans_biggest_ligand(imol);
   if (res.first) {
      r = scm_list_2(scm_from_int(res.second), residue_spec_to_scm(coot::residue_spec_t(res.first)));
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *new_molecule_sans_biggest_ligand_py(int imol) {

   PyObject *r = Py_False;
   std::pair<mmdb::Residue *, int> res = new_molecule_sans_biggest_ligand(imol);
   if (res.first) {
      r = PyList_New(2);
      PyList_SetItem(r, 0, PyLong_FromLong(res.second));
      PyList_SetItem(r, 1, residue_spec_to_py(coot::residue_spec_t(res.first)));
   }
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON


std::pair<mmdb::Residue *, int>
new_molecule_sans_biggest_ligand(int imol) {

   mmdb::Residue *res = NULL;
   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {

      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::Residue *r = coot::util::get_biggest_hetgroup(mol);
      if (r) {

	 // copy mol, create a new molecule and delete the residue
	 // with the spec of r from the new molecule.
	 //
	 res = r;
	 mmdb::Manager *new_mol = new mmdb::Manager;
	 new_mol->Copy(mol, mmdb::MMDBFCM_All);
	 atom_selection_container_t asc = make_asc(new_mol);
	 std::string label = "Copy_of_";
	 label += graphics_info_t::molecules[imol].name_;
	 imol_new = graphics_info_t::create_molecule();
	 graphics_info_t g;
	 g.molecules[imol_new].install_model(imol_new, asc, g.Geom_p(), label, 1);
	 update_go_to_atom_window_on_new_mol();
      }
   }
   return std::pair<mmdb::Residue *, int> (res, imol_new);
}



#ifdef USE_GUILE
bool
residues_torsions_match_scm(int imol_1, SCM res_1,
			    int imol_2, SCM res_2,
			    float tolerance) {

   bool r = false;

   graphics_info_t g;
   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {
	 coot::residue_spec_t rs1 = residue_spec_from_scm(res_1);
	 coot::residue_spec_t rs2 = residue_spec_from_scm(res_2);
	 if (!rs1.unset_p() && !rs2.unset_p()) {
	    mmdb::Residue *r_1 = g.molecules[imol_1].get_residue(rs1);
	    mmdb::Residue *r_2 = g.molecules[imol_2].get_residue(rs2);
	    if (r_1 && r_2) {
	       mmdb::Manager *mol1 = g.molecules[imol_1].atom_sel.mol;
	       mmdb::Manager *mol2 = g.molecules[imol_2].atom_sel.mol;
// 	       r = coot::compare_residue_torsions(mol1, r_1,
// 						  mol2, r_2,
// 						  tolerance,
// 						  g.Geom_p());
	    }
	 }
      }
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
bool
residues_torsions_match_py(int imol_1, PyObject *res_1,
			   int imol_2, PyObject *res_2,
			   float tolerance) {

   bool r = false;

   graphics_info_t g;
   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {
	 coot::residue_spec_t rs1 = residue_spec_from_py(res_1);
	 coot::residue_spec_t rs2 = residue_spec_from_py(res_2);
	 if (!rs1.unset_p() && !rs2.unset_p()) {
	    mmdb::Residue *r_1 = g.molecules[imol_1].get_residue(rs1);
	    mmdb::Residue *r_2 = g.molecules[imol_2].get_residue(rs2);
	    if (r_1 && r_2) {
	       mmdb::Manager *mol1 = g.molecules[imol_1].atom_sel.mol;
	       mmdb::Manager *mol2 = g.molecules[imol_2].atom_sel.mol;
// 	       r = coot::compare_residue_torsions(mol1, r_1,
// 						  mol2, r_2,
// 						  tolerance,
// 						  g.Geom_p());
	    }
	 }
      }
   }
   return r;
}
#endif // USE_PYTHON




#include "analysis/kolmogorov.hh"

#ifdef USE_GUILE
double kolmogorov_smirnov_scm(SCM l1, SCM l2) {

   double result = -1;
   if (scm_is_true(scm_list_p(l1)) && scm_is_true(scm_list_p(l2))) {
      SCM length_scm_1 = scm_length(l1);
      SCM length_scm_2 = scm_length(l2);
      unsigned int len_l1 = scm_to_int(length_scm_1);
      unsigned int len_l2 = scm_to_int(length_scm_2);
      std::vector<double> v1;
      std::vector<double> v2;
      for (unsigned int i=0; i<len_l1; i++) {
	 SCM item = scm_list_ref(l1, scm_from_int(i));
	 if (scm_is_true(scm_number_p(item)))
	    v1.push_back(scm_to_double(item));
      }
      for (unsigned int i=0; i<len_l2; i++) {
	 SCM item = scm_list_ref(l2, scm_from_int(i));
	 if (scm_is_true(scm_number_p(item)))
	    v2.push_back(scm_to_double(item));
      }
      result = nicholls::get_KS(v1, v2);
   }
   return result;
}
#endif

#include "analysis/stats.hh"

#ifdef USE_GUILE
double kolmogorov_smirnov_vs_normal_scm(SCM l1, double mean, double std_dev) {

   double result = -1;
   SCM result_scm = SCM_BOOL_F;
   if (scm_is_true(scm_list_p(l1))) {
      SCM length_scm_1 = scm_length(l1);
      unsigned int len_l1 = scm_to_int(length_scm_1);
      std::vector<double> v1;
      for (unsigned int i=0; i<len_l1; i++) {
	 SCM item = scm_list_ref(l1, scm_from_int(i));
	 if (scm_is_true(scm_number_p(item)))
	    v1.push_back(scm_to_double(item));
      }
      result = coot::stats::get_kolmogorov_smirnov_vs_normal(v1, mean, std_dev);
   }
   return result;
}
#endif // USE_GUILE


#ifdef USE_GUILE
SCM kullback_liebler_scm(SCM l1, SCM l2) {

   SCM result_scm = SCM_BOOL_F;
   if (scm_is_true(scm_list_p(l1)) && scm_is_true(scm_list_p(l2))) {
      SCM length_scm_1 = scm_length(l1);
      SCM length_scm_2 = scm_length(l2);
      unsigned int len_l1 = scm_to_int(length_scm_1);
      unsigned int len_l2 = scm_to_int(length_scm_2);
      std::vector<double> v1;
      std::vector<double> v2;
      for (unsigned int i=0; i<len_l1; i++) {
	 SCM item = scm_list_ref(l1, scm_from_int(i));
	 if (scm_is_true(scm_number_p(item)))
	    v1.push_back(scm_to_double(item));
      }
      for (unsigned int i=0; i<len_l2; i++) {
	 SCM item = scm_list_ref(l2, scm_from_int(i));
	 if (scm_is_true(scm_number_p(item)))
	    v2.push_back(scm_to_double(item));
      }
      std::pair<double, double> result = nicholls::get_KL(v1, v2);
      SCM r = scm_list_2(scm_from_double(result.first), scm_from_double(result.second));
   }
   return result_scm;
}
#endif


#ifdef USE_PYTHON
double kolmogorov_smirnov_py(PyObject *l1, PyObject *l2) {

   double result = -1;
   if (PyList_Check(l1) && PyList_Check(l2)) {
      unsigned int len_l1 = PyList_Size(l1);
      unsigned int len_l2 = PyList_Size(l2);
      std::vector<double> v1;
      std::vector<double> v2;
      for (unsigned int i=0; i<len_l1; i++) {
         PyObject *item = PyList_GetItem(l1, i);
         if (PyFloat_Check(item))
            v1.push_back(PyFloat_AsDouble(item));
      }
      for (unsigned int i=0; i<len_l2; i++) {
         PyObject *item = PyList_GetItem(l2, i);
         if (PyFloat_Check(item))
            v2.push_back(PyFloat_AsDouble(item));
      }
      result = nicholls::get_KS(v1, v2);
   }
   return result;
}
#endif

#ifdef USE_PYTHON
double kolmogorov_smirnov_vs_normal_py(PyObject *l1, double mean, double std_dev) {

   double result = -1;

   if (PyList_Check(l1)) {
      unsigned int len_l1 = PyList_Size(l1);
      std::vector<double> v1;
      for (unsigned int i=0; i<len_l1; i++) {
         PyObject *item = PyList_GetItem(l1, i);
         if (PyFloat_Check(item))
            v1.push_back(PyFloat_AsDouble(item));
      }
      result = coot::stats::get_kolmogorov_smirnov_vs_normal(v1, mean, std_dev);
   }
   return result;
}
#endif // USE_PYTHON


#ifdef USE_PYTHON
PyObject *kullback_liebler_py(PyObject *l1, PyObject *l2) {

   PyObject *result_py = Py_False;
   if (PyList_Check(l1) && PyList_Check(l2)) {
      unsigned int len_l1 = PyList_Size(l1);
      unsigned int len_l2 = PyList_Size(l2);
      std::vector<double> v1;
      std::vector<double> v2;
      for (unsigned int i=0; i<len_l1; i++) {
         PyObject *item = PyList_GetItem(l1, i);
         if (PyFloat_Check(item))
            v1.push_back(PyFloat_AsDouble(item));
      }
      for (unsigned int i=0; i<len_l2; i++) {
         PyObject *item = PyList_GetItem(l2, i);
         if (PyFloat_Check(item))
            v2.push_back(PyFloat_AsDouble(item));
      }
      std::pair<double, double> result = nicholls::get_KL(v1, v2);
      result_py = PyList_New(2);
      PyList_SetItem(result_py, 0, PyFloat_FromDouble(result.first));
      PyList_SetItem(result_py, 1, PyFloat_FromDouble(result.second));
   }
   if (PyBool_Check(result_py)) {
      Py_INCREF(result_py);
   }
   return result_py;
}
#endif


// Returning void ATM.  We should return an interesting object at some
// stage. Perhaps a coot::geometry_distortion_info_container_t?
//
double
print_residue_distortions(int imol, std::string chain_id, int res_no, std::string ins_code) {

   double total_distortion = 0.0;

   if (! is_valid_model_molecule(imol)) {
      std::cout << "Not a valid model molecule " << imol << std::endl;
   } else {
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (! residue_p) {
	 std::cout << "Residue not found in molecule " << imol << " "
		   << coot::residue_spec_t(chain_id, res_no, ins_code) << std::endl;
      } else {
	 bool with_nbcs = true;
	 coot::geometry_distortion_info_container_t gdc = g.geometric_distortions(imol, residue_p, with_nbcs);
	 int n_restraints_bonds    = 0;
	 int n_restraints_angles   = 0;
	 int n_restraints_torsions = 0;
	 int n_restraints_chirals  = 0;
	 int n_restraints_planes   = 0;
	 int n_restraints_nbcs     = 0;
	 double sum_penalties_bonds    = 0;
	 double sum_penalties_angles   = 0;
	 double sum_penalties_planes   = 0;
	 double sum_penalties_chirals  = 0;
	 double sum_penalties_torsions = 0;
	 double sum_penalties_nbcs     = 0;
	 std::vector<std::pair<std::string,double> > penalty_string_bonds;
	 std::vector<std::pair<std::string,double> > penalty_string_nbcs;
	 std::vector<std::pair<std::string,double> > penalty_string_angles;
	 std::cout << "Residue Distortion List: \n";
	 for (unsigned int i=0; i<gdc.geometry_distortion.size(); i++) {
	    const coot::simple_restraint &rest = gdc.geometry_distortion[i].restraint;
	    if (rest.restraint_type == coot::BOND_RESTRAINT) {
	       n_restraints_bonds++;
	       mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
	       mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
	       if (at_1 && at_2) {
		  clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		  clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		  double d = sqrt((p2-p1).lengthsq());
		  double distortion = d - rest.target_value;
		  double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
		  std::string s = std::string("bond ")
		     + std::string(at_1->name) + std::string(" to ") + std::string(at_2->name)
		     + std::string(" target_value: ") + coot::util::float_to_string_using_dec_pl(rest.target_value, 3)
		     + std::string(" d: ") + coot::util::float_to_string_using_dec_pl(d, 3)
		     + std::string(" sigma: ") + coot::util::float_to_string_using_dec_pl(rest.sigma, 3)
		     + std::string(" length-devi ") + coot::util::float_to_string_using_dec_pl(distortion, 3)
		     + std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
		  penalty_string_bonds.push_back(std::pair<std::string,double> (s, pen_score));
		  sum_penalties_bonds += pen_score;
	       }
	    }

	    if (rest.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	       n_restraints_bonds++;
	       mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
	       mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
	       if (at_1 && at_2) {
		  clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		  clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		  double d = sqrt((p2-p1).lengthsq());
		  double distortion = d - rest.target_value;
		  double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
		  if (distortion > 0) {
		     pen_score = 0;
		  } else {
		     // don't tell me about NBC that are not too close
		     std::string s = std::string("nbc  ")
			+ std::string(at_1->name) + std::string(" to ") + std::string(at_2->name)
			+ std::string(" target_value: ") + coot::util::float_to_string_using_dec_pl(rest.target_value, 3)
			+ std::string(" d: ") + coot::util::float_to_string_using_dec_pl(d, 3)
			+ std::string(" sigma: ") + coot::util::float_to_string_using_dec_pl(rest.sigma, 3)
			+ std::string(" length-devi ") + coot::util::float_to_string_using_dec_pl(distortion, 3)
			+ std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
		     penalty_string_nbcs.push_back(std::pair<std::string,double> (s, pen_score));
		     sum_penalties_nbcs += pen_score;
		  }
	       }
	    }

	    if (rest.restraint_type == coot::ANGLE_RESTRAINT) {
	       n_restraints_angles++;
	       mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
	       mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
	       mmdb::Atom *at_3 = residue_p->GetAtom(rest.atom_index_3);
	       if (at_1 && at_2 && at_3) {
		  clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		  clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		  clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
		  double angle_rad = clipper::Coord_orth::angle(p1, p2, p3);
		  double angle = clipper::Util::rad2d(angle_rad);
		  double distortion = angle - rest.target_value;
		  double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
		  std::string s = std::string("angle ")
		     + std::string(at_1->name) + std::string(" - ")
		     + std::string(at_2->name) + std::string(" - ")
		     + std::string(at_3->name)
		     + std::string("  target: ") + coot::util::float_to_string(rest.target_value)
		     + std::string(" model_angle: ") + coot::util::float_to_string(angle)
		     + std::string(" sigma: ") + coot::util::float_to_string(rest.sigma)
		     + std::string(" angle-devi ") + coot::util::float_to_string(distortion)
		     + std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
		  penalty_string_angles.push_back(std::pair<std::string,double> (s, pen_score));
		  sum_penalties_angles += pen_score;
	       }
	    }

	    if (rest.restraint_type == coot::TORSION_RESTRAINT) {
	       n_restraints_torsions++;
	       mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
	       mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
	       mmdb::Atom *at_3 = residue_p->GetAtom(rest.atom_index_3);
	       mmdb::Atom *at_4 = residue_p->GetAtom(rest.atom_index_4);
	       if (at_1 && at_2 && at_3 && at_4) {
		  clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		  clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		  clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
		  clipper::Coord_orth p4(at_3->x, at_3->y, at_3->z);
		  double torsion_rad = clipper::Coord_orth::torsion(p1, p2, p3, p4);
		  double torsion = clipper::Util::rad2d(torsion_rad);
		  double distortion = rest.torsion_distortion(torsion);
		  double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
		  std::string s = std::string("torsion ")
		     + std::string(at_1->name) + std::string(" - ")
		     + std::string(at_2->name) + std::string(" - ")
		     + std::string(at_3->name) + std::string(" - ")
		     + std::string(at_4->name)
		     + std::string("  target: ") + coot::util::float_to_string(rest.target_value)
		     + std::string(" model_torsion: ") + coot::util::float_to_string(torsion)
		     + std::string(" sigma: ") + coot::util::float_to_string(rest.sigma)
		     + std::string(" torsion-devi ") + coot::util::float_to_string(distortion)
		     + std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
		  penalty_string_angles.push_back(std::pair<std::string,double> (s, pen_score));
		  sum_penalties_angles += pen_score;
	       }
	    }

	    if (rest.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	       n_restraints_chirals++;
	       double chiral_volume_distortion_limit = 0.1; // arbitrons
	       chiral_volume_distortion_limit = 0;
	       if (gdc.geometry_distortion[i].distortion_score > chiral_volume_distortion_limit) {
		  mmdb::Atom *at_c = residue_p->GetAtom(rest.atom_index_centre);
		  mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
		  mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
		  mmdb::Atom *at_3 = residue_p->GetAtom(rest.atom_index_3);
		  if (at_c && at_1 && at_2 && at_3) {
		     std::cout << "   chiral volume penalty score for chiral vol centred at: "
			       << at_c->name << " with neighbours "
			       << at_1->name << " "
			       << at_2->name << " "
			       << at_3->name << " penalty-score: "
			       << gdc.geometry_distortion[i].distortion_score
			       << std::endl;
		  }
	       }
	       double pen_score = gdc.geometry_distortion[i].distortion_score;
	       sum_penalties_chirals += pen_score;
	    }

	    if (rest.restraint_type == coot::PLANE_RESTRAINT) {
	       n_restraints_planes++;
	       std::vector<mmdb::Atom *> plane_atoms;
	       for (unsigned int iat=0; iat<rest.plane_atom_index.size(); iat++) {
		  mmdb::Atom *at = residue_p->GetAtom(rest.plane_atom_index[iat].first);
		  if (at)
		     plane_atoms.push_back(at);
	       }
	       std::string penalty_string("   plane ");
	       for (unsigned int iat=0; iat<plane_atoms.size(); iat++) {
		  penalty_string += plane_atoms[iat]->name;
		  penalty_string += " ";
	       }
	       std::size_t len = penalty_string.length();
	       if (len < 88) // to match bonds
		  penalty_string += std::string(88-len, ' ');
	       double pen_score = gdc.geometry_distortion[i].distortion_score;
	       sum_penalties_planes += pen_score;
	       penalty_string += std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
	       std::cout << penalty_string << std::endl;
	    }
	 }

	 std::sort(penalty_string_bonds.begin(),  penalty_string_bonds.end(),  coot::util::sd_compare);
	 std::sort(penalty_string_angles.begin(), penalty_string_angles.end(), coot::util::sd_compare);

	 std::reverse(penalty_string_bonds.begin(),  penalty_string_bonds.end());
	 std::reverse(penalty_string_angles.begin(), penalty_string_angles.end());

	 // sorted list, line by line
	 for (unsigned int i=0; i<penalty_string_bonds.size(); i++)
	    std::cout << "   " << penalty_string_bonds[i].first << std::endl;
	 for (unsigned int i=0; i<penalty_string_angles.size(); i++)
	    std::cout << "   " << penalty_string_angles[i].first << std::endl;


	 // Summary:
	 double av_penalty_bond = 0;
	 double av_penalty_angle = 0;
	 double av_penalty_total = 0;
	 if (n_restraints_bonds > 0)
	    av_penalty_bond = sum_penalties_bonds/double(n_restraints_bonds);
	 if (n_restraints_angles > 0)
	    av_penalty_angle = sum_penalties_angles/double(n_restraints_angles);
	 if ((n_restraints_bonds+n_restraints_angles) > 0) {
	    av_penalty_total = (sum_penalties_bonds+sum_penalties_angles)/(n_restraints_bonds+n_restraints_angles);
	 }
	 total_distortion =
	    sum_penalties_bonds  +
	    sum_penalties_angles +
	    sum_penalties_torsions +
	    sum_penalties_chirals  +
	    sum_penalties_planes;

	 std::cout << "Residue Distortion Summary: \n   "
		   << n_restraints_bonds  << " bond restraints\n   "
		   << n_restraints_angles << " angle restraints\n"
		   << "   sum of bond  distortions penalties:  " << sum_penalties_bonds  << "\n"
		   << "   sum of angle distortions penalties:  " << sum_penalties_angles << "\n"
		   << "   average bond  distortion penalty:    " << av_penalty_bond  << "\n"
		   << "   average angle distortion penalty:    " << av_penalty_angle << "\n"
		   << "   total distortion penalty:            " << total_distortion
		   << "\n"
		   << "   average distortion penalty:          " << av_penalty_total
		   << std::endl;
      }
   }
   return total_distortion;

}

// Returning void ATM.  We should return an interesting object at some
// stage. Perhaps a coot::geometry_distortion_info_container_t?
//
void
display_residue_distortions(int imol, std::string chain_id, int res_no, std::string ins_code) {

   if (! is_valid_model_molecule(imol)) {
      std::cout << "Not a valid model molecule " << imol << std::endl;
   } else {
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (! residue_p) {
	 std::cout << "Residue not found in molecule " << imol << " "
		   << coot::residue_spec_t(chain_id, res_no, ins_code) << std::endl;
      } else {
	 bool with_nbcs = true;
	 coot::geometry_distortion_info_container_t gdc = g.geometric_distortions(imol, residue_p, with_nbcs);
         unsigned int n_angles = 0;
	 if (gdc.geometry_distortion.size()) {

	    std::string name = std::string("Ligand Distortion of ");
	    name += residue_p->GetChainID();
	    name += " ";
	    name += coot::util::int_to_string(residue_p->GetSeqNum());
	    name += " ";
	    name += residue_p->GetResName();
            gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
	    int new_obj = new_generic_object_number(name.c_str());
            meshed_generic_display_object &obj = g.generic_display_objects[new_obj];
            obj.attach_to_molecule(imol);
            std::cout << "in display_residue_distortions() with " << gdc.geometry_distortion.size()
                      << " distortions" << std::endl;

	    for (unsigned int i=0; i<gdc.geometry_distortion.size(); i++) {
	       coot::simple_restraint &rest = gdc.geometry_distortion[i].restraint;
	       if (rest.restraint_type == coot::BOND_RESTRAINT) {
		  mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
		  mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
		  if (at_1 && at_2) {
		     clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		     clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		     double d = sqrt((p2-p1).lengthsq());
		     double distortion = d - rest.target_value;
		     double pen_score = fabs(distortion/rest.sigma);
		     coot::colour_holder ch(distortion, 0.1, 5, true, "");
                     ch.scale_intensity(0.5);
                     const unsigned int n_slices = 16;
                     float line_radius = 0.104f;  // more than 0.103 (atoms)
                     std::pair<glm::vec3, glm::vec3> pos_pair(glm::vec3(coord_orth_to_glm(p1)),
                                                              glm::vec3(coord_orth_to_glm(p2)));
                     obj.add_cylinder(pos_pair, ch, line_radius, n_slices, true, true,
                                      meshed_generic_display_object::ROUNDED_CAP,
                                      meshed_generic_display_object::ROUNDED_CAP);
		  }
	       }

	       if (rest.restraint_type == coot::ANGLE_RESTRAINT) {
                  n_angles++;
                  // if (n_angles > 1) continue;
		  mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
		  mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
		  mmdb::Atom *at_3 = residue_p->GetAtom(rest.atom_index_3);
		  if (at_1 && at_2 && at_3) {
		     clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		     clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		     clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
		     double angle_rad = clipper::Coord_orth::angle(p1, p2, p3);
		     double angle = clipper::Util::rad2d(angle_rad);
		     double distortion = fabs(angle - rest.target_value);
		     coot::colour_holder ch(distortion, 0.1, 5, true, "");
                     ch.scale_intensity(0.6);

		     try {
			float radius = 0.66;
			float radius_inner = 0.1; // should match bond line_radius
			coot::arc_info_type arc_angle_info(at_1, at_2, at_3);
                        meshed_generic_display_object::arc_t arc(arc_angle_info, radius, radius_inner, ch);
                        obj.add_arc(arc);
		     }
		     catch (const std::runtime_error &rte) {
			std::cout << "WARNING:: " << rte.what() << std::endl;
		     }
		  }
	       }
	       if (rest.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
		  mmdb::Atom *at_c = residue_p->GetAtom(rest.atom_index_centre);
		  mmdb::Atom *at_1 = residue_p->GetAtom(rest.atom_index_1);
		  mmdb::Atom *at_2 = residue_p->GetAtom(rest.atom_index_2);
		  mmdb::Atom *at_3 = residue_p->GetAtom(rest.atom_index_3);
		  if (at_c && at_1 && at_2 && at_3) {
		     clipper::Coord_orth pc(at_c->x, at_c->y, at_c->z);
		     clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
		     clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
		     clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
		     clipper::Coord_orth bl_1 = 0.6 * pc + 0.4 * p1;
		     clipper::Coord_orth bl_2 = 0.6 * pc + 0.4 * p2;
		     clipper::Coord_orth bl_3 = 0.6 * pc + 0.4 * p3;
		     double distortion = sqrt(fabs(gdc.geometry_distortion[i].distortion_score));
		     coot::colour_holder ch(distortion, 0.1, 5, true, "");
                     ch.scale_intensity(0.5);
		     // to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
		     //    			bl_1.x(), bl_1.y(), bl_1.z(),
		     //    			bl_2.x(), bl_2.y(), bl_2.z());
		     // to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
		     //    			bl_1.x(), bl_1.y(), bl_1.z(),
		     //    			bl_3.x(), bl_3.y(), bl_3.z());
		     // to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
		     //    			bl_2.x(), bl_2.y(), bl_2.z(),
		     //    			bl_3.x(), bl_3.y(), bl_3.z());
                     float line_radius = 0.11f;
                     const unsigned int n_slices = 16;

                     std::pair<glm::vec3, glm::vec3> pos_pair(glm::vec3(coord_orth_to_glm(bl_1)),
                                                              glm::vec3(coord_orth_to_glm(bl_2)));

                     std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > v =
                        { std::make_pair(bl_1, bl_2), std::make_pair(bl_1, bl_2), std::make_pair(bl_2, bl_3) };

                     std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> >::const_iterator it;
                     for (it=v.begin(); it!=v.end(); ++it)
                        obj.add_cylinder(pos_pair, ch, line_radius, n_slices, true, true,
                                         meshed_generic_display_object::ROUNDED_CAP,
                                         meshed_generic_display_object::ROUNDED_CAP);

		     // return (if possible) the atom attached to
		     // at_c that is not at_1, at_2 or at_3.
		     mmdb::Atom *at_4th = coot::chiral_4th_atom(residue_p, at_c, at_1, at_2, at_3);
		     if (at_4th) {
			std::cout << "    " << coot::atom_spec_t(at_4th) << std::endl;
			clipper::Coord_orth p4(at_4th->x, at_4th->y, at_4th->z);
			clipper::Coord_orth bl_4 = 0.6 * pc + 0.4 * p4;
			// to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
			// 			   bl_1.x(), bl_1.y(), bl_1.z(),
			// 			   bl_4.x(), bl_4.y(), bl_4.z());
			// to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
			// 			   bl_2.x(), bl_2.y(), bl_2.z(),
			// 			   bl_4.x(), bl_4.y(), bl_4.z());
			// to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
			// 			   bl_3.x(), bl_3.y(), bl_3.z(),
			// 			   bl_4.x(), bl_4.y(), bl_4.z());

                        std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > v4 =
                           { std::make_pair(bl_1, bl_4), std::make_pair(bl_2, bl_4), std::make_pair(bl_3, bl_4) };

                        for (it=v4.begin(); it!=v4.end(); ++it)
                           obj.add_cylinder(pos_pair, ch, line_radius, n_slices, true, true,
                                            meshed_generic_display_object::ROUNDED_CAP,
                                            meshed_generic_display_object::ROUNDED_CAP);

		     } else {
			// make 4th tetrahedron point from the others
			clipper::Coord_orth neighb_sum = p1 + p2 + p3;
			clipper::Coord_orth neighb_average = 0.33333333 * neighb_sum;
			clipper::Coord_orth dir_unit(clipper::Coord_orth(pc - neighb_average).unit());
			clipper::Coord_orth p4(pc + 1.2 * dir_unit);
			clipper::Coord_orth bl_4 = 0.6 * pc + 0.4 * p4;
			// to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
			// 			   bl_1.x(), bl_1.y(), bl_1.z(),
			// 			   bl_4.x(), bl_4.y(), bl_4.z());
			// to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
			// 			   bl_2.x(), bl_2.y(), bl_2.z(),
			// 			   bl_4.x(), bl_4.y(), bl_4.z());
			// to_generic_object_add_line(new_obj, ch.hex().c_str(), 2,
			// 			   bl_3.x(), bl_3.y(), bl_3.z(),
			// 			   bl_4.x(), bl_4.y(), bl_4.z());

                        std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > v4 =
                           { std::make_pair(bl_1, bl_4), std::make_pair(bl_2, bl_4), std::make_pair(bl_3, bl_4) };

                        for (it=v4.begin(); it!=v4.end(); ++it)
                           obj.add_cylinder(pos_pair, ch, line_radius, n_slices, true, true,
                                            meshed_generic_display_object::ROUNDED_CAP,
                                            meshed_generic_display_object::ROUNDED_CAP);
		     }
		  }
	       }
	    }
            // attach_buffers() here or upstairs
            Material material;
            material.do_specularity = true;
            material.specular_strength = 0.9;
            obj.mesh.setup(material);
	    set_display_generic_object(new_obj, 1);
	    graphics_draw();
	 }
      }
   }
}





void
write_dictionary_from_residue(int imol, std::string chain_id, int res_no, std::string ins_code, std::string cif_file_name) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (! residue_p) {
	 std::cout << "Residue not found in molecule " << imol << " "
		   << coot::residue_spec_t(chain_id, res_no, ins_code) << std::endl;
      } else {
	 mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
	 if (mol) {
	    coot::dictionary_residue_restraints_t d(mol);
	    d.write_cif(cif_file_name);
	 }
	 delete mol;
      }
   }
}

void
add_dictionary_from_residue(int imol, std::string chain_id, int res_no, std::string ins_code) {

   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (! residue_p) {
	 std::cout << "Residue not found in molecule " << imol << " "
		   << coot::residue_spec_t(chain_id, res_no, ins_code) << std::endl;
      } else {
	 mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
	 if (mol) {
	    coot::dictionary_residue_restraints_t d(mol);
	    std::cout << "INFO:: replacing restraints for type \""
		      << d.residue_info.comp_id << "\"" << std::endl;
	    int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
	    g.Geom_p()->replace_monomer_restraints(d.residue_info.comp_id, imol_enc, d);

	    if (0) {
	       std::pair<bool, coot::dictionary_residue_restraints_t>
		  r = g.Geom_p()->get_monomer_restraints(d.residue_info.comp_id, imol);
	       if (! r.first) {
		  std::cout << "-------------------- problem retrieving restraints " << std::endl;
	       } else {
		  std::cout << "-------------------- got restraints " << std::endl;
		  for (unsigned int ib=0; ib<r.second.bond_restraint.size(); ib++) {
		     const coot::dict_bond_restraint_t &rest = r.second.bond_restraint[ib];
		     std::cout << ib << "   " << rest.atom_id_1_4c() << " " << rest.atom_id_2_4c() << " "
			       << rest.value_dist() << std::endl;
		  }
	       }
	    }
	 }
	 delete mol;
      }
   }
}


void display_residue_hydrogen_bond_atom_status_using_dictionary(int imol, std::string chain_id, int res_no,
								std::string ins_code) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (! residue_p) {
	 std::cout << "Residue not found in molecule " << imol << " "
		   << coot::residue_spec_t(chain_id, res_no, ins_code) << std::endl;
      } else {
	 coot::h_bonds hb;
	 mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	 int SelHnd_lig = mol->NewSelection();
	 mol->SelectAtoms(SelHnd_lig, 0, chain_id.c_str(),
			  res_no, ins_code.c_str(),
			  res_no, ins_code.c_str(),
			  "*", "*", "*", "*");

	 std::pair<bool, int> status = hb.check_hb_status(SelHnd_lig, mol, *g.Geom_p(), imol);
	 if (! status.first) {
	    std::cout << "WARNING:: ===================== no HB status on atoms of ligand! ======="
		      << "=========" << std::endl;
	 } else {
	    // happy path
	    std::string name = "HB Acceptor/Donor/Both/H for ";
	    name += residue_p->GetChainID();
	    name += " ";
	    name += coot::util::int_to_string(residue_p->GetSeqNum());
	    name += " ";
	    name += residue_p->GetInsCode();
	    meshed_generic_display_object features_obj;
	    mmdb::PPAtom residue_atoms = 0;
	    int n_residue_atoms;
	    mol->GetSelIndex(SelHnd_lig, residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       mmdb::Atom *at = residue_atoms[iat];
	       int hb_type = coot::HB_UNASSIGNED;
	       at->GetUDData(status.second, hb_type);
	       if (hb_type != coot::HB_UNASSIGNED) {
		  clipper::Coord_orth centre = coot::co(at);
		  meshed_generic_display_object::sphere_t sphere(centre, 0.5);
		  if (hb_type == coot::HB_DONOR) {
		     sphere.col = glm::vec4(0.2, 0.6, 0.7, 1.0);
		  }
		  if (hb_type == coot::HB_ACCEPTOR) {
		     sphere.col = glm::vec4(0.8, 0.2, 0.2, 1.0);
		  }
		  if (hb_type == coot::HB_BOTH) {
		     sphere.col = glm::vec4(0.8, 0.2, 0.8, 1.0);
		  }
		  if (hb_type == coot::HB_HYDROGEN) {
		     sphere.radius = 0.35;
		  }
		  if (hb_type == coot::HB_DONOR    ||
		      hb_type == coot::HB_ACCEPTOR ||
		      hb_type == coot::HB_BOTH     ||
		      hb_type == coot::HB_HYDROGEN) {
		     features_obj.add_sphere(sphere);
		  }
	       }
	    }
	    features_obj.mesh.set_draw_this_mesh(true);
	    g.generic_display_objects.push_back(features_obj);
	    graphics_draw();
	 }
	 mol->DeleteSelection(SelHnd_lig);
      }
   }
}

// create a new dictionary entry and molecule.
void
invert_chiral_centre(int imol,
		     std::string chain_id,
		     int res_no,
		     std::string ins_code,
		     std::string atom_name) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::pair<mmdb::Residue *, coot::dictionary_residue_restraints_t> rr =
	 g.molecules[imol].invert_chiral_centre(chain_id, res_no, ins_code, atom_name,
						*g.Geom_p());
      if (rr.first) {
	 // add new restraints to dictionary
	 std::string new_type = rr.second.residue_info.comp_id;
	 // THIS ONE IS COMPLICATED, FIXME
	 g.Geom_p()->replace_monomer_restraints(new_type, imol, rr.second);
	 // add a new molecule from this residue
	 mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(rr.first);
	 // I think we can delete residue now that we have made mol.
	 delete rr.first;
	 rr.first = NULL;
	 int g_mol = g.create_molecule();
	 atom_selection_container_t asc = make_asc(mol);
	 std::string label = "New Residue " + new_type;
	 g.molecules[g_mol].install_model(g_mol, asc, g.Geom_p(), label, 1);
	 graphics_draw();
      }
   }
}

bool
comprised_of_organic_set_p(const std::string &comp_id) {

   bool r = false;
   graphics_info_t g;
   int imol = 0; // dummy
   if (g.Geom_p()->have_dictionary_for_residue_type(comp_id, imol, g.cif_dictionary_read_number++)) {
      int imol = 0; // dummy;
      std::pair<bool, coot::dictionary_residue_restraints_t> rp =
	 g.Geom_p()->get_monomer_restraints(comp_id, imol);
      if (rp.first) {
         r = rp.second.comprised_of_organic_set();
      }
   }
   return r;
}

bool
enhanced_ligand_coot_p() {

   bool r = false;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   r = true;
#endif

   return r;
}



#ifdef USE_SQLITE3
#include <stdexcept>
#include <sqlite3.h>
#include "ligand/ligands-db.hh"
#endif // USE_SQLITE3

// we want to read in the built-in database to convert these scores to percentiles
// return -1 (test for negative) on failure
// metric_name examples: density_correlation coot_diff_map_KS_2 mogul_z_worst bumps_1
// range 0->100 or -1
//
// reverse order for mogul and bumps (for example) because low-is-good
//
double get_ligand_percentile(std::string metric_name, double metric_value, short int reverse_order) {

   double pc = -1;
#ifdef USE_SQLITE3
   std::string database_name = coot::package_data_dir() + "/data/ligands-2016.db";
   bool low_is_good = reverse_order;

   coot::ligand_metrics lm(database_name);
   std::pair<int, int> idx_pair = lm.get_index(metric_value, metric_name, low_is_good);
   if (idx_pair.second != 0) {
      pc = 100.0 * double(idx_pair.first)/double(idx_pair.second);
      std::cout << metric_name << "DEBUG:: percentile for " << metric_value << " is " << pc << "%"
		<< std::endl;
   }
#endif // USE_SQLITE3
   return pc;
}


// ---------------------------------------------------


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//                 coot built-in contact dots
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

#include "coot-utils/atom-overlaps.hh"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>
// #include <glm/ext.hpp> // 20240326-PE
#include <time.h>
#include "glarea_tick_function.hh"


void
coot_contact_dots_for_ligand_instancing_version(int imol, const coot::residue_spec_t &res_spec) {


#if 0 // 20210911-PE - goodbye ball clashes.
   auto setup_ball_clashes = [] () {

                                      std::vector<Instanced_Markup_Mesh_attrib_t> balls;
                                      balls.resize(c.clashes.size());
                                      std::string mesh_name = molecule_name_stub;
                                      mesh_name += " eeek clashes";
                                      Instanced_Markup_Mesh im_in(mesh_name);
                                      g.molecules[imol].instanced_meshes.push_back(im_in);
                                      Instanced_Markup_Mesh &clash_im = g.molecules[imol].instanced_meshes.back();
                                      clash_im.setup_octasphere(octasphere_subdivisions);
                                      clash_im.setup_instancing_buffers(c.clashes.size());
                                      const float point_size = ball_size;
                                      for (unsigned int i=0; i<c.clashes.size(); i++) {
                                         glm::vec3 position(c.clashes[i].first.x(), c.clashes[i].first.y(), c.clashes[i].first.z());
                                         Instanced_Markup_Mesh_attrib_t attribs(clash_col_glm, position, point_size);
                                         balls[i] = attribs;
                                         std::pair<glm::vec3, glm::vec3> pos_pair(glm::vec3(coord_orth_to_glm(c.clashes[i].first)),
                                                                                  glm::vec3(coord_orth_to_glm(c.clashes[i].second)));
                                      }
                                      clash_im.update_instancing_buffers(balls);
                                   };
#endif


   graphics_info_t g;
   unsigned int octasphere_subdivisions = 3; // make a member of graphics_info_t with an API

   mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
   mmdb::Residue *residue_p = coot::util::get_residue(res_spec, mol);
   if (!residue_p) {
      std::cout << "Can't find residue" << res_spec << std::endl;
   } else {

      // I don't like this part
      std::map<std::string, coot::colour_holder> colour_map;
      colour_map["blue"      ] = ::colour_values_from_colour_name("blue");
      colour_map["sky"       ] = ::colour_values_from_colour_name("sky");
      colour_map["sea"       ] = ::colour_values_from_colour_name("sea");
      colour_map["greentint" ] = ::colour_values_from_colour_name("greentint");
      colour_map["darkpurple"] = ::colour_values_from_colour_name("darkpurple");
      colour_map["green"     ] = ::colour_values_from_colour_name("green");
      colour_map["orange"    ] = ::colour_values_from_colour_name("orange");
      colour_map["orangered" ] = ::colour_values_from_colour_name("orangered");
      colour_map["yellow"    ] = ::colour_values_from_colour_name("yellow");
      colour_map["yellowtint"] = ::colour_values_from_colour_name("yellowtint");
      colour_map["red"       ] = ::colour_values_from_colour_name("red");
      colour_map["#55dd55"   ] = ::colour_values_from_colour_name("#55dd55");
      colour_map["hotpink"   ] = ::colour_values_from_colour_name("hotpink");
      colour_map["grey"      ] = ::colour_values_from_colour_name("grey");
      colour_map["magenta"   ] = ::colour_values_from_colour_name("magenta");
      colour_map["royalblue" ] = ::colour_values_from_colour_name("royalblue");

      auto colour_string_to_colour_holder = [colour_map] (const std::string &c) {
                                               return colour_map.find(c)->second;
                                            };
      std::vector<mmdb::Residue *> neighbs = coot::residues_near_residue(residue_p, mol, 5);
      coot::atom_overlaps_container_t overlaps(residue_p, neighbs, mol, imol, g.Geom_p(), 0.5, 0.25);
      float cdd = graphics_info_t::contact_dots_density;

      coot::atom_overlaps_dots_container_t c = overlaps.contact_dots_for_ligand(cdd * 2.0f);

      gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
      std::string molecule_name_stub = "Molecule " + coot::util::int_to_string(imol) + ": Ligand Contact Dots ";

      bool extra_annotation = graphics_info_t::get_exta_annotation_state();

      float ball_size = 0.07; // about the right size
      ball_size = 0.06; // 20211129-PE try this
      float tube_radius = ball_size;
      g.setup_cylinder_clashes(c, imol, tube_radius, extra_annotation);

      std::unordered_map<std::string, std::vector<coot::atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
      for (it=c.dots.begin(); it!=c.dots.end(); ++it) {
         float specular_strength = 0.5; //  default
         const std::string &type = it->first;
         // std::cout << "::::::::::: type  " << type << std::endl;
         // if (type != "big-overlap") continue;
         const std::vector<coot::atom_overlaps_dots_container_t::dot_t> &v = it->second;
         float point_size = ball_size;
         // if (type == "vdw-surface") point_size = 0.05;
         if (type == "vdw-surface") point_size = 0.03;
         if (type == "vdw-surface") specular_strength= 0.1; // dull, reduces zoomed out speckles
         std::string mesh_name = molecule_name_stub + type;
         Instanced_Markup_Mesh im_in(mesh_name);
         g.molecules[imol].instanced_meshes.push_back(im_in);
         Instanced_Markup_Mesh &im = g.molecules[imol].instanced_meshes.back();

         im.setup_octasphere(octasphere_subdivisions);
         im.setup_instancing_buffers(v.size());
         std::vector<Instanced_Markup_Mesh_attrib_t> balls;
         balls.resize(v.size());
         for (unsigned int i=0; i<v.size(); i++) {
            coot::colour_holder ch = colour_string_to_colour_holder(v[i].col);
            // std::cout << "colour convert " << v[i].col << " " << ch << std::endl;
            glm::vec4 colour(ch.red, ch.green, ch.blue, 1.0);
            glm::vec3 position(v[i].pos.x(), v[i].pos.y(), v[i].pos.z());
            Instanced_Markup_Mesh_attrib_t attribs(colour, position, point_size);
            attribs.specular_strength = specular_strength;
            balls[i] = attribs;
         }
         // std::cout << ":::::: im.name " << im.get_name() << std::endl;
         im.update_instancing_buffers(balls);
      }
   }
   graphics_draw();
}

void
coot_contact_dots_for_ligand_internal(int imol, const coot::residue_spec_t &res_spec) {

   graphics_info_t g;
   mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
   mmdb::Residue *residue_p = coot::util::get_residue(res_spec, mol);
   if (!residue_p) {
      std::cout << "Can't find residue" << res_spec << std::endl;
   } else {
      std::vector<mmdb::Residue *> neighbs = coot::residues_near_residue(residue_p, mol, 5);
      coot::atom_overlaps_container_t overlaps(residue_p, neighbs, mol, imol, g.Geom_p(), 0.5, 0.25);
      coot::atom_overlaps_dots_container_t c = overlaps.contact_dots_for_ligand();
      std::cout << "------------- coot_contact_dots_for_ligand_internal(): score " << c.score() << std::endl;

      // for quick colour lookups.
      std::map<std::string, coot::colour_holder> colour_map;
      colour_map["blue"      ] = coot::old_generic_display_object_t::colour_values_from_colour_name("blue");
      colour_map["sky"       ] = coot::old_generic_display_object_t::colour_values_from_colour_name("sky");
      colour_map["sea"       ] = coot::old_generic_display_object_t::colour_values_from_colour_name("sea");
      colour_map["greentint" ] = coot::old_generic_display_object_t::colour_values_from_colour_name("greentint");
      colour_map["darkpurple"] = coot::old_generic_display_object_t::colour_values_from_colour_name("darkpurple");
      colour_map["green"     ] = coot::old_generic_display_object_t::colour_values_from_colour_name("green");
      colour_map["orange"    ] = coot::old_generic_display_object_t::colour_values_from_colour_name("orange");
      colour_map["orangered" ] = coot::old_generic_display_object_t::colour_values_from_colour_name("orangered");
      colour_map["yellow"    ] = coot::old_generic_display_object_t::colour_values_from_colour_name("yellow");
      colour_map["yellowtint"] = coot::old_generic_display_object_t::colour_values_from_colour_name("yellowtint");
      colour_map["red"       ] = coot::old_generic_display_object_t::colour_values_from_colour_name("red");
      colour_map["#55dd55"   ] = coot::old_generic_display_object_t::colour_values_from_colour_name("#55dd55");
      colour_map["hotpink"   ] = coot::old_generic_display_object_t::colour_values_from_colour_name("hotpink");
      colour_map["grey"      ] = coot::old_generic_display_object_t::colour_values_from_colour_name("grey");
      colour_map["magenta"   ] = coot::old_generic_display_object_t::colour_values_from_colour_name("magenta");
      colour_map["royalblue" ] = coot::old_generic_display_object_t::colour_values_from_colour_name("royalblue");

      gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
      Material material;
      material.specular_strength *= 0.5;
      std::unordered_map<std::string, std::vector<coot::atom_overlaps_dots_container_t::dot_t> >::const_iterator it;
      for (it=c.dots.begin(); it!=c.dots.end(); ++it) {
	 const std::string &type = it->first;
	 const std::vector<coot::atom_overlaps_dots_container_t::dot_t> &v = it->second;
	 std::string obj_name = "Molecule ";
	 obj_name += coot::util::int_to_string(imol) + ": " + type;
	 int obj_index = g.generic_object_index(obj_name);
         if (obj_index == -1)
            obj_index = new_generic_object_number_for_molecule(obj_name, imol);
         else
            g.generic_display_objects[obj_index].clear();
	 int point_size = 2;
	 for (unsigned int i=0; i<v.size(); i++) {
	    const std::string &colour_string = v[i].col;
            coot::colour_holder ch = colour_map[colour_string];
            if (type == "vdw-surface") {
               point_size = 1;
               if (false) {
                  ch = colour_values_from_colour_name("#ffc3a0");
                  ch.scale_intensity(0.5);
                  material.shininess = 5.0;
                  material.specular_strength = 0.2;
               }
            }
	    to_generic_object_add_point_internal(obj_index, colour_string, ch, point_size, v[i].pos);
	 }
         // now setup() that mesh
         // g.generic_display_objects[obj_index].mesh.setup(&g.shader_for_moleculestotriangles, material); // 20210910-PE
         g.generic_display_objects[obj_index].mesh.setup(material);
	 if (type != "vdw-surface")
	    set_display_generic_object_simple(obj_index, 1); // a function with no redraw
      }
      std::string clashes_name = "Molecule " + coot::util::int_to_string(imol) + ":";
      clashes_name += " eeeee ligand_internal clashes";
      int clashes_obj_index = g.generic_object_index(clashes_name);
      if (clashes_obj_index == -1)
         clashes_obj_index = new_generic_object_number_for_molecule(clashes_name, imol);
      else
         g.generic_display_objects[clashes_obj_index].clear();
      meshed_generic_display_object &obj = g.generic_display_objects[clashes_obj_index];
      coot::colour_holder clash_col = colour_values_from_colour_name("#ff59b4");
      float line_radius = 0.062f;
      const unsigned int n_slices = 16;
      for (unsigned int i=0; i<c.clashes.size(); i++) {
         std::pair<glm::vec3, glm::vec3> pos_pair(glm::vec3(coord_orth_to_glm(c.clashes[i].first)),
                                                  glm::vec3(coord_orth_to_glm(c.clashes[i].second)));
         obj.add_cylinder(pos_pair, clash_col, line_radius, n_slices, true, true,
                          meshed_generic_display_object::ROUNDED_CAP,
                          meshed_generic_display_object::ROUNDED_CAP);
      }
      // obj.mesh.setup(&g.shader_for_moleculestotriangles, material); 20210910-PE
      obj.mesh.setup(material);
      set_display_generic_object(clashes_obj_index, 1);
   }
}

#ifdef USE_PYTHON
void
coot_contact_dots_for_ligand_py(int imol, PyObject *ligand_spec_py) {

   coot::residue_spec_t res_spec = residue_spec_from_py(ligand_spec_py);
   if (is_valid_model_molecule(imol)) {
      coot_contact_dots_for_ligand_instancing_version(imol, res_spec);
   }
}
#endif

#ifdef USE_PYTHON
void switch_HIS_protonation_py(int imol, PyObject *residue_spec_py) {

   coot::residue_spec_t res_spec = residue_spec_from_py(residue_spec_py);
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].switch_HIS_protonation(res_spec);
   }
   graphics_draw();
}
#endif // USE_PYTHON

#ifdef USE_GUILE
void switch_HIS_protonation_scm(int imol, SCM residue_spec_scm) {

   coot::residue_spec_t res_spec = residue_spec_from_scm(residue_spec_scm);
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].switch_HIS_protonation(res_spec);
   }
   graphics_draw();
}
#endif // USE_GUILE


void coot_add_hydrogen_atoms(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].reduce(g.Geom_p());
      graphics_draw();
   }
}

void add_hydrogen_atoms(int imol) {
   coot_add_hydrogen_atoms(imol);
}

void coot_reduce(int imol) {
   coot_add_hydrogen_atoms(imol);
}



#ifdef USE_GUILE
void
coot_contact_dots_for_ligand_scm(int imol, SCM ligand_spec_scm) {

   coot::residue_spec_t res_spec = residue_spec_from_scm(ligand_spec_scm);
   if (is_valid_model_molecule(imol)) {
      coot_contact_dots_for_ligand_instancing_version(imol, res_spec);
   }
}
#endif

//! \brief set if all atom contact should ignore water-water interactions (default off)
void set_all_atom_contact_dots_ignore_water(short int state) {

   graphics_info_t::all_atom_contact_dots_ignore_water_flag = state;

}

//! \brief set if all atom contact should make the vdw surface (default off/false)
void set_all_atom_contact_dots_do_vdw_surface(short int state) {

   graphics_info_t::all_atom_contact_dots_do_vdw_surface = state; // short int to bool conversion
}



// all-atom contact dots.  This is not the place for this definition (not a ligand function)
//
void coot_all_atom_contact_dots_old(int imol) {

   std::cout << "coot_all_atom_contact_dots_old()" << std::endl;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      // spike-length ball-radius
      bool ignore_waters = g.all_atom_contact_dots_ignore_water_flag;
      coot::atom_overlaps_container_t overlaps(mol, g.Geom_p(), ignore_waters, 0.5, 0.25); // spike length, ball radius
      // dot density
      float dd = graphics_info_t::contact_dots_density;
      coot::atom_overlaps_dots_container_t c = overlaps.all_atom_contact_dots(dd, true);

      std::unordered_map<std::string, std::vector<coot::atom_overlaps_dots_container_t::dot_t> >::const_iterator it;

      // for quick colour lookups.
      std::map<std::string, coot::colour_holder> colour_map;
      colour_map["blue"      ] = colour_values_from_colour_name("blue");
      colour_map["sky"       ] = colour_values_from_colour_name("sky");
      colour_map["sea"       ] = colour_values_from_colour_name("sea");
      colour_map["greentint" ] = colour_values_from_colour_name("greentint");
      colour_map["green"     ] = colour_values_from_colour_name("green");
      colour_map["darkpurple"] = colour_values_from_colour_name("darkpurple");
      colour_map["orange"    ] = colour_values_from_colour_name("orange");
      colour_map["orangered" ] = colour_values_from_colour_name("orangered");
      colour_map["yellow"    ] = colour_values_from_colour_name("yellow");
      colour_map["yellowtint"] = colour_values_from_colour_name("yellowtint");
      colour_map["red"       ] = colour_values_from_colour_name("red");
      colour_map["#55dd55"   ] = colour_values_from_colour_name("#55dd55");
      colour_map["hotpink"   ] = colour_values_from_colour_name("hotpink");
      colour_map["grey"      ] = colour_values_from_colour_name("grey");
      colour_map["magenta"   ] = colour_values_from_colour_name("magenta");

      Material material;
      material.specular_strength *= 0.4;
      gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
      for (it=c.dots.begin(); it!=c.dots.end(); ++it) {
	 const std::string &type = it->first;
	 const std::vector<coot::atom_overlaps_dots_container_t::dot_t> &v = it->second;
	 std::string obj_name = "Molecule ";
	 obj_name += coot::util::int_to_string(imol) + ": " + type;
	 int obj_index = g.generic_object_index(obj_name);
         if (obj_index == -1)
            obj_index = new_generic_object_number_for_molecule(obj_name, imol);
         else
            g.generic_display_objects[obj_index].clear();
	 std::string col = "#445566";
	 int point_size = 2;
	 if (type == "vdw-surface") point_size = 1;
	 if (type == "vdw-surface") continue; // for testing - if this line is included - catatonia on proteins
	 for (unsigned int i=0; i<v.size(); i++) {
	    const std::string &colour_string = v[i].col;
            coot::colour_holder ch = colour_map[colour_string];
	    to_generic_object_add_point_internal(obj_index, colour_string, ch, point_size, v[i].pos);
	 }
         // g.generic_display_objects[obj_index].mesh.setup(&g.shader_for_moleculestotriangles, material);
         g.generic_display_objects[obj_index].mesh.setup(material);
	 if (type != "vdw-surface")
	    set_display_generic_object_simple(obj_index, 1); // should be a function with no redraw
      }
      std::string clashes_name = "Molecule " + coot::util::int_to_string(imol) + ":";
      clashes_name += " clashes";
      int clashes_obj_index = g.generic_object_index(clashes_name);
      if (clashes_obj_index == -1)
         clashes_obj_index = new_generic_object_number_for_molecule(clashes_name, imol);
      else
         g.generic_display_objects[clashes_obj_index].clear();
      meshed_generic_display_object &obj = g.generic_display_objects[clashes_obj_index];
      coot::colour_holder clash_col = colour_values_from_colour_name("#ff59b4");
      float line_radius = 0.062f;
      const unsigned int n_slices = 16;
      for (unsigned int i=0; i<c.clashes.size(); i++) {
         obj.add_cylinder(std::pair<glm::vec3, glm::vec3>(glm::vec3(coord_orth_to_glm(c.clashes[i].first)),
                                                          glm::vec3(coord_orth_to_glm(c.clashes[i].second))),
                          clash_col, line_radius, n_slices, true, true,
                          meshed_generic_display_object::ROUNDED_CAP,
                          meshed_generic_display_object::ROUNDED_CAP);

      }
      // obj.mesh.setup(&g.shader_for_moleculestotriangles, material); 20210910-PE
      obj.mesh.setup(material);
      set_display_generic_object_simple(clashes_obj_index, 1);
      graphics_draw();
   }
}

void coot_all_atom_contact_dots_instanced(mmdb::Manager *mol, int imol) {

   unsigned int octasphere_subdivisions = 1; // make a member of graphics_info_t with an API

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.coot_all_atom_contact_dots_instanced(mol, imol);
   }

}

void set_contact_dots_density(float density) {

   graphics_info_t::contact_dots_density = density;
}

//! \brief set the number of subdivisions of contact dot density sphere (1=low (default), 3=high)
void set_contact_dot_sphere_n_subdivisions(unsigned int n_subdivisions) {

   graphics_info_t::contact_dot_sphere_subdivisions = n_subdivisions;
}



void coot_all_atom_contact_dots_instanced(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      coot_all_atom_contact_dots_instanced(mol, imol);
   }
   graphics_draw();
}

void coot_all_atom_contact_dots(int imol) {

   coot_all_atom_contact_dots_instanced(imol);
}

//! \brief as above, but use the active molecule
void coot_all_atom_contact_dots_active_molecule() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot_all_atom_contact_dots_instanced(imol);
   }
}


#ifdef USE_GUILE
SCM linked_residues_scm(SCM residue_centre_scm, int imol, float close_dist_max) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {

      coot::residue_spec_t residue_spec = residue_spec_from_scm(residue_centre_scm);
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(residue_spec);
      if (residue_p) {
         std::vector<mmdb::Residue *> v = coot::simple_residue_tree(residue_p, mol, close_dist_max);
         r = SCM_EOL;
         for (std::size_t i=0; i<v.size(); i++) {
            mmdb::Residue *residue_p = v[i];
            if (residue_p) {
               coot::residue_spec_t spec(residue_p);
               SCM item_scm = residue_spec_to_scm(spec);
               r = scm_cons(item_scm, r);
            }
         }
      }
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *linked_residues_py(PyObject *residue_centre_py, int imol, float close_dist_max) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {

      coot::residue_spec_t residue_spec = residue_spec_from_py(residue_centre_py);
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(residue_spec);
      if (residue_p) {
         std::vector<mmdb::Residue *> v = coot::simple_residue_tree(residue_p, mol, close_dist_max);
         r = PyList_New(v.size());
         for (std::size_t i=0; i<v.size(); i++) {
            mmdb::Residue *residue_p = v[i];
            if (residue_p) {
               coot::residue_spec_t spec(residue_p);
               PyObject *item = residue_spec_to_py(spec);
               PyList_SetItem(r, i, item);
            } else {
               // or put the valid residue_p in a vector first
               PyList_SetItem(r, i, Py_False); // Hmm..
            }
         }
      }
   }

   if (PyBool_Check(r))
      Py_INCREF(r);
   return r;

}
#endif



/*! \brief make chain masked maps

   needs to return a list of values
 */
int make_masked_maps_split_by_chain(int imol, int imol_map) {

   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         graphics_info_t g;
         coot::ligand lig;
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         lig.set_map_atom_mask_radius(3.3);
         lig.import_map_from(g.molecules[imol_map].xmap);
         float contour_level = g.molecules[imol_map].get_contour_level();
         std::vector<std::pair<std::string, clipper::Xmap<float> > > maps = lig.make_masked_maps_split_by_chain(mol);
         std::cout << "INFO:: made " << maps.size() << " masked maps" << std::endl;
         bool is_em_flag = g.molecules[imol_map].is_EM_map();
         for(unsigned int i=0; i<maps.size(); i++) {
            std::string map_name =maps[i].first;
            int imol_new_map = graphics_info_t::create_molecule();
            g.molecules[imol_new_map].install_new_map(maps[i].second, map_name, is_em_flag);
            g.molecules[imol_new_map].set_contour_level(contour_level);
         }
         graphics_draw();
      } else {
         std::cout << "WARNING:: molecule " << imol_map << " is not a valid map molecule"
                   << std::endl;
      }
   } else {
      std::cout << "WARNING:: molecule " << imol_map << " is not a valid model molecule"
                << std::endl;
   }
   return 0;
}

void smiles_to_simple_3d(const std::string &smiles) {

   // 20230605-PE this was done via scripting previously. Hmm

   graphics_info_t g;
   short int lang = coot::STATE_PYTHON;
   std::string tlc = "LIG";
   std::vector<coot::command_arg_t> args = { coot::command_arg_t(tlc), smiles };
   std::string sc = g.state_command("generator_3d_import", "new_molecule_by_smiles_string", args, lang);
   PyObject *r1 = safe_python_command_with_return("import generator_3d_import");
   std::cout << "smiles_to_simple_3d(): r1: " << r1 << std::endl;
   std::cout << "smiles_to_simple_3d(): calling this: " << sc << std::endl;
   PyObject *r2 = safe_python_command_with_return(sc);

}

#include "curl-utils.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "utils/subprocess.hpp"

/* ------------------------------------------------------------------------- */
/*                      Acedrg for dictionary                                */
/* ------------------------------------------------------------------------- */

int get_pdbe_cif_for_comp_id(const std::string &comp_id) {
   int status =-1;
   xdg_t xdg;
   std::filesystem::path user_data_dir = xdg.get_data_home();
   std::filesystem::path cif_file_path = user_data_dir / std::string(comp_id + ".cif");
   if (std::filesystem::exists(cif_file_path)) {
      std::cout << "INFO:: found " << cif_file_path << std::endl;
   } else {
      std::string url = "https://www.ebi.ac.uk/pdbe/static/files/pdbechem_v2/" + comp_id + ".cif";
      std::cout << "INFO:: downloading " << url << " "  << cif_file_path << std::endl;
      status = coot_get_url(url, cif_file_path.string());
   }
   return status;
}

bool graphics_info_t::acedrg_running = false; // set the static on startup

void run_acedrg_for_ccd_dict_async(const std::string &residue_name,
                                    const std::string &cif_file_path) {

   auto run_acedrg_func = [] (const std::string &residue_name,
                              const std::string &cif_file_path) {
      std::vector<std::string> cmd_list =
      {"acedrg", "-r", residue_name, "--noGeoOpt", "--coords", "-c", cif_file_path};
      try {
         graphics_info_t::acedrg_running = true;
         subprocess::OutBuffer obuf = subprocess::check_output(cmd_list);
         graphics_info_t::acedrg_running = false;
      }
      catch (const subprocess::CalledProcessError &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
         graphics_info_t::acedrg_running = false; // reset
      }
   };

   auto check_it = +[] () {
      std::cout << "check !!!!! " << graphics_info_t::acedrg_running << std::endl;

      if (graphics_info_t::acedrg_running) {
         return (gboolean)G_SOURCE_CONTINUE;
      } else{
         GtkWidget *w = widget_from_builder("acedrg_running_frame");
         if (w)
            gtk_widget_set_visible(w, FALSE);
         // where is the Acedrg output?
         std::string acedrg_output_file_name = "AcedrgOut.cif";
         if (std::filesystem::exists(acedrg_output_file_name)) {
            read_cif_dictionary(acedrg_output_file_name);
         }
         return (gboolean)G_SOURCE_REMOVE;
      }
   };

   if (!graphics_info_t::acedrg_running) {
      std::thread thread(run_acedrg_func, residue_name, cif_file_path);
      thread.detach();
      GSourceFunc f = GSourceFunc(check_it);
      g_timeout_add(1000, f, nullptr);
      GtkWidget *w = widget_from_builder("acedrg_running_frame");
      if (w)
         gtk_widget_set_visible(w, TRUE);
      else
         std::cout << "ERROR:: can't find acedrg_running_frame" << std::endl;
   } else {
      std::cout << "INFO:: acedrg is already running" << std::endl;
   }

}

void make_acedrg_dictionary_via_CCD_dictionary(int imol, const coot::residue_spec_t &spec) {

   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(spec);
      if (residue_p) {
         std::string residue_name = residue_p->GetResName();
         int status = get_pdbe_cif_for_comp_id(residue_name);
         if (status == CURLE_OK) {
            xdg_t xdg;
            std::filesystem::path user_data_dir = xdg.get_data_home();
            std::filesystem::path cif_file_path = user_data_dir / std::string(residue_name + ".cif");
            run_acedrg_for_ccd_dict_async(residue_name, cif_file_path.string());
         }
      }
   }

}


