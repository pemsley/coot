
#include "Python.h"
#include "graphics-info.h"

// static 
gint
graphics_info_t::drag_refine_refine_intermediate_atoms() {

   int retprog = -1;

#ifdef HAVE_GSL

   thread_for_refinement_loop_threaded();

   //    g.run_post_intermediate_atoms_moved_hook_maybe();  // put this somewhere else
   //   (after the refinement has finished.)                          

#endif // HAVE_GSL

   return retprog;
}

#include "ligand/backrub-rotamer.hh"

// return true if flip moving_atoms_asc was found
bool graphics_info_t::pepflip_intermediate_atoms_other_peptide() {

   bool status = false;
   if (moving_atoms_asc->mol) {

      mmdb::Atom *at_close = NULL;
      float min_dist_sqrd = 4.0;

      coot::Cartesian pt(graphics_info_t::RotationCentre());

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	 coot::Cartesian atom_pos(at->x, at->y, at->z);
	 coot::Cartesian diff = atom_pos - pt;
	 if (diff.amplitude_squared() < min_dist_sqrd) {
	    min_dist_sqrd = diff.amplitude_squared();
	    at_close = at;
	 }
      }

      if (at_close) {
	 mmdb::Residue *residue_p = at_close->residue;
	 if (residue_p) {
	    mmdb::Atom *at_other = 0;
	    std::string at_name_close(at_close->GetAtomName());
            const char *alt_conf = at_close->altLoc;
	    if (at_name_close == " N  ") { // PDBv3 FIXME
	       at_other = residue_p->GetAtom(" CA ", 0, alt_conf);
	    } else {
	       at_other = residue_p->GetAtom(" N ", 0, alt_conf);
	    }
	    status = pepflip_intermediate_atoms(at_other);
	 }
      } else {
         add_status_bar_text("No close atom");
      }
   }
   return status;
}


// return true if flip moving_atoms_asc was found
bool graphics_info_t::pepflip_intermediate_atoms() {

   bool status = false;
   if (moving_atoms_asc->mol) {

      mmdb::Atom *at_close = NULL;
      float min_dist_sqrd = 4.0;

      coot::Cartesian pt(graphics_info_t::RotationCentre());

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	 coot::Cartesian atom_pos(at->x, at->y, at->z);
	 coot::Cartesian diff = atom_pos - pt;
	 if (diff.amplitude_squared() < min_dist_sqrd) {
	    min_dist_sqrd = diff.amplitude_squared();
	    at_close = at;
	 }
      }

      if (at_close) {
	 status = pepflip_intermediate_atoms(at_close);
      }
   }
   return status;
}

// return true if flip moving_atoms_asc was found
bool graphics_info_t::pepflip_intermediate_atoms(mmdb::Atom *at_close) {

   // If we are on CA, C, O or any other atom except N
   // flip (this_res)-(next_res)
   // If we are on N
   // flip (prev_res)-(this_res)

   std::cout << "in pepflip_intermediate_atoms() with at_close " << coot::atom_spec_t(at_close)
             << std::endl;

   bool status = false;

   if (! at_close) {

      std::cout << "INFO:: No close atom" << std::endl;

   } else {

      mmdb::Residue *res_this = at_close->residue;
      std::string atom_name = at_close->name;
      const char *alt_conf = at_close->altLoc;

      // if N is at active atom then we want prev,this
      // otherwise we want this,next
      //
      mmdb::Residue *res_1 = NULL;
      mmdb::Residue *res_2 = NULL;
      if (atom_name == " N  ") {
         res_1 = moving_atoms_asc->get_previous(res_this);
         res_2 = res_this;
      } else {
         res_1 = res_this;
         res_2 = moving_atoms_asc->get_next(res_this);
      }

      if (res_1 && res_2) {
         mmdb::Atom *at_1_ca = res_1->GetAtom(" CA ", 0, alt_conf);
         mmdb::Atom *at_1_c  = res_1->GetAtom(" C  ", 0, alt_conf);
         mmdb::Atom *at_1_o  = res_1->GetAtom(" O  ", 0, alt_conf);
         mmdb::Atom *at_2_ca = res_2->GetAtom(" CA ", 0, alt_conf);
         mmdb::Atom *at_2_n  = res_2->GetAtom(" N  ", 0, alt_conf);
         mmdb::Atom *at_2_h  = res_2->GetAtom(" H  ", 0, alt_conf);

         if (at_1_ca && at_2_ca) {

            // tell the refinement to stop, wait for it to stop, move the atoms and then restart

            continue_threaded_refinement_loop = false;
            while (restraints_lock) {
               std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }

            clipper::Coord_orth base(at_1_ca->x, at_1_ca->y, at_1_ca->z);
            clipper::Coord_orth  top(at_2_ca->x, at_2_ca->y, at_2_ca->z);
            clipper::Coord_orth dir = top - base;
            coot::util::rotate_atom_about(dir, base, M_PI, at_1_c);
            coot::util::rotate_atom_about(dir, base, M_PI, at_1_o);
            coot::util::rotate_atom_about(dir, base, M_PI, at_2_n);
            coot::util::rotate_atom_about(dir, base, M_PI, at_2_h); // does null check

            // in the case where the refinement finishes before the
            // timeout function begins, we would like to force a redraw of the
            // bonds (so we can see that the atoms flipped). Hmm.. I am not
            // sure if this works/is-needed - difficult to test. Maybe
            // it should be added to refinement_of_last_restraints_needs_reset().

            // std::cout << "Atoms really moved - restarting refinement" << std::endl;
            threaded_refinement_loop_counter++;

            refinement_of_last_restraints_needs_reset();
            thread_for_refinement_loop_threaded();

            status = true;
         }
      }
   }
   graphics_draw();
   return status;
}

bool
graphics_info_t::backrub_rotamer_intermediate_atoms() {

   bool state = false;
   if (moving_atoms_asc->mol) {

      mmdb::Atom *at_close = NULL;
      float min_dist_sqrd = 4.0;

      coot::Cartesian pt(graphics_info_t::RotationCentre());

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	 coot::Cartesian atom_pos(at->x, at->y, at->z);
	 coot::Cartesian diff = atom_pos - pt;
	 if (diff.amplitude_squared() < min_dist_sqrd) {
	    min_dist_sqrd = diff.amplitude_squared();
	    at_close = at;
	 }
      }

      if (at_close) {

	 std::string chain_id = at_close->GetChainID();
	 int res_no = at_close->GetSeqNum();
	 std::string ins_code = at_close->GetInsCode();
	 std::string alt_conf = at_close->altLoc;
	 mmdb::Manager *mol = moving_atoms_asc->mol;
	 mmdb::Residue *this_res = at_close->residue;
	 mmdb::Residue *next_res = coot::util::get_following_residue(coot::residue_spec_t(this_res), mol);
	 mmdb::Residue *prev_res = coot::util::get_previous_residue(coot::residue_spec_t(this_res), mol);
	 int imol_map = Imol_Refinement_Map();
	 if (is_valid_map_molecule(imol_map)) {
	    if (this_res && prev_res && next_res) {
	       std::string monomer_type = this_res->GetResName();
	       std::pair<short int, coot::dictionary_residue_restraints_t> p =
		  Geom_p()->get_monomer_restraints(monomer_type, coot::protein_geometry::IMOL_ENC_ANY);
	       const coot::dictionary_residue_restraints_t &rest = p.second;

	       if (p.first) {
		  try {

		     // we can set the idx of the atoms in this_res, prev_res, next_res
		     // here if we need to speed up the update_moving_atoms_from_molecule_atoms()
		     // function.

		     coot::backrub br(chain_id, this_res, prev_res, next_res, alt_conf, mol,
				      &molecules[imol_map].xmap); // use a pointer for the map
		     std::pair<coot::minimol::molecule,float> m = br.search(rest);

		     continue_threaded_refinement_loop = false;
		     while (restraints_lock) {
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		     }

		     update_moving_atoms_from_molecule_atoms(m.first);
		     state = true;

		     refinement_of_last_restraints_needs_reset();
		     thread_for_refinement_loop_threaded();

		     // drag_refine_refine_intermediate_atoms();
		     // graphics_draw();
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "WARNING:: thrown " << rte.what() << std::endl;
		  }
	       } else {
		  std::string m = "Can't get all the residues needed for rotamer fit";
		  add_status_bar_text(m);
	       }
	    }
	 }
      }
   }
   return state;
}

// return true if the isomerisation was made
//
bool
graphics_info_t::cis_trans_conversion_intermediate_atoms() {

   // prefer to choose the peptide that contains a cis-peptide if you can find one on
   // either side of the given residue
   //
   // PDBv3 FIXME
   //
   auto get_this_and_next_residue_for_peptide_link = [] (atom_selection_container_t *asc,
                                                         mmdb::Atom *closest_atom) {

                                                        mmdb::Residue *residue_of_closest_atom = closest_atom->residue;
                                                        mmdb::Residue *res_this = residue_of_closest_atom;
                                                        mmdb::Residue *res_prev = asc->get_previous(residue_of_closest_atom);
                                                        mmdb::Residue *res_next = asc->get_next(residue_of_closest_atom);
                                                        std::pair<bool, double> omega_1 = coot::util::omega_torsion(res_prev, res_this, "");
                                                        std::pair<bool, double> omega_2 = coot::util::omega_torsion(res_this, res_next, "");

                                                        if (omega_1.first) {
                                                           if (coot::util::is_cis(omega_1.second)) {
                                                              return std::make_pair(res_prev, res_this);
                                                           } else {
                                                              if (omega_2.first) {
                                                                 if (coot::util::is_cis(omega_2.second)) {
                                                                    return std::make_pair(res_this, res_next);
                                                                 } else {
                                                                    // they are both not cis, so choose pair depending on
                                                                    // the closest atom
                                                                    std::string name_of_closest_atom(closest_atom->GetResName());
                                                                    if (name_of_closest_atom == " N  ") {
                                                                       return std::make_pair(res_this, res_next);
                                                                    } else {
                                                                       return std::make_pair(res_prev, res_this);
                                                                    }
                                                                 }
                                                              } else {
                                                                 // no next residue
                                                                 return std::make_pair(res_prev, res_this);
                                                              }
                                                           }
                                                        } else {
                                                           // no prev residue
                                                           return std::make_pair(res_this, res_next);
                                                        }
                                                        return std::pair<mmdb::Residue *, mmdb::Residue *>(res_this, res_next);
                                                     };

   bool state = false;
   if (moving_atoms_asc->mol) {

      mmdb::Atom *at_close = NULL;
      float min_dist_sqrd = 4.0;

      coot::Cartesian pt(graphics_info_t::RotationCentre());

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	 coot::Cartesian atom_pos(at->x, at->y, at->z);
	 coot::Cartesian diff = atom_pos - pt;
	 if (diff.amplitude_squared() < min_dist_sqrd) {
	    min_dist_sqrd = diff.amplitude_squared();
	    at_close = at;
	 }
      }

      if (at_close) {

         std::pair<mmdb::Residue *, mmdb::Residue *> res_pair = get_this_and_next_residue_for_peptide_link(moving_atoms_asc, at_close);
	 mmdb::Residue *res_1 = res_pair.first;
	 mmdb::Residue *res_2 = res_pair.second;

	 if (res_1 && res_2) {
	    mmdb::Atom *at_1_ca = res_1->GetAtom(" CA ");
	    mmdb::Atom *at_1_c  = res_1->GetAtom(" C  ");
	    mmdb::Atom *at_2_ca = res_2->GetAtom(" CA ");

	    if (at_1_ca && at_2_ca && at_1_c) {

	       mmdb::Manager *standard_residues_mol = standard_residues_asc.mol;
	       mmdb::Manager *mol = moving_atoms_asc->mol;

	       // tell the refinement to stop, wait for it to stop, move the atoms and then restart

	       continue_threaded_refinement_loop = false;
	       while (restraints_lock) {
		  std::this_thread::sleep_for(std::chrono::milliseconds(10));
	       }

               std::pair<bool, double> omega = coot::util::omega_torsion(res_1, res_2, "");
               if (omega.first) {
                  bool is_currently_cis = coot::util::is_cis(omega.second);

                  //  now we can manipulate the coordinates
                  coot::util::cis_trans_conversion(res_1, res_2, mol, standard_residues_mol);

                  if (is_currently_cis) { // well, before the conversion, I mean
                     last_restraints->add_trans_peptide_restraint(res_1, res_2);
                  } else {
                     // it was trans, now is cis and so we want to remove the trans-peptide restraint
                     last_restraints->remove_trans_peptide_restraint(res_1, res_2);
                  }
               }

               // deliberately outside previous test
               refinement_of_last_restraints_needs_reset(); // maybe not needed.
               thread_for_refinement_loop_threaded();
            }
	 }
      }
   }
   return state;
}


void
graphics_info_t::update_moving_atoms_from_molecule_atoms(const coot::minimol::molecule &mm) {

   if (moving_atoms_asc) {
      if (moving_atoms_asc->n_selected_atoms) {

	 int imod = 1;
	 mmdb::Model *model_p = moving_atoms_asc->mol->GetModel(imod);
	 if (! model_p) {
	    std::cout << "Null model in update_moving_atoms_from_molecule_atoms() " << std::endl;
	 } else {

	    for (unsigned int ifrag_mm=0; ifrag_mm<mm.fragments.size(); ifrag_mm++) {
	       const coot::minimol::fragment &frag = mm.fragments[ifrag_mm];
	       const std::string &frag_chain_id = frag.fragment_id;

	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  std::string moving_atom_chain_id = chain_p->GetChainID();
		  if (frag_chain_id == moving_atom_chain_id) {
		     for (int ires_mm=frag.min_res_no(); ires_mm<=frag.max_residue_number(); ires_mm++) {
			const coot::minimol::residue &residue_mm = frag[ires_mm];
			int nres = chain_p->GetNumberOfResidues();
			for (int ires=0; ires<nres; ires++) {
			   mmdb::Residue *residue_p = chain_p->GetResidue(ires);
			   if (residue_p->GetSeqNum() == residue_mm.seqnum) {
			      std::string ins_code(residue_p->GetInsCode());
			      if (ins_code == residue_mm.ins_code) {
				 for (unsigned int iatom_mm=0; iatom_mm<residue_mm.atoms.size(); iatom_mm++) {
				    const coot::minimol::atom &atom_mm = residue_mm.atoms[iatom_mm];
				    const std::string &atom_name_mm = atom_mm.name;
				    int n_atoms = residue_p->GetNumberOfAtoms();
				    for (int iat=0; iat<n_atoms; iat++) {
				       mmdb::Atom *at = residue_p->GetAtom(iat);
				       std::string atom_name(at->GetAtomName());
				       if (atom_name == atom_name_mm) {
					  std::string altLoc_mm = atom_mm.altLoc;
					  std::string at_altloc = at->altLoc;
					  if (at_altloc == altLoc_mm) {
					     at->x = atom_mm.pos.x();
					     at->y = atom_mm.pos.y();
					     at->z = atom_mm.pos.z();
					     break; // we found the right atom
					  }
				       }
				    }
				 }
				 break; // found the right residue
			      }
			   }
			}
		     }
		     break; // found the right chain
		  }
	       }
	    }
	 }
      }
   }
}



#ifdef USE_PYTHON
void graphics_info_t::register_post_intermediate_atoms_moved_hook(PyObject *function) {

   std::cout << "::::::::::: set post_intermediate_atoms_moved_hook to " << function << std::endl;
   post_intermediate_atoms_moved_hook = function;

}
#endif



#ifdef USE_PYTHON
void graphics_info_t::run_post_intermediate_atoms_moved_hook_maybe() {

   if (post_intermediate_atoms_moved_hook) {

      graphics_info_t g;
      PyObject *o = g.get_intermediate_atoms_bonds_representation();

      if (PyBool_Check(o)) {
	 // no useful representation, (must be False)
      } else {
	 PyObject *py_main = PyImport_AddModule("__main__");
	 if (py_main) {
	    // wraps PyDict_SetItemString (String because the key is a string (the variable name))
	    PyModule_AddObject(py_main, "intermediate_atoms_representation_internal", o);
	 }
      }
   }
}
#endif

// static
mmdb::Atom *
graphics_info_t::get_moving_atoms_active_atom(const coot::Cartesian &rc, float within_radius_limit) {
   mmdb::Atom *at_moving = 0;
   if (moving_atoms_asc) {
      if (moving_atoms_asc->n_selected_atoms > 0) {
         float min_dist_sqrd = 4.0;
         float best_dist_sqrd = min_dist_sqrd;
         for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
            mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
            coot::Cartesian atom_pos(at->x, at->y, at->z);
            coot::Cartesian diff = atom_pos - rc;
            if (diff.amplitude_squared() < best_dist_sqrd) {
               if (diff.amplitude_squared() < min_dist_sqrd) {
                  min_dist_sqrd = diff.amplitude_squared();
                  at_moving = at;
                  best_dist_sqrd = diff.amplitude_squared();
               }
            }
         }
      }
   }
   return at_moving;
};


#include "ligand/side-chain.hh"

//static
int
graphics_info_t::side_chain_flip_180_intermediate_atoms() {

   int status = 0;

   std::cout << "start side_chain_flip_180_intermediate_atoms()" << std::endl;

   if (moving_atoms_asc) {
      if (moving_atoms_asc->n_selected_atoms > 0) {
	 status = 1;
         coot::Cartesian rc = get_rotation_centre();
         float within_radius_limit = 9.0;
         mmdb::Atom *at = get_moving_atoms_active_atom(rc, within_radius_limit);
         if (at) {

            // do I need to get the refinement lock here?

            std::cout << "found at " << at << std::endl;
            coot::atom_spec_t atom_spec(at);
            coot::residue_spec_t spec(atom_spec);
            std::string alt_conf(at->altLoc);
            mmdb::Manager *mol = moving_atoms_asc->mol;
            coot::do_180_degree_side_chain_flip(spec, alt_conf, mol, Geom_p());
	    refinement_of_last_restraints_needs_reset();
	    thread_for_refinement_loop_threaded();
            graphics_draw();
         } else {
            std::cout << "Failed to find active atom in moving atoms " << std::endl;
         }
      } else {
         std::cout << "INFO:: no moving atoms - this should not happen" << std::endl;
      }
   } else {
      std::cout << "INFO:: no moving atoms " << std::endl;
   }
   return status;
}


// Note that this works different to all other button operattions!
// This is not a key-binding/active atom style action either
//
// This operation allows you to pick on intermediate atoms (to activate the flip)
//
// I am not sure that I like it or want it.
//
void
graphics_info_t::side_chain_flip_180_moving_atoms_residue(const coot::residue_spec_t &spec,
                                                          const std::string &alt_conf) {

   if (moving_atoms_displayed_p()) {

      coot::do_180_degree_side_chain_flip(spec, alt_conf, moving_atoms_asc->mol, Geom_p());
      graphics_draw();
      pick_pending_flag = 0;
      normal_cursor();
      model_fit_refine_unactive_togglebutton("model_refine_dialog_do_180_degree_sidechain_flip_togglebutton");
   }
}


#include "coot-utils/jed-flip.hh"

// static
int
graphics_info_t::jed_flip_intermediate_atoms() {

   int status = 0;

   if (moving_atoms_asc) {
      if (moving_atoms_asc->n_selected_atoms > 0) {
	 status = 1;

	 // get active atom
	 mmdb::Atom *active_atom = nullptr;
	 float min_dist_sqrd = 4.0;

	 coot::Cartesian pt(RotationCentre_x(),
			    RotationCentre_y(),
			    RotationCentre_z());

	 for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	    mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	    coot::Cartesian atom_pos(at->x, at->y, at->z);
	    coot::Cartesian diff = atom_pos - pt;
	    if (diff.amplitude_squared() < min_dist_sqrd) {
	       min_dist_sqrd = diff.amplitude_squared();
	       active_atom = at;
	    }
	 }

	 if (active_atom) {
	    mmdb::Residue *residue_p = active_atom->residue;
	    int imol = imol_moving_atoms;
	    bool invert_selection = false;
	    coot::util::jed_flip(imol, residue_p, active_atom, invert_selection, Geom_p());

	    // add_drag_refine_idle_function();
	    // drag_refine_refine_intermediate_atoms();

	    refinement_of_last_restraints_needs_reset();
	    thread_for_refinement_loop_threaded();
	 }
      }
   }
   graphics_draw();
   return status;

}

#include "coot-utils/atom-selection-container.hh"
#include "ideal/crankshaft.hh"

// static
int
graphics_info_t::crankshaft_peptide_rotation_optimization_intermediate_atoms() {

   // there is some repeated code here, consider factoring it out.

   int n_threads = coot::get_max_number_of_threads() - 1;
   if (n_threads < 1) n_threads = 1;

   int status = 0;
   if (moving_atoms_asc) {
      if (moving_atoms_asc->n_selected_atoms > 0) {
	 status = 1;

	 // get active atom
	 mmdb::Atom *active_atom = nullptr;
	 float min_dist_sqrd = 4.0;

	 coot::Cartesian pt(RotationCentre_x(),
			    RotationCentre_y(),
			    RotationCentre_z());

	 for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	    mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	    coot::Cartesian atom_pos(at->x, at->y, at->z);
	    coot::Cartesian diff = atom_pos - pt;
	    if (diff.amplitude_squared() < min_dist_sqrd) {
	       min_dist_sqrd = diff.amplitude_squared();
	       active_atom = at;
	    }
	 }

	 if (active_atom) {
	    mmdb::Residue *residue_p = active_atom->residue;
	    int imol = imol_moving_atoms;
	    coot::residue_spec_t residue_spec(residue_p);
	    unsigned int n_peptides = 3;
	    int n_samples = -1; // auto

	    graphics_info_t g;
	    int imol_map = g.Imol_Refinement_Map();
	    if (is_valid_map_molecule(imol_map)) {
	       const clipper::Xmap<float> &xmap = molecules[imol_map].xmap;
	       float w = g.geometry_vs_map_weight;
	       int n_solutions = 1;
	       std::vector<mmdb::Manager *> mols =
		  coot::crankshaft::crank_refine_and_score(residue_spec, n_peptides, xmap,
							   moving_atoms_asc->mol, w,
							   n_samples, n_solutions, &static_thread_pool, n_threads);
	       if (mols.size() == 1) { // not 0

		  // this feels super-dangerous, replacing the coordinates of the atoms
		  // of moving-atoms asc with those of the crankshaft mol.
		  // There are more crankshaft atoms than moving atoms (moving atoms
		  // don't include the fixed atoms used in refinement).
		  //
		  atom_selection_container_t asc = make_asc(mols[0]);
		  for (int iat=0; iat<moving_atoms_asc->n_selected_atoms; iat++) {
		     if (iat<asc.n_selected_atoms) {
			mmdb::Atom *at = moving_atoms_asc->atom_selection[iat];
			mmdb::Atom *asc_at = asc.atom_selection[iat];
			at->x = asc_at->x;
			at->y = asc_at->y;
			at->z = asc_at->z;
		     }
		  }
		  // add_drag_refine_idle_function(); // pre-threaded refinement
		  // drag_refine_refine_intermediate_atoms();
		  refinement_of_last_restraints_needs_reset();
		  thread_for_refinement_loop_threaded();
	       } else {
		  g.add_status_bar_text("Couldn't crankshaft this");
	       }
	    }
	 }
      }
   }
   graphics_draw();
   return status;
}



// static
void
graphics_info_t::rebond_molecule_corresponding_to_moving_atoms() {

   if (moving_atoms_asc) {

      // we were doing some refinement and decided to give up/reject, in that case, we neeed to redraw
      // the the "static" molecule from which the moving atoms were derived because, during refinement,
      // the "static" atoms were drawn without the atoms that correspond to the moving atoms.
      //
      if (is_valid_model_molecule(imol_moving_atoms)) {
	 std::set<int> empty_set;
	 graphics_info_t::molecules[graphics_info_t::imol_moving_atoms].make_bonds_type_checked(empty_set);
      }
   }
}

void
graphics_info_t::set_regenerate_bonds_needs_make_bonds_type_checked(bool state) {
   regenerate_bonds_needs_make_bonds_type_checked_flag = state;
}

bool
graphics_info_t::get_regenerate_bonds_needs_make_bonds_type_checked_state() {
   return regenerate_bonds_needs_make_bonds_type_checked_flag;
}

