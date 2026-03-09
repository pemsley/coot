/* ligand/ligand.cc
 *
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Author Paul Emsley
 * Copyright 2008, 2009 by The University of Oxford
 * Copyright 2012, 2015, 2016 by Medical Research Council
 * Author Paul Emsley
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
 * 02110-1301, USA.
 */

// This file was written as the USA started its (2003) attack on Iraq.
//

// #include <stdio.h> // for snprintf

#include <fstream>

#include <queue> // new fangled idea from Kevin (to do the flood
                 // filling and keeping the ligand centres via
                 // Coord_grids).

#include <algorithm>

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if defined _MSC_VER
// #define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
// #define S_IWUSR S_IWRITE
#define snprintf _snprintf
#else
#include <unistd.h>
#endif

// #if !defined(WINDOWS_MINGW) && !defined(_MSC_VER)
// #include <pwd.h>
// #include <sys/types.h>
// #endif

#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/core/xmap.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_interp.h"

#include "ligand.hh" // has mmdb-manager because mmdb::PPAtom is part
                     // of mask_map interface.

#include <mmdb2/mmdb_coormngr.h> // for GetMassCenter()

#include "utils/win-compat.hh"
#include "utils/coot-utils.hh"

#include "coords/mmdb-extras.hh"   // 220403
#include "coords/mmdb.hh"

#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "utils/logging.hh"

extern logging logger;


std::pair<coot::minimol::molecule, coot::minimol::molecule>
coot::make_mols_from_atom_selection_string(mmdb::Manager *mol,
                                           std::string atom_selection_string,
                                           bool fill_masking_molecule_flag) {
   int SelHnd = mol->NewSelection();
   mol->Select(SelHnd, mmdb::STYPE_ATOM, atom_selection_string.c_str(), mmdb::SKEY_NEW);
   mmdb::PPAtom atom_selection = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(SelHnd, atom_selection, n_selected_atoms);

   mmdb::Manager *mol_from_selected =
      coot::util::create_mmdbmanager_from_atom_selection(mol, SelHnd);
   // atom selection in mol gets inverted by this function:
   mmdb::Manager *mol_from_non_selected =
      coot::util::create_mmdbmanager_from_atom_selection(mol, SelHnd, 1);

   coot::minimol::molecule range_mol  = coot::minimol::molecule(mol_from_selected);
   coot::minimol::molecule masked_mol = coot::minimol::molecule(mol_from_non_selected);
   delete mol_from_selected;
   delete mol_from_non_selected;
   mol->DeleteSelection(SelHnd);
   return std::pair<coot::minimol::molecule, coot::minimol::molecule> (masked_mol, range_mol);
}

std::pair<coot::minimol::molecule, coot::minimol::molecule>
coot::make_mols_from_atom_selection(mmdb::Manager *mol,
                                    int udd_atom_selection_fitting_atoms,
                                    bool fill_masking_molecule_flag) {

   // the caller is control of the atom selection, we don't delete it here.

   mmdb::PPAtom atom_selection = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(udd_atom_selection_fitting_atoms, atom_selection, n_selected_atoms);

   mmdb::Manager *mol_from_selected =
      coot::util::create_mmdbmanager_from_atom_selection(mol, udd_atom_selection_fitting_atoms, false);
   // atom selection in mol gets inverted by this function:
   mmdb::Manager *mol_from_non_selected =
      coot::util::create_mmdbmanager_from_atom_selection(mol, udd_atom_selection_fitting_atoms, true);

   coot::minimol::molecule range_mol  = coot::minimol::molecule(mol_from_selected);
   coot::minimol::molecule masked_mol = coot::minimol::molecule(mol_from_non_selected);
   delete mol_from_selected;
   delete mol_from_non_selected;
   return std::pair<coot::minimol::molecule, coot::minimol::molecule> (masked_mol, range_mol);
}


coot::ligand::ligand() {

   n_clusters = 0;
   write_solutions = 1; // by default write out pdb files.
   write_orientation_solutions = 0; // by default, don't write orientations.
   write_raw_waters = 0; // default not, Eleanor wants them on
                         // however, so a function to turn them on is
                         // provided.

   water_molecule_volume = 11.0; // 11 captures the phosphate and the
                                 // ligand in the tutorial, 15 does
                                 // not!

   dont_test_rotations = 0;  // do test (the following) rotations by default.
   fit_fraction = 0.75;
   masked_map_val = 0.0;
   map_atom_mask_radius = 2.0;
   verbose_reporting = 0;
   do_cluster_size_check_flag = 1;
   do_chemically_sensible_test_flag = 1;
   do_sphericity_test_flag = 1;
   default_b_factor = 30.0;

   var_limit = 0.12;
   water_to_protein_distance_lim_min = 2.4;
   water_to_protein_distance_lim_max = 3.2;

   // We mask the map several times in water-flood mode.  We want to
   // calculate and keep the stats from the first masked map.  the
   // pair.first shows that we have sensible value in pair.second.
   //
   xmap_masked_stats.first = 0;

   // Size match test is on for ligands and off for rigid body
   // refinement (which calls the find_centre_by_ligand to make it so)
   //
   do_size_match_test = 1;

   // After we have found the initial site, we need to test the 4
   // different (180 degree) rotations of the ligand that matches the
   // eigen vectors:
   origin_rotations.push_back(clipper::Mat33<double>(1,0,0,0, 1,0,0,0, 1));
   origin_rotations.push_back(clipper::Mat33<double>(1,0,0,0,-1,0,0,0,-1)); // x rotation
   origin_rotations.push_back(clipper::Mat33<double>(-1,0,0,0,1,0,0,0,-1)); // y rotation
   origin_rotations.push_back(clipper::Mat33<double>(-1,0,0,0,-1,0,0,0,1)); // z rotation

   // To get the angle contributions, we need to find the vector
   // perpendicular to the vector that joins the rotation centre to
   // the atom position.  And we need that, in the plane corresponding
   // to the axis about which we are rotating, (i.e. we need the
   // vector in the XY plane for rotations about z (i.e. the z
   // component should be zero).
   //
   rotation_component[0] = clipper::RTop_orth(clipper::Mat33<double>(0,0,0,0,0,1,0,-1,0),
                                              clipper::Coord_frac(0,0,0));
   rotation_component[1] = clipper::RTop_orth(clipper::Mat33<double>(0,0,-1,0,0,0,1,0,0),
                                              clipper::Coord_frac(0,0,0));
   rotation_component[2] = clipper::RTop_orth(clipper::Mat33<double>(0,1,0,-1,0,0,0,0,0),
                                              clipper::Coord_frac(0,0,0));

   // gets set later, hopefully
   map_rms = -1.0;
}

// coot::ligand::~ligand() {

//    std::cout << "destructing ligands..." << std::endl;
//    if (fitted_ligand_vec.size() > 0) {
//       // need to delete the ligand because of mmdb usage
//       for (int i=0; i<fitted_ligand_vec.size(); i++) {
//          std::cout << "destructing ligand " << i << std::endl;
//          delete fitted_ligand_vec[i];
//       }
//    }
// }


// Set xmap and xmap_pristine from the input map:
//
void
coot::ligand::import_map_from(const clipper::Xmap<float> &map_in) {

   // xmap_pristine.init(map_in.spacegroup(), map_in.cell(), map_in.grid_sampling());
   xmap_pristine = map_in;
   xmap_cluster = xmap_pristine;
   xmap_masked = xmap_pristine;
   clipper::Xmap_base::Map_reference_index ix;
   float sum_mean = 0;
   float sum_sq = 0;
   int n = 0;
   for (ix = map_in.first(); !ix.last(); ix.next() )  { // iterator index.
      float v = map_in[ix];
      sum_mean += v;
      sum_sq += v * v;
      n++;
   }
   if (n > 0) {
      float mean = sum_mean/float(n);
      map_rms = sqrt(sum_sq/float(n) - mean*mean);
   }
   calculate_gradient_scale();
}



// Set xmap and xmap_pristine from the input map:
//
// Create a new form of the function that sets th e rms of the map,
// when then gets used in a kludgey way to set the gradient scale.
//
void
coot::ligand::import_map_from(const clipper::Xmap<float> &map_in, float rms_map_in) {

   xmap_pristine = map_in;
   xmap_cluster = xmap_pristine;
   xmap_masked  = xmap_pristine;
   map_rms = rms_map_in;
   calculate_gradient_scale();
}


short int
coot::ligand::map_fill_from_mtz(std::string mtz_file_name,
                                std::string f_col,
                                std::string phi_col,
                                std::string weight_col,
                                short int use_weights,
                                short int is_diff_map,
                                float map_sampling_rate) { // 1.5 default

   std::cout << "............................. map_fill_from_mtz " << mtz_file_name << std::endl;
   clipper::HKL_info myhkl;
   clipper::MTZdataset myset;
   clipper::MTZcrystal myxtl;

  std::cout << "reading mtz file " << mtz_file_name << std::endl;
  if (! coot::is_regular_file(mtz_file_name))
     return 0;

  clipper::CCP4MTZfile mtzin;
  mtzin.open_read( mtz_file_name );       // open new file
  mtzin.import_hkl_info( myhkl );         // read sg, cell, reso, hkls
  clipper::HKL_data< clipper::datatypes::F_sigF<float> >   f_sigf_data(myhkl, myxtl);
  clipper::HKL_data< clipper::datatypes::Phi_fom<float> > phi_fom_data(myhkl, myxtl);
  clipper::HKL_data< clipper::datatypes::F_phi<float> >       fphidata(myhkl, myxtl);


  if ( use_weights ) {
     clipper::String dataname = "/*/*/[" + f_col + " " + f_col + "]";
     std::cout << dataname << "\n";
     mtzin.import_hkl_data(  f_sigf_data, myset, myxtl, dataname );
     dataname = "/*/*/[" + phi_col + " " + weight_col + "]";
     std::cout << dataname << "\n";
     mtzin.import_hkl_data( phi_fom_data, myset, myxtl, dataname );
     mtzin.close_read();
     std::cout << "We should use the weights: " << weight_col << std::endl;

     fphidata.compute(f_sigf_data, phi_fom_data,
                      clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

  } else {
     clipper::String dataname = "/*/*/[" + f_col + " " + phi_col + "]";
     mtzin.import_hkl_data(     fphidata, myset, myxtl, dataname );
     mtzin.close_read();
  }
  std::cout << "Number of reflections: " << myhkl.num_reflections() << "\n";

  std::cout << "finding ASU unique map points..." << std::endl;
  clipper::Grid_sampling gs(myhkl.spacegroup(),
                            myhkl.cell(),
                            myhkl.resolution(),
                            map_sampling_rate);
  xmap_pristine.init(myhkl.spacegroup(), myhkl.cell(), gs);

  std::cout << "Grid..." << std::string(xmap_pristine.grid_sampling().format()).c_str() << "\n";

  std::cout << "doing fft..." << std::endl;
  xmap_pristine.fft_from( fphidata );                  // generate map
  std::cout << "done fft..." << std::endl;

  map_statistics(); //set map_rms
  xmap_cluster = xmap_pristine; // xmap_cluster gets scribbled on later.
  xmap_masked = xmap_pristine;
  calculate_gradient_scale();
  return 1;
}

// Manipulate the map.
// Where there are grid points close to atoms, set them to 0.0 density.
//
// The interface used by c-interface.cc's execute_ligand_search():
//
void
coot::ligand::mask_map(mmdb::Manager *mol, short int mask_waters_flag) {

   float atom_radius = map_atom_mask_radius; // Angstroems

   mmdb::PPAtom atoms;

   int n_atoms = make_selected_atoms(&atoms, mol); // modify atoms
   protein_atoms.init(mol);

   // xmap_pre_cluster = xmap; // copy to the map that (will be) masked not clustered.

   mmdb::realtype xmc, ymc, zmc; // filled by reference
   GetMassCenter (atoms, n_atoms, xmc, ymc, zmc);
   protein_centre = clipper::Coord_orth(xmc, ymc, zmc);
   logger.log(log_t::INFO, logging::ltw("Protein centre at: "), logging::ltw(protein_centre.format()));
   // std::cout << "INFO:: Protein centre at: " << protein_centre.format() << std::endl;

   // std::cout << "masking....";
   for(int i=0; i<n_atoms; i++) {
      clipper::Coord_orth co(atoms[i]->x, atoms[i]->y, atoms[i]->z);

      std::string res_name (atoms[i]->residue->name);
      if (mask_waters_flag) {
         mask_around_coord(co, atom_radius);  // mask xmap_masked
      } else {
         // don't mask waters
         if (res_name != "WAT" && res_name != "HOH" ) {
//             std::cout << "DEBUG:: masking around " << co.format() << " "
//                       << atom_radius << std::endl;
            mask_around_coord(co,atom_radius);  // masks xmap_cluster (implicitly).
         }
      }
   }
   xmap_masked = xmap_cluster;
   // std::cout << "masking done." << std::endl;
}

void
coot::ligand::mask_map(mmdb::Manager *mol,
                       int SelectionHandle,
                       short int invert_flag) {
   // if invert_flag is 0 put the map to 0
   // where the atoms are.  If 1, put map to
   // 0 where the atoms are not (e.g. ligand
   // density figure)

   mmdb::PPAtom atom_selection = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);
   logger.log(log_t::INFO, logging::ltw("Masking around "), logging::ltw(n_selected_atoms), logging::ltw(" atoms"));
   // std::cout << "INFO:: Masking around " << n_selected_atoms << " atoms" << std::endl;

   if (invert_flag == 0) {
      for (int i=0; i<n_selected_atoms; i++) {
         clipper::Coord_orth co(atom_selection[i]->x,
                                atom_selection[i]->y,
                                atom_selection[i]->z);
         mask_around_coord(co, map_atom_mask_radius); // mask xmap_cluster
      }

   } else {

      // e.g. for ligand figure: Let's make a real mask map that
      // contains only 1 and 0, first fill it with 1s and then mask
      // out the atoms (set those densities to 0).
      //
      // So now we have a map that contains 0 where we what density,
      // so we run through both maps looking at were the real masking
      // map is 1 and set those points to 0.0 in the xmap_masked.
      //
      // xmap_masked is the map that mask_around_coord without a map
      // pointer argument masks, of course.

      clipper::Xmap<int> mask(xmap_cluster.spacegroup(),
                              xmap_cluster.cell(),
                              xmap_cluster.grid_sampling());
      clipper::Xmap_base::Map_reference_index ix;
      for (ix=mask.first(); !ix.last(); ix.next()) {
         mask[ix] = 1;
      }
      // now mask the mask map
      for (int i=0; i<n_selected_atoms; i++) {
         mmdb::Atom *at = atom_selection[i];
         clipper::Coord_orth co = coot::co(at);
         mask_around_coord(co, map_atom_mask_radius, &mask);
      }

      // compare map and mask and reset where atom are not:
      for (ix=mask.first(); !ix.last(); ix.next()) {
         if (mask[ix] == 1) {
            xmap_cluster[ix] = 0.0;
         }
      }
   }
   xmap_masked = xmap_cluster;
}

// coot usage (coot has beforehand prepared a mol that does not have in it the
// residue of interest).
//
void
coot::ligand::mask_map(const minimol::molecule &mol, short int mask_waters_flag) {

   protein_atoms = mol;
   mask_map(mask_waters_flag);
}

void
coot::ligand::mask_map(bool mask_waters_flag) {

   xmap_cluster = xmap_pristine; // copy to the map that (will be) masked not clustered.
   float atom_radius = map_atom_mask_radius; // Angstroems (this is a guess).

   std::cout << "masking....";
   for (unsigned int ifrag=0; ifrag<protein_atoms.fragments.size(); ifrag++) {
      for (int ires=protein_atoms.fragments[ifrag].min_res_no();
           ires<=protein_atoms.fragments[ifrag].max_residue_number();
           ires++) {
         if (mask_waters_flag) {

            // no special test needed
            for (unsigned int iatom=0;
                 iatom<protein_atoms.fragments[ifrag][ires].atoms.size();
                 iatom++) {
               mask_around_coord(protein_atoms[ifrag][ires][iatom].pos, atom_radius);
            }

         } else {
            // we need to test if this residue is a water
            if (protein_atoms[ifrag][ires].name != "WAT" &&
                protein_atoms[ifrag][ires].name != "HOH") {
               for (unsigned int iatom=0;
                    iatom<protein_atoms.fragments[ifrag][ires].atoms.size();
                    iatom++) {
                  mask_around_coord(protein_atoms[ifrag][ires][iatom].pos, atom_radius);
               }
            }
         }
      }
   }
   xmap_masked = xmap_cluster;
   std::cout << "masking done\n";
   // output_map(xmap_masked, "post-protein-masking.map"); // debugging flood mode
}

// Masks xmap_cluster now.
//
// (Really?)
//
void
coot::ligand::mask_around_coord(const clipper::Coord_orth &co, float atom_radius) {
   clipper::Coord_frac cf = co.coord_frac(xmap_cluster.cell());

   clipper::Coord_frac box0(
                            cf.u() - atom_radius/xmap_cluster.cell().descr().a(),
                            cf.v() - atom_radius/xmap_cluster.cell().descr().b(),
                            cf.w() - atom_radius/xmap_cluster.cell().descr().c());

   clipper::Coord_frac box1(
                            cf.u() + atom_radius/xmap_cluster.cell().descr().a(),
                            cf.v() + atom_radius/xmap_cluster.cell().descr().b(),
                            cf.w() + atom_radius/xmap_cluster.cell().descr().c());

   clipper::Grid_map grid(box0.coord_grid(xmap_cluster.grid_sampling()),
                          box1.coord_grid(xmap_cluster.grid_sampling()));

   float atom_radius_sq = atom_radius * atom_radius;
   int nhit = 0;
   int nmiss = 0;

   clipper::Xmap_base::Map_reference_coord ix( xmap_cluster, grid.min() ), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
         for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
            if ( (iw.coord().coord_frac(xmap_cluster.grid_sampling()).coord_orth(xmap_cluster.cell()) - co).lengthsq() < atom_radius_sq) {
                // std::cout << "masked " << masked_map_val << " point at "
               // << iw.coord().coord_frac(xmap_masked.grid_sampling()).coord_orth(xmap_masked.cell()).format()
               // << " centre point: " << co.format() << " "
               //                           << (iw.coord().coord_frac(xmap_masked.grid_sampling()).coord_orth(xmap_masked.cell()) - co).lengthsq()
               //                           << std::endl;
               xmap_cluster[iw] = masked_map_val;
               nhit++;
            } else {
               nmiss++;
            }
         }
      }
   }
   //    std::cout << "nhit " << nhit << " nmiss " << nmiss << std::endl;
}

void
coot::ligand::mask_around_coord(const clipper::Coord_orth &co, float atom_radius,
                                clipper::Xmap<float> *xmap_p) {  // alter xmap_p

   clipper::Coord_frac cf = co.coord_frac(xmap_p->cell());
   clipper::Coord_frac box0(cf.u() - atom_radius/xmap_p->cell().descr().a(),
                            cf.v() - atom_radius/xmap_p->cell().descr().b(),
                            cf.w() - atom_radius/xmap_p->cell().descr().c());

   clipper::Coord_frac box1(cf.u() + atom_radius/xmap_p->cell().descr().a(),
                            cf.v() + atom_radius/xmap_p->cell().descr().b(),
                            cf.w() + atom_radius/xmap_p->cell().descr().c());

   clipper::Grid_map grid(box0.coord_grid(xmap_p->grid_sampling()),
                          box1.coord_grid(xmap_p->grid_sampling()));

   float atom_radius_sq = atom_radius * atom_radius;
   float masked_map_val = 0;

   clipper::Xmap_base::Map_reference_coord ix(*xmap_p, grid.min()), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
         for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
            if ( (iw.coord().coord_frac(xmap_p->grid_sampling()).coord_orth(xmap_p->cell()) - co).lengthsq() < atom_radius_sq) {
               (*xmap_p)[iw] = masked_map_val;
            }
         }
      }
   }
}

// Templete this.
void
coot::ligand::mask_around_coord(const clipper::Coord_orth &co, float atom_radius,
                                clipper::Xmap<int> *xmap_p) {  // alter xmap_p

   clipper::Coord_frac cf = co.coord_frac(xmap_p->cell());
   clipper::Coord_frac box0(cf.u() - atom_radius/xmap_p->cell().descr().a(),
                            cf.v() - atom_radius/xmap_p->cell().descr().b(),
                            cf.w() - atom_radius/xmap_p->cell().descr().c());

   clipper::Coord_frac box1(cf.u() + atom_radius/xmap_p->cell().descr().a(),
                            cf.v() + atom_radius/xmap_p->cell().descr().b(),
                            cf.w() + atom_radius/xmap_p->cell().descr().c());

   clipper::Grid_map grid(box0.coord_grid(xmap_p->grid_sampling()),
                          box1.coord_grid(xmap_p->grid_sampling()));

   float atom_radius_sq = atom_radius * atom_radius;
   int masked_map_val = 0;

   clipper::Xmap_base::Map_reference_coord ix(*xmap_p, grid.min()), iu, iv, iw;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
         for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
            if ( (iw.coord().coord_frac(xmap_p->grid_sampling()).coord_orth(xmap_p->cell()) - co).lengthsq() < atom_radius_sq) {
               (*xmap_p)[iw] = masked_map_val;
            }
         }
      }
   }
}

void
coot::ligand::set_masked_map_value(float v) {
   masked_map_val = v;
}

int
coot::ligand::make_selected_atoms(mmdb::PPAtom *atoms_p, mmdb::Manager *mol) {

   int selHnd = mol->NewSelection();
   mol->SelectAtoms (selHnd, 0, "*",
                     mmdb::ANY_RES, // starting resno, an int
                     "*", // any insertion code
                     mmdb::ANY_RES, // ending resno
                     "*", // ending insertion code
                     "*", // any residue name
                     "*", // atom name
                     "*", // elements
                     "*"  // alt loc.
                     );
   int nSelAtoms;
   mol->GetSelIndex(selHnd, *atoms_p, nSelAtoms);
   return nSelAtoms;
}


// Using the clipper::MMDB method.
int
coot::ligand::mask_by_atoms(std::string pdb_filename) {

   std::cout << "INFO:: Reading pdb file: " << pdb_filename << std::endl;

   // get mol from pdb_filename:
   atom_selection_container_t asc = get_atom_selection(pdb_filename, true, false, false);

   protein_atoms.init(asc.mol);
   bool mask_waters = 0;
   mask_map(mask_waters); // dont mask waters

   return 1;
}

void
coot::ligand::calculate_gradient_scale() {

   double t = (xmap_pristine.cell().a()/xmap_pristine.grid_sampling().nu() +
               xmap_pristine.cell().c()/xmap_pristine.grid_sampling().nv() +
               xmap_pristine.cell().b()/xmap_pristine.grid_sampling().nw())/3.0;

   gradient_scale = 0.3*t*t;  // somewhat arbitary.

   // std::cout << " DEBUG:: map_rms " << map_rms << std::endl;

   if (map_rms > 0.0)
      gradient_scale *= 0.25/map_rms;    // ditto.

   // std::cout << "The gradient scale will be: " << gradient_scale << std::endl;

}

clipper::Map_stats
coot::ligand::map_statistics() {

   clipper::Map_stats stats(xmap_pristine);

   std::cout << "Map stats:          mean: " << stats.mean() << " and std dev: "
             << stats.std_dev() << std::endl;

   map_rms = stats.std_dev();

   clipper::Map_stats stats_pristine(xmap_pristine);
   std::cout << "Pristine Map stats: mean: " << stats_pristine.mean() << " and std dev: "
             << stats_pristine.std_dev() << std::endl;

   std::cout << "Grid sampling: " << xmap_pristine.grid_sampling().format() << std::endl;
   std::cout << "Cell:          " << xmap_pristine.cell().format() << std::endl;

   return stats;

}



// We are passed the number of standard deviations from the mean
// (above which points are considered as being part of a cluster).
//
// What do we want to get? A list of clusters.
//
// Each cluster has a centre, with scores corresponding to the number
// of points and the density thereof and the list of points.
//
// So every time we meet a new point > cut_off, and when trace_along
// has finished (back here), we add that new cluster to the cluster
// list.
//
// Notice that we trash xmap.
void
coot::ligand::find_clusters_old(float z_cut_off_in) {

   clipper::Map_stats stats(xmap_pristine);
   clipper::Skeleton_basic::Neighbours neighb(xmap_pristine);

   // When the density has been masked very negative (eg. -2.0) the
   // mean is low (e.g.  -0.536255), so we can't use the mean
   // sensibly.  Let's ignore it and just use z sigma.
   //
   // cut_off = stats.mean() + z_cut_off_in*stats.std_dev();
   cut_off = z_cut_off_in*stats.std_dev();
   std::cout << "Using density cut-off: " << cut_off << " sigma ";
   std::cout << " (mean " << stats.mean() << " stdev: " << stats.std_dev()
             << ")" << std::endl;
   // xmap_cluster = xmap_pristine; // masked not clustered.
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = xmap_pristine.first(); !ix.last(); ix.next()) {
      if (xmap_pristine[ix] > cut_off) {
         n_clusters++;
         cluster.push_back(map_point_cluster());
         trace_along(ix.coord(), neighb, n_clusters);
      }
   }
   calculate_cluster_centres_and_eigens();
   std::vector <clipper::Coord_orth> sampled_protein_coords =
      make_sample_protein_coords();
   move_ligand_centres_close_to_protein(sampled_protein_coords);
   std::cout << "There were " << n_clusters << " clusters " << std::endl;
   std::sort(cluster.begin(), cluster.end(), compare_clusters);

   print_cluster_details();
}


// We are passed the number of standard deviations from the mean
// (above which points are considered as being part of a cluster).
//
// What do we want to get? A list of clusters.
//
// Each cluster has a centre, with scores corresponding to the number
// of points and the density thereof and the list of points.
//
// So every time we meet a new point > cut_off, and when trace_along
// has finished (back here), we add that new cluster to the cluster
// list.
//
// Notice that we trash xmap.
void
coot::ligand::find_clusters_int(float z_cut_off_in) {

   clipper::Map_stats stats(xmap_cluster);
   clipper::Skeleton_basic::Neighbours neighb(xmap_cluster);
   cut_off = z_cut_off_in*stats.std_dev();
   std::cout << "Using density cut-off: " << cut_off;
   std::cout << " (mean " << stats.mean() << " stdev: " << stats.std_dev()
             << ")" << std::endl;
   // xmap_cluster = xmap; // masked not clustered.
   clipper::Xmap<int> cluster_map;
   cluster_map.init(xmap_pristine.spacegroup(),
                    xmap_pristine.cell(),
                    xmap_pristine.grid_sampling());

   clipper::Xmap_base::Map_reference_index ix;
   for (ix = cluster_map.first(); !ix.last(); ix.next()) {
      cluster_map[ix] = 0;
   }

   int an_int = 1;
   for (ix = xmap_pristine.first(); !ix.last(); ix.next()) {
      if (xmap_pristine[ix] > cut_off) {
         cluster_map[ix] = an_int;
         an_int++;
      }
   }


   // OK, so now, every point that is part of a cluster has an individual int.
   //
   // Now we want to run over the map lots of times assigning each
   // point to its largest neighbour (if its neighbours are larger
   // than itself):

   int n_changed = 1; // non-zero start
   int n_neighbs = neighb.size();
   int v;

   clipper::Coord_grid c_g;
   clipper::Coord_grid c_g_start;
   while(n_changed) {
      n_changed = 0;
      for (ix = cluster_map.first(); !ix.last(); ix.next()) {
         if (cluster_map[ix]) {
            c_g_start = ix.coord();
            for (int i=0; i<n_neighbs; i++) {
               c_g = c_g_start + neighb[i];
               v = cluster_map.get_data(c_g);
               if (v > cluster_map[ix]) {
                  cluster_map[ix] = v;
                  n_changed++;
               }
            }
         }
      }
      std::cout << "nchanged this round was " << n_changed << std::endl;
   }

   // So now every cluster is marked with a characteristic int.  we
   // need to assign the cluster score and the Coord_grids (map_grid)
   //
   std::vector <int> cluster_index;
   for (ix = cluster_map.first(); !ix.last(); ix.next()) {
      if (cluster_map[ix]) {

      }
   }



   calculate_cluster_centres_and_eigens();
   std::vector <clipper::Coord_orth> sampled_protein_coords =
      make_sample_protein_coords();
   move_ligand_centres_close_to_protein(sampled_protein_coords);
   std::cout << "There were " << n_clusters << " clusters " << std::endl;
   std::sort(cluster.begin(), cluster.end(), compare_clusters);

   print_cluster_details();

}

void
coot::ligand::set_default_b_factor(float f) {
   default_b_factor = f;
}


// We are passed the number of standard deviations from the mean
// (above which points are considered as being part of a cluster).
//
// z_cut_off_in used to be the number of standard deviations of the map
// *after* clustering.
//
// Then chat with Lynn Ten Eyck who said that this was confusing and
// it ought to be the sigma level of the input map (pre-masked
// (i.e. the pristine map)).  So the actual level that we cut at
// should be: z_cut_off_in*map_rms
//
// What do we want to get? A list of clusters.
//
// Each cluster has a centre, with scores corresponding to the number
// of points and the density thereof and the list of points.
//
// So every time we meet a new point > cut_off, and when trace_along
// has finished (back here), we add that new cluster to the cluster
// list.
//
// Notice that we trash cluster_map
//
void
coot::ligand::find_clusters(float z_cut_off_in) {

   std::vector <clipper::Coord_orth> sampled_protein_coords =
      make_sample_protein_coords();

   find_clusters_internal(z_cut_off_in, sampled_protein_coords);

}


void
coot::ligand::find_clusters_internal(float z_cut_off_in,
                                     const std::vector <clipper::Coord_orth> &sampled_protein_coords) {

   std::cout << "INFO:: find_clusters map_rms is " << map_rms << std::endl;

   if (xmap_masked_stats.first == 0) {
      clipper::Map_stats stats(xmap_cluster);
      xmap_masked_stats.first = 1;
      xmap_masked_stats.second.first  = stats.mean();
      xmap_masked_stats.second.second = stats.std_dev();
   }

   clipper::Skeleton_basic::Neighbours neighb(xmap_cluster);

   cut_off = z_cut_off_in * map_rms; // map_rms is from the pristine map

   z_cut_off_in_save = z_cut_off_in; // used in water_fit when we call
                                     // this function again.
   std::cout << "INFO:: Using density cut-off: " << cut_off
             << " (" << z_cut_off_in << " sigma) ";
   std::cout << " (mean " << xmap_masked_stats.second.first
             << " stdev: " << xmap_masked_stats.second.second << ")" << std::endl;
   std::cout << "INFO:: Blobs with volume larger than " << water_molecule_volume
             << " A^3 are too big to be considered waters." << std::endl;
   std::cout << "INFO:: Using water to protein distance limits: "
             << water_to_protein_distance_lim_min << " "
             << water_to_protein_distance_lim_max << std::endl;


   clipper::Xmap_base::Map_reference_index ix;
   std::cout << "INFO:: Finding clusters...";
   std::cout.flush();
   clipper::Xmap<int> cluster_map;
   cluster_map.init(xmap_pristine.spacegroup(), xmap_pristine.cell(),
                    xmap_pristine.grid_sampling());
   for (ix = cluster_map.first(); !ix.last(); ix.next())
      cluster_map[ix] = 0;

   std::queue<clipper::Coord_grid> q;
   clipper::Coord_grid c_g;
   clipper::Coord_grid c_g_start;
   for (ix = xmap_cluster.first(); !ix.last(); ix.next()) {
      if (xmap_cluster[ix] > cut_off) {
         if (! cluster_map[ix]) {
            coot::map_point_cluster mpc;
            c_g_start = ix.coord();
//             mpc.map_grid.push_back(c_g_start);
//             mpc.score += xmap.get_data(c_g_start);
//             cluster_map.set_data(c_g_start, 1);
            q.push(c_g_start);

            //            std::cout << "DEBUG:: q size: " << q.size() << std::endl;
            while (q.size()) {
               c_g_start = q.front();
               q.pop();
               for (int i=0; i<neighb.size(); i++) {
                  c_g = c_g_start + neighb[i];
                  if (xmap_cluster.get_data(c_g) > cut_off) {
                     if (! cluster_map.get_data(c_g)) {
//                         std::cout << "DEBUG:: cluster " << n_clusters
//                                   << " pushing back " << c_g.format()
//                                   << std::endl;
                        cluster_map.set_data(c_g, 1);
                        mpc.map_grid.push_back(c_g);
                        mpc.score += xmap_cluster.get_data(c_g);
                        q.push(c_g);
                     }
                  }
               }
            }

//             std::cout << "pushing back cluster " << n_clusters
//                       << mpc.score << " " << mpc.map_grid.size() << std::endl;
            if (mpc.map_grid.size() > 0) {
               cluster.push_back(mpc);
               n_clusters++;
            }
         }
      }
   }
   std::cout << "done" << std::endl; // finding clusters.

   calculate_cluster_centres_and_eigens();
   move_ligand_centres_close_to_protein(sampled_protein_coords);
   //    std::cout << "There were " << n_clusters << " clusters " << std::endl;
   std::sort(cluster.begin(), cluster.end(), compare_clusters);

   if (verbose_reporting)
      print_cluster_details(true);
}

void
coot::ligand::cluster_from_point(clipper::Coord_orth pt,
                                 float z_cut_off_in) {

   // convert pt to a grid point first.
   clipper::Coord_frac a_cf   = pt.coord_frac(xmap_cluster.cell());
   clipper::Coord_map  a_cm   = a_cf.coord_map(xmap_cluster.grid_sampling());
   clipper::Coord_grid a_grid = a_cf.coord_grid(xmap_cluster.grid_sampling());
   clipper::Skeleton_basic::Neighbours neighb(xmap_cluster);
   float br = 3.0;
   float density_crit = z_cut_off_in * map_rms;

   // std::cout << " DEBUG:: cluster_from_point density_crit " << density_crit << std::endl;

   // Coord_map is like Coord_grid, but non-integer.

   // we need to find a starting point that is in the cluster (close
   // to pt).

   clipper::Coord_grid c_g_start;  // initially unset

   clipper::Coord_grid c_g = a_grid;

   //std::cout << " DEBUG:: cluster_from_point starting at a_grid: " << a_grid.format()
   // << std::endl;

   // scored grid values, in a set of neighbours, get them all and
   // take the top of the sorted list.
   std::vector<std::pair<clipper::Coord_grid, float> > sg;

   float d = xmap_cluster.get_data(a_grid);
   if (d > density_crit) {
      c_g_start = c_g;
   } else {

      bool found_site = 0;
      std::queue<clipper::Coord_grid> q;

      // if no points in the map are above density_crit, then this
      // could loop endlessly, so add a termination condition on
      // counts.
      //
      int count_max = 200000;
      int count = 0;

      // we don't want to add a neighbour if we've added this grid point
      // to the queue already.
      clipper::Xmap<short int> neighbour_map; // cant use bools
      neighbour_map.init(xmap_pristine.spacegroup(), xmap_pristine.cell(),
                         xmap_pristine.grid_sampling());
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = neighbour_map.first(); !ix.last(); ix.next())
         neighbour_map[ix] = 0;

      q.push(a_grid);

      while (! found_site && (count < count_max)) {
         count++;
         c_g = q.front();
         q.pop();
         neighbour_map.set_data(c_g, 1); // mark as examined.
         d = xmap_cluster.get_data(c_g);
         if (d > density_crit) {
            c_g_start = c_g;
            found_site = 1;
         } else {
            for (int i=0; i<neighb.size(); i++) {
               if (neighbour_map.get_data(c_g + neighb[i]) == 0)
                  q.push(c_g + neighb[i]);
            }
         }
      }
      if (! found_site) {
         std::cout << "Hopeless - no ligand density cluster found"
                   << std::endl;
         return;
      }
   }

   // now we have c_g_start

   clipper::Xmap<int> cluster_map;
   cluster_map.init(xmap_pristine.spacegroup(), xmap_pristine.cell(),
                    xmap_pristine.grid_sampling());
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = cluster_map.first(); !ix.last(); ix.next())
      cluster_map[ix] = 0;
   std::queue<clipper::Coord_grid> q2;
   q2.push(c_g_start);
   coot::map_point_cluster mpc;
   while (q2.size()) {
      c_g_start = q2.front();
      q2.pop();
      for (int i=0; i<neighb.size(); i++) {
         c_g = c_g_start + neighb[i];
         if (xmap_cluster.get_data(c_g) > density_crit) {
            if (! cluster_map.get_data(c_g)) {
               cluster_map.set_data(c_g, 1);
               mpc.map_grid.push_back(c_g);
               mpc.score += xmap_cluster.get_data(c_g);
               q2.push(c_g);
            }
         }
      }
   }
   if (mpc.map_grid.size() > 0) {
      cluster.push_back(mpc);
      n_clusters++; // class variable
   }

   calculate_cluster_centres_and_eigens();
   std::vector <clipper::Coord_orth> sampled_protein_coords =
      make_sample_protein_coords();
   move_ligand_centres_close_to_protein(sampled_protein_coords);

   if (verbose_reporting)
      print_cluster_details();
}


int
coot::ligand::n_grid_limit_for_water_cluster() const {
   // very crude.  Fixme.
   //
   float vol = xmap_pristine.cell().volume();
   float ngrid =
      xmap_pristine.grid_sampling().nu() *
      xmap_pristine.grid_sampling().nv() *
      xmap_pristine.grid_sampling().nw();

   float grid_vol = vol/ngrid;
   float n_grid_lim = 15.0/grid_vol; // 4/3 Pi r^3: with r=1.53

   return int(n_grid_lim);
}


// cluster size maxes out, starts new.
void
coot::ligand::find_clusters_water_flood(float z_cut_off_in,
                                        const std::vector <clipper::Coord_orth> &sampled_protein_coords) {
   if (xmap_masked_stats.first == 0) {
      clipper::Map_stats stats(xmap_cluster);
      xmap_masked_stats.first  = 1;
      xmap_masked_stats.second.first  = stats.mean();
      xmap_masked_stats.second.second = stats.std_dev();
      z_cut_off_in_save = z_cut_off_in; // used in water_fit from where we call
                                        // this function again.
   }

   clipper::Skeleton_basic::Neighbours neighb(xmap_cluster);
   cut_off = z_cut_off_in*xmap_masked_stats.second.second;
   std::cout << "Using density cut-off: " << cut_off;
   std::cout << " (mean " << xmap_masked_stats.second.first << " stdev: "
             << xmap_masked_stats.second.second << ")" << std::endl;
   clipper::Xmap_base::Map_reference_index ix;
   std::cout.flush();
   clipper::Xmap<int> cluster_map;
   cluster_map.init(xmap_pristine.spacegroup(), xmap_pristine.cell(),
                    xmap_pristine.grid_sampling());
   for (ix = cluster_map.first(); !ix.last(); ix.next())
      cluster_map[ix] = 0;

   std::queue<clipper::Coord_grid> q;
   clipper::Coord_grid c_g;
   clipper::Coord_grid c_g_start;
   int n_grids_crit = n_grid_limit_for_water_cluster();

   for (ix = xmap_cluster.first(); !ix.last(); ix.next()) {
      if (xmap_cluster[ix] > cut_off) {
         if (! cluster_map[ix]) {
            coot::map_point_cluster mpc;
            c_g_start = ix.coord();
            q.push(c_g_start);

            while (int(q.size()) && (int(mpc.map_grid.size()) < n_grids_crit)) {
               c_g_start = q.front();
               q.pop();
               for (int i=0; i<neighb.size(); i++) {
                  c_g = c_g_start + neighb[i];
                  if (xmap_cluster.get_data(c_g) > cut_off) {
                     if (! cluster_map.get_data(c_g)) {
                        cluster_map.set_data(c_g, 1);
                        mpc.map_grid.push_back(c_g);
                        mpc.score += xmap_cluster.get_data(c_g);
                        q.push(c_g);
                     }
                  }
               }
            }

            if (mpc.map_grid.size() > 0) {
               cluster.push_back(mpc);
               n_clusters++;
            }
         }
      }
   }

   calculate_cluster_centres_and_eigens();
   move_ligand_centres_close_to_protein(sampled_protein_coords);
   std::cout << "There were " << n_clusters << " clusters " << std::endl;
   std::sort(cluster.begin(), cluster.end(), compare_clusters);

   if (verbose_reporting)
      print_cluster_details();

}



// Fake up a cluster centre
//
// Use the installed ligand to provide a centre.  This is for use in
// "rigid body refinement", where the "ligand" (range of residues
// (typically 1)) is at the place from where we start.
//
void
coot::ligand::find_centre_by_ligand(short int do_size_match_flag) {

   if (initial_ligand.size() != 1) {
      std::cout << "initial_ligand size() is " << initial_ligand.size()
                << " we expected to be of size 1" << std::endl;
   } else {

      do_size_match_test = do_size_match_flag;
      find_centre_by_ligand_internal(0);
   }
}

// Is this function of any use really?  Each new ligand generates a
// new cluster against which the ligands (all of them) are tested.
//
void
coot::ligand::find_centres_by_ligands() {

   for(unsigned int ilig=0; ilig<initial_ligand.size(); ilig++) {
      find_centre_by_ligand_internal(ilig);
   }
}

void
coot::ligand::find_centre_by_ligand_internal(int ilig) {

   map_point_cluster mpc;
   std::vector<minimol::atom *> atoms = initial_ligand[ilig].select_atoms_serial();

   clipper::Coord_orth rp(0,0,0);
   for (unsigned int i=0; i<atoms.size(); i++) {
      rp += atoms[i]->pos;
   }
   clipper::Coord_orth mean_pos(0,0,0);
   if (atoms.size() > 0) {
      double scale = 1.0/double(atoms.size());
      mean_pos = clipper::Coord_orth(rp[0]*scale, rp[1]*scale, rp[2]*scale);
   }

   clipper::RTop_orth rtop(initial_ligand_eigenvectors[0], mean_pos);
   mpc.eigenvectors_and_centre = rtop;
   mpc.eigenvalues = initial_ligand_eigenvalues[ilig];
   cluster.push_back(mpc);

   // Now fake in some map_grids so that the ligand size in grid
   // unit approximately matches the ligand size in A^3 from the
   // ligand coordinates.
   //
   // This is no longer needed now that we have do_size_match_test:
   //
//    int ngrids = 7*atoms.size();
//    for(int i=0; i< ngrids; i++)
//       cluster[cluster.size()-1].map_grid.push_back(clipper::Coord_grid());
}

std::vector <clipper::Coord_orth>
coot::ligand::make_sample_protein_coords() const {

   std::vector <clipper::Coord_orth> sample;
   int atom_count = 0;
   int max_atom_loop_count = 4;

   // std::cout << "protein_mmdb has " << protein_atoms.fragments.size() << " fragments\n";
   for(unsigned int ifrag=0; ifrag<protein_atoms.fragments.size(); ifrag++) {
      // very strange.  If we put the limits in the for loop then the loop
      // doesn't get run (Mac).  putting maxresno and minresno here and the
      // loop runs and we find Marianne's error (20050114).  Code corrected.
//       std::cout << "protein_mmdb has " << protein_atoms[ifrag].residues.size()
//                 << " residues in fragment " << ifrag << " with min res no: "
//                 << protein_atoms[ifrag].min_res_no() << " and max res no "
//                 << protein_atoms[ifrag].max_residue_number() << "\n";
      int maxresno = protein_atoms[ifrag].max_residue_number();
      int minresno = protein_atoms[ifrag].min_res_no();
      //       std::cout << minresno << " " << maxresno << std::endl;
      for(int ires=minresno; ires<=maxresno; ires++) {
//          std::cout << "getting " << ires << " of "
//                    << protein_atoms[ifrag].max_residue_number() << "\n";
         for (unsigned int iat=0; iat<protein_atoms[ifrag][ires].atoms.size(); iat++) {
            atom_count++;
            if (atom_count == max_atom_loop_count) {
               atom_count = 0;
               sample.push_back(protein_atoms[ifrag][ires][iat].pos);
            }
         }
      }
   }
   // std::cout << "returning sample of size: " << sample.size() << std::endl;
   return sample;
}

// Move the centres (in eigenvectors_and_centre) and rotate the
// eigenvectors.
//
void
coot::ligand::move_ligand_centres_close_to_protein(const std::vector<clipper::Coord_orth> &sampled_protein_coords) {

   clipper::Coord_orth mv_centre;
   int n = sampled_protein_coords.size();
   if (n > 0) {
      for (unsigned int i=0; i<cluster.size() ; i++)
         //move_ligand_sites_close_to_protein(i);
         move_ligand_site_close_to_protein_using_shape(i, sampled_protein_coords);
   } // otherwise just leave them where they are in the asymmetric unit.

}

void
coot::ligand::calculate_cluster_centres_and_eigens() {

   clipper::Coord_orth running_centre(0.0,0.0,0.0);
   clipper::Coord_frac cf;
   clipper::Coord_orth co;
   double scale;
   clipper::Coord_orth diff_pt;
   clipper::Coord_orth mean_pos;
   double diff_x, diff_y, diff_z;

   for (unsigned int i=0; i<cluster.size(); i++) {
      running_centre = clipper::Coord_orth(0.0, 0.0, 0.0);
      for (unsigned int j=0; j<cluster[i].map_grid.size(); j++) {
         cf = cluster[i].map_grid[j].coord_frac(xmap_pristine.grid_sampling());
         co = cf.coord_orth(xmap_pristine.cell());
         running_centre += co;
      }
      scale = 1/double(cluster[i].map_grid.size()); // scale by number of
                                                    // points
      mean_pos =
         clipper::Coord_orth(running_centre.x() * scale,
                             running_centre.y() * scale,
                             running_centre.z() * scale);
      // cluster[i].centre = mean_pos;
      diff_x = 0; diff_y = 0; diff_z = 0;

      //
      for (unsigned int j=0; j<cluster[i].map_grid.size(); j++) {
         cf = cluster[i].map_grid[j].coord_frac(xmap_pristine.grid_sampling());
         co = cf.coord_orth(xmap_pristine.cell());
         diff_pt = co - mean_pos;
         diff_x += diff_pt.x() * diff_pt.x();
         diff_y += diff_pt.y() * diff_pt.y();
         diff_z += diff_pt.z() * diff_pt.z();
      }
      cluster[i].std_dev = clipper::Coord_orth(sqrt(diff_x*scale),
                                               sqrt(diff_y*scale),
                                               sqrt(diff_z*scale));
      // we need matrix, because Matrix (not Mat33) has eigen.
      clipper::Matrix<double> mat(3,3);
      for (int ii=0; ii<3; ii++)
         for (int jj=0; jj<3; jj++)
            mat(ii,jj) = 0.0;

      for (unsigned int j=0; j< cluster[i].map_grid.size(); j++) {
         cf = cluster[i].map_grid[j].coord_frac(xmap_pristine.grid_sampling());
         co = cf.coord_orth(xmap_pristine.cell());
         mat(0,0) += (co.x() - mean_pos.x()) * (co.x() - mean_pos.x());
         mat(0,1) += (co.x() - mean_pos.x()) * (co.y() - mean_pos.y());
         mat(0,2) += (co.x() - mean_pos.x()) * (co.z() - mean_pos.z());
         mat(1,0) += (co.y() - mean_pos.y()) * (co.x() - mean_pos.x());
         mat(1,1) += (co.y() - mean_pos.y()) * (co.y() - mean_pos.y());
         mat(1,2) += (co.y() - mean_pos.y()) * (co.z() - mean_pos.z());
         mat(2,0) += (co.z() - mean_pos.z()) * (co.x() - mean_pos.x());
         mat(2,1) += (co.z() - mean_pos.z()) * (co.y() - mean_pos.y());
         mat(2,2) += (co.z() - mean_pos.z()) * (co.z() - mean_pos.z());
      }
      std::vector<double> eigens = mat.eigen(true);
      // some jiggery pokery if the mat now has a negative determinant.
      clipper::Mat33<double> m33 = mat33(mat);
      double determinant = m33.det();
      if (determinant < 0) {
         // we need to swap the 1 and 2th rows [leave 0th row - the top]
         // std::cout << "DEBUG:: negative determinant in eigen matrix, swapping columns!\n";
         for (int q = 0; q < 3; q++ )
            clipper::Util::swap(m33(q, 1), m33(q, 2));
         clipper::Util::swap(eigens[1], eigens[2]);
      }
      cluster[i].eigenvectors_and_centre = clipper::RTop_orth(m33, mean_pos);
      cluster[i].eigenvalues = eigens;
   }

//    for (int i=0; i<cluster.size(); i++) {
//       std::cout << "cluster score, points, eigen trn " << i << " "
//                 << cluster[i].score << " " << cluster[i].map_grid.size()
//                 << " " << cluster[i].eigenvectors_and_centre.trn().format() << std::endl;
//    }

}

bool
coot::compare_clusters(const map_point_cluster &a,
                       const map_point_cluster &b) {

   return (a.score > b.score);
}


double
coot::map_point_cluster::volume(const clipper::Xmap<float> &xmap_ref) const {

   double cell_vol = xmap_ref.cell().volume();
   double n_grid_pts =
      xmap_ref.grid_sampling().nu() *
      xmap_ref.grid_sampling().nv() *
      xmap_ref.grid_sampling().nw();
   double grid_point_vol = cell_vol/n_grid_pts;
   return map_grid.size() * grid_point_vol;
}

void
coot::ligand::print_cluster_details(bool show_grid_points) const {

   int ncount = 0;
   int max_clusters = 10;
   if (cluster.size() < 10)
      max_clusters = 10;
   std::cout << "There are " << cluster.size() << " clusters\n";
   std::cout << "Here are the top " << max_clusters << " clusters:\n";
   for (unsigned int i=0; i<cluster.size(); i++) {
      ncount++;
      if (ncount == max_clusters) break;

      std::cout << "  Number: "  << i << " # grid points: "
                << cluster[i].map_grid.size() << " score: "
                << cluster[i].score << "     \n"
                << cluster[i].eigenvectors_and_centre.format() << "   "
                << cluster[i].std_dev.format() << " eigenvalues: "
                << cluster[i].eigenvalues[0] << " "
                << cluster[i].eigenvalues[1] << " "
                << cluster[i].eigenvalues[2] << " "
                << std::endl;

      if (show_grid_points) {
         clipper::Cell cell = xmap_pristine.cell();
         clipper::Grid_sampling gs = xmap_pristine.grid_sampling();
         for (unsigned int j=0; j<cluster[i].map_grid.size(); j++) {
            std::cout << "   "
                      << cluster[i].map_grid[j].format() << " "
                      << cluster[i].map_grid[j].coord_frac(gs).coord_orth(cell).format()
                      << std::endl;
         }
      }
   }
}

float
coot::ligand::get_cluster_volume(unsigned int iclust) const {

   float r = -1;
   if (iclust < cluster.size()) {
      const auto &c(cluster[iclust]);
      float vol = xmap_pristine.cell().volume();
      float ngrid = xmap_pristine.grid_sampling().nu() * xmap_pristine.grid_sampling().nv() * xmap_pristine.grid_sampling().nw();
      float grid_vol = vol/ngrid;
      float cluster_vol = grid_vol * c.map_grid.size();
      r = cluster_vol;
   }
   return r;
 }

void
coot::ligand::output_centres() {

   std::ofstream cen_out("centres.list");
   if (! cen_out) {
      // error
      std::cout << "Could not open " << "centres.list" << " for some reason\n";
   } else {
      for (unsigned int i=0; i<cluster.size(); i++) {
         cen_out << cluster[i].eigenvectors_and_centre.format()
                 << std::endl;
      }
   }
}

void
coot::ligand::trace_along(const clipper::Coord_grid &cg_start,
                          const clipper::Skeleton_basic::Neighbours &neighb,
                          int n_clusters) {

   cluster[n_clusters-1].score += xmap_cluster.get_data(cg_start);
   cluster[n_clusters-1].map_grid.push_back(cg_start);
   xmap_cluster.set_data(cg_start,0.0);

   clipper::Coord_grid c_g;
   for(int i=0; i< neighb.size(); i++) {
      c_g = cg_start + neighb[i];
      if (xmap_cluster.get_data(c_g) > cut_off) {
         trace_along(c_g, neighb, n_clusters);
      }
   }
}


void
coot::ligand::output_map(std::string filename) const {

   clipper::CCP4MAPfile mapout;
   mapout.open_write(filename);
   mapout.export_xmap(xmap_cluster);
   mapout.close_write();
}


void
coot::ligand::output_map(const clipper::Xmap<float> &xmap, const std::string &filename) const {

   clipper::CCP4MAPfile mapout;
   mapout.open_write(filename);
   mapout.export_xmap(xmap);
   mapout.close_write();
}


// Apply symmetry to the ligand centres so that they lie "next to"
// the protein that was the mask.  Pick the closest to the (precalculated)
// centre of the protein and move that cluster centre.
//
// Also adjust the eigenvectors according to the symmetry rotation that
// provided the smallest distance.
//
void
coot::ligand::move_ligand_sites_close_to_protein(int i) {
//
// Help code from Kevin:
// Spacegroup spgr;
// Cell cell;
// Coord_frac cell_shift;
// int i;
// Rtop_orth orthop = RTop_frac( spgr.symop(i).rot(), spgr.symop(i).trn() +
// cell_shift ).rtop_orth( cell );

   clipper::Coord_orth point(cluster[i].eigenvectors_and_centre.trn());

   clipper::Coord_orth s(999999999.9,9999999999.9,9999999999.9);
   clipper::Coord_orth t_point; // calculated for each symm/shift
   float min_dist = 999999999999.9;
   float t_dist;
   clipper::RTop_orth save_transformation(clipper::Mat33<double>(0,0,0,0,0,0,0,0,0),
                                          clipper::Coord_orth(0,0,0));

   int n = xmap_pristine.spacegroup().num_symops();
   clipper::Coord_frac cell_shift;
   for (int ii=0; ii<n; ii++) {
      for (int x_shift = -1; x_shift<2; x_shift++) {
         for (int y_shift = -1; y_shift<2; y_shift++) {
            for (int z_shift = -1; z_shift<2; z_shift++) {

               cell_shift = clipper::Coord_frac(x_shift, y_shift, z_shift);
               clipper::RTop_orth orthop = clipper::RTop_frac(xmap_pristine.spacegroup().symop(ii).rot(), xmap_pristine.spacegroup().symop(ii).trn() + cell_shift).rtop_orth(xmap_pristine.cell());
               t_point = point.transform(orthop);
               t_dist = clipper::Coord_orth::length(protein_centre,t_point);
               if (t_dist < min_dist) {
                  s = t_point;
                  min_dist = t_dist;
                  save_transformation = orthop;
               }
            }
         }
      }
   }

   cluster[i].eigenvectors_and_centre = clipper::RTop_orth(save_transformation * cluster[i].eigenvectors_and_centre);
   // cluster[i].eigenvectors_and_centre = r;

}


// This is the same function as above, except the testing function
// uses sampled_protein_coords (not just a simple distance check to
// the centre of the protein).
//
void
coot::ligand::move_ligand_site_close_to_protein_using_shape (int iclust,
                                                             const std::vector<clipper::Coord_orth> &sampled_protein_coords) {

   clipper::Coord_orth point(cluster[iclust].eigenvectors_and_centre.trn());


   clipper::Coord_orth s(999999999.9,9999999999.9,9999999999.9);
   clipper::Coord_orth t_point; // calculated for each symm/shift
   float min_dist = 999999999999.9;
   float t_dist;
   clipper::RTop_orth save_transformation(clipper::Mat33<double>(0,0,0,0,0,0,0,0,0),
                                          clipper::Coord_orth(0,0,0));  // was unset


   int n_sampled = sampled_protein_coords.size();

////BEGIN NEW CODE RESOLVES ISSUE WITH DISAPERING LIGANDS AND WATERS
   clipper::Vec3<double>  vcrd(0.0,0.0,0.0),resvec;
   clipper::Mat33<double> imat;
   int nn=0,mm=0,kk=0;
   for(int i=0;i<n_sampled;i++) {
      vcrd+=sampled_protein_coords[i];
   }
   vcrd = 1.0/((float)sampled_protein_coords.size())*vcrd;
   imat          = (xmap_pristine.cell()).matrix_frac();
   resvec = imat*vcrd;
   nn = floor(resvec[0]); mm = floor(resvec[1]); kk = floor(resvec[2]); //GEOM-CENTER INDICES
////END CODE

   if (n_sampled > 0) {
      int n = xmap_pristine.spacegroup().num_symops();
      clipper::Coord_frac cell_shift;
      for (int isym=0; isym<n; isym++) {
         for (int x_shift = -1; x_shift<2; x_shift++) {
            for (int y_shift = -1; y_shift<2; y_shift++) {
               for (int z_shift = -1; z_shift<2; z_shift++) {
                  cell_shift = clipper::Coord_frac(x_shift+nn, y_shift+mm, z_shift+kk);
                  clipper::RTop_orth orthop = clipper::RTop_frac(xmap_pristine.spacegroup().symop(isym).rot(), xmap_pristine.spacegroup().symop(isym).trn() + cell_shift).rtop_orth(xmap_pristine.cell());

                  t_point = point.transform(orthop);
                  t_dist = min_dist_to_protein(t_point, sampled_protein_coords);
                  if (t_dist < min_dist) {
                     s = t_point;
                     min_dist = t_dist;
                     save_transformation = orthop;
                  }
               }
            }
         }
      }
   }
   cluster[iclust].eigenvectors_and_centre = clipper::RTop_orth(save_transformation * cluster[iclust].eigenvectors_and_centre);

}

double
coot::ligand::min_dist_to_protein(const clipper::Coord_orth &point,
                                  const std::vector<clipper::Coord_orth> &sampled_protein_coords) const {

   double dist = 9999999.9;
   double this_dist;
   int n = sampled_protein_coords.size();
   if (n > 0) {
      for (int i=0; i<n; i++) {
         this_dist = clipper::Coord_orth::length(point, sampled_protein_coords[i]);
         if (this_dist < dist)
            dist = this_dist;
      }
   } else {
      dist = 0.0;
   }
   return dist;
}


// Turn the cluster centres into anisotropic atoms, so that we can see
// the direction eigenvalues as if it were anisotropy of an atom.
// This involves some tricky conversion of the eigenvectors to
// anitropic Us (the "diagonalized form"?).
void
coot::ligand::make_pseudo_atoms() {

//    clipper::MMDB mmdb(xmap_pristine.spacegroup(), xmap_pristine.cell());
//    clipper::NDBModel model;
//    clipper::NDBChain chain;
//    clipper::NDBResidue residue;
//    clipper::NDBAtom atom;
//    clipper::DBModel m1 = mmdb.add_model(model);
//    clipper::DBChain c1 = m1.add_chain(chain);

//    model.set_id("A");
//    chain.set_id("A");
//    residue.set_type("GLY");
//    residue.set_seqnum(1);
//    atom = clipper::NDBAtom::null();
//    atom.set_occupancy(1.0);
//    for (unsigned int i=0; i<cluster.size(); i++) {

//       if (cluster[i].map_grid.size() > 2) {
//          atom.set_coord_orth(clipper::Coord_orth(cluster[i].eigenvectors_and_centre.trn()));
//          atom.set_element(" C");
//          atom.set_type(" CA ");

//          // Matrix calculation stuff
//          clipper::Mat33<double> diag; // diagonal eigenvalues
//          diag = clipper::Mat33<double>::identity();
//          for (int j=0; j<3; j++)
//             diag(j,j) = 0.15*sqrt(cluster[i].eigenvalues[j]);
//          clipper::Mat33<double> e = cluster[i].eigenvectors_and_centre.rot();
//          clipper::Mat33<double> r = e * diag * e.transpose();
//          //       std::cout << "r"    << std::endl <<    r.format() << std::endl;
//          //       std::cout << "diag" << std::endl << diag.format() << std::endl;
//          clipper::U_aniso_orth u(r(0,0),r(1,1),r(2,2),r(0,1),r(0,2),r(1,2));
//          atom.set_u_aniso_orth(u);
//          clipper::DBResidue r1 = c1.add_residue(residue);
//          r1.set_seqnum(i+1);
//          r1.add_atom(atom);
//       }
//    }
//    mmdb.finalise_edit();
// //    mmdb.write_file("aniso-ligand-centres.pdb");
//    std::string file_name("aniso-ligand-centres.pdb");
//    mmdb.pcmmdbmanager()->WritePDBASCII( (char *)file_name.c_str() );

}

clipper::Mat33<double>
coot::ligand::mat33(const clipper::Matrix<double> &mat) const {

   clipper::Mat33<double> m;
   for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
         m(i,j) = mat(i,j);
      }
   }
   return m;
}

// Get the centre of initial ligand and the eigen vectors and values
//
void
coot::ligand::make_ligand_properties(int ilig) {

   std::vector<minimol::molecule> a;
   fitted_ligand_vec.push_back(a);

   if (ilig >= int(initial_ligand_model_centre.size())) {
      initial_ligand_model_centre.resize(ilig+1);
      initial_ligand_eigenvectors.resize(ilig+1);
      initial_ligand_eigenvalues.resize (ilig+1);
   }

   // Get Centre:
   //
   // select all atoms
   std::vector<minimol::atom *> atoms = initial_ligand[ilig].select_atoms_serial();

   if (atoms.size() == 0) {
      std::cout << "ERROR in ligand coordinates - none found" << std::endl;
      // exit(1); // No exit() calls from a library
      return;
      //
      // But wrong for interactive program.  To be fixed for that case.
   }
   //
   // Get ligand eigen values and vectors (for the class variables)
   //
   // First calculate the mean position (ignoring atomic mass)
   //
   clipper::Coord_orth running_centre(0.0, 0.0, 0.0);
   for (unsigned int i=0; i<atoms.size(); i++) {
      running_centre += atoms[i]->pos;
   }
   double scale = 1/double(atoms.size());
   clipper::Coord_orth mean_pos(running_centre.x() * scale,
                                running_centre.y() * scale,
                                running_centre.z() * scale);

//    std::cout << "inital position comparison: "
//              << initial_ligand_model_centre.format() << std::endl
//              << mean_pos.format() << std::endl;

   initial_ligand_model_centre[ilig] = mean_pos;

   clipper::Matrix<double> mat(3,3);
   for (int ii=0; ii<3; ii++)
      for (int jj=0; jj<3; jj++)
         mat(ii,jj) = 0.0;
   for (unsigned int i=0; i<atoms.size(); i++) {
      clipper::Coord_orth co = atoms[i]->pos;
               mat(0,0) += (co.x() - mean_pos.x()) * (co.x() - mean_pos.x());
         mat(0,1) += (co.x() - mean_pos.x()) * (co.y() - mean_pos.y());
         mat(0,2) += (co.x() - mean_pos.x()) * (co.z() - mean_pos.z());
         mat(1,0) += (co.y() - mean_pos.y()) * (co.x() - mean_pos.x());
         mat(1,1) += (co.y() - mean_pos.y()) * (co.y() - mean_pos.y());
         mat(1,2) += (co.y() - mean_pos.y()) * (co.z() - mean_pos.z());
         mat(2,0) += (co.z() - mean_pos.z()) * (co.x() - mean_pos.x());
         mat(2,1) += (co.z() - mean_pos.z()) * (co.y() - mean_pos.y());
         mat(2,2) += (co.z() - mean_pos.z()) * (co.z() - mean_pos.z());
   }
   // now using the class variable initial_ligand_eigenvectors

   initial_ligand_eigenvalues[ilig]  = mat.eigen(true); // overwrites itself.
   initial_ligand_eigenvectors[ilig] = mat33(mat);
   clipper::Mat33<double> m2 = mat33(mat);

   if (m2.det() < 0) {
      m2(0,2) = -m2(0,2);
      m2(1,2) = -m2(1,2);
      m2(2,2) = -m2(2,2);
      initial_ligand_eigenvectors[ilig] = m2;
   }

   if (verbose_reporting) {
      std::cout << "ligand eigen values: "
                << initial_ligand_eigenvalues[ilig][0] << "  "
                << initial_ligand_eigenvalues[ilig][1] << "  "
                << initial_ligand_eigenvalues[ilig][2] << "  " << std::endl;
      std::cout << " ligand eigen vectors: " << std::endl
                << initial_ligand_eigenvectors[ilig].format() << std::endl;
   }
//    std::cout << "DEBUG:: initial_ligand_eigenvectors[" << ilig << "] determinant "
//              << initial_ligand_eigenvectors[ilig].det() << std::endl;


}


// // This is the prefered method of getting the ligand info into the class.
// //
// //
// void
// coot::ligand::install_ligand(const clipper::MMDB &ligand_in) {

//    int ilig = initial_ligand.size();
//    initial_ligand.resize(ilig+1);

//    initial_ligand[ilig].init(ligand_in.pcmmdbmanager());
//    make_ligand_properties(ilig);
// }


void
coot::ligand::install_ligand(mmdb::Manager *mol) {

   int ilig = initial_ligand.size();
   initial_ligand.resize(ilig+1);

   initial_ligand[ilig].init(mol);
   make_ligand_properties(ilig);
}

//
void
coot::ligand::install_ligand(std::string ligand_pdb_filename) {

   std::cout << "Reading ligand pdb file: " << ligand_pdb_filename << std::endl;

   int ilig = initial_ligand.size();
   initial_ligand.resize(ilig+1);
   initial_ligand[ilig].read_file(ligand_pdb_filename);
   make_ligand_properties(ilig);
}

void
coot::ligand::install_ligand(const coot::minimol::molecule &ligand) {

   int ilig = initial_ligand.size();
   initial_ligand.resize(ilig+1);

   initial_ligand[ilig] = ligand;
   make_ligand_properties(ilig);
//    std::cout << "DEBUG:: initial_ligand now of size: " << initial_ligand.size()
//              << std::endl;
}


std::ostream&
coot::operator<<(std::ostream &s, const coot::ligand_score_card &lsc) {

   s << "[ligand-score: #" << lsc.ligand_no << " at-score: " <<  lsc.get_score()
     << " r-state: [" << lsc.correlation.first;

   if (lsc.correlation.first)
      s << " correl-score " << lsc.correlation.second;
   s << "]";

   s << " (atom-score: "
     << lsc.score_per_atom << ") many-atoms-fit: " << lsc.many_atoms_fit
     << " n-atoms: " << lsc.n_ligand_atoms << "]";
   return s;
}



// Externally accessible routine.
//
// Go down the cluster list, starting at biggest cluster fitting ligands
// until the cluter number is max_placements.
//
void
coot::ligand::fit_ligands_to_clusters(int max_n_clusters) {

   final_ligand.resize(max_n_clusters);
   save_ligand_score.resize(max_n_clusters);

   for (int iclust=0; iclust<int(cluster.size()) && iclust<max_n_clusters; iclust++) {
     fit_ligands_to_cluster(iclust);
   }
}

#include <atomic>
#include <thread>

#include "utils/split-indices.hh"

// We are given a particular site.  Fit all different ligand types to
// that site and write out the best one if it has a score above 0.0.
//
// Surely n_lig_max is the maximum *cluster* index?
//
void
coot::ligand::fit_ligands_to_cluster(int iclust) {

   // for debugging
   write_orientation_solutions = 0;
   bool debug = false;

   minimol::molecule ior_holder;

   coot::ligand_score_card this_scorecard;

   // Just like the helix algorithm, we search all eigenvectors.
   clipper::Mat33<double> no_rotation    (1, 0,  0, 0, 1,  0, 0, 0, 1);
   clipper::Mat33<double> y_axis_rotation(0, 0, -1, 0, 1,  0, 1, 0, 0);
   clipper::Mat33<double> x_axis_rotation(1, 0,  0, 0, 0, -1, 0, 1, 0);
   clipper::RTop_orth no_rotation_op(no_rotation, clipper::Coord_orth(0,0,0));
   clipper::RTop_orth  y_axis_op(y_axis_rotation, clipper::Coord_orth(0,0,0));
   clipper::RTop_orth  x_axis_op(x_axis_rotation, clipper::Coord_orth(0,0,0));

   std::vector<clipper::RTop_orth> eigen_orientations;
   eigen_orientations.push_back(no_rotation_op);
   eigen_orientations.push_back(y_axis_op);
   eigen_orientations.push_back(x_axis_op);

   std::atomic<bool> results_lock(false);
   auto get_results_lock = [&results_lock] () {
                              bool unlocked = false;
                              while (! results_lock.compare_exchange_weak(unlocked, true)) {
                                 std::this_thread::sleep_for(std::chrono::microseconds(1));
                                 unlocked = false;
                              }
                           };
   auto release_results_lock = [&results_lock] () {
                                results_lock = false;
                             };

   // I wasn't sure how to make these functions static, so we have a local versions.

   auto cluster_ligand_size_match = [] (const map_point_cluster &cluster, const minimol::molecule &ligand, float grid_vol) {
                                       float cluster_vol = grid_vol * cluster.map_grid.size();
                                       unsigned int n_atoms = ligand.count_atoms();
                                       float ligand_vol = (4.0 * 3.14159/3.0) * 1.78 * static_cast<float>(n_atoms);
                                       return ((ligand_vol/cluster_vol < 7.0) && (ligand_vol/cluster_vol > 0.8));
                                    };

   auto transform_ligand_atom = [] (const clipper::Coord_orth &position,
                                    const clipper::RTop_orth &cluster_rtop,
                                    const clipper::Mat33<double> &initial_ligand_eigenvector,
                                    const clipper::Coord_orth &initial_ligand_model_centre,
                                    const clipper::Mat33<double> &origin_rotation) {
                                   clipper::Coord_orth zero(0,0,0);
                                   clipper::RTop_orth ligand_op(initial_ligand_eigenvector, initial_ligand_model_centre);
                                   clipper::RTop_orth lopi(ligand_op.inverse());
                                   clipper::RTop_orth origin_rotation_op(origin_rotation, zero);
                                   clipper::Coord_orth a = position.transform(lopi);
                                   a = a.transform(origin_rotation_op);
                                   a = a.transform(cluster_rtop);
                                   return a;
                                };

   auto fit_ligand_copy = [transform_ligand_atom] (unsigned int ilig,
                                                   const map_point_cluster &cluster,
                                                   const minimol::molecule &ligand_in,
                                                   const clipper::Mat33<double> &initial_ligand_eigenvector,
                                                   const clipper::Coord_orth &initial_ligand_model_centre,
                                                   const clipper::RTop_orth &eigen_orientation,
                                                   const clipper::Xmap<float> &xmap_masked,
                                                   const clipper::Xmap<float> &xmap_pristine,
                                                   const clipper::RTop_orth rotation_component[3],
                                                   float gradient_scale) {

                             auto ligand = ligand_in;
                             std::vector<minimol::atom *> atoms_p = ligand.select_atoms_serial();

                             // First move the ligand to the site of the cluster:
                             for(unsigned int ii=0; ii<atoms_p.size(); ii++)
                                atoms_p[ii]->pos = transform_ligand_atom(atoms_p[ii]->pos, cluster.eigenvectors_and_centre,
                                                                         initial_ligand_eigenvector, initial_ligand_model_centre,
                                                                         eigen_orientation.rot());

                             if (false) {
                                std::cout << "initial_ligand_eigenvector:\n" << initial_ligand_eigenvector.format() << std::endl;
                                std::cout << "initial_ligand_model_centre:\n" << initial_ligand_model_centre.format() << std::endl;
                                std::cout << "eigen_orientation:\n" << eigen_orientation.format() << std::endl;
                                for(unsigned int ii=0; ii<atoms_p.size(); ii++)
                                   std::cout << "lf: ilig " << ilig
                                             << " post-transform_ligand_atom " << ii << atoms_p[ii]->name << " "
                                             << atoms_p[ii]->pos.format() << std::endl;
                             }

                             float fit_fraction = 0.1;
                             ligand_score_card lsc_pre = score_orientation(atoms_p, std::cref(xmap_pristine), fit_fraction);
                             rigid_body_refine_ligand(&atoms_p, std::cref(xmap_masked), std::cref(xmap_pristine),
                                                      rotation_component, gradient_scale); // ("rigid body") move atoms.
                             ligand_score_card lsc_post = score_orientation(atoms_p, std::cref(xmap_pristine), fit_fraction);
                             lsc_post.set_ligand_number(ilig);
                             if (false)
                                std::cout << "lf: ilig " << ilig << " pre-score " << lsc_pre.get_score() <<  " "
                                          <<  " post " << lsc_post.get_score()  << std::endl;
                             return std::make_pair(ligand, lsc_post);
                          };

   auto fit_ligand_to_cluster = [cluster_ligand_size_match,
                                 fit_ligand_copy,
                                 get_results_lock,
                                 release_results_lock] (unsigned int ilig, unsigned int iclust, // needed for debugging
                                                        const map_point_cluster &cluster,
                                                        const minimol::molecule &ligand,
                                                        const clipper::Mat33<double> &initial_ligand_eigenvector,
                                                        const clipper::Coord_orth &initial_ligand_model_centre,
                                                        const std::vector<clipper::Mat33<double> > &origin_rotations,
                                                        const std::vector<clipper::RTop_orth> &eigen_orientations,
                                                        const clipper::Xmap<float> &xmap_masked,
                                                        const clipper::Xmap<float> &xmap_pristine,
                                                        const clipper::RTop_orth rotation_component[3],
                                                        float gradient_scale,
                                                        float grid_vol,
                                                        bool do_size_match_test,
                                                        bool dont_test_rotations,
                                                        bool debug,
                                                        std::vector<std::pair<minimol::molecule, ligand_score_card> > &results) {

                                   if (!do_size_match_test || cluster_ligand_size_match(cluster, ligand, grid_vol)) {
                                      if (debug)
                                         std::cout << "ligand " << ilig << " passes the size match test " << "for cluster number "
                                                   << iclust << std::endl;

                                      int n_rot = origin_rotations.size();
                                      int n_eigen_oris = 3;
                                      if (dont_test_rotations) {
                                         n_rot = 1;  // the first one is the identity matrix
                                         n_eigen_oris = 1;
                                      }

                                      n_rot = 1; // all the rots give the same solution - I don't know why.

                                      for (int i_eigen_ori=0; i_eigen_ori<n_eigen_oris; i_eigen_ori++) {
                                         for (int ior=0; ior<n_rot; ior++) {
                                            std::pair<minimol::molecule, ligand_score_card> scored_ligand =
                                               fit_ligand_copy(ilig, cluster, ligand,
                                                               initial_ligand_eigenvector,
                                                               initial_ligand_model_centre,
                                                               eigen_orientations[i_eigen_ori],
                                                               std::cref(xmap_masked),
                                                               std::cref(xmap_pristine),
                                                               rotation_component, gradient_scale);
                                            get_results_lock();
                                            results.push_back(scored_ligand);
                                            release_results_lock();
                                         }
                                      }
                                   }
                                };

   float vol = xmap_pristine.cell().volume();
   float ngrid = xmap_pristine.grid_sampling().nu() * xmap_pristine.grid_sampling().nv() * xmap_pristine.grid_sampling().nw();
   float grid_vol = vol/ngrid;

   // FIXME
   // remember to add back the code for debugging solutions post-generation after multithreading

   unsigned int n_ligands = initial_ligand.size();
   // std::cout << "INFO:: fitting " << n_ligands << " ligands into cluster " << iclust << std::endl;

   if (false) {
      if (! initial_ligand.empty()) {
         std::vector<coot::minimol::atom *> atoms_p = initial_ligand[0].select_atoms_serial();
         for (unsigned int iat=0; iat<atoms_p.size(); iat++) {
            const auto &atom = *(atoms_p[iat]);
            std::cout << "   inital ligand 0 " << atom.name << " " << atom.pos.format() << std::endl;
         }
      }
   }

   std::vector<std::pair<minimol::molecule, ligand_score_card> > big_vector_of_results;

   // we split them into batches of 4 - that's fast enough for a nice speed up, but won't
   // create a massive (20Gb) memory load. I don't know how to find out where the memory
   // is allocated.
   //
   std::vector<std::vector<unsigned int> > indices;
   unsigned int n_per_batch = 10;
   unsigned int n_batches = n_ligands / n_per_batch + 1;
   coot::split_indices(&indices, n_ligands, n_batches);

   for (unsigned int ii=0; ii<indices.size(); ii++) {
      const std::vector<unsigned int> &ligand_indices = indices[ii];
      std::vector<std::thread> threads;
      for (unsigned int jj=0; jj<ligand_indices.size(); jj++) {
         const unsigned int &ilig = ligand_indices[jj];
         const minimol::molecule &ligand = initial_ligand[ilig];
         threads.push_back(std::thread(fit_ligand_to_cluster,
                                       ilig, iclust,
                                       std::cref(cluster[iclust]),
                                       std::cref(ligand),
                                       std::cref(initial_ligand_eigenvectors[ilig]),
                                       std::cref(initial_ligand_model_centre[ilig]),
                                       std::cref(origin_rotations),
                                       std::cref(eigen_orientations),
                                       std::cref(xmap_masked),
                                       std::cref(xmap_pristine),
                                       rotation_component, // make this a vector?
                                       gradient_scale,
                                       grid_vol,
                                       do_size_match_test,
                                       dont_test_rotations,
                                       debug,
                                       std::ref(big_vector_of_results)));
      }
      for (unsigned int ilig=0; ilig<ligand_indices.size(); ilig++)
         threads[ilig].join();
   }

   if (write_orientation_solutions) {
      // we no longer have access to i_eigen_ori or ior here, I think. Maybe they can be added to
      // ligand_score_card?
      // so that write_orientation_solution(mol) becomes a member function of ligand_score_card?
      //
      // old-function: write_orientation_solution(iclust, ilig, i_eigen_ori, ior, fitted_ligand_vec[ilig][iclust]);
   }

   final_ligand[iclust] = big_vector_of_results;

   sort_final_ligand(iclust);
}

void
coot::ligand::write_orientation_solution(unsigned int iclust,
                                         unsigned int ilig,
                                         unsigned int i_eigen_ori,
                                         unsigned int ior,
                                         const coot::minimol::molecule &mol) const {

   std::string ori_sol_file_name = "ori-sol-cluster:_";
   ori_sol_file_name += util::int_to_string(iclust);
   ori_sol_file_name += "-ligno:_";
   ori_sol_file_name += util::int_to_string(ilig);
   ori_sol_file_name += "-eigen:_";
   ori_sol_file_name += util::int_to_string(i_eigen_ori);
   ori_sol_file_name += "-ori:_";
   ori_sol_file_name += util::int_to_string(ior);
   ori_sol_file_name += ".pdb";
   fitted_ligand_vec[ilig][iclust].write_file(ori_sol_file_name, default_b_factor);

}



// sort the final_ligand for a given cluster, so that the best-scoring
// solution is in position 0.
//
// changes final_ligand vector
//
void
coot::ligand::sort_final_ligand(unsigned int iclust) {

   std::sort(final_ligand[iclust].begin(),
             final_ligand[iclust].end(),
             compare_scored_ligands);
   // lowest score is now in 0th position
   std::reverse(final_ligand[iclust].begin(),
                final_ligand[iclust].end());

   if (false)
      for (unsigned int isol=0; isol<final_ligand[iclust].size(); isol++)
         std::cout << "post reverse: solution " << isol << " of " << final_ligand[iclust].size()
                   << " " << final_ligand[iclust][isol].second << std::endl;
}

// static
bool
coot::ligand::compare_scored_ligands(const std::pair<coot::minimol::molecule, ligand_score_card> &sl_1,
                                     const std::pair<coot::minimol::molecule, ligand_score_card> &sl_2) {
   return (sl_1.second.get_score() < sl_2.second.get_score());
}
// static
bool
coot::ligand::compare_scored_ligands_using_correlation(const std::pair<coot::minimol::molecule, ligand_score_card> &sl_1,
                                                       const std::pair<coot::minimol::molecule, ligand_score_card> &sl_2) {


   if (sl_1.second.correlation.first && sl_2.second.correlation.first)
      return (sl_1.second.correlation < sl_2.second.correlation);
   return (sl_1.second.get_score() < sl_2.second.get_score());
}

unsigned int
coot::ligand::n_ligands_for_cluster(unsigned int iclust,
                                    float frac_limit_of_peak_score) const {

   unsigned int n = 0;
   float top_score = -1;

   if (final_ligand[iclust].size() > 0) {
      top_score = final_ligand[iclust][0].second.get_score();
      for (unsigned int i=0; i<final_ligand[iclust].size(); i++) {
         if (final_ligand[iclust][i].second.get_score() > frac_limit_of_peak_score * top_score)
            n++;
      }
   }
   std::cout << "debug:: n_ligands_for_cluster() top_score " << top_score << " and "
             << n << " are decent out of " << final_ligand[iclust].size()
             << std::endl;
   return n;
}

unsigned int
coot::ligand::n_ligands_for_cluster(unsigned int iclust) const {

   unsigned int n = 0;
   if (iclust < final_ligand.size())
      n = final_ligand[iclust].size();
   return n;
}




// generate correlation scores for the top n_sol solutions and re-sort
//
void
coot::ligand::score_and_resort_using_correlation(unsigned int iclust, unsigned int n_sol) {

   bool debug = false;

//    #pragma omp parallel for
//    for (unsigned int i=0; i<final_ligand[iclust].size(); i++) {
//       usleep(int(100000 * float(util::random())/float(RAND_MAX)));
//       std::cout << "   parallel i " << i << std::endl;
//    }

   unsigned int n_ligs = final_ligand[iclust].size();


   if (debug)
      std::cout << "score_and_resort_using_correlation iclust: " << iclust << " n_ligs " << n_ligs
                << " n_sol " << n_sol << std::endl;

//    #pragma omp parallel for
   for (unsigned int i=0; i<n_ligs; i++) {
      if (i < n_sol) {

         const minimol::molecule &lig_mol = final_ligand[iclust][i].first;
         mmdb::Manager *mol = lig_mol.pcmmdbmanager(); // d
         std::vector<residue_spec_t> specs;
         residue_spec_t spec(lig_mol[0].fragment_id,
                             lig_mol[0].min_res_no(), "");
         specs.push_back(spec);
         short int mode = 0; // all atoms
         std::vector<residue_spec_t> neighb_specs; // Dummy value currently.
                                                   // Don't count grid points of the spec residues
                                                   // that are part of other residues in
                                                   // the region.

         double c = util::map_to_model_correlation(mol, specs, neighb_specs,
                                                   mode, 1.5, xmap_pristine);
         if (debug)
            std::cout << "----- in get_correl() constructed spec for i "
                      << i << " " << spec
                      << " which has correlation " << c << std::endl;

         std::pair<bool, double> p(true, c);
         final_ligand[iclust][i].second.correlation = p;
         delete mol;
      }
   }

   std::sort(final_ligand[iclust].begin(),
             final_ligand[iclust].end(),
             compare_scored_ligands_using_correlation);
   std::reverse(final_ligand[iclust].begin(),
                final_ligand[iclust].end());


   if (debug) {
      std::cout << "INFO post-sort: ------------------ iclust: " << iclust
                <<  " size: " << final_ligand[iclust].size()
                << " solutions for cluster "
                << iclust << " ------------- " << std::endl;
      for (unsigned int isol=0; isol<final_ligand[iclust].size(); isol++)
         std::cout << "   post correl " << isol << " of "
                   << final_ligand[iclust].size() << " "
                   << final_ligand[iclust][isol].second << std::endl;
   }

}

double
coot::ligand::get_correl(const minimol::molecule &lig_mol) const {

   mmdb::Manager *mol = lig_mol.pcmmdbmanager();
   std::vector<residue_spec_t> specs;
   residue_spec_t spec(lig_mol[0].fragment_id,
                       lig_mol[0].min_res_no(), "");
   specs.push_back(spec);
   std::vector<residue_spec_t> neighb_specs; // Dummy (empty) value currently.
                                             // Don't count grid points of the spec residues
                                             // that are part of other residues in
                                             // the region.
   short int mode = 0; // all atoms
   double c = util::map_to_model_correlation(mol, specs, neighb_specs, mode, 1.5, xmap_pristine);
   if (0)
      std::cout << "----- in get_correl() constructed spec " << spec
                << " which has correlation " << c << std::endl;
   delete mol;
   return c;
}

// this should only be run post-sort (post-correlation sort)
void
coot::ligand::limit_solutions(unsigned int iclust,
                              float frac_max_correl_lim,
                              int max_n_solutions,
                              float tolerance,
                              bool filter_by_torsion_match) { // false

   bool debug = false;
   if (final_ligand[iclust].size()) {
      float min_correl = final_ligand[iclust][0].second.correlation.second * frac_max_correl_lim;
      if (debug)
         std::cout << "INFO:: ..... in limit_solutions() min_correl is " << min_correl << std::endl;
      final_ligand[iclust].erase(std::remove_if(final_ligand[iclust].begin(),
                                                final_ligand[iclust].end(),
                                                scored_ligand_eraser(min_correl)),
                                 final_ligand[iclust].end());
   }

   if (filter_by_torsion_match) {
      // make a vector of final solution mmdb::Residues:
      std::vector<std::pair<mmdb::Residue *, mmdb::Manager *>  >
         final_solution_residues(final_ligand[iclust].size());
      for (unsigned int i=0; i<final_ligand[iclust].size(); i++) {
         mmdb::Manager *mol = get_solution(iclust, i).pcmmdbmanager();
         if (mol) {
            mmdb::Residue *r = util::get_first_residue(mol);
            if (r) {
               std::pair<mmdb::Residue *, mmdb::Manager *> p(r, mol);
               final_solution_residues[i] = p;
            }
         }
      }

      // now we have the vector, let's test for similar
   }

   if (debug)
      for (unsigned int isol=0; isol<final_ligand[iclust].size(); isol++)
         std::cout << "limit solutions: " << isol << " of "
                   << final_ligand[iclust].size() << " "
                   << final_ligand[iclust][isol].second << std::endl;
}





// tinker with mmmol
void
coot::ligand::set_cell_and_symm(coot::minimol::molecule *mmmol) const {

   float a[6];
   a[0] = xmap_pristine.cell().descr().a();
   a[1] = xmap_pristine.cell().descr().b();
   a[2] = xmap_pristine.cell().descr().c();
   a[3] = clipper::Util::rad2d(xmap_pristine.cell().descr().alpha());
   a[4] = clipper::Util::rad2d(xmap_pristine.cell().descr().beta());
   a[5] = clipper::Util::rad2d(xmap_pristine.cell().descr().gamma());


   mmmol->set_cell(a);
   mmmol->set_spacegroup(xmap_pristine.spacegroup().symbol_hm().c_str());

}

short int
coot::ligand::similar_eigen_values(int iclust, int ilig) const {

   float fac = 0.3;
   short int i = 1;

   std::cout << "comparing eigens: " << std::endl;
   for (int ii=0; ii<3; ii++) {
      std::cout << initial_ligand_eigenvalues[ilig][ii] << " "
                << sqrt(cluster[iclust].eigenvalues[ii]) << std::endl;
   }

   for (int ii=0; ii<3; ii++) {
      if (initial_ligand_eigenvalues[ilig][ii] > (1.0+fac)*sqrt(cluster[iclust].eigenvalues[ii]) ||
          initial_ligand_eigenvalues[ilig][ii] < (1.0-fac)*sqrt(cluster[iclust].eigenvalues[ii])) {
         return 0;
      }
   }
   std::cout << std::endl;

   return i;
}

// ior is the orientation index, if it is negative, we are writing out
// the best ligand for that site, and we don't need the orientation
// index (ior) in the filename.
std::string
coot::ligand::ligand_filename(int n_count, int ior) const {

   std::string outfile("ligand-");
   if (ior >= 0) {
      outfile += util::int_to_string(ior);
   }
   outfile += ".pdb";
   if (ior < 0)
      outfile = "best-orientation-" + outfile;
   return outfile;
}


std::string
coot::ligand::get_first_residue_name(const coot::minimol::molecule &mol) const {

   std::string name("");

   for (unsigned int ifrag=0; ifrag<mol.fragments.size(); ifrag++) {
      for (int ires=mol.fragments[ifrag].min_res_no(); ires<=mol.fragments[ifrag].max_residue_number(); ires++) {
         name = mol[ifrag][ires].name;
         if (name != "")
            break;
      }
      if (name != "")
         break;
   }
   return name;
}

// a rough and ready calculator, so that we don't try to fit tris to a
// polysaccharide or vica versa.
//
// The volume of the ligand has to be at least a little bit similar to
// the cluster of electron density points.
//
// c.f. cluster_is_possible_water
//
short int
coot::ligand::cluster_ligand_size_match(int iclust, int ilig) {

   float vol = xmap_pristine.cell().volume();
   float ngrid =
      xmap_pristine.grid_sampling().nu() *
      xmap_pristine.grid_sampling().nv() *
      xmap_pristine.grid_sampling().nw();

   float grid_vol = vol/ngrid;
   float cluster_vol = grid_vol * cluster[iclust].map_grid.size();

//     std::cout << "INFO: grid_vol " << grid_vol << " cluster_vol "
//               << cluster_vol << std::endl;

   std::vector<minimol::atom *> atoms = initial_ligand[ilig].select_atoms_serial();
   int n_lig_atoms = 0;
   for (unsigned int i=0; i<atoms.size(); i++) {
      if (atoms[i]->element != " H") {
         n_lig_atoms++;
      }
   }

   // 1.2^3 = 1.78
   float ligand_vol = (4.0*3.14159/3.0)*1.78*float(n_lig_atoms);

//     std::cout << "cluster_ligand_size_match: "
//               << cluster[iclust].map_grid.size()
//               << " " << cluster_vol << " vs " << ligand_vol << std::endl;

   // Check whether ligand_vol and cluster_vol are "more or less" the
   // same.  Realise that ligand_vol is quite possible to be bigger
   // than cluster_vol for reasonable density and reasonably contour
   // level for a matching ligand.
   //
   if ( (ligand_vol/cluster_vol) < 7.0 &&  // limit for big ligand vs small density
        (ligand_vol/cluster_vol) > 0.8 )  // limit for small ligand vs bit density
      {
      return 1;
   } else {
      return 0;
   }
}

// Return the score.  The score is the sum of the densities at the
// atom centres of the atoms in the ligand after (rigid body)
// refinement.
//
// Note that we use xmap_pristine to score.
//
coot::ligand_score_card
coot::ligand::fit_ligand_copy(int iclust, int ilig, int ior) {

   // checkme - is fitted_ligand_vec[ilig] and
   // fitted_ligand_vec[ilig][iclust] properly dimensioned before we
   // get here?
   //
   //
   if (iclust >= int(fitted_ligand_vec[ilig].size())) {
      fitted_ligand_vec[ilig].resize(iclust + 1);
   }
   fitted_ligand_vec[ilig][iclust] = initial_ligand[ilig]; // Oh joy!  The relief is
                                                            // ... substantial.

   std::vector<minimol::atom *> atoms_p = fitted_ligand_vec[ilig][iclust].select_atoms_serial();

//    std::cout << "DEBUG:: There are "<< atoms_p.size() << " atoms in "
//              << "atoms_p in fit_ligand_copy" << std::endl;

   // First move the ligand to the site of the cluster:
   for(unsigned int ii=0; ii<atoms_p.size(); ii++) {
      atoms_p[ii]->pos = transform_ligand_atom(atoms_p[ii]->pos,ilig,iclust,ior);
   }

   // rigid_body_refine_ligand(&atoms_p, xmap_pristine); // ("rigid body") move atoms.
   // Joel suggests that we have another go at the masked map:
   rigid_body_refine_ligand(&atoms_p, xmap_masked, xmap_pristine, rotation_component, gradient_scale); // ("rigid body") move atoms.

   // update_fitted_ligand_vec(i, atoms); // no need now.  We will
   // reintroduce this when I convert to RTops to do the
   // manipulations.  We will need to convert the atoms for final
   // output.
   // ligand_score_card s = score_orientation(atoms_p, xmap_pristine); // JB says try masked map
   float ff = 0.1; // fit fraction
   ligand_score_card s = score_orientation(atoms_p, xmap_masked, ff);
   s.set_ligand_number(ilig);
   return s;
}

// Return the score.  The score is the sum of the densities at the
// atom centres of the atoms in the ligand after (rigid body)
// refinement.
//
// Note that we use xmap_pristine to score.
//
coot::ligand_score_card
coot::ligand::fit_ligand_copy(int iclust, int ilig, int ior, const clipper::RTop_orth &eigen_ori) {

   // checkme - is fitted_ligand_vec[ilig] and
   // fitted_ligand_vec[ilig][iclust] properly dimensioned before we
   // get here?
   //
   //
   if (iclust >= int(fitted_ligand_vec[ilig].size())) {
      fitted_ligand_vec[ilig].resize(iclust + 1);
   }
   fitted_ligand_vec[ilig][iclust] = initial_ligand[ilig]; // Oh joy!  The relief is
                                                            // ... substantial.

   std::vector<minimol::atom *> atoms_p = fitted_ligand_vec[ilig][iclust].select_atoms_serial();

   // First move the ligand to the site of the cluster:
   for(unsigned int ii=0; ii<atoms_p.size(); ii++) {
      atoms_p[ii]->pos = transform_ligand_atom(atoms_p[ii]->pos,ilig,iclust,ior, eigen_ori);
   }

   // rigid_body_refine_ligand(&atoms_p, xmap_pristine); // ("rigid body") move atoms.
   // Joel suggests that we have another go at the masked map:
   rigid_body_refine_ligand(&atoms_p, xmap_masked, xmap_pristine, rotation_component, gradient_scale); // ("rigid body") move atoms.

   // update_fitted_ligand_vec(i, atoms); // no need now.  We will
   // reintroduce this when I convert to RTops to do the
   // manipulations.  We will need to convert the atoms for final
   // output.
   // ligand_score_card s = score_orientation(atoms_p, xmap_pristine); // JB says try masked map
   float ff = 0.1;
   ligand_score_card s = score_orientation(atoms_p, xmap_masked, ff);
   // std::cout << "  ligand_score_card: " << s << std::endl;
   s.set_ligand_number(ilig);
   return s;
}


// Return the operator that tranforms the initial ligand to
// the position and orientation of the i_cluster'th cluster.
//
// Get your geometry hat on!
//
// ior is the ith origin rotation (see the constructor).
// i is the i'th cluster
//
clipper::Coord_orth
coot::ligand::transform_ligand_atom(const clipper::Coord_orth &a_in,
                                    int ilig, int iclust, int ior) const {

   return transform_ligand_atom(a_in, ilig,
                                cluster[iclust].eigenvectors_and_centre,
                                ior);
}

clipper::Coord_orth
coot::ligand::transform_ligand_atom(const clipper::Coord_orth &a_in,
                                    int ilig, const clipper::RTop_orth &cluster_rtop, int ior) const {

   clipper::Coord_orth a;

   clipper::RTop_orth ligand_op(initial_ligand_eigenvectors[ilig],
                                initial_ligand_model_centre[ilig]);

   clipper::RTop_orth lopi = clipper::RTop_orth (ligand_op.inverse());

   clipper::RTop_orth origin_rotate_op = clipper::RTop_orth(origin_rotations[ior],
                                                            clipper::Coord_orth(0,0,0));

//    std::cout << std::endl << std::endl << std::endl;
//    std::cout << "ligand_op"    << std::endl << ligand_op.format() << std::endl;
//    std::cout << "eigen_centre" << std::endl
//              << cluster[i].eigenvectors_and_centre.format() << std::endl;

   a = a_in.transform(lopi); // move to origin
   a = a.transform(origin_rotate_op); // rotate round origin
   a = a.transform(cluster_rtop);

   return a;
}





// Return the operator that tranforms the initial ligand to
// the position and orientation of the i_cluster'th cluster.
//
// Get your geometry hat on!
//
// ior is the ith origin rotation (see the constructor).
// i is the i'th cluster
//
clipper::Coord_orth
coot::ligand::transform_ligand_atom(const clipper::Coord_orth &a_in,
                                    int ilig, int iclust, int ior,
                                    const clipper::RTop_orth &eigen_ori) const {

   clipper::Coord_orth a;

   clipper::RTop_orth ligand_op(initial_ligand_eigenvectors[ilig],
                                initial_ligand_model_centre[ilig]);

   clipper::RTop_orth lopi = clipper::RTop_orth (ligand_op.inverse());

   clipper::RTop_orth origin_rotate_op = clipper::RTop_orth(origin_rotations[ior],
                                                            clipper::Coord_orth(0,0,0));

//    std::cout << std::endl << std::endl << std::endl;
//    std::cout << "ligand_op"    << std::endl << ligand_op.format() << std::endl;
//    std::cout << "eigen_centre" << std::endl
//              << cluster[i].eigenvectors_and_centre.format() << std::endl;

   a = a_in.transform(lopi); // move to origin
   a = a.transform(eigen_ori);
   a = a.transform(origin_rotate_op); // rotate round origin
   a = a.transform(cluster[iclust].eigenvectors_and_centre);

   return a;
}


// Construct a new mmdb(using space group and cell)
// add a empty model
//

// Things to add:
//
// Rotating the ligand so that eigenvalues match.
//
// Try all 4 different orientations and refine each
//
//   Refine with the system discussed with Kevin
//
// Apply symmetry to the ligand centres so that they lie "next to"
// The protein that was the mask.
//

// Let's concentrate on ligand issue (not solvent) tonight (20030323)
//

// Move the atoms by rigid body to fit the xmap.
//
// static
void
coot::ligand::rigid_body_refine_ligand(std::vector<minimol::atom *> *atoms_p,
                                       const clipper::Xmap<float> &xmap_fitting,
                                       const clipper::Xmap<float> &xmap_pristine,
                                       const clipper::RTop_orth rotation_component[3],
                                       float gradient_scale) {

   int n_atoms = atoms_p->size();
   int round_max = 300;
   int iround = 0;
   double move_by_length = 1.0; // just to past test initially
   double angle_sum = 1.0;      // as above
   bool debug = false;

   auto mean_ligand_position = [] (const std::vector<minimol::atom *> &atoms) {
                                  clipper::Coord_orth p(0,0,0);
                                  for (unsigned int ii=0; ii<atoms.size(); ii++)
                                     p += atoms[ii]->pos;
                                  double sc = 1.0/static_cast<double>(atoms.size());
                                  return clipper::Coord_orth(sc * p);
                               };

   auto apply_angles_to_ligand = [] (const clipper::Vec3<double> &angles ,
                                     const clipper::Coord_orth &mean_pos,
                                     const std::vector<minimol::atom *> *atoms_p) {

                                    double scale_factor = 1.0;
                                    unsigned int n_atoms = atoms_p->size();
                                    if (n_atoms > 2)
                                       scale_factor = 1.5/sqrt(static_cast<double>(n_atoms));

                                    double sin_t = sin(-scale_factor*angles[0]);
                                    double cos_t = cos(-scale_factor*angles[0]);
                                    clipper::Mat33<double> x_mat(1,0,0, 0,cos_t,-sin_t, 0,sin_t,cos_t);

                                    sin_t = sin(-scale_factor*angles[1]);
                                    cos_t = cos(-scale_factor*angles[1]);
                                    clipper::Mat33<double> y_mat(cos_t,0,sin_t, 0,1,0, -sin_t,0,cos_t);

                                    sin_t = sin(-scale_factor*angles[2]);
                                    cos_t = cos(-scale_factor*angles[2]);
                                    clipper::Mat33<double> z_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);

                                    clipper::Mat33<double> angle_mat = x_mat * y_mat * z_mat;
                                    clipper::RTop_orth angle_op(angle_mat, clipper::Coord_frac(0,0,0));

                                    for (unsigned int ii=0; ii<atoms_p->size(); ii++) {
                                       (*atoms_p)[ii]->pos -= mean_pos;
                                       (*atoms_p)[ii]->pos = (*atoms_p)[ii]->pos.transform(angle_op);
                                       (*atoms_p)[ii]->pos += mean_pos;
                                    }
                                 };


   while ((iround < round_max) &&
          ((move_by_length > 0.002) || (angle_sum > 0.002))) {

      clipper::Coord_orth midpoint(0,0,0);
      float ff = 0.1; // doesn't matter, we don't want warning messages

      if (debug)
         std::cout << "---------------------  start of the rigid body round " << iround << "\n"
                   << "   score is " << score_orientation(*atoms_p, xmap_fitting, ff) << std::endl;

      for (unsigned int ii=0; ii<atoms_p->size(); ii++) {
         midpoint += (*atoms_p)[ii]->pos;
      }

      double scale = 1.0/double(atoms_p->size());
      clipper::Coord_orth mean_pos = clipper::Coord_orth(midpoint.x() * scale,
                                                         midpoint.y() * scale,
                                                         midpoint.z() * scale);

      if (debug)
         std::cout << "   iround " << iround << " real mean pos mid point in rigid_body: point 0"
                   << mean_pos.format() << std::endl;

      // Get the average gradient in orthogonal x, y, z directions:
      //
      clipper::Grad_map<float> grad;
      clipper::Grad_frac<float> grad_frac;
      clipper::Grad_orth<float> grad_orth;
      std::vector<clipper::Grad_orth<float> > grad_vec(n_atoms);
      float dv, sum_dx = 0, sum_dy = 0, sum_dz = 0;
      for (int ii=0; ii<n_atoms; ii++) {
         const clipper::Coord_orth &atom_pos = (*atoms_p)[ii]->pos;
         clipper::Coord_frac atom_pos_frc = atom_pos.coord_frac(xmap_pristine.cell());
         clipper::Coord_map  atom_pos_map = atom_pos_frc.coord_map(xmap_pristine.grid_sampling());
         clipper::Interp_cubic::interp_grad(xmap_fitting, atom_pos_map, dv, grad);
         grad_frac = grad.grad_frac(xmap_pristine.grid_sampling());
         grad_orth = grad_frac.grad_orth(xmap_pristine.cell());
         // std::cout << "gradients: " << grad_orth.format() << std::endl;
         sum_dx += grad_orth.dx();
         sum_dy += grad_orth.dy();
         sum_dz += grad_orth.dz();
         grad_vec[ii] = grad_orth;
      }
      double datfrac = 1.0/double (n_atoms);
      clipper::Grad_orth<float> av_grad(sum_dx * datfrac,
                                        sum_dy * datfrac,
                                        sum_dz * datfrac);

      clipper::Coord_orth moved_by = clipper::Coord_orth(gradient_scale*av_grad.dx(),
                                                         gradient_scale*av_grad.dy(),
                                                         gradient_scale*av_grad.dz());

      if (debug)
         std::cout << "   iround " << iround << " moving by " << moved_by.format() << std::endl;

      clipper::Vec3<double> angles;
      bool angles_are_valid = false;
      if (atoms_p->size() > 1) {
         angles = get_rigid_body_angle_components(*atoms_p, mean_pos, grad_vec,
                                                  rotation_component, gradient_scale);
         angles_are_valid = true;
      }

      move_by_length = sqrt(moved_by.lengthsq());

      if (debug) {
         std::cout << "   iround " << iround << " " << move_by_length
                   << " " << moved_by.x()
                   << " " << moved_by.y()
                   << " " << moved_by.z() << "\n";
         std::cout << "      mean pos in rigid body refine ligand: point 1 "
                   << mean_pos.format() << std::endl;
         std::cout << "      moved by: " <<  moved_by.format() << std::endl;
         std::cout << "      mean by function: "
                   << mean_ligand_position(*atoms_p).format() << std::endl;
      }

      mean_pos += moved_by;
      for (int ii=0; ii<n_atoms; ii++)
         (*atoms_p)[ii]->pos += moved_by;

      // Consider moving the generation (not application) of the
      // angles above the application of the translations.
      //
      if (debug)
         std::cout << "   iround " << iround << " Now to apply the angles: "
                   << clipper::Util::rad2d(angles[0]) << " "
                   << clipper::Util::rad2d(angles[1]) << " "
                   << clipper::Util::rad2d(angles[2]) << std::endl;

      if (angles_are_valid) {
         apply_angles_to_ligand(angles,mean_pos, atoms_p);

         // set angle_sum for next round
         angle_sum = 0.0;
         angle_sum += fabs(clipper::Util::rad2d(angles[0]));
         angle_sum += fabs(clipper::Util::rad2d(angles[1]));
         angle_sum += fabs(clipper::Util::rad2d(angles[2]));
         if (debug)
            std::cout << "   iround " << iround << " moved: " << move_by_length
                      << " angle_sum: " << angle_sum << std::endl;
      }

      iround++;

   } // irounds

}

clipper::Coord_orth
coot::ligand::mean_ligand_position(const std::vector<minimol::atom *> &atoms) const {

   clipper::Coord_orth p(0,0,0);
   for (unsigned int ii=0; ii<atoms.size(); ii++) {
           p += atoms[ii]->pos;
   }

   double scale = 1/double(atoms.size());
   clipper::Coord_orth a(p.x() * scale,
                            p.y() * scale,
                            p.z() * scale);

   return a;
}

// Return the rigid body angles \alpha_x, \alpha_y, \alpha_z in
// radians.
//
// We will need to turn them into a RTop_orth later.
//
clipper::Vec3<double>
coot::ligand::get_rigid_body_angle_components(const std::vector<minimol::atom *> &atoms,
                                              const clipper::Coord_orth &mean_pos,
                                              const std::vector<clipper::Grad_orth<float> > &grad_vec,
                                              const clipper::RTop_orth rotation_component[3],
                                              float gradient_scale) {

   clipper::Vec3<double> a;  // return this
   bool debug = false;

   if (debug)
      std::cout << "DEBUG:: rigid_body: gradient_scale is " << gradient_scale << std::endl;

   // For each atom, we multiply the vector between the atom position
   // and the centre of rotation (V) by rotation_component to get the
   // vector pependicular to this vector in the appropriate plane
   // (e.g. XY plane for rotation round Z) - let's call that new vector
   // Vp.  We get the unit vector of Vp: Vpu.  We form the dot product
   // of the vector of gradients with Vpu.  This gives us a distance
   // (+ve or -ve) d_Vpu.
   //
   // We want an angle however, and that is atan2(d_Vpu,|V|).
   // Average that angle for contributions from all the atoms.
   // Do this for \alpha_x, \alpha_y, \alpha_z.

   double sum_alpha[3];
   double sum_cos_theta[3];
   double sum_grad[3];
   double Vp_rms_sum[3];

   for (int ir=0; ir<3; ir++) {
      sum_grad[ir]      = 0.0;
      sum_alpha[ir]     = 0.0;
      Vp_rms_sum[ir]    = 0.0;
      sum_cos_theta[ir] = 0.0; // debugging
   }

   clipper::Coord_orth V, Vp[3];
   double dot_prod[3];
   for (unsigned int ii=0; ii<atoms.size(); ii++) {
      V = atoms[ii]->pos - mean_pos; // vector from centre to atom
      clipper::Coord_orth grad(grad_vec[ii].dx(),
                               grad_vec[ii].dy(),
                               grad_vec[ii].dz());
      for (int ir=0; ir<3; ir++) {  // rotation axis
         Vp[ir] = V.transform(rotation_component[ir]);
         Vp_rms_sum[ir] += Vp[ir].lengthsq();
         dot_prod[ir] = clipper::Coord_orth::dot(grad,Vp[ir]);
         sum_grad[ir] += dot_prod[ir];
         if (debug)
            std::cout << "   iat: " << ii << " V(ec) " << V.format()
                      << " grad: " << grad.format() << " * " << Vp[ir].format() << " is "
                      << dot_prod[ir] << " now sum_grad[" << ir << "] = " << sum_grad[ir] << std::endl;
      }
   }

   double Vp_av_len[3];
   for (int ir=0; ir<3; ir++) {  // rotation axis
      Vp_av_len[ir] = sqrt(Vp_rms_sum[ir]/double(atoms.size()));
   }

   for (int ir=0; ir<3; ir++) {

      // a[ir] = gradient_scale * sum_grad[ir]/Vp_av_len[ir];
      a[ir] = gradient_scale * 0.1 * sum_grad[ir]/(Vp_av_len[ir] * sqrt(static_cast<double>(atoms.size())));

      if (debug) {
         std::cout << "  a[" << ir << "] is " << sum_grad[ir] << "/" << Vp_av_len[ir] << "/sqrt("
                   << atoms.size() << ") = " << a[ir] << "     " << a[ir] * 57.3 << " degrees " << std::endl;
          std::cout << "Vp_av_len[" << ir << "] is " << Vp_av_len[ir] << "    and sum_grad["
                    << ir << "] is " << sum_grad[ir] << "  ";
          std::cout << "  a[" << ir << "] is " << a[ir]*57.3 << " degrees " << std::endl;
      }
   }

   return a;
}

void
coot::ligand::apply_angles_to_ligand(const clipper::Vec3<double> &angles ,
                                     const std::vector<minimol::atom *> *atoms_p,
                                     const clipper::Coord_orth &mean_pos) {

   double sin_t;
   double cos_t;
   double scale_factor = 1.0;
   unsigned int n_atoms = atoms_p->size();
   if (n_atoms > 2)
      scale_factor = 1.5/double(sqrt(n_atoms));

   sin_t = sin(-scale_factor*angles[0]);
   cos_t = cos(-scale_factor*angles[0]);
   clipper::Mat33<double> x_mat(1,0,0, 0,cos_t,-sin_t, 0,sin_t,cos_t);

   sin_t = sin(-scale_factor*angles[1]);
   cos_t = cos(-scale_factor*angles[1]);
   clipper::Mat33<double> y_mat(cos_t,0,sin_t, 0,1,0, -sin_t,0,cos_t);

   sin_t = sin(-scale_factor*angles[2]);
   cos_t = cos(-scale_factor*angles[2]);
   clipper::Mat33<double> z_mat(cos_t,-sin_t,0, sin_t,cos_t,0, 0,0,1);

   clipper::Mat33<double> angle_mat = x_mat * y_mat * z_mat;
   clipper::RTop_orth angle_op(angle_mat, clipper::Coord_frac(0,0,0));

   for (unsigned int ii=0; ii<atoms_p->size(); ii++) {
      (*atoms_p)[ii]->pos -= mean_pos;
      (*atoms_p)[ii]->pos = (*atoms_p)[ii]->pos.transform(angle_op);
      (*atoms_p)[ii]->pos += mean_pos;
   }
}

// return the score
//
// static
coot::ligand_score_card
coot::ligand::score_orientation(const std::vector<minimol::atom *> &atoms,
                                const clipper::Xmap<float> &xmap_fitting,
                                float fit_fraction,
                                bool use_linear_interpolation) {

   coot::ligand_score_card score_card;
   int n_positive_atoms = 0;

   int n_non_hydrogens = 0;

   for (unsigned int ii=0; ii<atoms.size(); ii++) {
      const clipper::Coord_orth &atom_pos = atoms[ii]->pos;
      clipper::Coord_frac atom_pos_frc = atom_pos.coord_frac(xmap_fitting.cell());
      if (!atoms[ii]->is_hydrogen_p()) {
         // float dv = xmap_fitting.interp<clipper::Interp_cubic>(atom_pos_frc);
         float dv = 0;
         if (use_linear_interpolation)
            dv = xmap_fitting.interp<clipper::Interp_linear>(atom_pos_frc);
         else
            dv = xmap_fitting.interp<clipper::Interp_cubic>(atom_pos_frc); // faster and accurate enough
         score_card.add(dv * atoms[ii]->occupancy);
         n_non_hydrogens++;
         if (dv > 0)
            n_positive_atoms++;
      }
   }

   if (atoms.size() > 0) {
      // fit_fraction is initially 0.75, but can be changed by an public member
      // function.
      if (n_non_hydrogens > 0) {
         score_card.set_n_ligand_atoms(n_non_hydrogens);
         if (0)
            std::cout << "fit fraction test: is " << n_positive_atoms << "/"
                      << n_non_hydrogens << " ("
                      << float(n_positive_atoms)/float(n_non_hydrogens)
                      << ") < " << fit_fraction << std::endl;
         if (float(n_positive_atoms)/float(n_non_hydrogens) >= fit_fraction ) { // arbitary
            score_card.many_atoms_fit = 1; // consider using a member function
            score_card.score_per_atom = score_card.get_score()/float(n_non_hydrogens);
         } else {
            if (false) // too noisy
               std::cout << "WARNING:: badly fitting atoms, failing fit_fraction test "
                         << n_positive_atoms << " / " << n_non_hydrogens << " vs " << fit_fraction
                         << std::endl;
         }
      } else {
         // Pathalogical case.  No non-hydrogens in ligand.  This code
         // should never realistically be run...
         score_card.many_atoms_fit = 0;
         score_card.score_per_atom = -1.0;
      }

//       std::cout << "for score card score: "  << score_card.score << std::endl;
//       for (int i=0; i< atoms.size(); i++) {

//          clipper::Coord_orth atom_pos(atoms[i]->pos.x(),
//                                       atoms[i]->pos.y(),
//                                       atoms[i]->pos.z());
//          clipper::Coord_frac atom_pos_frc = atom_pos.coord_frac(xmap.cell());
//          dv = xmap_fitting.interp<clipper::Interp_cubic>(atom_pos_frc);

//          std::cout << i << " " << " " << atom_pos.format() << " " << dv << std::endl;
//       }
   }
   return score_card;
}

float
coot::ligand::score_position(const clipper::Coord_orth &atom_pos,
                             const clipper::Xmap<float> &xmap_fitting) const {

   clipper::Coord_frac atom_pos_frc = atom_pos.coord_frac(xmap_fitting.cell());
   return xmap_fitting.interp<clipper::Interp_cubic>(atom_pos_frc);
}



short int
coot::ligand::masking_molecule_has_atoms() const {

   short int i;

   if (protein_atoms.is_empty())
      i = 0;
   else
      i = 1;

   return i;
}


// This list of peaks above a particular level (and one day
// optionally below the same but negative level):
//
std::vector<std::pair<clipper::Coord_orth, float> >
coot::ligand::cluster_centres() const {

   std::vector<std::pair<clipper::Coord_orth, float> > c;

   clipper::Coord_orth pos;
   float d;
   for (unsigned int ic=0; ic<cluster.size(); ic++) {

      clipper::Coord_orth pos_l(cluster[ic].eigenvectors_and_centre.trn());
      d = density_at_point(pos_l, xmap_cluster);
      std::pair<clipper::Coord_orth, float> p(pos_l, d);
      c.push_back(p);
   }

   return c;
}


void
coot::ligand::import_reference_coords(mmdb::Manager *mol) {

}



/* Flip the ligand (usually active residue) around its eigen
   vectors to the next flip number.  Immediate replacement (like
   flip peptide). We need to undo the current flip number first
   though (if flip_number is not 1). */
coot::minimol::molecule
coot::ligand::flip_ligand(short int flip_number) const {

   // the unflips rotate (back) 180 degrees
   int unflip_number = flip_number -1;
   if (unflip_number == -1)
      unflip_number = 3;

   coot::minimol::molecule m = initial_ligand[0];

   std::vector<minimol::atom *> atoms_p = m.select_atoms_serial();

   // std::cout << "DEBUG:: flip_number is "  << flip_number << std::endl;
   for (unsigned int ii=0; ii<atoms_p.size(); ii++) {
      clipper::Coord_orth lig_centre = initial_ligand_model_centre[0];
      clipper::Mat33<double> lig_eigenvectors = initial_ligand_eigenvectors[0];
      clipper::RTop_orth cl_op(lig_eigenvectors, lig_centre);
      atoms_p[ii]->pos = transform_ligand_atom(atoms_p[ii]->pos,0,cl_op,unflip_number);
      atoms_p[ii]->pos = transform_ligand_atom(atoms_p[ii]->pos,0,cl_op,flip_number);
   }

   return m;
}



std::vector<std::pair<std::string, clipper::Xmap<float> > >
coot::ligand::make_masked_maps_split_by_chain(mmdb::Manager *mol) {

   std::vector<std::pair<std::string, clipper::Xmap<float> > > v;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      if (n_chains > 0) {
         clipper::Xmap<float> save_map = xmap_masked;
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int sel_hnd = mol->NewSelection(); // d
            std::string chain_id(chain_p->GetChainID());
            // is there a better way to select the atoms of this chain?
            std::string atom_selection_str = "/1/" +  chain_id;
            mol->Select(sel_hnd, mmdb::STYPE_ATOM,
               atom_selection_str.c_str(), mmdb::SKEY_NEW);
            mask_map(mol, sel_hnd, true); // change xmap_cluster
            mol->DeleteSelection(sel_hnd);
            std::string name = "Masked Map for Chain " + chain_id;
            v.push_back(std::pair<std::string, clipper::Xmap<float> > (name, xmap_masked));
            xmap_cluster = save_map;
         }
      }
   }
   return v;

}
