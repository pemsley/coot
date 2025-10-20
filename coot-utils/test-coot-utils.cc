/* coot-utils/test-utils.cc
 *
 * Copyright 2005, 2006 by Paul Emsley, The University of York
 * Copyright 2014 by Medical Research Council
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

#include <iostream>
#include <algorithm>
#include <chrono>

#include "clipper/core/rotation.h"

#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"
#include "lsq-improve.hh"
#include "helix-analysis.hh"
#include "glyco-tree.hh"


class testing_data {
public:
   static coot::protein_geometry geom;
   testing_data() {
      if (geom.size() == 0)
         geom.init_standard();
   }
};

coot::protein_geometry testing_data::geom;


namespace coot {
   class SortableChainsManager : public mmdb::Manager {
   public:
      void SortChains();
   };
}

void
coot::SortableChainsManager::SortChains() {

   for (int imod=1; imod<=GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      std::vector<std::pair<mmdb::Chain *, std::string> > chain_ids(nchains);
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         std::string chain_id = chain_p->GetChainID();
         chain_ids[ichain] = std::pair<mmdb::Chain *, std::string> (chain_p, chain_id);
      }
      // now chain_ids is full
      std::sort(chain_ids.begin(), chain_ids.end(), sort_chains_util);
      for (int ichain=0; ichain<nchains; ichain++) {
         // model_p->Chain[ichain] = chain_ids[ichain].first;
      }
   }
   PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   FinishStructEdit();
}

void test_euler_angles() {

   clipper::Euler_ccp4 e(M_PI/2.0, 0, 0);
   clipper::Rotation r(e);

   std::cout << "Rotation from euler angles: \n"
             << r.matrix().format() << std::endl;

}


void split_test(const std::string &r) {

   std::vector<std::string> v = coot::util::split_string(r, " ");
   std::cout << ":" << r << ": -> ";
   for (unsigned int i=0; i<v.size(); i++) {
      std::cout << ":" << v[i] << ": ";
   }
   std::cout << "\n";
}

bool
test_quaternion_matrix(clipper::Mat33<double> m){

   coot::util::quaternion q1(m);
   clipper::Mat33<double> m2 = q1.matrix();
   coot::util::quaternion q2(m2);
   bool match = q1.is_similar_p(q2);

   std::cout << "   " << q1 << "\n" << "   " << q2 << "   ";
   std::cout << "match: " << match << std::endl;
   return match;
}

bool test_quaternion_quaternion(coot::util::quaternion q) {

   clipper::Mat33<double> m2 = q.matrix();
   coot::util::quaternion q2(m2);
   bool match = q.is_similar_p(q2);


   std::cout << std::endl;
   std::cout << m2.format() << std::endl;
   std::cout << "   " << q << "\n" << "   " << q2 << "   ";
   std::cout << "match: " << match << std::endl;
   return match;

}

void test_sort_chains() {

   std::string file_name = "test-sort-chains.pdb";
   mmdb::Manager *mol = new mmdb::Manager;
   mol->ReadCoorFile(file_name.c_str());
   coot::sort_chains(mol);
   mol->WritePDBASCII("test-sort-chains-sorted.pdb");
}

void test_lsq_improve() {

   std::cout << "========================= lsq-improve ==================" << std::endl;
   mmdb::Manager *mol_1 = new mmdb::Manager;
   mmdb::Manager *mol_2 = new mmdb::Manager;

   std::string ref_pdb = "tutorial-modern.pdb";
   std::string mov_pdb = "1py3-matched-A6-13.pdb";

   mmdb::ERROR_CODE err_1 = mol_1->ReadCoorFile(ref_pdb.c_str());
   mmdb::ERROR_CODE err_2 = mol_2->ReadCoorFile(mov_pdb.c_str());

   if (err_1) {
      std::cout << "There was an error reading " << ref_pdb << ".\n";
   }
   if (err_2) {
      std::cout << "There was an error reading " << mov_pdb << ".\n";
   }

   for (unsigned int i=0; i<1; i++) {

      try {
         coot::lsq_improve lsq_imp(mol_1, "//A/1-50", mol_2, "//A/4-50");
         lsq_imp.improve();
         clipper::RTop_orth rtop = lsq_imp.rtop_of_moving();
         std::cout << "rtop:\n" << rtop.format() << std::endl;
         coot::util::transform_mol(mol_2, rtop);
         // mol_2->WritePDBASCII("lsq-improved.pdb");

      }
      catch (const std::runtime_error &rte) {
         std::cout << "lsq_improve ERROR::" << rte.what() << std::endl;
      }
   }

}

int test_string_manipulation() {


   /*
   std::string s;
   s = "AVasdfasdfC";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
             << std::endl;
   s = "   AVC";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
             << std::endl;
   s = " AVC ";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
             << std::endl;
   s = "C";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
             << std::endl;
   s = "";
   std::cout << s << " cuts to :" << coot::util::remove_leading_spaces(s) << ":"
             << std::endl;
   */

   std::string a("ABCDefgh");
   std::cout << a << " downcased: " << coot::util::downcase(a) << std::endl;
   std::cout << a << "   upcased: " << coot::util::upcase(a) << std::endl;

   std::string s("Cottage");
   std::string r("tag");
   std::cout  << "removing :" << r << ": from :" << s << ": gives :"
              << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "tage";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
             << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "e";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
             << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
             << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "ball";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
             << coot::util::remove_string(s, r) <<  ":" << std::endl;
   r = "Cottage";
   std::cout << "removing :" << r << ": from :" << s << ": gives :"
             << coot::util::remove_string(s, r) <<  ":" << std::endl;

   r = "Columns";
   split_test(r);
   r = "Columns ";
   split_test(r);
   r = " Columns ";
   split_test(r);
   r = "Columns     of   letters  ";
   split_test(r);
   r = " Columns     of   letters  ";
   split_test(r);

   return 0;

}

int test_matrices() {


   clipper::Mat33<double> m1 (1,0,0, 0,1,0, 0,0,1);
   test_quaternion_matrix(m1);
   clipper::Mat33<double> m2 (0,1,0, 1,0,0, 0,0,-1);
   test_quaternion_matrix(m2);
   // this one from quat-convert.scm:
   clipper::Mat33<double> m3( 0.0347695872187614, 0.773433089256287,   0.632923781871796,
                              0.774806916713715,  0.379149734973907,  -0.505885183811188,
                             -0.631241261959076,  0.507983148097992,  -0.586078405380249);
        // -> (-0.557 -0.694704 -0.0007537 0.454928)
   test_quaternion_matrix(m3);

   coot::util::quaternion q1(1,0,0,0);
   coot::util::quaternion q2(0,1,0,0);
   coot::util::quaternion q3(0,0,1,0);
   coot::util::quaternion q4(0,0,0,1);
   coot::util::quaternion q5(-0.557, -0.694704, -0.0007537, 0.454928);
   test_quaternion_quaternion(q1);
   test_quaternion_quaternion(q2);
   test_quaternion_quaternion(q3);
   test_quaternion_quaternion(q4);
   test_quaternion_quaternion(q5);

   return 0;
}

int test_glyco_tree() {

   testing_data t;
   int dynamic_add_status_1 = t.geom.try_dynamic_add("NAG", 1);
   int dynamic_add_status_2 = t.geom.try_dynamic_add("MAN", 1);
   int dynamic_add_status_3 = t.geom.try_dynamic_add("BMA", 1);
   int dynamic_add_status_4 = t.geom.try_dynamic_add("GAL", 1);
   // int dynamic_add_status_4 = t.geom.try_dynamic_add("GLB", 1); minimal

   mmdb::Manager *mol = new mmdb::Manager;
   // std::string file_name = "3u2s.pdb";
   // coot::residue_spec_t spec("G", 560, "");

   std::string file_name = "sweet2-test-1.pdb";
   coot::residue_spec_t spec("", 1, "");

   mol->ReadCoorFile(file_name.c_str());
   mmdb::Residue *r = coot::util::get_residue(spec, mol);
   if (! r) {
      std::cout << "No residue " << spec << std::endl;
   } else {
      coot::glyco_tree_t gt(r, mol, &t.geom);
   }


   delete mol;
   return 0;
}

int test_helix_analysis() {

   mmdb::Manager *mol = new mmdb::Manager;
   // std::string file_name = "theor-helix-down-z.pdb";
   std::string file_name = "helix-just-off-z.pdb";
   // file_name = "../src/pdb2qc1-sans-sans-ASN.pdb";
   coot::residue_spec_t spec("A", 10, "");
   // coot::residue_spec_t spec("B", 201, "");

   mol->ReadCoorFile(file_name.c_str());
   mmdb::Residue *r = coot::util::get_residue(spec, mol);
   if (! r) {
      std::cout << "No residue " << spec << std::endl;
   } else {
      coot::helix_params_container_t h;
      // h.make(mol, "B", 201, 207);
      h.make(mol, "A", 10, 15);
   }
   delete mol;
   return 0;
}

#include <fstream>

int test_qq_plot() {

   int status = 0;

   std::string file_name = "random-300.tab";

   std::ifstream f(file_name.c_str());

   if (f) {

      std::vector<double> data;
      std::string line;
      while (std::getline(f, line)) {
         std::vector<std::string> bits = coot::util::split_string_no_blanks(line, " ");
         for (unsigned int ibit=0; ibit<bits.size(); ibit++) {
            try {
               double v = coot::util::string_to_float(bits[ibit]);
               data.push_back(v);
            }
            catch (const std::runtime_error &rte) {
               std::cout << "   " << rte.what() << std::endl;
            }
         }
      }

      coot::util::qq_plot_t qq(data);
      std::vector<std::pair<double, double> > qqd = qq.qq_norm();

      for (unsigned int i=0; i<qqd.size(); i++) {
         std::cout << "plot " << i << " " << "   " << qqd[i].first << "   "
                   << qqd[i].second << std::endl;
      }
   }
   return status;
}

#include "coot-least-squares.hh"
int test_least_squares_fit() {

   int status = 0;
   std::vector<std::pair<double, double> > data(3);
   data[0] = std::pair<double, double> (0,-0.1);
   data[1] = std::pair<double, double> (1,2);
   data[2] = std::pair<double, double> (2,4);

   coot::least_squares_fit lsq(data);

   std::cout << "  lsq m " << lsq.m() << std::endl;
   std::cout << "  lsq c " << lsq.c() << std::endl;

   return status;
}

#include "atom-overlaps.hh"

int test_atom_overlaps() {

   int status = 0;


//    testing_data t;
//    t.geom.try_dynamic_add("MG",  1);
//    t.geom.try_dynamic_add("824", 1);
//    t.geom.init_refmac_mon_lib("824-acedrg.cif", 1);

   coot::protein_geometry geom;
   geom.init_standard();
   geom.try_dynamic_add("824", 1);
   geom.try_dynamic_add("MG", 1);

   mmdb::Manager *mol = new mmdb::Manager;
   std::string file_name = "1x8b-H.pdb";
   coot::residue_spec_t spec("A", 901, "");

   int read_status = mol->ReadCoorFile(file_name.c_str());

   if (read_status == mmdb::Error_NoError) {
     mmdb::Residue *residue_p = coot::util::get_residue(spec, mol);
     if (residue_p) {
        int imol = 0;
        std::vector<mmdb::Residue *> neighbs = coot::residues_near_residue(residue_p, mol, 5);
        coot::atom_overlaps_container_t overlaps(residue_p, neighbs, mol, imol, &geom, 0.5, 0.25);
        coot::atom_overlaps_dots_container_t c = overlaps.contact_dots_for_ligand();
     } else {
       std::cout << "Can't find residue" << spec << std::endl;
     }
   } else {
     std::cout << "Failed to read " << file_name << std::endl;
   }
   return status;
}

int test_all_atom_overlaps() {

   int status = 0;
   coot::protein_geometry geom;
   geom.init_standard();
   geom.set_verbose(false);
   geom.try_dynamic_add("824", 1);
   geom.try_dynamic_add("MG", 1);

   // for 5hcj
   geom.try_dynamic_add("LMT", 1);
   geom.try_dynamic_add("MBR", 1);
   geom.try_dynamic_add("CL", 1);
   geom.try_dynamic_add("NA", 1);

   mmdb::Manager *mol = new mmdb::Manager;
   std::string file_name = "1x8b-all-H-no-water.pdb";

   // file_name = "5hcj-with-coot-Hs.pdb";

   int read_status = mol->ReadCoorFile(file_name.c_str());

   if (read_status == mmdb::Error_NoError) {
      // spike length and probe radius (which are not used in overlaps)

      bool ignore_water = false;
      coot::atom_overlaps_container_t overlaps(mol, &geom, ignore_water, 0.5, 0.25);

      overlaps.make_all_atom_overlaps();
      std::vector<coot::atom_overlap_t> olv = overlaps.overlaps;
      std::cout << "Found " << olv.size() << " atom overlaps" << std::endl;
      for (std::size_t ii=0; ii<olv.size(); ii++) {
         const coot::atom_overlap_t &o = olv[ii];
         std::cout << "Overlap " << ii << " "
                   << coot::atom_spec_t(o.atom_1) << " "
                   << coot::atom_spec_t(o.atom_2) << " overlap-vol "
                   << o.overlap_volume << " r_1 "
                   << o.r_1 << " r_2 " << o.r_2 << std::endl;
      }
   }

   delete mol;
   return status;
}

// #include "cp.hh"

int test_cp() {
   mmdb::Manager *mol = new mmdb::Manager;
   std::string file_name = "MAN.pdb";
   mol->ReadCoorFile(file_name.c_str());
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (! model_p) {
      std::cout << "Null model" << std::endl;
   } else {
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         if (! chain_p) {
            std::cout << "Null chain" << std::endl;
         } else {
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (! residue_p) {
                  std::cout << "Null residue" << std::endl;
               } else {
                  /*
                  coot::cp_t cp;
                  double a = cp.amplitude(residue_p);
                  std::vector<double> t(6);
                  t[0] = -62.576; t[1] =  38.474; t[2] =  16.080;
                  t[3] = -56.662; t[4] =  32.146; t[5] =  28.071;
                  // cp.amplitude(t);
                  */
               }
            }
         }
      }
   }
   return 1;
}

#include "reduce.hh"

int test_reduce() {

   mmdb::Manager *mol = new mmdb::Manager;
   std::string file_name = "1x8b.pdb";
   mol->ReadCoorFile(file_name.c_str());
   // doing this 100 times takes 6s - might be quicker if I don't keep adding Hs
   // to the same residues :-)
   int imol = 0; // dummy
   coot::reduce r(mol, imol);
   bool go_nuclear = false;
   r.add_hydrogen_atoms(go_nuclear);
   mol->WritePDBASCII("reduced.pdb");
   delete mol;
   return 1;
}

int test_glyco_link_by_geometry() {

   std::vector<std::string> file_names;
   file_names.push_back("beta1-6-example.pdb");

   for (std::size_t i=0; i<file_names.size(); i++) {
      const std::string &file_name = file_names[i];
      if (coot::file_exists(file_name)) {
         mmdb::Manager *mol = new mmdb::Manager;
         mol->ReadCoorFile(file_name.c_str());

         // find linked carbohydrates and test axial vsl equatorial
      }
   }

   return 1;

}

#include <clipper/ccp4/ccp4_map_io.h>
#include "coot-map-utils.hh"

int test_soi(int argc, char **argv) {

   std::string file_name = "test.map";
   if (argc == 2)
      file_name = argv[1];

   try {
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(file_name);
      file.import_xmap(xmap);

      coot::util::soi_variance sv(xmap);

      /* we can't do this yet
      clipper::Xmap<float> outmap = sv.proc(0.66);
      clipper::CCP4MAPfile outmapfile;
      outmapfile.open_write("soi.map");
      outmapfile.export_xmap(outmap);
      outmapfile.close_write();
      */
      sv.proc(0.66);

   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << file_name << std::endl;
   }

   return 1;
}

#include <iomanip>
#include "coot-map-heavy.hh"

#define HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#include "utils/split-indices.hh"
#include "utils/ctpl.h"

void density_for_atoms_multithread(int thread_index,
                                   const atom_selection_container_t &asc,
                                   const clipper::RTop<> &rtop_og,
                                   const std::pair<unsigned int, unsigned int> &atom_index_range,
                                   const clipper::NXmap<float> &nxmap,
                                   float *dv,
                                   std::atomic<unsigned int> &done_count_for_threads) {

   for (unsigned int i=atom_index_range.first; i<atom_index_range.second; i++) {
      mmdb::Atom *at = asc.atom_selection[i];
      clipper::Coord_orth pt = coot::co(at);
      clipper::Coord_map cm_try_2(rtop_og * pt);
      float dn = coot::util::density_at_point_by_cubic_interp(nxmap, cm_try_2);
      *dv += dn;
   }

   done_count_for_threads++;
}

#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

int test_nxmap_simple(int argc, char **argv) {

   int status = 0;

   if (argc <= 2) {
      std::cout << "Usage: " << argv[0] << " map_file_name pdb_file_name" << std::endl;
   } else {
      // Happy path
      std::string map_file_name = argv[1];
      std::string pdb_file_name = argv[2];
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);
      clipper::NXmap<float> nxmap = coot::util::make_nxmap(xmap, asc);

      clipper::CCP4MAPfile mapout;
      mapout.open_write("nxmap.map");
      mapout.set_cell(xmap.cell());
      mapout.export_nxmap(nxmap);
      mapout.close_write();
   }
   return 1;
}

int test_nxmap(int argc, char **argv) {

   int status = 0;

   if (argc <= 2) {
      std::cout << "Usage: " << argv[0] << " map_file_name pdb_file_name" << std::endl;
   } else {
      // Happy path
      std::string map_file_name = argv[1];
      std::string pdb_file_name = argv[2];
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);

      clipper::NXmap<float> nxmap = coot::util::make_nxmap(xmap, asc);
      int n_atoms_max = asc.n_selected_atoms;

      std::cout << "debug: xmap  grid " <<  xmap.grid_sampling().format() << std::endl;
      std::cout << "debug: nxmap grid " << nxmap.grid().format() << std::endl;
      clipper::RTop<> rtop_og = nxmap.operator_orth_grid();
      clipper::RTop<> rtop_go = nxmap.operator_grid_orth();
      std::cout << "operators\n" << rtop_og.format() << std::endl;
      std::cout << rtop_go.format() << std::endl;

      float min_x =  999;
      float max_x = -999;
      float min_y =  999;
      float max_y = -999;
      float min_z =  999;
      float max_z = -999;

      clipper::NXmap_base::Map_reference_index ix;
      for (ix = nxmap.first(); !ix.last(); ix.next() )  { // iterator index.
         clipper::Coord_grid cg = ix.coord();
         clipper::Coord_map  cm = cg.coord_map();
         clipper::Coord_orth pt = nxmap.coord_orth(cm);
         // std::cout << "    " << pt.format() << std::endl;
         if (pt.x() < min_x) min_x = pt.x();
         if (pt.x() > max_x) max_x = pt.x();
         if (pt.y() < min_y) min_y = pt.y();
         if (pt.y() > max_y) max_y = pt.y();
         if (pt.z() < min_z) min_z = pt.z();
         if (pt.z() > max_z) max_z = pt.z();
      }

      std::cout << "nx grid extents: x " << min_x << " " << max_x << std::endl;
      std::cout << "nx grid extents: y " << min_y << " " << max_y << std::endl;
      std::cout << "nx grid extents: z " << min_z << " " << max_z << std::endl;

      auto tp_0 = std::chrono::high_resolution_clock::now();
      for(int i=0; i<n_atoms_max; i++) {
         mmdb::Atom *at = asc.atom_selection[i];
         clipper::Coord_orth pt = coot::co(at);
         float dx = coot::util::density_at_point(xmap, pt);
      }
      auto tp_1 = std::chrono::high_resolution_clock::now();

      for(int i=0; i<n_atoms_max; i++) {
         // std::cout << "atom i " << i << std::endl;
         mmdb::Atom *at = asc.atom_selection[i];
         clipper::Coord_orth pt = coot::co(at);

         clipper::Coord_map cm_try_2(rtop_og * pt);
         float dn = coot::util::density_at_point_by_cubic_interp(nxmap, cm_try_2);
      }
      auto tp_2 = std::chrono::high_resolution_clock::now();

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
      // let's do that with the thread pool

      unsigned int n_threads = 3; // 3 is faster than 4.
      ctpl::thread_pool thread_pool(n_threads);
      std::atomic<unsigned int> done_count_for_threads(0);
      std::vector<float> dv(n_threads, 0.0);
      std::vector<std::pair<unsigned int, unsigned int> > ranges =
         coot::atom_index_ranges(n_atoms_max, n_threads);
      auto tp_3 = std::chrono::high_resolution_clock::now();
      for (std::size_t i=0; i<ranges.size(); i++) {
         thread_pool.push(density_for_atoms_multithread,
                          std::cref(asc),
                          std::cref(rtop_og),
                          std::cref(ranges[i]),
                          std::cref(nxmap),
                          &dv[i],
                          std::ref(done_count_for_threads));
      }
      auto tp_4 = std::chrono::high_resolution_clock::now();
      while (done_count_for_threads < ranges.size()) {
         std::this_thread::sleep_for(std::chrono::microseconds(1));
      }
      auto tp_5 = std::chrono::high_resolution_clock::now();
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

      auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
      std::cout << "Timings:: d10  " << std::setw(4) << d10 << " microseconds" << std::endl;
      std::cout << "Timings:: d21  " << std::setw(4) << d21 << " microseconds" << std::endl;
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
      auto d43 = std::chrono::duration_cast<std::chrono::microseconds>(tp_4 - tp_3).count();
      auto d54 = std::chrono::duration_cast<std::chrono::microseconds>(tp_5 - tp_4).count();
      std::cout << "Timings:: d43  " << std::setw(4) << d43 << " microseconds" << std::endl;
      std::cout << "Timings:: d54  " << std::setw(4) << d54 << " microseconds" << std::endl;
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   }

   return status;

}

void
nxmap_fft_to(const clipper::Xmap<float> &xmap,
             clipper::HKL_data<clipper::data32::F_phi> &fphidata) {

   float mg = coot::util::max_gridding(xmap);
   clipper::Resolution reso(2.0 * mg);
   clipper::HKL_info myhkl(xmap.spacegroup(), xmap.cell(), reso, true);
   // clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(myhkl);
   fphidata.init(myhkl, xmap.cell());

   // clipper::Xmap<float> xmap = somehow nxmap-to-xmap

   xmap.fft_to(fphidata);
}

int
test_nxmap_edcalc(int argc, char **argv) {

   if (argc <= 2) {
      std::cout << "Usage: " << argv[0] << " map_file_name pdb_file_name" << std::endl;
   } else {
      std::string map_file_name = argv[1];
      std::string pdb_file_name = argv[2];
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);

      // now a few residues
      coot::residue_spec_t spec_1("A", 10, "");
      coot::residue_spec_t spec_2("A", 11, "");
      coot::residue_spec_t spec_3("A", 12, "");
      coot::residue_spec_t spec_4("A", 13, "");
      coot::residue_spec_t spec_5("A", 14, "");

      int few_residues_selection_handle = asc.mol->NewSelection();

      spec_1.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_NEW);
      spec_2.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      spec_3.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      spec_4.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      spec_5.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);

      std::cout << "Making nxmap... " << std::endl;
      clipper::NXmap<float> nxmap_ref =
         coot::util::make_nxmap(xmap, asc.mol, few_residues_selection_handle);
      std::cout << "ED calc..." << std::endl;
      clipper::NXmap<float> nxmap_edcalc =
         coot::util::make_edcalc_map(nxmap_ref,  // for metrics
                                     asc.mol, few_residues_selection_handle);

      asc.mol->DeleteSelection(few_residues_selection_handle);

      clipper::CCP4MAPfile mapout;
      mapout.open_write("edcalc.map");
      mapout.set_cell(xmap.cell());
      mapout.export_nxmap(nxmap_edcalc);
      mapout.close_write();

      clipper::CCP4MAPfile mapout_ref;
      mapout.open_write("edcalc_ref.map");
      mapout.set_cell(xmap.cell());
      mapout.export_nxmap(nxmap_ref);
      mapout.close_write();

      clipper::HKL_data<clipper::data32::F_phi> fphi_calc;
      clipper::HKL_data<clipper::data32::F_phi> fphi_ref;

      std::cout << "fft-calc" << std::endl;
      nxmap_fft_to(xmap, fphi_calc);

      // how do I make an xmap from an nxmap for ref? Too hard for now.

      std::cout << "fft-ref" << std::endl;
      nxmap_fft_to(xmap, fphi_ref);

      std::cout << "done ffts" << std::endl;
      int n_f_calc = fphi_calc.data_size();
      int n_f_ref  = fphi_ref.data_size();

      // std::cout << "info:: n_f_calc " << fphi_calc.debug() << std::endl;
      // std::cout << "info:: n_f_ref  " << n_f_ref  << std::endl;

      // fphi_calc.debug();
      // fphi_ref.debug();

      int count = 0;
      clipper::HKL_info::HKL_reference_index hri;

      hri = fphi_calc.first();
      if (hri.last()) {
         std::cout << "booo... first is last " << std::endl;
      }
      for (hri = fphi_calc.first(); !hri.last(); hri.next()) {
         std::cout << "   " << fphi_calc[hri].f() << " " << fphi_ref[hri].f()
                   << std::endl;

         count++;
         if (count == 10)
            break;
      }

      // bricks

      std::cout << "---------- testing bricks " << std::endl;
      int t1 = coot::get_brick_id_inner(0,0,0, 2,2,2);
      int t2 = coot::get_brick_id_inner(1,0,0, 6,4,8);
      int t3 = coot::get_brick_id_inner(0,1,0, 6,4,8);
      int t4 = coot::get_brick_id_inner(1,1,0, 6,4,8);
      int t5 = coot::get_brick_id_inner(3,2,2, 6,4,8);

      std::cout << "test: t1, t2, t3, t4 t5 "
                << t1 << " " << t2 << " " << t3 << " " << t4 << " " << t5
                << std::endl;


      float atom_max_radius = 3.0;
      std::vector<std::vector<int> > bricks =
         coot::molecule_to_bricks(asc.mol, asc.SelectionHandle, atom_max_radius);

      std::cout << "found " << bricks.size() << " bricks for "
                << asc.n_selected_atoms << " atoms " << std::endl;

      // this looks fine
      if (false) {
         for (std::size_t ii=0; ii<bricks.size(); ii++) {
            std::cout << "--- brick " << ii << std::endl;
            for (std::size_t jj=0; jj<bricks[ii].size(); jj++) {
               std::cout << "   " << bricks[ii][jj] << " "
                         << coot::atom_spec_t(asc.atom_selection[bricks[ii][jj]]) << "\n";
            }
         }
      }
   }
   return 1;
}

int
test_xmap_edcalc(int argc, char **argv) {

   if (argc <= 2) {
      std::cout << "Usage: " << argv[0] << " map_file_name pdb_file_name" << std::endl;
   } else {
      std::string map_file_name = argv[1];
      std::string pdb_file_name = argv[2];
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);

      // now a few residues
      coot::residue_spec_t spec_1("A", 10, "");
      coot::residue_spec_t spec_2("A", 11, "");
      coot::residue_spec_t spec_3("A", 12, "");
      coot::residue_spec_t spec_4("A", 13, "");
      coot::residue_spec_t spec_5("A", 14, "");

      int few_residues_selection_handle = asc.mol->NewSelection();

      spec_1.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_NEW);
      spec_2.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      spec_3.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      spec_4.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      spec_5.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);

      clipper::Xmap<float> calc_atom_map(mmdb::Manager *mol,
                                         int atom_selection_handle,
                                         const clipper::Cell &cell,
                                         const clipper::Spacegroup &space_group,
                                         const clipper::Grid_sampling &sampling);

      asc.mol->DeleteSelection(few_residues_selection_handle);

   }

   return 0;
}


int
test_string_split() {

   int status = 1;

   std::cout << "test_string_split()" << std::endl;

   std::string s = "All Loop candidates ";
   std::vector<std::string> v = coot::util::split_string(s, " ");
   std::cout << "0: \"" << v[0] << "\"" << std::endl;
   std::cout << "1: \"" << v[1] << "\"" << std::endl;
   std::cout << "2: \"" << v[2] << "\"" << std::endl;

   return status;

}

#include "bonded-atoms.hh"

int test_bonded_atoms(int argc, char **argv) {

   int status = 1;
   if (argc > 2) {
      std::string pdb_file_name = argv[2];
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);

      auto tp_0 = std::chrono::high_resolution_clock::now();
      std::vector<std::vector<unsigned int> > b = coot::make_bonds(asc.mol,
                                                                   asc.n_selected_atoms,
                                                                   asc.UDDAtomIndexHandle);

      auto tp_1 = std::chrono::high_resolution_clock::now();

      // Timings: 2844 us make_bonds()
      //           193 us find_1_4_connections()
      //
      std::vector<std::vector<unsigned int> > connections_1_4 = coot::find_1_4_connections(b);

      for (std::size_t i=0; i<connections_1_4.size(); i++) {
         const std::vector<unsigned int> &v1 = connections_1_4[i];
         mmdb::Atom *at_i = asc.atom_selection[i];
         for (std::size_t j=0; j<v1.size(); j++) {
            mmdb::Atom *at_j = asc.atom_selection[v1[j]];
            if (false)
               std::cout << " 1-4: " << coot::atom_spec_t(at_i) << " "
                         << coot::atom_spec_t(at_j) << std::endl;
         }
      }
      auto tp_2 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
      std::cout << "Timings:: d10  " << std::setw(4) << d10 << " microseconds" << std::endl;
      std::cout << "Timings:: d21  " << std::setw(4) << d21 << " microseconds" << std::endl;
   }


   return status;
}

#include "helix-like.hh"

int test_helix_like(int argc, char **argv) {

    if (argc == 2) {
       std::string file_name = argv[1];
       atom_selection_container_t asc = get_atom_selection(file_name, false, true, true);
       if (asc.read_success) {
          std::vector<std::string> ch_ids;
          int imod = 1;
          mmdb::Model *model_p = asc.mol->GetModel(imod);
          if (model_p) {
             int n_chains = model_p->GetNumberOfChains();
             for (int ichain=0; ichain<n_chains; ichain++) {
                mmdb::Chain *chain_p = model_p->GetChain(ichain);
                ch_ids.push_back(chain_p->GetChainID());
             }
          }

          for (unsigned int ich=0; ich<ch_ids.size(); ich++) {
             int residue_selection_handle = asc.mol->NewSelection();
             asc.mol->Select (residue_selection_handle, mmdb::STYPE_RESIDUE, 0,
                    ch_ids[ich].c_str(),
                    mmdb::ANY_RES, "*",  // starting res
                    mmdb::ANY_RES, "*",  // ending res
                    "*",  // residue name
                    "*",  // Residue must contain this atom name?
                    "*",  // Residue must contain this Element?
                    "*",  // altLocs
                    mmdb::SKEY_NEW // selection key
                    );
             std::cout << "test like_a_helix()" << std::endl;
             coot::like_a_helix(asc.mol, residue_selection_handle);
             asc.mol->DeleteSelection(residue_selection_handle);
          }
       }
    }
    return 0;
}

#include "cablam-markup.hh"

int
test_cablam(int argc, char **argv) {

   std::string cablam_log_file_name = "../src/cablam.log";
   std::string file_name = argv[1];
   atom_selection_container_t asc = get_atom_selection(file_name, false, true, true);
   if (asc.read_success) {
      // parse this log file and call the above function for each cablam outlier residue
      std::vector<coot::cablam_markup_t> v =
      coot::make_cablam_markups(asc.mol, cablam_log_file_name);
   }
   return 0;
}


#include "merge-atom-selections.hh"

int
test_merge_fragments(int argc, char **argv) {

   if (argc == 2) {
      std::string file_name = argv[1];
      atom_selection_container_t asc = get_atom_selection(file_name, false, true, true);
      if (asc.read_success) {
         std::vector<std::string> ch_ids;
         int imod = 1;
         mmdb::Model *model_p = asc.mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               ch_ids.push_back(chain_p->GetChainID());
            }
         }
         coot::merge_atom_selections(asc.mol);
         asc.mol->WritePDBASCII("merged.pdb");
      }
   }
   return 1;
}

#include <fstream>

int
test_make_a_difference_map(int argc, char **argv) {

   // This function is for all map and all model - not a model selection
   // that is quite a bit more tricky.

   // B-factor scaling:
   // <f^2> = k <Fo^2> * exp(-2B * r)
   //    where r = 1/d^2 where d is in A.
   // -2B * r + ln(k) = log(<f^2>/<Fo^2>)

   int status = 1;

   if (argc <= 2) {
      std::cout << "Usage: " << argv[0] << " map_file_name pdb_file_name" << std::endl;
   } else {

      // ----- 1 ------------------ Parse

      std::string map_file_name = argv[1];
      std::string pdb_file_name = argv[2];
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);

      // Matching structure factors on an atom selection is a different function.
      //
      // now a few residues
      // coot::residue_spec_t spec_1("A", 10, "");
      // coot::residue_spec_t spec_2("A", 11, "");
      // int few_residues_selection_handle = asc.mol->NewSelection();
      // spec_1.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_NEW);
      // spec_2.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      // spec_3.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      // spec_4.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);
      // spec_5.select_atoms(asc.mol, few_residues_selection_handle, mmdb::SKEY_OR);

      // ----- 2 ------------------ Make Calc Map

      auto tp_0 = std::chrono::high_resolution_clock::now();
      clipper::Xmap<float> xmap_calc = coot::util::calc_atom_map(asc.mol,
                                                                 asc.SelectionHandle,
                                                                 xmap.cell(),
                                                                 xmap.spacegroup(),
                                                                 xmap.grid_sampling());

      auto tp_1 = std::chrono::high_resolution_clock::now();
      if (true) {
         clipper::CCP4MAPfile outmapfile;
         outmapfile.open_write("test-calc.map");
         outmapfile.export_xmap(xmap_calc);
         outmapfile.close_write();
      }


      if (false) {
         xmap_calc = xmap;
         clipper::Xmap_base::Map_reference_index ix;
         for (ix = xmap_calc.first(); !ix.last(); ix.next() ) {
            xmap_calc[ix] = xmap[ix] * 2.0 + 0.11;
         }
      }

      // ----- 3 ------------------ Get Amp vs Reso data for B-factor from map

      std::vector<coot::amplitude_vs_resolution_point> pts_ref  = coot::util::amplitude_vs_resolution(xmap,      15);
      std::vector<coot::amplitude_vs_resolution_point> pts_calc = coot::util::amplitude_vs_resolution(xmap_calc, 15);

      if (true) { // debugging B-factor estimation
         for (unsigned int ii=0; ii< pts_ref.size(); ii++) {
            std::cout << "pts_ref " << pts_ref[ii].average << " " << pts_ref[ii].resolution_recip << std::endl;
         }
         for (unsigned int ii=0; ii< pts_ref.size(); ii++) {
            std::cout << "pts_calc " << pts_calc[ii].average << " " << pts_calc[ii].resolution_recip << std::endl;
         }
      }


      // ----- 4 ------------------ Get B-factors from maps

      // asc.mol->DeleteSelection(few_residues_selection_handle);

      // clipper::HKL_data<clipper::data32::F_phi> fphi_calc;
      // clipper::HKL_data<clipper::data32::F_phi> fphi_ref;

      clipper::Resolution reso(1.06);
      float reso_min_for_scaling = 0.05; // 4.47A or so is the minimum resolution
                                         // for B-factor scaling
                                         // If there are only data below that
                                         // then use simple (1 parameter) scaling.
      double reso_max = 0.29; // inv res sq., calculate this from an analysis of the map f values
      reso_max = 0.89; // 1pwg

      std::pair<bool, float> reso_low_invresolsq(true, 0.05);
      std::pair<bool, float> reso_high_invresolsq(true, reso_max);

      float b1 = coot::util::b_factor(pts_ref,  reso_low_invresolsq, reso_high_invresolsq);
      float b2 = coot::util::b_factor(pts_calc, reso_low_invresolsq, reso_high_invresolsq);

      std::cout << "b1 " << b1 << " b2 " << b2 << std::endl;

      // ----- 5 ------------------ Calculate Structure Factors

      clipper::HKL_info myhkl(xmap.spacegroup(), xmap.cell(), reso, true);

      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphi_calc(myhkl);
      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphi_ref(myhkl);

      std::cout << "fft-calc" << std::endl;
      auto tp_2 = std::chrono::high_resolution_clock::now();
      xmap_calc.fft_to(fphi_calc);

      auto tp_3 = std::chrono::high_resolution_clock::now();
      std::cout << "fft-ref" << std::endl;
      xmap.fft_to(fphi_ref);
      auto tp_4 = std::chrono::high_resolution_clock::now();

      int count = 0;
      clipper::HKL_info::HKL_reference_index hri;
      hri = fphi_calc.first();
      if (hri.last()) {
         std::cout << "booo... first is last " << std::endl;
      }
      if (true) { // just check that the SFS contain data
         for (hri = fphi_ref.first(); !hri.last(); hri.next()) {
            std::cout << "   " << hri.hkl().format() << " " << fphi_ref[hri].f() << " " << fphi_ref[hri].phi()
                      << std::endl;
            count++;
            if (count == 10)
               break;
         }
         count = 0;
         for (hri = fphi_calc.first(); !hri.last(); hri.next()) {
            std::cout << "   " << hri.hkl().format() << " " << fphi_calc[hri].f() << " " << fphi_calc[hri].phi()
                      << std::endl;
            count++;
            if (count == 10)
               break;
         }
      }

      // ----- 6 ------------------ Scale and get stats for R-factor

      const unsigned int n_bins = 20;

      // if these are data from an mtz say, we'd need to check for isnan, but it's from a map
      // so there is no nan data.

      std::vector<double> sum_fo(n_bins, 0.0);
      std::vector<double> sum_fc(n_bins, 0.0);
      std::vector<unsigned int> counts(n_bins, 0);

      for (hri = fphi_calc.first(); !hri.last(); hri.next()) {
         if (hri.hkl() != clipper::HKL(0,0,0)) {
            float irs = hri.invresolsq();
            if (irs < reso_max) {
               int bin_idx = static_cast<int>(static_cast<float>(n_bins) * irs/reso_max);
               if (bin_idx == n_bins) bin_idx--;
               if (irs >= reso_min_for_scaling) {
                  sum_fc[bin_idx] += fabs(fphi_calc[hri].f());
                  counts[bin_idx]++;
               }
            }
         }
      }

      for (hri = fphi_ref.first(); !hri.last(); hri.next()) {
         if (hri.hkl() != clipper::HKL(0,0,0)) {
            float irs = hri.invresolsq();
            if (irs < reso_max) {
               int bin_no = static_cast<int>(n_bins * irs/reso_max);
               if (bin_no == n_bins) bin_no--;
               if (irs >= reso_min_for_scaling)
                  sum_fo[bin_no] += fabs(fphi_ref[hri].f());
            }
         }
      }

      for (std::size_t i=0; i<n_bins; i++) {
         float bin_reso = reso_max * (static_cast<float>(i)/static_cast<float>(n_bins) + 0.5);
         if (counts[i] > 0) {
            std::cout << i << " " << bin_reso << " "
                      << sum_fo[i]/static_cast<float>(counts[i]) << " "
                      << sum_fc[i]/static_cast<float>(counts[i]) << " "
                      << sum_fo[i]/sum_fc[i] << "\n";
         }
      }
   }
   return status;
}

int test_map_molecule_centre(int argc, char **argv) {

   if (argc > 1) {
      std::string fn = argv[1];
      try {
         clipper::CCP4MAPfile file;
         clipper::Xmap<float> xmap;
         file.open_read(fn);
         file.import_xmap(xmap);
         coot::util::map_molecule_centre(xmap);
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "Failed to open " << fn << std::endl;
      }
   }
   return 0;
}

int
test_CO_orientations(int argc, char **argv) {
   int status = 0;
   coot::util::multi_parse_prosmart_log_and_gen_CO_plot();
   return status;
}

int test_fsc(int argc, char **argv) {
   int status = 0;
   if (argc > 2) {
      std::string fn1 = argv[1];
      std::string fn2 = argv[2];
      try {
         std::cout << "idx resolution count ff-sum f1-sum f2-sum fsc" << std::endl;
         clipper::CCP4MAPfile file;
         clipper::Xmap<float> xmap_1;
         clipper::Xmap<float> xmap_2;
         std::cout << "# reading map 1" << std::endl;
         file.open_read(fn1);
         file.import_xmap(xmap_1);
         file.close_read();
         std::cout << "# reading map 2" << std::endl;
         file.open_read(fn2);
         file.import_xmap(xmap_2);
         file.close_read();
         coot::util::fsc(xmap_1, xmap_2);
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "Failed to open " << fn1 << std::endl;
      }
   }
   return status;
}

int test_flip(int argc, char **argv) {
   int status = 0;
   if (argc > 1) {
           std::string fn = argv[1];
           try {
              clipper::CCP4MAPfile file;
              clipper::Xmap<float> xmap;
              std::cout << "# reading map" << std::endl;
              file.open_read(fn);
              file.import_xmap(xmap);
              file.close_read();
              coot::util::flip_hand(&xmap);
              clipper::CCP4MAPfile outmapfile;
              outmapfile.open_write("flipped-hand.map");
              outmapfile.export_xmap(xmap);
              outmapfile.close_write();
           }
           catch (const clipper::Message_base &exc) {
                   std::cout << "Failed to flip" << fn << std::endl;
           }
   }
   return status;
}

int test_interface_residues(int argc, char **argv) {

   if (argc > 1) {
      std::string pdb_file_name = argv[1]; // 6lzg
      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);
      if (asc.read_success) {
         float min_dist = 3.8;
         std::pair<std::set<mmdb::Residue *>, std::set<mmdb::Residue *> > ir = coot::interface_residues(asc.mol, "A", "B", min_dist);

         std::set<mmdb::Residue *>::const_iterator it;
         for (it=ir.first.begin(); it!=ir.first.end(); ++it) {
            mmdb::Residue *r = *it;
            std::cout << "   B " << coot::residue_spec_t(r) << std::endl;
         }
         for (it=ir.second.begin(); it!=ir.second.end(); ++it) {
            mmdb::Residue *r = *it;
            std::cout << "   B " << coot::residue_spec_t(r) << std::endl;
         }
      }
   }
   return 0;

}

#include "fib-sphere.hh"
void
test_fibonacci() {

   unsigned int n_samples = 450;
   std::vector<clipper::Coord_orth> pts = coot::fibonacci_sphere(n_samples);
   for (unsigned int i=0; i<n_samples; i++)
      std::cout << "spherical " << pts[i].x() << " " << pts[i].y() << " " << pts[i].z() << "\n";

}

#include "polar-atoms.hh"

void
test_polar_atom_analysis(int argc, char **argv) {

   if (argc > 1) {
      std::string pdb_file_name = argv[1]; // 6lzg
      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);
      if (asc.read_success) {
         coot::buried_unsatisfied_polar_atoms(asc.mol);
      }
   }
}

#include "fragment-container.hh"

void
test_fragment_maker(int argc, char **argv) {
   if (argc > 1) {
      std::string pdb_file_name = argv[1]; // 6lzg
      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);
      if (asc.read_success) {
         coot::fragment_container_t fragments = coot::make_fragments(asc.mol);
         fragments.print_fragments();
      }
   }
}

#include <clipper/core/atomsf.h>

void
test_correlation_of_residue_runs(int argc, char **argv) {

   bool is_cryo_em = true;
   bool exclude_NOC = true;

   if (argc > 5) {
      std::string pdb_file_name = argv[1];
      std::string chain_id      = argv[2];
      std::string map_file_name = argv[3];
      std::string atom_mask_radius_str = argv[4];
      std::string NOC_mask_radius_str  = argv[5];
      float atom_mask_radius = 3.0;
      float NOC_mask_radius = 3.0;

      try {
         atom_mask_radius = coot::util::string_to_float(atom_mask_radius_str);
         NOC_mask_radius  = coot::util::string_to_float(NOC_mask_radius_str);
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR::" << rte.what() << std::endl;
         return;
      }

      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);
      if (asc.read_success) {

         std::cout << "pdb read success " << pdb_file_name << std::endl;

         clipper::CCP4MAPfile file;
         clipper::Xmap<float> xmap;
         std::cout << "# reading map" << std::endl;
         file.open_read(map_file_name);
         file.import_xmap(xmap);
         file.close_read();

#if 0 // I am compiling this with old coot for some reason.
         if (is_cryo_em)
            clipper::ScatteringFactors::selectScattteringFactorsType(clipper::SF_ELECTRON);

         // if exclude_NOC is set then the second part of the pair (for side chains) is also filled.

         unsigned int n_residue_per_residue_range = 11;
         std::pair<std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>,
                   std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t> > residue_stats =
            coot::util::map_to_model_correlation_stats_per_residue_run(asc.mol, chain_id, xmap,
                                                                       n_residue_per_residue_range, exclude_NOC,
                                                                       atom_mask_radius, NOC_mask_radius);
         std::cout << "INFO:: We got " << residue_stats.first.size()  << " residue all-atom correlations"   << std::endl;
         std::cout << "INFO:: We got " << residue_stats.second.size() << " residue side-chain correlations" << std::endl;

         std::map<coot::residue_spec_t, coot::util::density_correlation_stats_info_t>::const_iterator it;
         for (it=residue_stats.first.begin(); it!=residue_stats.first.end(); ++it) {
            const coot::residue_spec_t &rs(it->first);
            const coot::util::density_correlation_stats_info_t &stats(it->second);
            std::cout << "   all-atom-stats " << rs << " " << stats.correlation() << " from " << stats.n << " points ";
            std::cout << std::endl;
         }
         for (it=residue_stats.second.begin(); it!=residue_stats.second.end(); ++it) {
            const coot::residue_spec_t &rs(it->first);
            const coot::util::density_correlation_stats_info_t &stats(it->second);
            std::cout << "   side-chain-stats " << rs << " " << stats.correlation() << " from " << stats.n << " points ";
            std::cout << std::endl;
         }
#endif
      } else {
         std::cout << "Failed to read " << pdb_file_name << std::endl;
      }
   } else {
      std::cout << "Usage: " << argv[0] << "pdb_file_name chain_id map_file_name atom-mask-radius NOC-mask-radius"
                << std::endl;
   }
}

#include "merge-C-and-N-terminii.hh"

void
test_merge_C_and_N_terminii(int argc, char **argv) {

   if (argc > 1) {
      std::string pdb_file_name = argv[1];

      std::cout << "Getting atoms... " << std::endl;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);
      if (asc.read_success) {
         std::cout << "pdb read success " << pdb_file_name << std::endl;
         coot::merge_C_and_N_terminii_0_gap(asc.mol);
         asc.mol->WritePDBASCII("C-N-merged.pdb");
      }
   }

}

void
test_compare_structure_factors(int argc, char **argv) {

   if (argc> 2) {
      std::string map_file_name_1 = argv[1];
      std::string map_file_name_2 = argv[2];
      clipper::CCP4MAPfile file_1;
      clipper::CCP4MAPfile file_2;
      clipper::Xmap<float> xmap_1;
      clipper::Xmap<float> xmap_2;
      std::cout << "# reading map " << map_file_name_1 << std::endl;
      file_1.open_read(map_file_name_1);
      std::cout << "# reading map " << map_file_name_2 << std::endl;
      file_2.open_read(map_file_name_2);
      std::cout << "# importing map " << map_file_name_1 << std::endl;
      file_1.import_xmap(xmap_1);
      std::cout << "# importing map " << map_file_name_2 << std::endl;
      file_2.import_xmap(xmap_2);

      coot::util::compare_structure_factors(xmap_1, xmap_2);
   }

}

void
test_ncs_chain_match(int argc, char **argv) {

   if (argc > 1) {
      std::string pdb_file_name(argv[1]);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, true, false);
      int imod = 1;
      std::vector<std::vector<mmdb::Chain *> > v = coot::ncs_related_chains(asc.mol, imod);
      std::cout << "------------------------------------ ncs --------------------------" << std::endl;
      unsigned int n_chains = 0;
      for (unsigned int iv1=0; iv1<v.size(); ++iv1) {
         // std::cout << "New Set .... " << std::endl;
         std::vector<mmdb::Chain *> &vv = v[iv1];
         for (unsigned int iv2=0; iv2<vv.size(); ++iv2) {
            mmdb::Chain *c = vv[iv2];
            std::string chain_id(c->GetChainID());
            std::cout << " " << chain_id;
            n_chains++;
         }
         std::cout << std::endl;
      }
      std::cout << "Found " << n_chains << " total chains" << std::endl;
   }
}

void
test_dictionary_conformers(int argc, char **argv) {

   std::string monomer_type = "TYR";
   int imol = 0;
   coot::protein_geometry geom;
   geom.init_standard();
   geom.set_verbose(false);
   std::pair<bool, coot::dictionary_residue_restraints_t> r = geom.get_monomer_restraints(monomer_type, imol);
   if (r.first) {
      bool delete_clash_confs = true;
      std::cout << "------------------- calling get_dictionary_conformers() " << std::endl;
      std::vector<mmdb::Residue *> confs = coot::util::get_dictionary_conformers(r.second, delete_clash_confs);
      for (unsigned int i=0; i<confs.size(); i++) {
         mmdb::Residue *res = confs[i];
         mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(res);
         std::string fn = "conf-" + std::to_string(i) + ".pdb";
         mol->WritePDBASCII(fn.c_str());
         delete mol;
      }
      for (unsigned int i=0; i<confs.size(); i++)
         delete confs[i];
   }
}


void test_partition_map_by_chain(int argc, char **argv) {

   if (argc > 2) {
      std::string map_file_name_1 = argv[1];
      clipper::CCP4MAPfile file_1;
      clipper::Xmap<float> xmap_1;
      std::cout << "# reading map " << map_file_name_1 << std::endl;
      file_1.open_read(map_file_name_1);
      file_1.import_xmap(xmap_1);


      clipper::Cell cell = xmap_1.cell();
      clipper::Spacegroup sg = xmap_1.spacegroup();
      clipper::Grid_sampling gs = xmap_1.grid_sampling();

      std::cout << "Cell:" << cell.format() << std::endl;
      std::cout << "Spacegroup:" << sg.symbol_hm() << std::endl;
      std::cout << "Grid Sampling:" << gs.format() << std::endl;

      std::string pdb_file_name = argv[2];
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, false);

      if (asc.read_success) {
         std::string status_string;
         std::vector<std::pair<std::string, clipper::Xmap<float> > > maps =
            coot::util::partition_map_by_chain(xmap_1, asc.mol, &status_string);
         for (size_t i=0; i < maps.size(); i++) {
            std::string map_file_name = "partitioned-" + maps[i].first + ".map";
            const auto &xmap = maps[i].second;
            clipper::CCP4MAPfile outmapfile;
            outmapfile.open_write(map_file_name);
            outmapfile.export_xmap(xmap);
            outmapfile.close_write();
         }
      }
   }
}


void test_make_mask_map(int argc, char **argv) {

   if (argc > 2) {
      std::string map_file_name_1 = argv[1];
      clipper::CCP4MAPfile file_1;
      clipper::Xmap<float> xmap;
      std::cout << "# reading map " << map_file_name_1 << std::endl;
      file_1.open_read(map_file_name_1);
      file_1.import_xmap(xmap);

      clipper::Cell cell = xmap.cell();
      clipper::Spacegroup sg = xmap.spacegroup();
      clipper::Grid_sampling gs = xmap.grid_sampling();

      std::cout << "Cell: " << cell.format() << std::endl;
      std::cout << "Spacegroup: " << sg.symbol_hm() << std::endl;
      std::cout << "Grid Sampling: " << gs.format() << std::endl;

      std::string pdb_file_name = argv[2];
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, false);

      if (asc.read_success) {
         mmdb::Manager *mol = asc.mol;
         int selection_handle = mol->NewSelection();
         std::string selection_string = "/";
         asc.mol->Select(selection_handle, mmdb::STYPE_ATOM, selection_string.c_str(), mmdb::SKEY_NEW);
         float radius = 4.5f;
         float smooth = 1.0f;
         clipper::Xmap<float> mask_xmap =
            coot::util::make_map_mask(sg, cell, gs, asc.mol, selection_handle, radius, smooth);
         clipper::CCP4MAPfile outmapfile;
         std::string map_file_name = "A-chain-mask.map";
         outmapfile.open_write(map_file_name);
         outmapfile.export_xmap(mask_xmap);
         outmapfile.close_write();
         mol->DeleteSelection(selection_handle);
      }
   }
}

#include "cremer-pople.hh"

void
test_cremer_pople(int argc, char **argv) {

   auto is_member = [] (const std::string &rn, const std::vector<std::string> &res_names) {
      return std::find(res_names.begin(), res_names.end(), rn) != res_names.end();
   };

   if (argc > 1) {
      std::string pdb_file_name(argv[1]);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false);
      if (asc.read_success) {
         mmdb::Manager *mol = asc.mol;
         std::vector<std::string> res_names = {"NAG", "MAN", "BMA", "GAL", "FUC"};
         for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int n_res = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<n_res; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        std::string rn = residue_p->GetResName();
                        if (is_member(rn, res_names)) {
                           coot::cremer_pople_t cpi(residue_p);
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

int test_ribose_torsions(int argc, char **argv) {

   std::string pdb_file_name = "8b0x.cif";
   if (argc > 1)
      pdb_file_name = argv[1];

   int status = 0;
   atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, true);
   std::vector<std::string> base_names = {"A", "G", "T", "U"};
   std::string alt_conf = "";
   if (asc.mol) {
      int imod = 1;
      mmdb::Model *model_p = asc.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  std::string res_name = residue_p->GetResName();
                  if (std::find(base_names.begin(), base_names.end(), res_name) != base_names.end()) {
                     mmdb::Atom *at_O4_prime = nullptr;
                     mmdb::Atom *at_C1_prime = nullptr;
                     mmdb::Atom *at_C2_prime = nullptr;
                     mmdb::Atom *at_C3_prime = nullptr;
                     mmdb::Atom *at_C4_prime = nullptr;
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (! at->isTer()) {
                           std::string atom_name = at->GetAtomName();
                           if (atom_name == " O4'") at_O4_prime = at;
                           if (atom_name == " C1'") at_C1_prime = at;
                           if (atom_name == " C2'") at_C2_prime = at;
                           if (atom_name == " C3'") at_C3_prime = at;
                           if (atom_name == " C4'") at_C4_prime = at;
                        }
                     }
                     if (at_O4_prime && at_C1_prime && at_C2_prime && at_C3_prime && at_C4_prime) {
                        coot::atom_quad aq_1(at_O4_prime, at_C1_prime, at_C2_prime, at_C3_prime);
                        coot::atom_quad aq_2(at_C1_prime, at_C2_prime, at_C3_prime, at_C4_prime);
                        coot::atom_quad aq_3(at_C2_prime, at_C3_prime, at_C4_prime, at_O4_prime);
                        coot::atom_quad aq_4(at_C3_prime, at_C4_prime, at_O4_prime, at_C1_prime);
                        coot::atom_quad aq_5(at_C4_prime, at_O4_prime, at_C1_prime, at_C2_prime);

                        double t1 = aq_1.torsion();
                        double t2 = aq_2.torsion();
                        double t3 = aq_3.torsion();
                        double t4 = aq_4.torsion();
                        double t5 = aq_5.torsion();

                        coot::pucker_analysis_info_t pi(residue_p, alt_conf);

                        std::cout << "puckered-atom" << pi.puckered_atom() << " torsions: "
                                  << std::setw(9) << t1 << " " << std::setw(9) << t2 << " "
                                  << std::setw(9) << t3 << " " << std::setw(9) << t4 << " "
                                  << std::setw(9) << t5 << std::endl;
                     }
                  }
               }
            }
         }
      }
   }
   return status;
}


int main(int argc, char **argv) {

   mmdb::InitMatType();

   if (true)
      test_ribose_torsions(argc, argv);

   if (false)
      test_cremer_pople(argc, argv);

   if (false)
      test_make_mask_map(argc, argv);

   if (false)
      test_partition_map_by_chain(argc, argv);

   if (false)
      test_dictionary_conformers(argc, argv);

   if (false)
      test_ncs_chain_match(argc, argv);

   if (false)
      test_all_atom_overlaps();

   if (false)
      test_glyco_link_by_geometry();

   if (false)
      test_string_manipulation();

   if (false)
      test_sort_chains();

   if (false)
      test_euler_angles();

   if (false)
      test_lsq_improve();

   if (false)
      test_glyco_tree();

   if (false)
      test_helix_analysis();

   if (false)
      test_qq_plot();

   if (false)
      test_least_squares_fit();

   if (false)
      test_atom_overlaps();

   if (false)
      for (unsigned int i=0; i<2; i++) {
         test_all_atom_overlaps();
      }

   if (false)
      test_reduce();

   if (false)
       test_cp();

   if (false)
      test_soi(argc, argv);

   if (false)
      test_bonded_atoms(argc, argv);

   if (false)
      test_string_split();

   if (false)
      test_xmap_edcalc(argc, argv);

   if (false)
      test_nxmap_edcalc(argc, argv);

   if (false)
      test_helix_like(argc, argv);

   if (false)
      test_merge_fragments(argc, argv);

   if (false)
      test_make_a_difference_map(argc, argv);

   if (false)
      test_nxmap_simple(argc, argv);

   if (false)
      test_CO_orientations(argc, argv);

   if (false)
      test_cablam(argc, argv);

   if (false)
      test_map_molecule_centre(argc, argv);

   if (false)
      test_fsc(argc, argv);

   if (false)
      test_flip(argc, argv);

   if (false)
      test_interface_residues(argc, argv);

   if (false)
      test_fibonacci();

   if (false)
      test_polar_atom_analysis(argc, argv);

   if (false)
      test_fragment_maker(argc, argv);

   if (false)
      test_correlation_of_residue_runs(argc, argv);

   if (false)
      test_merge_C_and_N_terminii(argc, argv);

   if (false)
      test_compare_structure_factors(argc, argv);

   return 0;
}
