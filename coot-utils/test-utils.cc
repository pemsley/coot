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

#include "clipper/core/rotation.h"

#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"
#include "lsq-improve.hh"
#include "helix-analysis.hh"


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
      catch (std::runtime_error rte) {
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
	std::vector<mmdb::Residue *> neighbs = coot::residues_near_residue(residue_p, mol, 5);
	coot::atom_overlaps_container_t overlaps(residue_p, neighbs, mol, &geom, 0.5, 0.25);
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
      bool waters = true; // include waters
      coot::atom_overlaps_container_t overlaps(mol, &geom, waters, 0.5, 0.25);
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
	    mmdb::Residue *residue_p;
	    mmdb::Atom *at;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
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
   r.add_hydrogen_atoms();
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


int main(int argv, char **argc) {

   if (1)
      test_all_atom_overlaps();

   if (0)
      test_glyco_link_by_geometry();

   if (0)
      test_string_manipulation();

   if (0)
      test_sort_chains();

   if (0)
      test_euler_angles();

   if (0)
      test_lsq_improve();

   if (0)
      test_glyco_tree();

   if (0)
      test_helix_analysis();

   if (0)
      test_qq_plot();

   if (0)
      test_least_squares_fit();
   
   if (0)
      test_atom_overlaps();

   if (0)
      for (unsigned int i=0; i<2; i++) {
	 test_all_atom_overlaps();
      }
   
   if (0)
      test_reduce();

//    if (true)
//       test_cp();
   
   return 0;
}

