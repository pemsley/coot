
#include "testing.hh"

#ifdef BUILT_IN_TESTING

#include <iostream>
#include <sstream>

#include "mmdb-extras.h"
#include "mmdb.h"
#include "primitive-chi-angles.hh"

#include "wligand.hh"

std::string greg_test(const std::string &file_name) {

   std::string d = getenv("HOME");
   d += "/data/greg-data/";
   d += file_name;
   return d;
} 

std::string stringify(double x) {
   std::ostringstream o;
   if (!(o << x))
      throw std::runtime_error("stringify(double)");
   return o.str();
}


#ifdef STRINGIFY_SIZE_T
std::string stringify(size_t x) {
   std::ostringstream o;
   if (!(o << x))
      throw std::runtime_error("stringify(size_t)");
   return o.str();
}
#endif

std::string stringify(int i) {
   std::ostringstream o;
   if (!(o << i))
      throw std::runtime_error("stringify(int)");
   return o.str();
}

std::string stringify(unsigned int i) {
   std::ostringstream o;
   if (!(o << i))
      throw std::runtime_error("stringify(unsigned int)");
   return o.str();
}

int test_internal() {

   int status = 1;

   try { 
      status = test_alt_conf_rotamers();
   }
   catch (std::runtime_error mess) {
      std::cout << "Failed test_alt_conf_rotamers(). " << mess.what() << std::endl;
      status = 0; 
   }
   if (status) {
      try {
	 status = test_wiggly_ligands();
      }
      catch (std::runtime_error mess) {
	 std::cout << "Failed test_wiggly_ligands(). " << mess.what() << std::endl;
	 status = 0; 
      }
   }

   return status; 
}

int test_alt_conf_rotamers() {

   int status = 1;

   std::string filename = greg_test("tutorial-modern.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(filename);
   bool ifound = 0;

   int imod = 1;
   if (atom_sel.read_success > 0) { 
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id = chain_p->GetChainID();
	 if (chain_id == "B") {
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int resno = residue_p->GetSeqNum();
	       if (resno == 72) {

		  ifound = 1;
		  coot::primitive_chi_angles prim_chis(residue_p);
		  std::vector<coot::alt_confed_chi_angles> chis = prim_chis.get_chi_angles();
		  if (chis.size() != 2) {
		     std::string mess = "chis.size() is ";
		     mess += stringify(int(chis.size()));
		     throw std::runtime_error(mess);
		  }

		  int n_rots_found = 0;
		  for (unsigned int i_rot=0; i_rot<2; i_rot++) {
// 		     std::cout << "DEBUG:: chi: " << chis[i_rot].alt_conf << " " 
// 			       << chis[i_rot].chi_angles[0].first  << " "
// 			       << chis[i_rot].chi_angles[0].second << std::endl;
		     
		     if (chis[i_rot].alt_conf == "A") {
			float chi = chis[i_rot].chi_angles[0].second;
			if (chi > 60.0)
			   if (chi < 61.0)
			      n_rots_found++;
		     }
		     if (chis[i_rot].alt_conf == "B") {
			float chi = chis[i_rot].chi_angles[0].second;
			if (chi > -75.0)
			   if (chi < -74.0)
			      n_rots_found++;
		     }
		  }
		  if (n_rots_found != 2) {
		     std::string mess = " found only ";
		     mess += stringify(n_rots_found);
		     mess += " rotamers ";
		     throw std::runtime_error(mess);
		  }
	       }

	       if (resno == 80) {
		  coot::primitive_chi_angles prim_chis(residue_p);
		  std::vector<coot::alt_confed_chi_angles> chis = prim_chis.get_chi_angles();
		  if (chis.size() != 1) {
		     std::string mess = "chis.size() is ";
		     mess += stringify(int(chis.size()));
		     mess += " for ";
		     mess += stringify(79);
		     throw std::runtime_error(mess);
		  }
		  int n_rots_found = 0;
		  if (chis[0].alt_conf == "") {
		     float chi_1 = chis[0].chi_angles[0].second;
		     float chi_2 = chis[0].chi_angles[1].second;
		     if (chi_1 > -58.0)
			if (chi_1 < -57.0)
			   if (chi_2 > 95.0)
			      if (chi_2 < 96.0)
			   n_rots_found++;
		  }
		  if (n_rots_found != 1) {
		     std::string mess = " Oops found ";
		     mess += stringify(n_rots_found);
		     mess += " rotamers ";
		     mess += " for ";
		     mess += stringify(79);
		     throw std::runtime_error(mess);
		  }
	       } 
	    }
	 }
      }
   }

   if (! ifound) {
      std::string mess = "file not found ";
      mess += filename;
      throw std::runtime_error(mess);
   } 

   return status; 
}

int test_wiggly_ligands () {

   // This is not a proper test.  It just exercises some functions.  I
   // had set wiggly_ligand_n_samples to 100 and checked that the
   // distribution of chi angles looked reasonable (it did).
   // Analysing the chi angles can be done another day.
   //
   // BUA looks great.  So why does 3GP's wiggly ligand look terrible?
   // Not sure, but I think it is because there are torsions going on
   // in the ribose (torsionable bonds (not const - they are filtered
   // out) that shouldn't be - and the ring is breaking).
   
   int r = 1;
   coot::protein_geometry geom;
   std::string cif_file_name = greg_test("libcheck_BUA.cif");
   int geom_stat = geom.init_refmac_mon_lib(cif_file_name, 0);
   if (geom_stat == 0) {
      std::string m = "Critical cif dictionary reading failure.";
      std::cout << m << std::endl;
      throw std::runtime_error(m);
   }
   coot::wligand wlig;
   coot::minimol::molecule mmol;
   mmol.read_file(greg_test("monomer-BUA.pdb"));
   int wiggly_ligand_n_samples = 10;
   try { 
      bool optim_geom = 0;
      bool fill_vec = 1;
      std::vector<coot::minimol::molecule> ms = 
	 wlig.install_simple_wiggly_ligands(&geom, mmol, wiggly_ligand_n_samples,
					    optim_geom, fill_vec);

      // jump out if no returned molecules
      if (ms.size() != wiggly_ligand_n_samples) {
	 std::cout << "FAIL: ms.size() != wiggly_ligand_n_samples " << ms.size() << " "
		   << wiggly_ligand_n_samples << std::endl;
	 return 0;
      }

      //
      for (unsigned int imol=0; imol<ms.size(); imol++) {
	 std::string file_name = "test-wiggly-ligand-";
	 file_name += stringify(imol);
	 file_name += ".pdb";
	 ms[imol].write_file(file_name, 10.0);
      }
   }
   catch (std::runtime_error mess) {
      std::cout << mess.what() << std::endl;
   } 
   return r;
}


#endif // BUILT_IN_TESTING
