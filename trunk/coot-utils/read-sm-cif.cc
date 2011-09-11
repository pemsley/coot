
#include <iostream>
#include <string>
#include <stdexcept>

#include <mmdb/mmdb_manager.h>
#include <clipper/core/clipper_util.h>
#include <clipper/core/spacegroup.h>
#include "clipper/core/clipper_instance.h" // tidy up space group cache

#include "coot-utils.hh"

#include "read-sm-cif.hh"

// This can throw a std::runtime_error.
// 
clipper::Cell
coot::smcif::get_cell(PCMMCIFData data) const {

   pstr cell_a = NULL;
   pstr cell_b = NULL;
   pstr cell_c = NULL;
   pstr cell_alpha = NULL;
   pstr cell_beta  = NULL;
   pstr cell_gamma = NULL;
   
   int ierr = 0;
   ierr += data->GetString (cell_a,     "" ,"_cell_length_a");
   ierr += data->GetString (cell_b,     "" ,"_cell_length_b");
   ierr += data->GetString (cell_c,     "" ,"_cell_length_c");
   ierr += data->GetString (cell_alpha, "" ,"_cell_angle_alpha");
   ierr += data->GetString (cell_beta,  "" ,"_cell_angle_beta");
   ierr += data->GetString (cell_gamma, "" ,"_cell_angle_gamma");

   clipper::Cell cell;

   if (! ierr) {
      if (0)
	 std::cout << "make cell from " 
		   << cell_a << " " 
		   << cell_b << " " 
		   << cell_c << " " 
		   << cell_alpha << " " 
		   << cell_beta  << " " 
		   << cell_gamma << " " 
		   << std::endl;
      std::vector<std::string> a_v     = coot::util::split_string_no_blanks(cell_a, "(");
      std::vector<std::string> b_v     = coot::util::split_string_no_blanks(cell_b, "(");
      std::vector<std::string> c_v     = coot::util::split_string_no_blanks(cell_c, "(");
      std::vector<std::string> alpha_v = coot::util::split_string_no_blanks(cell_alpha, "(");
      std::vector<std::string> beta_v  = coot::util::split_string_no_blanks(cell_beta,  "(");
      std::vector<std::string> gamma_v = coot::util::split_string_no_blanks(cell_gamma, "(");

      double a     = coot::util::string_to_float(a_v[0]);
      double b     = coot::util::string_to_float(b_v[0]);
      double c     = coot::util::string_to_float(c_v[0]);
      double alpha = coot::util::string_to_float(alpha_v[0]);
      double beta  = coot::util::string_to_float( beta_v[0]);
      double gamma = coot::util::string_to_float(gamma_v[0]);
      clipper::Cell_descr cell_descr(a,b,c,
				     clipper::Util::d2rad(alpha),
				     clipper::Util::d2rad(beta),
				     clipper::Util::d2rad(gamma));
      clipper::Cell cell(cell_descr);
   }
   // Oh dear, we are returning an empty cell, maybe sometimes
   return cell;
}


// can throw an clipper::Message_base exception
// 
clipper::Spacegroup
coot::smcif::get_space_group(const std::vector<std::string> &symm_strings) const {

   std::string symmetry_ops;
   for (unsigned int isym=0; isym<symm_strings.size(); isym++) { 
      symmetry_ops += symm_strings[isym];
      symmetry_ops += " ; ";
   }
   clipper::Spacegroup space_group;
   space_group.init(clipper::Spgr_descr(symmetry_ops, clipper::Spgr_descr::Symops));
   std::cout << "DEBUG:: space group initialised with symbol \""
	     << space_group.symbol_hm() << "\"" << std::endl;
   return space_group;
}

std::vector<CAtom *>
coot::smcif::read_coordinates(PCMMCIFData data, const clipper::Cell &cell, const clipper::Spacegroup &spg) const {

   std::vector<CAtom *> atom_vec;
   const char *loopTagsAtom[6] = { "_atom_site_label",
			       "_atom_site_type_symbol",
			       "_atom_site_fract_x",
			       "_atom_site_fract_y",
			       "_atom_site_fract_z", ""};
   const char *loopTagsAniso[2] = { "_atom_site_aniso_label", ""};


   int ierr = 0;
   pstr S = NULL;
   PCMMCIFLoop loop = data->FindLoop((pstr *) loopTagsAtom);
   if (loop) {
      int ll = loop->GetLoopLength();
      if (ll >= 0) {

	 char *symbol = NULL;
	 char *label  = NULL;
	 realtype xf,yf,zf, occ, tf;
	 realtype x,y,z;
	 int ierr_tot = 0;
	 
	    
	 for (unsigned int il=0; il<ll; il++) {

	    label  = loop->GetString(loopTagsAtom[0], il, ierr);
	    ierr_tot += ierr;
	    symbol = loop->GetString(loopTagsAtom[1], il, ierr);
	    ierr_tot += ierr;
	    loop->GetReal(xf, loopTagsAtom[2], il, ierr);
	    ierr_tot += ierr;
	    loop->GetReal(yf, loopTagsAtom[3], il, ierr);
	    ierr_tot += ierr;
	    loop->GetReal(zf, loopTagsAtom[4], il, ierr);
	    ierr_tot += ierr;

	    occ = 1; // hack
	    tf = 10.0;

	    // can we get a real value for occ?
	    int ierr_tf = 0;
	    loop->GetReal(tf, "_atom_site_U_iso_or_equiv", il, ierr_tf);
	    if (! ierr_tf) {
	       tf *= 8 * M_PI * M_PI; // PDB-scaled
	    } 

	    int ierr_occ = 0;
	    loop->GetReal(occ, "_atom_site_occupancy", il, ierr_occ);

	    
	    if (ierr_tot == 0) {
	       CAtom *at = new CAtom;
	       clipper::Coord_frac cf(xf,yf,zf);
	       std::cout << " found atom: " << label << " " << symbol << " "
			 << cf.format() << std::endl;
	       clipper::Coord_orth co = cf.coord_orth(cell);
	       at->SetCoordinates(co.x(), co.y(),co.z(), occ, tf);
	       // label -> 4c atom name conversion? 
	       at->SetAtomName(label);
	       at->SetElementName(symbol);
	       at->Het = 1; // all SM cifs atoms are HETATMs :)
	       atom_vec.push_back(at);
	    } else {
	       if (0) 
		  std::cout << "reject atom at loop count " << il << std::endl;
	    } 
	 }
      }
   }

   // Aniso atoms
   //
   std::vector<coot::simple_sm_u> u_aniso_vec;
   
   loop = data->FindLoop((pstr *) loopTagsAniso);
   if (loop) {
      int ll = loop->GetLoopLength();
      char *label  = NULL;
      realtype u11, u22, u33, u12, u13, u23;
      int ierr_tot = 0;
      for (unsigned int il=0; il<ll; il++) {
	 label  = loop->GetString(loopTagsAniso[0], il, ierr);
	 ierr_tot += ierr;
	 loop->GetReal(u11, "_atom_site_aniso_U_11", il, ierr);
	 ierr_tot += ierr;
	 loop->GetReal(u22, "_atom_site_aniso_U_22", il, ierr);
	 ierr_tot += ierr;
	 loop->GetReal(u33, "_atom_site_aniso_U_33", il, ierr);
	 ierr_tot += ierr;
	 loop->GetReal(u12, "_atom_site_aniso_U_12", il, ierr);
	 ierr_tot += ierr;
	 loop->GetReal(u13, "_atom_site_aniso_U_13", il, ierr);
	 ierr_tot += ierr;
	 loop->GetReal(u23, "_atom_site_aniso_U_23", il, ierr);
	 ierr_tot += ierr;

	 if (! ierr_tot) {
	    // label -> atom name conversion here?
	    coot::simple_sm_u smu(label, u11, u22, u33, u12, u13, u23);
	    u_aniso_vec.push_back(smu);
	 }
      }

      // now put those aniso Us into the atom_vec;
      //
      double a = cell.a();
      double b = cell.b();
      double c = cell.c();
      for (unsigned int ianiso=0; ianiso<u_aniso_vec.size(); ianiso++) { 
	 for (unsigned int iat=0; iat<atom_vec.size(); iat++) {
	    CAtom *at = atom_vec[iat];
	    if (u_aniso_vec[ianiso].label == std::string(at->GetAtomName())) {
	       clipper::U_aniso_frac caf(u11/(a*a), u22/(b*b), u33/(c*c),
					 u12/(a*b), u13/(a*c), u23/(b*c));
	       clipper::U_aniso_orth cao = caf.u_aniso_orth(cell);
	       at->u11 = cao(0,0);
	       at->u22 = cao(1,1);
	       at->u33 = cao(2,2);
	       at->u12 = cao(0,1);
	       at->u13 = cao(0,2);
	       at->u23 = cao(1,2);
	       at->WhatIsSet |= ASET_Anis_tFac; // is anisotropic
	    } 
	 }
      }
   } 
   return atom_vec;
}


CMMDBManager *
coot::smcif::read_sm_cif(const std::string &file_name) const {

   CMMDBManager *mol = NULL;
   pstr S = NULL;
   PCMMCIFData data = new CMMCIFData();
   data->SetFlag (CIFFL_SuggestCategories);
   int ierr = data->ReadMMCIFData (file_name.c_str());
   if (!ierr) {

// testing      
//       int ierr = data->GetString (S, "" ,"_chemical_formula_sum");
//       if (! ierr) { 
// 	 printf("chemical-formula-sum: %s\n", S);
//       } else {
// 	 printf("error getting chemical-formula-sum string.\n");
//       } 

      ierr = data->GetString (S, "", "_[local]_cod_chemical_formula_sum_orig");
      if (!ierr)
	 printf("_[local]_cod_chemical_formula_sum_orig: %s\n", S);

      try { 
	 clipper::Cell cell = get_cell(data);
	 std::cout << "got cell: " << cell.format() << std::endl;

	 std::vector<std::string> symm_strings;
	 const char *loopTag1[2] = { "_symmetry_equiv_pos_as_xyz",
				     ""};

	 PCMMCIFLoop loop = data->FindLoop((pstr *) loopTag1);
	 if (loop) {
	    int ll = loop->GetLoopLength();
	    // std::cout << "loop length: " << ll << std::endl;
	    if (ll > 0) { 
	       for (unsigned int il=0; il<ll; il++) {

		  S = loop->GetString(loopTag1[0], il, ierr);
		  if (! ierr) {
		     // std::cout << "symmetry: " << S << std::endl;
		     symm_strings.push_back(S);
		  } else {
		     std::cout << "error symmetry-equiv-pos-as-xyz string.\n";
		  } 
	       }
	    } 
	    if (symm_strings.size()) {
	       try { 
		  clipper::Spacegroup spg = get_space_group(symm_strings);

		  std::vector<CAtom *> atoms = read_coordinates(data, cell, spg);
		  std::cout << "read " << atoms.size() << " atoms" << std::endl;

		  if (atoms.size()) {

		     mol = new CMMDBManager;
		     CModel *model_p = new CModel;
		     CChain *chain_p = new CChain;
		     CResidue *residue_p = new CResidue;
		     chain_p->SetChainID("");
		     residue_p->seqNum = 1;
		     residue_p->SetResName("XXX");
		     for (unsigned int iat=0; iat<atoms.size(); iat++)
			residue_p->AddAtom(atoms[iat]);
		     chain_p->AddResidue(residue_p);
		     model_p->AddChain(chain_p);
		     mol->AddModel(model_p);
		  } 
	       } 
	       catch (clipper::Message_base exc) {
		  std::cout << "Oops, trouble.  No such spacegroup\n";
	       }
	    } else {
	       std::cout << "ERROR:: no symm strings" << std::endl;
	    } 
	 } else {
	    std::cout << "No symmetry loop" << std::endl;
	 }
	 
      }

      catch (std::runtime_error rte) {
	 std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }
      
   delete data;
   delete S;
   
   return mol;
} 

// int main(int argc, char **argv) {

//    if (argc > 1) {
//       InitMatType(); // delete me when not stand-alone
//       std::string file_name = argv[1];
//       coot::smcif smcif;
//       CMMDBManager *mol = smcif.read_sm_cif(file_name);
//    }
//    clipper::ClipperInstantiator::instance().destroy();
//    return 0;
// } 

