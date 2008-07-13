
#include "testing.hh"

#ifdef BUILT_IN_TESTING

#include <iostream>
#include <sstream>

#include "clipper/core/ramachandran.h"

#include "coot-coord-utils.hh"
#include "coot-rama.hh"
#include "primitive-chi-angles.hh"

#include "wligand.hh"
#include "simple-restraint.hh"

// a shorthand so that the push back line doesn't get too long:
typedef std::pair<int(*)(), std::string> named_func;

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

bool close_double(double a, double b) {
   double small_diff = 0.001;
   return (fabs(a-b) < small_diff);
}

bool close_pair_test(std::pair<double, double> a, std::pair<double, double> b) {
   return (close_double(a.first, b.first) && close_double(a.second, b.second));
}

std::ostream& operator<< (std::ostream &s, std::pair<double, double> b) {
   s << "[" << b.first << "," << b.second << "]";
   return s;
} 

int test_internal() {

   int status = 1;
   std::vector<named_func> functions;
   functions.push_back(named_func(test_alt_conf_rotamers, "test_alt_conf_rotamers"));
   functions.push_back(named_func(test_wiggly_ligands, "test_wiggly_ligands"));
   functions.push_back(named_func(test_ramachandran_probabilities, "test_ramachandran_probabilities"));

   for (unsigned int i_func=0; i_func<functions.size(); i_func++) {
      std::cout << "Entering test: " << functions[i_func].second << std::endl;
      try { 
	 status = functions[i_func].first();
	 if (status == 0) 
	    break;
      }
      catch (std::runtime_error mess) {
	 std::cout << "Failed " << functions[i_func].second << " " << mess.what() << std::endl;
	 status = 0;
	 break;
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
   atom_sel.clear_up();

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

// Return a new mol, and a residue selection.  Delete the residue
// selection and mol when you are done with it.
// 
residue_selection_t
test_ramachandran_probabilities_refine_fragment(atom_selection_container_t atom_sel,
						PCResidue *SelResidues,
						int nSelResidues,
						const std::string &chain_id,
						int resno_mid,
						coot::protein_geometry geom,
						bool enable_rama_refinement) {

   
   // now refine a bit of structure:
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   short int have_flanking_residue_at_start = 0;
   short int have_flanking_residue_at_end = 0;
   short int have_disulfide_residues = 0;  // other residues are included in the
   std::string altconf = "";
   short int in_alt_conf_split_flag = 0;
   char *chn = (char *) chain_id.c_str(); // mmdb thing.  Needs updating on new mmdb?
	       
   std::pair<CMMDBManager *, int> residues_mol_pair = 
      coot::util::create_mmdbmanager_from_res_selection(atom_sel.mol,
							SelResidues, nSelResidues, 
							have_flanking_residue_at_start,
							have_flanking_residue_at_end,
							altconf,
							chain_id,
							in_alt_conf_split_flag);
   
   coot::restraints_container_t restraints(resno_mid-1,
					   resno_mid+1,
					   have_flanking_residue_at_start,
					   have_flanking_residue_at_end,
					   have_disulfide_residues,
					   altconf,
					   chn,
					   residues_mol_pair.first,
					   fixed_atom_specs);

   short int do_rama_restraints = 0;
   short int do_residue_internal_torsions = 0;
   short int do_link_torsions = 1;
   float rama_plot_restraint_weight = 1.0;
	       
   coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
   short int do_peptide_torsion_restraints = coot::restraints_container_t::NO_LINK_TORSION;
   if (enable_rama_refinement) { 
      do_peptide_torsion_restraints = coot::restraints_container_t::LINK_TORSION_RAMACHANDRAN_GOODNESS;
      do_rama_restraints = 1;
      flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
   } 

   std::cout << " ===================== enable_rama_refinement: " << enable_rama_refinement
	     << " " << do_rama_restraints << " ========================" << std::endl;
   
   coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
   int nrestraints = 
      restraints.make_restraints(geom, flags,
				 do_residue_internal_torsions,
				 do_link_torsions,
				 rama_plot_restraint_weight,
				 do_rama_restraints,
				 pseudos);
   restraints.minimize(flags);

   int post_refine_selHnd = residues_mol_pair.first->NewSelection();
   int post_refine_nSelResidues; 
   PCResidue *post_refine_SelResidues = NULL;
   residues_mol_pair.first->Select(post_refine_selHnd, STYPE_RESIDUE, 0,
				   chn,
				   resno_mid-1, "",
				   resno_mid+1, "",
				   "*",  // residue name
				   "*",  // Residue must contain this atom name?
				   "*",  // Residue must contain this Element?
				   "",   // altLocs
				   SKEY_NEW // selection key
				   );
   residues_mol_pair.first->GetSelIndex(post_refine_selHnd,
					post_refine_SelResidues,
					post_refine_nSelResidues);

   residue_selection_t res_sel;
   res_sel.mol = residues_mol_pair.first;
   res_sel.SelectionHandle = post_refine_selHnd;
   res_sel.nSelResidues = post_refine_nSelResidues;
   res_sel.SelResidues = post_refine_SelResidues;
   return res_sel;
} 

int test_ramachandran_probabilities() {

   int r = 0;

   std::string file_name = greg_test("crashes_on_cootaneering.pdb");
   atom_selection_container_t atom_sel = get_atom_selection(file_name);

   if (! atom_sel.read_success)
      throw std::runtime_error(file_name + ": file not found.");


   std::string chain_id = "A";
   char *chn = (char *) chain_id.c_str(); // mmdb thing.  Needs updating on new mmdb?
   std::vector<int> resnos;
   resnos.push_back(12);  // fail
   resnos.push_back(14);  // phi=-57.7411  psi=-32.0451
   resnos.push_back(15);  // phi=-53.7037  psi=-47.837
   resnos.push_back(16);  // phi=-57.4047  psi=-42.6739

   coot::protein_geometry geom;
   geom.init_standard();
   int n_correct = 0; 
   for (int i=0; i<resnos.size(); i++) { 
      int selHnd = atom_sel.mol->NewSelection();
      int nSelResidues; 
      PCResidue *SelResidues = NULL;
      atom_sel.mol->Select(selHnd, STYPE_RESIDUE, 0,
			   chn,
			   resnos[i]-1, "",
			   resnos[i]+1, "",
			   "*",  // residue name
			   "*",  // Residue must contain this atom name?
			   "*",  // Residue must contain this Element?
			   "",   // altLocs
			   SKEY_NEW // selection key
			   );
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

      try { 
	 coot::util::phi_psi_t angles = coot::util::ramachandran_angles(SelResidues, nSelResidues);
	 for (int ires=0; ires<3; ires++) 
	    geom.try_dynamic_add(SelResidues[ires]->GetResName(), ires);

	 if (i > 0) { 
	    std::pair<double,double> angles_pair(angles.phi(), angles.psi());
	    std::pair<double, double> expected;
	    if (i==1)
	       expected = std::pair<double, double> (-57.7411, -32.0451);
	    if (i==2)
	       expected = std::pair<double, double> (-53.7037, -47.837);
	    if (i==3)
	       expected = std::pair<double, double> (-57.4047, -42.6739);
	 
	    bool close = close_pair_test(angles_pair, expected);
	    if (! close) {
	       std::cout << " Fail on Ramachandran angle test " << angles << " should be "
			 << expected << std::endl;
	       r = 0;
	       break;
	    } else {

	       // get the probability
	       std::string residue_type = SelResidues[1]->GetResName();

	       clipper::Ramachandran::TYPE rama_type = clipper::Ramachandran::NonGlyPro;
	       if (residue_type == "GLY")
		  rama_type = clipper::Ramachandran::Gly;
	       if (residue_type == "PRO")
		  rama_type = clipper::Ramachandran::Pro;
	       clipper::Ramachandran rama(rama_type);

	       double prob = rama.probability(clipper::Util::d2rad(angles.phi()),
					      clipper::Util::d2rad(angles.psi()));

	       int enable_rama_refinement = 0;
	       residue_selection_t refined_res_sel =
		  test_ramachandran_probabilities_refine_fragment(atom_sel, SelResidues, nSelResidues,
								  chain_id, resnos[i], geom,
								  enable_rama_refinement);
	       

	       if (0) { 
		  // Let's look at the refined structure. Write them out as pdb files ;-/
		  std::string tmp_file_name = "rama-test-";
		  tmp_file_name += coot::util::int_to_string(i);
		  tmp_file_name += ".pdb";
		  refined_res_sel.mol->WritePDBASCII((char *)tmp_file_name.c_str());
	       }

	       coot::util::phi_psi_t post_refine_angles =
		  coot::util::ramachandran_angles(refined_res_sel.SelResidues, refined_res_sel.nSelResidues);
	       refined_res_sel.mol->DeleteSelection(refined_res_sel.SelectionHandle);
	       delete refined_res_sel.mol;
	       refined_res_sel.mol = 0;
	       
	       double post_refine_prob =
		  rama.probability(clipper::Util::d2rad(post_refine_angles.phi()),
				   clipper::Util::d2rad(post_refine_angles.psi()));

	       // now with Ramachandran refinement:
	       //
	       enable_rama_refinement = 1;
	       residue_selection_t rama_refined_res_sel =
		  test_ramachandran_probabilities_refine_fragment(atom_sel, SelResidues, nSelResidues,
								  chain_id, resnos[i], geom,
								  enable_rama_refinement);
	       coot::util::phi_psi_t rama_refine_angles =
		  coot::util::ramachandran_angles(rama_refined_res_sel.SelResidues,
						  rama_refined_res_sel.nSelResidues);
	       rama_refined_res_sel.mol->DeleteSelection(rama_refined_res_sel.SelectionHandle);
	       delete rama_refined_res_sel.mol;
	       rama_refined_res_sel.mol = 0;
	       
	       double rama_refine_prob =
		  rama.probability(clipper::Util::d2rad(rama_refine_angles.phi()),
				   clipper::Util::d2rad(rama_refine_angles.psi()));
	       std::cout << "--------------------------------------\n";
	       std::cout << "Pre-refine         Rama probability residue " << resnos[i] << ": "
			 << prob << std::endl;
	       std::cout << "Post-simple refine Rama probability residue " << resnos[i] << ": "
			 << post_refine_prob << std::endl;
	       std::cout << "Post-Rama   refine Rama probability residue " << resnos[i] << ": "
			 << rama_refine_prob << std::endl;
	       std::cout << "--------------------------------------\n";

	       // 5% better probability needed.
	       // 
	       if (rama_refine_prob > (post_refine_prob*1.05))
		  n_correct++;
	    }
	 }
      }
      catch (std::runtime_error mess) {
	 if (i==0) { 
	    // std::cout << "resno set " << i << " was correct " << std::endl;
	    n_correct++;
	 }
      } 

      atom_sel.mol->DeleteSelection(selHnd);
   }

   if (n_correct != 4) {
      std::cout << "Failed to get 4 rama angles improvements " << std::endl;
      r = 0;
   } else {
      r = 1;  // return success.
      // std::cout << "n_correct is 4" << std::endl;
   } 

   return r;
} 


#endif // BUILT_IN_TESTING
