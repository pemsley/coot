
// The problem with carbohydrate-building - as it stands.
//
// The set of reference structures has pyranose-BETA1-4-pyranose.  And
// we paste a MAN onto the extension poart of that to give us initial
// coordinates for the MAN. The link consists of part pre-built base
// residue and part of this LSQ fitted MAN - and because the C1, O4 C4
// are not in the same positions as the reference
// pyranose-BETA1-4-pyranose, this damages the bonds and angle of the
// glycosidic linkage.
//
// Instead, we need to build the MAN by torsions from the base
// residue.
//
// What do we need to do this?
//
// 1: For each link type, how to build the core 6 atoms 1: For each
// carbohydrate how to build the approproate decorations.

#include <iomanip>
#include <fstream>
#include "glyco-z.hh"

std::ostream& operator<<(std::ostream &o, const atom_by_torsion_t &abt) {

   o << "atom " << abt.atom_name << " " 
     << abt.element << " based-on "
     << abt.prior_atom_1.first << " "
     << abt.prior_atom_1.second << " "
     << abt.prior_atom_2.first << " "
     << abt.prior_atom_2.second << " "
     << abt.prior_atom_3.first << " "
     << abt.prior_atom_3.second << " "
     << "bond-length: ";
   o.setf(std::ios::fixed, std::ios::floatfield);
   o.precision(5);
   o << abt.bond_length
     << " angle: " << abt.angle
     << " tors: " << std::setw(10) << abt.torsion;
   return o;
} 


clipper::Coord_orth
atom_by_torsion_t::pos(CResidue *base_residue_p, CResidue *ext_residue_p) const {

   CAtom *at_1 = NULL;
   CAtom *at_2 = NULL;
   CAtom *at_3 = NULL;

   if (0) { 
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      base_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) { 
	 CAtom *at = residue_atoms[iat];
	 std::cout << "   " << coot::atom_spec_t(at) << " vs "
		   << prior_atom_1.first << " " << prior_atom_1.second << std::endl;
      }
   }
      
   if (prior_atom_1.first)
      at_1 = base_residue_p->GetAtom(prior_atom_1.second.c_str());
   else 
      at_1 = ext_residue_p->GetAtom(prior_atom_1.second.c_str());
   if (prior_atom_2.first)
      at_2 = base_residue_p->GetAtom(prior_atom_2.second.c_str());
   else 
      at_2 = ext_residue_p->GetAtom(prior_atom_2.second.c_str());
   if (prior_atom_3.first)
      at_3 = base_residue_p->GetAtom(prior_atom_3.second.c_str());
   else 
      at_3 = ext_residue_p->GetAtom(prior_atom_3.second.c_str());

   if (at_1 && at_2 && at_3) {
      clipper::Coord_orth p1 = coot::co(at_1);
      clipper::Coord_orth p2 = coot::co(at_2);
      clipper::Coord_orth p3 = coot::co(at_3);
      clipper::Coord_orth new_pos = clipper::Coord_orth(p3,p2,p1,
							bond_length,
							clipper::Util::d2rad(angle),
							clipper::Util::d2rad(torsion));
      return new_pos;
   } else {
      std::string m = "missing atom";
      if (! at_1) m += " at_1";
      if (! at_2) m += " at_2";
      if (! at_3) m += " at_3";
      throw std::runtime_error(m);
   }
}
   
CResidue *link_by_torsion_t::make_residue(CResidue *base_residue_p) const {

   CResidue *r = NULL;
   if (geom_atom_torsions.size()) {
      r = new CResidue;
      r->SetResName(new_residue_type.c_str());
      r->seqNum = new_res_no;; 
      for (unsigned int i=0; i<geom_atom_torsions.size(); i++) {
	 const atom_by_torsion_t &gat = geom_atom_torsions[i];
	 clipper::Coord_orth p = geom_atom_torsions[i].pos(base_residue_p, r);
	 CAtom *atom = new CAtom(r); // does an add atom
	 atom->SetAtomName(gat.atom_name.c_str());
	 atom->SetElementName(gat.element.c_str());
	 atom->SetCoordinates(p.x(), p.y(), p.z(), 1.0, 30);
	 if (0) 
	    std::cout << "   " << gat.atom_name << " " << p.format()  << std::endl;
      }
   } 
   return r;
} 


atom_by_torsion_t get_atom_by_torsion(const atom_by_torsion_base_t &names,
				      CResidue *residue_1_p,  // reference/lower
				      CResidue *residue_2_p   // extension residue
				      ) {
   atom_by_torsion_t abt;
   if (0) 
      std::cout << "get "
		<< names.prior_atom_1.first << " " << names.prior_atom_1.second << " "
		<< names.prior_atom_2.first << " " << names.prior_atom_2.second << " "
		<< names.prior_atom_3.first << " " << names.prior_atom_3.second << " "
		<< std::endl;
   PPCAtom residue_atoms_1 = 0;
   PPCAtom residue_atoms_2 = 0;
   int n_residue_atoms_1;
   int n_residue_atoms_2;
   residue_1_p->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   residue_2_p->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

   CAtom *p_new = residue_2_p->GetAtom(names.atom_name.c_str());
   if (p_new) { 
      CAtom *p_1 = NULL;
      CAtom *p_2 = NULL;
      CAtom *p_3 = NULL;

      // I could just GetAtom() here.
      for (unsigned int iat=0; iat<n_residue_atoms_1; iat++) { 
	 CAtom *at = residue_atoms_1[iat];
	 std::string nb_name = coot::util::remove_whitespace(at->name);
	 if (names.prior_atom_1.first)
	    if (names.prior_atom_1.second == nb_name)
	       p_1 = at;
	 if (names.prior_atom_2.first)
	    if (names.prior_atom_2.second == nb_name)
	       p_2 = at;
	 if (names.prior_atom_3.first)
	    if (names.prior_atom_3.second == nb_name)
	       p_3 = at;
      }
      for (unsigned int iat=0; iat<n_residue_atoms_2; iat++) { 
	 CAtom *at = residue_atoms_2[iat];
	 std::string nb_name = coot::util::remove_whitespace(at->name);
	 if (! names.prior_atom_1.first)
	    if (names.prior_atom_1.second == nb_name)
	       p_1 = at;
	 if (! names.prior_atom_2.first)
	    if (names.prior_atom_2.second == nb_name)
	       p_2 = at;
	 if (! names.prior_atom_3.first)
	    if (names.prior_atom_3.second == nb_name)
	       p_3 = at;
      }

      if (p_1 && p_2 && p_3) {
	 coot::atom_quad q(p_new, p_1, p_2, p_3);
	 if (0) { 
	    std::cout << "check quad: " << q << " " << q.get_atom_name_quad() << std::endl;
	    std::cout << "   check quad "
		      << coot::atom_spec_t(q.atom_1) << " "
		      << coot::atom_spec_t(q.atom_2) << " "
		      << coot::atom_spec_t(q.atom_3) << " "
		      << coot::atom_spec_t(q.atom_4) << std::endl;
	    }
	 clipper::Coord_orth pos_n = coot::co(p_new);
	 clipper::Coord_orth pos_1 = coot::co(p_1);
	 clipper::Coord_orth pos_2 = coot::co(p_2);
	 clipper::Coord_orth pos_3 = coot::co(p_3);
	 if (0) { 
	    std::cout << "... " << coot::atom_spec_t(p_new) << " has pos " << pos_n.format() << std::endl;
	    std::cout << "... " << coot::atom_spec_t(p_1)   << " has pos " << pos_1.format() << std::endl;
	    }
	 double bl = sqrt((pos_n - pos_1).lengthsq());
	 double a = q.angle_2();
	 double t = q.torsion();
	 abt = atom_by_torsion_t(names, bl, a, t);
      }
   }
   return abt;
}

   


// Alpha 1-6: BMA-MAN
// C1: O6p C6p C5p
// C2: C1  O6p C6p
// O2: C2  C1  O6p 
// C3: C2  C1  O6p 
// O3: C3  C2  C1  
// C4: C3  C2  C1  
// O4: C4  C3  C2
// C5: C4  C3  C2
// O5: C5  C4  C3
// C6: C5  C4  C3
// O6: C6  C5  C4


link_by_torsion_base_t pyranose_link_1_6_to_core() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   ats.push_back(atom_by_torsion_base_t("C1", "C", BS(true,  "O6"), BS(true,  "C6"), BS(true,  "C5")));
   ats.push_back(atom_by_torsion_base_t("C2", "C", BS(false, "C1"), BS(true,  "O6"), BS(true,  "C6")));
   ats.push_back(atom_by_torsion_base_t("C3", "C", BS(false, "C2"), BS(false, "C1"), BS(true,  "O6")));
   ats.push_back(atom_by_torsion_base_t("C4", "C", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("C5", "C", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("O5", "O", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

link_by_torsion_base_t pyranose_link_1_4_to_core() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   ats.push_back(atom_by_torsion_base_t("C1", "C", BS(true,  "O4"), BS(true,  "C4"), BS(true,  "C3")));
   ats.push_back(atom_by_torsion_base_t("C2", "C", BS(false, "C1"), BS(true,  "O4"), BS(true,  "C4")));
   ats.push_back(atom_by_torsion_base_t("C3", "C", BS(false, "C2"), BS(false, "C1"), BS(true,  "O4")));
   
   ats.push_back(atom_by_torsion_base_t("C4", "C", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("C5", "C", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("O5", "O", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

link_by_torsion_base_t pyranose_link_1_2_to_core() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   ats.push_back(atom_by_torsion_base_t("C1", "C", BS(true,  "O2"), BS(true,  "C2"), BS(true,  "C1")));
   ats.push_back(atom_by_torsion_base_t("C2", "C", BS(false, "C1"), BS(true,  "O2"), BS(true,  "C2")));
   ats.push_back(atom_by_torsion_base_t("C3", "C", BS(false, "C2"), BS(false, "C1"), BS(true,  "O2")));
   
   ats.push_back(atom_by_torsion_base_t("C4", "C", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("C5", "C", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("O5", "O", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

link_by_torsion_base_t pyranose_link_1_3_to_core() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   ats.push_back(atom_by_torsion_base_t("C1", "C", BS(true,  "O3"), BS(true,  "C3"), BS(true,  "C2")));
   ats.push_back(atom_by_torsion_base_t("C2", "C", BS(false, "C1"), BS(true,  "O3"), BS(true,  "C3")));
   ats.push_back(atom_by_torsion_base_t("C3", "C", BS(false, "C2"), BS(false, "C1"), BS(true,  "O3")));
   
   ats.push_back(atom_by_torsion_base_t("C4", "C", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("C5", "C", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("O5", "O", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

link_by_torsion_base_t pyranose_link_2_3_to_core() {

   // different - the extending residue keeps its oxygen
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   ats.push_back(atom_by_torsion_base_t("O3", "C", BS(true,  "C2"), BS(true,  "C3"), BS(true,  "C4")));
   ats.push_back(atom_by_torsion_base_t("C3", "C", BS(false, "O3"), BS(true,  "C2"), BS(true,  "C3")));
   ats.push_back(atom_by_torsion_base_t("C2", "C", BS(false, "C3"), BS(false, "O3"), BS(true,  "C2")));

   ats.push_back(atom_by_torsion_base_t("C1", "C", BS(false, "C2"), BS(false, "C3"), BS(false, "O3")));
   ats.push_back(atom_by_torsion_base_t("C4", "C", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("C5", "C", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("O5", "O", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

link_by_torsion_base_t mannose_decorations() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   ats.push_back(atom_by_torsion_base_t("O2", "O", BS(false, "C2"), BS(false, "C1"), BS(true,  "O6")));
   ats.push_back(atom_by_torsion_base_t("O3", "O", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("O4", "O", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("C6", "C", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   ats.push_back(atom_by_torsion_base_t("O6", "O", BS(false, "C6"), BS(false, "C5"), BS(false, "C4")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

// When the link type is xxx1-Y, then we can simply add up the torsions.
// When the link type is xxx2-3, then we have an O3 (which would otherwise be "decoration")
//   that is part of the link atoms.  In that case, one of the O3s should be omitted.
//   The add() function should check if this atom exists already and if so not add the new
//   one.
// 
link_by_torsion_base_t add_link_by_torsions(const link_by_torsion_base_t &v1,
					    const link_by_torsion_base_t &v2) {

   link_by_torsion_base_t r = v1;
   for (unsigned int i=0; i<v2.atom_torsions.size(); i++)
      r.add(v2.atom_torsions[i]);
   return r;
} 

link_by_torsion_base_t
get_names_for_link_type(const std::string &link_type) {
   
   link_by_torsion_base_t r;
   if (link_type == "ALPHA1-6") r = pyranose_link_1_6_to_core();
   if (link_type == "ALPHA1-2") r = pyranose_link_1_2_to_core();
   if (link_type == "ALPHA1-3") r = pyranose_link_1_3_to_core();
   if (link_type == "ALPHA2-3") r = pyranose_link_2_3_to_core();
   if (link_type == "BETA1-4")  r = pyranose_link_1_4_to_core();
   return r;
}
   

std::pair<CResidue *, CResidue *> get_residue_pair(CMMDBManager *mol) {

   std::pair<CResidue *, CResidue *> r(NULL, NULL);
   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 if (r.first) { 
	    r.second = residue_p;
	    break;
	 } else {
	    r.first = residue_p;
	 } 
      }
      if (r.first && r.second)
	 break;
   }
   return r;
}

void write_tree(CResidue *ref_res_p, CResidue *ext_res_p) {

   // Just a test function
   
   link_by_torsion_base_t a16 = add_link_by_torsions(pyranose_link_1_6_to_core(), mannose_decorations());

   for (unsigned int i=0; i<a16.atom_torsions.size(); i++) {
      atom_by_torsion_t abt = get_atom_by_torsion(a16.atom_torsions[i], ref_res_p, ext_res_p);
      if (abt.filled()) {
	 std::cout << abt << std::endl;
      }
   }
} 

void
link_by_torsion_t::init(CResidue *ref_res_p, CResidue *ext_res_p) {

   if (0) 
      std::cout << "in link_by_torsion_t::init() have " << atom_torsions.size()
		<< "  atom torsions" << std::endl;
   for (unsigned int i=0; i<atom_torsions.size(); i++) {
      atom_by_torsion_t abt = get_atom_by_torsion(atom_torsions[i], ref_res_p, ext_res_p);
      if (! abt.filled()) {
	 std::cout << "Missing atom! " << abt << std::endl;
      } else {
	 add(abt);
      } 
   }
}

void
link_by_torsion_t::print() const {

   std::cout << "--- link_by_torsion_t " << geom_atom_torsions.size()
	     << " atom torsions" << std::endl;
   for (unsigned int i=0; i<geom_atom_torsions.size(); i++)
      std::cout << "   " << std::setw(2) << i << " " << geom_atom_torsions[i] << std::endl;
} 
 
void
link_by_torsion_t::write(const std::string &file_name) const {

   std::ofstream f(file_name.c_str());
   if (f)
      for (unsigned int i=0; i<geom_atom_torsions.size(); i++)
	 f << "  "  << " " << geom_atom_torsions[i] << "\n";
}

#include <string>

// read what is written by write() function
void
link_by_torsion_t::read(const std::string &file_name) {

   std::ifstream f(file_name.c_str());

   if (f) {
      std::vector<std::string> lines;
      std::string line;
      while (std::getline(f, line)) { 
	 lines.push_back(line);
      }

      if (lines.size()) {
	 for (unsigned int i=0; i<lines.size(); i++) { 
	    std::vector<std::string> bits = coot::util::split_string_no_blanks(lines[i]);
	    if (bits.size() == 16) {
	       if (bits[0] == "atom") {
		  // std::cout << "parse line " << lines[i] << std::endl;
		  try {
		     std::vector<std::pair<bool, std::string> > a(4);
		     std::string new_atom_name = bits[1];
		     std::string new_atom_ele  = bits[2];
		     a[1].first  = coot::util::string_to_int(bits[4]);
		     a[1].second = bits[5];
		     a[2].first  = coot::util::string_to_int(bits[6]);
		     a[2].second = bits[7];
		     a[3].first  = coot::util::string_to_int(bits[8]);
		     a[3].second = bits[9];
		     double bl_fs      = coot::util::string_to_float(bits[11]);
		     double angle_fs   = coot::util::string_to_float(bits[13]);
		     double torsion_fs = coot::util::string_to_float(bits[15]);

		     atom_by_torsion_base_t ba(new_atom_name, new_atom_ele,
					       a[1], a[2], a[3]);
		     atom_by_torsion_t abt(ba, bl_fs, angle_fs, torsion_fs);
		     add(abt);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "Failed to parse: " << line << std::endl;
		  }
	       }
	    }
	 }
      }
   }
   if (geom_atom_torsions.size() == 0) 
      std::cout << "After read()ing, we have " << geom_atom_torsions.size()
		<< " atom torsions" << std::endl;
}



int old_main(int argc, char **argv) {

   if (argc > 2) {
      std::string file_name = argv[1];
      std::string link_type = argv[2]; // "ALPHA1-6";
      std::string new_residue_type = "MAN";
      
      CMMDBManager *mol = new CMMDBManager;
      int status = mol->ReadPDBASCII(file_name.c_str());
      if (status != Error_NoError) {
	 std::cout << "ERROR:: on reading " << file_name << std::endl;
      } else {
	 std::pair<CResidue *, CResidue *> p = get_residue_pair(mol);
	 if (p.first && p.second) {

	    link_by_torsion_base_t a16_to_core = get_names_for_link_type(link_type);
	    if (! a16_to_core.filled()) {
	       std::cout << "link by torsion core not filled with atom torsions " << std::endl;
	    } else { 
// 	       int new_res_no = 1;
// 	       link_by_torsion_t l("", new_residue_type, a16_to_core, p.first, p.second);
// 	       if (! l.filled()) {
// 		  std::cout << "link by torsion not filled with atom torsions " << std::endl;
// 	       } else { 
// 		  l.print();
// 		  std::string file_name = "atoms-by-torsion-";
// 		  file_name += link_type;
// 		  file_name += "-";
// 		  file_name += new_residue_type;
// 		  file_name += ".tab";
// 		  l.write(file_name);
// 		  l.read(file_name);
// 		  try { 
// 		     CResidue *r = l.make_residue(p.first);
// 		     if (r) {
// 			CMMDBManager *mol = coot::util::create_mmdbmanager_from_residue(r);
// 			mol->WritePDBASCII("test-out.pdb");
// 		     }
// 		  }
// 		  catch (const std::runtime_error &rte) {
// 		     std::cout << rte.what() << std::endl;
// 		  }
// 	       }
	    }
	 }
      }
      delete mol;
   }
   return 0;
} 


int main(int argc, char **argv) {

   if (argc > 2) {

      if (std::string(argv[1]) == "test") {

	 // use the link-by-torsion reference files to make a new residue (i.e. to test
	 // the function)

	 if (argc > 3) {
	    std::string link_type = argv[2];
	    std::string new_residue_type = argv[3];
	    
	    std::string file_name = "link-by-torsion-to-pyranose-core-" + link_type + ".tab";
	    if (coot::file_exists(file_name)) {
	       link_by_torsion_t l(file_name);
	       l.new_residue_type = new_residue_type;

	       // Get a base residue
	       CMMDBManager *mol = new CMMDBManager;
	       std::string pdb_file_name =
		  "pdb-templates/pyranose-pyranose-via-" + link_type + ".pdb";
	       
	       mol->ReadPDBASCII(pdb_file_name.c_str());
	       CResidue *base_residue_p = coot::util::get_first_residue(mol);
	       CResidue *r = l.make_residue(base_residue_p);
	       if (r) {
		  CMMDBManager *mol = coot::util::create_mmdbmanager_from_residue(r);
		  mol->WritePDBASCII("new-residue.pdb");
	       }
	    }
	 }

      } else {

	 // make the link-by-torsion reference files

	 std::string file_name = argv[1];
	 std::string link_type = argv[2]; // e.g. "ALPHA1-6";
	 std::string new_residue_type = "MAN";
      
	 CMMDBManager *mol = new CMMDBManager;
	 int status = mol->ReadPDBASCII(file_name.c_str());
	 if (status != Error_NoError) {
	    std::cout << "ERROR:: on reading " << file_name << std::endl;
	 } else {
	    std::pair<CResidue *, CResidue *> p = get_residue_pair(mol);
	    if (p.first && p.second) {

	       // generate link torsions, write link torsions
	       link_by_torsion_base_t link_to_core = get_names_for_link_type(link_type);
	       if (! link_to_core.filled()) {
		  std::cout << "ERROR:: " << link_type << std::endl;
	       } else {
		  link_by_torsion_t l(link_to_core, p.first, p.second);
		  if (l.filled()) {
		     std::string file_name = "link-by-torsion-to-pyranose-core-" + link_type + ".tab";
		     l.write(file_name);
		  } 
	       } 
	    }
	 }
      }
   }
   return 0;
}
   

 
