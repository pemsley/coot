/* coot-utils/glyco-torsions.cc
 * 
 * Copyright 2011, 2012 by The University of Oxford
 * Copyright 2014 by Medical Research Council
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

// The problem with carbohydrate-building - as it stands.
//
// The set of reference structures has pyranose-BETA1-4-pyranose.  And
// we paste a MAN onto the extension part of that to give us initial
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
#include "utils/coot-utils.hh"
#include "glyco-torsions.hh"

std::ostream&
coot::operator<<(std::ostream &o, const atom_by_torsion_t &abt) {

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


// ext_residue_p is the residue being added.
//
clipper::Coord_orth
coot::atom_by_torsion_t::pos(mmdb::Residue *base_residue_p, mmdb::Residue *ext_residue_p) const {

   mmdb::Atom *at_1 = NULL;
   mmdb::Atom *at_2 = NULL;
   mmdb::Atom *at_3 = NULL;

   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   base_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   if (false) {
      std::cout << "positioning atom " << atom_name << std::endl;
      for (int iat=0; iat<n_residue_atoms; iat++) { 
	 mmdb::Atom *at = residue_atoms[iat];
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
      unsigned int n_missing = 0;
      if (! at_1) n_missing++;
      if (! at_2) n_missing++;
      if (! at_3) n_missing++;
      std::string m = "missing atom";
      if (n_missing>1)
	 m += "s";
      m += " in atom_by_torsion_t::pos() when positioning ";
      m += atom_name;
      m += " : ";
      if (! at_1) m += " at_1 " + prior_atom_1.second;
      if (! at_2) m += " at_2 " + prior_atom_2.second;
      if (! at_3) m += " at_3 " + prior_atom_3.second;
      m += " of ";
      m += util::int_to_string(n_residue_atoms);
      m += " base atoms";
      throw std::runtime_error(m);
   }
}

coot::link_by_torsion_t::link_by_torsion_t(const std::string &link_type,
					   const std::string &new_residue_comp_id_in) {

   new_residue_type = new_residue_comp_id_in;
   new_res_no = 1; // FIXME

   std::string fn = link_type_to_file_name(link_type, new_residue_comp_id_in);
   read(fn);

   // now handle the decorations
   // 
   std::string decor_file_name = comp_id_to_decoration_file_name(new_residue_comp_id_in);
   if (! file_exists(decor_file_name)) {
      std::cout << "No file " << decor_file_name << std::endl;
   } else { 
      coot::link_by_torsion_t decor(decor_file_name);
      if (! decor.filled()) {
	 std::cout << "Decorations not filled from " << decor_file_name
		   << std::endl;
      } else {
	 add(decor);
      }
   }
}

std::string
coot::link_by_torsion_t::link_type_to_file_name(const std::string &link_type,
						const std::string &new_res_comp_id) const {

   // try to get new_res_comp_id-specific template, if not, fall back to generic
   //
   std::string p = package_data_dir();
   std::string f = "link-by-torsion-to-" + new_res_comp_id + "-core-" + link_type + ".tab";
   std::string ff = util::append_dir_file(p,f);

   std::cout << "......... checking for " << ff << std::endl;

   if (file_exists(ff)) {
      return ff;
   } else {
      f = "link-by-torsion-to-pyranose-core-" + link_type + ".tab";
      ff = util::append_dir_file(p,f);
      std::cout << "..that failed - trying  " << ff << std::endl;
   }
   return ff;
} 

std::string
coot::link_by_torsion_t::link_type_to_file_name(const std::string &link_type) const {

   std::string p = package_data_dir();
   // std::string d = util::append_dir_dir(p, "pdb-templates"); not yet
   std::string f = "link-by-torsion-to-pyranose-core-" + link_type + ".tab";
   std::string ff = util::append_dir_file(p,f);

   return ff;
} 

std::string
coot::link_by_torsion_t::comp_id_to_decoration_file_name(const std::string &comp_id) const {

   std::string p = package_data_dir();
   std::string f = new_residue_type + "-decorations.tab";
   std::string ff = util::append_dir_file(p,f);

   return ff;
} 
   
mmdb::Residue *
coot::link_by_torsion_t::make_residue(mmdb::Residue *base_residue_p) const {

   mmdb::Residue *r = NULL;
   if (geom_atom_torsions.size()) {
      r = new mmdb::Residue;
      r->SetResName(new_residue_type.c_str());
      r->seqNum = new_res_no;
      for (unsigned int i=0; i<geom_atom_torsions.size(); i++) {
	 const atom_by_torsion_t &gat = geom_atom_torsions[i];
	 // std::cout << "in make_residue() i: " << i << " " << gat << std::endl;
	 std::string f = gat.filled_atom_name(); // FIXME PDBv3 the function call is not needed
	 clipper::Coord_orth p = gat.pos(base_residue_p, r);
	 mmdb::Atom *atom = new mmdb::Atom(r); // does an add atom
	 atom->SetAtomName(f.c_str());
	 atom->SetElementName(gat.element.c_str());
	 atom->SetCoordinates(p.x(), p.y(), p.z(), 1.0, b_factor);
	 atom->Het = 1;
	 if (false)
	    std::cout << "   " << gat.atom_name << " " << p.format()  << std::endl;
      }
   } 
   return r;
} 


coot::atom_by_torsion_t::atom_by_torsion_t(const atom_by_torsion_base_t &names,
					   mmdb::Residue *residue_1_p,  // reference/lower
					   mmdb::Residue *residue_2_p   // extension residue
					   ) {
   
   if (false)
      std::cout << "get "
		<< names.prior_atom_1.first << " " << names.prior_atom_1.second << " "
		<< names.prior_atom_2.first << " " << names.prior_atom_2.second << " "
		<< names.prior_atom_3.first << " " << names.prior_atom_3.second << " "
		<< std::endl;
   mmdb::PPAtom residue_atoms_1 = 0;
   mmdb::PPAtom residue_atoms_2 = 0;
   int n_residue_atoms_1;
   int n_residue_atoms_2;
   residue_1_p->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   residue_2_p->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

   mmdb::Atom *p_new = residue_2_p->GetAtom(names.atom_name.c_str());
   if (p_new) { 
      mmdb::Atom *p_1 = NULL;
      mmdb::Atom *p_2 = NULL;
      mmdb::Atom *p_3 = NULL;

      // I could just GetAtom() here.
      for (int iat=0; iat<n_residue_atoms_1; iat++) { 
	 mmdb::Atom *at = residue_atoms_1[iat];
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
      for (int iat=0; iat<n_residue_atoms_2; iat++) { 
	 mmdb::Atom *at = residue_atoms_2[iat];
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
	    std::cout << "... " << coot::atom_spec_t(p_new) << " has pos " << pos_n.format()
		      << std::endl;
	    std::cout << "... " << coot::atom_spec_t(p_1)   << " has pos " << pos_1.format()
		      << std::endl;
	 }
	 double bl = sqrt((pos_n - pos_1).lengthsq());
	 double a = q.angle_2();
	 double t = q.torsion();
	 // bleugh, FIXME with an init() function.
	 *this = atom_by_torsion_t(names, bl, a, t);
      }
   }
}

coot::link_by_torsion_base_t coot::asn_pyranose_link_to_core() {

   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;

   ats.push_back(atom_by_torsion_base_t("C1", "C", BS(true,  "ND2"), BS(true,  "CG" ), BS(true,  "CB" )));
   ats.push_back(atom_by_torsion_base_t("C2", "C", BS(false, "C1" ), BS(true,  "ND2"), BS(true,  "CG" )));
   ats.push_back(atom_by_torsion_base_t("C3", "C", BS(false, "C2" ), BS(false, "C1" ), BS(true,  "ND2")));
   ats.push_back(atom_by_torsion_base_t("C4", "C", BS(false, "C3" ), BS(false, "C2" ), BS(false, "C1" )));
   ats.push_back(atom_by_torsion_base_t("C5", "C", BS(false, "C4" ), BS(false, "C3" ), BS(false, "C2" )));
   ats.push_back(atom_by_torsion_base_t("O5", "O", BS(false, "C5" ), BS(false, "C4" ), BS(false, "C3" )));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
} 

coot::link_by_torsion_base_t coot::ser_pyranose_link_to_core() {

   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;

   ats.push_back(atom_by_torsion_base_t("C1", "C", BS(true,  "OG"),  BS(true,  "CB" ), BS(true,  "CA" )));
   ats.push_back(atom_by_torsion_base_t("C2", "C", BS(false, "C1" ), BS(true,  "OG" ), BS(true,  "CB" )));
   ats.push_back(atom_by_torsion_base_t("C3", "C", BS(false, "C2" ), BS(false, "C1" ), BS(true,  "OG" )));
   ats.push_back(atom_by_torsion_base_t("C4", "C", BS(false, "C3" ), BS(false, "C2" ), BS(false, "C1" )));
   ats.push_back(atom_by_torsion_base_t("C5", "C", BS(false, "C4" ), BS(false, "C3" ), BS(false, "C2" )));
   ats.push_back(atom_by_torsion_base_t("O5", "O", BS(false, "C5" ), BS(false, "C4" ), BS(false, "C3" )));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}


coot::link_by_torsion_base_t coot::pyranose_link_1_6_to_core() {
   
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

coot::link_by_torsion_base_t coot::pyranose_link_1_4_to_core() {
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

coot::link_by_torsion_base_t coot::pyranose_link_1_2_to_core() {
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

coot::link_by_torsion_base_t coot::pyranose_link_1_3_to_core() {
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

coot::link_by_torsion_base_t coot::pyranose_link_2_3_to_core() {

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

coot::link_by_torsion_base_t coot::mannose_decorations() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   // I don't like this depending on O6 of the previous residue
   ats.push_back(atom_by_torsion_base_t("O2", "O", BS(false, "C2"), BS(false, "C1"), BS(true,  "O6")));
   ats.push_back(atom_by_torsion_base_t("O3", "O", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("O4", "O", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("C6", "C", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   ats.push_back(atom_by_torsion_base_t("O6", "O", BS(false, "C6"), BS(false, "C5"), BS(false, "C4")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

// FIXME
coot::link_by_torsion_base_t coot::glucose_decorations() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   // I don't like this depending on O6 of the previous residue
   ats.push_back(atom_by_torsion_base_t("O2", "O", BS(false, "C2"), BS(false, "C1"), BS(true,  "O6")));
   ats.push_back(atom_by_torsion_base_t("O3", "O", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("O4", "O", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("C6", "C", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   ats.push_back(atom_by_torsion_base_t("O6", "O", BS(false, "C6"), BS(false, "C5"), BS(false, "C4")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}


coot::link_by_torsion_base_t coot::fucose_decorations() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   ats.push_back(atom_by_torsion_base_t("O2", "O", BS(false, "C2"), BS(false, "C1"), BS(true,  "O3")));
   ats.push_back(atom_by_torsion_base_t("O3", "O", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("O4", "O", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("C6", "C", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   ats.push_back(atom_by_torsion_base_t("O6", "O", BS(false, "C6"), BS(false, "C5"), BS(false, "C4")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

coot::link_by_torsion_base_t coot::galactose_decorations() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   // O4 of previous residue, from an 1-4 link.
   ats.push_back(atom_by_torsion_base_t("O2", "O", BS(false, "C2"), BS(false, "C1"), BS(true,  "O4")));
   ats.push_back(atom_by_torsion_base_t("O3", "O", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("O4", "O", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("C6", "C", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   ats.push_back(atom_by_torsion_base_t("O6", "O", BS(false, "C6"), BS(false, "C5"), BS(false, "C4")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

// FIXME
coot::link_by_torsion_base_t coot::NAG_decorations() {
   link_by_torsion_base_t l;
   std::vector<atom_by_torsion_base_t> ats;
   ats.push_back(atom_by_torsion_base_t("N2", "N", BS(false, "C2"), BS(false, "C1"), BS(false,  "O5")));
   ats.push_back(atom_by_torsion_base_t("C7", "C", BS(false, "N2"), BS(false, "C2"), BS(false,  "C1")));
   ats.push_back(atom_by_torsion_base_t("C8", "C", BS(false, "C7"), BS(false, "N2"), BS(false,  "C2")));
   ats.push_back(atom_by_torsion_base_t("O7", "O", BS(false, "C7"), BS(false, "N2"), BS(false,  "C2")));
   
   ats.push_back(atom_by_torsion_base_t("O3", "O", BS(false, "C3"), BS(false, "C2"), BS(false, "C1")));
   ats.push_back(atom_by_torsion_base_t("O4", "O", BS(false, "C4"), BS(false, "C3"), BS(false, "C2")));
   ats.push_back(atom_by_torsion_base_t("C6", "C", BS(false, "C5"), BS(false, "C4"), BS(false, "C3")));
   ats.push_back(atom_by_torsion_base_t("O6", "O", BS(false, "C6"), BS(false, "C5"), BS(false, "C4")));
   for (unsigned int i=0; i<ats.size(); i++) l.add(ats[i]);
   return l;
}

coot::link_by_torsion_base_t
coot::get_decorations(const std::string &new_comp_id) {
   
   if (new_comp_id == "MAN") return mannose_decorations();
   if (new_comp_id == "BMA") return mannose_decorations();
   if (new_comp_id == "GLC") return glucose_decorations();
   if (new_comp_id == "FUC") return fucose_decorations(); // was glucose_decorations()!
   if (new_comp_id == "FUL") return fucose_decorations();
   if (new_comp_id == "GAL") return galactose_decorations();
   if (new_comp_id == "NAG") return NAG_decorations();
   link_by_torsion_base_t empty;
   return empty;
} 


// When the link type is xxx1-Y, then we can simply add up the torsions.
// When the link type is xxx2-3, then we have an O3 (which would otherwise be "decoration")
//   that is part of the link atoms.  In that case, one of the O3s should be omitted.
//   The add() function should check if this atom exists already and if so not add the new
//   one.
// 
// coot::link_by_torsion_base_t coot::add_link_by_torsions(const coot::link_by_torsion_base_t &v1,
// 							const coot::link_by_torsion_base_t &v2) {

//    link_by_torsion_base_t r = v1;
//    for (unsigned int i=0; i<v2.atom_torsions.size(); i++)
//       r.add(v2.atom_torsions[i]);
//    return r;
// } 

coot::link_by_torsion_base_t
coot::get_names_for_link_type(const std::string &link_type) {

   // std::cout << "here in get_names_for_link_type() " << link_type << std::endl;
   
   link_by_torsion_base_t r;
   if (link_type == "ALPHA1-6") r = pyranose_link_1_6_to_core();
   if (link_type == "ALPHA1-2") r = pyranose_link_1_2_to_core();
   if (link_type == "ALPHA1-3") r = pyranose_link_1_3_to_core();
   if (link_type == "ALPHA2-3") r = pyranose_link_2_3_to_core();
   if (link_type == "BETA1-2")  r = pyranose_link_1_2_to_core();
   if (link_type == "BETA1-3")  r = pyranose_link_1_3_to_core();
   if (link_type == "BETA1-4")  r = pyranose_link_1_4_to_core();
   if (link_type == "BETA1-6")  r = pyranose_link_1_6_to_core();
   if (link_type == "NAG-ASN")  r = asn_pyranose_link_to_core();
   if (link_type == "NAG-SER")  r = ser_pyranose_link_to_core();
   return r;
}
   
// static
std::pair<mmdb::Residue *, mmdb::Residue *>
coot::link_by_torsion_t::get_residue_pair(mmdb::Manager *mol) {

   std::pair<mmdb::Residue *, mmdb::Residue *> r(NULL, NULL);
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::Residue *residue_p;
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

// void coot::write_tree(mmdb::Residue *ref_res_p, mmdb::Residue *ext_res_p) {

//    // Just a test function
   
//    link_by_torsion_base_t a16 = add_link_by_torsions(pyranose_link_1_6_to_core(), mannose_decorations());

//    for (unsigned int i=0; i<a16.atom_torsions.size(); i++) {
//       atom_by_torsion_t abt = get_atom_by_torsion(a16.atom_torsions[i], ref_res_p, ext_res_p);
//       if (abt.filled()) {
// 	 std::cout << abt << std::endl;
//       }
//    }
// } 

void
coot::link_by_torsion_t::init(mmdb::Residue *ref_res_p, mmdb::Residue *ext_res_p) {

   b_factor = 31; // default
   if (0) 
      std::cout << "in link_by_torsion_t::init() have " << atom_torsions.size()
		<< "  atom torsions" << std::endl;
   for (unsigned int i=0; i<atom_torsions.size(); i++) {
      atom_by_torsion_t abt(atom_torsions[i], ref_res_p, ext_res_p);
      if (! abt.filled()) {
	 std::cout << "Missing atom! " << abt << std::endl;
      } else {
	 add(abt);
      } 
   }
}

void
coot::link_by_torsion_t::set_temperature_factor(float b) {
   b_factor = b;
} 


void
coot::link_by_torsion_t::print() const {

   for (unsigned int i=0; i<geom_atom_torsions.size(); i++)
      std::cout << "   " << std::setw(2) << i << " " << geom_atom_torsions[i] << std::endl;
} 
 
void
coot::link_by_torsion_t::write(const std::string &file_name) const {

   std::ofstream f(file_name.c_str());
   if (f)
      for (unsigned int i=0; i<geom_atom_torsions.size(); i++)
	 f << "  "  << " " << geom_atom_torsions[i] << "\n";
}

// read what is written by write() function
void
coot::link_by_torsion_t::read(const std::string &file_name) {

   if (! file_exists(file_name)) {
      std::cout << "ERROR:: file not found " << file_name << std::endl;
      return;
   } else {
      std::cout << "reading " << file_name << std::endl;
   } 

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
   // std::cout << "finished read() with " << geom_atom_torsions.size() << " torsions" << std::endl;
   if (geom_atom_torsions.size() == 0) 
      std::cout << "After read()ing, we have " << geom_atom_torsions.size()
		<< " atom torsions" << std::endl;
}


