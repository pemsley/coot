/* geometry/protein-geometry.cc
 * 
 * Copyright 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
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

#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>  // needed for sort? Yes.
#include <stdexcept>  // Thow execption.

#include "utils/win-compat.hh"
#include "mini-mol/atom-quads.hh"
#include "geometry/protein-geometry.hh"
#include "utils/coot-utils.hh"

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define DATADIR "C:/coot/share"
#define PKGDATADIR DATADIR
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

#include "clipper/core/clipper_util.h"

#include "compat/coot-sysdep.h"

#include "lbg-graph.hh"


// static
unsigned int
coot::chem_link::make_hash_code(const std::string &comp_id_1, const std::string &comp_id_2, const std::string &group_1, const std::string &group_2) {

   unsigned int hash = 0;
   unsigned int hash_c1 = 0;
   unsigned int hash_c2 = 0;
   unsigned int hash_g1 = 0;
   unsigned int hash_g2 = 0;

   std::string local_group_1 = group_1;
   std::string local_group_2 = group_2;

   if (local_group_1 == "L-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "L-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "P-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "P-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "M-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "M-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "D-pyranose")   local_group_1 = "pyranose";
   if (local_group_2 == "D-pyranose")   local_group_2 = "pyranose";
   if (local_group_1 == "D-SACCHARIDE") local_group_1 = "pyranose";  // CCD annotation for MAN, etc
   if (local_group_2 == "D-SACCHARIDE") local_group_2 = "pyranose";
   if (local_group_1 == "SACCHARIDE")   local_group_1 = "pyranose";  // CCD annotation for FUC
   if (local_group_2 == "SACCHARIDE")   local_group_2 = "pyranose";

   if (local_group_1 == "RNA") local_group_1 = "DNA/RNA";
   if (local_group_2 == "RNA") local_group_2 = "DNA/RNA";

   for (unsigned int i = 0; i < comp_id_1.length(); i++) {
     unsigned int chr = comp_id_1[i];
     hash_c1  = ((hash_c1 << 5) - hash_c1) + chr;
     hash_c1 |= 0; // Convert to 32bit integer
   }
   for (unsigned int i = 0; i < comp_id_2.length(); i++) {
     unsigned int chr = comp_id_2[i];
     hash_c2  = ((hash_c2 << 5) - hash_c2) + chr;
     hash_c2 |= 0;
   }
   for (unsigned int i = 0; i < local_group_1.length(); i++) {
     unsigned int chr = local_group_1[i];
     hash_g1  = ((hash_g1 << 5) - hash_g1) + chr;
     hash_g1 |= 0;
   }
   for (unsigned int i = 0; i < local_group_2.length(); i++) {
     unsigned int chr = local_group_2[i];
     hash_g2  = ((hash_g2 << 5) - hash_g2) + chr;
     hash_g2 |= 0;
   }

   // hash = hash_c1 + 2 * hash_c2 + 4 * hash_g1 + hash_g2;

   hash = hash_g1 + 8 * hash_g2;

   // hash |= 0;

   if (false)
      std::cout << "debug:: in make_hash_code \"" << group_1 << "\" \"" << group_2 << "\" -> \""
		<< local_group_1 << "\" \"" << local_group_2 << "\"" << " return hash " << hash
                << std::endl;

   return hash;
}

std::pair<bool, bool>
coot::chem_link::matches_comp_ids_and_groups_hashed(const std::string &comp_id_1,
                                                    const std::string &group_1,
                                                    const std::string &comp_id_2,
                                                    const std::string &group_2) const {

   bool match = false; // initially
   bool order_switch = false;

   std::string local_group_1 = group_1;
   std::string local_group_2 = group_2;

   // chem_links specify "peptide" or "pyranose", but comp_groups are "L-peptide"/"D-pyranose".
   // So allow them to match.
   // 201201013 (Friday) allow M-peptides to match too.
   if (local_group_1 == "L-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "L-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "M-peptide")  local_group_1 = "peptide";
   if (local_group_2 == "M-peptide")  local_group_2 = "peptide";
   if (local_group_1 == "D-pyranose") local_group_1 = "pyranose";
   if (local_group_2 == "D-pyranose") local_group_2 = "pyranose";
   if (local_group_1 == "D-SACCHARIDE") local_group_1 = "pyranose";  // CCD annotation for MAN, etc
   if (local_group_1 == "SACCHARIDE") local_group_1 = "pyranose";    // CCD annotation for FUC
   if (local_group_2 == "D-SACCHARIDE") local_group_2 = "pyranose";
   if (local_group_2 == "SACCHARIDE") local_group_2 = "pyranose";

   if (local_group_1 == "RNA") local_group_1 = "DNA/RNA";
   if (local_group_2 == "RNA") local_group_2 = "DNA/RNA";

   unsigned int hash_test = make_hash_code(comp_id_1, comp_id_2, local_group_1, local_group_2);

   if (hash_test == hash_code) {
      match = true;
      order_switch = false;
   }

   return std::pair<bool, bool>(match, order_switch);
}


std::pair<bool, bool>
coot::chem_link::matches_comp_ids_and_groups_hashed(unsigned int hash_test_forwards,
                                                    unsigned int hash_test_backwards) const {

   bool match = false; // initially
   bool order_switch = false;

   if (hash_test_forwards == hash_code)
      match = true;

   if (hash_test_backwards == hash_code) {
      match = true;
      order_switch = true;
   }

   return std::pair<bool, bool>(match, order_switch);
}



//  return (match,order_switch) pair
bool
coot::chem_link::matches_comp_ids_and_groups(const std::string &comp_id_1,
					     const std::string &group_1,
					     const std::string &comp_id_2,
					     const std::string &group_2) const {

   bool debug = false;

   if (debug) {

#if 1
      std::cout << "   ------ DEBUG:: in matches_comp_ids_and_groups() "
		<< id << " chem_link_name: " << chem_link_name << ": input comp_ids "
		<< comp_id_1 << " and " << comp_id_2 << "          vs ref link-comp_id-1 :"
		<< chem_link_comp_id_1 << ": ref link-comp_id-2 :"
		<< chem_link_comp_id_2 << ":" << std::endl; 
      std::cout << "           for chem_link_comp_name " << chem_link_name << ": input groups "
		<< group_1 << " and " << group_2 << "             vs ref link-group-1 :"
		<< chem_link_group_comp_1 << ": ref link-group-2 :"
		<< chem_link_group_comp_2 << ":" << std::endl;
#endif

      std::cout << "   chem-link self group    1: \"" << chem_link_group_comp_1 << "\" self group 2   \"" << chem_link_group_comp_2 << "\"" << std::endl;
      std::cout << "   chem-link input group   1: \"" << group_1                << "\" input group 2  \"" << group_2                << "\"" << std::endl;
      std::cout << "   chem-link self comp_id  1: \"" << chem_link_comp_id_1    << "\" self comp_id 2 \"" << chem_link_comp_id_2    << "\"" << std::endl;
      std::cout << "   chem-link input comp_id 1: \"" << comp_id_1              << "\" input group 2  \"" << comp_id_2              << "\"" << std::endl;

   }

   unsigned int hash_test = make_hash_code(comp_id_1, comp_id_2, group_1, group_2);

   if (debug)
      std::cout << "   checking:: hash comparison this: " << hash_code << " test: "
		<< hash_test << " " << *this << "\n";

   bool match = false; // initially
   bool order_switch = false;

   std::string local_group_1 = group_1; 
   std::string local_group_2 = group_2;

   // chem_links specify "peptide" or "pyranose", but comp_groups are "L-peptide"/"D-pyranose".
   // So allow them to match.
   // 201201013 (Friday) allow M-peptides to match too.
   // If you are thinkg of adding/changing this, change the one in make_hash_code() too.
   if (local_group_1 == "L-peptide")    local_group_1 = "peptide";
   if (local_group_2 == "L-peptide")    local_group_2 = "peptide";
   if (local_group_1 == "P-peptide")    local_group_1 = "peptide";
   if (local_group_2 == "P-peptide")    local_group_2 = "peptide";
   if (local_group_1 == "M-peptide")    local_group_1 = "peptide";
   if (local_group_2 == "M-peptide")    local_group_2 = "peptide";
   if (local_group_1 == "D-pyranose")   local_group_1 = "pyranose";
   if (local_group_2 == "D-pyranose")   local_group_2 = "pyranose";
   if (local_group_1 == "D-SACCHARIDE") local_group_1 = "pyranose";  // CCD annotation for MAN, etc
   if (local_group_2 == "D-SACCHARIDE") local_group_2 = "pyranose";
   if (local_group_1 ==   "SACCHARIDE") local_group_1 = "pyranose";  // CCD annotation for FUC
   if (local_group_2 ==   "SACCHARIDE") local_group_2 = "pyranose";

   // if (local_group_1 == "RNA") local_group_1 = "DNA/RNA";
   // if (local_group_2 == "RNA") local_group_2 = "DNA/RNA";

   std::string self_group_1 = chem_link_group_comp_1;
   std::string self_group_2 = chem_link_group_comp_2;

   if (self_group_1 == "RNA") self_group_1 = "RNA/DNA"; // to match the input type ! (which was converted to make the hash)
   if (self_group_2 == "RNA") self_group_2 = "RNA/DNA";

   if (debug)
      std::cout << "     sighx... check match: group_1: dict \""
		<< self_group_1 << "\" vs model \"" << local_group_1 << "\" -and- dict group_2: \""
		<< self_group_2 << "\" vs model \"" << local_group_2 << "\"\n";

   if (local_group_2 == "SACCHARIDE") local_group_2 = "pyranose";

   if (debug) {
      std::cout << "   ------ DEBUG:: in matches_comp_ids_and_groups() "
		<< id << " chem_link_name " << chem_link_name << ": input comp_ids "
		<< comp_id_1 << " and " << comp_id_2 << " vs ref link-comp_id-1 :"
		<< chem_link_comp_id_1 << ": ref link-comp_id-2:"
		<< chem_link_comp_id_2 << ":" << std::endl;
      std::cout << "         for chem_link_comp_name " << chem_link_name << ": input groups "
		<< local_group_1 << " and " << local_group_2 << " vs ref link-group-1 :"
		<< chem_link_group_comp_1 << ": ref link-group-2 :"
		<< chem_link_group_comp_2 << ":" << std::endl;
   }

   if (((self_group_1 == "") || (self_group_1 == local_group_1)) &&
       ((self_group_2 == "") || (self_group_2 == local_group_2))) {
      if (((chem_link_comp_id_1 == "") || (chem_link_comp_id_1 == comp_id_1)) &&
	  ((chem_link_comp_id_2 == "") || (chem_link_comp_id_2 == comp_id_2))) {
	 match = true;
      }
   }

   if (((chem_link_group_comp_1 == "DNA/RNA") && (local_group_1 == "RNA") && 
	(chem_link_group_comp_2 == "DNA/RNA") && (local_group_2 == "RNA")))
      match = true;

#if 0 // 20220917-PE I don't know what this was supposed to do -
      // it looks suspicious.
   if (((chem_link_group_comp_1 == "DNA/RNA") && (local_group_1 == "DNA") && 
	(chem_link_group_comp_1 == "DNA/RNA") && (local_group_2 == "DNA")))
      match = true;
#endif

   return match;
}

std::ostream& coot::operator<<(std::ostream &s, coot::chem_link lnk) {

   std::string p1 = lnk.chem_link_comp_id_1;
   std::string p2 = lnk.chem_link_comp_id_2;
   std::string id = lnk.id;
   int l1 = id.length();
   if (l1 < 5) {
      int d = 5 -l1;
      id.append(d, ' ');
   }
   s << "[chem_link: id: " << id
     << " [comp_id1: \"" << p1 << "\" group_1: \"" << lnk.chem_link_group_comp_1
     << "\" mod_1: \"" << lnk.chem_link_mod_id_1 << "\"] to "
     << " [comp_id2: \"" << p2 << "\" group_2: \"" << lnk.chem_link_group_comp_2
     << "\" mod_2: \"" << lnk.chem_link_mod_id_2 << "\"] " << lnk.chem_link_name << "]";
   return s; 
}

std::ostream& coot::operator<<(std::ostream &s, coot::list_chem_mod mod) {

   s << "[list_chem_mod: id: " << mod.id << " " 
     << "name: " << mod.name 
     << " comp_id :" << mod.comp_id 
     << ": group_id: " << mod.group_id
     << "]";
   return s;
} 



double
coot::dict_chiral_restraint_t::assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds,
							   const std::vector <dict_angle_restraint_t> &angles) {

   double vol = -1;
   double a = -1, b = -1, c = -1;
   double alpha = -1, beta = -1, gamma = -1;
   std::string mmdb_centre_atom =  atom_id_mmdb_expand(local_atom_id_centre);
   std::string mmdb_local_atom_id_1 = atom_id_mmdb_expand(local_atom_id_1);
   std::string mmdb_local_atom_id_2 = atom_id_mmdb_expand(local_atom_id_2);
   std::string mmdb_local_atom_id_3 = atom_id_mmdb_expand(local_atom_id_3);
   
   // local_atom_id_centre to local_atom_id_1 bond length
   for (unsigned int i=0; i<bonds.size(); i++) {
      try { 
	 if (bonds[i].atom_id_1_4c() == mmdb_centre_atom) {
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_1)) { 
	       a = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_2)) { 
	       b = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_3)) { 
	       c = bonds[i].value_dist();
	    }
	 }
	 if (bonds[i].atom_id_2_4c() == atom_id_mmdb_expand(local_atom_id_centre)) {
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_1)) { 
	       a = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_2)) { 
	       b = bonds[i].value_dist();
	    }
	    if (bonds[i].atom_id_1_4c() == atom_id_mmdb_expand(local_atom_id_3)) { 
	       c = bonds[i].value_dist();
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 // do nothing, it's not really an error if the dictionary
	 // doesn't have target geometry (the bonding description came
	 // from a Chemical Component Dictionary entry for example).
      } 
   }

   for (unsigned int i=0; i<angles.size(); i++) {
      if (angles[i].atom_id_2_4c() == mmdb_centre_atom) {
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_2 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_3) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_2 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_3))  {
	    alpha = clipper::Util::d2rad(angles[i].angle());
	 }
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_3) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_3))  {
	    beta = clipper::Util::d2rad(angles[i].angle());
	 }
	 if ((angles[i].atom_id_1_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_3_4c() == mmdb_local_atom_id_2) ||
	     (angles[i].atom_id_3_4c() == mmdb_local_atom_id_1 &&
	      angles[i].atom_id_1_4c() == mmdb_local_atom_id_2))  {
	    gamma = clipper::Util::d2rad(angles[i].angle());
	 }
      }
   }

   
   if (a > 0 && b > 0 && c > 0) {
//       std::cout << "DEBUG:: found all distances in chiral restraint" << std::endl;
      if (alpha > 0 && beta > 0 && gamma > 0) {
// 	 std::cout << "DEBUG:: found all angles in chiral restraint" << std::endl;
	 vol = assign_chiral_volume_target_internal(a, b, c, alpha, beta, gamma);
      } else {
// 	 std::cout << "DEBUG:: failed to find all angles in chiral restraint"
// 		   << alpha << " " << beta << " " << gamma << std::endl;
      }
   } else {
//       std::cout << "DEBUG:: failed to find all distances in chiral restraint"
// 		<< a << " " << b << " " << c << std::endl;
   }
   return vol;
}

double
coot::dict_link_chiral_restraint_t::assign_chiral_volume_target(const std::vector <dict_bond_restraint_t> &bonds_1,
								const std::vector <dict_angle_restraint_t> &angles_1,
								const std::vector <dict_bond_restraint_t> &bonds_2,
								const std::vector <dict_angle_restraint_t> &angles_2,
								const std::vector <dict_link_bond_restraint_t> &link_bonds,
								const std::vector <dict_link_angle_restraint_t> &link_angles) {

   double d = 0;

   return d;
} 


// angles in radians.
double
coot::dict_chiral_restraint_t::assign_chiral_volume_target_internal(double a, double b, double c,
								    double alpha, double beta, double gamma) {

   // set target_volume
   // from: abc ( 1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma) + 2(cos(alpha) + cos(beta) + cos(gamma)))^0.5

   double cos_alpha = cos(alpha);
   double cos_beta  = cos(beta);
   double cos_gamma = cos(gamma);
   
   double cos_2_alpha = cos_alpha * cos_alpha;
   double cos_2_beta  = cos_beta  * cos_beta;
   double cos_2_gamma = cos_gamma * cos_gamma;


   double param = 1-cos_2_alpha-cos_2_beta-cos_2_gamma + 2*cos_alpha*cos_beta*cos_gamma;
   if (false) {
      std::cout << "assign_chiral_volume_target_internal() input a, b, c, alpha, beta, gamma " << a << " "
		<< b << " " << c << " "
		<< clipper::Util::rad2d(alpha) << " "
		<< clipper::Util::rad2d(beta) << " "
		<< clipper::Util::rad2d(gamma) << " " << std::endl;

      double a_bit = 1-cos_2_alpha-cos_2_beta-cos_2_gamma;
      double b_bit = 2 * cos_alpha * cos_beta * cos_gamma;
      double c_bit = a_bit + b_bit;
      std::cout << "    bits: " << a_bit << " " << b_bit << " " << c_bit << std::endl;
      std::cout << "    param: " << param << std::endl;
   }
   if (param < 0) param = 0;
   target_volume_ = volume_sign * a*b*c*sqrt(param);

   // volume_sigma_ = 0.2;  // seems reasonable give target voluemes of about 2.6
   // volume_sigma_ = 0.3; // test
   // volume_sigma_ = 0.18; // meh, 0.3 seemed to give "too many" chiral errors when
                         // refining large numbers of atoms in poor/EM maps.
   // 20210228-PE OK, molprobity is giving "too many" CB outliers - let's tighten this up a bit.
   // 20210315-PE OK, a bit more (was 0.13)
   volume_sigma_ = 0.1;

   if (false)
      std::cout << "DEBUG:: assign_chiral_volume_target_internal() target_volume chiral: "
		<< target_volume_ << std::endl;

   // give the appropriate dictionary, this can return a target_volume_ of nan
   return target_volume_;
}


// return "" on failure.
// no order switch is considered.
// 
std::string
coot::protein_geometry::find_glycosidic_linkage_type(mmdb::Residue *first, mmdb::Residue *second) const {

   // Fixup needed for PDBv3

   bool debug = false;
   double critical_dist = 2.4; // A, less than that and Coot should
			       // try to make the bond.
                               // 20170505: changed to 2.4, was 3.0.
                               // Needed to stop the beta1-6 linked MAN
                               // on a BMA linkinking to the NAG to which
                               // the BMA is bonded (A605-602 in 3whe).
                              
   mmdb::PPAtom res_selection_1 = NULL;
   mmdb::PPAtom res_selection_2 = NULL;
   int i_no_res_atoms_1;
   int i_no_res_atoms_2;
   double d;
   std::vector<coot::glycosidic_distance> close;
   
   first->GetAtomTable( res_selection_1, i_no_res_atoms_1);
   second->GetAtomTable(res_selection_2, i_no_res_atoms_2);

   for (int i1=0; i1<i_no_res_atoms_1; i1++) { 
      clipper::Coord_orth a1(res_selection_1[i1]->x,
			     res_selection_1[i1]->y,
			     res_selection_1[i1]->z);
      for (int i2=0; i2<i_no_res_atoms_2; i2++) {
	 clipper::Coord_orth a2(res_selection_2[i2]->x,
				res_selection_2[i2]->y,
				res_selection_2[i2]->z);
	 d = (a1-a2).lengthsq();
	 if (d < critical_dist*critical_dist) {
	    close.push_back(coot::glycosidic_distance(res_selection_1[i1],
						      res_selection_2[i2],
						      sqrt(d)));
	 }
      }
   }

   std::sort(close.begin(), close.end());

   // if you consider to uncomment this to debug a repulsion instead
   // the forming of a glycosidic bond, consider the residue numbering
   // of the residues involved: the "residue 1" should have the O4 and
   // the "residue 2" (+1 residue number) should have the C1.
   // 
   if (debug) {
      std::cout << "DEBUG:: find_glycosidic_linkage_type() "
		<< "number of sorted distances: "
		<< close.size() << std::endl;
      for (unsigned int i=0; i<close.size(); i++) {
	 std::cout << "#### glyco close: " << close[i].distance << "  "
		   << close[i].at1->GetChainID() << " " 
		   << close[i].at1->GetSeqNum() << " " 
		   << close[i].at1->GetAtomName() << " " 
		   << " to "
		   << close[i].at2->GetChainID() << " " 
		   << close[i].at2->GetSeqNum() << " " 
		   << close[i].at2->GetAtomName() << " " 
		   << std::endl;
      }
   }

   std::string link_type("");

   // glyco_chiral constructor can throw an exception
   try { 
   
      float smallest_link_dist = 99999.9;
      for (unsigned int i=0; i<close.size(); i++) {
	 std::string name_1(close[i].at1->name);
	 std::string name_2(close[i].at2->name);


	 // First test the NAG-ASN link (that order - as per dictionary)
	 //
	 if (name_1 == " C1 ")
	    if (name_2 == " ND2")
	       if (close[i].distance < smallest_link_dist) {
		  smallest_link_dist = close[i].distance;
		  link_type = "NAG-ASN";
	       }
      
	 if (name_1 == " O4 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-4");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-4";
		  }
	       }
      
	 if (name_1 == " O2 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-2");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-2";
		  }
	       }

	 if (name_1 == " O3 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-3");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-3";
		  }
	       }
      

	 // The BETA2-3 link should never happen :-)
	 // There are no biosynthetic pathways to make an BETA2-3 link for a SIA.
	 // (SIA BETA2-3 would be axial if it existed)
	 //
	 if (name_1 == " C2 " )
	    if (name_2 == " O3 ")
	       if (std::string(close[i].at1->GetResName()) == "SIA") {
		  if (close[i].distance < smallest_link_dist) {
		     coot::atom_quad glyco_chiral_quad(first, second, "ALPHA2-3");
		     std::cout << "   glyco_chiral ALPHA2-3 "
			       << close[i].at1->GetResName() << " "
			       << close[i].at2->GetResName() << " "
			       << glyco_chiral_quad.chiral_volume() << std::endl;
		     if (glyco_chiral_quad.chiral_volume() > 0.0) {
			smallest_link_dist = close[i].distance;
			link_type = "ALPHA2-3";
		     }
		  }
	       }

	 // 20180111 Add ALPHA2-6 links for SIA
	 if (name_1 == " C2 " )
	    if (name_2 == " O6 ")
	       if (std::string(close[i].at1->GetResName()) == "SIA") {
		  if (close[i].distance < smallest_link_dist) {
		     coot::atom_quad glyco_chiral_quad(first, second, "ALPHA2-6");
		     std::cout << "   glyco_chiral ALPHA2-6 "
			       << close[i].at1->GetResName() << " "
			       << close[i].at2->GetResName() << " "
			       << glyco_chiral_quad.chiral_volume() << std::endl;
		     if (glyco_chiral_quad.chiral_volume() > 0.0) {
			smallest_link_dist = close[i].distance;
			link_type = "ALPHA2-6";
		     }
		  }
	       }

	 if (name_1 == " O6 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "BETA1-6");
		  if (glyco_chiral_quad.chiral_volume() > 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "BETA1-6";
		  }
	       }
	       
	 if (name_1 == " O2 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-2");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-2";
		  }
	       }
	       
	 if (name_1 == " O3 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-3");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-3";
		  }
	       }

	 if (name_1 == " C2 " )
	    if (name_2 == " O3 ")
	       if (std::string(close[i].at1->GetResName()) == "SIA") { 
		  if (close[i].distance < smallest_link_dist) {
		     coot::atom_quad glyco_chiral_quad(first, second, "ALPHA2-3");
		     std::cout << "   glyco_chiral ALPHA2-3 "
			       << close[i].at1->GetResName() << " "
			       << close[i].at2->GetResName() << " "
			       << glyco_chiral_quad.chiral_volume() << std::endl;
		     if (glyco_chiral_quad.chiral_volume() < 0.0) { 
			smallest_link_dist = close[i].distance;
			link_type = "ALPHA2-3";
		     }
		  }
	       }
      
	 if (name_1 == " O4 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-4");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-4";
		  }
	       }
      
	 if (name_1 == " O6 " )
	    if (name_2 == " C1 ")
	       if (close[i].distance < smallest_link_dist) {
		  coot::atom_quad glyco_chiral_quad(first, second, "ALPHA1-6");
		  if (glyco_chiral_quad.chiral_volume() < 0.0) { 
		     smallest_link_dist = close[i].distance;
		     link_type = "ALPHA1-6";
		  }
	       }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING::" << rte.what() << std::endl;
   }

   if (debug)
      std::cout << "   debug:: find_glycosidic_linkage_type() for "
		<< first->GetChainID() << " " << first->GetSeqNum() << " " << first->GetInsCode()
		<< first->GetResName() << ","
		<< second->GetChainID() << " " << second->GetSeqNum() << " " << second->GetInsCode()
		<< second->GetResName() 
		<< " returns \"" << link_type << "\""
		<< std::endl;
   
   return link_type;
}

std::string
coot::protein_geometry::find_glycosidic_linkage_type(mmdb::Residue *first,
						     mmdb::Residue *second,
						     mmdb::Manager *mol) const {
   
   // First check that the residues are LINK - or are sequential.
   // If so, then we can find the link as above.
   
   std::string link_type;
   bool are_linked = false;
   bool are_sequential = false;

   // Test for sequential/tandem
   //
   std::string chain_id_1 =  first->GetChainID();
   std::string chain_id_2 = second->GetChainID();
   int resno_1 =  first->GetSeqNum();
   int resno_2 = second->GetSeqNum();
   if (chain_id_1 == chain_id_2) {
      if (resno_1 == (resno_2+1)) are_sequential = true;
      if (resno_2 == (resno_1+1)) are_sequential = true;
   } 

   if (! are_sequential) {
      // Test for link in molecule
      std::string ins_code_1 =  first->GetInsCode();
      std::string ins_code_2 = second->GetInsCode();
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) { 
	 mmdb::LinkContainer *links = model_p->GetLinks();
	 int n_links = model_p->GetNumberOfLinks();
	 for (int ilink=1; ilink<=n_links; ilink++) { 
	    mmdb::Link *link = model_p->GetLink(ilink);
	    if (link) {
	       are_linked = are_linked_in_order(first, second, link);
	       if (! are_linked)
		  are_linked = are_linked_in_order(first, second, link);
	       if (are_linked)
		  break;
	    }
	 }
      }
   } 


   if (are_linked || are_sequential)
      link_type = find_glycosidic_linkage_type(first, second);

   return link_type;
}


bool
coot::protein_geometry::are_linked_in_order(mmdb::Residue *first,
					    mmdb::Residue *second,
					    mmdb::Link *link) const {

   bool linked = false;
   // c.f link_atoms() - but we can't use that because that's a
   // coot-utils function and geometry is more primitive.
   std::string link_chain_id_1(link->chainID1);
   std::string link_chain_id_2(link->chainID2);
   std::string chain_id_1 =  first->GetChainID();
   std::string chain_id_2 = second->GetChainID();
   int resno_1 =  first->GetSeqNum();
   int resno_2 = second->GetSeqNum();
   if (link_chain_id_1 == chain_id_1) { 
      if (link_chain_id_2 == chain_id_2) {
	 int link_reso_1 = link->seqNum1;
	 int link_reso_2 = link->seqNum2;
	 if (link_reso_1 == resno_1) { 
	    if (link_reso_2 == resno_2) {
	       std::string link_ins_code_1 = link->insCode1;
	       std::string link_ins_code_2 = link->insCode2;
	       std::string ins_code_1 =  first->GetInsCode();
	       std::string ins_code_2 = second->GetInsCode();
	       if (link_ins_code_1 == ins_code_1) { 
		  if (link_ins_code_2 == ins_code_2) {
		     linked = true;
		  }
	       }
	    }
	 }
      }
   }

   return linked;
} 
					    

std::pair<std::string, bool>
coot::protein_geometry::find_glycosidic_linkage_type_with_order_switch(mmdb::Residue *first, mmdb::Residue *second) const {

   std::pair<std::string, bool> r("", false);

   std::string l = find_glycosidic_linkage_type(first, second);

   if (l == "") { 
      l = find_glycosidic_linkage_type(second,first);
      if (l != "") {
	 r.first = l;
	 r.second = true;
      } 
   } else {
      r.first = l;
      r.second = false;
   } 
   return r;
} 




void
coot::protein_geometry::assign_link_chiral_volume_targets() {

   for (unsigned int idict=0; idict<dict_link_res_restraints.size(); idict++) {
      if (dict_link_res_restraints[idict].has_unassigned_chiral_volumes()) {
	 dict_link_res_restraints[idict].assign_link_chiral_volume_targets();
      }
   }
}


void
coot::protein_geometry::print_chem_links() const {

   // simple vector
   // for (unsigned int i_chem_link=0; i_chem_link<chem_link_vec.size(); i_chem_link++)
   // std::cout<< i_chem_link << " " << chem_link_vec[i_chem_link] << "\n";

   std::map<unsigned int, std::vector<chem_link> >::const_iterator it;

   for (it=chem_link_map.begin(); it!=chem_link_map.end(); ++it) {
      const std::vector<chem_link> &v = it->second;
      std::vector<chem_link>::const_iterator itv;
      for (itv=v.begin(); itv!=v.end(); ++itv) {
         const chem_link &cl = *itv;
         std::cout << "     " << it->first << " " << cl << "\n";
      }
   }

}



bool
coot::protein_geometry::linkable_residue_types_p(const std::string &this_res_type,
						 const std::string &env_res_type) {

   std::pair<short int, coot::dictionary_residue_restraints_t> r1 = get_monomer_restraints(this_res_type, IMOL_ENC_ANY);
   std::pair<short int, coot::dictionary_residue_restraints_t> r2 = get_monomer_restraints(env_res_type, IMOL_ENC_ANY);

   bool r = 0;
   if (r1.first) {
      if (r1.second.residue_info.group != "non-polymer")
	 r = 1;
   }
   if (r2.first) {
      if (r2.second.residue_info.group != "non-polymer")
	 r = 1;
   }
   return r;
} 

