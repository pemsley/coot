/* ideal/pepflip.cc
 * 
 * Copyright 2002, 2003 The University of York
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

#include <string>
#include "coot-utils/coot-coord-utils.hh"
#include "pepflip.hh"

// mmdb-style interface
// 
// Given a mol and a residue number and chain of the first residue
// in the peptide (i.e. the residue with the C and O atoms), flip
// the C and O atoms of this peptide and the N of the next one (if
// exists) round a line joining this Ca to the next one.
//
// Return status is 0 if the flip did not happen (because, for
// example, either or both of the Ca's could not be found).
//
// mol is manipulated.
int
coot::pepflip(CMMDBManager *mol,
	      const std::string &chain_id,
	      int resno, 
	      const std::string &ins_code, 
	      const std::string &altconf) {

   int status = pepflip_standard(mol, chain_id, resno, ins_code, altconf);
   if (status)
      return status;
   else
      return pepflip_internal_to_residue(mol, chain_id, resno, ins_code, altconf);
}

   
int
coot::pepflip_standard(CMMDBManager *mol,
	      const std::string &chain_id,
	      int resno, 
	      const std::string &ins_code, 
	      const std::string &altconf) {


   // We need to find the 2 Ca's.  If we do not we return 0.
   //
   // Then we need to find the coordinates of the C, O of this atom
   // and the N of the next one.  Put them into a vector of
   // Coord_orths.  If size of the resulting Coord_orths vector is 3,
   // return 0.  Life is too short to mess about with indexing and
   // names due to missing atoms.

   // We then move the Ca of the resno Ca to the origin, and the other
   // atoms correspondingly.
   // 
   // The RTop around the vector...?
   //
   int status = 0; // fail initially

   std::vector<const char *> second_atoms;
   std::vector<const char *> first_atoms(2);
   first_atoms[0] = " C  ";
   first_atoms[1] = " O  ";
   second_atoms.push_back(" N  ");
   second_atoms.push_back(" H  ");

   CAtom *ca1 = NULL;
   CAtom *ca2 = NULL;

   CResidue *first_res = coot::util::get_residue(chain_id, resno, ins_code, mol);
   std::vector<CAtom *> flipping_atoms;
   if (first_res) {
      coot::residue_spec_t rs(first_res);
      CResidue *second_res = coot::util::get_following_residue(rs, mol);
      if (second_res) {
	 PCAtom *first_residue_atoms = NULL;
	 PCAtom *second_residue_atoms = NULL;
	 int n_first_residue_atoms; 
	 int n_second_residue_atoms;
	 first_res->GetAtomTable(first_residue_atoms, n_first_residue_atoms);
	 second_res->GetAtomTable(second_residue_atoms, n_second_residue_atoms);
	 for (unsigned int iat=0; iat<n_first_residue_atoms; iat++) {
	    std::string atom_name(first_residue_atoms[iat]->name);
	    std::string alt_conf_atom(first_residue_atoms[iat]->altLoc);
	    if (alt_conf_atom == altconf || alt_conf_atom == "") {
	       if (atom_name == " CA " ) {
		  ca1 = first_residue_atoms[iat];
	       }
	       for (unsigned int i=0; i<first_atoms.size(); i++) {
		  if (atom_name == first_atoms[i]) {
		     flipping_atoms.push_back(first_residue_atoms[iat]);
		  } 
	       }
	    } 
	 }
	 for (unsigned int iat=0; iat<n_second_residue_atoms; iat++) {
	    std::string atom_name(second_residue_atoms[iat]->name);
	    std::string alt_conf_atom(second_residue_atoms[iat]->altLoc);
	    if (alt_conf_atom == altconf || alt_conf_atom == "") {
	       if (atom_name == " CA " ) {
		  ca2 = second_residue_atoms[iat];
	       }
	       for (unsigned int i=0; i<second_atoms.size(); i++) {
		  if (atom_name == second_atoms[i]) {
		     flipping_atoms.push_back(second_residue_atoms[iat]);
		  } 
	       }
	    }
	 }
	 if (! ca1) {
	    std::cout << "WARNING:: No first CA atom found" << std::endl;
	 } else { 
	    if (! ca2) {
	       std::cout << "WARNING:: No second CA atom found" << std::endl;
	    } else {
	       status = 1;
	       std::vector<clipper::Coord_orth> cas(2);
	       cas[0] = clipper::Coord_orth(ca1->x, ca1->y, ca1->z);
	       cas[1] = clipper::Coord_orth(ca2->x, ca2->y, ca2->z);
	       std::vector<clipper::Coord_orth> v = 
		  flip_internal(cas, flipping_atoms);
	       for (unsigned int i=0; i<v.size(); i++) {
		  flipping_atoms[i]->x = v[i].x();
		  flipping_atoms[i]->y = v[i].y();
		  flipping_atoms[i]->z = v[i].z();
	       } 
	    } 
	 } 
      }
   }

   return status; // 1 success, 0 failure
}

   
int
coot::pepflip_internal_to_residue(CMMDBManager *mol,
				  const std::string &chain_id,
				  int resno, 
				  const std::string &ins_code, 
				  const std::string &altconf) {

   int status = 0;
   CResidue *residue_p = coot::util::get_residue(chain_id, resno, ins_code, mol);
   if (residue_p) {
      CAtom *c_at  = NULL;
      CAtom *o_at  = NULL;
      CAtom *ca_at = NULL;
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 std::string atom_name(at->name);
	 std::string atom_alt_conf(at->altLoc);
	 // PDBv3 FIXME
	 if (atom_alt_conf == altconf) { 
	    if (atom_name == " CA ")
	       ca_at = at;
	    if (atom_name == " C  ")
	       c_at = at;
	    if (atom_name == " O  ")
	       o_at = at;
	 }
	 if (c_at && o_at && ca_at) {
	    clipper::Coord_orth p1(ca_at->x, ca_at->y, ca_at->z);
	    clipper::Coord_orth p2(c_at->x,   c_at->y,  c_at->z);
	    clipper::Coord_orth p3(o_at->x,   o_at->y,  o_at->z);
	    clipper::Coord_orth p3_new = util::rotate_round_vector(p2-p1, p3, p1, M_PI);
	    o_at->x = p3_new.x();
	    o_at->y = p3_new.y();
	    o_at->z = p3_new.z();
	    status = true;
	 }
      }
   }
   return status;
} 




std::vector<clipper::Coord_orth> 
coot::flip_internal(const std::vector<clipper::Coord_orth> &ca_in,
		    const std::vector<CAtom *> &atoms) {

   std::vector<clipper::Coord_orth> atoms_orth(atoms.size()); // returned thing
   std::vector<clipper::Coord_orth> cas = ca_in;

   clipper::Coord_orth trans = cas[0];

   cas[0] -= trans;
   cas[1] -= trans;

   for (unsigned int i=0;i<atoms.size(); i++) {
      atoms_orth[i] = clipper::Coord_orth(atoms[i]->x, atoms[i]->y, atoms[i]->z);
      atoms_orth[i] -= trans;
   }

   clipper::Coord_orth ca_vec_unit(cas[1].unit());

   // Polar coordinates:

   // ca_vec_unit[0] = ( l )    ( sin omega cos phi )
   // ca_vec_unit[1] = ( m )  = ( sin omega sin phi )
   // ca_vec_unit[2] = ( n )    ( cos omega )

   double l = ca_vec_unit[0];
   double m = ca_vec_unit[1];
   double n = ca_vec_unit[2];

   double ll = l*l;
   double mm = m*m;
   double nn = n*n;
   
   // The Rotation matrix applying omega and phi and 180 around k.
   // 
   // cos k = -1,    sin k = 0:
   // 
   // ( l**2-(m**2+n**2)   2lm                 2nl              )
   // ( 2lm                m**2-(l**2+n**2)    2mn              )
   // ( 2nl                2mn                 n**2-(l**2+m**2) )
   //
   // (Amore documentation) Thanks for that pointer EJD :).

   clipper::Mat33<double> r(ll-(mm+nn),   2.0*l*m,          2.0*n*l, 
			    2.0*l*m,      mm-(ll+nn),       2.0*m*n,          
			    2.0*n*l,      2.0*m*n,          nn-(ll+mm) );

   clipper::RTop_orth rtop(r, clipper::Coord_orth(0,0,0));

   for (int i=0;i<atoms_orth.size(); i++) { 
      atoms_orth[i] = atoms_orth[i].transform(rtop);
      atoms_orth[i] += trans;
   }

   return atoms_orth;
}
