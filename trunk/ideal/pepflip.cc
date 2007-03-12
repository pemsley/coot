/* src/pepflip.cc
 * 
 * Copyright 2002, 2003 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <string>
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
// Typically, one would copy one's mol (and save it) before calling
// this.
// 
int
coot::pepflip(CMMDBManager *mol, int resno, 
	      const std::string &altconf,
	      const std::string &chain_id) {

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


   // So, let's get the Ca's:
   //
   int SelHnd_ca1 = mol->NewSelection();
   int SelHnd_ca2 = mol->NewSelection();
   PPCAtom SelAtom_ca1;
   PPCAtom SelAtom_ca2;
   int nSelAtoms_ca1;
   int nSelAtoms_ca2;
   short int found_ca = 0;
   int status = 0;
   std::string altloc_composite;
   altloc_composite = "[" ;
   altloc_composite += altconf;
   altloc_composite += "";
   altloc_composite += "]";
   

   mol->SelectAtoms(SelHnd_ca1,
		    0,
		    (char *) chain_id.c_str(),
		    resno, "*",
		    resno, "*",
		    "*", // rnames
		    " CA ",
		    "*", // elements
		    "*"
		    );
   
   mol->GetSelIndex(SelHnd_ca1, SelAtom_ca1, nSelAtoms_ca1);

   mol->SelectAtoms(SelHnd_ca2,
		    0,
		    (char *)chain_id.c_str(),
		    resno+1, "*",
		    resno+1, "*",
		    "*", // rnames
		    " CA ",
		    "*", // elements
		    "*" // altconf
		    );
   
   mol->GetSelIndex(SelHnd_ca1, SelAtom_ca1, nSelAtoms_ca1);
   mol->GetSelIndex(SelHnd_ca2, SelAtom_ca2, nSelAtoms_ca2);

   if (nSelAtoms_ca1 > 0)
      if (nSelAtoms_ca2 > 0) 
	 found_ca = 1;
      

   if (!found_ca) {
      std::cout << "WARNING: Failed pepflip: failed to find CA atoms for some reason."
		<< " Sorry\n";
      return 0;
   } else { 
      std::cout << "DEBUG:: found CAs\n";
   } 

   std::vector<clipper::Coord_orth> ca(2);
   ca[0] = clipper::Coord_orth(SelAtom_ca1[0]->x,
			       SelAtom_ca1[0]->y,
			       SelAtom_ca1[0]->z);
   ca[1] = clipper::Coord_orth(SelAtom_ca2[0]->x,
			       SelAtom_ca2[0]->y,
			       SelAtom_ca2[0]->z);

   std::cout << "   CA 0 " << ca[0].format() << std::endl;
   std::cout << "   CA 1 " << ca[1].format() << std::endl;

   mol->DeleteSelection(SelHnd_ca1);
   mol->DeleteSelection(SelHnd_ca2);

   // Now the other atoms:

   
   std::vector<char *> first_atoms(3);
   std::vector<int> selHnd_first(3);
   std::vector<PPCAtom> SelAtom_first(3);
   std::vector<int> nSelAtoms_first(3);
   PPCAtom  SelAtom_second = NULL;
   int nSelAtoms_second;
   int SelHnd_second  = mol->NewSelection();

   first_atoms[0] = "  C ";
   first_atoms[1] = "  O ";
   first_atoms[2] = "  H ";
   
   char *second_atom = "  N ";
   std::vector<CAtom *> flipped_atom(1);
   
   for (int iat=0; iat<3; iat++) { 
      SelAtom_first[iat] = NULL;
      selHnd_first[iat] = mol->NewSelection();
      mol->SelectAtoms(selHnd_first[iat],
		       0,
		       (char *)chain_id.c_str(),
		       resno, "*",
		       resno, "*",
		       "*", // rnames
		       first_atoms[iat],
		       "*", // elements
		       "*" // altLocs 
		       );
      mol->GetSelIndex(selHnd_first[iat], SelAtom_first[iat], nSelAtoms_first[iat]);
      if (nSelAtoms_first[iat] > 0) { 
	 flipped_atom[0] = SelAtom_first[iat][0];
	 std::cout << " start coords for iat: " << iat << " "
		   << flipped_atom[0]->x << " " 
		   << flipped_atom[0]->y << " " 
		   << flipped_atom[0]->z << " " 
		   << std::endl;
	 std::vector<clipper::Coord_orth> c = flip_internal(ca, flipped_atom);
	 SelAtom_first[iat][0]->x = c[0].x();
	 SelAtom_first[iat][0]->y = c[0].y();
	 SelAtom_first[iat][0]->z = c[0].z();
	 status = 1;
	 std::cout << iat << " assinging coords" << c[0].format() << std::endl;
      }
      mol->DeleteSelection(selHnd_first[iat]);
   }

   mol->SelectAtoms(SelHnd_second,
		    0,
		    (char *)chain_id.c_str(),
		    resno+1, "*",
		    resno+1, "*",
		    "*", // rnames
		    second_atom,
		    "*", // elements
		    "*" // altLocs 
		    );
   mol->GetSelIndex(SelHnd_second, SelAtom_second, nSelAtoms_second);
   if (nSelAtoms_second > 0) { 
      flipped_atom[0] = SelAtom_second[0];
      std::cout << " start coords for second " 
		<< flipped_atom[0]->x << " " 
		<< flipped_atom[0]->y << " " 
		<< flipped_atom[0]->z << " " 
		<< std::endl;
      std::vector<clipper::Coord_orth> c = flip_internal(ca, flipped_atom);
      SelAtom_second[0]->x = c[0].x();
      SelAtom_second[0]->y = c[0].y();
      SelAtom_second[0]->z = c[0].z();
      std::cout << "assinging coords for second " << c[0].format() << std::endl;
   }
   mol->DeleteSelection(SelHnd_second);

   return status; // 1 success, 0 failure
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
