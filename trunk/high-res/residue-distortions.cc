/* high-res/residue-distortions.cc
 * 
 * Copyright 2004 by Paul Emsley, The University of York
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

#include "coot-atom-graph.hh"

 
// sc_nodes: a vector of sidechain nodes, starting at the C-beta.
//
// This removes symmetry-effects from the atom positions, then goes on
// to find distortion for each residue type.  The Ca is added to
// side_chain_pos (it was not in sc_nodes).


// OK, another go at that, were we put the residue type on the inside.
//
// i_node_ca indexes the atom_info vector: 
void
coot::atom_graph::score_all_side_chain_types(int i_node_ca,
					     const std::vector<coot::node_info> &sc_nodes,
					     coot::sequence_assignment::side_chain_score_t *scs) const {

   double d = -1;
   std::vector<clipper::RTop_orth> transformations;
   std::vector<clipper::Coord_orth> side_chain_pos(sc_nodes.size()+1);

   if (atom_info[i_node_ca].size() > 0) { 
      side_chain_pos[0] = atom_info[i_node_ca][0].atom.pos;

      for (int i=0; i<sc_nodes.size(); i++) {
	 side_chain_pos[i+1] = atoms[sc_nodes[i].index].pos;
	 if (sc_nodes[i].symm_trans_needed_flag) {
	    transformations.push_back(sc_nodes[i].rtop);
	    side_chain_pos[i+1] = get_transformed_atom(side_chain_pos[i+1], transformations);
	 }
      }

      int max_resno_in_chain = 100; // FIXME
      for (int residue_type=0; residue_type<20; residue_type++) { 

	 d = distortion_score_side_chain(residue_type, side_chain_pos);
	 scs->add_score(atom_info[i_node_ca][0].chain_number,
			chain_id(atom_info[i_node_ca][0].chain_number),
			atom_info[i_node_ca][0].residue_number,
			max_resno_in_chain,
			residue_type, d);
      }
   }

}

// 
//
// Actually, we don't want to do this (remove/apply symmetry) for
// every residue type.  We should do it once (here) - and do the
// residue type looping inside (here).
//
// This is old, redundant code.
//
// Return -1 on (indexing) failure.
// 
double
coot::atom_graph::distortion_score_side_chain_old(int i_node_ca,
					      const std::string &residue_type,
					      const std::vector<coot::node_info> &sc_nodes) const {

   double d = -1;
   std::vector<clipper::RTop_orth> transformations;
   std::vector<clipper::Coord_orth> side_chain_pos(sc_nodes.size()+1);

   if (atom_info[i_node_ca].size() > 0) { 
      side_chain_pos[0] = atom_info[i_node_ca][0].atom.pos;

      for (int i=0; i<sc_nodes.size(); i++) {
	 side_chain_pos[i+1] = atoms[sc_nodes[i].index].pos;
	 if (sc_nodes[i].symm_trans_needed_flag) {
	    transformations.push_back(sc_nodes[i].rtop);
	    side_chain_pos[i+1] = get_transformed_atom(side_chain_pos[i+1], transformations);
	 }
      }

      // Doesn't compile any more...
      // d = distortion_score_side_chain(residue_type, side_chain_pos);
   }
   return d;
}

double
//residue type is a  coot::sequence_assignment::side_chain_name_index 
coot::atom_graph::distortion_score_side_chain(int residue_type,
					      const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (residue_type == coot::sequence_assignment::SER) {
      d = distortion_score_ser(c);
   }
   if (residue_type == coot::sequence_assignment::ALA) {
      d = distortion_score_ala(c);
   }
   if (residue_type == coot::sequence_assignment::CYS) {
      d = distortion_score_cys(c);
   }
   if (residue_type == coot::sequence_assignment::LEU) {
      d = distortion_score_leu(c);
   }
   if (residue_type == coot::sequence_assignment::ASP) {
      d = distortion_score_asp(c);
   }
   if (residue_type == coot::sequence_assignment::ASN) {
      d = distortion_score_asn(c);
   }
   if (residue_type == coot::sequence_assignment::VAL) {
      d = distortion_score_val(c);
   }
   if (residue_type == coot::sequence_assignment::ILE) {
      d = distortion_score_ile(c);
   }
   if (residue_type == coot::sequence_assignment::THR) {
      d = distortion_score_thr(c);
   }
   if (residue_type == coot::sequence_assignment::GLU) {
      d = distortion_score_glu(c);
   }
   if (residue_type == coot::sequence_assignment::GLN) {
      d = distortion_score_gln(c);
   }
   if (residue_type == coot::sequence_assignment::LYS) {
      d = distortion_score_lys(c);
   }
   if (residue_type == coot::sequence_assignment::PHE) {
      d = distortion_score_phe(c);
   }
   if (residue_type == coot::sequence_assignment::TYR) {
      d = distortion_score_tyr(c);
   }
   if (residue_type == coot::sequence_assignment::PRO) {
      d = distortion_score_pro(c);
   }

   return d;

}

double
coot::atom_graph::distortion_score_ser(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 3) {
      std::cout << "ERROR in distortion_score_ser: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]);
      double len2 = clipper::Coord_orth::length(c[1], c[2]);
      double len3 = clipper::Coord_orth::length(c[0], c[2]);

      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.417);
      d += squared(len3 - 2.43);
   }

   return d;
}

double
coot::atom_graph::distortion_score_ala(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 2) {
      std::cout << "ERROR in distortion_score_ala: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]);
      d  = squared(len1 - 1.521);
   }
   return d;
}

double
coot::atom_graph::distortion_score_cys(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 3) {
      std::cout << "ERROR in distortion_score_cys: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]);
      double len2 = clipper::Coord_orth::length(c[1], c[2]);
      double len3 = clipper::Coord_orth::length(c[0], c[2]);

      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.808);
      d += squared(len3 - 2.81);
   }

   return d;
}



double
coot::atom_graph::distortion_score_leu(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 5) {
      std::cout << "ERROR in distortion_score_leu: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG  angle
      double len4 = clipper::Coord_orth::length(c[2], c[3]); // CG  GD1 bond
      double len5 = clipper::Coord_orth::length(c[2], c[4]); // CG  GD2 bond
      double len6 = clipper::Coord_orth::length(c[1], c[4]); // CB  GD2 angle
      double len7 = clipper::Coord_orth::length(c[1], c[3]); // CB  GD2 angle
      

      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.53);
      d += squared(len3 - 2.6);

      d += squared(len4 - 1.52);
      d += squared(len5 - 1.52);

      d += squared(len6 - 2.51);
      d += squared(len7 - 2.51);

   }

   return d;
}

double
coot::atom_graph::distortion_score_asp(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 5) {
      std::cout << "ERROR in distortion_score_asp: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG  angle
      double len4 = clipper::Coord_orth::length(c[2], c[3]); // CG  OD1 bond
      double len5 = clipper::Coord_orth::length(c[2], c[4]); // CG  OD2 bond
      double len6 = clipper::Coord_orth::length(c[1], c[4]); // CB  OD2 bond
      double len7 = clipper::Coord_orth::length(c[1], c[3]); // CB  OD2 bond
      

      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.52);
      d += squared(len3 - 2.53);

      d += squared(len4 - 1.25);
      d += squared(len5 - 1.25);

      d += squared(len6 - 2.4);
      d += squared(len7 - 2.4);

   }

   return d;
}


double
coot::atom_graph::distortion_score_asn(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 5) {
      std::cout << "ERROR in distortion_score_asn: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG  angle
      double len4 = clipper::Coord_orth::length(c[2], c[3]); // CG  OD1 bond
      double len5 = clipper::Coord_orth::length(c[2], c[4]); // CG  OD2 bond
      double len6 = clipper::Coord_orth::length(c[1], c[4]); // CB  OD2 bond
      double len7 = clipper::Coord_orth::length(c[1], c[3]); // CB  OD2 bond
      

      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.52);
      d += squared(len3 - 2.53);

      d += squared(len4 - 1.28); // average of CG->OD1 CG->ND2
      d += squared(len5 - 1.28);

      d += squared(len6 - 2.4);
      d += squared(len7 - 2.4);

   }

   return d;
}


double
coot::atom_graph::distortion_score_val(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 4) {
      std::cout << "ERROR in distortion_score_val: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG1 bond
      double len3 = clipper::Coord_orth::length(c[1], c[3]); // CB  CG2 bond
      double len4 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG1 angle
      double len5 = clipper::Coord_orth::length(c[0], c[3]); // CA  CG2 angle
      

      d  = squared(len1 - 1.54);
      d += squared(len2 - 1.521);
      d += squared(len3 - 1.521);

      d += squared(len4 - 2.52); 
      d += squared(len5 - 2.52);

   }

   return d;
}

double
coot::atom_graph::distortion_score_ile(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 5) {
      std::cout << "ERROR in distortion_score_val: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG1 bond
      double len3 = clipper::Coord_orth::length(c[1], c[3]); // CB  CG2 bond
      double len4 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG1 angle
      double len5 = clipper::Coord_orth::length(c[0], c[3]); // CA  CG2 angle

      double len6 = clipper::Coord_orth::length(c[2], c[4]); // CG1 CD1 bond
      double len7 = clipper::Coord_orth::length(c[1], c[4]); // CB  CD1 angle
      

      d  = squared(len1 - 1.54);
      d += squared(len2 - 1.521);
      d += squared(len3 - 1.521);

      d += squared(len4 - 2.52); 
      d += squared(len5 - 2.52);

      d += squared(len6 - 1.513);
      d += squared(len7 - 2.58);

   }

   return d;
}




double
coot::atom_graph::distortion_score_thr(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 4) {
      std::cout << "ERROR in distortion_score_thr: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG1 bond
      double len3 = clipper::Coord_orth::length(c[1], c[3]); // CB  CG2 bond
      double len4 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG1 angle
      double len5 = clipper::Coord_orth::length(c[0], c[3]); // CA  CG2 angle
      

      d  = squared(len1 - 1.54);
      d += squared(len2 - 1.476); // average of CB->OG1/CG2
      d += squared(len3 - 1.476);

      d += squared(len4 - 2.47); // average again 
      d += squared(len5 - 2.47);

   }

   return d;
}


double
coot::atom_graph::distortion_score_glu(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 6) {
      std::cout << "ERROR in distortion_score_glu: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG  angle
      double len4 = clipper::Coord_orth::length(c[2], c[3]); // CG  GD  bond
      double len5 = clipper::Coord_orth::length(c[1], c[3]); // CB  CD  angle
      double len6 = clipper::Coord_orth::length(c[3], c[4]); // CD  OE1 bond
      double len7 = clipper::Coord_orth::length(c[3], c[5]); // CD  OE2 bond

      double len8 = clipper::Coord_orth::length(c[2], c[4]); // CG  OE1 angle
      double len9 = clipper::Coord_orth::length(c[2], c[5]); // CG  OE2 angle
      
      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.52);
      d += squared(len3 - 2.53);

      d += squared(len4 - 1.516);
      d += squared(len5 - 2.53);

      d += squared(len6 - 1.25);
      d += squared(len7 - 1.25);

      d += squared(len8 - 2.4);
      d += squared(len9 - 2.4);
   }

   return d;
}

double
coot::atom_graph::distortion_score_gln(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 6) {
      std::cout << "ERROR in distortion_score_gln: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG  angle
      double len4 = clipper::Coord_orth::length(c[2], c[3]); // CG  GD  bond
      double len5 = clipper::Coord_orth::length(c[1], c[3]); // CB  CD  angle
      double len6 = clipper::Coord_orth::length(c[3], c[4]); // CD  OE1 bond
      double len7 = clipper::Coord_orth::length(c[3], c[5]); // CD  OE2 bond

      double len8 = clipper::Coord_orth::length(c[2], c[4]); // CG  OE1 angle
      double len9 = clipper::Coord_orth::length(c[2], c[5]); // CG  OE2 angle
      
      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.52);
      d += squared(len3 - 2.53);

      d += squared(len4 - 1.516);
      d += squared(len5 - 2.53);

      d += squared(len6 - 1.25); // FIXME, needs averaging
      d += squared(len7 - 1.25);

      d += squared(len8 - 2.4);  // FIXME, needs averaging
      d += squared(len9 - 2.4);
   }

   return d;
}


double
coot::atom_graph::distortion_score_lys(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 6) {
      std::cout << "ERROR in distortion_score_lys: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[2], c[3]); // CG  GD  bond
      double len4 = clipper::Coord_orth::length(c[3], c[4]); // CD  CE  bond
      double len5 = clipper::Coord_orth::length(c[4], c[5]); // CE  NZ  bond

      double len6 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG  angle
      double len7 = clipper::Coord_orth::length(c[1], c[3]); // CB  CD  angle
      double len8 = clipper::Coord_orth::length(c[2], c[4]); // CG  CE  angle
      double len9 = clipper::Coord_orth::length(c[3], c[5]); // CD  NZ  angle
      
      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.52);
      d += squared(len3 - 1.52);

      d += squared(len4 - 1.52);
      d += squared(len5 - 1.49);

      d += squared(len6 - 2.56);
      d += squared(len7 - 2.51);

      d += squared(len8 - 2.51);
      d += squared(len9 - 2.49);
   }

   return d;
}

double
coot::atom_graph::distortion_score_pro(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 5) {
      std::cout << "ERROR in distortion_score_pro: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[2], c[3]); // CG  GD  bond
      double len4 = clipper::Coord_orth::length(c[3], c[4]); // CD   N  bond

      double len5 = clipper::Coord_orth::length(c[0], c[2]); // CA  CG  angle
      double len6 = clipper::Coord_orth::length(c[1], c[3]); // CB  CD  angle

      
      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.52);
      d += squared(len3 - 1.52);

      d += squared(len4 - 2.52);
      d += squared(len5 - 2.49);
   }
   return d;
}


double
coot::atom_graph::distortion_score_phe(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 8) {
      std::cout << "ERROR in distortion_score_phe: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[2], c[3]); // CG  GD  bond
      double len4 = clipper::Coord_orth::length(c[2], c[4]); //   bond
      double len5 = clipper::Coord_orth::length(c[3], c[5]); //   bond
      double len6 = clipper::Coord_orth::length(c[4], c[6]); //   bond
      double len7 = clipper::Coord_orth::length(c[5], c[7]); //   bond
      double len8 = clipper::Coord_orth::length(c[6], c[7]); //   bond

      double len9 =  clipper::Coord_orth::length(c[0], c[2]); // 
      double len10 = clipper::Coord_orth::length(c[1], c[4]); // 
      double len11 = clipper::Coord_orth::length(c[1], c[3]); // 
      double len12 = clipper::Coord_orth::length(c[2], c[6]); // 
      double len13 = clipper::Coord_orth::length(c[2], c[5]); // 
      double len14 = clipper::Coord_orth::length(c[6], c[5]); //

      double torsion1 = clipper::Coord_orth::torsion(c[2],c[4],c[6],c[7]);
      double torsion2 = clipper::Coord_orth::torsion(c[2],c[3],c[5],c[7]);

      double abs_cv1 = abs_chiral_vol(c[2], c[1], c[3], c[4]);

      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.52);
      d += squared(len3 - 1.52);

      d += squared(len4 - 2.52);
      d += squared(len5 - 2.49);
   }
   return d;
}

double
coot::atom_graph::distortion_score_tyr(const std::vector <clipper::Coord_orth> &c) const {

   double d = -1;
   if (c.size() != 9) {
      std::cout << "ERROR in distortion_score_tyr: c length: " << c.size() << "\n";
   } else { 
      double len1 = clipper::Coord_orth::length(c[0], c[1]); // CA  CB  bond
      double len2 = clipper::Coord_orth::length(c[1], c[2]); // CB  CG  bond
      double len3 = clipper::Coord_orth::length(c[2], c[3]); // CG  GD  bond
      double len4 = clipper::Coord_orth::length(c[2], c[4]); //   bond
      double len5 = clipper::Coord_orth::length(c[3], c[5]); //   bond
      double len6 = clipper::Coord_orth::length(c[4], c[6]); //   bond
      double len7 = clipper::Coord_orth::length(c[5], c[7]); //   bond
      double len8 = clipper::Coord_orth::length(c[6], c[7]); //   bond

      double len9 =  clipper::Coord_orth::length(c[0], c[2]); // 
      double len10 = clipper::Coord_orth::length(c[1], c[4]); // 
      double len11 = clipper::Coord_orth::length(c[1], c[3]); // 
      double len12 = clipper::Coord_orth::length(c[2], c[6]); // 
      double len13 = clipper::Coord_orth::length(c[2], c[5]); // 
      double len14 = clipper::Coord_orth::length(c[6], c[5]); //

      double torsion1 = clipper::Coord_orth::torsion(c[2],c[4],c[6],c[7]);
      double torsion2 = clipper::Coord_orth::torsion(c[2],c[3],c[5],c[7]);

      double abs_cv1 = abs_chiral_vol(c[2], c[1], c[3], c[4]);

      d  = squared(len1 - 1.53);
      d += squared(len2 - 1.52);
      d += squared(len3 - 1.52);

      d += squared(len4 - 2.52);
      d += squared(len5 - 2.49);
   }
   return d;
}

