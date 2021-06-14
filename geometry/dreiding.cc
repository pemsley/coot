/* geometry/dreiding.cc
 * 
 * Copyright 2013 by Medical Research Council
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

#include "protein-geometry.hh"


// can thow a std::runtime_error
double
coot::protein_geometry::dreiding_torsion_energy(const std::string &comp_id,
						int imol_enc,
						mmdb::Atom *atom_0,
						mmdb::Atom *atom_1,
						mmdb::Atom *atom_2,
						mmdb::Atom *atom_3) const {

   double d = 0;

   if (!(atom_0 && atom_1 && atom_2 && atom_3)) {
      throw std::runtime_error("Null atom in dreiding_torsion_energy");
   } else {
      // happy path
      int indx = get_monomer_restraints_index(comp_id, imol_enc, true);
      if (indx != -1) {
	 const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[indx].second;
	 std::vector<std::string> name(4);
	 std::vector<std::string> energy_type(4);
	 std::vector<int> sp_hybrid(4);
	 name[0] = atom_0->name;
	 name[1] = atom_1->name;
	 name[2] = atom_2->name;
	 name[3] = atom_3->name;
	 for (unsigned int i=0; i<4; i++) { 
	    energy_type[i] = restraints.type_energy(name[i]);
	    std::map<std::string, energy_lib_atom>::const_iterator atom_map_it =
	       energy_lib.atom_map.find(energy_type[i]);
	    if (atom_map_it == energy_lib.atom_map.end()) {
	       std::string m = "No energy lib for type ";
	       m += energy_type[i];
	       throw std::runtime_error(m);
	    }
	    sp_hybrid[i] = atom_map_it->second.sp_hybridisation;
	 }
	 clipper::Coord_orth p0(atom_0->x, atom_0->y, atom_0->z);
	 clipper::Coord_orth p1(atom_1->x, atom_1->y, atom_1->z);
	 clipper::Coord_orth p2(atom_2->x, atom_2->y, atom_2->z);
	 clipper::Coord_orth p3(atom_3->x, atom_3->y, atom_3->z);
	 // double phi = clipper::Coord_orth::torsion(p0,p1,p2,p3);
	 // d = dreiding_torsion_energy(phi, sp_hybrid[1], sp_hybrid[2], "dummy", false, false);
      }
   }
   return d;
}

double
coot::protein_geometry::dreiding_torsion_energy(const std::string &comp_id,
						int imol_enc,
						const atom_quad &quad) const {

   return dreiding_torsion_energy(comp_id, imol_enc,
				  quad.atom_1, quad.atom_2, quad.atom_3, quad.atom_4);
}

// double
// coot::protein_geometry::dreiding_torsion_energy(double phi, dreiding_torsion_energy_t dr) const {
//    return dreiding_torsion_energy(phi, dr.Vjk, dr.phi0_jk, dr.n_jk);
// } 

coot::protein_geometry::dreiding_torsion_energy_t
coot::protein_geometry::dreiding_torsion_energy_params(const std::string &comp_id,
						       int imol_enc,
						       const atom_quad &quad) const {

   dreiding_torsion_energy_t dr(0,0,0);

   if (! quad.filled_p()) {
      throw std::runtime_error("Null atom in dreiding_torsion_energy params");
   } else {
      // happy path
      int indx = get_monomer_restraints_index(comp_id, imol_enc, true);
      if (indx != -1) {
	 const coot::dictionary_residue_restraints_t &restraints = dict_res_restraints[indx].second;
	 std::vector<std::string> name(4);
	 std::vector<std::string> energy_type(4);
	 std::vector<int> sp_hybrid(4);
	 name[0] = quad.atom_1->name;
	 name[1] = quad.atom_2->name;
	 name[2] = quad.atom_3->name;
	 name[3] = quad.atom_4->name;
	 for (unsigned int i=0; i<4; i++) { 
	    energy_type[i] = restraints.type_energy(name[i]);
	    std::map<std::string, energy_lib_atom>::const_iterator atom_map_it =
	       energy_lib.atom_map.find(energy_type[i]);
	    if (atom_map_it == energy_lib.atom_map.end()) {
	       std::string m = "No energy lib for type ";
	       m += energy_type[i];
	       throw std::runtime_error(m);
	    }
	    sp_hybrid[i] = atom_map_it->second.sp_hybridisation;
	 }
      }
   }
   return dr;
} 



coot::protein_geometry::dreiding_torsion_energy_t 
coot::protein_geometry::dreiding_torsion_energy_params(double phi, int sp_a1, int sp_a2,
						       const std::string &bond_order,
						       bool at_1_deloc_or_arom,
						       bool at_2_deloc_or_arom) const {
   dreiding_torsion_energy_t dr(0,0,0);

   // single bond, sp3s
   if (sp_a1 == 3 && sp_a2 == 3) {
      dr.Vjk = 2.0;
      dr.n_jk = 3;
      dr.phi0_jk = M_PI; // radians!
   }
   // single bond sp2 and sp3
   if  ((sp_a1 == 2 && sp_a2 == 3) || (sp_a1 == 3 && sp_a2 == 2)) {
      dr.Vjk = 1.0;
      dr.n_jk = 6;
      dr.phi0_jk = 0.0;
   }
   // double bond sp2-sp2
   if (sp_a1 == 2 && sp_a2 == 2) {
      dr.Vjk = 45.0;
      dr.n_jk = 2;
      dr.phi0_jk = M_PI;
   }
   // resonance bond (deloc or arom)

   // single bond between two sp2 atoms - or deloc atoms
   //   exception: exocyclic single bond involving 2 aromatics

   if (bond_order == "aromatic" || bond_order == "deloc") {
      
   }

   // g: anything with sp1
   if (sp_a1 == 1 || sp_a2 == 1)
      dr.Vjk = 0;

   return dr;
}


// phi is in radians!
// 
double
coot::protein_geometry::dreiding_torsion_energy(double phi,
						double Vjk, double phi0_jk, double n_jk) const {
   double E = 0.5*Vjk*(1.0-cos(n_jk*(phi-phi0_jk)));
   return E;
}
