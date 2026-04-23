/*
 * analysis/daca.cc
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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


#include <map>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "compat/coot-sysdep.h"
#include "utils/coot-utils.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/helix-like.hh"
#include "daca.hh"

coot::daca::box_index_t::box_index_t(const clipper::Coord_orth &pos) {
   box_width = 1.0;
   idx_x = floor(pos.x()/box_width);
   idx_y = floor(pos.y()/box_width);
   idx_z = floor(pos.z()/box_width);
}

// the reverse of the above - make a point in the middle of the box
clipper::Coord_orth
coot::daca::box_index_t::coord_orth() const {
   double x = static_cast<double>(idx_x) * box_width + 0.5 * box_width;
   double y = static_cast<double>(idx_y) * box_width + 0.5 * box_width;
   double z = static_cast<double>(idx_z) * box_width + 0.5 * box_width;
   return clipper::Coord_orth(x,y,z);
}

float
coot::daca::box_index_t::d_squared() const {
   clipper::Coord_orth pt = coord_orth();
   return (pt.x() * pt.x() + pt.y() * pt.y() + pt.z() * pt.z());
}

float
coot::daca::box_index_t::d() const {
   return sqrtf(d_squared());
}

bool
coot::daca::box_index_t::operator<(const coot::daca::box_index_t &other) const {
   if (idx_x < other.idx_x) return true;
   if (idx_x > other.idx_x) return false;
   if (idx_y < other.idx_y) return true;
   if (idx_y > other.idx_y) return false;
   if (idx_z < other.idx_z) return true;
   if (idx_z > other.idx_z) return false;
   return false;
}


float
coot::daca::gompertz_scale(const float &dist_squared) {

   float box_width = 1.0; // should match the box width of the box_index_t.
                          // Should be transfered?

   std::map<float, float>::const_iterator it = envelope_distance_map.find(dist_squared);
   if (it != envelope_distance_map.end()) {
      return it->second;
   } else {
      float x = sqrtf(dist_squared);
      const float m = 9.2 * box_width; // 9.2 not 8.0 - for better tailing off of the function
      const float a = 1.0;
      const float b = 7.0;
      const float c = 1.0;
      float g = a * exp(-b * exp(-c * (m - x)));
      envelope_distance_map[dist_squared] = g;
      return g;
   }

}


// I can't get where this should go.
#if 0
std::ostream &
coot::daca::operator<<(std::ostream &s, const coot::daca::box_index_t &bi) {
   s << "[box " << bi.x << " " << bi.y << " " << bi.z << "]";
   return s;
}
#endif

void
coot::daca::fill_reference_fragments() {

   std::string pkg_data_dir = coot::package_data_dir();
   std::string fn = util::append_dir_file(pkg_data_dir, "standard-residues.pdb");
   if (file_exists(fn)) {
      atom_selection_container_t asc = get_atom_selection(fn, false, false, false);
      if (asc.read_success) {
         // make reference fragments for each of the residues
         // std::cout << "Now do things with residues in standard_residues\n";
         int imod = 1;
         mmdb::Model *model_p = asc.mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     std::string res_name(residue_p->GetResName());
                     std::vector<std::vector<std::string> > atom_name_sets =
                        atom_names_for_fragments(res_name);

                     for(unsigned int iset=0; iset<atom_name_sets.size(); iset++) {
                        std::vector<clipper::Coord_orth> v; // this has special order
                        const std::vector<std::string> &us = atom_name_sets[iset];
                        std::vector<std::string>::const_iterator it;
                        for (it=us.begin(); it!=us.end(); ++it) {
                           const std::string &frag_atom_name = *it;
                           int n_residue_atoms;
                           mmdb::PPAtom residue_atoms = 0;
                           residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
                           for (int iat=0; iat<n_residue_atoms; iat++) {
                              mmdb::Atom *at = residue_atoms[iat];
                              std::string this_atom_name(at->GetAtomName());
                              if (this_atom_name == frag_atom_name) {
                                 clipper::Coord_orth pos = co(at);
                                 v.push_back(pos);
                                 break;
                              }
                           }
                        }
                        if (v.size() == us.size()) {
                           // add another set of coordinates to reference_fragments[res_name]
                           clipper::Coord_orth sum(0,0,0); // for calculating the centre of the fragment
                           for (unsigned int ii=0; ii<v.size(); ii++)
                              sum += v[ii];
                           double m = 1.0/static_cast<double>(v.size());
                           clipper::Coord_orth fragment_centre(sum * m);
                           for (unsigned int ii=0; ii<v.size(); ii++)
                              v[ii] -= fragment_centre;
                           reference_fragments[res_name].push_back(v);
                           if (false)
                              std::cout << " filling " << residue_p << " " << res_name << " "
                                        << v.size() << std::endl;
                        }
                     }
                  }
               }
            }
         }
      }
   } else {
      std::cout << "File not found " << fn << std::endl;
   }

   if (false) { // debugging
      std::map<std::string, std::vector<std::vector<clipper::Coord_orth> > >::const_iterator it;
      for (it =reference_fragments.begin(); it!=reference_fragments.end(); it++){
         std::cout << "Reference Residue type " << it->first << "\n";
         std::vector<std::vector<clipper::Coord_orth> >::const_iterator itvv;
         for (itvv=it->second.begin(); itvv!=it->second.end(); itvv++) {
            const std::vector<clipper::Coord_orth> &v = *itvv;
            std::vector<clipper::Coord_orth>::const_iterator itv;
            for (itv=v.begin(); itv!=v.end(); itv++)
               std::cout << itv->format() << " " ;
            std::cout << std::endl;
         }
      }
   }

}

std::vector<std::pair<mmdb::Atom *, std::string> >
coot::daca::make_typed_atoms(mmdb::Model *model_p, const coot::protein_geometry &geom) const {

   std::vector<std::pair<mmdb::Atom *, std::string> > v;

   std::map<std::string, dictionary_residue_restraints_t> dictionary_map;

   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               std::string res_type = residue_p->GetResName();
               std::map<std::string, dictionary_residue_restraints_t>::const_iterator it;
               it = dictionary_map.find(res_type);
               if (it == dictionary_map.end()) {
                  std::pair<bool, dictionary_residue_restraints_t> restraints =
                     geom.get_monomer_restraints(res_type, protein_geometry::IMOL_ENC_ANY);
                  if (restraints.first) {
                     dictionary_map[res_type] = restraints.second;
                  }
               }
            }
         }
      }
   }

   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               int n_atoms = residue_p->GetNumberOfAtoms();
               if (n_atoms > 0) {
                  std::string res_type = residue_p->GetResName();
                  std::map<std::string, dictionary_residue_restraints_t>::const_iterator it;
                  it = dictionary_map.find(res_type);
                  if (it != dictionary_map.end()) {
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (at) {
                           std::string atom_name(at->GetAtomName());
                           const std::string type = it->second.type_energy(atom_name);

                           // check  for atom name being "N" here //  C is Correct type
                           if (atom_name == " N  ") {
                              std::pair<mmdb::Atom *, std::string> p(at, "NH1");
                              v.push_back(p);
                           } else {
                              if (! type.empty()) {
                                 std::pair<mmdb::Atom *, std::string> p(at, type);
                                 v.push_back(p);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return v;
}

std::vector<std::pair<mmdb::Atom *, std::string> >
coot::daca::make_symmetry_typed_atoms(mmdb::Manager *mol,
                                      mmdb::Model *model_p,
                                      const protein_geometry &geom,
                                      float expansion_radius,
                                      std::vector<mmdb::Atom *> *symm_atom_store_p) const {

   std::vector<std::pair<mmdb::Atom *, std::string> > result;

   int n_symops = mol->GetNumberOfSymOps();
   if (n_symops == 0) {
      std::cout << "INFO:: no symmetry operators — symmetry contacts not included" << std::endl;
      return result;
   }

   // Build dictionary map for typing (same as make_typed_atoms)
   std::map<std::string, dictionary_residue_restraints_t> dictionary_map;
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               std::string res_type = residue_p->GetResName();
               if (dictionary_map.find(res_type) == dictionary_map.end()) {
                  auto restraints = geom.get_monomer_restraints(res_type, protein_geometry::IMOL_ENC_ANY);
                  if (restraints.first)
                     dictionary_map[res_type] = restraints.second;
               }
            }
         }
      }
   }

   // Compute bounding box of ASU atoms
   float xmin = 1e9f, xmax = -1e9f;
   float ymin = 1e9f, ymax = -1e9f;
   float zmin = 1e9f, zmax = -1e9f;
   int n_asu_atoms;
   mmdb::PPAtom asu_atoms = nullptr;
   mol->GetAtomTable(asu_atoms, n_asu_atoms);
   for (int i=0; i<n_asu_atoms; i++) {
      mmdb::Atom *at = asu_atoms[i];
      if (! at) continue;
      if (at->x < xmin) xmin = at->x;
      if (at->x > xmax) xmax = at->x;
      if (at->y < ymin) ymin = at->y;
      if (at->y > ymax) ymax = at->y;
      if (at->z < zmin) zmin = at->z;
      if (at->z > zmax) zmax = at->z;
   }
   // Expand bounding box by contact search distance
   xmin -= expansion_radius; xmax += expansion_radius;
   ymin -= expansion_radius; ymax += expansion_radius;
   zmin -= expansion_radius; zmax += expansion_radius;

   // ASU centre and radius for quick symop filtering
   float cx = 0.5f * (xmin + xmax + 2.0f * expansion_radius);
   float cy = 0.5f * (ymin + ymax + 2.0f * expansion_radius);
   float cz = 0.5f * (zmin + zmax + 2.0f * expansion_radius);
   // not needed after bbox expansion — cx/cy/cz already centred on original bbox

   // Determine cell translation range
   mmdb::realtype a[6], vol;
   int orthcode;
   mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
   float min_cell = static_cast<float>(std::min({a[0], a[1], a[2]}));
   if (min_cell < 1.0f) {
      std::cout << "WARNING:: unit cell dimension too small (" << min_cell
                << ") — symmetry contacts not included" << std::endl;
      return result;
   }
   int ishift = 1 + static_cast<int>(expansion_radius / min_cell);

   unsigned int n_symm_atoms = 0;
   for (int isym=0; isym<n_symops; isym++) {
      for (int xs=-ishift; xs<=ishift; xs++) {
         for (int ys=-ishift; ys<=ishift; ys++) {
            for (int zs=-ishift; zs<=ishift; zs++) {
               // Skip identity operation
               if (isym == 0 && xs == 0 && ys == 0 && zs == 0) continue;

               mmdb::mat44 mat;
               int err = mol->GetTMatrix(mat, isym, xs, ys, zs);
               if (err != 0) continue;

               // For each atom in the model, transform and check bounding box
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int nres = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (! residue_p) continue;
                     std::string res_type = residue_p->GetResName();
                     auto it_dict = dictionary_map.find(res_type);
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (! at) continue;

                        // Transform the atom position
                        float tx = mat[0][0]*at->x + mat[0][1]*at->y + mat[0][2]*at->z + mat[0][3];
                        float ty = mat[1][0]*at->x + mat[1][1]*at->y + mat[1][2]*at->z + mat[1][3];
                        float tz = mat[2][0]*at->x + mat[2][1]*at->y + mat[2][2]*at->z + mat[2][3];

                        // Quick bounding box test
                        if (tx < xmin || tx > xmax) continue;
                        if (ty < ymin || ty > ymax) continue;
                        if (tz < zmin || tz > zmax) continue;

                        // This atom is close enough — create a copy with transformed coords
                        mmdb::Atom *new_at = new mmdb::Atom;
                        new_at->Copy(at);
                        new_at->x = tx;
                        new_at->y = ty;
                        new_at->z = tz;
                        new_at->SetResidue(nullptr); // not part of the ASU
                        symm_atom_store_p->push_back(new_at);

                        // Type the atom
                        if (it_dict != dictionary_map.end()) {
                           std::string atom_name(at->GetAtomName());
                           if (atom_name == " N  ") {
                              result.push_back(std::make_pair(new_at, std::string("NH1")));
                           } else {
                              std::string type = it_dict->second.type_energy(atom_name);
                              if (! type.empty())
                                 result.push_back(std::make_pair(new_at, type));
                           }
                        }
                        n_symm_atoms++;
                     }
                  }
               }
            }
         }
      }
   }

   std::cout << "INFO:: added " << n_symm_atoms << " symmetry-related atoms from "
             << n_symops << " symmetry operators (ishift=" << ishift << ")" << std::endl;
   return result;
}

std::vector<std::vector<std::string> >
coot::daca::atom_names_for_fragments(const std::string &res_name) const {
   std::vector<std::vector<std::string> > v;

   std::vector<std::string> all{" CA ", " C  ", " O  "};
   v.push_back(all);

   if (res_name == "GLY") {
      std::vector<std::string> s{" N  ", " C  ", " CA "};
      v.push_back(s);
   }
   if (res_name == "ALA") {
      std::vector<std::string> s{" N  ", " C  ", " CA ", " CB "};
      v.push_back(s);
   }
   if (res_name == "CYS") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " SG "};
      v.push_back(s1);
      v.push_back(s2);
   }
   if (res_name == "ASP") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " OD1", " OD2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "GLU") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD "};
      std::vector<std::string> s4{" CG ", " CD ", " OE1", " OE2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
   }
   if (res_name == "PHE") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "HIS") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " ND1", " CE1", " NE2", " CD2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "ILE") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG1", " CG2"};
      std::vector<std::string> s3{" CB ", " CG1", " CD1"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "LYS") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD "};
      std::vector<std::string> s4{" CG ", " CD ", " CE "};
      std::vector<std::string> s5{" CD ", " CE ", " NZ "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
      v.push_back(s5);
   }
   if (res_name == "LEU") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD1", " CD2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "MET") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " SD "};
      std::vector<std::string> s4{" CG ", " SD ", " CE "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
   }
   if (res_name == "MSE") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " SE "};
      std::vector<std::string> s4{" CG ", " SE ", " CE "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
   }
   if (res_name == "ASN") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " OD1", " ND2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "PRO") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB ", " CG ", " CD "};
      v.push_back(s1);
   }
   if (res_name == "GLN") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD "};
      std::vector<std::string> s4{" CG ", " CD ", " OE1", " NE2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
   }
   if (res_name == "ARG") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD "};
      std::vector<std::string> s4{" CG ", " CD ", " NE "};
      std::vector<std::string> s5{" NE ", " CZ ", " NH1", " NH2"}; // + CD?
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
      v.push_back(s4);
      v.push_back(s5);
   }
   if (res_name == "SER") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " OG "};
      v.push_back(s1);
      v.push_back(s2);
   }
   if (res_name == "THR") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " OG1", " CG2"};
      v.push_back(s1);
      v.push_back(s2);
   }
   if (res_name == "VAL") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG1", " CG2"};
      v.push_back(s1);
      v.push_back(s2);
   }
   if (res_name == "TRP") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD1", " CD2", " NE1", " CE2", " CE3", " CZ2", " CZ3", " CH2"};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   if (res_name == "TYR") {
      std::vector<std::string> s1{" N  ", " C  ", " CA ", " CB "};
      std::vector<std::string> s2{" CA ", " CB ", " CG "};
      std::vector<std::string> s3{" CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " OH "};
      v.push_back(s1);
      v.push_back(s2);
      v.push_back(s3);
   }
   return v;
}

std::vector<std::vector<mmdb::Atom *> >
coot::daca::get_daca_fragments(mmdb::Residue *reference_residue_p) const {

   std::vector<std::vector<mmdb::Atom *> > v;
   std::string res_name(reference_residue_p->GetResName());
   std::vector<std::vector<std::string> > atom_name_vec_vec = atom_names_for_fragments(res_name);
   std::vector<std::vector<std::string> >::const_iterator it_1;
   std::vector<std::string>::const_iterator it_2;
   // Find the first alt conf in this residue, in case we need it as a fallback
   std::string first_alt_conf;
   {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      reference_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         std::string alt_loc(residue_atoms[iat]->altLoc);
         if (! alt_loc.empty()) {
            first_alt_conf = alt_loc;
            break;
         }
      }
   }

   for (it_1=atom_name_vec_vec.begin(); it_1!=atom_name_vec_vec.end(); it_1++){
      const std::vector<std::string> &atom_names = *it_1;
      std::vector<mmdb::Atom *> atom_vec;
      for (it_2=atom_names.begin(); it_2!=atom_names.end(); it_2++) {
         const std::string &atom_name = *it_2;
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms;
         reference_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            std::string this_atom_name(at->GetAtomName());
            if (atom_name == this_atom_name) {
               std::string alt_loc(at->altLoc);
               if (alt_loc.empty()) {
                  atom_vec.push_back(at);
                  break;
               }
            }
         }
      }
      // If we didn't find all atoms with blank alt conf, try the first alt conf
      if (atom_vec.size() != atom_names.size() && ! first_alt_conf.empty()) {
         atom_vec.clear();
         for (it_2=atom_names.begin(); it_2!=atom_names.end(); it_2++) {
            const std::string &atom_name = *it_2;
            mmdb::PPAtom residue_atoms = 0;
            int n_residue_atoms;
            reference_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               std::string this_atom_name(at->GetAtomName());
               if (atom_name == this_atom_name) {
                  std::string alt_loc(at->altLoc);
                  if (alt_loc.empty() || alt_loc == first_alt_conf) {
                     atom_vec.push_back(at);
                     break;
                  }
               }
            }
         }
      }
      if (atom_names.size() == atom_vec.size()) {
         v.push_back(atom_vec);
      }
   }
   return v;
}

std::pair<bool, clipper::RTop_orth>
coot::daca::get_frag_to_reference_rtop(const std::string &res_name,
                                       const unsigned int &frag_idx,
                                       const std::vector<mmdb::Atom *> &fragment_atoms) const {
   clipper::RTop_orth rtop;
   bool status = true;
   std::map<std::string, std::vector<std::vector<clipper::Coord_orth> > >::const_iterator it;
   it = reference_fragments.find(res_name);
   if (it != reference_fragments.end()) {
      if (frag_idx < it->second.size()) {
         const std::vector<clipper::Coord_orth> &ref_atom_positions = it->second[frag_idx];
         // convert fragment_atoms to coordinates
         std::vector<clipper::Coord_orth> residue_fragment_atoms;
         for (unsigned int i=0; i<fragment_atoms.size(); i++) {
            clipper::Coord_orth pos = co(fragment_atoms[i]);
            residue_fragment_atoms.push_back(pos);
         }
         if (ref_atom_positions.size() == residue_fragment_atoms.size()) {
            clipper::RTop_orth rtop_1(residue_fragment_atoms, ref_atom_positions);
            rtop = rtop_1;
         } else {
            std::cout << "size error in get_frag_to_reference_rtop()" << std::endl;
            status = false;
         }
      } else {
         std::cout << "index vector error in get_frag_to_reference_rtop() " << frag_idx << " res_name"
                   << std::endl;
         status = false;
      }

   } else {
      std::cout << "index residue error in get_frag_to_reference_rtop() " << frag_idx << " res_name"
                << std::endl;
      status = false;
   }
   return std::pair<bool, clipper::RTop_orth> (status, rtop);
}

void
coot::daca::add_to_box(mode_t mode,
                       const std::string &residue_type,
                       bool is_helical_flag,
                       unsigned int frag_index,
                       const box_index_t &box_index,
                       const std::string &atom_type,
                       unsigned int counts) {

   std::string box_key = residue_type + "-non-helical";
   if (is_helical_flag) box_key = residue_type + "-helical";

   if (mode == REFERENCE) {

      std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it =
         boxes.find(box_key);
      if (it == boxes.end()) {
         std::cout << "error in boxes " << box_key << std::endl;
      } else {
         boxes[box_key][frag_index][atom_type][box_index] += counts;
      }
   }

   if (mode == ANALYSIS) {
      if (frag_index >= boxes_for_testing[residue_type].size())
         boxes_for_testing[box_key].resize(6);
      boxes_for_testing[box_key][frag_index][atom_type][box_index] += counts;
   }

}

int
coot::daca::get_reference_counts(const std::string &residue_type,
                                 bool is_helical_flag,
                                 unsigned int frag_index,
                                 const box_index_t &box_index,
                                 const std::string &atom_type) const {

   int score = -1; // not found
   std::string box_key = residue_type + "-non-helical";
   if (is_helical_flag) box_key = residue_type + "-helical";

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it =
      boxes.find(box_key);
   if (it == boxes.end()) // should never happen
       return score;
   const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &frag_boxes = it->second;
   std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_typed_box =
      frag_boxes[frag_index].find(atom_type);
   if (it_typed_box != frag_boxes[frag_index].end()) {
      std::map<box_index_t, unsigned int>::const_iterator it_box = it_typed_box->second.find(box_index);
      if (it_box != it_typed_box->second.end()) {
         score = it_box->second;
      } else {
         std::cout << "Miss " << box_key << " " << frag_index << " " << atom_type << " "
                   << std::setw(2) << box_index.idx_x << " "
                   << std::setw(2) << box_index.idx_y << " "
                   << std::setw(2) << box_index.idx_z << " "
                   << std::endl;
      }
   } else {
      std::cout << "Miss:: " << box_key << " atom type " << atom_type << std::endl;
   }

   return score;

}

void
coot::daca::debug_boxes(const std::string &prefix) const {

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it;
   for (it=boxes.begin(); it!=boxes.end(); it++) {
      const std::string residue_type = it->first;
      std::cout << "========== debug_boxes(): " << prefix << " Residue Type " << residue_type << std::endl;
      const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &frag_boxes = it->second;
      std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > >::const_iterator it_v;
      for (unsigned int ifrag=0; ifrag<frag_boxes.size(); ifrag++) {
         const std::map<std::string, std::map<box_index_t, unsigned int> > &typed_boxes = frag_boxes[ifrag];
         std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_typed_box;
         for (it_typed_box=typed_boxes.begin(); it_typed_box!=typed_boxes.end(); it_typed_box++) {
            std::string atom_type = it_typed_box->first;
            if (residue_type.substr(0,3) == "ARG") {
               if (ifrag == 0) {
                  std::cout << "========== debug_boxes(): " << prefix << " Residue Type " << residue_type << " frag index "
                            << ifrag << " atom_type " << atom_type << std::endl;
                  if (true) {
                     std::map<box_index_t, unsigned int>::const_iterator it_box;
                     for (it_box=it_typed_box->second.begin(); it_box!=it_typed_box->second.end(); it_box++) {
                        const box_index_t &bi = it_box->first;
                        unsigned int count = it_box->second;
                        std::cout << " "
                                  << std::setw(2) << bi.idx_x << " " << std::setw(2) << bi.idx_y << " " << std::setw(2) << bi.idx_z << " "
                                  << std::setw(3) << count << std::endl;
                     }
                  }
               }
            }
         }
      }
   }
}


void
coot::daca::write_tables(const std::string &dir) const {

   std::cout << "write_tables(): write " << boxes.size() << " boxes " << std::endl;
   coot::util::create_directory(dir);

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it;
   for (it=boxes.begin(); it!=boxes.end(); ++it) {
      const std::string &residue_type = it->first;
      std::cout << "============= write_tables(): Residue Type " << residue_type << std::endl;
      const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &frag_boxes = it->second;
      for (unsigned int i=0; i<frag_boxes.size(); i++) {
         const std::map<std::string, std::map<box_index_t, unsigned int> > &typed_boxes = frag_boxes[i];
         std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_typed_box;
         for (it_typed_box=typed_boxes.begin(); it_typed_box!=typed_boxes.end(); ++it_typed_box) {
            std::string atom_type = it_typed_box->first;
            if (false)
               std::cout << "----------------- write_tables(): Residue Type " << residue_type << " " << i << " atom type "
                         << atom_type << std::endl;
            std::string box_file_name = residue_type + "-" + util::int_to_string(i) + "-" + atom_type + ".table";
            std::string full_box_file_name = coot::util::append_dir_file(dir, box_file_name);
            std::ofstream f(full_box_file_name.c_str());
            if (f) {
               std::map<box_index_t, unsigned int>::const_iterator it_box;
               for (it_box=it_typed_box->second.begin(); it_box!=it_typed_box->second.end(); ++it_box) {
                  const box_index_t &bi = it_box->first;
                  unsigned int count = it_box->second;
                  f << " "
                    << std::setw(2) << bi.idx_x << " " << std::setw(2) << bi.idx_y << " " << std::setw(2) << bi.idx_z << " "
                    << std::setw(3) << count << "\n";
               }
               f.close();
            }
         }
      }
   }
}

void
coot::daca::read_many_tables(const std::vector<std::string> &dirs) {
   presize_boxes();
   for (unsigned int i=0; i<dirs.size(); i++) {
      std::cout << "read tables directory " << dirs[i] << std::endl;
      read_tables(dirs[i]);
   }
}

void
coot::daca::read_tables(const std::string &dir) {

   if (! boxes_have_been_resized)
      presize_boxes();

   std::string glob_pattern = "*.table";
   std::vector<std::string> files = coot::util::glob_files(dir, glob_pattern);
   for (unsigned int i=0; i<files.size(); i++) {
      std::string file_name = files[i];
      // std::cout << "read table file " << file_name << std::endl;

      std::pair<std::string, std::string> z_parts = coot::util::split_string_on_last_slash(file_name);
      std::vector<std::string> fn_parts = coot::util::split_string(z_parts.second, "-");

      if (false) {
         std::cout << "fn_parts: " << std::endl;
         for (unsigned int i=0; i<fn_parts.size(); i++)
            std::cout << fn_parts[i] << " ";
         std::cout << std::endl;
      }

      if (fn_parts.size() == 4 || fn_parts.size() == 5) {
         try {
            std::string res_name = fn_parts[0];
            std::string ss_type = "helical";
            int ss_type_index = 0;
            unsigned int frag_string_index = 2;
            unsigned int atom_type_index = 3;
            bool is_helical_flag = true;
            if (fn_parts[1] == "non") {
               ss_type = "non-helical";
               ss_type_index = 1;
               frag_string_index = 3;
               atom_type_index = 4;
               is_helical_flag = false;
            }
            std::string frag_string = fn_parts[frag_string_index];
            int frag_index = coot::util::string_to_int(frag_string);
            const std::string &at_raw = fn_parts[atom_type_index];
            unsigned int l = at_raw.size();
            std::string atom_type = at_raw.substr(0,l-6);
            if (false)
               std::cout << " decoded: " << res_name << " " << ss_type << " " << frag_index
                         << " " << atom_type << std::endl;

            std::string line;
            std::vector<std::string> lines;
            std::ifstream f(files[i].c_str());
            while (std::getline(f, line)) {
               lines.push_back(line);
            }
            for (unsigned int j=0; j<lines.size(); j++) {
               const std::string &line = lines[j];
               std::vector<std::string> parts = coot::util::split_string_on_whitespace_no_blanks(line);
               if (parts.size() == 4) {
                  // .. x y z count
                  try {
                     int x = coot::util::string_to_int(parts[0]);
                     int y = coot::util::string_to_int(parts[1]);
                     int z = coot::util::string_to_int(parts[2]);
                     int c = coot::util::string_to_int(parts[3]);
                     box_index_t bi(x,y,z);
                     add_to_box(REFERENCE, res_name, is_helical_flag, frag_index, bi, atom_type, c);
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << "failed to parse " << line << " from " << files[i] << " " << rte.what() << std::endl;
                  }
               }
            }
         }
         catch (const std::runtime_error &rte) { }
      }
   }
}


void
coot::daca::fill_helix_flags(mmdb::Model *model_p, mmdb::Manager *mol) {

   std::vector<std::string> ch_ids;
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         ch_ids.push_back(chain_p->GetChainID());
      }
   }

   for (unsigned int ich=0; ich<ch_ids.size(); ich++) {
      int residue_selection_handle = mol->NewSelection();
      mol->Select (residue_selection_handle, mmdb::STYPE_RESIDUE, 0,
         ch_ids[ich].c_str(),
         mmdb::ANY_RES, "*",  // starting res
         mmdb::ANY_RES, "*",  // ending res
         "*",  // residue name
         "*",  // Residue must contain this atom name?
         "*",  // Residue must contain this Element?
         "*",  // altLocs
         mmdb::SKEY_NEW // selection key
      );

      std::vector<mmdb::Residue *> helical_residues_in_chain = like_a_helix(mol, residue_selection_handle);
      for (unsigned int i=0; i<helical_residues_in_chain.size(); i++)
         helical_residues.push_back(helical_residues_in_chain[i]);

      mol->DeleteSelection(residue_selection_handle);
   }
}

bool
coot::daca::atom_is_close_to_a_residue_atom(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const {
   float d_close = 1.7 + 1.7 + 1.5; // or so
   float dd_close = d_close * d_close;
   bool status = false;
   int n_residue_atoms;
   mmdb::PPAtom residue_atoms = 0;
   reference_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *ref_at = residue_atoms[iat];
      float dd =
         (at->x - ref_at->x) * (at->x - ref_at->x) +
         (at->y - ref_at->y) * (at->y - ref_at->y) +
         (at->z - ref_at->z) * (at->z - ref_at->z);
      if (dd < dd_close) {
         status = true;
         break;
      }
   }
   return status;
}


bool
coot::daca::atom_is_neighbour_mainchain(mmdb::Atom *at, mmdb::Residue *reference_residue_p) const {
   bool status = false;
   if (! at->residue) return false; // symmetry atom — not a neighbour
   int idx_res_1 = reference_residue_p->index;
   int idx_res_2 = at->residue->index;
   int idx_delta = abs(idx_res_2 - idx_res_1);
   if (idx_delta < 2) {
      std::string atom_name(at->GetAtomName());
      if (atom_name == " N  ") { status = true; }
      if (atom_name == " CA ") { status = true; }
      if (atom_name == " C  ") { status = true; }
      if (atom_name == " O  ") { status = true; }
   }
   return status;
}

void
coot::daca::presize_boxes(mode_t mode) {

   std::vector<std::string> residue_types = { "GLY", "ALA", "CYS", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU",
                                              "MET", "MSE", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP",
                                              "TYR"};

   if (mode == REFERENCE) {
      boxes_have_been_resized = true;
      for (auto type : residue_types) {
         const std::vector<std::string> h_types = {"-helical", "-non-helical"};
         for (auto h : h_types) {
            std::string key = type + h;
            boxes[key].resize(6);
         }
      }
   }
   if (mode == ANALYSIS) {
      for (auto type : residue_types) {
         const std::vector<std::string> h_types = {"-helical", "-non-helical"};
         for (auto h : h_types) {
            std::string key = type + h;
            boxes_for_testing[key].resize(6);
         }
      }
   }
}

#include "geometry/main-chain.hh"

float
coot::daca::calculate_daca(mmdb::Residue *reference_residue_p,
                           const std::vector<std::pair<mmdb::Atom *, std::string> > &typed_atoms,
                           coot::daca::mode_t mode) {

   bool print_scores = true;

   // brain-dead distance search (sad face)

   // has fill_helix_flags() been called before now?

   double d_crit = 8.0; // or something
   double dd_crit = d_crit * d_crit;

   float reference_counts = 0.0f;

   std::string res_name(reference_residue_p->GetResName());
   int reference_residue_seqnum = reference_residue_p->GetSeqNum();
   std::vector<std::vector<mmdb::Atom *> > fragments = get_daca_fragments(reference_residue_p);
   if (false)
      std::cout << "debug:: fragments.size() " << fragments.size() << " "
                << residue_spec_t(reference_residue_p)
                << " " << reference_residue_p->GetResName() << std::endl;
   for (unsigned int ifrag=0; ifrag<fragments.size(); ifrag++) {
      const std::vector<mmdb::Atom *> &atom_vec(fragments[ifrag]);
      std::vector<clipper::Coord_orth> reference_positions_vec;
      clipper::Coord_orth sum(0,0,0); // for calculating the centre of the fragment
      std::vector<mmdb::Atom *>::const_iterator it;
      for (it=atom_vec.begin(); it!=atom_vec.end(); it++) {
         clipper::Coord_orth pos = co(*it);
         reference_positions_vec.push_back(pos);
         sum += pos;
      }
      if (reference_positions_vec.size() > 2) {
         if (reference_positions_vec.size() == atom_vec.size()) {

            double m = 1.0/static_cast<double>(reference_positions_vec.size());
            clipper::Coord_orth frag_centre(sum * m);
            // Get the RTop that transforms the fragment to a reference
            // fragment at the origin.
            std::pair<bool, clipper::RTop_orth> frag_to_reference_rtop_pair =
               get_frag_to_reference_rtop(res_name, ifrag, atom_vec);
            if (! frag_to_reference_rtop_pair.first) continue;
            const clipper::RTop_orth &frag_to_reference_rtop = frag_to_reference_rtop_pair.second;
            for (unsigned int ita=0; ita<typed_atoms.size(); ita++) {
               mmdb::Atom *at = typed_atoms[ita].first;
               const std::string &atom_type = typed_atoms[ita].second;

               // don't consider atoms in this residue, of course
               if (at->residue == reference_residue_p)
                  continue;

               // don't consider peptide neighbour mainchain
               // (symmetry atoms have residue set to nullptr — skip this check for them)
               if (at->residue) {
                  int res_no_delta = at->residue->GetSeqNum() - reference_residue_seqnum;
                  if (std::abs(res_no_delta) < 2)
                     if (at->residue->chain == reference_residue_p->chain)
                        if (is_main_chain_p(at))
                           continue;
               }

               double dd =
                  (at->x - frag_centre.x()) * (at->x - frag_centre.x()) +
                  (at->y - frag_centre.y()) * (at->y - frag_centre.y()) +
                  (at->z - frag_centre.z()) * (at->z - frag_centre.z());
               if (dd < dd_crit) {
                  // Good, found something
                  if (atom_is_close_to_a_residue_atom(at, reference_residue_p)) {
                     if (! atom_is_neighbour_mainchain(at, reference_residue_p)) {
                        clipper::Coord_orth at_pos = co(at);
                        clipper::Coord_orth transformed_pos = frag_to_reference_rtop * at_pos;
                        box_index_t box_index(transformed_pos);
                        bool helical_flag = false;
                        if (std::find(helical_residues.begin(), helical_residues.end(), reference_residue_p) != helical_residues.end())
                        helical_flag = true;
                        if (mode == REFERENCE)
                           add_to_box(mode, res_name, helical_flag, ifrag, box_index, typed_atoms[ita].second);
                        if (mode == ANALYSIS) {
                           // Vriend & Sander V(d) envelope scoring:
                           // For each neighbouring box within radius R, compute
                           // V(d) = 1 - d^2/R^2, weighted by the reference count.
                           // R = 1.5 Angstroms (the box_width is 1.0).
                           std::string box_key = res_name + "-non-helical";
                           if (helical_flag) box_key = res_name + "-helical";
                           const float R = 1.5f;
                           const float R_sq = R * R;
                           float atom_score = 0.0f;
                           auto it_bk = boxes.find(box_key);
                           if (it_bk != boxes.end()) {
                              const auto &frag_boxes = it_bk->second;
                              auto it_at = frag_boxes[ifrag].find(atom_type);
                              if (it_at != frag_boxes[ifrag].end()) {
                                 const auto &box_map = it_at->second;
                                 for (int dx=-1; dx<=1; dx++) {
                                    for (int dy=-1; dy<=1; dy++) {
                                       for (int dz=-1; dz<=1; dz++) {
                                          box_index_t nb(box_index.idx_x + dx,
                                                         box_index.idx_y + dy,
                                                         box_index.idx_z + dz);
                                          auto it_nb = box_map.find(nb);
                                          if (it_nb != box_map.end()) {
                                             clipper::Coord_orth box_centre = nb.coord_orth();
                                             double d_sq =
                                                (transformed_pos.x() - box_centre.x()) * (transformed_pos.x() - box_centre.x()) +
                                                (transformed_pos.y() - box_centre.y()) * (transformed_pos.y() - box_centre.y()) +
                                                (transformed_pos.z() - box_centre.z()) * (transformed_pos.z() - box_centre.z());
                                             if (d_sq < R_sq) {
                                                float vd = 1.0f - static_cast<float>(d_sq) / R_sq;
                                                atom_score += vd * static_cast<float>(it_nb->second);
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                           reference_counts += atom_score;
                        }
                     }
                  }
               }
            }
         } else {
            std::cout << "OOps in atom set vs reference set size test " << std::endl;
         }
      } else {
         std::cout << "ERROR:: in calculate_daca(): This can't happen. reference positions size "
                   << reference_positions_vec.size() << " " << residue_spec_t(reference_residue_p)
                   << std::endl;
      }
   }
   return reference_counts;
}

void
coot::daca::write_tables_using_reference_structures_from_dir(const std::string &dir_name,
                                                             const std::string &output_tables_dir) {

   protein_geometry geom;
   geom.init_standard();
   std::vector<std::string> files = util::glob_files(dir_name, "*.pdb");
   // also search one level deeper (e.g. top2018_.../1d/1d1i/*.pdb)
   std::vector<std::string> files_sub = util::glob_files(dir_name, "*/*.pdb");
   files.insert(files.end(), files_sub.begin(), files_sub.end());

   std::cout << "in write_tables_using_reference_structures_from_dir() " << dir_name
             << " " << output_tables_dir << " found " << files.size() << " PDB files"
             << std::endl;

   presize_boxes(REFERENCE);
   for (unsigned int i=0; i<files.size(); i++) {
      std::string fn = files[i];
      atom_selection_container_t asc = get_atom_selection(fn, false, true, false);
      if (asc.read_success) {

         std::cout << "write_tables()... read pdb file " << fn << std::endl;

         if (false) { // bring this back when the consolidated tables are in  place.
            bool side_chain_only = false;
            std::vector<std::pair<mmdb::Residue *, float> > se = solvent_exposure(asc.mol, side_chain_only);
            for (unsigned int ii=0; ii<se.size(); ii++) {
               std::string rn(se[ii].first->GetResName());
               std::cout << "se " << fn << " " << coot::residue_spec_t(se[ii].first)
                         << " " << rn
                         << " " << se[ii].second << std::endl;
            }
         }

         mmdb::Model *model_p = asc.mol->GetModel(1);
         if (model_p) {
            fill_helix_flags(model_p, asc.mol);
            std::vector<std::pair<mmdb::Atom *, std::string> > ta = make_typed_atoms(model_p, geom);
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     std::string res_name(residue_p->GetResName());
                     if (res_name == "HOH") continue;
                     if (! util::is_standard_amino_acid_name(res_name)) continue;
                     calculate_daca(residue_p, ta, REFERENCE);
                  }
               }
            }
         }
      }
   }

   // debug_boxes("done-write-tables-using-reference-structures");

   write_tables(output_tables_dir);

}

void
coot::daca::score_molecule(const std::string &pdb_file_name) {

   std::cout << "score_molecule() " << pdb_file_name << std::endl;
   if (! coot::file_exists(pdb_file_name)) {
      std::cout << "No such file " << pdb_file_name << std::endl;
      return;
   }

   atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, false, false);
   if (! asc.read_success) return;
   mmdb::Model *model_p = asc.mol->GetModel(1);
   if (! model_p) return;

   std::cout << "INFO:: scoring " << pdb_file_name << std::endl;

   protein_geometry geom;
   geom.init_standard();
   presize_boxes(ANALYSIS);
   fill_helix_flags(model_p, asc.mol);

   std::vector<std::pair<mmdb::Residue *, float> > se = solvent_exposure(asc.mol);
   std::map<mmdb::Residue *, float> se_as_map;
   for (unsigned int i=0; i<se.size(); i++)
      se_as_map[se[i].first] = se[i].second;

   std::vector<std::pair<mmdb::Atom *, std::string> > ta = make_typed_atoms(model_p, geom);

   // Add symmetry-related atoms as contact partners
   std::vector<mmdb::Atom *> symm_atom_store;
   std::vector<std::pair<mmdb::Atom *, std::string> > symm_ta =
      make_symmetry_typed_atoms(asc.mol, model_p, geom, 10.0f, &symm_atom_store);
   ta.insert(ta.end(), symm_ta.begin(), symm_ta.end());

   // First pass: collect per-residue scores
   float se_half = 50.0f; // solvent exposure half-weight parameter
   struct residue_score_t {
      int res_number;
      std::string res_type;
      std::string ss_type;
      float raw_score;
      float se_score;
      float z_score;     // per-residue-type z-score
      float sef;         // solvent exposure factor: 1/(1 + se/se_half)
      float smoothed;    // sliding window average of z_score * sef
   };
   std::vector<residue_score_t> results;

   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) {
         mmdb::Residue *residue_p = chain_p->GetResidue(ires);
         if (! residue_p) continue;
         std::string res_name(residue_p->GetResName());
         if (res_name == "HOH") continue;
         if (! util::is_standard_amino_acid_name(res_name)) continue;

         residue_score_t rs;
         rs.res_number = residue_p->GetSeqNum();
         rs.res_type = res_name;
         rs.raw_score = calculate_daca(residue_p, ta, ANALYSIS);
         rs.se_score = -1.0f;
         auto it = se_as_map.find(residue_p);
         if (it != se_as_map.end())
            rs.se_score = it->second;
         rs.ss_type = "helix";
         if (std::find(helical_residues.begin(), helical_residues.end(), residue_p) == helical_residues.end())
            rs.ss_type = "non-helical";
         rs.z_score = 0.0f;
         rs.sef = (rs.se_score >= 0.0f) ? 1.0f / (1.0f + rs.se_score / se_half) : 0.0f;
         rs.smoothed = 0.0f;
         results.push_back(rs);
      }
   }

   // Per-residue-type normalization: compute mean and stddev, then z-scores
   std::map<std::string, std::vector<unsigned int> > type_indices;
   for (unsigned int i=0; i<results.size(); i++)
      type_indices[results[i].res_type].push_back(i);

   for (auto &ti : type_indices) {
      const std::vector<unsigned int> &indices = ti.second;
      if (indices.size() < 2) continue;
      double sum = 0.0;
      for (auto idx : indices)
         sum += results[idx].raw_score;
      double mean = sum / static_cast<double>(indices.size());
      double sum_sq = 0.0;
      for (auto idx : indices)
         sum_sq += (results[idx].raw_score - mean) * (results[idx].raw_score - mean);
      double stddev = std::sqrt(sum_sq / static_cast<double>(indices.size()));
      if (stddev > 1e-6) {
         for (auto idx : indices)
            results[idx].z_score = static_cast<float>((results[idx].raw_score - mean) / stddev);
      }
   }

   // Sliding window average of z_score * sef (window = 7, centred)
   // Skip residues with negative solvent exposure (missing data)
   int half_w = 3;
   for (int i=0; i<static_cast<int>(results.size()); i++) {
      float wsum = 0.0f;
      int count = 0;
      for (int j=i-half_w; j<=i+half_w; j++) {
         if (j >= 0 && j < static_cast<int>(results.size())) {
            if (results[j].se_score >= 0.0f) {
               wsum += results[j].z_score * results[j].sef;
               count++;
            }
         }
      }
      if (count > 0)
         results[i].smoothed = wsum / static_cast<float>(count);
   }

   // Output
   float cumulative = 0.0f;
   for (unsigned int i=0; i<results.size(); i++) {
      const residue_score_t &rs = results[i];
      cumulative += rs.raw_score;
      std::cout << "residue_number " << rs.res_number
                << " type " << rs.res_type
                << " SS-type " << rs.ss_type
                << " score " << std::fixed << std::setprecision(1) << rs.raw_score
                << " z_score " << std::fixed << std::setprecision(2) << rs.z_score
                << " sef " << std::fixed << std::setprecision(3) << rs.sef
                << " smoothed " << std::fixed << std::setprecision(2) << rs.smoothed
                << " daca_sum_score " << std::fixed << std::setprecision(1) << cumulative
                << " solvent_exposure " << rs.se_score
                << "\n";
   }

   // Clean up symmetry atoms
   for (unsigned int i=0; i<symm_atom_store.size(); i++)
      delete symm_atom_store[i];
}

void
coot::daca::make_data_for_figure_2(const std::string &pdb_dir) {

   protein_geometry geom;
   geom.init_standard();

   std::vector<std::string> files = util::glob_files(pdb_dir, "*.pdb");
   std::vector<std::string> files_sub = util::glob_files(pdb_dir, "*/*.pdb");
   files.insert(files.end(), files_sub.begin(), files_sub.end());
   std::vector<std::string> ent_files = util::glob_files(pdb_dir, "*.ent");
   files.insert(files.end(), ent_files.begin(), ent_files.end());
   std::vector<std::string> ent_files_sub = util::glob_files(pdb_dir, "*/*.ent");
   files.insert(files.end(), ent_files_sub.begin(), ent_files_sub.end());

   std::cout << "fig2: found " << files.size() << " PDB files in " << pdb_dir << std::endl;

   // header line
   std::cout << "fig2-data: res_type ss_type score solvent_exposure" << std::endl;

   for (unsigned int ifile=0; ifile<files.size(); ifile++) {
      const std::string &fn = files[ifile];
      atom_selection_container_t asc = get_atom_selection(fn, false, false, false);
      if (! asc.read_success) continue;
      mmdb::Model *model_p = asc.mol->GetModel(1);
      if (! model_p) continue;

      fill_helix_flags(model_p, asc.mol);
      std::vector<std::pair<mmdb::Atom *, std::string> > ta = make_typed_atoms(model_p, geom);

      std::vector<std::pair<mmdb::Residue *, float> > se = solvent_exposure(asc.mol);
      std::map<mmdb::Residue *, float> se_as_map;
      for (unsigned int i=0; i<se.size(); i++)
         se_as_map[se[i].first] = se[i].second;

      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (! residue_p) continue;
            std::string res_name(residue_p->GetResName());
            if (res_name == "HOH") continue;
            if (! util::is_standard_amino_acid_name(res_name)) continue;

            float daca_score = calculate_daca(residue_p, ta, ANALYSIS);

            float se_score = -1.0f;
            auto it = se_as_map.find(residue_p);
            if (it != se_as_map.end())
               se_score = it->second;
            if (se_score < 0.0f) continue;

            std::string ss_type = "helix";
            if (std::find(helical_residues.begin(), helical_residues.end(), residue_p) == helical_residues.end())
               ss_type = "non-helical";

            std::cout << "fig2-data: " << res_name << " " << ss_type << " "
                      << std::fixed << std::setprecision(1) << daca_score << " "
                      << std::fixed << std::setprecision(1) << se_score
                      << std::endl;
         }
      }

      if ((ifile + 1) % 10 == 0)
         std::cerr << "fig2: processed " << (ifile + 1) << "/" << files.size() << std::endl;
   }
}

std::pair<int, int>
coot::daca::self_test(const std::string &pdb_file_name) {

   int n_total = 0;
   int n_misses = 0;

   if (! coot::file_exists(pdb_file_name)) {
      std::cout << "self_test(): No such file " << pdb_file_name << std::endl;
      return std::make_pair(-1, -1);
   }

   atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, false);
   if (! asc.read_success) {
      std::cout << "self_test(): Failed to read " << pdb_file_name << std::endl;
      return std::make_pair(-1, -1);
   }

   mmdb::Model *model_p = asc.mol->GetModel(1);
   if (! model_p) {
      std::cout << "self_test(): No model in " << pdb_file_name << std::endl;
      return std::make_pair(-1, -1);
   }

   protein_geometry geom;
   geom.init_standard();

   fill_helix_flags(model_p, asc.mol);
   std::vector<std::pair<mmdb::Atom *, std::string> > ta = make_typed_atoms(model_p, geom);

   // REFERENCE pass — populate boxes from this structure
   presize_boxes(REFERENCE);
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) {
         mmdb::Residue *residue_p = chain_p->GetResidue(ires);
         if (residue_p) {
            std::string res_name(residue_p->GetResName());
            if (res_name == "HOH") continue;
            if (! util::is_standard_amino_acid_name(res_name)) continue;
            calculate_daca(residue_p, ta, REFERENCE);
         }
      }
   }

   // ANALYSIS pass — replay the same residue loop.
   // For each contact, get_reference_counts() must return > 0 because
   // the reference distribution was built from this same structure.
   // We replicate the inner logic of calculate_daca() here so that
   // we can count hits and misses directly.

   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) {
         mmdb::Residue *residue_p = chain_p->GetResidue(ires);
         if (! residue_p) continue;
         std::string res_name(residue_p->GetResName());
         if (res_name == "HOH") continue;
         if (! util::is_standard_amino_acid_name(res_name)) continue;

         int seqnum = residue_p->GetSeqNum();
         std::vector<std::vector<mmdb::Atom *> > fragments = get_daca_fragments(residue_p);
         for (unsigned int ifrag=0; ifrag<fragments.size(); ifrag++) {
            const std::vector<mmdb::Atom *> &atom_vec(fragments[ifrag]);
            if (atom_vec.size() < 3) continue;
            std::pair<bool, clipper::RTop_orth> rtop_pair =
               get_frag_to_reference_rtop(res_name, ifrag, atom_vec);
            if (! rtop_pair.first) continue;
            const clipper::RTop_orth &rtop = rtop_pair.second;

            // fragment centre
            clipper::Coord_orth sum(0,0,0);
            for (unsigned int ia=0; ia<atom_vec.size(); ia++)
               sum += co(atom_vec[ia]);
            double m = 1.0/static_cast<double>(atom_vec.size());
            clipper::Coord_orth frag_centre(sum * m);

            bool helical_flag = (std::find(helical_residues.begin(),
                                           helical_residues.end(),
                                           residue_p) != helical_residues.end());

            for (unsigned int ita=0; ita<ta.size(); ita++) {
               mmdb::Atom *at = ta[ita].first;
               if (at->residue == residue_p) continue;
               int res_no_delta = at->residue->GetSeqNum() - seqnum;
               if (std::abs(res_no_delta) < 2)
                  if (at->residue->chain == residue_p->chain)
                     if (is_main_chain_p(at))
                        continue;
               double dd =
                  (at->x - frag_centre.x()) * (at->x - frag_centre.x()) +
                  (at->y - frag_centre.y()) * (at->y - frag_centre.y()) +
                  (at->z - frag_centre.z()) * (at->z - frag_centre.z());
               if (dd < 64.0) { // 8.0^2
                  if (atom_is_close_to_a_residue_atom(at, residue_p)) {
                     if (! atom_is_neighbour_mainchain(at, residue_p)) {
                        clipper::Coord_orth transformed_pos = rtop * co(at);
                        box_index_t box_index(transformed_pos);
                        int counts = get_reference_counts(res_name, helical_flag,
                                                          ifrag, box_index,
                                                          ta[ita].second);
                        n_total++;
                        if (counts <= 0)
                           n_misses++;
                     }
                  }
               }
            }
         }
      }
   }

   std::cout << "self_test(): n_total_contacts " << n_total << " n_misses " << n_misses << std::endl;
   return std::make_pair(n_total, n_misses);
}

std::pair<int, int>
coot::daca::self_test_perturbed(const std::string &pdb_file_name, float perturbation) {

   int n_total = 0;
   int n_misses = 0;

   if (! coot::file_exists(pdb_file_name)) {
      std::cout << "self_test_perturbed(): No such file " << pdb_file_name << std::endl;
      return std::make_pair(-1, -1);
   }

   atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, false);
   if (! asc.read_success) {
      std::cout << "self_test_perturbed(): Failed to read " << pdb_file_name << std::endl;
      return std::make_pair(-1, -1);
   }

   mmdb::Model *model_p = asc.mol->GetModel(1);
   if (! model_p) {
      std::cout << "self_test_perturbed(): No model in " << pdb_file_name << std::endl;
      return std::make_pair(-1, -1);
   }

   protein_geometry geom;
   geom.init_standard();

   fill_helix_flags(model_p, asc.mol);
   std::vector<std::pair<mmdb::Atom *, std::string> > ta = make_typed_atoms(model_p, geom);

   // REFERENCE pass — populate boxes from the unperturbed structure
   presize_boxes(REFERENCE);
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) {
         mmdb::Residue *residue_p = chain_p->GetResidue(ires);
         if (residue_p) {
            std::string res_name(residue_p->GetResName());
            if (res_name == "HOH") continue;
            if (! util::is_standard_amino_acid_name(res_name)) continue;
            calculate_daca(residue_p, ta, REFERENCE);
         }
      }
   }

   // Perturb all atom coordinates
   srand(42); // fixed seed for reproducibility
   int n_atoms;
   mmdb::PPAtom all_atoms = nullptr;
   asc.mol->GetAtomTable(all_atoms, n_atoms);
   for (int i=0; i<n_atoms; i++) {
      mmdb::Atom *at = all_atoms[i];
      if (at) {
         float dx = perturbation * (2.0f * static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 1.0f);
         float dy = perturbation * (2.0f * static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 1.0f);
         float dz = perturbation * (2.0f * static_cast<float>(rand()) / static_cast<float>(RAND_MAX) - 1.0f);
         at->x += dx;
         at->y += dy;
         at->z += dz;
      }
   }

   // ANALYSIS pass on the perturbed structure
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) {
         mmdb::Residue *residue_p = chain_p->GetResidue(ires);
         if (! residue_p) continue;
         std::string res_name(residue_p->GetResName());
         if (res_name == "HOH") continue;
         if (! util::is_standard_amino_acid_name(res_name)) continue;

         int seqnum = residue_p->GetSeqNum();
         std::vector<std::vector<mmdb::Atom *> > fragments = get_daca_fragments(residue_p);
         for (unsigned int ifrag=0; ifrag<fragments.size(); ifrag++) {
            const std::vector<mmdb::Atom *> &atom_vec(fragments[ifrag]);
            if (atom_vec.size() < 3) continue;
            std::pair<bool, clipper::RTop_orth> rtop_pair =
               get_frag_to_reference_rtop(res_name, ifrag, atom_vec);
            if (! rtop_pair.first) continue;
            const clipper::RTop_orth &rtop = rtop_pair.second;

            // fragment centre
            clipper::Coord_orth sum(0,0,0);
            for (unsigned int ia=0; ia<atom_vec.size(); ia++)
               sum += co(atom_vec[ia]);
            double m = 1.0/static_cast<double>(atom_vec.size());
            clipper::Coord_orth frag_centre(sum * m);

            bool helical_flag = (std::find(helical_residues.begin(),
                                           helical_residues.end(),
                                           residue_p) != helical_residues.end());

            for (unsigned int ita=0; ita<ta.size(); ita++) {
               mmdb::Atom *at = ta[ita].first;
               if (at->residue == residue_p) continue;
               int res_no_delta = at->residue->GetSeqNum() - seqnum;
               if (std::abs(res_no_delta) < 2)
                  if (at->residue->chain == residue_p->chain)
                     if (is_main_chain_p(at))
                        continue;
               double dd =
                  (at->x - frag_centre.x()) * (at->x - frag_centre.x()) +
                  (at->y - frag_centre.y()) * (at->y - frag_centre.y()) +
                  (at->z - frag_centre.z()) * (at->z - frag_centre.z());
               if (dd < 64.0) { // 8.0^2
                  if (atom_is_close_to_a_residue_atom(at, residue_p)) {
                     if (! atom_is_neighbour_mainchain(at, residue_p)) {
                        clipper::Coord_orth transformed_pos = rtop * co(at);
                        box_index_t box_index(transformed_pos);
                        int counts = get_reference_counts(res_name, helical_flag,
                                                          ifrag, box_index,
                                                          ta[ita].second);
                        n_total++;
                        if (counts <= 0)
                           n_misses++;
                     }
                  }
               }
            }
         }
      }
   }

   std::cout << "self_test_perturbed(): perturbation " << perturbation
             << " n_total_contacts " << n_total << " n_misses " << n_misses << std::endl;
   return std::make_pair(n_total, n_misses);
}

void
coot::daca::compare_boxes() const {

   unsigned int n_daca = 0;
   unsigned int n_hits = 0;
   unsigned int sum = 0;

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it;
   for (it =boxes_for_testing.begin(); it!=boxes_for_testing.end(); it++) {
      const std::string &res_name_with_ss(it->first);
      const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &v(it->second);
      for (unsigned int idx_frag=0; idx_frag<v.size(); idx_frag++) {
         const std::map<std::string, std::map<box_index_t, unsigned int> > &m1(v[idx_frag]);
         std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_1;
         for (it_1=m1.begin(); it_1!=m1.end(); it_1++) {
            const std::string &atom_type = it_1->first;
            const std::map<box_index_t, unsigned int> &m2(it_1->second);
            std::map<box_index_t, unsigned int>::const_iterator it_2;
            for (it_2=m2.begin(); it_2!=m2.end(); it_2++) {
               const box_index_t &bi = it_2->first;
               const unsigned int &count_analysis = it_2->second;
               n_daca++;

               // does there exist a reference count for that?
               std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::const_iterator it_ref;
               it_ref = boxes.find(res_name_with_ss);
               if (it_ref == boxes.end()) {
                  std::cout << "Failed to find reference for type " << res_name_with_ss << std::endl;
               } else {
                  const std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &v_ref(it_ref->second);
                  if (! v_ref.empty()) {
                     const std::map<std::string, std::map<box_index_t, unsigned int> > &m1_ref(v_ref[idx_frag]); // vector sized correctly?
                     std::map<std::string, std::map<box_index_t, unsigned int> >::const_iterator it_1_ref;
                     it_1_ref = m1_ref.find(atom_type);
                     if (it_1_ref == m1_ref.end()) {
                        std::cout << "Failed to find reference for type " << res_name_with_ss
                                  << " frag-index " << idx_frag << " atom-type " << atom_type
                                  << " we have map size " << m1_ref.size() << std::endl;
                     } else {
                        const std::map<box_index_t, unsigned int> &m2_ref(it_1_ref->second);
                        std::map<box_index_t, unsigned int>::const_iterator it_2_ref;
                        it_2_ref = m2_ref.find(bi);
                        if (it_2_ref == m2_ref.end()) {
                           std::cout << "Failed to find reference for " << res_name_with_ss << " "
                                     << idx_frag << " " << atom_type << " box_index "
                                     << bi.idx_x << " " << bi.idx_y << " " << bi.idx_z << std::endl;

                        } else {
                           int count_ref = it_2_ref->second;
                           if (false)
                              std::cout << res_name_with_ss << " " << idx_frag << " " << atom_type << " "
                                        << bi.idx_x << " " << bi.idx_y << " " << bi.idx_z << " "
                                        << count_ref << "\n";
                           sum += count_ref;
                           // std::cout << "sum " << sum << "\n";
                           n_hits++;
                        }
                     }
                  } else {
                     std::cout << "v_ref is empty for " << it_ref->first << std::endl;
                  }
               }
            }
         }
      }
   }
   std::cout << "compare_boxes() n_daca " << n_daca << " n_hits " << n_hits
             << " sum " << sum << std::endl;
}


void
coot::daca::normalize() {

   // iterators not const because we want to modify the contents of the boxes
   //
   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::iterator it;
   for (it =boxes.begin(); it!=boxes.end(); it++) {
      const std::string &res_name_with_ss(it->first);
      std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &v(it->second);
      for (unsigned int idx_frag=0; idx_frag<v.size(); idx_frag++) {
         std::map<std::string, std::map<box_index_t, unsigned int> > &m1(v[idx_frag]);
         std::map<std::string, std::map<box_index_t, unsigned int> >::iterator it_1;
         for (it_1=m1.begin(); it_1!=m1.end(); it_1++) {
            const std::string &atom_type = it_1->first;
            std::map<box_index_t, unsigned int> &m2(it_1->second);
            unsigned int n_count_sum = 0;
            std::map<box_index_t, unsigned int>::const_iterator it_2;
            for (it_2=m2.begin(); it_2!=m2.end(); it_2++) {
               const box_index_t &bi = it_2->first;
               const unsigned int &counts = it_2->second;
               n_count_sum += counts;
            }

            if (false)
               std::cout << "normalize " << res_name_with_ss << " "
                         << "frag-index " << idx_frag << " "
                         << "atom_type " << atom_type << " "
                         << n_count_sum << std::endl;
            float scale_factor = 1000000.0f / static_cast<float>(n_count_sum);
            std::map<box_index_t, unsigned int>::iterator it_counts;
            for (it_counts=m2.begin(); it_counts!=m2.end(); it_counts++) {
               unsigned int counts = it_counts->second;
               it_counts->second = static_cast<int>(scale_factor * static_cast<float>(counts));
            }
         }
      }
   }
}

void
coot::daca::normalize_v2() {

   // For each (res_type+SS, frag_index, atom_type), normalize the spatial
   // distribution so that counts sum to a fixed value (1,000,000).
   // This converts raw accumulated counts into a probability-like distribution,
   // so that structures contributing different numbers of residues to the
   // reference database don't dominate.

   const float target_sum = 1000000.0f;
   unsigned int n_distributions = 0;

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::iterator it;
   for (it=boxes.begin(); it!=boxes.end(); it++) {
      std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &v(it->second);
      for (unsigned int idx_frag=0; idx_frag<v.size(); idx_frag++) {
         std::map<std::string, std::map<box_index_t, unsigned int> > &m1(v[idx_frag]);
         std::map<std::string, std::map<box_index_t, unsigned int> >::iterator it_at;
         for (it_at=m1.begin(); it_at!=m1.end(); it_at++) {
            std::map<box_index_t, unsigned int> &m2(it_at->second);
            unsigned int n_count_sum = 0;
            for (auto it_box=m2.begin(); it_box!=m2.end(); it_box++)
               n_count_sum += it_box->second;
            if (n_count_sum > 0) {
               float scale_factor = target_sum / static_cast<float>(n_count_sum);
               for (auto it_box=m2.begin(); it_box!=m2.end(); it_box++)
                  it_box->second = static_cast<unsigned int>(scale_factor * static_cast<float>(it_box->second));
               n_distributions++;
            }
         }
      }
   }
   std::cout << "normalize_v2(): normalized " << n_distributions << " distributions" << std::endl;
}

// 2aaa get the pdb redo model, build on this and check with privateer.

// as in the verb
void
coot::daca::envelope() {

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::iterator it;
   for (it =boxes.begin(); it!=boxes.end(); ++it) {
      std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &v(it->second);
      for (unsigned int idx_frag=0; idx_frag<v.size(); idx_frag++) {
         std::map<std::string, std::map<box_index_t, unsigned int> > &m1(v[idx_frag]);
         std::map<std::string, std::map<box_index_t, unsigned int> >::iterator it_1;
         for (it_1=m1.begin(); it_1!=m1.end(); ++it_1) {
            std::map<box_index_t, unsigned int> &m2(it_1->second);
            std::map<box_index_t, unsigned int>::iterator it_2;
            for (it_2=m2.begin(); it_2!=m2.end(); ++it_2) {
               const box_index_t &bi = it_2->first;
               const unsigned int &counts = it_2->second;
               unsigned int c = counts;
               float dd = bi.d_squared();
               float scale_factor = gompertz_scale(dd);
               it_2->second = static_cast<int>(scale_factor * static_cast<float>(counts));
               if (false) // checking that the scaling is sane
                  std::cout << "d " << sqrt(dd) << " scale " << scale_factor << " " << c
                            << " " << it_2->second << std::endl;
            }
         }
      }
   }
}

void
coot::daca::smooth() {

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >::iterator it;

   std::map<std::string, std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > >
      copy_boxes = boxes;

   for (it=boxes.begin(); it!=boxes.end(); ++it) {
      const std::string &res_name_with_ss(it->first);
      std::vector<std::map<std::string, std::map<box_index_t, unsigned int> > > &v(it->second);
      for (unsigned int idx_frag=0; idx_frag<v.size(); idx_frag++) {
         std::map<std::string, std::map<box_index_t, unsigned int> > &m1(v[idx_frag]);
         std::map<std::string, std::map<box_index_t, unsigned int> >::iterator it_1;
         for (it_1=m1.begin(); it_1!=m1.end(); ++it_1) {
            const std::string &atom_type = it_1->first;
            std::map<box_index_t, unsigned int> &m2(it_1->second);
            std::map<box_index_t, unsigned int>::const_iterator it_2;
            for (it_2=m2.begin(); it_2!=m2.end(); ++it_2) {
               const box_index_t &bi = it_2->first;
               const unsigned int &counts = it_2->second;
               const int box_idx_min = -8;
               const int box_idx_max =  7;
               for (int delta_x = -1; delta_x <= 1; delta_x++) {
                  for (int delta_y = -1; delta_y <= 1; delta_y++) {
                     for (int delta_z = -1; delta_z <= 1; delta_z++) {
                        if ((delta_x == 0) && (delta_y == 0) && (delta_z == 0)) {
                           // do nothing
                        } else {
                           int idx_neighb_x = bi.idx_x + delta_x;
                           int idx_neighb_y = bi.idx_y + delta_y;
                           int idx_neighb_z = bi.idx_z + delta_z;
                           if (idx_neighb_x >= box_idx_min) {
                              if (idx_neighb_x <= box_idx_max) {
                                 if (idx_neighb_y >= box_idx_min) {
                                    if (idx_neighb_y <= box_idx_max) {
                                       if (idx_neighb_z >= box_idx_min) {
                                          if (idx_neighb_z <= box_idx_max) {
                                             box_index_t neighb_box_index(idx_neighb_x, idx_neighb_y, idx_neighb_z);
                                             // int contrib = static_cast<int>(0.125 * static_cast<float>(counts));
                                             int contrib = counts;
                                             copy_boxes[res_name_with_ss][idx_frag][atom_type][neighb_box_index] += contrib;
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   boxes = copy_boxes;
}

double
coot::daca::get_radius(const std::string &ele) const {

   double radius = 1.70;
   if (ele == " H")
      radius = 1.20;
   if (ele == " N")
      radius = 1.55;
   if (ele == " O")
      radius = 1.52;
   if (ele == " S")
      radius = 1.8;
   // PDBv3
   if (ele == "H")
      radius = 1.20;
   if (ele == "N")
      radius = 1.55;
   if (ele == "O")
      radius = 1.52;
   if (ele == "S")
      radius = 1.8;
   return radius;
}


#include "coot-utils/fib-sphere.hh"

// eposure of the side chain atoms only are considered
//
std::vector<std::pair<mmdb::Residue *, float> >
coot::daca::solvent_exposure(mmdb::Manager *mol, bool side_chain_only) const {

   std::vector<std::pair<mmdb::Residue *, float> > v; // return the residue count map, not this
   if (! mol) return v;

   float max_dist = 2 * (1.7 + 1.4);
   mmdb::PPAtom atom_selection = 0;
   int n_atoms;

   int SelHnd = mol->NewSelection(); // d
   mol->SelectAtoms(SelHnd, 1, "*",
                    mmdb::ANY_RES, "*",
                    mmdb::ANY_RES, "*",
                    "*","*","!H","*", mmdb::SKEY_NEW);

   std::map<int, std::set<int> > contact_map;

   // fill contact_map
   mol->GetSelIndex(SelHnd, atom_selection, n_atoms);
   if (n_atoms) {
      mmdb::Contact *pscontact = NULL; // d
      int n_contacts;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
	 for (int j=0; j<4; j++)
	    my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mol->SeekContacts(atom_selection, n_atoms,
			atom_selection, n_atoms,
			0, max_dist,
			0, // allow contacts from atoms in same residue
			pscontact, n_contacts,
			0, &my_matt, i_contact_group); // makes reverses also
      if (n_contacts > 0) {

	 if (pscontact) {
            // we could do the selection wthout waters, but also filter out waters this way
            std::vector<bool> is_water(n_atoms, false);
            for (int iat=0; iat<n_atoms; iat++) {
               std::string rn(atom_selection[iat]->residue->GetResName());
               if (rn == "HOH")
                  is_water[iat] = true;
            }
            std::vector<float> radius(n_atoms);
            for (int iat=0; iat<n_atoms; iat++) {
               std::string ele(atom_selection[iat]->element);
               radius[iat] = get_radius(ele); // could be more clever, use atom type.
            }
	    for (int i=0; i<n_contacts; i++) {
               // check for not a water here. (both ways)
               if (is_water[pscontact[i].id1]) continue;
               if (is_water[pscontact[i].id2]) continue;

               mmdb::Atom *at = atom_selection[pscontact[i].id1];
               if (is_main_chain_p(at)) continue;
               contact_map[pscontact[i].id1].insert(pscontact[i].id2);
            }
            delete [] pscontact;


            // OK, contact_map is filled.

            // so now let's count the dots around each of the atoms

            // for faster execution, need a function that takes 2 atoms and a squared dist
            // (so that there is no conversion to doubles)
            // bool within_distance_criterion(mmdb::Atom *at_1, mmdb::Atom *at_2, float dist_sq);

            auto dot_count = [] (int atom_index,
                                 const std::set<int> &neighbour_atoms,
                                 float atom_radius,
                                 mmdb::PPAtom atom_selection,
                                 const std::vector<clipper::Coord_orth> &unit_sphere_points) {
                                int count = 0;
                                const float water_radius = 1.8;
                                float radius = atom_radius + water_radius;
                                double dd_crit = radius * radius;
                                mmdb::Atom *at = atom_selection[atom_index];
                                clipper::Coord_orth atom_position = co(at);
                                std::set<int>::const_iterator it;
                                for (unsigned int i=0; i<unit_sphere_points.size(); i++) {
                                   bool inside_another_atom = false;
                                   clipper::Coord_orth pt_on_sphere(radius * unit_sphere_points[i]);
                                   clipper::Coord_orth pt = atom_position + pt_on_sphere;
                                   for (it=neighbour_atoms.begin(); it!=neighbour_atoms.end(); ++it) {
                                      mmdb::Atom *at_neigb = atom_selection[*it];
                                      clipper::Coord_orth pt_neighb = co(at_neigb);
                                      double dd = (pt-pt_neighb).lengthsq();
                                      if (false)
                                         std::cout << i << " pt-on-sphere " << pt_on_sphere.format()
                                                   << " pt " << pt.format()
                                                   << " d " << sqrt(dd) << " " << dd << std::endl;
                                      if (dd < dd_crit) {
                                         inside_another_atom = true;
                                         break;
                                      }
                                   }
                                   if (! inside_another_atom)
                                      count += 1;
                                }
                                if (false) {
                                   std::string rn(at->residue->GetResName());
                                   if (rn == "HOH") {
                                      std::cout << "HOH " << residue_spec_t(at->residue)
                                                << " with n-neighbs " << neighbour_atoms.size()
                                                << " returning " << count << std::endl;
                                   }
                                }
                                return count;
                             };

            const unsigned int n_sphere_points = 106; // so that each point covers about 1A^2
            std::vector<clipper::Coord_orth> unit_sphere_points =
               coot::fibonacci_sphere(n_sphere_points);

            std::map<mmdb::Residue *, int> residue_count_map;
            std::map<int, std::set<int> >::const_iterator it;
            for (it=contact_map.begin(); it!=contact_map.end(); ++it) {
               int atom_index = it->first;
               mmdb::Atom *at = atom_selection[atom_index];
               const std::set<int> &neighbours = it->second;
               int n_dots_for_atom = 0;
               if (neighbours.size() == 0)
                  n_dots_for_atom = n_sphere_points; // nothing can block the sphere points
               else
                  n_dots_for_atom = dot_count(atom_index, neighbours, radius[atom_index],
                                              atom_selection, unit_sphere_points);

               residue_count_map[at->residue] += n_dots_for_atom;
            }

            {
               if (false) {
                  std::cout << "contact map:"  << std::endl;
                  std::map<mmdb::Residue *, int>::const_iterator it_rc;
                  for (it_rc=residue_count_map.begin(); it_rc!=residue_count_map.end(); ++it_rc) {
                     std::string rn = it_rc->first->GetResName();
                     std::cout << "    " << residue_spec_t(it->first) << " " << rn << " "
                               << it_rc->second << std::endl;
                  }
               }

               std::map<mmdb::Residue *, int>::const_iterator it_rc;
               for (it_rc=residue_count_map.begin(); it_rc!=residue_count_map.end(); ++it_rc) {
                  std::pair<mmdb::Residue *, float> p(it_rc->first, it_rc->second);
                  v.push_back(p);
               }
            }
         }
      }
   }

   return v;
}

std::vector<std::pair<mmdb::Residue *, float> >
coot::daca::solvent_exposure_old_version_v2(mmdb::Manager *mol,
                                            bool side_chain_only) const {

   std::vector<std::pair<mmdb::Residue *, float> > v;
   if (! mol) return v;

   mmdb::PPAtom atom_selection = 0;
   int n_atoms;

   int SelHnd = mol->NewSelection(); // d
   mol->SelectAtoms(SelHnd, 1, "*",
                    mmdb::ANY_RES, "*",
                    mmdb::ANY_RES, "*",
                    "*","*","!H","*", mmdb::SKEY_NEW);

   std::map<mmdb::Residue *, std::set<mmdb::Atom *> > residue_neighbouring_atoms;

   mol->GetSelIndex(SelHnd, atom_selection, n_atoms);
   if (n_atoms) {

      float max_dist = 5.7;
      mmdb::Contact *pscontact = NULL; // d
      int n_contacts;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
	 for (int j=0; j<4; j++)
	    my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mol->SeekContacts(atom_selection, n_atoms,
			atom_selection, n_atoms,
			0, max_dist,
			0, // in same residue
			pscontact, n_contacts,
			0, &my_matt, i_contact_group); // makes reverses also
      if (n_contacts > 0) {

	 if (pscontact) {
	    for (int i=0; i<n_contacts; i++) {
               mmdb::Atom *at_1 = atom_selection[pscontact[i].id1];
               mmdb::Atom *at_2 = atom_selection[pscontact[i].id2];

               mmdb::Residue *residue_p_1 = at_1->GetResidue();
               mmdb::Residue *residue_p_2 = at_2->GetResidue();
               if (residue_p_2 == residue_p_1) continue;
               std::string res_name_1(residue_p_1->GetResName());
               std::string res_name_2(at_2->residue->GetResName());
               if (res_name_1 == "HOH") continue;
               if (res_name_2 == "HOH") continue;
               if (! util::is_standard_amino_acid_name(res_name_1)) continue;

               if (! side_chain_only)
                  residue_neighbouring_atoms[residue_p_1].insert(at_2);
               else
                  if (!is_main_chain_p(at_1))
                     residue_neighbouring_atoms[residue_p_1].insert(at_2);
            }
            delete [] pscontact;
         }
      }

      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               std::map<mmdb::Residue *, std::set<mmdb::Atom *> >::const_iterator it =
                  residue_neighbouring_atoms.find(residue_p);
               if (it != residue_neighbouring_atoms.end()) {
                  std::pair<mmdb::Residue *, float> p(residue_p, it->second.size());
                  v.push_back(p);
               }
            }
         }
      }
   }
   mol->DeleteSelection(SelHnd);
   return v;
}


std::vector<std::pair<mmdb::Atom *, float> >
coot::daca::solvent_exposure_old_version(int SelHnd_in, mmdb::Manager *mol) const {

   std::vector<std::pair<mmdb::Atom *, float> > v;
   if (mol) {

      double dot_density = 0.5;
      //
      double phi_step = 5.0 * (M_PI/180.0);
      double theta_step = 5.0 * (M_PI/180.0);
      phi_step   /= dot_density;
      theta_step /= dot_density;

      double water_radius = 1.4;
      double fudge = 1.0;
      mmdb::PPAtom atoms = 0;
      int n_atoms;
      mol->GetSelIndex(SelHnd_in, atoms, n_atoms);
      std::vector<double> radius(n_atoms);

      for (int iat=0; iat<n_atoms; iat++) {
	 std::string ele(atoms[iat]->element);
	 radius[iat] = get_radius(ele);
      }

      mmdb::PPAtom atoms_all = 0;
      int n_atoms_all;
      int SelHnd_all = mol->NewSelection();
      mol->SelectAtoms(SelHnd_all, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
      mol->GetSelIndex(SelHnd_all, atoms_all, n_atoms_all);

      for (int iatom=0; iatom<n_atoms; iatom++) {
	 if (! atoms[iatom]->isTer()) {
	    clipper::Coord_orth centre(atoms[iatom]->x,
				       atoms[iatom]->y,
				       atoms[iatom]->z);
	    bool even = 1;
	    int n_points = 0;
	    int n_sa = 0;
	    for (double theta=0; theta<M_PI; theta+=theta_step) {
	       double phi_step_inner = phi_step + 0.1 * pow(theta-0.5*M_PI, 2);
	       for (double phi=0; phi<2*M_PI; phi+=phi_step_inner) {
		  if (even) {
		     double r = fudge * (radius[iatom] + water_radius);
		     clipper::Coord_orth pt(r*cos(phi)*sin(theta),
					    r*sin(phi)*sin(theta),
					    r*cos(theta));
		     pt += centre;
		     n_points++;

		     // now, is pt closer to (a water centre around)
		     // another atom?

		     bool is_solvent_accessible = true;
		     for (int i_all=0; i_all<n_atoms_all; i_all++) {
			// don't exclude from self
			mmdb::Atom *other_at = atoms_all[i_all];
			std::string other_res_name = other_at->GetResName();
			if (other_res_name != "HOH") {
			   if (atoms[iatom] != other_at) {
			      std::string other_ele = other_at->element;
			      if (other_ele != " H") {
				 double other_atom_r = fudge * (get_radius(other_ele) + water_radius);
				 double other_atom_r_sq = other_atom_r * other_atom_r;
				 clipper::Coord_orth pt_other(other_at->x, other_at->y, other_at->z);
				 if ((pt-pt_other).lengthsq() < other_atom_r_sq) {
				    is_solvent_accessible = 0;
				    break;
				 }
			      }
			   }
			}
		     }
		     if (is_solvent_accessible)
			n_sa++;
		  }
		  even = 1 - even;
	       }
	    }

	    double exposure_frac = double(n_sa)/double(n_points);
	    if (0)
	       std::cout << "Atom " << atoms[iatom]->name << " has exposure " << n_sa << "/" << n_points
			 << " = " << exposure_frac << std::endl;
	    std::pair<mmdb::Atom *, float> p(atoms[iatom], exposure_frac);
	    v.push_back(p);
	 }
      }
      mol->DeleteSelection(SelHnd_all); // presumably this was missing before... 20101230
   }
   return v;
}


void
coot::daca::cook() {

   // smooth();
   // envelope();
   normalize_v2();
}
