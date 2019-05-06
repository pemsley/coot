/* ligand/dipole.cc
 * 
 * Copyright 2009 by The University of Oxford
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
 * 02110-1301, USA.
 */

// no dependency on coords files
// #include "coords/mmdb-extras.h" // this and next
// #include "coords/mmdb.h"        // for nice printing of atom pointer

#include <stdexcept>
#include "utils/coot-utils.hh" // for int_to_string()

#include "dipole.hh"

coot::dipole::dipole() {
   dipole_is_good_flag = 0;
} 

// Thow an exception on failure to make a dipole
// 
coot::dipole::dipole(const coot::dictionary_residue_restraints_t &rest,
		     mmdb::Residue *residue_p) {
   
   
   std::vector<std::pair<coot::dictionary_residue_restraints_t, mmdb::Residue *> > dict_res_pairs;
   std::pair<coot::dictionary_residue_restraints_t, mmdb::Residue *> p(rest, residue_p);
   dict_res_pairs.push_back(p);
   init(dict_res_pairs);
}

coot::dipole::dipole(std::vector<std::pair<coot::dictionary_residue_restraints_t, mmdb::Residue *> > dict_res_pairs) {
   init(dict_res_pairs);
}

// This needs not to fail with a residue of 1 atom, because we need
// fill_charged_atoms to work for the partially charged surface.
// 
void
coot::dipole::init(std::vector<std::pair<coot::dictionary_residue_restraints_t, mmdb::Residue *> > dict_res_pairs) {

   dipole_is_good_flag = 0;

   // Find the centre of the residue
   double sum_x = 0; 
   double sum_y = 0; 
   double sum_z = 0;
   int n_points = 0;

   for (unsigned int ires=0; ires<dict_res_pairs.size(); ires++) {
      mmdb::Residue *residue_p = dict_res_pairs[ires].second;
      coot::dictionary_residue_restraints_t rest = dict_res_pairs[ires].first;

      mmdb::PPAtom SelAtoms;
      int nSelAtoms;
      residue_p->GetAtomTable(SelAtoms, nSelAtoms);

      for (int i_res_at=0; i_res_at<nSelAtoms; i_res_at++) {
	 sum_x += SelAtoms[i_res_at]->x;
	 sum_y += SelAtoms[i_res_at]->y;
	 sum_z += SelAtoms[i_res_at]->z;
	 n_points++;
      }
   }

   if (n_points == 0) {
      std::string mess = "No atoms in ";
      mess += coot::util::int_to_string(dict_res_pairs.size());
      mess += " residue";
      if (dict_res_pairs.size() != 1)
	 mess += "s";
      mess += " ";
      for (unsigned int i=0; i<dict_res_pairs.size(); i++) {
         mmdb::Residue *residue_p = dict_res_pairs[i].second;
         mess += residue_p->GetChainID();
         mess += " ";
         mess += residue_p->GetSeqNum();
         mess += ", ";
      }
      throw std::runtime_error(mess);
   }


   double multiplier = 1.0/double(n_points);

   // class member:
   residue_centre = clipper::Coord_orth(sum_x*multiplier,
					sum_y*multiplier,
					sum_z*multiplier);
   
   std::vector<std::pair<mmdb::Atom *, float> > charged_ats = charged_atoms(dict_res_pairs);
   std::vector<std::pair<float, clipper::Coord_orth> > charged_points(charged_ats.size());
   for (unsigned int i=0; i<charged_ats.size(); i++) {
      clipper::Coord_orth p(charged_ats[i].first->x,
			    charged_ats[i].first->y,
			    charged_ats[i].first->z);
      charged_points[i] =
	 std::pair<float, clipper::Coord_orth>(charged_ats[i].second, p);
   }

   clipper::Coord_orth dip(0,0,0);
   for (unsigned int ii=0; ii<charged_points.size(); ii++) { 
      clipper::Coord_orth scaled(charged_points[ii].second - residue_centre);
	 scaled = charged_points[ii].first * scaled;
	 dip += scaled;
	 dipole_is_good_flag = 1;
   }

   if (! dipole_is_good_flag) {
      std::string mess = "Dipole is not good for ";
      mess += coot::util::int_to_string(dict_res_pairs.size());
      mess += " residue";
      if (dict_res_pairs.size() != 1)
	 mess += "s";
      mess += " ";
      for (unsigned int i=0; i<dict_res_pairs.size(); i++) {
         mmdb::Residue *residue_p = dict_res_pairs[i].second;
         mess += residue_p->GetChainID();
         mess += " ";
         mess += coot::util::int_to_string(residue_p->GetSeqNum());
         mess += " ";
         mess += residue_p->GetResName();
         mess += ", ";
      }
      throw std::runtime_error(mess);
   }
      
   dipole_ = dip;
} 

void
coot::dipole::fill_charged_atoms(mmdb::Residue *residue_p,
				 const coot::dictionary_residue_restraints_t &rest) {

   std::vector<std::pair<mmdb::Atom *, float> > v = charged_atoms(residue_p, rest);
   for (unsigned int i=0; i<v.size(); i++) {
      v[i].first->charge = v[i].second;
   }
}

std::vector<std::pair<mmdb::Atom *, float> >
coot::dipole::charged_atoms(mmdb::Residue *residue_p,
			    const coot::dictionary_residue_restraints_t &rest) const {

   std::vector<std::pair<mmdb::Atom *, float> > charged_ats;
   mmdb::PPAtom SelAtoms;
   int nSelAtoms;
   residue_p->GetAtomTable(SelAtoms, nSelAtoms);
   int n_dict_atom = rest.atom_info.size();
	 	 
   for (int i_res_at=0; i_res_at<nSelAtoms; i_res_at++) {
      bool found_match = 0;
      mmdb::Atom *at = SelAtoms[i_res_at];
      std::string atom_name = at->name;
      for (int j=0; j<n_dict_atom; j++) {
	 if (rest.atom_info[j].partial_charge.first) {
	    if (atom_name == rest.atom_info[j].atom_id_4c) {
	       std::pair<mmdb::Atom *, float> p(at, rest.atom_info[j].partial_charge.second);
	       charged_ats.push_back(p);
	       break; 
	    }
	 } else {
	    std::cout << "    no partial charge for "
		      << rest.atom_info[j].atom_id << std::endl;
	 }
      }
   }
   return charged_ats;
}



std::vector<std::pair<mmdb::Atom *, float> >
coot::dipole::charged_atoms(std::vector<std::pair<coot::dictionary_residue_restraints_t, mmdb::Residue *> > dict_res_pairs) const {

   std::vector<std::pair<mmdb::Atom *, float> > charged_ats;

   for (unsigned int i=0; i<dict_res_pairs.size(); i++) {
      mmdb::Residue *residue_p = dict_res_pairs[i].second;
      coot::dictionary_residue_restraints_t rest = dict_res_pairs[i].first;
      std::vector<std::pair<mmdb::Atom *, float> > residue_charged_ats = charged_atoms(residue_p, rest);
      for (unsigned int j=0; j<residue_charged_ats.size(); j++) {
         charged_ats.push_back(residue_charged_ats[j]);
      }
   }
   return charged_ats;
}


clipper::Coord_orth
coot::dipole::get_unit_dipole() const {

   clipper::Coord_orth d = dipole_;
   double l = d.lengthsq();
   double sc = 1/sqrt(l);
   return sc * d;
} 
		    

std::ostream&
coot::operator<<(std::ostream &s, const coot::dipole &d) {

   s << "[" << d.dipole_.x() << " " << d.dipole_.y() << " " << d.dipole_.z()
     << "]";

   return s;
}
