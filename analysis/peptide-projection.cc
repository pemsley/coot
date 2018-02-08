/* analysis/peptide-projection.cc
 * 
 * Copyright 2016 by Medical Research Council
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

#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coot-utils/coot-coord-utils.hh"

#include "stats.hh"

mmdb::Manager *get_manager(const std::string &file_name) {

   mmdb::ERROR_CODE err;
   mmdb::Manager *mol = new mmdb::Manager;
   err = mol->ReadCoorFile(file_name.c_str());
   
   if (err) {
      std::cout << "There was an error reading " << file_name << ". \n";
      std::cout << "ERROR " << err << " READ: "
		<< mmdb::GetErrorDescription(err) << std::endl;
      int  error_count;
      char error_buf[500];
      mol->GetInputBuffer(error_buf, error_count);
      if (error_count >= 0) { 
	 std::cout << "         LINE #" << error_count << "\n     "
		   << error_buf << std::endl << std::endl;
      } else {
	 if (error_count == -1) { 
	    std::cout << "       CIF ITEM: " << error_buf
		      << std::endl << std::endl;
	 }
      }
      delete mol;
      mol = NULL;
   }

   return mol;
}

bool
is_trans_peptide_by_distance(const clipper::Coord_orth &p1,
			     const clipper::Coord_orth &p2) {

   double d = clipper::Coord_orth::length(p1, p2);

   bool r = false;
   if (d < 3.9) {
      if (d > 3.6) {
	 r = true;
      }
   }
   return r;
}
			     

std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
get_trans_CAs(mmdb::Manager *mol) {

   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v;

   // for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) { 
      mmdb::Chain *chain_p;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_this;
	 mmdb::Residue *residue_next;
	 mmdb::Atom *at_this;
	 mmdb::Atom *at_next;
	 for (int ires=0; ires<(nres-1); ires++) { 
	    residue_this = chain_p->GetResidue(ires);
	    residue_next = chain_p->GetResidue(ires+1);
	    int n_atoms_this = residue_this->GetNumberOfAtoms();

	    // PDBv3 FIXME
	    at_this = residue_this->GetAtom(" CA ");
	    at_next = residue_next->GetAtom(" CA ");

	    if (at_this) { 
	       if (at_next) {
		  std::string alt_conf_a = at_this->altLoc;
		  std::string alt_conf_b = at_next->altLoc;
		  if (alt_conf_a == "") { 
		     if (alt_conf_b == "") {
			if ((at_this->GetSeqNum() +1) == at_next->GetSeqNum()) { 
			   clipper::Coord_orth pos_1 = coot::co(at_this);
			   clipper::Coord_orth pos_2 = coot::co(at_next);
			   if (is_trans_peptide_by_distance(pos_1, pos_2)) {
			      std::pair<mmdb::Atom *, mmdb::Atom *> p(at_this, at_next);
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
   
   return v;
}

void analyse_distances(const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &v) {

   std::cout << "# analyse_distances() " << v.size() << std::endl;
   for (unsigned int i=0; i<v.size(); i++) {
      clipper::Coord_orth ca_1_pos = coot::co(v[i].first);
      clipper::Coord_orth ca_2_pos = coot::co(v[i].second);
      mmdb::Atom *c_1 = v[i].first->residue->GetAtom(" C  ");
      mmdb::Atom *o_1 = v[i].first->residue->GetAtom(" O  ");
      mmdb::Atom *n_2 = v[i].second->residue->GetAtom(" N  ");

      if (c_1) {
	 if (o_1) {
	    if (n_2) {
	       clipper::Coord_orth c_1_pos = coot::co(c_1);
	       clipper::Coord_orth o_1_pos = coot::co(o_1);
	       clipper::Coord_orth n_2_pos = coot::co(n_2);
	       double ca_ca_dist = clipper::Coord_orth::length(ca_1_pos, ca_2_pos);
	       clipper::Coord_orth ca_vec = ca_2_pos - ca_1_pos;
	       clipper::Coord_orth ca_c_vec = c_1_pos - ca_1_pos;
	       clipper::Coord_orth ca_n_vec = n_2_pos - ca_1_pos;
	       clipper::Coord_orth ca_2_n_vec = n_2_pos - ca_2_pos;
	       double dp_c = clipper::Coord_orth::dot(ca_vec, ca_c_vec);
	       double dp_n = clipper::Coord_orth::dot(ca_vec, ca_n_vec);
	       double dp_n_2 = clipper::Coord_orth::dot(ca_vec, ca_2_n_vec);
	       double projection_length_c = dp_c/ca_ca_dist;
	       double projection_length_n = dp_n/ca_ca_dist;
	       double projection_length_n_2 = dp_n_2/ca_ca_dist;
	       std::cout << "   "
			 << projection_length_c << " " 
			 << projection_length_n << " " 
			 << projection_length_n_2 << " " 
			 << "\n";
	    }
	 }
      }
   }
}

clipper::Coord_orth
mid_point(const clipper::Coord_orth &pt_1,
	  const clipper::Coord_orth &pt_2) {

   return clipper::Coord_orth(0.5*(pt_1.x() + pt_2.x()),
			      0.5*(pt_1.y() + pt_2.y()),
			      0.5*(pt_1.z() + pt_2.z()));
}


clipper::Coord_orth
fraction_point(const clipper::Coord_orth &pt_1,
	       const clipper::Coord_orth &pt_2,
	       double fr) {

   double m = 1.0 - fr;
   return clipper::Coord_orth(fr*pt_1.x() + m*pt_2.x(),
			      fr*pt_1.y() + m*pt_2.y(),
			      fr*pt_1.z() + m*pt_2.z());
}

double analyse_mid_points_fp(const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &v,
			     double fraction_point_1,
			     double fraction_point_2) {

   double mean;
   std::vector<double> ds;
   ds.reserve(100);
   // std::cout << "# analyse_distances() " << v.size() << std::endl;
   for (unsigned int i=0; i<v.size(); i++) { 
      clipper::Coord_orth ca_1_pos = coot::co(v[i].first);
      clipper::Coord_orth ca_2_pos = coot::co(v[i].second);
      mmdb::Atom *c_1 = v[i].first->residue->GetAtom(" C  ");
      mmdb::Atom *o_1 = v[i].first->residue->GetAtom(" O  ");
      mmdb::Atom *n_2 = v[i].second->residue->GetAtom(" N  ");

      if (c_1) {
	 if (o_1) {
	    if (n_2) {
	       clipper::Coord_orth c_1_pos = coot::co(c_1);
	       clipper::Coord_orth o_1_pos = coot::co(o_1);
	       clipper::Coord_orth n_2_pos = coot::co(n_2);
	       double ca_ca_dist = clipper::Coord_orth::length(ca_1_pos, ca_2_pos);

	       clipper::Coord_orth mp_ca = fraction_point(ca_1_pos, ca_2_pos, fraction_point_1);
	       clipper::Coord_orth mp_nc = fraction_point(c_1_pos,   n_2_pos, fraction_point_2);

	       double d = sqrt((mp_ca-mp_nc).lengthsq());

	       // std::cout << " " << d << std::endl;

	       ds.push_back(d);
	    }
	 }
      }
   }

   coot::stats::single s(ds);
   return s.mean();
}

void analyse_mid_points(const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &v) {

   double m = analyse_mid_points_fp(v, 0.5, 0.485);
   // std::cout << "m " << m << std::endl;

   int n_samples = 30;
   for (int i= -n_samples; i<n_samples; i++) {
      double fraction_point_i = 0.5 + double(i) / double(n_samples) * 0.2;
      for (int j= -n_samples; j<n_samples; j++) {
	 double fraction_point_j = 0.5 + double(j) / double(n_samples) * 0.2;
	 double m = analyse_mid_points_fp(v, fraction_point_i, fraction_point_j);
	 std::cout << " " << fraction_point_i << " " << fraction_point_j << " " << m << "\n";
      }
   }
}



int main(int argc, char **argv) {

   int status = 0;

   if (argc > 1) {
      std::string file_name = argv[1];
      mmdb::Manager *mol = get_manager(file_name);
      if(mol) {
	 std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v = get_trans_CAs(mol);
	 // analyse_distances(v);
	 analyse_mid_points(v);
      }
   }
   return status;
}
