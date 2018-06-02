/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
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

// #define ANALYSE_REFINEMENT_TIMING

// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL

#include <string.h> // for strcmp

#ifdef ANALYSE_REFINEMENT_TIMING
#include <sys/time.h> // for gettimeofday()
#endif // ANALYSE_REFINEMENT_TIMING

#include <fstream>
#include <algorithm> // for sort
#include <stdexcept>
#include <iomanip>
#ifdef HAVE_CXX_THREAD
#include <thread>
#include <chrono>
#endif // HAVE_CXX_THREAD

#include "simple-restraint.hh"

double
coot::geometry_distortion_info_container_t::print() const {

   double total_distortion = 0.0;
   
   int n_restraints_bonds    = 0;
   int n_restraints_angles   = 0;
   int n_restraints_torsions = 0;
   int n_restraints_chirals  = 0;
   int n_restraints_planes   = 0;
   double sum_penalties_bonds    = 0;
   double sum_penalties_angles   = 0;
   double sum_penalties_planes   = 0;
   double sum_penalties_chirals  = 0;
   double sum_penalties_torsions = 0;
   std::vector<std::pair<std::string,double> > penalty_string_bonds;
   std::vector<std::pair<std::string,double> > penalty_string_angles;
   std::vector<std::pair<std::string,double> > penalty_string_torsions;
   std::vector<std::pair<std::string,double> > penalty_string_planes;
   std::cout << "Residue Distortion List: \n";
   for (unsigned int i=0; i< geometry_distortion.size(); i++) { 
      const coot::simple_restraint &rest = geometry_distortion[i].restraint;
      if (rest.restraint_type == coot::BOND_RESTRAINT) {
	 n_restraints_bonds++;
	 mmdb::Atom *at_1 = atom[rest.atom_index_1];
	 mmdb::Atom *at_2 = atom[rest.atom_index_2];
	 if (at_1 && at_2) {
	    clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
	    double d = sqrt((p2-p1).lengthsq());
	    double distortion = d - rest.target_value;
	    double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
	    std::string chain_res_1 = std::string(at_1->GetChainID()) + " " + coot::util::int_to_string(at_1->GetSeqNum());
	    std::string chain_res_2 = std::string(at_2->GetChainID()) + " " + coot::util::int_to_string(at_2->GetSeqNum());
	    std::string s = std::string("bond ")
	       + chain_res_1 + " "
	       + std::string(at_1->name) + std::string(" to ")
	       + chain_res_2 + " "
	       + std::string(at_2->name)
	       + std::string(" target_value: ") + coot::util::float_to_string_using_dec_pl(rest.target_value, 3)
	       + std::string(" d: ") + coot::util::float_to_string_using_dec_pl(d, 3)
	       + std::string(" sigma: ") + coot::util::float_to_string_using_dec_pl(rest.sigma, 3)
	       + std::string(" length-devi ") + coot::util::float_to_string_using_dec_pl(distortion, 3)
	       + std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
	    penalty_string_bonds.push_back(std::pair<std::string,double> (s, pen_score));
	    sum_penalties_bonds += pen_score;
	 }
      }

      if (rest.restraint_type == coot::ANGLE_RESTRAINT) {
	 n_restraints_angles++;
	 mmdb::Atom *at_1 = atom[rest.atom_index_1];
	 mmdb::Atom *at_2 = atom[rest.atom_index_2];
	 mmdb::Atom *at_3 = atom[rest.atom_index_3];
	 if (at_1 && at_2 && at_3) {
	    clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
	    clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
	    double angle_rad = clipper::Coord_orth::angle(p1, p2, p3);
	    double angle = clipper::Util::rad2d(angle_rad);
	    double distortion = angle - rest.target_value;
	    double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
	    std::string chain_res_1 = std::string(at_1->GetChainID()) + " " + coot::util::int_to_string(at_1->GetSeqNum());
	    std::string chain_res_2 = std::string(at_2->GetChainID()) + " " + coot::util::int_to_string(at_2->GetSeqNum());
	    std::string chain_res_3 = std::string(at_3->GetChainID()) + " " + coot::util::int_to_string(at_3->GetSeqNum());
	    std::string s = std::string("angle ")
	       + chain_res_1 + " "
	       + std::string(at_1->name)
	       + std::string(" - ")
	       + chain_res_2 + " "
	       + std::string(at_2->name)
	       + std::string(" - ")
	       + chain_res_3 + " "
	       + std::string(at_3->name)
	       + std::string("  target: ") + coot::util::float_to_string(rest.target_value)
	       + std::string(" model_angle: ") + coot::util::float_to_string(angle)
	       + std::string(" sigma: ") + coot::util::float_to_string(rest.sigma)
	       + std::string(" angle-devi ") + coot::util::float_to_string(distortion)
	       + std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
	    penalty_string_angles.push_back(std::pair<std::string,double> (s, pen_score));
	    sum_penalties_angles += pen_score;
	 }
      }

      if (rest.restraint_type == TORSION_RESTRAINT) {
	 n_restraints_torsions++;
	 mmdb::Atom *at_1 = atom[rest.atom_index_1];
	 mmdb::Atom *at_2 = atom[rest.atom_index_2];
	 mmdb::Atom *at_3 = atom[rest.atom_index_3];
	 mmdb::Atom *at_4 = atom[rest.atom_index_4];
	 if (at_1 && at_2 && at_3 && at_4) {
	    clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
	    clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
	    clipper::Coord_orth p4(at_4->x, at_4->y, at_4->z);
	    double torsion_rad = clipper::Coord_orth::torsion(p1, p2, p3, p4);
	    double torsion = clipper::Util::rad2d(torsion_rad);
	    double distortion = rest.torsion_distortion(torsion);
	    double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
	    std::string chain_res_1 = std::string(at_1->GetChainID()) + " " + coot::util::int_to_string(at_1->GetSeqNum());
	    std::string chain_res_2 = std::string(at_2->GetChainID()) + " " + coot::util::int_to_string(at_2->GetSeqNum());
	    std::string chain_res_3 = std::string(at_3->GetChainID()) + " " + coot::util::int_to_string(at_3->GetSeqNum());
	    std::string chain_res_4 = std::string(at_4->GetChainID()) + " " + coot::util::int_to_string(at_4->GetSeqNum());
	    std::string s = std::string("torsion ")
	       + chain_res_1 + " "
	       + std::string(at_1->name) + std::string(" - ")
	       + chain_res_2 + " "
	       + std::string(at_2->name) + std::string(" - ")
	       + chain_res_3 + " "
	       + std::string(at_3->name) + std::string(" - ")
	       + chain_res_4 + " "
	       + std::string(at_4->name)
	       + std::string("  target: ") + coot::util::float_to_string(rest.target_value)
	       + std::string(" model_torsion: ") + coot::util::float_to_string(torsion)
	       + std::string(" sigma: ") + coot::util::float_to_string(rest.sigma)
	       + std::string(" torsion-devi ") + coot::util::float_to_string(distortion)
	       + std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
	    penalty_string_torsions.push_back(std::pair<std::string,double> (s, pen_score));
	    sum_penalties_angles += pen_score;
	 }
      }

      if (rest.restraint_type == TRANS_PEPTIDE_RESTRAINT) {
	 mmdb::Atom *at_1 = atom[rest.atom_index_1];
	 mmdb::Atom *at_2 = atom[rest.atom_index_2];
	 mmdb::Atom *at_3 = atom[rest.atom_index_3];
	 mmdb::Atom *at_4 = atom[rest.atom_index_4];
	 if (at_1 && at_2 && at_3 && at_4) {
	    clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
	    clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
	    clipper::Coord_orth p4(at_4->x, at_4->y, at_4->z);
	    double torsion_rad = clipper::Coord_orth::torsion(p1, p2, p3, p4);
	    double torsion = clipper::Util::rad2d(torsion_rad);
	    double distortion = rest.torsion_distortion(torsion);
	    // std::cout << "---------- restraint is trans_peptide_restraint with per "
	    // << rest.periodicity << " model " << torsion
	    // << " torsion_distortion: " << distortion << std::endl;
	    std::string chain_res_1 = std::string(at_1->GetChainID()) + " " + coot::util::int_to_string(at_1->GetSeqNum());
	    std::string chain_res_2 = std::string(at_2->GetChainID()) + " " + coot::util::int_to_string(at_2->GetSeqNum());
	    std::string chain_res_3 = std::string(at_3->GetChainID()) + " " + coot::util::int_to_string(at_3->GetSeqNum());
	    std::string chain_res_4 = std::string(at_4->GetChainID()) + " " + coot::util::int_to_string(at_4->GetSeqNum());
	    std::string s = std::string("trans-peptide ")
	       + chain_res_1 + " "
	       + std::string(at_1->name) + std::string(" - ")
	       + chain_res_2 + " "
	       + std::string(at_2->name) + std::string(" - ")
	       + chain_res_3 + " "
	       + std::string(at_3->name) + std::string(" - ")
	       + chain_res_4 + " "
	       + std::string(at_4->name)
	       + std::string("  target: ")        + coot::util::float_to_string(rest.target_value)
	       + std::string(" model_torsion: ")  + coot::util::float_to_string(torsion)
	       + std::string(" sigma: ")          + coot::util::float_to_string(rest.sigma)
	       + std::string(" penalty-score:  ") + coot::util::float_to_string_using_dec_pl(distortion, 3);
	    penalty_string_torsions.push_back(std::pair<std::string,double> (s, distortion));
	    sum_penalties_angles += distortion;
	 }
      }
	    
      if (rest.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	 n_restraints_chirals++;
	 double chiral_limit = 2.0;  // (based on histogram of CVs of A chain of tutorial.)
	 // chiral_limit = 0.0;
	 if (geometry_distortion[i].distortion_score > chiral_limit) {
	    mmdb::Atom *at_c = atom[rest.atom_index_centre];
	    mmdb::Atom *at_1 = atom[rest.atom_index_1];
	    mmdb::Atom *at_2 = atom[rest.atom_index_2];
	    mmdb::Atom *at_3 = atom[rest.atom_index_3];
	    if (at_c && at_1 && at_2 && at_3) {
	       std::cout << "   chiral volume problem centred at: "
			 << at_c->GetChainID() << " "
			 << at_c->GetSeqNum() << " "
			 << at_c->GetResName() << " "
			 << at_c->name 
			 << "  with neighbours "
			 << at_1->name << " "
			 << at_2->name << " "
			 << at_3->name << " "
			 << geometry_distortion[i].distortion_score
			 << std::endl;
	    }
	 }
	 double pen_score = geometry_distortion[i].distortion_score;
	 sum_penalties_chirals += pen_score;
      }

      if (rest.restraint_type == coot::PLANE_RESTRAINT) {
	 n_restraints_planes++;
	 std::vector<mmdb::Atom *> plane_atoms;
	 for (unsigned int iat=0; iat<rest.plane_atom_index.size(); iat++) { 
	    mmdb::Atom *at = atom[rest.plane_atom_index[iat].first];
	    if (at)
	       plane_atoms.push_back(at);
	 }
	 std::string penalty_string("plane ");
	 for (unsigned int iat=0; iat<plane_atoms.size(); iat++) {
	    penalty_string += std::string(plane_atoms[iat]->GetChainID()) + " ";
	    penalty_string += coot::util::int_to_string(plane_atoms[iat]->GetSeqNum()) + " ";
	    penalty_string += plane_atoms[iat]->name;
	    penalty_string += " ";
	 }
	 std::size_t len = penalty_string.length();
	 if (len < 88) // to match bonds
	    penalty_string += std::string(88-len, ' ');
	 double pen_score = geometry_distortion[i].distortion_score;
	 sum_penalties_planes += pen_score;
	 penalty_string += std::string(" penalty-score:  ") + coot::util::float_to_string(pen_score);
	 // std::cout << penalty_string << std::endl;
	 penalty_string_planes.push_back(std::pair<std::string, double > (penalty_string, pen_score));
      }
   }
	 
   std::sort(penalty_string_bonds.begin(),    penalty_string_bonds.end(),    coot::util::sd_compare);
   std::sort(penalty_string_angles.begin(),   penalty_string_angles.end(),   coot::util::sd_compare);
   std::sort(penalty_string_torsions.begin(), penalty_string_torsions.end(), coot::util::sd_compare);
   std::sort(penalty_string_planes.begin(),   penalty_string_planes.end(),   coot::util::sd_compare);

   std::reverse(penalty_string_bonds.begin(),    penalty_string_bonds.end());
   std::reverse(penalty_string_angles.begin(),   penalty_string_angles.end());
   std::reverse(penalty_string_torsions.begin(), penalty_string_torsions.end());
   std::reverse(penalty_string_planes.begin(),   penalty_string_planes.end());

   // sorted list, line by line
   for (unsigned int i=0; i<penalty_string_planes.size(); i++)
      std::cout << "   " << penalty_string_planes[i].first << std::endl;
   for (unsigned int i=0; i<penalty_string_bonds.size(); i++)
      std::cout << "   " << penalty_string_bonds[i].first << std::endl;
   for (unsigned int i=0; i<penalty_string_angles.size(); i++)
      std::cout << "   " << penalty_string_angles[i].first << std::endl;
   for (unsigned int i=0; i<penalty_string_torsions.size(); i++)
      std::cout << "   " << penalty_string_torsions[i].first << std::endl;

	 
   // Summary:
   double av_penalty_bond = 0;
   double av_penalty_angle = 0;
   double av_penalty_total = 0;
   if (n_restraints_bonds > 0)
      av_penalty_bond = sum_penalties_bonds/double(n_restraints_bonds);
   if (n_restraints_angles > 0)
      av_penalty_angle = sum_penalties_angles/double(n_restraints_angles);
   if ((n_restraints_bonds+n_restraints_angles) > 0) { 
      av_penalty_total = (sum_penalties_bonds+sum_penalties_angles)/(n_restraints_bonds+n_restraints_angles);
   }
   total_distortion =
      sum_penalties_bonds  +
      sum_penalties_angles +
      sum_penalties_torsions +
      sum_penalties_chirals  + 
      sum_penalties_planes;
	 
   std::cout << "Residue Distortion Summary: \n   "
	     << n_restraints_bonds  << " bond restraints\n   "
	     << n_restraints_angles << " angle restraints\n"
	     << "   sum of bond  distortions penalties:  " << sum_penalties_bonds  << "\n"
	     << "   sum of angle distortions penalties:  " << sum_penalties_angles << "\n"
	     << "   average bond  distortion penalty:    " << av_penalty_bond  << "\n"
	     << "   average angle distortion penalty:    " << av_penalty_angle << "\n"
	     << "   total distortion penalty:            " << total_distortion
	     << "\n"
	     << "   average distortion penalty:          " << av_penalty_total
	     << std::endl;

   return total_distortion;
}


double
coot::geometry_distortion_info_container_t::distortion() const {

   // why do some of these have their distortion set already, and others calculated now?

   double total_distortion = 0.0;
   for (unsigned int i=0; i< geometry_distortion.size(); i++) {
      const coot::simple_restraint &rest = geometry_distortion[i].restraint;
      if (rest.restraint_type == coot::BOND_RESTRAINT) {
	 mmdb::Atom *at_1 = atom[rest.atom_index_1];
	 mmdb::Atom *at_2 = atom[rest.atom_index_2];
	 if (at_1 && at_2) {
	    clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
	    double d = sqrt((p2-p1).lengthsq());
	    double distortion = d - rest.target_value;
	    double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
	    total_distortion += pen_score;
	 }
      }

      if (rest.restraint_type == coot::ANGLE_RESTRAINT) {
	 mmdb::Atom *at_1 = atom[rest.atom_index_1];
	 mmdb::Atom *at_2 = atom[rest.atom_index_2];
	 mmdb::Atom *at_3 = atom[rest.atom_index_3];
	 if (at_1 && at_2 && at_3) {
	    clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
	    clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
	    double angle_rad = clipper::Coord_orth::angle(p1, p2, p3);
	    double angle = clipper::Util::rad2d(angle_rad);
	    double distortion = angle - rest.target_value;
	    double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
	    total_distortion += pen_score;
	 }
      }

      if (rest.restraint_type == TORSION_RESTRAINT) {
	 mmdb::Atom *at_1 = atom[rest.atom_index_1];
	 mmdb::Atom *at_2 = atom[rest.atom_index_2];
	 mmdb::Atom *at_3 = atom[rest.atom_index_3];
	 mmdb::Atom *at_4 = atom[rest.atom_index_4];
	 if (at_1 && at_2 && at_3 && at_4) {
	    clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
	    clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
	    clipper::Coord_orth p4(at_4->x, at_4->y, at_4->z);
	    double torsion_rad = clipper::Coord_orth::torsion(p1, p2, p3, p4);
	    double torsion = clipper::Util::rad2d(torsion_rad);
	    double distortion = rest.torsion_distortion(torsion);
	    double pen_score = distortion*distortion/(rest.sigma*rest.sigma);
	    total_distortion += pen_score;
	 }
      }

      if (rest.restraint_type == TRANS_PEPTIDE_RESTRAINT) {
	 mmdb::Atom *at_1 = atom[rest.atom_index_1];
	 mmdb::Atom *at_2 = atom[rest.atom_index_2];
	 mmdb::Atom *at_3 = atom[rest.atom_index_3];
	 mmdb::Atom *at_4 = atom[rest.atom_index_4];
	 if (at_1 && at_2 && at_3 && at_4) {
	    clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
	    clipper::Coord_orth p3(at_3->x, at_3->y, at_3->z);
	    clipper::Coord_orth p4(at_4->x, at_4->y, at_4->z);
	    double torsion_rad = clipper::Coord_orth::torsion(p1, p2, p3, p4);
	    double torsion = clipper::Util::rad2d(torsion_rad);
	    double pen_score = rest.torsion_distortion(torsion);
	    total_distortion += pen_score;
	 }
      }

      if (rest.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	 total_distortion += geometry_distortion[i].distortion_score;
      }

      if (rest.restraint_type == coot::PLANE_RESTRAINT) {
	 total_distortion += geometry_distortion[i].distortion_score;
      }
   }
   return total_distortion;
}



coot::geometry_distortion_info_container_t
coot::restraints_container_t::geometric_distortions(coot::restraint_usage_Flags flags) {

   restraints_usage_flag = flags;
   setup_gsl_vector_variables();  //initial positions in x array
   coot::geometry_distortion_info_container_t dv = distortion_vector(x);
   return dv;
}

coot::geometry_distortion_info_container_t
coot::restraints_container_t::geometric_distortions() {

   // we don't want to do this if it has already been done. Hmmm.
   //
   // that's because this can be called when we are part-way through a refinement
   //
   if (!x)
      setup_gsl_vector_variables();  //initial positions in x array

   coot::geometry_distortion_info_container_t dv = distortion_vector(x);
   return dv;
} 


std::ostream&
coot::operator<<(std::ostream &s, geometry_distortion_info_container_t gdic) {

   s << "[ chain :" << gdic.chain_id << ": residues " << gdic.min_resno << " to " << gdic.max_resno
     << " residues: \n" ;
   for (unsigned int ires=0; ires<gdic.geometry_distortion.size(); ires++)
      s << "      " << gdic.geometry_distortion[ires] << "\n";
   s << "]\n";
   return s;
} 

std::ostream&
coot::operator<<(std::ostream &s, geometry_distortion_info_t gdi) {

   if (gdi.set) {
      s << gdi.restraint << " " << gdi.residue_spec << " distortion: " << gdi.distortion_score;
   } else {
      s << "{geometry_distortion_info-unset}";
   }
   return s;
}

coot::omega_distortion_info_container_t
coot::restraints_container_t::omega_trans_distortions(const coot::protein_geometry &geom,
						      bool mark_cis_peptides_as_bad_flag) {

   // restraints_usage_flag = flags;  // not used?
   setup_gsl_vector_variables();  //initial positions in x array
   std::string chain_id("");

   if (n_atoms > 0)
      chain_id = atom[0]->GetChainID();
   else
      chain_id = "blank"; // shouldn't happen.

   mmdb::Chain *chain_p = atom[0]->GetChain();
   
   // I think there will be need for some sort of scaling thing here.
   std::pair<short int, int> minr = coot::util::min_resno_in_chain(chain_p);
   std::pair<short int, int> maxr = coot::util::max_resno_in_chain(chain_p);
   int min_resno = 0;
   int max_resno = 0;
   if (minr.first && maxr.first) {
      min_resno = minr.second;
      max_resno = maxr.second;
   }
   double scale = 1.2; // how much to scale the omega difference by so that
		       // it looks good on the graph.

   if (!mark_cis_peptides_as_bad_flag)
      scale = 2.5;

   coot::omega_distortion_info_container_t dc(chain_id, min_resno, max_resno);

   // we need to pick out the CA and C of this and N and Ca of next.

   int n_chain_residues = chain_p->GetNumberOfResidues();
   for (int ires_serial=0; ires_serial<(n_chain_residues-1); ires_serial++) {
   
      mmdb::Residue *first  = chain_p->GetResidue(ires_serial);
      mmdb::Residue *second = chain_p->GetResidue(ires_serial+1);

      if (first && second) {

	 if (! chain_p->isSolventChain()) { 

	    std::pair<std::string, bool> lt = find_link_type_complicado(first, second, geom);

	    if (lt.first == "TRANS" || lt.first == "PTRANS" ||
		lt.first == "CIS"   || lt.first == "PCIS") {
	       
	       // So we have atoms selected in both residues, lets look for those atoms:

	       mmdb::Atom *at;
	       clipper::Coord_orth ca_first, c_first, n_next, ca_next;
	       short int got_ca_first = 0, got_c_first = 0, got_n_next = 0, got_ca_next = 0;
	       mmdb::PPAtom res_selection = NULL;
	       int i_no_res_atoms;

	       first->GetAtomTable(res_selection, i_no_res_atoms);
	       if (i_no_res_atoms > 0) {
		  for (int iresatom=0; iresatom<i_no_res_atoms; iresatom++) {
		     at = res_selection[iresatom];
		     std::string atom_name(at->name);
		     if (atom_name == " CA ") {
			ca_first = clipper::Coord_orth(at->x, at->y, at->z);
			got_ca_first = 1;
		     }
		     if (atom_name == " C  ") {
			c_first = clipper::Coord_orth(at->x, at->y, at->z);
			got_c_first = 1;
		     }
		  }
	       }
	       second->GetAtomTable(res_selection, i_no_res_atoms);
	       if (i_no_res_atoms > 0) {
		  for (int iresatom=0; iresatom<i_no_res_atoms; iresatom++) {
		     at = res_selection[iresatom];
		     std::string atom_name(at->name);
		     if (atom_name == " CA ") {
			ca_next = clipper::Coord_orth(at->x, at->y, at->z);
			got_ca_next = 1;
		     }
		     if (atom_name == " N  ") {
			n_next = clipper::Coord_orth(at->x, at->y, at->z);
			got_n_next = 1;
		     }
		  }
	       }

	       if (got_ca_first && got_c_first && got_n_next && got_ca_next) {
		  // the omega angle belongs to the second residue
		  double tors = clipper::Coord_orth::torsion(ca_first, c_first, n_next, ca_next);
		  double torsion = clipper::Util::rad2d(tors);
		  torsion = (torsion > 0.0) ? torsion : 360.0 + torsion;
		  std::string info = chain_id;
		  info += " ";
		  info += coot::util::int_to_string(second->GetSeqNum());
		  info += " ";
		  info += second->name;
		  info += " Omega: ";
		  info += coot::util::float_to_string(torsion);
		  double distortion = fabs(180.0 - torsion);

		  // Add comment on Cis peptide, if it is:
		  if (torsion < 90.0 || torsion > 270.0) 
		     info += "  (Cis)";

		  if (!mark_cis_peptides_as_bad_flag)
		     // consider the cases: torsion=1 and torsion=359
		     if (distortion > 90.0) {
			distortion = torsion;
			if (distortion > 180.0) {
			   distortion -= 360;
			   distortion = fabs(distortion);
			}
		     }
		  distortion *= scale;
		  omega_distortion_info_t odi(second->GetSeqNum(), distortion, info);
		  dc.omega_distortions.push_back(odi);
	       } else {
		  std::cout << "INFO:: failed to get all atoms for omega torsion "
			    << "chain " << chain_id << " residues "
			    << first->GetSeqNum() << " to " << second->GetSeqNum()
			    << std::endl;
	       }
	    }
	 }
      }
   }
   return dc;
}

// return value in distortion
void
coot::distortion_score_single_thread(const gsl_vector *v, void *params,
				     int idx_start, int idx_end, double *distortion) {
   // first extract the object from params
   //
   coot::restraints_container_t *restraints = static_cast<coot::restraints_container_t *>(params);

   double d = 0;
   for (int i=idx_start; i<idx_end; i++) {

      if (restraints->restraints_usage_flag & coot::NON_BONDED_MASK) { // 16:
	 if ( (*restraints)[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	    d = coot::distortion_score_non_bonded_contact( (*restraints)[i], v);
	    // std::cout << "dsm: nbc  single-thread " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::BONDS_MASK) { // 1: bonds
	 if ( (*restraints)[i].restraint_type == coot::BOND_RESTRAINT) {
	    d = coot::distortion_score_bond((*restraints)[i], v);
	    // std::cout << "dsm: bond  single-thread " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::ANGLES_MASK) { // 2: angles
	 if ( (*restraints)[i].restraint_type == coot::ANGLE_RESTRAINT) {
	    d = coot::distortion_score_angle((*restraints)[i], v);
	    // std::cout << "dsm: angle single-thread " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & TRANS_PEPTIDE_MASK) {
	 if ( (*restraints)[i].restraint_type == TRANS_PEPTIDE_RESTRAINT) {
	    double d =  coot::distortion_score_trans_peptide(restraints->at(i), v);
	    // std::cout << "dsm: trans_peptide single-thread " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::TORSIONS_MASK) { // 4: torsions
	 if ( (*restraints)[i].restraint_type == coot::TORSION_RESTRAINT) {
	    try {
	       double d =  coot::distortion_score_torsion(i, restraints->at(i), v);
	       // std::cout << "dsm: torsion single-thread " << d << std::endl;
	       *distortion += d;
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "ERROR::" << rte.what() << std::endl;
	    }
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::PLANES_MASK) { // 8: planes
	 if ( (*restraints)[i].restraint_type == coot::PLANE_RESTRAINT) {
	    d =  coot::distortion_score_plane((*restraints)[i], v);
	    // std::cout << "dsm: plane single-thread " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::PARALLEL_PLANES_MASK) { // 128
	 if ( (*restraints)[i].restraint_type == coot::PARALLEL_PLANES_RESTRAINT) {
	    d =  coot::distortion_score_parallel_planes((*restraints)[i], v);
	    // std::cout << "dsm: paralelplane single-thread " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) {
   	 if ( (*restraints)[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
   	    d = coot::distortion_score_chiral_volume( (*restraints)[i], v);
   	    *distortion += d;
	    // std::cout << "dsm: chiral single-trhead " << d << std::endl;
	    continue;
   	 }
      }

      if (restraints->restraints_usage_flag & coot::RAMA_PLOT_MASK) {
   	 if ( (*restraints)[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) {
	    if (restraints->rama_type == restraints_container_t::RAMA_TYPE_ZO) {
	       d = coot::distortion_score_rama( (*restraints)[i], v, restraints->ZO_Rama(), restraints->get_rama_plot_weight());
	    } else {
	       d = coot::distortion_score_rama( (*restraints)[i], v, restraints->LogRama());
	    }
	    // std::cout << "dsm: rama single-thread " << d << std::endl;
   	    *distortion += d; // positive is bad...  negative is good.
	    continue;
   	 }
      }

      if (restraints->restraints_usage_flag & coot::START_POS_RESTRAINT_MASK) {
	 if ( (*restraints)[i].restraint_type == coot::START_POS_RESTRAINT) {
	    *distortion = coot::distortion_score_start_pos((*restraints)[i], params, v);
	    // std::cout << "dsm: start-pos single-thread " << d << std::endl;
	 }
      }

      if (restraints->restraints_usage_flag & coot::GEMAN_MCCLURE_DISTANCE_MASK) {
	 if ( (*restraints)[i].restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	    d = coot::distortion_score_geman_mcclure_distance((*restraints)[i], v,
	                                                        restraints->geman_mcclure_alpha);
	    *distortion += d;
	    // std::cout << "dsm: geman-mcclure single-thread idx " << i << " " << d << std::endl;
	 }
      }
   }
}

#ifdef HAVE_CXX_THREAD
// return value in distortion
void
coot::distortion_score_multithread(int thread_id, const gsl_vector *v, void *params,
				   int idx_start, int idx_end, double *distortion,
				   std::atomic<unsigned int> &done_count_for_threads) {

   // first extract the object from params
   //
   coot::restraints_container_t *restraints = static_cast<coot::restraints_container_t *>(params);

   double d = 0;
   for (int i=idx_start; i<idx_end; i++) {

      if (restraints->restraints_usage_flag & coot::NON_BONDED_MASK) { // 16:
	 if ( (*restraints)[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	    d = coot::distortion_score_non_bonded_contact( (*restraints)[i], v);
	    // std::cout << "dsm: nbc  thread_idx " << thread_id << " idx " << i << " " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::BONDS_MASK) { // 1: bonds
	 if ( (*restraints)[i].restraint_type == coot::BOND_RESTRAINT) {
	    d = coot::distortion_score_bond((*restraints)[i], v);
	    // std::cout << "dsm: bond  thread_idx " << thread_id << " idx " << i << " " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::ANGLES_MASK) { // 2: angles
	 if ( (*restraints)[i].restraint_type == coot::ANGLE_RESTRAINT) {
	    d = coot::distortion_score_angle((*restraints)[i], v);
	    // std::cout << "dsm: angle thread_idx " << thread_id << " idx " << i << " " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & TRANS_PEPTIDE_MASK) {
	 if ( (*restraints)[i].restraint_type == TRANS_PEPTIDE_RESTRAINT) {
	       double d =  coot::distortion_score_trans_peptide((*restraints)[i], v);
	    *distortion += d;
	    // std::cout << "dsm: trans-peptide " << thread_id << " idx " << i << " " << d << std::endl;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::TORSIONS_MASK) { // 4: torsions
	 if ( (*restraints)[i].restraint_type == coot::TORSION_RESTRAINT) {
	    double d =  coot::distortion_score_torsion(i, restraints->at(i), v);
	    // std::cout << "dsm: torsion " << thread_id << " idx " << i << " " << d << std::endl;
	    *distortion += d;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::PLANES_MASK) { // 8: planes
	 if ( (*restraints)[i].restraint_type == coot::PLANE_RESTRAINT) {
	    d =  coot::distortion_score_plane((*restraints)[i], v);
	    *distortion += d;
	    // std::cout << "dsm: plane " << thread_id << " idx " << i << " " << d << std::endl;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::PARALLEL_PLANES_MASK) { // 128
	 if ( (*restraints)[i].restraint_type == coot::PARALLEL_PLANES_RESTRAINT) {
	    d =  coot::distortion_score_parallel_planes((*restraints)[i], v);
	    *distortion += d;
	    // std::cout << "dsm: parallel plane " << thread_id << " idx " << i << " " << d << std::endl;
	    continue;
	 }
      }

      if (restraints->restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) { 
   	 if ( (*restraints)[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) { 
   	    d = coot::distortion_score_chiral_volume( (*restraints)[i], v);
   	    *distortion += d;
	    // std::cout << "dsm: chiral " << thread_id << " idx " << i << " " << d << std::endl;
	    continue;
   	 }
      }

      if (restraints->restraints_usage_flag & coot::RAMA_PLOT_MASK) {
   	 if ( (*restraints)[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) {
	    // std::cout << "dsm: rama " << thread_id << " idx " << i << " " << d << std::endl;
	    if (restraints->rama_type == restraints_container_t::RAMA_TYPE_ZO) {
	       d = coot::distortion_score_rama( (*restraints)[i], v, restraints->ZO_Rama(), restraints->get_rama_plot_weight());
	    } else {
	       d = coot::distortion_score_rama( (*restraints)[i], v, restraints->LogRama());
	    }
   	    *distortion += d; // positive is bad...  negative is good.
	    continue;
   	 }
      }

      if (restraints->restraints_usage_flag & coot::GEMAN_MCCLURE_DISTANCE_MASK) {
   	 if ( (*restraints)[i].restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	    d = coot::distortion_score_geman_mcclure_distance((*restraints)[i], v,
							      restraints->geman_mcclure_alpha);
	    *distortion += d;
	    // std::cout << "dsm: geman-mclure " << thread_id << " idx " << i << " " << d << std::endl;
	 }
      }

      if (restraints->restraints_usage_flag & coot::START_POS_RESTRAINT_MASK) {
	 if ( (*restraints)[i].restraint_type == coot::START_POS_RESTRAINT) {
	    d = coot::distortion_score_start_pos((*restraints)[i], params, v);
	    *distortion += d;
	    // std::cout << "dsm: start-pos " << thread_id << " idx " << i << " " << d << std::endl;
	 }
      }
   }
   done_count_for_threads++; // atomic operation
}
#endif // HAVE_CXX_THREAD

// Return the distortion score.
//
double coot::distortion_score(const gsl_vector *v, void *params) {

#ifdef ANALYSE_REFINEMENT_TIMING
#endif // ANALYSE_REFINEMENT_TIMING

   // so we are comparing the geometry of the value in the gsl_vector
   // v and the ideal values.
   //

   // first extract the object from params
   //
   coot::restraints_container_t *restraints_p = static_cast<coot::restraints_container_t *>(params);

   double distortion = 0;

   // distortion += starting_structure_diff_score(v, params);

   // tmp debugging stuff
   double nbc_diff = 0.0;
   double d;
   int restraints_size = restraints_p->size();

#ifdef HAVE_CXX_THREAD

   if (false) // debug
      std::cout << "here in coot::distortion_score() thread_pool is "
		<< restraints_p->thread_pool_p << std::endl;

   if (restraints_p->thread_pool_p) {

      int n_per_thread = restraints_size/restraints_p->n_threads;
      std::atomic<unsigned int> done_count_for_threads(0);

      double distortions[restraints_p->n_threads];
      for (unsigned int i_thread=0; i_thread<restraints_p->n_threads; i_thread++) {
	 int idx_start = i_thread * n_per_thread;
	 int idx_end   = idx_start + n_per_thread;
	 // for the last thread, set the end restraint index
	 if (i_thread == (restraints_p->n_threads - 1))
	    idx_end = restraints_size; // for loop uses iat_start and tests for < iat_end
	 distortions[i_thread] = 0;

	 if (false)
	    std::cout << "pushing thread doing index range " << idx_start << " " << idx_end
		      << std::endl;
	 restraints_p->thread_pool_p->push(distortion_score_multithread,
					   v, params, idx_start, idx_end, &distortions[i_thread],
					   std::ref(done_count_for_threads));
	 if (false)
	    std::cout << "    now thread pool info: size: "
		      << restraints_p->thread_pool_p->size() << " idle: "
		      << restraints_p->thread_pool_p->n_idle() 
		      << std::endl;
      }

      while (done_count_for_threads != restraints_p->n_threads) {
	 std::this_thread::sleep_for(std::chrono::microseconds(1));
      }
	 
      for (unsigned int i_thread=0; i_thread<restraints_p->n_threads; i_thread++)
	 distortion += distortions[i_thread];

      // double distortion_single;
      // distortion_score_multithread(999, v, params, 0, restraints_size, &distortion_single);
      // std::cout << "distortion  multi: " << distortion << std::endl;
      // std::cout << "distortion single: " << distortion_single << std::endl;

   } else {
      // "return" value is passed pointer
      distortion_score_single_thread(v, params, 0, restraints_size, &distortion);
   }

#else

   // "return" value is passed pointer
   distortion_score_single_thread(v, params, 0, restraints_size, &distortion);

#endif // HAVE_CXX_THREAD

//     std::cout << "nbc_diff   distortion: " << nbc_diff << std::endl;
//     std::cout << "post-terms distortion: " << distortion << std::endl;

   if (restraints_p->include_map_terms())
      // multi-thread this too:
      distortion += coot::electron_density_score(v, params); // good map fit: low score
   
#ifdef ANALYSE_REFINEMENT_TIMING
#endif // ANALYSE_REFINEMENT_TIMING

   return distortion; 
}

coot::geometry_distortion_info_container_t
coot::restraints_container_t::distortion_vector(const gsl_vector *v) const {

   std::string chainid;
   if (n_atoms > 0)
      chainid = atom[0]->GetChainID();
   else
      chainid = "blank";

   coot::geometry_distortion_info_container_t distortion_vec_container(atom, n_atoms, chainid);
   double distortion = 0.0;

   for (unsigned int i=0; i<restraints_vec.size(); i++) {
      int atom_index = -1; // initial unset, the atom index used to spec the residue.
      const simple_restraint &rest = restraints_vec[i];
      std::vector<int> atom_indices;
      if (restraints_usage_flag & coot::BONDS_MASK) 
	 if (restraints_vec[i].restraint_type == coot::BOND_RESTRAINT) { 
	    distortion = coot::distortion_score_bond(rest, v);
	    atom_index = restraints_vec[i].atom_index_1;
	    atom_indices.push_back(rest.atom_index_1);
	    atom_indices.push_back(rest.atom_index_2);
	 } 

      if (restraints_usage_flag & coot::ANGLES_MASK) 
	 if (restraints_vec[i].restraint_type == coot::ANGLE_RESTRAINT) { 
	    distortion = coot::distortion_score_angle(rest, v);
	    atom_index = restraints_vec[i].atom_index_1;
	    atom_indices.push_back(rest.atom_index_1);
	    atom_indices.push_back(rest.atom_index_2);
	    atom_indices.push_back(rest.atom_index_3);
	 }

      if (restraints_usage_flag & coot::TORSIONS_MASK) {
	 if (restraints_vec[i].restraint_type == coot::TORSION_RESTRAINT) {
	    // distortion_score_torsion can throw a std::runtime_error
	    try {
	       distortion = coot::distortion_score_torsion(i, rest, v);
	       atom_index = restraints_vec[i].atom_index_1;
	       atom_indices.push_back(rest.atom_index_1);
	       atom_indices.push_back(rest.atom_index_2);
	       atom_indices.push_back(rest.atom_index_3);
	       atom_indices.push_back(rest.atom_index_4);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "ERROR::" << rte.what() << std::endl;
	    }
	 }
      }

      if (restraints_usage_flag & coot::PLANES_MASK) 
	 if (restraints_vec[i].restraint_type == coot::PLANE_RESTRAINT) { 
	    distortion = coot::distortion_score_plane(restraints_vec[i], v);
	    atom_index = restraints_vec[i].plane_atom_index[0].first;
	    for (unsigned int ii=0; ii<rest.plane_atom_index.size(); ii++) { 
	       if (rest.plane_atom_index[ii].first >= 0)
		  atom_indices.push_back(rest.plane_atom_index[ii].first);
	    }
	 }

      if (restraints_usage_flag & coot::PARALLEL_PLANES_MASK) 
	 if (restraints_vec[i].restraint_type == coot::PARALLEL_PLANES_RESTRAINT) { 
	    distortion = coot::distortion_score_parallel_planes(restraints_vec[i], v);
	    atom_index = restraints_vec[i].plane_atom_index[0].first;
	 } 
      if (restraints_usage_flag & coot::NON_BONDED_MASK)
	 if (restraints_vec[i].restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	    distortion = coot::distortion_score_non_bonded_contact(restraints_vec[i], v);
	    atom_index = restraints_vec[i].atom_index_1;
	    atom_indices.push_back(rest.atom_index_1);
	    atom_indices.push_back(rest.atom_index_2);
	    // std::cout << " NBC i " << i << " " << distortion << std::endl;
	 }
      if (restraints_usage_flag & coot::GEMAN_MCCLURE_DISTANCE_MASK)
	 if (restraints_vec[i].restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	    distortion = coot::distortion_score_geman_mcclure_distance(restraints_vec[i], v,
								       geman_mcclure_alpha);
	    atom_index = restraints_vec[i].atom_index_1;
	    atom_indices.push_back(rest.atom_index_1);
	    atom_indices.push_back(rest.atom_index_2);
	 }

      if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK)
   	 if (restraints_vec[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
   	    distortion = coot::distortion_score_chiral_volume(restraints_vec[i], v);
	    atom_index = restraints_vec[i].atom_index_centre;
	    atom_indices.push_back(rest.atom_index_centre);
	    atom_indices.push_back(rest.atom_index_1);
	    atom_indices.push_back(rest.atom_index_2);
	    atom_indices.push_back(rest.atom_index_3);
	 } 

      if (restraints_usage_flag & coot::RAMA_PLOT_MASK) 
    	 if (restraints_vec[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) { 
	    if (rama_type == restraints_container_t::RAMA_TYPE_ZO) {
	       distortion = coot::distortion_score_rama(restraints_vec[i], v, ZO_Rama(), get_rama_plot_weight());
	    } else {
	       distortion = coot::distortion_score_rama(restraints_vec[i], v, lograma);
	    }
	    atom_index = restraints_vec[i].atom_index_1;
	    atom_indices.push_back(rest.atom_index_1);
	    atom_indices.push_back(rest.atom_index_2);
	    atom_indices.push_back(rest.atom_index_3);
	    atom_indices.push_back(rest.atom_index_4);
	    atom_indices.push_back(rest.atom_index_5);
	 }

      if (atom_index != -1) {
	 coot::residue_spec_t rs(atom[atom_index]->GetResidue());
	 coot::geometry_distortion_info_t gdi(distortion, rest, rs);
	 gdi.atom_indices = atom_indices;
	 distortion_vec_container.geometry_distortion.push_back(gdi);
      }
   }


   // Find max_resno and min_resno.
   // 
   // Notice that we only use BOND_RESTRAINT to do this (atom_index_1 is
   // not defined in plane restraint).
   int max_resno =  -9999999;
   int min_resno = 999999999;

   for (unsigned int i=0; i<distortion_vec_container.geometry_distortion.size(); i++) {
      if (false)
	 std::cout << "distortion_vector() restraint " << i << " of "
		   << distortion_vec_container.geometry_distortion.size() << " "
		   << restraints_vec[i].restraint_type << std::endl;
      if (restraints_usage_flag & coot::BONDS_MASK)
	 if (restraints_vec[i].restraint_type == coot::BOND_RESTRAINT) {
	   const int &idx1 = distortion_vec_container.geometry_distortion[i].atom_indices[0];
	   const int &idx2 = distortion_vec_container.geometry_distortion[i].atom_indices[1];

	   // std::cout << "idx1 " << idx1 << std::endl;
	   // std::cout << "idx2 " << idx2 << std::endl;

	   int this_resno1 = distortion_vec_container.atom[idx1]->GetSeqNum();
	   int this_resno2 = distortion_vec_container.atom[idx2]->GetSeqNum();
	   if (this_resno1 < min_resno)
	      min_resno = this_resno1;
	   if (this_resno2 < min_resno)
	      min_resno = this_resno2;
	   if (this_resno1 > max_resno)
	      max_resno = this_resno1;
	   if (this_resno2 > max_resno)
	      max_resno = this_resno2;
	 }
   }
   distortion_vec_container.set_min_max(min_resno, max_resno);

   return distortion_vec_container;
}


double
coot::distortion_score_bond(const coot::simple_restraint &bond_restraint,
			    const gsl_vector *v) {

   int idx = 3*(bond_restraint.atom_index_1 - 0); 
   clipper::Coord_orth a1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(bond_restraint.atom_index_2 - 0); 
   clipper::Coord_orth a2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   
   double weight = 1.0/(bond_restraint.sigma * bond_restraint.sigma);
   double bit = (clipper::Coord_orth::length(a1,a2) - bond_restraint.target_value);
   
   return weight * bit *bit;
}

double
coot::distortion_score_geman_mcclure_distance(const coot::simple_restraint &restraint,
					      const gsl_vector *v,
					      const double &alpha) {

   int idx = 3*restraint.atom_index_1;
   clipper::Coord_orth a1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*restraint.atom_index_2;
   clipper::Coord_orth a2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   
   // Let z = (boi - bi)/sigma
   // so S_i = z^2/(1 + alpha * z^2)
   //
   double bl = clipper::Coord_orth::length(a1,a2);
   double bit = bl - restraint.target_value;
   double z = bit/restraint.sigma;
   double distortion = z*z/(1+alpha*z*z);

   if (false)
      std::cout << "distortion_score_geman_mcclure_distance: " << bl
		<< " sigma " << restraint.sigma
		<< " target " << restraint.target_value
		<< " alpha " << alpha
		<< " z " << z
		<< " distortion " << distortion << std::endl;

   // return z * z; // least squares
   return distortion;
}

double
coot::distortion_score_angle(const coot::simple_restraint &angle_restraint,
			     const gsl_vector *v) {


   int idx = 3*(angle_restraint.atom_index_1); 
   clipper::Coord_orth a1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(angle_restraint.atom_index_2); 
   clipper::Coord_orth a2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(angle_restraint.atom_index_3); 
   clipper::Coord_orth a3(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   clipper::Coord_orth d1 = a1 - a2; 
   clipper::Coord_orth d2 = a3 - a2; 
   double len1 = clipper::Coord_orth::length(a1,a2);
   double len2 = clipper::Coord_orth::length(a3,a2);

   // len1 = len1 > 0.01 ? len1 : 0.01; 
   // len2 = len2 > 0.01 ? len2 : 0.01;
   if (len1 < 0.01) {
      len1 = 0.01;
      d1 = clipper::Coord_orth(0.01, 0.01, 0.01);
   } 
   if (len2 < 0.01) {
      len2 = 0.01;
      d2 = clipper::Coord_orth(0.01, 0.01, 0.01);
   } 

   double cos_theta = clipper::Coord_orth::dot(d1,d2)/(len1*len2);
   if (cos_theta < -1) cos_theta = -1.0;
   if (cos_theta >  1) cos_theta =  1.0;
   double theta = acos(cos_theta);
   double bit = clipper::Util::rad2d(theta) - angle_restraint.target_value;
   double weight = 1/(angle_restraint.sigma * angle_restraint.sigma);
   if (0)
      std::cout << "actual: " << clipper::Util::rad2d(theta)
		<< " cos_theta " << cos_theta
		<< " target: "
		<< angle_restraint.target_value
		<< " adding distortion score: "
		<< weight * pow(bit, 2.0) << std::endl;
   return weight * bit * bit;
}


//
// Return the distortion score from a single torsion restraint.
// 
// can throw a std::runtime_error if there is a problem calculating the torsion.
// 
double
coot::distortion_score_torsion(unsigned int idx_restraint,
			       const coot::simple_restraint &torsion_restraint,
			       const gsl_vector *v) {

   // First calculate the torsion:
   // theta = arctan(E/G);
   // where E = a.(bxc) and G = -a.c + (a.b)(b.c)

   int idx; 

   idx = 3*(torsion_restraint.atom_index_1);
   clipper::Coord_orth P1(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));
   idx = 3*(torsion_restraint.atom_index_2);
   clipper::Coord_orth P2(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));
   idx = 3*(torsion_restraint.atom_index_3);
   clipper::Coord_orth P3(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));
   idx = 3*(torsion_restraint.atom_index_4);
   clipper::Coord_orth P4(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));

//    P1 = clipper::Coord_orth( 1.0,  0.0, 1.0);
//    P2 = clipper::Coord_orth( 0.0, -1.0, 1.0);
//    P3 = clipper::Coord_orth( 0.0,  0.0, 0.0);
//    P4 = clipper::Coord_orth(-1.0, -1.0, 1.0);
//    P4 = clipper::Coord_orth( 1.0,  1.0, 1.0);

   clipper::Coord_orth a = P2 - P1; 
   clipper::Coord_orth b = P3 - P2; 
   clipper::Coord_orth c = P4 - P3;

   // b*b * [ a.(bxc)/b ]
   double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
      sqrt( b.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
      + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

   double theta = clipper::Util::rad2d(atan2(E,G));
//    double clipper_theta = 
//       clipper::Util::rad2d(clipper::Coord_orth::torsion(P1, P2, P3, P4));

   if (clipper::Util::isnan(theta)) {
      std::string mess = "WARNING: distortion_score_torsion() observed torsion theta is a NAN!";
      throw std::runtime_error(mess);
   }

   // instabilty when the P2-P3-P4 or P1-P2-p3 line is linear. Give up with the derivatives
   // similar escape in the derivatives
   double al = sqrt(clipper::Coord_orth::dot(a,a));
   double bl = sqrt(clipper::Coord_orth::dot(b,b));
   double cl = sqrt(clipper::Coord_orth::dot(c,c));
   double cos_a1 = clipper::Coord_orth::dot(a,b)/(al*bl);
   double cos_a2 = clipper::Coord_orth::dot(b,c)/(bl*cl);

   // I don't like this escape at all
   //
   if (cos_a1 > 0.9 || cos_a2> 0.9) {
      return 0;
   }

   if (theta < 0.0) theta += 360.0;

   //if (torsion_restraint.periodicity == 1) {
   // }
   // 
   //    double diff = theta - torsion_restraint.target_value;

   double diff = 99999.9; 
   double tdiff; 
   double trial_target; 
   int per = torsion_restraint.periodicity;
   for(int i=0; i<per; i++) { 
      // trial_target = torsion_restraint.target_value + double(i)*360.0/double(per);  ??
      trial_target = torsion_restraint.target_value + double(i)*360.0/double(per); 
      if (trial_target >= 360.0) trial_target -= 360.0; 
      tdiff = theta - trial_target;
      if (tdiff < -180) tdiff += 360;
      if (tdiff >  180) tdiff -= 360;
      if (fabs(tdiff) < fabs(diff)) {
	 diff = tdiff;
      }
   }

   if (false)
      std::cout << "DEBUG:: atom index: " << torsion_restraint.atom_index_1
		<< " target " << torsion_restraint.target_value
		<< ", theta: " << theta << " \tdiff: " << diff << std::endl;

   if (diff >= 99999.0) { 
      std::cout << "Error in periodicity (" << per << ") check" << std::endl;
      std::cout << "target_value: " << torsion_restraint.target_value
		<< ", theta: " << theta << std::endl;
   }
   if (diff < -180.0) { 
      diff += 360.; 
   } else { 
      if (diff > 180.0) { 
	 diff -= 360.0; 
      }
   }


   if (false) { // debug
      double pen = diff*diff/(torsion_restraint.sigma * torsion_restraint.sigma);
      std::cout << "distortion_torsion " << idx_restraint << " theta (calc): " << theta
		<< " periodicity " << torsion_restraint.periodicity
		<< " target "      << torsion_restraint.target_value
		<< " diff: " << diff << " distortion " << pen << std::endl ;
      std::cout << "   in distortion_torsion: sigma = " << torsion_restraint.sigma
		<< ", weight=" << pow(torsion_restraint.sigma,-2.0)
		<< " and diff is " << diff << std::endl;
   }
   
   return diff*diff/(torsion_restraint.sigma * torsion_restraint.sigma);

}

double
coot::distortion_score_plane(const coot::simple_restraint &plane_restraint,
			     const gsl_vector *v) {

   coot::plane_distortion_info_t info =
      distortion_score_plane_internal(plane_restraint, v);

   return info.distortion_score;

}

double
coot::distortion_score_parallel_planes(const simple_restraint &ppr,
				       const gsl_vector *v) {

   double score = 0;
	 
   plane_distortion_info_t info =
      distortion_score_2_planes(ppr.plane_atom_index, ppr.atom_index_other_plane, ppr.sigma, v);
   if (0)
      std::cout << "parallel plane-combined abcd "
		<< info.abcd[0] << " "
		<< info.abcd[1] << " "
		<< info.abcd[2] << " "
		<< info.abcd[3] << std::endl;
   return info.distortion_score;
} 


double 
coot::distortion_score_chiral_volume(const coot::simple_restraint &chiral_restraint,
				     const gsl_vector *v) {

   
   int idx = 3*(chiral_restraint.atom_index_centre);
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   idx = 3*(chiral_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));

   clipper::Coord_orth a = a1 - centre;
   clipper::Coord_orth b = a2 - centre;
   clipper::Coord_orth c = a3 - centre;

   double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
   // double volume_sign = chiral_restraint.chiral_volume_sign;
   double distortion = cv - chiral_restraint.target_chiral_volume;

   if (0) {
      std::cout << "atom indices: "
		<< chiral_restraint.atom_index_centre << " "
		<< chiral_restraint.atom_index_1 << " " 
		<< chiral_restraint.atom_index_2 << " " 
		<< chiral_restraint.atom_index_3 << " ";
      std::cout << "DEBUG:: (distortion) chiral volume target "
		<< chiral_restraint.target_chiral_volume
		<< " chiral actual " << cv << " diff: " << distortion;
   }
 	     

   distortion *= distortion;
   distortion /= chiral_restraint.sigma * chiral_restraint.sigma;

   if (0)
      std::cout << " score: " << distortion << "\n";

   return distortion;
}

double
coot::distortion_score_rama(const coot::simple_restraint &rama_restraint,
			    const gsl_vector *v,
			    const LogRamachandran &lograma) {

   double distortion = 0;
   // First calculate the torsions:
   // theta = arctan(E/G); 
   // where E = a.(bxc) and G = -a.c + (a.b)(b.c)

   int idx; 

   idx = 3*(rama_restraint.atom_index_1);
   clipper::Coord_orth P1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_2); 
   clipper::Coord_orth P2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_3); 
   clipper::Coord_orth P3(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_4); 
   clipper::Coord_orth P4(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_5); 
   clipper::Coord_orth P5(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));

//    P1 = clipper::Coord_orth(1.0, 0.0, 1.0);
//    P2 = clipper::Coord_orth(0.0, -1.0, 1.0);
//    P3 = clipper::Coord_orth(0.0, 0.0, 0.0);
//    P4 = clipper::Coord_orth(-1.0, -1.0, 1.0);
//    P4 = clipper::Coord_orth(1.0, 1.0, 1.0);

   clipper::Coord_orth a = P2 - P1; 
   clipper::Coord_orth b = P3 - P2; 
   clipper::Coord_orth c = P4 - P3;
   clipper::Coord_orth d = P5 - P4;

   // Old (6 atom) wrong:
   // TRANS    psi      1 N (P1)    1 CA (P2)     1 C  (P3)    2 N (P4)
   // TRANS    phi      1 C (P3)    2  N (P4)     2 CA (P5)    2 C (P6)
   //
   // New assignements:
   // TRANS    phi    (1st C) (2nd N ) (2nd CA) (2nd C)
   // TRANS    psi    (2nd N) (2nd CA) (2nd C ) (3nd N)
   //
   // So Rama_atoms in this order:
   //   0       1        2      3         4
   //  P1      P2       P3     P4        P5
   // (1st C) (2nd N) (2nd CA) (2nd C) (3rd N)

   // ---------- phi ------------------
   // b*b * [ a.(bxc)/b ]
   double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
      sqrt( b.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
      + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

   double phi = clipper::Util::rad2d(atan2(E,G));
   if (phi < 180.0)
      phi += 360.0;
   if (phi > 180.0)
      phi -= 360.0;

   // ---------- psi ------------------
   // b*b * [ a.(bxc)/b ]
   double H = clipper::Coord_orth::dot(b, clipper::Coord_orth::cross(c,d)) *
      sqrt( c.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double I = - clipper::Coord_orth::dot(b,d)*c.lengthsq()
      + clipper::Coord_orth::dot(b,c)*clipper::Coord_orth::dot(c,d);

   double psi = clipper::Util::rad2d(atan2(H,I));
   if (psi < 180.0)
      psi += 360.0;
   if (psi > 180.0)
      psi -= 360.0;

   double lr = lograma.interp(clipper::Util::d2rad(phi), clipper::Util::d2rad(psi));
   double R = 10.0 * lr;
   // std::cout << "rama (lograma) distortion for " << phi << " " << psi << " is " << R << std::endl;

   if ( clipper::Util::isnan(phi) ) {
      std::cout << "WARNING: observed torsion phi is a NAN!" << std::endl;
      std::cout << "         debug-info: " << E << "/" << G << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_1 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_2 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_3 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_4 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_5 << std::endl;
      std::cout << "         debug-info: P1: " << P1.format() << std::endl;
      std::cout << "         debug-info: P2: " << P2.format() << std::endl;
      std::cout << "         debug-info: P3: " << P3.format() << std::endl;
      std::cout << "         debug-info: P4: " << P4.format() << std::endl;
      std::cout << "         debug-info: P5: " << P5.format() << std::endl;
      std::cout << "         debug-info: a: " << a.format() << std::endl;
      std::cout << "         debug-info: b: " << b.format() << std::endl;
      std::cout << "         debug-info: c: " << c.format() << std::endl;
      std::cout << "         debug-info: d: " << d.format() << std::endl;

      for (unsigned int i=0; i<15; i++) {
	 std::cout << "           in distortion_score_rama() " << i << " "
		   << gsl_vector_get(v, 3*i  ) << " "
		   << gsl_vector_get(v, 3*i+1) << " "
		   << gsl_vector_get(v, 3*i+2) << " " << std::endl;
      }
   }
   if ( clipper::Util::isnan(psi) ) {
      std::cout << "WARNING: observed torsion psi is a NAN!" << std::endl;
      std::cout << "         debug-info: " << H << "/" << I << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_1 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_2 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_3 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_4 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_5 << std::endl;
      std::cout << "         debug-info: P1: " << P1.format() << std::endl;
      std::cout << "         debug-info: P2: " << P2.format() << std::endl;
      std::cout << "         debug-info: P3: " << P3.format() << std::endl;
      std::cout << "         debug-info: P4: " << P4.format() << std::endl;
      std::cout << "         debug-info: P5: " << P5.format() << std::endl;
      std::cout << "         debug-info: a: " << a.format() << std::endl;
      std::cout << "         debug-info: b: " << b.format() << std::endl;
      std::cout << "         debug-info: c: " << c.format() << std::endl;
      std::cout << "         debug-info: d: " << d.format() << std::endl;
   }

   return R;
}

double 
coot::distortion_score_rama(const coot::simple_restraint &rama_restraint,
			    const gsl_vector *v,
			    const zo::rama_table_set &rama,
			    float rama_plot_weight) {
			    // const LogRamachandran &lograma) { // debugging

   double distortion = 0;
   // First calculate the torsions:
   // theta = arctan(E/G);
   // where E = a.(bxc) and G = -a.c + (a.b)(b.c)

   int idx;

   idx = 3*(rama_restraint.atom_index_1);
   clipper::Coord_orth P1(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_2);
   clipper::Coord_orth P2(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_3);
   clipper::Coord_orth P3(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_4);
   clipper::Coord_orth P4(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));
   idx = 3*(rama_restraint.atom_index_5);
   clipper::Coord_orth P5(gsl_vector_get(v,idx),
			  gsl_vector_get(v,idx+1),
			  gsl_vector_get(v,idx+2));

//    P1 = clipper::Coord_orth(1.0, 0.0, 1.0); 
//    P2 = clipper::Coord_orth(0.0, -1.0, 1.0); 
//    P3 = clipper::Coord_orth(0.0, 0.0, 0.0); 
//    P4 = clipper::Coord_orth(-1.0, -1.0, 1.0); 
//    P4 = clipper::Coord_orth(1.0, 1.0, 1.0); 

   clipper::Coord_orth a = P2 - P1;
   clipper::Coord_orth b = P3 - P2;
   clipper::Coord_orth c = P4 - P3;
   clipper::Coord_orth d = P5 - P4;

   // Old (6 atom) wrong:
   // TRANS    psi      1 N (P1)    1 CA (P2)     1 C  (P3)    2 N (P4)
   // TRANS    phi      1 C (P3)    2  N (P4)     2 CA (P5)    2 C (P6)
   //
   // New assignements:
   // TRANS    phi    (1st C) (2nd N ) (2nd CA) (2nd C) 
   // TRANS    psi    (2nd N) (2nd CA) (2nd C ) (3nd N)
   // 
   // So Rama_atoms in this order:
   //   0       1        2      3         4
   //  P1      P2       P3     P4        P5
   // (1st C) (2nd N) (2nd CA) (2nd C) (3rd N)

   // ---------- phi ------------------
   // b*b * [ a.(bxc)/b ]
   double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
      sqrt( b.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
      + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

   double phi = clipper::Util::rad2d(atan2(E,G));
   if (phi < 180.0)
      phi += 360.0;
   if (phi > 180.0)
      phi -= 360.0;

   // ---------- psi ------------------
   // b*b * [ a.(bxc)/b ]
   double H = clipper::Coord_orth::dot(b, clipper::Coord_orth::cross(c,d)) *
      sqrt( c.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double I = - clipper::Coord_orth::dot(b,d)*c.lengthsq()
      + clipper::Coord_orth::dot(b,c)*clipper::Coord_orth::dot(c,d);

   double psi = clipper::Util::rad2d(atan2(H,I));
   if (psi < 180.0)
      psi += 360.0;
   if (psi > 180.0)
      psi -= 360.0;

   // double lr_kdc = lograma.interp(clipper::Util::d2rad(phi), clipper::Util::d2rad(psi));
   // double R = 10.0 * lr;

   std::string residue_type = rama_restraint.rama_plot_residue_type;
   double lr = rama.value(residue_type, clipper::Util::d2rad(phi), clipper::Util::d2rad(psi));
   double R = -rama_plot_weight * lr;

   // std::cout << "zo-rama-distortion for " << phi << " " << psi << " is " << lr << " kdc: "
   // << lr_kdc << "\n";

   if ( clipper::Util::isnan(phi) ) {
      std::cout << "WARNING: observed torsion phi is a NAN!" << std::endl;
      std::cout << "         debug-info: " << E << "/" << G << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_1 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_2 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_3 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_4 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_5 << std::endl;
      std::cout << "         debug-info: P1: " << P1.format() << std::endl;
      std::cout << "         debug-info: P2: " << P2.format() << std::endl;
      std::cout << "         debug-info: P3: " << P3.format() << std::endl;
      std::cout << "         debug-info: P4: " << P4.format() << std::endl;
      std::cout << "         debug-info: P5: " << P5.format() << std::endl;
      std::cout << "         debug-info: a: " << a.format() << std::endl;
      std::cout << "         debug-info: b: " << b.format() << std::endl;
      std::cout << "         debug-info: c: " << c.format() << std::endl;
      std::cout << "         debug-info: d: " << d.format() << std::endl;

      for (unsigned int i=0; i<15; i++) { 
	 std::cout << "           in distortion_score_rama() " << i << " "
		   << gsl_vector_get(v, 3*i  ) << " " 
		   << gsl_vector_get(v, 3*i+1) << " " 
		   << gsl_vector_get(v, 3*i+2) << " " << std::endl;
      }
   } 
   if ( clipper::Util::isnan(psi) ) {
      std::cout << "WARNING: observed torsion psi is a NAN!" << std::endl;
      std::cout << "         debug-info: " << H << "/" << I << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_1 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_2 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_3 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_4 << std::endl;
      std::cout << "         debug-info: atom indices: " << rama_restraint.atom_index_5 << std::endl;
      std::cout << "         debug-info: P1: " << P1.format() << std::endl;
      std::cout << "         debug-info: P2: " << P2.format() << std::endl;
      std::cout << "         debug-info: P3: " << P3.format() << std::endl;
      std::cout << "         debug-info: P4: " << P4.format() << std::endl;
      std::cout << "         debug-info: P5: " << P5.format() << std::endl;
      std::cout << "         debug-info: a: " << a.format() << std::endl;
      std::cout << "         debug-info: b: " << b.format() << std::endl;
      std::cout << "         debug-info: c: " << c.format() << std::endl;
      std::cout << "         debug-info: d: " << d.format() << std::endl;
   } 

   return R;
}

double
coot::distortion_score_non_bonded_contact(const coot::simple_restraint &nbc_restraint,
					  const gsl_vector *v) {

   const int &idx_1 = 3*(nbc_restraint.atom_index_1); 
   const int &idx_2 = 3*(nbc_restraint.atom_index_2);
   
//    clipper::Coord_orth a1(gsl_vector_get(v, idx_1  ), 
// 			     gsl_vector_get(v, idx_1+1), 
// 			     gsl_vector_get(v, idx_1+2));
//    clipper::Coord_orth a2(gsl_vector_get(v, idx_2 ),
// 			     gsl_vector_get(v, idx_2+1), 
// 			     gsl_vector_get(v, idx_2+2));

//    double dist_sq = (a1-a2).lengthsq();

   double dist_sq = 0.0;

   double delta = gsl_vector_get(v, idx_1) - gsl_vector_get(v, idx_2);
   dist_sq += delta * delta;
   delta = gsl_vector_get(v, idx_1+1) - gsl_vector_get(v, idx_2+1);
   dist_sq += delta * delta;
   delta = gsl_vector_get(v, idx_1+2) - gsl_vector_get(v, idx_2+2);
   dist_sq += delta * delta;

   double r = 0.0;

   if (nbc_restraint.fixed_atom_flags[0] && nbc_restraint.fixed_atom_flags[1])
      return 0.0;

   if (false)
      std::cout << "in distortion_score_non_bonded_contact: " << idx_1 << " " << idx_2
		<< " comparing model: " << sqrt(dist_sq) << " min_dist: "
		<< nbc_restraint.target_value
		<< " with sigma " << nbc_restraint.sigma << std::endl;

   if (dist_sq < nbc_restraint.target_value * nbc_restraint.target_value) {
      double weight = 1.0/(nbc_restraint.sigma * nbc_restraint.sigma);
      double dist = sqrt(dist_sq);
      double bit = dist - nbc_restraint.target_value;
      r = weight * bit * bit;
   }
   return r;
}

coot::plane_distortion_info_t
coot::distortion_score_plane_internal(const coot::simple_restraint &plane_restraint,
				      const gsl_vector *v) {


   coot::plane_distortion_info_t info;
   double sum_devi = 0; // return this (after weighting)

   // Recall the equation of a plane: ax + by + cz + d = 0:
   // 
   // First find the centres of the sets of atoms: x_cen, y_cen,
   // z_cen.  We move the plane down so that it crosses the origin and
   // there for d=0 (we'll add it back later).  This makes the problem
   // into 3 equations, 3 unknowns, an eigenvalue problem, with the
   // smallest eigenvalue corresponding to the best fit plane.
   //
   // Google for: least squares plane:
   // http://www.infogoaround.org/JBook/LSQ_Plane.html
   //
   //
   //
   // So, first get the centres:
   double sum_x = 0, sum_y = 0, sum_z = 0;
   int idx;
   int n_atoms = plane_restraint.plane_atom_index.size();
   double dn_atoms = double(n_atoms); 

   if (n_atoms > 0) { 
      for (int i=0; i<n_atoms; i++) {

	 idx = 3*(plane_restraint.plane_atom_index[i].first);

// 	 std::cout << "atom_index vs n_atom " << plane_restraint.atom_index[i]
// 		   << " " << n_atoms << std::endl;
	 if (plane_restraint.plane_atom_index[i].first < 0) {
	    std::cout << "trapped bad plane restraint! " << plane_restraint.plane_atom_index[i].first
		      << std::endl;
	    // return info;
	 } else {
	    sum_x += gsl_vector_get(v,idx);
	    sum_y += gsl_vector_get(v,idx+1);
	    sum_z += gsl_vector_get(v,idx+2);
	 }
      }
      double x_cen = sum_x/dn_atoms;
      double y_cen = sum_y/dn_atoms;
      double z_cen = sum_z/dn_atoms;

      clipper::Matrix<double> mat(3,3);

      for (int i=0; i<n_atoms; i++) {
	 idx = 3*(plane_restraint.plane_atom_index[i].first);
	 // std::cout << "plane restraint adding plane deviations for atom index "
	 // << plane_restraint.atom_index[i] << std::endl;
	 if (plane_restraint.plane_atom_index[i].first < 0) {
	 } else { 
	    mat(0,0) += (gsl_vector_get(v,idx  ) - x_cen) * (gsl_vector_get(v,idx  ) - x_cen);
	    mat(1,1) += (gsl_vector_get(v,idx+1) - y_cen) * (gsl_vector_get(v,idx+1) - y_cen);
	    mat(2,2) += (gsl_vector_get(v,idx+2) - z_cen) * (gsl_vector_get(v,idx+2) - z_cen);
	    mat(0,1) += (gsl_vector_get(v,idx  ) - x_cen) * (gsl_vector_get(v,idx+1) - y_cen);
	    mat(0,2) += (gsl_vector_get(v,idx  ) - x_cen) * (gsl_vector_get(v,idx+2) - z_cen);
	    mat(1,2) += (gsl_vector_get(v,idx+1) - y_cen) * (gsl_vector_get(v,idx+2) - z_cen);
	 }
      }
      mat(1,0) = mat(0,1);
      mat(2,0) = mat(0,2);
      mat(2,1) = mat(1,2);

      if (false) { // debug
	 std::cout << "mat pre  eigens:\n";
	 for (unsigned int ii=0; ii<3; ii++) { 
	    for (unsigned int jj=0; jj<3; jj++) { 
	       std::cout << "mat(" << ii << "," << jj << ") = " << mat(ii,jj) << "   ";
	    }
	    std::cout << "\n";
	 }
      }

      std::vector<double> eigens = mat.eigen(true);

      // Let's now extract the values of a,b,c normalize them
      std::vector<double> abcd(4);
      abcd[0] = mat(0,0);
      abcd[1] = mat(1,0);
      abcd[2] = mat(2,0);

      if (0) { // debug
	 std::cout << "mat post eigens:\n";
	 for (unsigned int ii=0; ii<3; ii++) { 
	    for (unsigned int jj=0; jj<3; jj++) { 
	       std::cout << "mat(" << ii << "," << jj << ") = " << mat(ii,jj) << "   ";
	    }
	    std::cout << "\n";
	 }
      }
      
      double sqsum = 1e-20;

      for (int i=0; i<3; i++)
	 sqsum += abcd[i] * abcd[i];
      for (int i=0; i<3; i++)
	 abcd[i] /= sqsum;


      //set d, recall di = Axi+Byi+Czi-D, so xi = x_cen, yi = y_cen, zi = z_cen:
      abcd[3] = abcd[0]*x_cen + abcd[1]*y_cen + abcd[2]*z_cen;
      info.abcd = abcd;

      double val;
      for (int i=0; i<n_atoms; i++) {
	 idx = 3*(plane_restraint.plane_atom_index[i].first);
	 if (idx < 0) {
	 } else { 
	    // if (! plane_restraint.fixed_atom_flags[i] ) {
	    if (true) { // should fixed atoms contribute to the distortion of the plane?
	                // yes.
	       val = 
		  abcd[0]*gsl_vector_get(v,idx  ) +
		  abcd[1]*gsl_vector_get(v,idx+1) +
		  abcd[2]*gsl_vector_get(v,idx+2) -
		  abcd[3];
	       double r = val/plane_restraint.plane_atom_index[i].second; // .second is the weight
	       sum_devi += r*r;
	    }
	 }
      }

      if (false) {
	 std::cout << plane_restraint << " " << sum_devi << std::endl;
      }
   }

   if (false)
      std::cout << "DEBUG:: distortion_score_plane_internal() returning "
		<< sum_devi << " for " << plane_restraint.plane_atom_index.size() << " atoms"
		<< std::endl;
   
    // info.distortion_score = sum_devi / (plane_restraint.sigma * plane_restraint.sigma);
    // now individually weighted
    info.distortion_score = sum_devi;
   
   return info;
}

// Like the above, except use atom indices, not a restraint.
// 
// I copied and edited, rather than abstracted&refactored for speed reasons.
//
// 
coot::plane_distortion_info_t
coot::distortion_score_2_planes(const std::vector<std::pair<int, double> > &atom_index_set_1,
				const std::vector<std::pair<int, double> > &atom_index_set_2,
				const double &restraint_sigma,
				const gsl_vector *v) {

   coot::plane_distortion_info_t info;
   double sum_devi = 0; // return this (after weighting)

   // Recall the equation of a plane: ax + by + cz + d = 0:
   // 
   // First find the centres of the sets of atoms: x_cen, y_cen,
   // z_cen.  We move the plane down so that it crosses the origin and
   // there for d=0 (we'll add it back later).  This makes the problem
   // into 3 equations, 3 unknowns, an eigenvalue problem, with the
   // smallest eigenvalue corresponding to the best fit plane.
   //
   // Google for: least squares plane:
   // http://www.infogoaround.org/JBook/LSQ_Plane.html
   //
   //
   //
   // So, first get the centres:
   double sum_x = 0, sum_y = 0, sum_z = 0;
   int idx;

   if (atom_index_set_1.size() > 2 && atom_index_set_2.size() > 2) {

      double f_1 = 1.0/double(atom_index_set_1.size());
      double f_2 = 1.0/double(atom_index_set_2.size());
      
      // plane 1
      // 
      for (unsigned int i=0; i<atom_index_set_1.size(); i++) {
	 idx = 3*(atom_index_set_1[i].first);
	 if (atom_index_set_1[i].first < 0) {
	    std::cout << "trapped bad par plane restraint! " << atom_index_set_1[i].first
		      << std::endl;
	 } else {
	    sum_x += gsl_vector_get(v,idx);
	    sum_y += gsl_vector_get(v,idx+1);
	    sum_z += gsl_vector_get(v,idx+2);
	 }
      }
      double x_cen_1 = sum_x * f_1;
      double y_cen_1 = sum_y * f_1;
      double z_cen_1 = sum_z * f_1;
      info.centre_1 = clipper::Coord_orth(x_cen_1, y_cen_1, z_cen_1);
      
      
      // plane 2
      // 
      sum_x = 0; sum_y = 0; sum_z = 0;
      for (unsigned int i=0; i<atom_index_set_2.size(); i++) {
	 idx = 3*(atom_index_set_2[i].first);
	 if (atom_index_set_2[i].first < 0) {
	    std::cout << "trapped bad par plane restraint! " << atom_index_set_2[i].first
		      << std::endl;
	 } else {
	    sum_x += gsl_vector_get(v,idx);
	    sum_y += gsl_vector_get(v,idx+1);
	    sum_z += gsl_vector_get(v,idx+2);
	 }
      }
      double x_cen_2 = sum_x * f_2;
      double y_cen_2 = sum_y * f_2;
      double z_cen_2 = sum_z * f_2;
      info.centre_2 = clipper::Coord_orth(x_cen_2, y_cen_2, z_cen_2);

      clipper::Matrix<double> mat(3,3);
      
      for (unsigned int i=0; i<atom_index_set_1.size(); i++) {
	 idx = 3*(atom_index_set_1[i].first);
	 if (atom_index_set_1[i].first < 0) {
	 } else { 
	    mat(0,0) += (gsl_vector_get(v,idx  ) - x_cen_1) * (gsl_vector_get(v,idx  ) - x_cen_1);
	    mat(1,1) += (gsl_vector_get(v,idx+1) - y_cen_1) * (gsl_vector_get(v,idx+1) - y_cen_1);
	    mat(2,2) += (gsl_vector_get(v,idx+2) - z_cen_1) * (gsl_vector_get(v,idx+2) - z_cen_1);
	    mat(0,1) += (gsl_vector_get(v,idx  ) - x_cen_1) * (gsl_vector_get(v,idx+1) - y_cen_1);
	    mat(0,2) += (gsl_vector_get(v,idx  ) - x_cen_1) * (gsl_vector_get(v,idx+2) - z_cen_1);
	    mat(1,2) += (gsl_vector_get(v,idx+1) - y_cen_1) * (gsl_vector_get(v,idx+2) - z_cen_1);
	 }
      }
      for (unsigned int i=0; i<atom_index_set_2.size(); i++) {
	 idx = 3*(atom_index_set_2[i].first);
	 if (atom_index_set_2[i].first < 0) {
	 } else { 
	    mat(0,0) += (gsl_vector_get(v,idx  ) - x_cen_2) * (gsl_vector_get(v,idx  ) - x_cen_2);
	    mat(1,1) += (gsl_vector_get(v,idx+1) - y_cen_2) * (gsl_vector_get(v,idx+1) - y_cen_2);
	    mat(2,2) += (gsl_vector_get(v,idx+2) - z_cen_2) * (gsl_vector_get(v,idx+2) - z_cen_2);
	    mat(0,1) += (gsl_vector_get(v,idx  ) - x_cen_2) * (gsl_vector_get(v,idx+1) - y_cen_2);
	    mat(0,2) += (gsl_vector_get(v,idx  ) - x_cen_2) * (gsl_vector_get(v,idx+2) - z_cen_2);
	    mat(1,2) += (gsl_vector_get(v,idx+1) - y_cen_2) * (gsl_vector_get(v,idx+2) - z_cen_2);
	 }
      }
      mat(1,0) = mat(0,1);
      mat(2,0) = mat(0,2);
      mat(2,1) = mat(1,2);

      if (0) { // debug
	 std::cout << "mat pre  eigens:\n";
	 for (unsigned int ii=0; ii<3; ii++) { 
	    for (unsigned int jj=0; jj<3; jj++) { 
	       std::cout << "mat(" << ii << "," << jj << ") = " << mat(ii,jj) << "   ";
	    }
	    std::cout << "\n";
	 }
      }
      
      std::vector<double> eigens = mat.eigen(true);

      if (0) 
	 std::cout << "we get eigen values: "
		   << eigens[0] << "  "
		   << eigens[1] << "  "
		   << eigens[2] << std::endl;

      
      // Let's now extract the values of a,b,c normalize them
      std::vector<double> abcd(4);
      abcd[0] = mat(0,0);
      abcd[1] = mat(1,0);
      abcd[2] = mat(2,0);

      if (0) { // debug
	 std::cout << "mat post eigens:\n";
	 for (unsigned int ii=0; ii<3; ii++) { 
	    for (unsigned int jj=0; jj<3; jj++) { 
	       std::cout << "mat(" << ii << "," << jj << ") = " << mat(ii,jj) << "   ";
	    }
	    std::cout << "\n";
	 }
      }
      
      double sqsum = 1e-20;

      for (int i=0; i<3; i++)
	 sqsum += abcd[i] * abcd[i];
      for (int i=0; i<3; i++)
	 abcd[i] /= sqsum;

      abcd[3] = 0;  // we move the combined plane to the origin
      info.abcd = abcd;

      double val;
      for (unsigned int i=0; i<atom_index_set_1.size(); i++) {
	 idx = 3*(atom_index_set_1[i].first);
	 if (idx < 0) {
	 } else { 
	    val = 
	       abcd[0]*(gsl_vector_get(v,idx  ) - x_cen_1) +
	       abcd[1]*(gsl_vector_get(v,idx+1) - y_cen_1) +
	       abcd[2]*(gsl_vector_get(v,idx+2) - z_cen_1);
	    sum_devi += val * val;
	 }
      }
      for (unsigned int i=0; i<atom_index_set_2.size(); i++) {
	 idx = 3*(atom_index_set_2[i].first);
	 if (idx < 0) {
	 } else { 
	    val = 
	       abcd[0]*(gsl_vector_get(v,idx  ) - x_cen_2) +
	       abcd[1]*(gsl_vector_get(v,idx+1) - y_cen_2) +
	       abcd[2]*(gsl_vector_get(v,idx+2) - z_cen_2);
	    sum_devi += val * val;
	 }
      }
   }
   info.distortion_score = sum_devi / (restraint_sigma * restraint_sigma);
   return info;
}

#endif // HAVE_GSL
