/* coot-utils/coot-coord-utils-mmdb-free.cc
 *
 * The mmdb-free functions of coot-coord-utils.cc, split out so they can be
 * linked into builds that do not use mmdb2 (e.g. the gemmi backend). These are
 * exactly the coot::util definitions whose bodies touch no mmdb type; they are
 * declared (as before) in coot-coord-utils.hh and stay in namespace coot::util.
 * Moved verbatim from coot-coord-utils.cc.
 *
 * Copyright 2006-2016 as per coot-coord-utils.cc; Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#include <algorithm>
#include <stdexcept>
#include <string.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>

#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"

#include "compat/coot-sysdep.h"
#include "clipper/mmdb/clipper_mmdb.h"
#include "geometry/main-chain.hh"
#include "geometry/mol-utils.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "utils/logging.hh"
extern logging logger;

std::vector<std::string>
coot::util::standard_residue_types() {

   std::vector<std::string> v;
   v.push_back("ALA"); v.push_back("ARG"); v.push_back("ASP");
   v.push_back("ASN"); v.push_back("CYS"); v.push_back("SER");
   v.push_back("PRO"); v.push_back("PHE"); v.push_back("GLY");
   v.push_back("GLU"); v.push_back("GLN"); v.push_back("ILE");
   v.push_back("LEU"); v.push_back("TYR"); v.push_back("TRP");
   v.push_back("HIS"); v.push_back("LYS"); v.push_back("MET");
   v.push_back("VAL"); v.push_back("THR"); v.push_back("MSE"); 
   return v;
}

std::vector<std::string>
coot::util::PDB_standard_residue_types() {

   std::vector<std::string> v = coot::util::standard_residue_types();
   v.push_back("Td"); v.push_back("Tr"); v.push_back("T");
   v.push_back("Gd"); v.push_back("Gr"); v.push_back("G");
   v.push_back("Ad"); v.push_back("Ar"); v.push_back("A");

   v.push_back("DG"); v.push_back("DC"); v.push_back("DA");
   v.push_back("DU"); v.push_back("DT"); v.push_back("DI");

   v.push_back("UNK"); v.push_back("N");

   return v;
}

std::ostream&
coot::operator<<(std::ostream&  s, const coot::lsq_range_match_info_t &m) {

   s << "LSQ Match: (" << m.model_number_reference << ") " << m.reference_chain_id << " "
     << m.to_reference_start_resno << "-" << m.to_reference_end_resno
     << " to ("
     << m.model_number_matcher << ") " << m.matcher_chain_id << " "
     << m.from_matcher_start_resno << "-" << m.from_matcher_end_resno
     << " type: " << m.match_type_flag;
   return s;
}

float
coot::util::interquartile_range(const std::vector<float> &v_in) {

   float iqr = 0;
   std::vector<float> v = v_in;

   std::sort(v.begin(), v.end());
   unsigned int n = v.size();
   int q_1 = int(0.25 * n);
   int q_3 = int(0.75 * n);
   float v_1 = v[q_1];
   float v_3 = v[q_3];
   iqr = v_3 - v_1;
   return iqr;
}

coot::util::stats_data::stats_data(const std::vector<float> &v) {

   mean = 0;
   sd = 0;
   iqr = 0;
   double sum = 0;
   double sum_sq = 0;
   for (unsigned int i=0; i<v.size(); i++) {
      sum += v[i];
      sum_sq += v[i] * v[i];
   }
   if (v.size() > 0) {
      mean = sum/double(v.size());
      double var = sum_sq/double(v.size()) - mean * mean;
      if (var < 0) var = 0;
      sd = sqrt(var);
      iqr = interquartile_range(v);
   }
}

coot::util::stats_data::stats_data(const std::vector<double> &v) {

   mean = 0;
   sd = 0;
   iqr = 0;
   double sum = 0;
   double sum_sq = 0;
   for (unsigned int i=0; i<v.size(); i++) {
      sum += v[i];
      sum_sq += (v[i] * v[i]);
   }
   if (v.size() > 0) { 
      mean = sum/double(v.size());
      double var = sum_sq/double(v.size()) - mean * mean;
      if (var < 0) var = 0;
      sd = sqrt(var);
      std::vector<float> vf(v.size());
      for (unsigned int i=0; i<v.size(); i++) vf[i] = v[i];
      iqr = interquartile_range(vf);
   }
}

std::vector<std::pair<double, double> >
coot::util::qq_plot_t::qq_norm() {
   
   std::vector<std::pair<double, double> > v;

   std::sort(data.begin(), data.end());
   size_t stride = 1;
   std::vector<double> sorted_data(data.size());
   for (unsigned int i=0; i<data.size(); i++) { 
      sorted_data[i] = data[i];
   }
   stats_data sd = stats_data(data);

   // debugging
   std::vector<double> save_gs;
   std::vector<double> save_qs;
   

   // frac goes between 0 and 1.
   for (double frac=0.01; frac<1; frac+= 0.01) {

      double g = gsl_cdf_gaussian_Pinv(frac, sd.sd);
      double q = gsl_stats_quantile_from_sorted_data(sorted_data.data(), stride,
                                                     data.size(), frac);

      // mean correction (gs would otherwise have mean 0)
      double g_mc = g + sd.mean;
      if (0) 
         std::cout << "debug:: g " << g << " from frac " << frac
                   << " and sd " << sd.sd << std::endl;

      std::pair<double, double> p(g_mc, q);
      v.push_back(p);
      save_gs.push_back(g_mc);
      save_qs.push_back(q);
   }

   stats_data gs_data(save_gs);
   stats_data qs_data(save_qs);
   
   std::cout << "debug:: gs: mean " << gs_data.mean << " sd " <<  gs_data.sd  << std::endl;
   std::cout << "debug:: qs: mean " << qs_data.mean << " sd " <<  qs_data.sd  << std::endl;
   std::cout << "debug:: sd: mean " << sd.mean << " sd " <<  sd.sd  << std::endl;
   
   return v;
}

clipper::Mat33<double>
coot::util::quaternion::matrix() const {

   clipper::Mat33<double> mat;

   mat(0,0) = 1.0 - 2.0 * (q1 * q1 + q2 * q2);
   mat(0,1) = 2.0 * (q0 * q1 - q2 * q3);
   mat(0,2) = 2.0 * (q2 * q0 + q1 * q3);
   
   mat(1,0) = 2.0 * (q0 * q1 + q2 * q3);
   mat(1,1)= 1.0 - 2.0 * (q2 * q2 + q0 * q0);
   mat(1,2) = 2.0 * (q1 * q2 - q0 * q3);
   
   mat(2,0) = 2.0 * (q2 * q0 - q1 * q3);
   mat(2,1) = 2.0 * (q1 * q2 + q0 * q3);
   mat(2,2) = 1.0 - 2.0 * (q1 * q1 + q0 * q0);
   
   return mat;
}

coot::util::quaternion::quaternion(const clipper::Mat33<double> &m) {

   float pw = 1 + m(0,0) + m(1,1) + m(2,2);
   float px = 1 + m(0,0) - m(1,1) - m(2,2);
   float py = 1 - m(0,0) + m(1,1) - m(2,2); 
   float pz = 1 - m(0,0) - m(1,1) + m(2,2);

   float pr1 = sqrt( (pw>0) ? pw : 0) / 2.0;
   float pr2 = sqrt( (px>0) ? px : 0) / 2.0;
   float pr3 = sqrt( (py>0) ? py : 0) / 2.0;
   float pr4 = sqrt( (pz>0) ? pz : 0) / 2.0;

   q0 = convert_sign(pr2, m(2,1) - m(1,2));
   q1 = convert_sign(pr3, m(0,2) - m(2,0));
   q2 = convert_sign(pr4, m(1,0) - m(0,1));
   q3 = pr1;
   
}

float 
coot::util::quaternion::convert_sign(const float &x, const float &y) const {

   if ((x > 0) && (y > 0)) return  x;
   if ((x < 0) && (y > 0)) return -x; 
   if ((x > 0) && (y < 0)) return -x; 
   return  x; 
}

void
coot::util::quaternion::normalize() {

   double sum_sq = 0.0;
   sum_sq += q0*q0;
   sum_sq += q1*q1;
   sum_sq += q2*q2;
   sum_sq += q3*q3;
   if (sum_sq > 0.0) {
      double f = sqrt(1.0/sum_sq);
      q0 *= f;
      q1 *= f;
      q2 *= f;
      q3 *= f;
   }
} 

coot::util::quaternion
coot::util::quaternion::rotate(double angle, const clipper::Coord_orth &vec) const {

   coot::util::quaternion q(0,0,0,1);

   std::cout << "rotate() just a stub - fill me later!" << std::endl;

   return q;
}

coot::util::quaternion
coot::util::quaternion::inverse() const {

   coot::util::quaternion q(q0, q1, q2, -q3);
   return q;
}

clipper::RTop_orth
coot::util::quaternion::centroid_rtop(const std::vector<std::pair<clipper::RTop_orth,float> > &rtops) {

   if (rtops.size() == 0) { 
      return clipper::RTop_orth(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), clipper::Vec3<double>(0,0,0));
   } else {
      clipper::Coord_orth sum_trn(0,0,0);
      for (unsigned int i=0; i<rtops.size(); i++) { 
         quaternion q(rtops[i].first.rot());
         q0 += rtops[i].second * q.q0;
         q1 += rtops[i].second * q.q1;
         q2 += rtops[i].second * q.q2;
         q3 += rtops[i].second * q.q3;
         sum_trn += rtops[i].first.trn();
      }
      normalize();
      clipper::Mat33<double> m = matrix();
      double inv_n = 1.0/double(rtops.size());
      clipper::Coord_orth t(sum_trn.x() * inv_n, sum_trn.y() * inv_n, sum_trn.z() * inv_n);
      return clipper::RTop_orth(m, t);
   }
}

clipper::RTop_orth
coot::util::quaternion::centroid_rtop(const std::vector<std::pair<clipper::RTop_orth,float> > &rtops,
                                      bool robust_filter) {

   if (! robust_filter) { 
      return centroid_rtop(rtops);
   } else {
      if (rtops.size() < 2) {
         return centroid_rtop(rtops);
      } else { 
         clipper::Coord_orth sum_trn(0,0,0);
         for (unsigned int i=0; i<rtops.size(); i++) { 
            quaternion q(rtops[i].first.rot());
            q0 += rtops[i].second * q.q0;
            q1 += rtops[i].second * q.q1;
            q2 += rtops[i].second * q.q2;
            q3 += rtops[i].second * q.q3;
            sum_trn += rtops[i].first.trn();
         }
         normalize();

         double inv_n = 1.0/double(rtops.size());
         clipper::Coord_orth t(sum_trn.x() * inv_n, sum_trn.y() * inv_n, sum_trn.z() * inv_n);

         double sum_rotation_distance_sq = 0.0;
         double sum_translation_distance_sq = 0.0;
         std::vector<w_rtop_orth> deviance(rtops.size());
         for (unsigned int i=0; i<rtops.size(); i++) {
            // rotation
            quaternion q(rtops[i].first.rot());
            double d0 = q0 - rtops[i].second * q.q0;
            double d1 = q1 - rtops[i].second * q.q1;
            double d2 = q2 - rtops[i].second * q.q2;
            double d3 = q3 - rtops[i].second * q.q3;
            double d = d0*d0 + d1*d1 + d2*d2 + d3*d3;
            sum_rotation_distance_sq += d;
         
            // translation
            clipper::Coord_orth wpt(rtops[i].second * clipper::Coord_orth(rtops[i].first.trn()));
            double dt = (t-wpt).lengthsq();
            sum_translation_distance_sq += dt;
            if (0)
               std::cout << "for irtop " << i << " added rotation distance_sq " << d 
                         << " and translation distance "<< dt << std::endl;
            deviance[i].rtop   = rtops[i].first;
            deviance[i].weight = rtops[i].second;
            deviance[i].deviance = d * 1.888 + dt;
         }

         std::sort(deviance.begin(), deviance.end(), deviance_sorter);
         for (unsigned int i=0; i<deviance.size(); i++)
            if (0) 
               std::cout << "        deviance " << i << " " << deviance[i].weight << " "
                         << deviance[i].deviance << std::endl;
         std::vector<float> iqr_data(deviance.size());
         for (unsigned int i=0; i<deviance.size(); i++)
            iqr_data[i] = deviance[i].deviance;
         stats_data sd(iqr_data);
      
         clipper::Coord_orth sum_trn_filtered_dev(0,0,0);
         int n = 0;

         for (unsigned int i=0; i<deviance.size(); i++) {
            if (deviance[i].deviance < sd.mean + 0.5 * sd.iqr) {
               n++;
               quaternion q(rtops[i].first.rot());
               q0 += deviance[i].weight * q.q0;
               q1 += deviance[i].weight * q.q1;
               q2 += deviance[i].weight * q.q2;
               q3 += deviance[i].weight * q.q3;
               sum_trn_filtered_dev += deviance[i].weight * deviance[i].rtop.trn();
            }
         }

         //       std::cout << "rejecting deviances more than " << sd.mean + 0.5 * sd.iqr
         //                 << " leaves " << n << " from " << rtops.size() 
         //                 << std::endl;
      
         if (n > 0) {

            normalize();
            double inv_n_local = 1.0/double(n);
            clipper::Mat33<double> m = matrix();
            clipper::Coord_orth td(sum_trn_filtered_dev.x() * inv_n_local,
                                   sum_trn_filtered_dev.y() * inv_n_local,
                                   sum_trn_filtered_dev.z() * inv_n_local);
            return clipper::RTop_orth(m, td);
         } else { 
       
            // unfiltered
            clipper::Mat33<double> m = matrix();
            return clipper::RTop_orth(m, t);
         }
      }
   } 
}

bool
coot::util::quaternion::deviance_sorter(const w_rtop_orth &a,
                                        const w_rtop_orth &b) {
   return a.deviance < a.deviance;
}

void
coot::util::quaternion::test_quaternion() {

   // currently quaternions are tested in testcootutils
} 

std::ostream&  coot::util::operator<<(std::ostream& s, const coot::util::quaternion &q) {

   s << "(" << q.q0 << ", " << q.q1 << ", " << q.q2 << ", " << q.q3 << ")";
   return s;
}

std::ofstream& coot::util::operator<<(std::ofstream &s, const coot::util::quaternion &q) {

   // s << "(" << q.q0 << ", " << q.q1 << ", " << q.q2 << ", " << q.q3 << ")";
   return s;
} 

std::string
coot::util::single_letter_to_3_letter_code(char code) {

   if (code == 'G') return std::string("GLY");
   if (code == 'A') return std::string("ALA");
   if (code == 'V') return std::string("VAL");
   if (code == 'S') return std::string("SER");
   if (code == 'N') return std::string("ASN");
   if (code == 'P') return std::string("PRO");
   if (code == 'D') return std::string("ASP");
   if (code == 'C') return std::string("CYS");
   if (code == 'Q') return std::string("GLN");
   if (code == 'E') return std::string("GLU");
   if (code == 'H') return std::string("HIS");
   if (code == 'I') return std::string("ILE");
   if (code == 'L') return std::string("LEU");
   if (code == 'K') return std::string("LYS");
   if (code == 'M') return std::string("MET");
   if (code == 'F') return std::string("PHE");
   if (code == 'T') return std::string("THR");
   if (code == 'W') return std::string("TRP");
   if (code == 'Y') return std::string("TYR");
   if (code == 'R') return std::string("ARG");

   return std::string("");
}

std::string
coot::util::single_letter_to_3_letter_code(const std::string &code) {

   if (code == "G") return std::string("GLY");
   if (code == "A") return std::string("ALA");
   if (code == "V") return std::string("VAL");
   if (code == "S") return std::string("SER");
   if (code == "N") return std::string("ASN");
   if (code == "P") return std::string("PRO");
   if (code == "D") return std::string("ASP");
   if (code == "C") return std::string("CYS");
   if (code == "Q") return std::string("GLN");
   if (code == "E") return std::string("GLU");
   if (code == "H") return std::string("HIS");
   if (code == "I") return std::string("ILE");
   if (code == "L") return std::string("LEU");
   if (code == "K") return std::string("LYS");
   if (code == "M") return std::string("MET");
   if (code == "F") return std::string("PHE");
   if (code == "T") return std::string("THR");
   if (code == "W") return std::string("TRP");
   if (code == "Y") return std::string("TYR");
   if (code == "R") return std::string("ARG");

   return std::string("");
}

std::string
coot::graph_match_info_t::invent_new_name(const std::string &name_in,
                                          const std::string &ele,
                                          const std::vector<std::string> &residue_atom_name) const {
   std::string a("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
   bool found = 0;
   std::string new_name("XXXX");
   for (unsigned int i=0; i<a.size(); i++) { 
      for (unsigned int j=0; j<a.size(); j++) {
         std::string test_atom_name = "";
         if (ele.length() == 1) { 
            test_atom_name = " ";
            test_atom_name += ele;
         } else {
            test_atom_name = ele;
         }
         test_atom_name += a[i];
         test_atom_name += a[j];
         if (std::find(residue_atom_name.begin(), residue_atom_name.end(), test_atom_name)
             == residue_atom_name.end()) {
            found = 1;
            new_name = test_atom_name;
         }
         if (found)
            break;
      }
      if (found)
         break;
   }
   return new_name;
} 

bool
coot::util::chain_id_residue_vec_helper_t::operator<(const chain_id_residue_vec_helper_t &c) const { 

   return (chain_id < c.chain_id);
}

bool
coot::util::is_cis(const double &omega_torsion) {

   // -90 to +90
   bool is_cis_flag = false;
   if ((omega_torsion < 1.57) && (omega_torsion > -1.57))
      is_cis_flag = true;
   return is_cis_flag;
}

std::pair<double, double>
coot::lsq_plane_deviation(const std::vector<clipper::Coord_orth> &v,
                          const clipper::Coord_orth &pt) {

   coot::lsq_plane_info_t lpd(v);
   double val = lpd.plane_deviation(pt);
   double rms = lpd.plane_atoms_rms();
   return std::pair<double, double> (val, rms);
}

coot::lsq_plane_info_t::lsq_plane_info_t(const std::vector<clipper::Coord_orth> &v) {

   int n_atoms = v.size();
   clipper::Coord_orth sum(0,0,0);
   for (int i=0; i<n_atoms; i++)
      sum += v[i];
   double factor = 1/double(n_atoms);
   clipper::Coord_orth midpoint(sum.x()*factor, sum.y()*factor, sum.z()*factor);
   centre_ = midpoint;

   clipper::Matrix<double> mat(3,3);
   for (int i=0; i<n_atoms; i++) {
      mat(0,0) += (v[i].x() - midpoint.x()) * (v[i].x() - midpoint.x());
      mat(1,1) += (v[i].y() - midpoint.y()) * (v[i].y() - midpoint.y());
      mat(2,2) += (v[i].z() - midpoint.z()) * (v[i].z() - midpoint.z());
      mat(0,1) += (v[i].x() - midpoint.x()) * (v[i].y() - midpoint.y());
      mat(0,2) += (v[i].x() - midpoint.x()) * (v[i].z() - midpoint.z());
      mat(1,2) += (v[i].y() - midpoint.y()) * (v[i].z() - midpoint.z());
   }
   mat(1,0) = mat(0,1);
   mat(2,0) = mat(0,2);
   mat(2,1) = mat(1,2);

   if (0) { 
      std::cout << "  mat for eigens: " << std::endl;
      std::cout << "     " << mat(0,0) << "   " << mat(0,1) << "   " << mat(0,2) << std::endl;
      std::cout << "     " << mat(1,0) << "   " << mat(1,1) << "   " << mat(1,2) << std::endl;
      std::cout << "     " << mat(2,0) << "   " << mat(2,1) << "   " << mat(2,2) << std::endl;
   }
   std::vector<double> eigens = mat.eigen(true);
   // Let's now extract the values of a,b,c normalize them
   abcd.resize(4);
   
   abcd[0] = mat(0,0);
   abcd[1] = mat(1,0);
   abcd[2] = mat(2,0);

   if (0) 
      std::cout << " abcd - pre-values "
                << abcd[0] << " "
                << abcd[1] << " "
                << abcd[2] << " "
                << std::endl;
   
   double sqsum = 1e-20;
   
   for (int i=0; i<3; i++)
      sqsum += abcd[i] * abcd[i];
   for (int i=0; i<3; i++)
      abcd[i] /= sqsum;
   
   // set D, recall di = Axi+Byi+Czi-D, so when
   // xi = x_cen, yi = y_cen, zi = z_cen, d is 0,
   // so we can set D.
   // 
   abcd[3] = abcd[0]*midpoint.x() + abcd[1]*midpoint.y() + abcd[2]*midpoint.z();

   if (0) 
      std::cout << " abcd "
                << abcd[0] << " "
                << abcd[1] << " "
                << abcd[2] << " "
                << abcd[3] << std::endl;

   double var = 0;
   for (unsigned int i_plane_at=0; i_plane_at<v.size(); i_plane_at++) {
      double d =
         abcd[0]*v[i_plane_at].x() +
         abcd[1]*v[i_plane_at].y() +
         abcd[2]*v[i_plane_at].z() - abcd[3];
      var += d*d;
   }
   rms = 0;
   if (v.size() > 0)
      rms = sqrt(var/double(v.size()));

}

bool
coot::compare_atom_specs_user_float(const coot::atom_spec_t &a1, const coot::atom_spec_t &a2) {

   return a1.float_user_data < a2.float_user_data ? 1 : 0;

} 

std::string
coot::util::interesting_things_list(const std::vector<atom_spec_t> &v) {

#ifdef USE_GUILE
   // e.g. (list) for empty v
   // (list (list "button label" imol-no chain-id resno atom-name)
   //       (list "button label" imol-no chain-id resno atom-name)
   // )

   std::string r = " (list ";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].chain_id;
      atom_str += "\" ";
      atom_str += int_to_string(v[i].res_no);
      atom_str += " \"";
      atom_str += v[i].ins_code;
      atom_str += "\" \"";
      atom_str += v[i].atom_name;
      atom_str += "\" \"";
      atom_str += v[i].alt_conf;
      atom_str += " \"";

      std::string button_label("Clash gap: ");
      button_label += float_to_string(v[i].float_user_data);
      button_label += " : ";
      button_label += v[i].chain_id;
      button_label += " ";
      button_label += int_to_string(v[i].res_no);
      button_label += " ";
      if (v[i].ins_code != "") {
         button_label += v[i].ins_code;
          button_label += " ";
      }
      button_label += v[i].atom_name;
      if (v[i].alt_conf != "") {
         button_label += ",";
    button_label += v[i].alt_conf;
         button_label += " ";
      }

      std::string s = "(list ";
      s += single_quote(button_label);
      s += " ";
      s += int_to_string(v[i].int_user_data);
      s += " ";
      s += atom_str;
      s += ")\n";
      
      r += s;
   }

   r += ")";
   return r;
#else
#ifdef USE_PYTHON
// BL says:: we want to have [] lists in python, separated by commas (,)
   // e.g. [] for empty v
   // [["button label",imol-no,chain-id,resno,atom-name],
   //  ["button label",imol-no,chain-id,resno,atom-name]
   // ]

   std::string r = "[";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].chain_id;
      atom_str += "\",";
      atom_str += int_to_string(v[i].res_no);
      atom_str += ",\"";
      atom_str += v[i].ins_code;
      atom_str += "\",\"";
      atom_str += v[i].atom_name;
      atom_str += "\",\"";
      atom_str += v[i].alt_conf;
      atom_str += " \"";

      std::string button_label("Clash gap: ");
      button_label += float_to_string(v[i].float_user_data);
      button_label += " : ";
      button_label += v[i].chain_id;
      button_label += " ";
      button_label += int_to_string(v[i].res_no);
      button_label += " ";
      if (v[i].ins_code != "") {
         button_label += v[i].ins_code;
         button_label += " ";
      }
      button_label += v[i].atom_name;
      if (v[i].alt_conf != "") {
         button_label += ",";
         button_label += v[i].alt_conf;
         button_label += " ";
      }

      std::string s = "[";
      s += single_quote(button_label);
      s += ",";
      s += int_to_string(v[i].int_user_data);
      s += ",";
      s += atom_str;
      s += "],\n";

      r += s;
   }

   r += "]";
   return r;
#else    
   return "";
#endif // PYTHON
   
#endif // GUILE
}

std::string
coot::util::interesting_things_list_py(const std::vector<atom_spec_t> &v) {

#ifdef USE_PYTHON
   // BL says:: we want to have [] lists in python, separated by commas (,)
   // e.g. [] for empty v
   // [["button label",imol-no,chain-id,resno,atom-name],
   //  ["button label",imol-no,chain-id,resno,atom-name]
   // ]

   std::string r = "[";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].chain_id;
      atom_str += "\",";
      atom_str += int_to_string(v[i].res_no);
      atom_str += ",\"";
      atom_str += v[i].ins_code;
      atom_str += "\",\"";
      atom_str += v[i].atom_name;
      atom_str += "\",\"";
      atom_str += v[i].alt_conf;
      atom_str += " \"";

      std::string button_label("Clash gap: ");
      button_label += float_to_string(v[i].float_user_data);
      button_label += " : ";
      button_label += v[i].chain_id;
      button_label += " ";
      button_label += int_to_string(v[i].res_no);
      button_label += " ";
      if (v[i].ins_code != "") {
         button_label += v[i].ins_code;
         button_label += " ";
      }
      button_label += v[i].atom_name;
      if (v[i].alt_conf != "") {
         button_label += ",";
         button_label += v[i].alt_conf;
         button_label += " ";
      }

      std::string s = "[";
      s += single_quote(button_label);
      s += ",";
      s += int_to_string(v[i].int_user_data);
      s += ",";
      s += atom_str;
      s += "],\n";

      r += s;
   }

   r += "]";
   return r;
#else
   return "";
#endif // PYTHON
}

std::string
coot::util::interesting_things_list_with_fix(const std::vector<coot::util::atom_spec_and_button_info_t> &v,
                                             const std::string &error_type) {

#ifdef USE_GUILE
   // e.g. (list) for empty v
   // (list (list "button label" imol-no chain-id resno atom-name)
   //       (list "button label" imol-no chain-id resno atom-name)
   // )
   //
   // if we have a fix, the callback function is not "" and then the
   // returned thing becomes:
   // 
   // (list (list "button label" imol-no chain-id resno atom-name callback-func)
   //       (list "button label" imol-no chain-id resno atom-name callback-func)
   // )
   //
   // where callback-func is e.g. (lambda() (do-180-degree-side-chain-flip 0 "A" 45 "" ""))

   std::string r = " (list ";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].as.chain_id;
      atom_str += "\" ";
      atom_str += int_to_string(v[i].as.res_no);
      atom_str += " \"";
      atom_str += v[i].as.ins_code;
      atom_str += "\" \"";
      atom_str += v[i].as.atom_name;
      atom_str += "\" \"";
      atom_str += v[i].as.alt_conf;
      atom_str += " \"";

      std::string button_label = v[i].button_label;

      std::string s = "(list ";
      s += single_quote(button_label);
      s += " ";
      s += int_to_string(v[i].as.int_user_data);
      s += " ";
      s += atom_str;

      if (v[i].callback_func != "") {
         s += " ";
         s +=  v[i].callback_func;
      }
      
      s += ")\n";
      
      r += s;
   }

   r += ")";
   return r;
#else
#ifdef USE_PYTHON
// BL says:: here again we need a [] list in python 
   std::string r = "[";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].as.chain_id;
      atom_str += "\",";
      atom_str += int_to_string(v[i].as.res_no);
      atom_str += ",\"";
      atom_str += v[i].as.ins_code;
      atom_str += "\",\"";
      atom_str += v[i].as.atom_name;
      atom_str += "\",\"";
      atom_str += v[i].as.alt_conf;
      atom_str += " \"";

      std::string button_label = v[i].button_label;

      std::string s = "[";
      s += single_quote(button_label);
      s += ",";
      s += int_to_string(v[i].as.int_user_data);
      s += ",";
      s += atom_str;

      if (v[i].callback_func != "") {
         s += ",";
         s +=  v[i].callback_func;
      }

      if (i<(v.size()-1)) {
         s += "],\n";
      } else {
         s += "]\n";
      }

      r += s;
   }

   r += "]";
   return r;
#else
   return "";
#endif // PYTHON
#endif // GUILE
}

std::string
coot::util::interesting_things_list_with_fix_py(const std::vector<coot::util::atom_spec_and_button_info_t> &v,
                                                const std::string &error_type) {
#ifdef USE_PYTHON
// BL says:: here again we need a [] list in python 
   std::string r = "[";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].as.chain_id;
      atom_str += "\",";
      atom_str += int_to_string(v[i].as.res_no);
      atom_str += ",\"";
      atom_str += v[i].as.ins_code;
      atom_str += "\",\"";
      atom_str += v[i].as.atom_name;
      atom_str += "\",\"";
      atom_str += v[i].as.alt_conf;
      atom_str += " \"";

      std::string button_label = v[i].button_label;

      std::string s = "[";
      s += single_quote(button_label);
      s += ",";
      s += int_to_string(v[i].as.int_user_data);
      s += ",";
      s += atom_str;

      if (v[i].callback_func != "") {
         s += ",";
         s +=  v[i].callback_func;
      }

      if (i<(v.size()-1)) {
         s += "],\n";
      } else {
         s += "]\n";
      }

      r += s;
   }

   r += "]";
   return r;
#else
   return "";
#endif // PYTHON
}

bool
coot::compare_atom_specs_user_float_in_pair(const std::pair<coot::atom_spec_t, std::string> &a,
                                            const std::pair<coot::atom_spec_t, std::string> &b) {

   return b.first.float_user_data < a.first.float_user_data ? 1 : 0;
}

clipper::Coord_orth
coot::util::rotate_around_vector(const clipper::Coord_orth &direction,
                                 const clipper::Coord_orth &position,
                                 const clipper::Coord_orth &origin_shift,
                                 double angle) {
   
   clipper::Coord_orth unit_vec = clipper::Coord_orth(direction.unit());
   
   double l = unit_vec[0];
   double m = unit_vec[1];
   double n = unit_vec[2];

   double ll = l*l;
   double mm = m*m;
   double nn = n*n;
   double cosk = cos(angle);
   double sink = sin(angle);
   double I_cosk = 1.0 - cosk;
   
   // The Rotation matrix angle w about vector with direction cosines l,m,n.
   // 
   // ( l**2+(m**2+n**2)cos k     lm(1-cos k)-nsin k        nl(1-cos k)+msin k   )
   // ( lm(1-cos k)+nsin k        m**2+(l**2+n**2)cos k     mn(1-cos k)-lsin k   )
   // ( nl(1-cos k)-msin k        mn(1-cos k)+lsin k        n*2+(l**2+m**2)cos k )
   //
   // (Amore documentation) Thanks for that pointer EJD :).
   
   clipper::Mat33<double> r( ll+(mm+nn)*cosk,    l*m*I_cosk-n*sink,  n*l*I_cosk+m*sink,
                             l*m*I_cosk+n*sink,  mm+(ll+nn)*cosk,    m*n*I_cosk-l*sink,
                             n*l*I_cosk-m*sink,  m*n*I_cosk+l*sink,  nn+(ll+mm)*cosk );
   clipper::RTop_orth rtop(r, clipper::Coord_orth(0,0,0));
   return origin_shift + (position-origin_shift).transform(rtop);
}

std::string
coot::util::cis_peptide_info_t::string() const {

   std::string s;
   s += chain_id_1;
   s += " ";
   s += int_to_string(resno_1);
   if (! ins_code_1.empty()) { 
      s += " ";
      s += ins_code_1;
   }
   s += " - ";
   s += chain_id_2;
   s += " ";
   s += int_to_string(resno_2);
   if (! ins_code_2.empty()) { 
      s += " ";
      s += ins_code_2;
   }
   return s;
}

clipper::Coord_orth
coot::util::average_position(std::vector<clipper::Coord_orth> &pts) {
   
   if (pts.size() > 0) {
      double xsum = 0.0;
      double ysum = 0.0;
      double zsum = 0.0;
      for (unsigned int i=0; i<pts.size(); i++) { 
         xsum += pts[i].x();
         ysum += pts[i].y();
         zsum += pts[i].z();
      }
      double denom=1.0/double(pts.size());
      return clipper::Coord_orth(denom*xsum, denom*ysum, denom*zsum);
   } else {
      return clipper::Coord_orth(0,0,0);
   }
}

double
coot::util::min_dist_to_points(const clipper::Coord_orth &pt,
                               const std::vector<clipper::Coord_orth> &others) {

   double best_dist = 9999999.9;
   for (unsigned int i=0; i<others.size(); i++) {
      double d = (pt - others[i]).lengthsq();
      if (d<best_dist) {
         best_dist = d;
      }
   }
   return sqrt(best_dist);
}

clipper::Coord_frac
coot::util::shift_to_origin(const std::vector<clipper::Coord_orth> &protein_coords,
                            clipper::Cell cell,
                            clipper::Spacegroup spacegroup) {

   clipper::Coord_orth median_pos = median_position(protein_coords);
   clipper::Coord_frac mpf = median_pos.coord_frac(cell);
   clipper::Coord_frac rf (round(-mpf.u()), round(-mpf.v()), round(-mpf.w()));
   return rf;
}

bool
coot::util::is_000_shift(const clipper::Coord_frac &cf_shift) {

   // This is for testing the results of the above function

   if (std::abs(cf_shift.u()) > 0.1) {
      return false;
   } else {
      if (std::abs(cf_shift.v()) > 0.1) {
         return false;
      } else {
         if (std::abs(cf_shift.w()) > 0.1) {
            return false;
         }
      }
   }
   return true;
}

clipper::Coord_orth
coot::util::median_position(const std::vector<clipper::Coord_orth> &pts) {

   if (pts.size() == 0 ) {
      std::string message = "No atoms in molecule - no mediain position";
      throw std::runtime_error(message);
   }
   
   std::vector<float> pts_x;
   std::vector<float> pts_y;
   std::vector<float> pts_z;
   for (unsigned int i=0; i<pts.size(); i++) {
      pts_x.push_back(pts[i].x());
      pts_y.push_back(pts[i].y());
      pts_z.push_back(pts[i].z());
   }
   std::sort(pts_x.begin(), pts_x.end());
   std::sort(pts_y.begin(), pts_y.end());
   std::sort(pts_z.begin(), pts_z.end());
   unsigned int mid_index = pts_x.size()/2;
   return clipper::Coord_orth(pts_x[mid_index], pts_y[mid_index], pts_z[mid_index]);
}

clipper::Coord_orth
coot::util::translate_close_to_origin(const clipper::Coord_orth pos,
                                      const clipper::Cell &cell) {

   clipper::Coord_frac cf = pos.coord_frac(cell);
   clipper::Coord_frac cfi(round(-cf.u()), round(-cf.v()), round(-cf.w()));
   return pos + cfi.coord_orth(cell);
} 

bool 
coot::position_residue_by_internal_coordinates::move_moving_residue() {

   bool status = 0;

   return status;
} 

