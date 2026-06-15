/*
 * MoleculesToTriangles/CXXClasses/DishyBase-gemmi.hh
 *
 * gemmi-native twin of DishyBase.h (mmdb->gemmi migration). Lives alongside the
 * original; types are in namespace coot::m2t so both can coexist in one build.
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 */

#ifndef DISHY_BASE_GEMMI_HH
#define DISHY_BASE_GEMMI_HH

#include <vector>
#include <string>
#include <utility>
#include <clipper/core/clipper_types.h>
#include <gemmi/model.hpp>
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

namespace coot {
   namespace m2t {

      class DishyPlaneLSQ_t {
         std::vector<float> abcd;
      public:
         FCXXCoord centre_;
         DishyPlaneLSQ_t(const std::vector<FCXXCoord> &v) {
            std::size_t n_atoms = v.size();
            if (n_atoms > 0) {
               FCXXCoord sum;
               for (std::size_t i=0; i<n_atoms; i++)
                  sum += v[i];
               centre_ = sum/float(v.size());

               clipper::Matrix<double> mat(3,3);
               for (std::size_t i=0; i<n_atoms; i++) {
                  mat(0,0) += (v[i].x() - centre_.x()) * (v[i].x() - centre_.x());
                  mat(1,1) += (v[i].y() - centre_.y()) * (v[i].y() - centre_.y());
                  mat(2,2) += (v[i].z() - centre_.z()) * (v[i].z() - centre_.z());
                  mat(0,1) += (v[i].x() - centre_.x()) * (v[i].y() - centre_.y());
                  mat(0,2) += (v[i].x() - centre_.x()) * (v[i].z() - centre_.z());
                  mat(1,2) += (v[i].y() - centre_.y()) * (v[i].z() - centre_.z());
               }
               mat(1,0) = mat(0,1);
               mat(2,0) = mat(0,2);
               mat(2,1) = mat(1,2);

               std::vector<double> eigens = mat.eigen(true);
               abcd.resize(4);
               abcd[0] = mat(0,0);
               abcd[1] = mat(1,0);
               abcd[2] = mat(2,0);

               double sqsum = 1e-20;
               for (int i=0; i<3; i++)
                  sqsum += abcd[i] * abcd[i];
               for (int i=0; i<3; i++)
                  abcd[i] /= sqsum;

               abcd[3] = abcd[0]*centre_.x() + abcd[1]*centre_.y() + abcd[2]*centre_.z();
            }
         }
         FCXXCoord normal() const {
            return FCXXCoord(-abcd[0], -abcd[1], -abcd[2]);
         }
      };

      class DishyBase_t {
      public:
         DishyBase_t(const FCXXCoord &centre_in, const FCXXCoord &normal_in, const float &rad_in,
                     const std::vector<gemmi::CRA> &ribose_atoms_in, const FCXXCoord &ribose_centre_in) :
            centre(centre_in), normal(normal_in), radius(rad_in),
            ribose_atoms(ribose_atoms_in), ribose_centre(ribose_centre_in) { idx = 0; }

         FCXXCoord centre;
         FCXXCoord normal;
         double radius;
         // CRA (not bare Atom*) so ribose atoms carry chain/residue context for colouring.
         std::vector<gemmi::CRA> ribose_atoms; // (O4', C1', C2', C3', C4')
         FCXXCoord ribose_centre;
         int idx;
         static std::vector<std::pair<int, int> > bondingPattern;
      };

      // one of these for every segment
      class DishyBaseContainer_t {
         void init();
      public:
         DishyBaseContainer_t() { init(); }
         std::vector<DishyBase_t> bases;
         bool index_order;
         std::vector<std::string> cytidine_base_names;
         std::vector<std::string> uracil_base_names;
         std::vector<std::string> adenine_base_names;
         std::vector<std::string> guanine_base_names;
         std::vector<std::string> thymine_base_names;
         void add(const DishyBase_t &db_in) { bases.push_back(db_in); }
      };
   }
}

#endif // DISHY_BASE_GEMMI_HH
