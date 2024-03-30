
#ifndef DISHY_BASE_H
#define DISHY_BASE_H

#include <vector>
#include <clipper/core/clipper_types.h>
#include <mmdb2/mmdb_manager.h>
#include "MoleculesToTriangles/CXXSurface/CXXCoord.h"

class DishyPlaneLSQ_t {
   std::vector<float> abcd;
public:
   FCXXCoord centre_;
   DishyPlaneLSQ_t(const std::vector<FCXXCoord> &v) {
      std::size_t n_atoms = v.size();
      if (n_atoms > 0) {
	 FCXXCoord sum;
	 for (int i=0; i<n_atoms; i++)
	    sum += v[i];
	 centre_ = sum/float(v.size());

	 clipper::Matrix<double> mat(3,3);
	 for (int i=0; i<n_atoms; i++) {
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
	 // Let's now extract the values of a,b,c normalize them
	 abcd.resize(4);
   
	 abcd[0] = mat(0,0);
	 abcd[1] = mat(1,0);
	 abcd[2] = mat(2,0);

	 double sqsum = 1e-20;
   
	 for (int i=0; i<3; i++)
	    sqsum += abcd[i] * abcd[i];
	 for (int i=0; i<3; i++)
	    abcd[i] /= sqsum;
   
	 // set D, recall di = Axi+Byi+Czi-D, so when
	 // xi = x_cen, yi = y_cen, zi = z_cen, d is 0,
	 // so we can set D.
	 // 
	 abcd[3] = abcd[0]*centre_.x() + abcd[1]*centre_.y() + abcd[2]*centre_.z();

	 double var = 0;
	 for (unsigned int i_plane_at=0; i_plane_at<v.size(); i_plane_at++) {
	    double d =
	       abcd[0]*v[i_plane_at].x() +
	       abcd[1]*v[i_plane_at].y() +
	       abcd[2]*v[i_plane_at].z() - abcd[3];
	    var += d*d;
	 }
	 float rms = 0;
	 if (v.size() > 0)
	    rms = sqrt(var/double(v.size()));
	 // std::cout << "debug:: n_atoms " << n_atoms << " rms: " << rms << std::endl;
      }
   }
   FCXXCoord normal() const {
         return FCXXCoord(abcd[0], abcd[1], abcd[2]);
   }
};

class DishyBase_t {

public:

   // color info not passed
   //
   DishyBase_t(const FCXXCoord &centre_in, const FCXXCoord &normal_in, const float &rad_in,
	       const std::vector<mmdb::Atom *> &ribose_atoms_in, const FCXXCoord &ribose_centre_in) :
      centre(centre_in), normal(normal_in), radius(rad_in), ribose_atoms(ribose_atoms_in), ribose_centre(ribose_centre_in) { idx = 0; }

   // ribose_atoms is guaranteed to be size 5 or 0 (fail, don't
   // do any nucleic acid things)
   //
   // Draw a stick between 0->1, 1-2, 2->3, 3->4, 4->0
   // and triangles fan from ribose_centre to 0,1,2,3,4,0
   // Draw a stick from ribose_atoms[1] to 1/3 of the way to
   // centre.
   std::vector<mmdb::Atom *> ribose_atoms; // in a particular order: (O4', C1', C2', C3', C4')
   FCXXCoord ribose_centre;
   int idx;
   FCXXCoord normal;
   FCXXCoord centre;
   double radius;
    static std::vector<std::pair<int, int> >bondingPattern;
};

// one of these for every segment
//
class DishyBaseContainer_t {

   void init();
public:
   DishyBaseContainer_t() { init(); }
   std::vector<DishyBase_t> bases;
   bool index_order; // so that 3'-5' order can be used as well as 5'-3'
                     // i.e. 3'-5' can use reverse coloring
   std::vector<std::string> cytidine_base_names;
   std::vector<std::string> uracil_base_names;
   std::vector<std::string> adenine_base_names;
   std::vector<std::string> guanine_base_names;
   std::vector<std::string> thymine_base_names;
   void add(const DishyBase_t &db_in) {
      bases.push_back(db_in);
   }

};


#endif // DISHY_BASE_H
