// -*-c++-*-
/* coords/Cartesian.cc
 * 
 * Copyright 2002, 2003, 2004, 2005 by The University of York
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

#include <math.h>
#include <vector>

#include "cos-sin.h"

#include "Cartesian.h"


coot::Cartesian
coot::Cartesian::operator*(const float &f) const {
   return coot::Cartesian(x_*f, y_*f, z_*f);
} 


void
coot::Cartesian::operator+=(const coot::Cartesian &in) {

   x_ += in.x_;
   y_ += in.y_;
   z_ += in.z_;
}

void
coot::Cartesian::operator-=(const coot::Cartesian &in) {

   x_ -= in.x_;
   y_ -= in.y_;
   z_ -= in.z_;
}

void
coot::Cartesian::operator*=(float s) {

   x_ *= s;
   y_ *= s;
   z_ *= s;
}

void
coot::Cartesian::operator/=(float s) {

   x_ /= s;
   y_ /= s;
   z_ /= s;
}


coot::Cartesian
coot::Cartesian::by_scalar(float scale) { 

   return coot::Cartesian(x_*scale, y_*scale, z_*scale); 
}

void
coot::Cartesian::invert_z(void) {
  z_ = -z_;
}

// Really consider doubles here.
// 
float
coot::Cartesian::distance_to_line(const coot::Cartesian &front, const coot::Cartesian &back) const {

   cos_sin table;

   coot::Cartesian line_vector = back - front; 
   coot::Cartesian front_to_point = (*this) - front;

   float lva = line_vector.amplitude();
   if (lva < 0.0001) { // some arbitary small number
      // (really, we are checking that front and back are not the same)
      //
      std::cout << "There is no vector between " << front << " and " << back << std::endl;
      std::cout << "So arbitarily returning 1.0" << std::endl;
      return 1.0;
   }
	 

   float front_to_point_amp = front_to_point.amplitude();
   float cos_theta_f = cos_angle_btwn_vecs(line_vector, front_to_point);
   float sin_theta_f = table(cos_theta_f);
   float d_f = sin_theta_f * front_to_point_amp;


   // now do the same for the back vector
   // and combine the two d's weightedly.

   coot::Cartesian back_to_point = (*this) - back;
   float cos_theta_b = cos_angle_btwn_vecs(line_vector, back_to_point);
   float sin_theta_b = table(cos_theta_b);
   float d_b = sin_theta_b * back_to_point.amplitude();

   //cout << "vectors: " << line_vector << " " << front_to_point << std::endl;
   //cout << "coses: " << cos_theta_f << " " << cos_theta_b << std::endl;
   
   //cout << "distances: " << d_f << " (" << sin_theta_f << ")"
   //	<< " and " << d_b  << " (" << sin_theta_b << ")"<< std::endl;

   // cout << std::endl;

   float weighted_d = (sin_theta_f*d_f+sin_theta_b*d_b)/(sin_theta_b+sin_theta_f);

   float click_front_weight = 0.25 * front_to_point_amp/lva;
   
   return weighted_d + click_front_weight;
}


// "Reduce" vector to length one
// 
short int 
coot::Cartesian::normalize() { 

   float a = amplitude(); 

   if ( ! (a > 0.0 ) ) { 
      // 
     std::cout << "ERROR in length of vector in normalize()" << std::endl; 

     return 0;

   } else { 
      
      float ar = 1/a; 

      x_ *= ar;
      y_ *= ar;
      z_ *= ar;
   }
   return 1;
} 

coot::Cartesian::Cartesian() {

   x_ = 0.0;
   y_ = 0.0;
   z_ = 0.0;
}

// coot::Cartesian
// coot::Cartesian::operator=(const coot::Cartesian &in) {

//    x_ = in.x_;
//    y_ = in.y_;
//    z_ = in.z_;

//    return *this;
// }

std::ostream& coot::operator<<(std::ostream&s, coot::Cartesian pt) {

   s << "(" << pt.x() << "," << pt.y() << "," << pt.z() << ")";
   return s;

}

std::ofstream& coot::operator<<(std::ofstream&s, coot::Cartesian pt) {

   s << "(" << pt.x() << "," << pt.y() << "," << pt.z() << ")";
   return s;

}

std::ostream& coot::operator<<(std::ostream& s, coot::CartesianPair pair) {

   s << pair.getStart() << "->" << pair.getFinish();

   return s; 

} 

std::ofstream& coot::operator<<(std::ofstream& s, coot::CartesianPair pair) {

   // IRIX CC cannot cope with these all being on the same line.
   // 
   s << pair.getStart();
   s << "->";
   s << pair.getFinish(); 

   return s; 

} 


// not a member function
float
coot::dot_product(const coot::Cartesian &v1, const coot::Cartesian &v2) {
   
   return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

// A friend, not a member function.
//
coot::Cartesian
coot::cross_product (const coot::Cartesian &v1, const coot::Cartesian &v2) {


   // Recall: a x b = (a2b3 - a3b2)i + (a3b1 - a1b3)j + (a1b2 - a2b1)k
   // Where a = a1 i + a2 j + a3 k and
   //       b = b1 i + b2 j + b3 k

      coot::Cartesian out(v1.y()*v2.z()-v1.z()*v2.y(),v1.z()*v2.x()-v1.x()*v2.z(),v1.x()*v2.y()-v1.y()*v2.x());
   return out;
}

coot::Cartesian
coot::Cartesian::CrossProduct(const coot::Cartesian &v1, const coot::Cartesian &v2) {
   return cross_product(v1, v2);
}



double
coot::cos_angle_btwn_vecs(const coot::Cartesian &v1, const coot::Cartesian &v2) {

   // take the dot product
   double a = dot_product(v1, v2);

   double b = a/(v1.amplitude()*v2.amplitude());

   // We once got an impossible cosine due to inaccuraices (I guess)
   // in this calculation, so clamp the values.

   if (b>1.0) {
      b = 1.0;
   } else {
      if ( b< -1.0 ) {
	 b = -1.0; 
      }
   }
   return b;
}

float coot::Cartesian::amplitude(void) const {

   return sqrt(x_*x_ + y_*y_ + z_*z_);
}

// static
float
coot::Cartesian::lengthsq(const Cartesian &atom_1_pos, const Cartesian &atom_2_pos) {
   return
      (atom_1_pos.x() - atom_2_pos.x()) * (atom_1_pos.x() - atom_2_pos.x()) +
      (atom_1_pos.y() - atom_2_pos.y()) * (atom_1_pos.y() - atom_2_pos.y()) +
      (atom_1_pos.z() - atom_2_pos.z()) * (atom_1_pos.z() - atom_2_pos.z());
}


coot::Cartesian
coot::Cartesian::mid_point(const coot::Cartesian &other) const {

   coot::Cartesian out((x_+other.x_)*0.5, (y_+other.y_)*0.5, (z_+other.z_)*0.5);

   return out;
}

std::vector<coot::Cartesian>
coot::Cartesian::third_points(const coot::Cartesian &other) const {

   // 
   std::vector<coot::Cartesian> out; 

   out.push_back(coot::Cartesian((x_+2*other.x_)*0.3333, (y_+2*other.y_)*0.3333,(z_+2*other.z_)*0.3333)); 
   out.push_back(coot::Cartesian((2*x_+other.x_)*0.3333, (2*y_+other.y_)*0.3333,(2*z_+other.z_)*0.3333)); 

   return out; 

} 

// This is a poor name for this function.
int
coot::Cartesian::within_box(const coot::Cartesian &front , const coot::Cartesian &back) const {

   coot::Cartesian a = back - front;      // front->back  vector
   coot::Cartesian b = (*this) - front;   // front->point vector
   coot::Cartesian c = back - (*this);    // point->back  vector

   if ( dot_product(a,b) >= 0.0 &&
	dot_product(a,c) >= 0.0 &&
	a.amplitude() >= b.amplitude()) {

      return 1; // TRUE
   } else {
      return 0;
   }
}


 
coot::CartesianPair::CartesianPair(const coot::Cartesian &start_,
			     const coot::Cartesian &finish_) {

   start = start_;
   finish = finish_;
} 

   
coot::CartesianPair::CartesianPair(void) {

} 
   
surface_face_data
on_a_face(const coot::Cartesian &a, const coot::Cartesian &b) {

   float d = 0.001;

   if ( fabsf (a.x() - b.x()) < d ) {
      return X_FACE;
   } else {
      if ( fabsf (a.y() - b.y()) < d ) {
	 return Y_FACE;
      } else {
	 if ( fabsf (a.z() - b.z()) < d ) {
	    return Z_FACE;
	 } else {
	    return NO_FACE; // not on a face
	 }
      }
   }
}


short int
is_an_in_triangle(surface_face_data face,
		  const coot::Cartesian &a, const coot::Cartesian &b) {


   switch (face) {

   case X_FACE:
      return (a.get_x() > b.get_x()) ? 1 : 0;
      break;
   case Y_FACE: 
      return (a.get_y() > b.get_y()) ? 1 : 0;
      break;
   case Z_FACE: 
      return (a.get_z() > b.get_z()) ? 1 : 0;
      break;
   case NO_FACE:  // add this case to make the compiler happy
      return 0;
   }

   // shouldn't happen.
   return 0;

}

float
coot::Cartesian::LineLength(const coot::Cartesian &a,
		      const coot::Cartesian &b) {

   coot::Cartesian c = b - a;
   return c.amplitude();
} 

double
coot::Cartesian::Angle(const coot::Cartesian &a,
		 const coot::Cartesian &b,
		 const coot::Cartesian &c) {

   coot::Cartesian vec1 = b - a;
   coot::Cartesian vec2 = b - c;
   return acos(cos_angle_btwn_vecs(vec1, vec2));

}

// c.f. BuildCas::angle_torsion_score()
double
coot::Cartesian::DihedralAngle(const coot::Cartesian &a,
			 const coot::Cartesian &b,
			 const coot::Cartesian &c,
			 const coot::Cartesian &d) {

   
   coot::Cartesian v1 = b - a;
   coot::Cartesian v2 = c - b;
   coot::Cartesian v3 = d - c;
   coot::Cartesian v3b= c - d; // backwards

   coot::Cartesian v_prod_1 = cross_product(v1, v2); 
   coot::Cartesian v_prod_2 = cross_product(v2, v3); 

   float dot_prod  = dot_product(v_prod_1, v_prod_2); 
   
   float cos_tor = dot_prod/(v_prod_1.amplitude()*v_prod_2.amplitude()); 

   // sign
   coot::Cartesian cross_1_2 = cross_product(v_prod_1,v_prod_2); 
   
   float sign = dot_product(cross_1_2, v2); 

   // float tor = sign < 0 ? acos(cos_tor) : -acos(cos_tor); // fix to clipper convention
   float tor = sign < 0 ? -acos(cos_tor) : acos(cos_tor); 

   return tor;
} 


coot::Cartesian
coot::Cartesian::GetCartFrom3Carts_intermediate(const coot::Cartesian &Atom1,
			     const coot::Cartesian &Atom2,
			     const coot::Cartesian &Atom3,
			     double blength,
			     double angle1,
			     double angle2,
			     int chiral) { // optional arg

   return position_by_torsion(Atom1, Atom2, Atom3, angle1, angle2, blength);
}

coot::Cartesian
coot::Cartesian::GetCartFrom3Carts(const coot::Cartesian &Atom1, double blength, 
			     const coot::Cartesian &Atom2, double angle1, 
			     const coot::Cartesian &Atom3, double angle2, 
			     int chiral) {

   return position_by_torsion(Atom3, Atom2, Atom1, angle1, angle2, blength);

} 



coot::Cartesian
coot::Cartesian::position_by_torsion(const coot::Cartesian &Atom_1, const coot::Cartesian &Atom_2, const coot::Cartesian &Atom_3,
			       float theta_2, float torsion, float dist) { 

   coot::Cartesian a1a2 = Atom_2 - Atom_1; 
   coot::Cartesian a2a3 = Atom_3 - Atom_2; 
   short int istat1 = 0;
   short int istat2 = 0;

   coot::Cartesian z_r_normal = a2a3; 
   istat1 = z_r_normal.normalize(); 
   if (! istat1)
      std::cout << "ERROR vector a2a3 is 0\n";
   coot::Cartesian y_r_normal = cross_product(a1a2, a2a3);  
   istat2 = y_r_normal.normalize(); 
   if (! istat2)
      std::cout << "ERROR yr is 0\n";
   coot::Cartesian x_r_normal = cross_product(y_r_normal, z_r_normal); 
   x_r_normal.normalize(); 

   float z_r = dist*sin(theta_2 - M_PI_2); 

   float l = dist*cos(theta_2 - M_PI_2); // consider when theta_2 is less than PI_BY_2
                                          // do we need a fabs here? 
                                          // Checked.  It seems not.
   
   float x_r = l*cos(torsion); 
   float y_r = l*sin(torsion); 

   coot::Cartesian x_r_vec = x_r_normal.by_scalar(x_r); 
   coot::Cartesian y_r_vec = y_r_normal.by_scalar(y_r); 
   coot::Cartesian z_r_vec = z_r_normal.by_scalar(z_r); 

   coot::Cartesian sum_bits = x_r_vec + y_r_vec + z_r_vec; 

   coot::Cartesian Atom_4 = Atom_3 + sum_bits; 

   if ( (istat1 == 0) || (istat2 == 0))
      Atom_4 = coot::Cartesian(-999.9, 0.0, 0.0);

   return Atom_4; 
}



coot::Cartesian
coot::Cartesian::rotate_about_vector(const coot::Cartesian &direction,
				     const coot::Cartesian &origin,
				     double angle) const {

   clipper::Coord_orth dir_c(direction.x(), direction.y(), direction.z());
   clipper::Coord_orth pos_c(x(),  y(),  z());
   clipper::Coord_orth ori_c(origin.x(),    origin.y(),    origin.z());
   clipper::Coord_orth unit_vec(dir_c.unit());

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
   clipper::Coord_orth tf = (pos_c-ori_c).transform(rtop);

   return origin + Cartesian(tf.x(), tf.y(), tf.z());

}
