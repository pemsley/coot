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

#ifndef CARTESIAN_H
#define CARTESIAN_H

#include <fstream>
#include <vector>

#include "clipper/core/coords.h"

// #include "cos-sin.h"

enum surface_face_data { NO_FACE = 0, X_FACE = 1, Y_FACE = 2, Z_FACE = 3 };

namespace coot { 

   class Cartesian { 

      float x_, y_, z_; 

   public: 

      // inlines (of disaproved of funtions). 
      //
      float get_x() const { return x_;};
      float get_y() const { return y_;};
      float get_z() const { return z_;};
   
      // actually, I don't like the gets part of the variable name. Make const ref also.

      const float &x() const { return x_;};
      const float &y() const { return y_;};
      const float &z() const { return z_;};

      Cartesian(float xi, float yi, float zi) { x_ = xi; y_ = yi; z_ = zi; }
      Cartesian();

      Cartesian(const clipper::Coord_orth &pt) {
	 x_ = pt.x();
	 y_ = pt.y();
	 z_ = pt.z();
      } 

      void set_them(float a, float b, float c) {   // tmp function
	 x_ = a; y_ = b; z_ = c; }


      float amplitude(void) const;
      float amplitude_squared(void) const { return (x_*x_ + y_*y_ + z_*z_); }
      float length(void) const { return amplitude(); }
      short int  normalize();  // return success status 0: fails, 1: OK

      Cartesian operator+(const Cartesian &in1) const { return Cartesian(x_+in1.x_, y_+in1.y_, z_+in1.z_); }
      Cartesian operator-(const Cartesian &in1) const { return Cartesian(x_-in1.x_, y_-in1.y_, z_-in1.z_); }
      Cartesian operator*(const float &f) const;

      void operator+=(const Cartesian &);
      void operator-=(const Cartesian &);
      void operator*=(float scale);
      void operator/=(float scale);

      Cartesian by_scalar(float scale);
      void invert_z(void);

      // 20221105-PE interesting
      Cartesian& operator=(const Cartesian &in) {
	 x_ = in.x_;
	 y_ = in.y_;
	 z_ = in.z_;
	 return *this;
      }
   
      void unit_vector_yourself() {
	 float l = (*this).amplitude();
	 (*this) /= l;
      }

      Cartesian unit() const {
	 float l = amplitude();
	 return Cartesian(x()/l, y()/l, z()/l);
      }

      // rotate this point about direction, with origin at origin.
      //
      // angle in radians
      Cartesian rotate_about_vector(const coot::Cartesian &direction,
				     const coot::Cartesian &origin,
				     double angle) const;


      float distance_to_line(const Cartesian &front, const Cartesian &back) const;
      int  within_box (const Cartesian &front, const Cartesian &back) const;

      Cartesian mid_point(const Cartesian &other) const; 

      // The following functions needed for use in mgtree
      // 
      static float LineLength(const Cartesian &a,
			      const Cartesian &b);
      static double Angle(const Cartesian &a,
			  const Cartesian &b,
			  const Cartesian &c);  // radians
      static double DihedralAngle(const Cartesian &a,
				  const Cartesian &b,
				  const Cartesian &c,
				  const Cartesian &d); // radians

      /*    static Cartesian GetCartFrom3Carts(const Cartesian &Atom1, */
      /* 				      const Cartesian &Atom2, */
      /* 				      const Cartesian &Atom3, */
      /* 				      double blength, */
      /* 				      double angle1, */
      /* 				      double angle2, */
      /* 				      int chiral=0); */

      static Cartesian GetCartFrom3Carts(const Cartesian &Atom1, double blength, 
					 const Cartesian &Atom2, double angle1, 
					 const Cartesian &Atom3, double angle2, 
					 int chiral=0);

      static Cartesian GetCartFrom3Carts_intermediate(const Cartesian &Atom1,
						      const Cartesian &Atom2,
						      const Cartesian &Atom3,
						      double blength,
						      double angle1,
						      double angle2,
						      int chiral=0); 

      static Cartesian position_by_torsion(const Cartesian &Atom_1, 
					   const Cartesian &Atom_2, 
					   const Cartesian &Atom_3,
					   float theta_2, float torsion, float dist); 

      static Cartesian CrossProduct(const Cartesian &Atom_1, 
				    const Cartesian &Atom_2);

      static float lengthsq(const Cartesian &c1, const Cartesian &c2);

      std::vector<Cartesian> third_points(const Cartesian &other) const;

   

      // should be static members?
      friend double  cos_angle_btwn_vecs(const Cartesian &a, const Cartesian &b);
      friend float           dot_product(const Cartesian &a, const Cartesian &b);
      friend Cartesian     cross_product(const Cartesian &a, const Cartesian &b);
      friend surface_face_data on_a_face(const Cartesian &a, const Cartesian &b);

      friend std::ostream&  operator<<(std::ostream&, Cartesian pt);
      friend std::ofstream& operator<<(std::ofstream&, Cartesian pt);

   };
   std::ostream&  operator<<(std::ostream&, Cartesian pt);
   std::ofstream& operator<<(std::ofstream&, Cartesian pt);

        
   class CartesianPair { 

      Cartesian start, finish;

   public:
      CartesianPair(const Cartesian &start_, const Cartesian &finish_);
      CartesianPair(void);
      void extentsBox(Cartesian centre, float dist) { 
	 start.set_them(centre.get_x() - dist, 
			centre.get_y() - dist, 
			centre.get_z() - dist); 
	 finish.set_them(centre.get_x() + dist, 
			 centre.get_y() + dist, 
			 centre.get_z() + dist);
      }

      friend std::ostream&  operator<<(std::ostream&,  CartesianPair);
      friend std::ofstream& operator<<(std::ofstream&, CartesianPair);

      // We don't approve of get functions, but the alternative is that
      // Cartesian, CartesianPair and Bond_lines get merged - and I 
      // don't want to do that.
      // 
      const Cartesian &getStart()  const { return start; }
      const Cartesian &getFinish() const { return finish;}
      float amplitude() const;
   };
   std::ostream&  operator<<(std::ostream&,  CartesianPair);
   std::ofstream& operator<<(std::ofstream&, CartesianPair);

   short int is_an_in_triangle(surface_face_data face,  const Cartesian &b,
			       const Cartesian &c); 
   // declared as friends above, now declare for real
   double  cos_angle_btwn_vecs(const Cartesian &a, const Cartesian &b);
   float           dot_product(const Cartesian &a, const Cartesian &b);
   Cartesian     cross_product(const Cartesian &a, const Cartesian &b);


   class CartesianPairInfo { 

   public:

      CartesianPairInfo() {
	 data = 0;
	 size = 0;
      }
      CartesianPair *data;
      int size;
   };


} // namespace coot

#endif // CARTESIAN_H
