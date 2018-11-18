
#ifndef ZO_RAMA_HH
#define ZO_RAMA_HH

#include <utility>
#include <vector>
#include <stdexcept>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>

#define PNG_SKIP_SETJMP_CHECK true
#define PNG_DEBUG 3
#include <png.h>

#include "utils/coot-utils.hh"
#include "phi-psi.hh"

namespace zo {

   typedef double realtype;

   class rama_coeffs {
   public:
      double A_cc;
      double A_cs;
      double A_sc;
      double A_ss;
      int order_phi, order_psi;

      rama_coeffs() {}
      rama_coeffs(int order_phi_in, int order_psi_in,
		  const double &A_1, const double &A_2,
		  const double &A_3, const double &A_4) :
	 A_cc(A_1), A_cs(A_2), A_sc(A_3), A_ss(A_4),
	 order_phi(order_phi_in), order_psi(order_psi_in) { }

      realtype value(const realtype &rama_phi, const realtype &rama_psi) const {

	 // things here can be sped up because
	 //
	 // cos(2a) = cos^2(a) - sin^2(a)
	 // sin(2a) = 2cos(a)sin(a)
	 //
	 // or more generally
	 //
	 // cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
	 // sin(a+b) = cos(a)sin(b) + sin(a)cos(b)

	 double v;
	 v  = A_cc * cosf(order_phi * rama_phi) * cosf(order_psi * rama_psi);
	 v += A_cs * cosf(order_phi * rama_phi) * sinf(order_psi * rama_psi);
	 v += A_sc * sinf(order_phi * rama_phi) * cosf(order_psi * rama_psi);
	 v += A_ss * sinf(order_phi * rama_phi) * sinf(order_psi * rama_psi);

	 return v;
      }

      std::pair<realtype, realtype> df(const realtype &rama_phi, const realtype &rama_psi) const {

	 realtype d1 = 0, d2 = 0;
	 d1  = A_cc * order_phi * -sinf(order_phi * rama_phi) * cosf(order_psi * rama_psi);
	 d1 += A_cs * order_phi * -sinf(order_phi * rama_phi) * sinf(order_psi * rama_psi);
	 d1 += A_sc * order_phi *  cosf(order_phi * rama_phi) * cosf(order_psi * rama_psi);
	 d1 += A_ss * order_phi *  cosf(order_phi * rama_phi) * sinf(order_psi * rama_psi);

	 d2  = A_cc * cosf(order_phi * rama_phi) * order_psi * -sinf(order_psi * rama_psi);
	 d2 += A_cs * cosf(order_phi * rama_phi) * order_psi *  cosf(order_psi * rama_psi);
	 d2 += A_sc * sinf(order_phi * rama_phi) * order_psi * -sinf(order_psi * rama_psi);
	 d2 += A_ss * sinf(order_phi * rama_phi) * order_psi *  cosf(order_psi * rama_psi);

	 return std::pair<realtype, realtype> (d1, d2);
      }
      std::pair<realtype,realtype> df_numerical(const realtype &rama_phi,
						const realtype &rama_psi) const {
	 realtype bit = 0.01;
	 realtype rc_1 = value(rama_phi-bit, rama_psi);
	 realtype rc_2 = value(rama_phi+bit, rama_psi);
	 realtype rc_3 = value(rama_phi, rama_psi-bit);
	 realtype rc_4 = value(rama_phi, rama_psi+bit);
	 realtype d1 = (rc_2-rc_1)/(2.0 * bit);
	 realtype d2 = (rc_4-rc_3)/(2.0 * bit);
	 return std::pair<realtype, realtype> (d1, d2);
      }
   };


   const int N_COEFFS = 5;

   // can throw a runtime_error exception.
   class rama_table {

   public:
      std::vector<rama_coeffs> rama_vec;

      // this can throw a runtime_error
      rama_table() {
	 // old style single table in file
	 // std::string fn = "all-non-pre-pro.tab";
	 // init(fn);
      }

      void make_a_png(int n_pixels, const std::string &file_name) const;

      void test_analytical_derivs() const;

      realtype value(const realtype &phi_in, const realtype &psi_in) const {
	 realtype sum = 0;
	 for (std::size_t i=0; i<rama_vec.size(); i++) {
	    sum += rama_vec[i].value(phi_in, psi_in);
	 }
	 return sum;
      }

      std::pair<realtype,realtype> df(const realtype &phi_in, const realtype &psi_in) const {
	 realtype sum_1 = 0;
	 realtype sum_2 = 0;
	 for (std::size_t i=0; i<rama_vec.size(); i++) {
	    std::pair<realtype,realtype> v = rama_vec[i].df(phi_in, psi_in);
	    sum_1 += v.first;
	    sum_2 += v.second;
	 }
	 return std::pair<realtype,realtype> (sum_1, sum_2);
      }

      std::pair<realtype,realtype> df_numerical(const realtype &phi_in, const realtype &psi_in) const {
	 realtype sum_1 = 0;
	 realtype sum_2 = 0;
	 for (std::size_t i=0; i<rama_vec.size(); i++) {
	    std::pair<realtype,realtype> v = rama_vec[i].df_numerical(phi_in, psi_in);
	    sum_1 += v.first;
	    sum_2 += v.first;
	 }
	 return std::pair<realtype,realtype> (sum_1, sum_2);

      }

      // this can throw a runtime_error.
      //
      void read(const std::string &file_name);

      // this can throw a runtime_error
      void init(const std::string &local_file_name) {
	 std::string dir1 = coot::package_data_dir();
	 std::string dir2 = coot::util::append_dir_dir(dir1, "data");
	 std::string dir3 = coot::util::append_dir_dir(dir2, "rama");
	 std::string dir4 = coot::util::append_dir_dir(dir3, "zo-tables");

	 std::string full = coot::util::append_dir_file(dir4, local_file_name);
	 read(full);
      }

   };

   // void make_a_png(const rama_coeffs &rama_coeff_table[N_COEFFS][N_COEFFS],
   // int width, int height, const std::string &file_name);

   void write_png_file(int width, int height, png_bytep *row_pointers,
		       const std::string &file_name);

   // rama_table_set is a container for rama tables.
   // each rama table is for a particular residue type (e.g. pre-PRO, GLY, PRO, LIV etc)
   // this rama_table_set is a container for these different types.
   // we will select the right table based on the residue type
   // 
   class rama_table_set {
   public:
      // enum type_t { ALL_NON_PRE_PRO, PRE_PRO, LIV };
      std::map<std::string, rama_table> table_map;

      rama_table_set() { init(); }

      std::string get_residue_type(const std::string &this_residue_type,
				   const std::string &next_residue_type) const;

      // fill the table_map with tables
      //
      // should this be public?
      //
      void init();

      // residue type eg "ALL!nP" "ALLnP" "GLY!nP"  "GLYnP" "PRO!nP"
      //
      std::pair<realtype,realtype> df(const std::string &residue_type,
				      const realtype &phi, const realtype &psi) const;

      std::pair<realtype,realtype> df(const coot::phi_psi_t &pp,
				      const std::string &residue_type) const {
	 return df(residue_type, pp.phi, pp.psi);
      }

      // residue type eg "ALL!nP" "ALLnP" "GLY!nP"  "GLYnP" "PRO!nP"
      //
      realtype value(const std::string &residue_type, const realtype &phi, const realtype &psi) const;
      realtype value(const coot::phi_psi_t &pp, const std::string &residue_type) const {
	 return value(residue_type, pp.phi, pp.psi);
      }

   };
}

#endif /// ZO_RAMA_HH
