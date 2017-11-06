
#ifndef ZO_RAMA_HH
#define ZO_RAMA_HH

#include <utility>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>

#define PNG_DEBUG 3
#include <png.h>

#include "coot/utils/coot-utils.hh"

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
      std::pair<realtype,realtype> df_numerical(const realtype &rama_phi, const realtype &rama_psi) const {
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
	 std::string fn = "all-non-pre-pro.tab";
	 init(fn);
      }

      // this can throw a runtime_error
      rama_table(const std::string &fn) {
	 init(fn);
      }

      void make_a_png(int n_pixels, const std::string &file_name);

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
      void read(const std::string &file_name) {
	 std::ifstream f(file_name);
	 std::string line;

	 std::cout << "INFO:: reading file " << file_name << std::endl;

	 if (f) {
	    while(std::getline(f, line)) {
	       std::vector<std::string> bits = coot::util::split_string_no_blanks(line);
	       if (bits.size() == 7) {
		  std::cout << "line: " << line << std::endl;
		  int idx_1 = coot::util::string_to_int(bits[0]);
		  int idx_2 = coot::util::string_to_int(bits[1]);
		  double A_cc = coot::util::string_to_double(bits[3]);
		  double A_cs = coot::util::string_to_double(bits[4]);
		  double A_sc = coot::util::string_to_double(bits[5]);
		  double A_ss = coot::util::string_to_double(bits[6]);
		  zo::rama_coeffs v(idx_1, idx_2, A_cc, A_cs, A_sc, A_ss);
		  rama_vec.push_back(v);
	       }
	    }
	 } else {
	    std::cout << "Warning:: file not found: " << file_name << std::endl;
	    throw std::runtime_error("Can't init zo-rama");
	 }
      }

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
      enum type_t { ALL_NON_PRE_PRO, PRE_PRO, LIV };
      std::map<type_t, rama_table> table_map;

      // actually write this at some stage
      type_t get_residue_type(const std::string &this_residue_type,
			      const std::string &next_residue_type) const;

      // fill the table_map with tables
      //
      void init() {
	 try {
	    table_map[ALL_NON_PRE_PRO] = rama_table("all-non-pre-pro.tab");
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "ERROR:: " << rte.what() << std::endl;
	 }
      }
      std::pair<realtype,realtype> df(const realtype &phi, const realtype &psi) const {
	 std::map<type_t, rama_table>::const_iterator it = table_map.find(ALL_NON_PRE_PRO);
	 return it->second.df(phi,psi);
      }

      realtype value(const realtype &phi, const realtype &psi) const {
	 std::map<type_t, rama_table>::const_iterator it = table_map.find(ALL_NON_PRE_PRO);
	 if (it != table_map.end())
	    return it->second.value(phi,psi);
	 else
	    return 0.0;
      }
   };

}

#endif /// ZO_RAMA_HH
