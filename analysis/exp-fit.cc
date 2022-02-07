
#include <iostream>

#include "exp-fit.hh"


// explicit
coot::exponential_fit_with_offset::exponential_fit_with_offset(const std::vector<std::pair<double, double> > &A_data) {

   auto compute_Sk_for_k = [&A_data] (unsigned int k, const std::vector<double> &S) {
                              if (k == 0) return 0.0;
                              double S_k_part_A = S[k-1];
                              double S_k_part_B = 0.5 * (A_data[k].second + A_data[k-1].second) * (A_data[k].first - A_data[k-1].first);
                              return S_k_part_A + S_k_part_B;
                           };

   auto theta_k = [] (const double &c_2, const double &x_k) {
                     return std::exp(c_2 * x_k);
                   };

   class mat_2x2 {
   public:
      mat_2x2() { identity_self(); }
      mat_2x2(const double &x_00, const double &x_01, const double &x_10, const double &x_11) {
         init();
         m[0][0] = x_00;
         m[0][1] = x_01;
         m[1][0] = x_10;
         m[1][1] = x_11;
      }
      std::vector<std::vector<double> > m;
      void init() {
         m.resize(2);
         m[0].resize(2);
         m[1].resize(2);
      }
      void identity_self() {
         init();
         m[0][0] = 1.0;
         m[0][1] = 0.0;
         m[1][0] = 0.0;
         m[1][1] = 1.0;
      }
      mat_2x2 operator*(const double &d) const {
         return mat_2x2(d * m[0][0], d * m[0][1], d * m[1][0], d * m[1][1]);
      }
      std::vector<double> operator*(const std::vector<double> &v) const {
         std::vector<double> r(2);
         r[0] = m[0][0] * v[0] + m[0][1] * v[1];
         r[1] = m[1][0] * v[0] + m[1][1] * v[1];
         return r;
      }
      mat_2x2 inverse() const {
         double det = m[0][0]*m[1][1] - m[0][1]*m[1][0];
         // this might be the wrong way around.
         mat_2x2 mm(m[1][1], -m[0][1], -m[1][0], m[0][0]);
         mat_2x2 mmi = mm * (1.0/det);
         return mmi;
      }
      mat_2x2 operator *(const mat_2x2 &o) const {
         mat_2x2 r(m[0][0] * o.m[0][0] + m[0][1] * o.m[1][0],
                   m[0][0] * o.m[0][1] + m[0][1] * o.m[1][1],
                   m[1][0] * o.m[0][0] + m[1][1] * o.m[1][0],
                   m[1][0] * o.m[0][1] + m[1][1] * o.m[1][1]);
         return r;
      }
      void print() const {
         std::cout << "    [ " << m[0][0] << "  " << m[0][1] << " ]" << std::endl;
         std::cout << "    [ " << m[1][0] << "  " << m[1][1] << " ]" << std::endl;
      }
   };

   // from the reference:
   //
   // Sigma (x_k - x_1)^2                 Sigma_A
   // Sigma (x_k - x_1)S_k                Sigma_B
   // Sigma S_k^2                         Sigma_C
   // Sigma (y_k - y_1)(x_k - x_1)        Sigma_D
   // Sigma (y_k - y_1)S_k                Sigma_E
   // Sigma theta(k)                      Sigma_F
   // Sigma theta(k)^2                    Sigma_G
   // Sigma y_k                           Sigma_H
   // Sigma y_k * theta_k                 Sigma_I

   double Sigma_A = 0.0;
   double Sigma_B = 0.0;
   double Sigma_C = 0.0;
   double Sigma_D = 0.0;
   double Sigma_E = 0.0;
   double Sigma_F = 0.0;
   double Sigma_G = 0.0;
   double Sigma_H = 0.0;
   double Sigma_I = 0.0;

   // the indexing is "off-by-one" (by design)

   std::vector<double> S_k(A_data.size(), 0.0);

   for (unsigned int k=0; k<A_data.size(); k++) {
      double S = compute_Sk_for_k(k, S_k);
      S_k[k] = S;
   }

   for (unsigned int k=0; k<A_data.size(); k++) {
      const auto &data_k = A_data[k];
      const auto &data_0 = A_data[0];
      double  x_d = data_k.first  - data_0.first;
      double  y_d = data_k.second - data_0.second;
      Sigma_A += x_d * x_d;
      Sigma_B += x_d * S_k[k];
      Sigma_C += S_k[k] * S_k[k];
      Sigma_D += x_d * y_d;
      Sigma_E += y_d * S_k[k];
   }

   mat_2x2 test_m(2, 0, 0, 2);
   mat_2x2 test_m_inverse = test_m.inverse();

   mat_2x2 m(Sigma_A, Sigma_B, Sigma_B, Sigma_C);
   mat_2x2 m_inverse = m.inverse();

   std::vector<double> A1_c2 = m.inverse() * std::vector<double>{Sigma_D, Sigma_E};
   const double &c_2 = A1_c2[1];

   for (unsigned int k=0; k<A_data.size(); k++) {
      const auto &data_k = A_data[k];

      double t = theta_k(c_2, data_k.first);
      Sigma_F += t;
      Sigma_G += t * t;
      Sigma_H += data_k.second;
      Sigma_I += t * data_k.second;
   }

   unsigned int n = A_data.size();
   mat_2x2 m2(static_cast<double>(n), Sigma_F, Sigma_F, Sigma_G);
   mat_2x2 m2_inverse = m2.inverse();
   std::vector<double> a_2_b_2 = m2.inverse() * std::vector<double>{Sigma_H, Sigma_I};

   a = a_2_b_2[0];
   b = a_2_b_2[1];
   c = c_2;

   // std::cout << "debug:: a: " << a << "  b: " << b << "  c: " << c << std::endl;

}
