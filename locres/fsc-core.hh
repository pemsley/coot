
#ifndef FSC_CORE_HH
#define FSC_CORE_HH

// The shell sums behind a Fourier Shell Correlation.
//
// This file uses only the standard library, so that coot-utils (the FSC of
// whole maps) and libcoot-locres (the FSC of small windows of maps) can share
// it without either depending on the other.

#include <vector>
#include <complex>
#include <cmath>

namespace coot {

   //! \brief Shell sums for a Fourier Shell Correlation
   //!
   //! Add pairs of structure factors, one from each map, to the resolution
   //! shell that they belong to, then ask for the correlation in each shell:
   //!
   //! FSC = sum Re(F1 conj(F2)) / sqrt(sum |F1|^2 * sum |F2|^2)
   class fsc_shells_t {
      std::vector<double> sum_f1_f2_;   // sum of Re(F1 conj(F2))
      std::vector<double> sum_f1_sq_;   // sum of |F1|^2
      std::vector<double> sum_f2_sq_;   // sum of |F2|^2
      std::vector<unsigned int> counts_;
   public:
      //! \param n_shells the number of resolution shells
      explicit fsc_shells_t(unsigned int n_shells) :
         sum_f1_f2_(n_shells, 0.0), sum_f1_sq_(n_shells, 0.0), sum_f2_sq_(n_shells, 0.0), counts_(n_shells, 0) {}

      //! \brief Add a pair of structure factors to a shell
      //!
      //! \param shell_index the shell, which must be less than the number of shells
      //! \param f_1 the structure factor from the first map
      //! \param f_2 the structure factor for the same reflection from the second map
      //! \param n_times how many times to count this pair: 2 if it stands for its Friedel mate as well
      void add(unsigned int shell_index, const std::complex<float> &f_1, const std::complex<float> &f_2,
               unsigned int n_times = 1) {
         double a_1 = f_1.real();
         double b_1 = f_1.imag();
         double a_2 = f_2.real();
         double b_2 = f_2.imag();
         double n = static_cast<double>(n_times);
         sum_f1_f2_[shell_index] += n * (a_1 * a_2 + b_1 * b_2);
         sum_f1_sq_[shell_index] += n * (a_1 * a_1 + b_1 * b_1);
         sum_f2_sq_[shell_index] += n * (a_2 * a_2 + b_2 * b_2);
         counts_[shell_index] += n_times;
      }

      //! \return the number of pairs that have been added to this shell
      unsigned int count(unsigned int shell_index) const { return counts_[shell_index]; }
      //! \return the sum of Re(F1 conj(F2)) for this shell
      double sum_f1_f2(unsigned int shell_index) const { return sum_f1_f2_[shell_index]; }
      //! \return the sum of |F1|^2 for this shell
      double sum_f1_sq(unsigned int shell_index) const { return sum_f1_sq_[shell_index]; }
      //! \return the sum of |F2|^2 for this shell
      double sum_f2_sq(unsigned int shell_index) const { return sum_f2_sq_[shell_index]; }

      //! \brief The correlation in a shell
      //!
      //! \param shell_index the shell, which must be less than the number of shells
      //! \return the FSC for this shell, or 0 if the shell is empty or if either
      //!         map has no signal in it, so the result is never NaN
      double fsc(unsigned int shell_index) const {
         double d = sum_f1_sq_[shell_index] * sum_f2_sq_[shell_index];
         if (d > 0.0)
            return sum_f1_f2_[shell_index] / std::sqrt(d);
         return 0.0;
      }
   };
}

#endif // FSC_CORE_HH
