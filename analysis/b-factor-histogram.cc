
#include <iostream>
#include <cmath>

#include "b-factor-histogram.hh"

// set n_atoms, n_bins, b_max
coot::b_factor_histogram::b_factor_histogram(mmdb::Manager *mol) {

   init();
   b_max = -1.0;
   n_atoms = 0;

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       int n_atoms_in_residue = residue_p->GetNumberOfAtoms();
	       for (int iat=0; iat<n_atoms_in_residue; iat++) {
		  mmdb::Atom *at = residue_p->GetAtom(iat);
		  const float &b = at->tempFactor;
		  if (b >= 0.0) {
		     n_atoms++;
		     if (b > b_max) {
			b_max = b;
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (n_atoms > 0) {
      n_bins = get_n_bins(); // use b_max and n_atoms
   }

   b_vector.resize(n_bins);

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       int n_atoms_in_residue = residue_p->GetNumberOfAtoms();
	       for (int iat=0; iat<n_atoms_in_residue; iat++) {
		  mmdb::Atom *at = residue_p->GetAtom(iat);
		  const float &b = at->tempFactor;
		  if (b >= 0.0) {
		     int bin_idx = b_to_bin(b);
		     b_vector[bin_idx].push_back(b);
		  }
	       }
	    }
	 }
      }
   }
}


coot::b_factor_histogram::b_factor_histogram(mmdb::Manager *mol, int atom_selection_handle) {

   init();
   b_max = -1.0;
   n_atoms = 0;

   mmdb::Atom **atom_selection = 0;
   int n_selection_atoms;
   mol->GetSelIndex(atom_selection_handle, atom_selection, n_selection_atoms);
   for (int i=0; i<n_selection_atoms; i++) {
      mmdb::Atom *at = atom_selection[i];
      const float &b = at->tempFactor;
      if (b >= 0.0) {
	 n_atoms++;
	 if (b > b_max) {
	    b_max = b;
	 }
      }
   }

   if (n_atoms > 0) {
      n_bins = get_n_bins(); // use b_max and n_atoms
   }
   b_vector.resize(n_bins);
   for (int i=0; i<n_selection_atoms; i++) {
      mmdb::Atom *at = atom_selection[i];
      const float &b = at->tempFactor;
      if (b >= 0.0) {
	 int bin_idx = b_to_bin(b);
	 b_vector[bin_idx].push_back(b);
      }
   }

}

int
coot::b_factor_histogram::b_to_bin(const float &b) const {

   int ibin = 0;
   double f = b/static_cast<double>(b_max);
   ibin = static_cast<int>(f*n_bins);

   // std::cout << "b_to_bin() ibin "  << ibin << " for f " << f << " b_max " << b_max
   // << std::endl;

   if (ibin >= n_bins)
      ibin = n_bins - 1;
   return ibin;
}


int
coot::b_factor_histogram::get_n_bins() const {

    float f_n_bins = 40.0f;
    int nb = static_cast<int>(f_n_bins);
    return nb;

}


std::vector<std::pair<double, double> >
coot::b_factor_histogram::get_data() const {

   std::vector<std::pair<double, double> > v(n_bins);

   for (std::size_t i=0; i<b_vector.size(); i++) {
      int freq = b_vector[i].size();
      double frac = static_cast<double>(i)/static_cast<double>(n_bins);
      double b_val = b_max * frac;
      std::pair<double, double> p(b_val, freq);
      v.push_back(p);
   }
   return v;
}

std::vector<std::pair<double, double> >
coot::b_factor_histogram::get_model() const {

   int n_model_points = 100; // smooth enough?

   std::vector<std::pair<double, double> > v(n_model_points+1); // so that we capture the biggest b.

   if (false) {
      std::cout << "in get_model() with alpha_estimate " << alpha_estimate << std::endl;
      std::cout << "in get_model() with beta_estimate "  << beta_estimate << std::endl;
   }

   double sf = 1.0;
   double sum = 0.0;
   for (int i=0; i<=n_model_points; i++) {
      double frac = static_cast<double>(i)/static_cast<double>(n_model_points);
      double b = b_max * frac;
      double b_model = ig(b + shift_estimate);
      sum += b_model;
   }
   sf = 1.0/sum;

   double sf_2 = 0.5 * static_cast<double>(b_max) / static_cast<double>(get_n_bins());

   for (int i=0; i<=n_model_points; i++) {
      double frac = static_cast<double>(i)/static_cast<double>(n_model_points);
      double b = b_max * frac;
      double b_model = ig(b + shift_estimate);
      std::pair<double, double> p(b, b_model * sf * sf_2 * n_atoms);
      // std::cout << "   " << b << " " << b_model << std::endl;
      v.push_back(p);
   }
   return v;
}

double
coot::b_factor_histogram::ig(const double &x) const {

   if (x <= 0.0) {
      return 0.0;
   } else {
      double g_a    = Gamma(alpha_estimate);
      double part_1 = (std::pow(beta_estimate, alpha_estimate) * std::pow(x, -alpha_estimate-1.0))/g_a;
      double part_2 = std::exp(-beta_estimate/x);

      if (false)
	 std::cout << "  ig: for x " << x << " g_a: " << g_a << " " << part_1 << " " << part_2 << " returning "
		   << part_1 * part_2 << std::endl;

      return part_1 * part_2;
   }
}

double
coot::b_factor_histogram::Gamma(const double &b) const {

   return std::tgamma(b);
}

void
coot::b_factor_histogram::init() {

   alpha_estimate = 0.0;
   beta_estimate  = 0.0;
   shift_estimate = 0.0;

}

void
coot::b_factor_histogram::model() {

   // This is pure invese gamma modelling - not shifted inverse gamma.

   // model this bf distribution with parameters alpha, beta
   // bf is shifted from b-factor:
   // bf = b-factor + shift
   // f(bf) = beta^alpha/Gamma(alpha) * (1/)^(alpha+1) exp(-beta/bf)
   // where Gamma(alpha) is
   // Gamma(alpha) = (n-1)!

   // \mathcal{IG}(x|\alpha,\beta) = \frac{\beta^\alpha x^{-\alpha-1}}{|\Gamma{\alpha}}exp(\frac{-\beta}{x})

   // Gaussian approximation to the moments, mean (mu) and variance (nu)
   // \mu = \frac{\beta}{\alpha-1}
   // \nu = \frac{\beta^2}{(\alpha-1)^2(\alpha-2))
   // \alpha(hat) = \frac{\mu^2}{\nu} + 2
   // \beta(hat)  = \mu(\frac{\mu^2}{\nu}+1)

   double sum = 0.0;
   double sum_sq = 0.0;
   int n = 0;
   for (std::size_t ibin=0; ibin<b_vector.size(); ibin++) {
      for (std::size_t i=0; i<b_vector[ibin].size(); i++) {
	 const float &b = b_vector[ibin][i];
	 sum += b;
	 sum_sq += b*b;
	 n++;
      }
   }
   double n_d = static_cast<double>(n);
   double mu = sum/n_d;
   double nu = sum_sq/n_d - mu * mu;
   if (nu < 0.0) nu = 0;

   alpha_estimate = mu * mu / nu + 2.0;
   beta_estimate  = mu * (mu*mu/nu+1);

   // optimize_estimates(); too complicated to make work at the moment.

}

#include <fstream>

#include "kolmogorov.hh"
#include "utils/coot-utils.hh"

std::vector<double>
coot::b_factor_histogram::select_from_model() const {

   unsigned int n_points = 200;
   std::vector<double> v;

   for (;;) {
      double b_phi = b_max * static_cast<double>(util::random())/static_cast<double>(RAND_MAX);
      double pr  = ig(b_phi + shift_estimate);
      double r   = static_cast<double>(util::random())/static_cast<double>(RAND_MAX);
      bool state = false;
      if (pr > r) state = true;
      // std::cout << " select_from_model " << b_phi << " " << pr << " " << state << std::endl;
      if (pr > r) {
	 v.push_back(b_phi);
	 if (v.size() == n_points)
	    break;
      }
   }

   std::ofstream f("bfm.tab");
   for (std::size_t i=0; i<v.size(); i++) {
      f << i << " " << v[i] << "\n";
   }
   f.close();

   return v;
}

void
coot::b_factor_histogram::optimize_estimates() {

   // this doesn't work

   double alpha_orig = alpha_estimate;
   std::vector<double> data_1;
   std::vector<double> data_2;
   for (std::size_t i=0; i<b_vector.size(); i++) {
      const std::vector<float> &bv = b_vector[i];
      for (std::size_t j=0; j<bv.size(); j++) {
	 data_1.push_back(bv[j]);
      }
   }

   std::cout << "alpha_orig " << alpha_orig << std::endl;

   double l = 0.99;
   for (std::size_t i=0; i<20; i++) {
      double f = static_cast<double>(i) * 0.05;
      double m_factor = (1.0-l) + 2.0 * l * f;
      alpha_estimate = m_factor * alpha_orig;

      std::vector<std::pair<double, double> > m = get_model();
      // how do I turn a probability model into a set of data?
      // data_2 = some function of m
      data_2 = select_from_model();
      std::pair<double, double> kl_divergences = nicholls::get_KL(data_1, data_2);

      std::cout << "f " << f << " l " << l << " alpha " << alpha_estimate << " k-l div: "
		<< kl_divergences.first << " " << kl_divergences.second
		<< std::endl;
   }

   // restore alpha_estimate
   alpha_estimate = alpha_orig;
}
