
#include <string.h> // memcpy
#include "simple-restraint.hh"
#include <gsl/gsl_siman.h>

/* some functions to test - in one dimension */
double E1(void *xp) {
   const double &x = *static_cast<double *>(xp);
   return exp(-pow((x-1.0),2.0))*sin(8.0*x);
}

double M1(void *xp, void *yp) {
   double x = *static_cast<double *>(xp);
   double y = *static_cast<double *>(yp);
   return fabs(x - y);
}

void S1(const gsl_rng *r, void *xp, double step_size) {
  double old_x = *static_cast<double *>(xp);
  double u = gsl_rng_uniform(r);
  double new_x = u * 2 * step_size - step_size + old_x;
  memcpy(xp, &new_x, sizeof(new_x));
}

void P1(void *xp) {
   printf("%12g", *static_cast<double *>(xp));
}

   /* how many points do we try before stepping */
#define N_TRIES 200
   /* how many iterations for each T? */
#define ITERS_FIXED_T 1000
   /* max step size in random walk */
#define STEP_SIZE 1.0
   /* Boltzmann constant */
#define K 1.0
   /* initial temperature */
#define T_INITIAL 0.008
   /* damping factor for temperature */
#define MU_T 1.003
#define T_MIN 2.0e-6

coot::restraints_container_t *sa_restraints;

void
coot::restraints_container_t::simulated_annealing() {
   
   gsl_rng_env_setup();
   const gsl_rng_type *T = gsl_rng_default;
   gsl_rng *r = gsl_rng_alloc(T);
   unsigned int sod = sizeof(double);

  // E1:  returns the energy of the configuration
  // S1:  modify the configuration using a random step taken from the generator r
  //      upto a max. dist of step_size
  // M1:  the distance between the 2 configurations
  // P1:  print the contents
  // not  used:
  // C1:  copy the configuration from source to dest
  // CC1: make a new copy of the configuration
  // D1:  destroy the configuration

   auto E1_m = [] (void *xp) {
                  return sa_restraints->get_distortion_score();
               };

   auto S1_m = [] (const gsl_rng *r, void *xp, double step_size) {
                  gsl_vector *x = static_cast<gsl_vector *>(xp);
                  // double new_x = u * 2 * step_size - step_size + old_x;
                  // now copy new_x into xp
                  unsigned int n_atoms = sa_restraints->get_n_atoms();
                  for (unsigned int i=0; i<n_atoms; i++) {
                     double u1 = gsl_rng_uniform(r);
                     double u2 = gsl_rng_uniform(r);
                     double u3 = gsl_rng_uniform(r);
                     *gsl_vector_ptr(x, i*3  ) += 0.1 * 2.0 * (u1 - 0.5);
                     *gsl_vector_ptr(x, i*3+1) += 0.1 * 2.0 * (u2 - 0.5);
                     *gsl_vector_ptr(x, i*3+2) += 0.1 * 2.0 * (u3 - 0.5);
                  }
               };

   auto M1_m = [] (void *xp, void *yp) {
                  gsl_vector *x = static_cast<gsl_vector *>(xp);
                  gsl_vector *y = static_cast<gsl_vector *>(xp);
                  unsigned int n_atoms = sa_restraints->get_n_atoms();
                  double sum_dist = 0;
                  for (unsigned int i=0; i<n_atoms; i++) {
                     double d1 = gsl_vector_get(y, i*3  ) - gsl_vector_get(x, i*3  );
                     double d2 = gsl_vector_get(y, i*3  ) - gsl_vector_get(x, i*3  );
                     double d3 = gsl_vector_get(y, i*3  ) - gsl_vector_get(x, i*3  );
                     double d_sqrd = d1 * d1 + d2 * d2 + d3 * d3;
                     sum_dist += std::sqrt(d_sqrd);
                  }
                  return sum_dist;
               };

   sa_restraints = this;
   gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE, K, T_INITIAL, MU_T, T_MIN};
   gsl_siman_solve(r, x, E1_m, S1_m, M1_m, P1, NULL, NULL, NULL, sod, params);
   gsl_rng_free (r);

}
