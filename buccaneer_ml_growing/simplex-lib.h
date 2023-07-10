/*! \file simplex-lib.h simplex optimiser library */
/* Copyright 2003-2006 Kevin Cowtan & University of York all rights reserved */


#ifndef SIMPLEX_LIB
#define SIMPLEX_LIB

#include <vector>


/*! Abstract base class for zero-th order function. */
class Target_fn_order_zero {
 public:
  Target_fn_order_zero() {}
  virtual ~Target_fn_order_zero() {}
  virtual int num_params() const = 0;
  virtual double operator() ( const std::vector<double>& args ) const = 0;
};


/*! Simplex optimiser. */
class Optimiser_simplex {
 public:
  enum TYPE { NORMAL, GRADIENT };
  Optimiser_simplex( double tolerance = 0.001, int max_cycles = 50, TYPE type =
NORMAL );
  std::vector<double> operator() ( const Target_fn_order_zero& target_fn, const
std::vector<std::vector<double> >& args ) const;
  void debug() { debug_mode = true; }
 private:
  double tolerance_, max_cycles_;
  TYPE type_;
  bool debug_mode;
  std::vector<double> params_;
};


#endif
