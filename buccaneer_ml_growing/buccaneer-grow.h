/*! \file buccaneer-grow.h buccaneer library */
/* (C) 2002-2006 Kevin Cowtan & University of York all rights reserved */


#include "buccaneer-prot.h"
#include "simplex-lib.h"

#ifndef BUCCANEER_GROW_H
#define BUCCANEER_GROW_H
#pragma once

//! Class for growing Ca chains using density
class Ca_grow {
 public:
  //! constructor: take maximum number of rasidues to add in either direction
  Ca_grow( int n_grow = 25 );
  //! grow the chain given a map and target
  bool operator() ( clipper::MiniMol& mol, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const bool& rama_grow ) const;

  //! Grow a chain at both ends
  static void grow( Ca_chain& chain, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const double& cutoff, const int& ngrow, const bool& rama_grow );
  //! Add a new Ca-group to the C-terminus of a chain by best fit to density.
  static Ca_group next_ca_group( const Ca_chain& chain, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const bool& rama_grow );
  //! Add a new Ca-group to the N-terminus of a chain by best fit to density.
  static Ca_group prev_ca_group( const Ca_chain& chain, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const bool& rama_grow );

  static void set_cpus( int cpus ) { ncpu = cpus; }
 private:
  struct Rama_ang { float phi,psi; };     //!< ramachandran angles
  struct Rama_ang1 { Rama_ang r1; };
  struct Rama_ang2 { Rama_ang r1, r2; };
  int ngrow;
  static const int max_conf1, max_conf2;
  static int ncpu;
};


//! class for growing for Ca groups
class Grow_threaded : public clipper::Thread_base {
 public:
  Grow_threaded() {}
  Grow_threaded( const std::vector<Ca_chain>& chains, const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const double& cutoff, const int& n_grow, const bool& rama_grow );
  void grow( const int& chn );
  const std::vector<Ca_chain>& result() const { return chains_; }
  //! run single or multi-threaded
  bool operator() ( int nthread = 0 );
  //! merge results from multiple threads
  void merge( const Grow_threaded& other );
 private:
  void Run();        //!< the thread 'Run' method
  static int count;  //!< Thread control parameter
  // all data required for calculation is stored in the class
  std::vector<Ca_chain> chains_;
  const clipper::Xmap<float>* xmap_;
  const LLK_map_target* llktarget_;
  std::vector<bool> done;
  double cutoff_;
  int ngrow;
  clipper::Ramachandran rama1, rama2;
  bool rama_grow_;
};


////! class for refining grown Ca groups
class Target_fn_refine_n_terminal_build : Target_fn_order_zero
{
 public:
  Target_fn_refine_n_terminal_build() {}
  Target_fn_refine_n_terminal_build( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const double& rot_step );
  ~Target_fn_refine_n_terminal_build() {}
  int num_params() const { return 4; }
  //! evaluate target function for EulerXYZr offset from rot_
  double operator() ( const std::vector<double>& args ) const;
  //! refine rotation
  Ca_group refine( const Ca_chain& chain, Ca_group& ca1, Ca_group& ca2 );
 private:
  const clipper::Xmap<float>* xmap_;
  const LLK_map_target* llktarget_;
  const clipper::Ramachandran* rama1_;
  const clipper::Ramachandran* rama2_;
  double rot_step_;
  const Ca_chain* chain_;
  double omega1_;
  double omega2_;
};


//! class for refining grown Ca groups
class Target_fn_refine_c_terminal_build : Target_fn_order_zero
{
 public:
  Target_fn_refine_c_terminal_build() {}
  Target_fn_refine_c_terminal_build( const clipper::Xmap<float>& xmap, const LLK_map_target& llktarget, const clipper::Ramachandran& rama1, const clipper::Ramachandran& rama2, const double& rot_step );
  ~Target_fn_refine_c_terminal_build() {}
  int num_params() const { return 4; }
  //! evaluate target function for EulerXYZr offset from rot_
  double operator() ( const std::vector<double>& args ) const;
  //! refine rotation
  Ca_group refine( const Ca_chain& chain, Ca_group& ca1, Ca_group& ca2 );
 private:
  const clipper::Xmap<float>* xmap_;
  const LLK_map_target* llktarget_;
  const clipper::Ramachandran* rama1_;
  const clipper::Ramachandran* rama2_;
  double rot_step_;
  const Ca_chain* chain_;
  double omega1_;
  double omega2_;
};

#endif
