// log-Ramachandran target for Ramachandran refinement
/*
 * 
 * Copyright 2008 The University of York
 * Author: Kevin Cowtan
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

#include <clipper/clipper.h>

using clipper::ftype;
using clipper::Prob_phi_2d;


  //! Ramachandran plot class
  /*! This class provides a reference LogRamachandran plot for Gly, Pro,
    other, and combinations of those types of residues. The source
    data comes from the best residues from the 'top500' best-determined
    structures list of D. C. and J. S. Richardson,
    http://kinemage.biochem.duke.edu/index.html

    The LogRamachandran plot is normalised by probability weighted
    Z-scoring, and then scaled by "scale" in the ctor. To combined it
    with least-squares geometry targets of the form (x-<x>)^2/sigma^2,
    use scale = 2.0 and simply add the result to the geometric
    residual.

    Setting the smooth parameter to true will fill in smoothly varying
    values for empty areas of the plot. */
  class LogRamachandran : private Prob_phi_2d
  {
  public:
    class Lgrad { public: ftype logp, DlogpDphi, DlogpDpsi; };
    //! enumeration of built-in Ramachandran tables
    enum TYPE { Gly, Pro, NonGlyPro, NonGly, All };
    //! null constructor
    LogRamachandran() {}
    //! constructor: from standard plot
    LogRamachandran( TYPE type, ftype scale = 1.0, bool fill = false );
    //! initialise: from standard plot
    void init( TYPE type, ftype scale = 1.0, bool fill = false );
     //! get log-probability for a particular pair of angles
    ftype interp( const ftype& phi, const ftype& psi ) const;
    //! get log-probability for a particular pair of angles
    Lgrad interp_grad( const ftype& phi, const ftype& psi ) const;
  };
