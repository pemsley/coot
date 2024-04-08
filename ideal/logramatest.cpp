/*
 * ideal/lograma-test.cc
 * 
 * Copyright 2008 by The University of York
 * Author: Kevin Cowtan
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */

#include "lograma.cpp"


int main() {
   // LogRamachandran lograma( LogRamachandran::All, 2.0, true ); // KDC
   LogRamachandran lograma;

   lograma.init( LogRamachandran::All, 2.0, true );

  for ( double phi = -150; phi < -0; phi += 1 ) {
    double phir = clipper::Util::d2rad( phi );
    double psir = clipper::Util::d2rad( 0.9 * phi );
    std::cout << phir << " " << psir << " " << lograma.interp( phir, psir ) << "\n";
  }

  for ( double phi = -180; phi < 179; phi += 50 )
    for ( double psi = -180; psi < 179; psi += 50 ) {
      double phir = clipper::Util::d2rad( phi );
      double psir = clipper::Util::d2rad( psi );
      std::cout << lograma.interp( phir, psir ) << "\t";
      LogRamachandran::Lgrad lgrd = lograma.interp_grad( phir, psir );
      std::cout << lgrd.logp << "\n";
      double d1 = (lograma.interp(phir+0.01,psir)-lograma.interp(phir-0.01,psir))/0.02;
      double d2 = (lograma.interp(phir,psir+0.01)-lograma.interp(phir,psir-0.01))/0.02;
      std::cout << "dphi " << lgrd.DlogpDphi << "\t" << d1 << "\n";
      std::cout << "dpsi " << lgrd.DlogpDpsi << "\t" << d2 << "\n";
    }

}
