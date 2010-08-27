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
