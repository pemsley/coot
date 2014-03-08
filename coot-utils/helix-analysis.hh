
#include <stdexcept>
#include <vector>

#include <clipper/core/clipper_types.h>
#include <clipper/core/coords.h>

#include "mini-mol/atom-quads.hh"

namespace coot {

   // All angles internal in radians here.
   // 
   class helix_params_t {
      clipper::Mat33<double> A;
      clipper::Coord_orth B;
      void calc_A();
      void calc_B();
   public:
      helix_params_t(int resno_start_in, atom_quad quad_in, double t) {
	 resno_start = resno_start_in;
	 quad = quad_in;
	 torsion = clipper::Util::d2rad(t);
	 calc_A();
	 calc_B();
      }
      int resno_start; // start of this quad
      atom_quad quad;
      double torsion;
   }; 

   class helix_params_container_t {
      // get the set of atoms starting from the given residue serial number.
      atom_quad get_quad(const std::string &atom_name, CChain *chain_p, int res_serial_no);
      CMMDBManager *mol;
	 
   public:
      helix_params_container_t() {}
      std::vector<helix_params_t> params;
      void make(CMMDBManager *mol, const std::string chain_id, int resno_start, int resno_end);
   };

}
