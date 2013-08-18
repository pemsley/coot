
#include "mini-mol/atom-quads.hh"
#include <clipper/core/coords.h>


namespace coot {

   // Given a residue selection: 
   // 
   // Find pairs (like Ramachandran)
   // Return diagnostics:
   //   a vector of cablam torsion angles (and scores when we
   //   have a table to look them up).

   class cablam_pseudo_torsion_info {
   public:
      torsion_atom_quad quad;
   }; 


   class cablam {
      clipper::Coord_orth get_closest_CA_CA_approach(const coot::torsion_atom_quad &quad) const;
 public:
      cablam(PCResidue *residues, int n_sel_residues);
   };

}
