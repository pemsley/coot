
#include "phi-psi-prob.hh"

coot::phi_psi_prob_t::phi_psi_prob_t(const coot::util::phi_psi_t &pp, const coot::Cartesian &pos, const ramachandrans_container_t &rc) {

   is_allowed_flag = true;
   phi_psi = pp;
   position = pos;
   const clipper::Ramachandran *rama = &rc.rama;

   if (phi_psi.residue_name() == "PRO") rama = &rc.rama_pro;
   if (phi_psi.residue_name() == "GLY") rama = &rc.rama_gly;

   // Use the right version of clipper!

   // if (phi_psi.residue_name() == "ILE" || phi_psi.residue_name() == "VAL" ) rama = &rc.rama_ileval;
   // if (phi_psi.is_pre_pro())
   // if (phi_psi.residue_name() != "GLY")
   // rama = &rc.rama_pre_pro;

   probability = rama->probability(clipper::Util::d2rad(phi_psi.phi()),
                                   clipper::Util::d2rad(phi_psi.psi()));

   if (probability < 0.002) // from src/rama-plot.hh
      is_allowed_flag = false;

}
