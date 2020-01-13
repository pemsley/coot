
#include "richardson-rotamer.hh"

std::vector<float>
coot::richardson_rotamer::probabilities() const {

   std::string rt = Residue_Type();
   if (rt == "MSE")
      rt = "MET";

   std::vector<simple_rotamer> rots = get_rotamers(rt, Probability_limit());

   std::vector<float> p(rots.size());
   for(unsigned int i=0; i<rots.size(); i++)
      p[i] = rots[i].Probability_rich();
   
   return p;
}
