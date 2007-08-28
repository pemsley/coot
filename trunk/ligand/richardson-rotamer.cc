
#include "richardson-rotamer.hh"

std::vector<float>
coot::richardson_rotamer::probabilities() const {

   std::vector<float> p;
   std::string rt = Residue_Type();
   if (rt == "MSE")
      rt = "MET";
   std::vector<coot::simple_rotamer> rots = rotamers(rt, Probability_limit());

   for(unsigned int i=0; i<rots.size(); i++)
      p.push_back(rots[i].Probability_rich());
   
   return p;
}
