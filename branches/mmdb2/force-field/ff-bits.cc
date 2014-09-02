
#include "use-rdkit.hh"

#include <ForceField/UFF/BondStretch.h>
#include <ForceField/UFF/AngleBend.h>


int main(int argc, char **argv) {

   ForceFields::UFF::AtomicParams p1,p2;
   double restLen,forceConstant;

   // sp3 carbon:
   p1.r1 = .757;
   p1.Z1 = 1.912;
   p1.GMP_Xi = 5.343;
  
   // sp3 - sp3: checks basics
   restLen=ForceFields::UFF::Utils::calcBondRestLength(1.0,&p1,&p1);
   TEST_ASSERT(RDKit::feq(restLen,1.514));
  
   forceConstant=ForceFields::UFF::Utils::calcBondForceConstant(restLen,&p1,&p1);
   TEST_ASSERT(RDKit::feq(forceConstant,699.5918));
  
   // sp2 carbon:
   p2.r1 = .732;
   p2.Z1 = 1.912;
   p2.GMP_Xi = 5.343;
   // sp2 - sp2: checks rBO
   restLen=ForceFields::UFF::Utils::calcBondRestLength(2.0,&p2,&p2);
   TEST_ASSERT(RDKit::feq(restLen,1.32883,1e-5));

   forceConstant=ForceFields::UFF::Utils::calcBondForceConstant(restLen,&p2,&p2);
   TEST_ASSERT(RDKit::feq(forceConstant,1034.69,1e-2));
  
   // sp3 nitrogen:
   p2.r1 = .700;
   p2.Z1 = 2.544;
   p2.GMP_Xi = 6.899;

   // Csp3 - Nsp3: checks rEN
   restLen=ForceFields::UFF::Utils::calcBondRestLength(1.0,&p1,&p2);
   TEST_ASSERT(RDKit::feq(restLen,1.462929,1e-5));

   forceConstant=ForceFields::UFF::Utils::calcBondForceConstant(restLen,&p1,&p2);
   TEST_ASSERT(RDKit::feq(forceConstant,1031.77,1e-2));


   // amide bond: check we can reproduce values from the UFF paper:
   // C_R:
   p1.r1 = .729;
   p1.Z1 = 1.912;
   p1.GMP_Xi = 5.343;
   // N_R:
   p2.r1 = .699;
   p2.Z1 = 2.544;
   p2.GMP_Xi = 6.899;

   restLen=ForceFields::UFF::Utils::calcBondRestLength(ForceFields::UFF::Params::amideBondOrder,&p1,&p2);
   TEST_ASSERT(RDKit::feq(restLen,1.368,1e-3)); // NOTE: the paper has 1.366

   forceConstant=ForceFields::UFF::Utils::calcBondForceConstant(restLen,&p1,&p2);
   TEST_ASSERT(RDKit::feq(forceConstant,1260.,1)); // NOTE: the paper has 1293

   std::cout << "restLen: " << restLen << " forceConstant: " << forceConstant << std::endl;

   return 0;

} 
