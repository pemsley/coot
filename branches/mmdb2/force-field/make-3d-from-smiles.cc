
#include "use-rdkit.hh"

#include <GraphMol/DistGeomHelpers/Embedder.h>

int
main(int argc, char **argv) {

   
   std::string smiles_string = "C1CCCCC1";

   if (argc > 1)
      smiles_string = argv[1];
      

   RDKit::RWMol *mol=new RDKit::RWMol();
   mol = RDKit::SmilesToMol(smiles_string);
   RDDepict::compute2DCoords(*mol);

   double vdwThresh=10.0;
   int confId = -1;
   bool ignoreInterfragInteractions=true;
   int maxIters = 200;

   int cid = RDKit::DGeomHelpers::EmbedMolecule(*mol);
   if(cid<0) {
      std::cout << "embedding failed." << std::endl;
   } else {

      ForceFields::ForceField *ff =
	 RDKit::UFF::constructForceField(*mol,
					 vdwThresh, confId,
					 ignoreInterfragInteractions);
      ff->initialize();
      int res=ff->minimize(maxIters);
      delete ff;

      std::cout << RDKit::MolToMolBlock(*mol, true, -1) << std::endl;
   }

   return 0;
}

