
#ifndef LIGAND_SCORED_MOLECULE_HH
#define LIGAND_SCORED_MOLECULE_HH

namespace coot {
   // Trivial class so that we can pass the (best orientation) ligand
   // back to fit_ligands_to_clusters()
   // 
   class scored_molecule {
   public:
      minimol::molecule mol;
      ligand_score_card score_card;
   };

}

#endif // LIGAND_SCORED_MOLECULE_HH
