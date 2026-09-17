/*
 * api/pdbqt-charges.hh
 *
 * Per-(residue, atom) partial charges for native PDBQT receptor charge
 * assignment. The values are the Gasteiger charges that AutoDockTools/Meeko
 * prepare_receptor writes, extracted (median per residue-atom) from the AutoDock
 * Vina example receptor files, which are distributed under the Apache License 2.0.
 */

#ifndef PDBQT_CHARGES_HH
#define PDBQT_CHARGES_HH

#include <map>
#include <string>

namespace coot {

   //! @return partial charges for the 20 standard amino acids: residue name ->
   //! (atom name -> charge). Non-polar hydrogens are merged into their parent
   //! heavy atom (as in the source files), so a residue present here should not
   //! have its non-polar hydrogens merged again.
   const std::map<std::string, std::map<std::string, float> > &pdbqt_reference_charges();

}

#endif // PDBQT_CHARGES_HH
