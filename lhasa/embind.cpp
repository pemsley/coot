#include "embind.hpp"
#include "lhasa.hpp"
#include "../layla/ligand_editor_canvas.hpp"
#include "../layla/utils.hpp"

EMSCRIPTEN_BINDINGS(lhasa) {
  function("remove_non_polar_hydrogens", &coot::layla::remove_non_polar_hydrogens);
  // function("rdkit_mol_from_smiles", &lhasa::rdkit_mol_from_smiles);
  // function("rdkit_mol_to_smiles", &lhasa::rdkit_mol_to_smiles);
  class_<CootLigandEditorCanvas>("LhasaCanvas")
    // remove the semicolon when adding functions
    .constructor<>();
    // .function("demo", &LhasaDemo::demo);
}